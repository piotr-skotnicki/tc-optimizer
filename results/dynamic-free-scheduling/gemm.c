/**
 * This version is stamped on Apr. 14, 2015
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* gemm.c: this file is part of PolyBench/C */

#include <omp.h>
#include <math.h>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

#define S2_I(i,k,j) C[i][j] += alpha * A[i][k] * B[k][j];
#define S2(i,k,j) S2_I((i),(k),(j))
#define S1_I(i,j) C[i][j] *= beta;
#define S1(i,j) S1_I((i),(j))

#include <isl/ctx.h>
#include <isl/space.h>
#include <isl/point.h>
#include <isl/val.h>
#include <isl/set.h>
#include <isl/map.h>

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "gemm.h"

/* Variable declaration/allocation. */
DATA_TYPE alpha;
DATA_TYPE beta;
POLYBENCH_2D_ARRAY_DECL(C,DATA_TYPE,NI,NJ,ni,nj);
POLYBENCH_2D_ARRAY_DECL(A,DATA_TYPE,NI,NK,ni,nk);
POLYBENCH_2D_ARRAY_DECL(B,DATA_TYPE,NK,NJ,nk,nj);

/* Array initialization. */
static
void init_array(int ni, int nj, int nk,
		DATA_TYPE *alpha,
		DATA_TYPE *beta,
		DATA_TYPE POLYBENCH_2D(C,NI,NJ,ni,nj),
		DATA_TYPE POLYBENCH_2D(A,NI,NK,ni,nk),
		DATA_TYPE POLYBENCH_2D(B,NK,NJ,nk,nj))
{
  int i, j;

  *alpha = 1.5;
  *beta = 1.2;
  for (i = 0; i < ni; i++)
    for (j = 0; j < nj; j++)
      C[i][j] = (DATA_TYPE) (i*j % ni) / ni;
  for (i = 0; i < ni; i++)
    for (j = 0; j < nk; j++)
      A[i][j] = (DATA_TYPE) (i*(j+1) % nk) / nk;
  for (i = 0; i < nk; i++)
    for (j = 0; j < nj; j++)
      B[i][j] = (DATA_TYPE) (i*(j+2) % nj) / nj;
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int ni, int nj,
		 DATA_TYPE POLYBENCH_2D(C,NI,NJ,ni,nj))
{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("C");
  for (i = 0; i < ni; i++)
    for (j = 0; j < nj; j++) {
	if ((i * ni + j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
	fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, C[i][j]);
    }
  POLYBENCH_DUMP_END("C");
  POLYBENCH_DUMP_FINISH;
}

int create_task(isl_point* point, void* user) {
  isl_val* p0 = isl_point_get_coordinate_val(point, isl_dim_set, 0);
  int ii0 = isl_val_get_num_si(p0);
  isl_val_free(p0);
  isl_val* p1 = isl_point_get_coordinate_val(point, isl_dim_set, 1);
  int ii1 = isl_val_get_num_si(p1);
  isl_val_free(p1);
  isl_val* p2 = isl_point_get_coordinate_val(point, isl_dim_set, 2);
  int ii2 = isl_val_get_num_si(p2);
  isl_val_free(p2);
  isl_val* p3 = isl_point_get_coordinate_val(point, isl_dim_set, 3);
  int ii3 = isl_val_get_num_si(p3);
  isl_val_free(p3);
  isl_point_free(point);
  #pragma omp task
  {
    if (ii0 >= 0 && ii1 == 1 && ii2 >= 0 && ii3 >= 0) {
      for (int c0 = 256 * ii0; c0 <= min(NI - 1, 256 * ii0 + 255); c0 += 1)
        for (int c2 = 256 * ii2; c2 <= min(NK - 1, 256 * ii2 + 255); c2 += 1)
          for (int c3 = 256 * ii3; c3 <= min(NJ - 1, 256 * ii3 + 255); c3 += 1)
            S2(c0, c2, c3);
    } else if (ii0 >= 0 && ii1 == 0 && ii2 >= 0 && ii3 == 0)
      for (int c0 = 256 * ii0; c0 <= min(NI - 1, 256 * ii0 + 255); c0 += 1)
        for (int c2 = 256 * ii2; c2 <= min(NJ - 1, 256 * ii2 + 255); c2 += 1)
          S1(c0, c2);
  }
  return 0;
}

/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_gemm()
{
#pragma scop
  isl_ctx* ctx = isl_ctx_alloc();
  isl_map* rtile = isl_map_read_from_str(ctx, "[NI, NK, NJ] -> { [i0, 1, i2, i3] -> [i0, 1, o2, i3] : i0 >= 0 and 256i0 <= -1 + NI and i2 >= 0 and i3 >= 0 and 256i3 <= -1 + NJ and o2 >= 1 + i2 and 256o2 <= -1 + NK; [i0, 0, i2, 0] -> [i0, 1, o2, i2] : i0 >= 0 and 256i0 <= -1 + NI and i2 >= 0 and 256i2 <= -1 + NJ and o2 >= 0 and 256o2 <= -1 + NK }");
  isl_set* ii_set = isl_set_read_from_str(ctx, "[NJ, NI, NK] -> { [i0, 1, i2, i3] : i0 >= 0 and 256i0 <= -1 + NI and i2 >= 0 and 256i2 <= -1 + NK and i3 >= 0 and 256i3 <= -1 + NJ; [i0, 0, i2, 0] : i0 >= 0 and 256i0 <= -1 + NI and i2 >= 0 and 256i2 <= -1 + NJ }");
  rtile = isl_map_fix_si(rtile, isl_dim_param, 0, NI);
  ii_set = isl_set_fix_si(ii_set, isl_dim_param, 0, NJ);
  rtile = isl_map_fix_si(rtile, isl_dim_param, 1, NK);
  ii_set = isl_set_fix_si(ii_set, isl_dim_param, 1, NI);
  rtile = isl_map_fix_si(rtile, isl_dim_param, 2, NJ);
  ii_set = isl_set_fix_si(ii_set, isl_dim_param, 2, NK);

  isl_set* uds = isl_set_subtract(isl_map_domain(isl_map_copy(rtile)), isl_map_range(isl_map_copy(rtile)));
  
  #pragma omp parallel
  #pragma omp single
  {
    while (!isl_set_is_empty(uds)) {
      isl_set_foreach_point(uds, &create_task, NULL);
      #pragma omp taskwait

      ii_set = isl_set_subtract(ii_set, uds);
      rtile = isl_map_intersect_domain(rtile, isl_set_copy(ii_set));
      uds = isl_set_subtract(isl_map_domain(isl_map_copy(rtile)), isl_map_range(isl_map_copy(rtile)));
    }
    
    isl_set_foreach_point(ii_set, &create_task, NULL);
    #pragma omp taskwait
  }
  
  isl_set_free(uds);
  isl_set_free(ii_set);
  isl_map_free(rtile);
  isl_ctx_free(ctx);
#pragma endscop
}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int ni = NI;
  int nj = NJ;
  int nk = NK;

  /* Initialize array(s). */
  init_array (ni, nj, nk, &alpha, &beta,
	      POLYBENCH_ARRAY(C),
	      POLYBENCH_ARRAY(A),
	      POLYBENCH_ARRAY(B));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_gemm ();

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(ni, nj,  POLYBENCH_ARRAY(C)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(C);
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(B);

  return 0;
}
