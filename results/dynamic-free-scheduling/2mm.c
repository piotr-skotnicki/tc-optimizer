/**
 * This version is stamped on Apr. 14, 2015
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* 2mm.c: this file is part of PolyBench/C */

#include <omp.h>
#include <math.h>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

#define S3_I(i,j) D[i][j] *= beta;
#define S3(i,j) S3_I((i),(j))
#define S2_I(i,j,k) tmp[i][j] += alpha * A[i][k] * B[k][j];
#define S2(i,j,k) S2_I((i),(j),(k))
#define S1_I(i,j) tmp[i][j] = SCALAR_VAL(0.0);
#define S1(i,j) S1_I((i),(j))
#define S4_I(i,j,k) D[i][j] += tmp[i][k] * C[k][j];
#define S4(i,j,k) S4_I((i),(j),(k))

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
#include "2mm.h"

/* Variable declaration/allocation. */
DATA_TYPE alpha;
DATA_TYPE beta;
POLYBENCH_2D_ARRAY_DECL(tmp,DATA_TYPE,NI,NJ,ni,nj);
POLYBENCH_2D_ARRAY_DECL(A,DATA_TYPE,NI,NK,ni,nk);
POLYBENCH_2D_ARRAY_DECL(B,DATA_TYPE,NK,NJ,nk,nj);
POLYBENCH_2D_ARRAY_DECL(C,DATA_TYPE,NJ,NL,nj,nl);
POLYBENCH_2D_ARRAY_DECL(D,DATA_TYPE,NI,NL,ni,nl);

/* Array initialization. */
static
void init_array(int ni, int nj, int nk, int nl,
		DATA_TYPE *alpha,
		DATA_TYPE *beta,
		DATA_TYPE POLYBENCH_2D(A,NI,NK,ni,nk),
		DATA_TYPE POLYBENCH_2D(B,NK,NJ,nk,nj),
		DATA_TYPE POLYBENCH_2D(C,NJ,NL,nj,nl),
		DATA_TYPE POLYBENCH_2D(D,NI,NL,ni,nl))
{
  int i, j;

  *alpha = 1.5;
  *beta = 1.2;
  for (i = 0; i < ni; i++)
    for (j = 0; j < nk; j++)
      A[i][j] = (DATA_TYPE) (i*j % ni) / ni;
  for (i = 0; i < nk; i++)
    for (j = 0; j < nj; j++)
      B[i][j] = (DATA_TYPE) (i*(j+1) % nj) / nj;
  for (i = 0; i < nj; i++)
    for (j = 0; j < nl; j++)
      C[i][j] = (DATA_TYPE) (i*(j+3) % nl) / nl;
  for (i = 0; i < ni; i++)
    for (j = 0; j < nl; j++)
      D[i][j] = (DATA_TYPE) (i*(j+2) % nk) / nk;
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int ni, int nl,
		 DATA_TYPE POLYBENCH_2D(D,NI,NL,ni,nl))
{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("D");
  for (i = 0; i < ni; i++)
    for (j = 0; j < nl; j++) {
	if ((i * ni + j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
	fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, D[i][j]);
    }
  POLYBENCH_DUMP_END("D");
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
  isl_val* p4 = isl_point_get_coordinate_val(point, isl_dim_set, 4);
  int ii4 = isl_val_get_num_si(p4);
  isl_val_free(p4);
  isl_point_free(point);
  #pragma omp task
  {
    if (ii0 == 1 && ii1 >= 0 && ii2 >= 0 && ii3 == 1 && ii4 >= 0) {
      for (int c1 = 256 * ii1; c1 <= min(NI - 1, 256 * ii1 + 255); c1 += 1)
        for (int c2 = 256 * ii2; c2 <= min(NL - 1, 256 * ii2 + 255); c2 += 1)
          for (int c4 = 256 * ii4; c4 <= min(NJ - 1, 256 * ii4 + 255); c4 += 1)
            S4(c1, c2, c4);
    } else if (ii0 == 0 && ii1 >= 0 && ii2 >= 0 && ii3 == 1 && ii4 >= 0) {
      for (int c1 = 256 * ii1; c1 <= min(NI - 1, 256 * ii1 + 255); c1 += 1)
        for (int c2 = 256 * ii2; c2 <= min(NJ - 1, 256 * ii2 + 255); c2 += 1)
          for (int c4 = 256 * ii4; c4 <= min(NK - 1, 256 * ii4 + 255); c4 += 1)
            S2(c1, c2, c4);
    } else if (ii0 == 1 && ii1 >= 0 && ii2 >= 0 && ii3 == 0 && ii4 == 0) {
      for (int c1 = 256 * ii1; c1 <= min(NI - 1, 256 * ii1 + 255); c1 += 1)
        for (int c2 = 256 * ii2; c2 <= min(NL - 1, 256 * ii2 + 255); c2 += 1)
          S3(c1, c2);
    } else if (ii0 == 0 && ii1 >= 0 && ii2 >= 0 && ii3 == 0 && ii4 == 0)
      for (int c1 = 256 * ii1; c1 <= min(NI - 1, 256 * ii1 + 255); c1 += 1)
        for (int c2 = 256 * ii2; c2 <= min(NJ - 1, 256 * ii2 + 255); c2 += 1)
          S1(c1, c2);
  }
  return 0;
}

/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_2mm()
{
#pragma scop
  isl_ctx* ctx = isl_ctx_alloc();
  isl_map* rtile = isl_map_read_from_str(ctx, "[NI, NL, NJ, NK] -> { [1, i1, i2, 1, i4] -> [1, i1, i2, 1, o4] : i1 >= 0 and 256i1 <= -1 + NI and i2 >= 0 and 256i2 <= -1 + NL and i4 >= 0 and o4 >= 1 + i4 and 256o4 <= -1 + NJ; [0, i1, i2, 1, i4] -> [0, i1, i2, 1, o4] : i1 >= 0 and 256i1 <= -1 + NI and i2 >= 0 and 256i2 <= -1 + NJ and i4 >= 0 and o4 >= 1 + i4 and 256o4 <= -1 + NK; [0, i1, i2, 1, i4] -> [1, i1, o2, 1, i2] : i1 >= 0 and 256i1 <= -1 + NI and i2 >= 0 and 256i2 <= -1 + NJ and i4 >= 0 and 256i4 <= -1 + NK and o2 >= 0 and 256o2 <= -1 + NL; [1, i1, i2, 0, 0] -> [1, i1, i2, 1, o4] : i1 >= 0 and 256i1 <= -1 + NI and i2 >= 0 and 256i2 <= -1 + NL and o4 >= 0 and 256o4 <= -1 + NJ; [0, i1, i2, 0, 0] -> [0, i1, i2, 1, o4] : i1 >= 0 and 256i1 <= -1 + NI and i2 >= 0 and 256i2 <= -1 + NJ and o4 >= 0 and 256o4 <= -1 + NK; [0, i1, i2, 0, 0] -> [1, i1, o2, 1, i2] : i1 >= 0 and 256i1 <= -1 + NI and i2 >= 0 and 256i2 <= -1 + NJ and o2 >= 0 and 256o2 <= -1 + NL }");
  isl_set* ii_set = isl_set_read_from_str(ctx, "[NJ, NI, NL, NK] -> { [1, i1, i2, 1, i4] : i1 >= 0 and 256i1 <= -1 + NI and i2 >= 0 and 256i2 <= -1 + NL and i4 >= 0 and 256i4 <= -1 + NJ; [0, i1, i2, 1, i4] : i1 >= 0 and 256i1 <= -1 + NI and i2 >= 0 and 256i2 <= -1 + NJ and i4 >= 0 and 256i4 <= -1 + NK; [1, i1, i2, 0, 0] : i1 >= 0 and 256i1 <= -1 + NI and i2 >= 0 and 256i2 <= -1 + NL; [0, i1, i2, 0, 0] : i1 >= 0 and 256i1 <= -1 + NI and i2 >= 0 and 256i2 <= -1 + NJ }");
  rtile = isl_map_fix_si(rtile, isl_dim_param, 0, NI);
  ii_set = isl_set_fix_si(ii_set, isl_dim_param, 0, NJ);
  rtile = isl_map_fix_si(rtile, isl_dim_param, 1, NL);
  ii_set = isl_set_fix_si(ii_set, isl_dim_param, 1, NI);
  rtile = isl_map_fix_si(rtile, isl_dim_param, 2, NJ);
  ii_set = isl_set_fix_si(ii_set, isl_dim_param, 2, NL);
  rtile = isl_map_fix_si(rtile, isl_dim_param, 3, NK);
  ii_set = isl_set_fix_si(ii_set, isl_dim_param, 3, NK);

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
  int nl = NL;

  /* Initialize array(s). */
  init_array (ni, nj, nk, nl, &alpha, &beta,
	      POLYBENCH_ARRAY(A),
	      POLYBENCH_ARRAY(B),
	      POLYBENCH_ARRAY(C),
	      POLYBENCH_ARRAY(D));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_2mm ();

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(ni, nl,  POLYBENCH_ARRAY(D)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(tmp);
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(B);
  POLYBENCH_FREE_ARRAY(C);
  POLYBENCH_FREE_ARRAY(D);

  return 0;
}

