/**
 * This version is stamped on Apr. 14, 2015
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* gesummv.c: this file is part of PolyBench/C */

#include <omp.h>
#include <math.h>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

#define S3_I(i,j) tmp[i] = A[i][j] * x[j] + tmp[i];
#define S3(i,j) S3_I((i),(j))
#define S2_I(i) y[i] = SCALAR_VAL(0.0);
#define S2(i) S2_I((i))
#define S1_I(i) tmp[i] = SCALAR_VAL(0.0);
#define S1(i) S1_I((i))
#define S5_I(i) y[i] = alpha * tmp[i] + beta * y[i];
#define S5(i) S5_I((i))
#define S4_I(i,j) y[i] = B[i][j] * x[j] + y[i];
#define S4(i,j) S4_I((i),(j))

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
#include "gesummv.h"

/* Variable declaration/allocation. */
DATA_TYPE alpha;
DATA_TYPE beta;
POLYBENCH_2D_ARRAY_DECL(A, DATA_TYPE, N, N, n, n);
POLYBENCH_2D_ARRAY_DECL(B, DATA_TYPE, N, N, n, n);
POLYBENCH_1D_ARRAY_DECL(tmp, DATA_TYPE, N, n);
POLYBENCH_1D_ARRAY_DECL(x, DATA_TYPE, N, n);
POLYBENCH_1D_ARRAY_DECL(y, DATA_TYPE, N, n);
  
/* Array initialization. */
static
void init_array(int n,
		DATA_TYPE *alpha,
		DATA_TYPE *beta,
		DATA_TYPE POLYBENCH_2D(A,N,N,n,n),
		DATA_TYPE POLYBENCH_2D(B,N,N,n,n),
		DATA_TYPE POLYBENCH_1D(x,N,n))
{
  int i, j;

  *alpha = 1.5;
  *beta = 1.2;
  for (i = 0; i < n; i++)
    {
      x[i] = (DATA_TYPE)( i % n) / n;
      for (j = 0; j < n; j++) {
	A[i][j] = (DATA_TYPE) (i*j % n) / n;
	B[i][j] = (DATA_TYPE) (i*j % n) / n;
      }
    }
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_1D(y,N,n))

{
  int i;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("y");
  for (i = 0; i < n; i++) {
    if (i % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
    fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, y[i]);
  }
  POLYBENCH_DUMP_END("y");
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
    if (ii0 >= 0 && ii1 == 2 && ii2 >= 0 && ii3 == 1) {
      for (int c0 = 512 * ii0; c0 <= min(N - 1, 512 * ii0 + 511); c0 += 1)
        for (int c2 = 512 * ii2; c2 <= min(N - 1, 512 * ii2 + 511); c2 += 1)
          S4(c0, c2);
    } else if (ii0 >= 0 && ii1 == 2 && ii2 >= 0 && ii3 == 0) {
      for (int c0 = 512 * ii0; c0 <= min(N - 1, 512 * ii0 + 511); c0 += 1)
        for (int c2 = 512 * ii2; c2 <= min(N - 1, 512 * ii2 + 511); c2 += 1)
          S3(c0, c2);
    } else if (ii0 >= 0 && ii1 == 3 && ii2 == 0 && ii3 == 0) {
      for (int c0 = 512 * ii0; c0 <= min(N - 1, 512 * ii0 + 511); c0 += 1)
        S5(c0);
    } else if (ii0 >= 0 && ii1 == 1 && ii2 == 0 && ii3 == 0) {
      for (int c0 = 512 * ii0; c0 <= min(N - 1, 512 * ii0 + 511); c0 += 1)
        S2(c0);
    } else if (ii0 >= 0 && ii1 == 0 && ii2 == 0 && ii3 == 0)
      for (int c0 = 512 * ii0; c0 <= min(N - 1, 512 * ii0 + 511); c0 += 1)
        S1(c0);
  }
  return 0;
}

/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_gesummv()
{
#pragma scop
  isl_ctx* ctx = isl_ctx_alloc();
  isl_map* rtile = isl_map_read_from_str(ctx, "[N] -> { [i0, 2, i2, 1] -> [i0, 2, o2, 1] : i0 >= 0 and 512i0 <= -1 + N and i2 >= 0 and o2 >= 1 + i2 and 512o2 <= -1 + N; [i0, 2, i2, 0] -> [i0, 2, o2, 0] : i0 >= 0 and 512i0 <= -1 + N and i2 >= 0 and o2 >= 1 + i2 and 512o2 <= -1 + N; [i0, 1, 0, 0] -> [i0, 2, o2, 1] : i0 >= 0 and 512i0 <= -1 + N and o2 >= 0 and 512o2 <= -1 + N; [i0, 0, 0, 0] -> [i0, 2, o2, 0] : i0 >= 0 and 512i0 <= -1 + N and o2 >= 0 and 512o2 <= -1 + N; [i0, 2, i2, 1] -> [i0, 3, 0, 0] : i0 >= 0 and 512i0 <= -1 + N and i2 >= 0 and 512i2 <= -1 + N; [i0, 2, i2, 0] -> [i0, 3, 0, 0] : i0 >= 0 and 512i0 <= -1 + N and i2 >= 0 and 512i2 <= -1 + N; [i0, 1, 0, 0] -> [i0, 3, 0, 0] : i0 >= 0 and 512i0 <= -1 + N; [i0, 0, 0, 0] -> [i0, 3, 0, 0] : i0 >= 0 and 512i0 <= -1 + N }");
  isl_set* ii_set = isl_set_read_from_str(ctx, "[N] -> { [i0, 2, i2, 1] : i0 >= 0 and 512i0 <= -1 + N and i2 >= 0 and 512i2 <= -1 + N; [i0, 2, i2, 0] : i0 >= 0 and 512i0 <= -1 + N and i2 >= 0 and 512i2 <= -1 + N; [i0, 3, 0, 0] : i0 >= 0 and 512i0 <= -1 + N; [i0, 1, 0, 0] : i0 >= 0 and 512i0 <= -1 + N; [i0, 0, 0, 0] : i0 >= 0 and 512i0 <= -1 + N }");
  rtile = isl_map_fix_si(rtile, isl_dim_param, 0, N);
  ii_set = isl_set_fix_si(ii_set, isl_dim_param, 0, N);

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
  int n = N;

  /* Initialize array(s). */
  init_array (n, &alpha, &beta,
	      POLYBENCH_ARRAY(A),
	      POLYBENCH_ARRAY(B),
	      POLYBENCH_ARRAY(x));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_gesummv ();

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(y)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(B);
  POLYBENCH_FREE_ARRAY(tmp);
  POLYBENCH_FREE_ARRAY(x);
  POLYBENCH_FREE_ARRAY(y);

  return 0;
}
