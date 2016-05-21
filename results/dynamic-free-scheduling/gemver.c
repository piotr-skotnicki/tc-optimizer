/**
 * This version is stamped on Apr. 14, 2015
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* gemver.c: this file is part of PolyBench/C */

#include <omp.h>
#include <math.h>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

#define S3_I(i) x[i] = x[i] + z[i];
#define S3(i) S3_I((i))
#define S2_I(i,j) x[i] = x[i] + beta * A[j][i] * y[j];
#define S2(i,j) S2_I((i),(j))
#define S1_I(i,j) A[i][j] = A[i][j] + u1[i] * v1[j] + u2[i] * v2[j];
#define S1(i,j) S1_I((i),(j))
#define S4_I(i,j) w[i] = w[i] + alpha * A[i][j] * x[j];
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
#include "gemver.h"

/* Variable declaration/allocation. */
DATA_TYPE alpha;
DATA_TYPE beta;
POLYBENCH_2D_ARRAY_DECL(A, DATA_TYPE, N, N, n, n);
POLYBENCH_1D_ARRAY_DECL(u1, DATA_TYPE, N, n);
POLYBENCH_1D_ARRAY_DECL(v1, DATA_TYPE, N, n);
POLYBENCH_1D_ARRAY_DECL(u2, DATA_TYPE, N, n);
POLYBENCH_1D_ARRAY_DECL(v2, DATA_TYPE, N, n);
POLYBENCH_1D_ARRAY_DECL(w, DATA_TYPE, N, n);
POLYBENCH_1D_ARRAY_DECL(x, DATA_TYPE, N, n);
POLYBENCH_1D_ARRAY_DECL(y, DATA_TYPE, N, n);
POLYBENCH_1D_ARRAY_DECL(z, DATA_TYPE, N, n);

/* Array initialization. */
static
void init_array (int n,
		 DATA_TYPE *alpha,
		 DATA_TYPE *beta,
		 DATA_TYPE POLYBENCH_2D(A,N,N,n,n),
		 DATA_TYPE POLYBENCH_1D(u1,N,n),
		 DATA_TYPE POLYBENCH_1D(v1,N,n),
		 DATA_TYPE POLYBENCH_1D(u2,N,n),
		 DATA_TYPE POLYBENCH_1D(v2,N,n),
		 DATA_TYPE POLYBENCH_1D(w,N,n),
		 DATA_TYPE POLYBENCH_1D(x,N,n),
		 DATA_TYPE POLYBENCH_1D(y,N,n),
		 DATA_TYPE POLYBENCH_1D(z,N,n))
{
  int i, j;

  *alpha = 1.5;
  *beta = 1.2;

  DATA_TYPE fn = (DATA_TYPE)n;

  for (i = 0; i < n; i++)
    {   
      u1[i] = i;
      u2[i] = ((i+1)/fn)/2.0;
      v1[i] = ((i+1)/fn)/4.0;
      v2[i] = ((i+1)/fn)/6.0;
      y[i] = ((i+1)/fn)/8.0;
      z[i] = ((i+1)/fn)/9.0;
      x[i] = 0.0;
      w[i] = 0.0;
      for (j = 0; j < n; j++)
        A[i][j] = (DATA_TYPE) (i*j % n) / n;
    }   
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_1D(w,N,n))
{
  int i;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("w");
  for (i = 0; i < n; i++) {
    if (i % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
    fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, w[i]);
  }
  POLYBENCH_DUMP_END("w");
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
  isl_point_free(point);
  #pragma omp task
  {
    if (ii0 == 3 && ii1 >= 0 && ii2 >= 0) {
      for (int c1 = 512 * ii1; c1 <= min(N - 1, 512 * ii1 + 511); c1 += 1)
        for (int c2 = 512 * ii2; c2 <= min(N - 1, 512 * ii2 + 511); c2 += 1)
          S4(c1, c2);
    } else if (ii0 == 1 && ii1 >= 0 && ii2 >= 0) {
      for (int c1 = 512 * ii1; c1 <= min(N - 1, 512 * ii1 + 511); c1 += 1)
        for (int c2 = 512 * ii2; c2 <= min(N - 1, 512 * ii2 + 511); c2 += 1)
          S2(c1, c2);
    } else if (ii0 == 0 && ii1 >= 0 && ii2 >= 0) {
      for (int c1 = 512 * ii1; c1 <= min(N - 1, 512 * ii1 + 511); c1 += 1)
        for (int c2 = 512 * ii2; c2 <= min(N - 1, 512 * ii2 + 511); c2 += 1)
          S1(c1, c2);
    } else if (ii0 == 2 && ii1 >= 0 && ii2 == 0)
      for (int c1 = 512 * ii1; c1 <= min(N - 1, 512 * ii1 + 511); c1 += 1)
        S3(c1);
  }
  return 0;
}

/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_gemver()
{
#pragma scop
  isl_ctx* ctx = isl_ctx_alloc();
  isl_map* rtile = isl_map_read_from_str(ctx, "[N, alpha, beta] -> { [1, i1, i2] -> [3, o1, i1] : i1 >= 0 and 512i1 <= -1 + N and i2 >= 0 and 512i2 <= -1 + N and o1 >= 0 and 512o1 <= -1 + N; [3, i1, i2] -> [3, i1, o2] : i1 >= 0 and 512i1 <= -1 + N and i2 >= 0 and o2 >= 1 + i2 and 512o2 <= -1 + N; [1, i1, i2] -> [1, i1, o2] : i1 >= 0 and 512i1 <= -1 + N and i2 >= 0 and o2 >= 1 + i2 and 512o2 <= -1 + N; [2, i1, 0] -> [3, o1, i1] : i1 >= 0 and 512i1 <= -1 + N and o1 >= 0 and 512o1 <= -1 + N; [0, i1, i2] -> [1, i2, i1] : i1 >= 0 and 512i1 <= -1 + N and i2 >= 0 and 512i2 <= -1 + N; [0, i1, i2] -> [3, i1, i2] : i1 >= 0 and 512i1 <= -1 + N and i2 >= 0 and 512i2 <= -1 + N; [1, i1, i2] -> [2, i1, 0] : i1 >= 0 and 512i1 <= -1 + N and i2 >= 0 and 512i2 <= -1 + N }");
  isl_set* ii_set = isl_set_read_from_str(ctx, "[N, alpha, beta] -> { [3, i1, i2] : i1 >= 0 and 512i1 <= -1 + N and i2 >= 0 and 512i2 <= -1 + N; [1, i1, i2] : i1 >= 0 and 512i1 <= -1 + N and i2 >= 0 and 512i2 <= -1 + N; [0, i1, i2] : i1 >= 0 and 512i1 <= -1 + N and i2 >= 0 and 512i2 <= -1 + N; [2, i1, 0] : i1 >= 0 and 512i1 <= -1 + N }");
  rtile = isl_map_fix_si(rtile, isl_dim_param, 0, N);
  ii_set = isl_set_fix_si(ii_set, isl_dim_param, 0, N);
  rtile = isl_map_fix_si(rtile, isl_dim_param, 1, alpha);
  ii_set = isl_set_fix_si(ii_set, isl_dim_param, 1, alpha);
  rtile = isl_map_fix_si(rtile, isl_dim_param, 2, beta);
  ii_set = isl_set_fix_si(ii_set, isl_dim_param, 2, beta);

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
	      POLYBENCH_ARRAY(u1),
	      POLYBENCH_ARRAY(v1),
	      POLYBENCH_ARRAY(u2),
	      POLYBENCH_ARRAY(v2),
	      POLYBENCH_ARRAY(w),
	      POLYBENCH_ARRAY(x),
	      POLYBENCH_ARRAY(y),
	      POLYBENCH_ARRAY(z));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_gemver ();

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(w)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(u1);
  POLYBENCH_FREE_ARRAY(v1);
  POLYBENCH_FREE_ARRAY(u2);
  POLYBENCH_FREE_ARRAY(v2);
  POLYBENCH_FREE_ARRAY(w);
  POLYBENCH_FREE_ARRAY(x);
  POLYBENCH_FREE_ARRAY(y);
  POLYBENCH_FREE_ARRAY(z);

  return 0;
}
