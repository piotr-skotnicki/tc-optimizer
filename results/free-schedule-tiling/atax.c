/**
 * This version is stamped on Apr. 14, 2015
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* atax.c: this file is part of PolyBench/C */

#include <omp.h>
#include <math.h>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

#define S3_I(i,j) tmp[i] = tmp[i] + A[i][j] * x[j];
#define S3(i,j) S3_I((i),(j))
#define S2_I(i) tmp[i] = SCALAR_VAL(0.0);
#define S2(i) S2_I((i))
#define S1_I(i) y[i] = 0;
#define S1(i) S1_I((i))
#define S4_I(i,j) y[j] = y[j] + A[i][j] * tmp[i];
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
#include "atax.h"

/* Variable declaration/allocation. */
POLYBENCH_2D_ARRAY_DECL(A, DATA_TYPE, M, N, m, n);
POLYBENCH_1D_ARRAY_DECL(x, DATA_TYPE, N, n);
POLYBENCH_1D_ARRAY_DECL(y, DATA_TYPE, N, n);
POLYBENCH_1D_ARRAY_DECL(tmp, DATA_TYPE, M, m);
  
/* Array initialization. */
static
void init_array (int m, int n,
		 DATA_TYPE POLYBENCH_2D(A,M,N,m,n),
		 DATA_TYPE POLYBENCH_1D(x,N,n))
{
  int i, j;
  DATA_TYPE fn;
  fn = (DATA_TYPE)n;

  for (i = 0; i < n; i++)
      x[i] = 1 + (i / fn);
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
      A[i][j] = (DATA_TYPE) ((i+j) % n) / (5*m);
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
    if (ii0 == 1 && ii1 >= 0 && ii2 == 2 && ii3 >= 0) {
      for (int c1 = 1024 * ii1; c1 <= min(M - 1, 1024 * ii1 + 1023); c1 += 1)
        for (int c3 = 1024 * ii3; c3 <= min(N - 1, 1024 * ii3 + 1023); c3 += 1)
          S4(c1, c3);
    } else if (ii0 == 1 && ii1 >= 0 && ii2 == 1 && ii3 >= 0) {
      for (int c1 = 1024 * ii1; c1 <= min(M - 1, 1024 * ii1 + 1023); c1 += 1)
        for (int c3 = 1024 * ii3; c3 <= min(N - 1, 1024 * ii3 + 1023); c3 += 1)
          S3(c1, c3);
    } else if (ii0 == 1 && ii1 >= 0 && ii2 == 0 && ii3 == 0) {
      for (int c1 = 1024 * ii1; c1 <= min(M - 1, 1024 * ii1 + 1023); c1 += 1)
        S2(c1);
    } else if (ii0 == 0 && ii1 >= 0 && ii2 == 0 && ii3 == 0)
      for (int c1 = 1024 * ii1; c1 <= min(N - 1, 1024 * ii1 + 1023); c1 += 1)
        S1(c1);
  }
  return 0;
}

/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_atax()
{
#pragma scop
  isl_ctx* ctx = isl_ctx_alloc();
  isl_map* rtile = isl_map_read_from_str(ctx, "[M, N] -> { [1, i1, 1, i3] -> [1, i1, 2, o3] : i1 >= 0 and 1024i1 <= -1 + M and i3 >= 0 and 1024i3 <= -1 + N and o3 >= 0 and 1024o3 <= -1 + N; [1, i1, 1, i3] -> [1, i1, 1, o3] : i1 >= 0 and 1024i1 <= -1 + M and i3 >= 0 and o3 >= 1 + i3 and 1024o3 <= -1 + N; [1, i1, 2, i3] -> [1, o1, 2, i3] : i1 >= 0 and i3 >= 0 and 1024i3 <= -1 + N and o1 >= 1 + i1 and 1024o1 <= -1 + M; [1, i1, 0, 0] -> [1, i1, 2, o3] : i1 >= 0 and 1024i1 <= -1 + M and o3 >= 0 and 1024o3 <= -1 + N; [1, i1, 0, 0] -> [1, i1, 1, o3] : i1 >= 0 and 1024i1 <= -1 + M and o3 >= 0 and 1024o3 <= -1 + N; [0, i1, 0, 0] -> [1, o1, 2, i1] : i1 >= 0 and 1024i1 <= -1 + N and o1 >= 0 and 1024o1 <= -1 + M }");
  isl_set* ii_set = isl_set_read_from_str(ctx, "[N, M] -> { [1, i1, 2, i3] : i1 >= 0 and 1024i1 <= -1 + M and i3 >= 0 and 1024i3 <= -1 + N; [1, i1, 1, i3] : i1 >= 0 and 1024i1 <= -1 + M and i3 >= 0 and 1024i3 <= -1 + N; [1, i1, 0, 0] : i1 >= 0 and 1024i1 <= -1 + M; [0, i1, 0, 0] : i1 >= 0 and 1024i1 <= -1 + N }");
  rtile = isl_map_fix_si(rtile, isl_dim_param, 0, M);
  ii_set = isl_set_fix_si(ii_set, isl_dim_param, 0, N);
  rtile = isl_map_fix_si(rtile, isl_dim_param, 1, N);
  ii_set = isl_set_fix_si(ii_set, isl_dim_param, 1, M);

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
  int m = M;
  int n = N;

  /* Initialize array(s). */
  init_array (m, n, POLYBENCH_ARRAY(A), POLYBENCH_ARRAY(x));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_atax ();

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(y)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(x);
  POLYBENCH_FREE_ARRAY(y);
  POLYBENCH_FREE_ARRAY(tmp);

  return 0;
}
