/**
 * This version is stamped on Apr. 14, 2015
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* covariance.c: this file is part of PolyBench/C */

#include <omp.h>
#include <math.h>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

#define S8_I(i,j) cov[j][i] = cov[i][j];
#define S8(i,j) S8_I((i),(j))
#define S3_I(j) mean[j] /= float_n;
#define S3(j) S3_I((j))
#define S2_I(j,i) mean[j] += data[i][j];
#define S2(j,i) S2_I((j),(i))
#define S1_I(j) mean[j] = SCALAR_VAL(0.0);
#define S1(j) S1_I((j))
#define S7_I(i,j) cov[i][j] /= (float_n - SCALAR_VAL(1.0));
#define S7(i,j) S7_I((i),(j))
#define S6_I(i,j,k) cov[i][j] += data[k][i] * data[k][j];
#define S6(i,j,k) S6_I((i),(j),(k))
#define S5_I(i,j) cov[i][j] = SCALAR_VAL(0.0);
#define S5(i,j) S5_I((i),(j))
#define S4_I(i,j) data[i][j] -= mean[j];
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
#include "covariance.h"

/* Variable declaration/allocation. */
DATA_TYPE float_n;
POLYBENCH_2D_ARRAY_DECL(data,DATA_TYPE,N,M,n,m);
POLYBENCH_2D_ARRAY_DECL(cov,DATA_TYPE,M,M,m,m);
POLYBENCH_1D_ARRAY_DECL(mean,DATA_TYPE,M,m);

/* Array initialization. */
static
void init_array (int m, int n,
		 DATA_TYPE *float_n,
		 DATA_TYPE POLYBENCH_2D(data,N,M,n,m))
{
  int i, j;

  *float_n = (DATA_TYPE)n;

  for (i = 0; i < N; i++)
    for (j = 0; j < M; j++)
      data[i][j] = ((DATA_TYPE) i*j) / M;
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int m,
		 DATA_TYPE POLYBENCH_2D(cov,M,M,m,m))

{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("cov");
  for (i = 0; i < m; i++)
    for (j = 0; j < m; j++) {
      if ((i * m + j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
      fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, cov[i][j]);
    }
  POLYBENCH_DUMP_END("cov");
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
    if (ii0 == 2 && ii1 >= 0 && ii2 >= ii1 && ii3 == 1 && ii4 >= 0) {
      for (int c1 = 128 * ii1; c1 <= min(M - 1, 128 * ii1 + 127); c1 += 1)
        for (int c2 = max(128 * ii2, c1); c2 <= min(M - 1, 128 * ii2 + 127); c2 += 1)
          for (int c4 = 128 * ii4; c4 <= min(N - 1, 128 * ii4 + 127); c4 += 1)
            S6(c1, c2, c4);
    } else if (ii0 == 2 && ii1 >= 0 && ii2 >= ii1 && ii3 == 3 && ii4 == 0) {
      for (int c1 = 128 * ii1; c1 <= min(M - 1, 128 * ii1 + 127); c1 += 1)
        for (int c2 = max(128 * ii2, c1); c2 <= min(M - 1, 128 * ii2 + 127); c2 += 1)
          S8(c1, c2);
    } else if (ii0 == 2 && ii1 >= 0 && ii2 >= ii1 && ii3 == 2 && ii4 == 0) {
      for (int c1 = 128 * ii1; c1 <= min(M - 1, 128 * ii1 + 127); c1 += 1)
        for (int c2 = max(128 * ii2, c1); c2 <= min(M - 1, 128 * ii2 + 127); c2 += 1)
          S7(c1, c2);
    } else if (ii0 == 0 && ii1 >= 0 && ii2 == 1 && ii3 >= 0 && ii4 == 0) {
      for (int c1 = 128 * ii1; c1 <= min(M - 1, 128 * ii1 + 127); c1 += 1)
        for (int c3 = 128 * ii3; c3 <= min(N - 1, 128 * ii3 + 127); c3 += 1)
          S2(c1, c3);
    } else if (ii0 == 2 && ii1 >= 0 && ii2 >= ii1 && ii3 == 0 && ii4 == 0) {
      for (int c1 = 128 * ii1; c1 <= min(M - 1, 128 * ii1 + 127); c1 += 1)
        for (int c2 = max(128 * ii2, c1); c2 <= min(M - 1, 128 * ii2 + 127); c2 += 1)
          S5(c1, c2);
    } else if (ii0 == 1 && ii1 >= 0 && ii2 >= 0 && ii3 == 0 && ii4 == 0) {
      for (int c1 = 128 * ii1; c1 <= min(N - 1, 128 * ii1 + 127); c1 += 1)
        for (int c2 = 128 * ii2; c2 <= min(M - 1, 128 * ii2 + 127); c2 += 1)
          S4(c1, c2);
    } else if (ii0 == 0 && ii1 >= 0 && ii2 == 2 && ii3 == 0 && ii4 == 0) {
      for (int c1 = 128 * ii1; c1 <= min(M - 1, 128 * ii1 + 127); c1 += 1)
        S3(c1);
    } else if (ii0 == 0 && ii1 >= 0 && ii2 == 0 && ii3 == 0 && ii4 == 0)
      for (int c1 = 128 * ii1; c1 <= min(M - 1, 128 * ii1 + 127); c1 += 1)
        S1(c1);
  }
  return 0;
}

/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_covariance()
{
#pragma scop
  isl_ctx* ctx = isl_ctx_alloc();
  isl_map* rtile = isl_map_read_from_str(ctx, "[M, N] -> { [2, i1, i2, 1, i4] -> [2, i1, i2, 1, o4] : i1 >= 0 and i2 >= i1 and 128i2 <= -1 + M and i4 >= 0 and o4 >= 1 + i4 and 128o4 <= -1 + N; [2, i1, i2, 0, 0] -> [2, i1, i2, 1, o4] : i1 >= 0 and i2 >= i1 and 128i2 <= -1 + M and o4 >= 0 and 128o4 <= -1 + N; [1, i1, i2, 0, 0] -> [2, o1, i2, 1, i1] : i1 >= 0 and 128i1 <= -1 + N and 128i2 <= -1 + M and o1 >= 0 and o1 <= i2; [1, i1, i2, 0, 0] -> [2, i2, o2, 1, i1] : i1 >= 0 and 128i1 <= -1 + N and i2 >= 0 and o2 >= i2 and 128o2 <= -1 + M; [2, i1, i2, 1, i4] -> [2, i1, i2, 3, 0] : i1 >= 0 and i2 >= i1 and 128i2 <= -1 + M and i4 >= 0 and 128i4 <= -1 + N; [2, i1, i2, 1, i4] -> [2, i1, i2, 2, 0] : i1 >= 0 and i2 >= i1 and 128i2 <= -1 + M and i4 >= 0 and 128i4 <= -1 + N; [0, i1, 1, i3, 0] -> [0, i1, 1, o3, 0] : i1 >= 0 and 128i1 <= -1 + M and i3 >= 0 and o3 >= 1 + i3 and 128o3 <= -1 + N; [0, i1, 1, i3, 0] -> [1, o1, i1, 0, 0] : i1 >= 0 and 128i1 <= -1 + M and i3 >= 0 and 128i3 <= -1 + N and o1 >= 0 and 128o1 <= -1 + N; [2, i1, i2, 2, 0] -> [2, i1, i2, 3, 0] : i1 >= 0 and i2 >= i1 and 128i2 <= -1 + M; [2, i1, i2, 0, 0] -> [2, i1, i2, 3, 0] : i1 >= 0 and i2 >= i1 and 128i2 <= -1 + M; [2, i1, i2, 0, 0] -> [2, i1, i2, 2, 0] : i1 >= 0 and i2 >= i1 and 128i2 <= -1 + M; [0, i1, 0, 0, 0] -> [0, i1, 1, o3, 0] : i1 >= 0 and 128i1 <= -1 + M and o3 >= 0 and 128o3 <= -1 + N; [0, i1, 1, i3, 0] -> [0, i1, 2, 0, 0] : i1 >= 0 and 128i1 <= -1 + M and i3 >= 0 and 128i3 <= -1 + N; [0, i1, 2, 0, 0] -> [1, o1, i1, 0, 0] : i1 >= 0 and 128i1 <= -1 + M and o1 >= 0 and 128o1 <= -1 + N; [0, i1, 0, 0, 0] -> [1, o1, i1, 0, 0] : i1 >= 0 and 128i1 <= -1 + M and o1 >= 0 and 128o1 <= -1 + N; [0, i1, 0, 0, 0] -> [0, i1, 2, 0, 0] : i1 >= 0 and 128i1 <= -1 + M }");
  isl_set* ii_set = isl_set_read_from_str(ctx, "[M, N] -> { [2, i1, i2, 1, i4] : i1 >= 0 and 128i1 <= -1 + M and i2 >= 0 and 128i2 <= -1 + M and i4 >= 0 and 128i4 <= -1 + N; [2, i1, i2, 3, 0] : i1 >= 0 and 128i1 <= -1 + M and i2 >= 0 and 128i2 <= -1 + M; [2, i1, i2, 2, 0] : i1 >= 0 and 128i1 <= -1 + M and i2 >= 0 and 128i2 <= -1 + M; [0, i1, 1, i3, 0] : i1 >= 0 and 128i1 <= -1 + M and i3 >= 0 and 128i3 <= -1 + N; [2, i1, i2, 0, 0] : i1 >= 0 and 128i1 <= -1 + M and i2 >= 0 and 128i2 <= -1 + M; [1, i1, i2, 0, 0] : i1 >= 0 and 128i1 <= -1 + N and i2 >= 0 and 128i2 <= -1 + M; [0, i1, 2, 0, 0] : i1 >= 0 and 128i1 <= -1 + M; [0, i1, 0, 0, 0] : i1 >= 0 and 128i1 <= -1 + M }");
  rtile = isl_map_fix_si(rtile, isl_dim_param, 0, M);
  ii_set = isl_set_fix_si(ii_set, isl_dim_param, 0, M);
  rtile = isl_map_fix_si(rtile, isl_dim_param, 1, N);
  ii_set = isl_set_fix_si(ii_set, isl_dim_param, 1, N);

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
  int m = M;

  /* Initialize array(s). */
  init_array (m, n, &float_n, POLYBENCH_ARRAY(data));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_covariance ();

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(m, POLYBENCH_ARRAY(cov)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(data);
  POLYBENCH_FREE_ARRAY(cov);
  POLYBENCH_FREE_ARRAY(mean);

  return 0;
}
