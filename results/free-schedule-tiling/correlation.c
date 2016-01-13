/**
 * This version is stamped on Apr. 14, 2015
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* correlation.c: this file is part of PolyBench/C */

#include <omp.h>
#include <math.h>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

#define S9_I(i,j) data[i][j] -= mean[j];
#define S9(i,j) S9_I((i),(j))
#define S8_I(j) stddev[j] = stddev[j] <= eps ? SCALAR_VAL(1.0) : stddev[j];
#define S8(j) S8_I((j))
#define S3_I(j) mean[j] /= float_n;
#define S3(j) S3_I((j))
#define S2_I(j,i) mean[j] += data[i][j];
#define S2(j,i) S2_I((j),(i))
#define S1_I(j) mean[j] = SCALAR_VAL(0.0);
#define S1(j) S1_I((j))
#define S7_I(j) stddev[j] = SQRT_FUN(stddev[j]);
#define S7(j) S7_I((j))
#define S6_I(j) stddev[j] /= float_n;
#define S6(j) S6_I((j))
#define S5_I(j,i) stddev[j] += (data[i][j] - mean[j]) * (data[i][j] - mean[j]);
#define S5(j,i) S5_I((j),(i))
#define S4_I(j) stddev[j] = SCALAR_VAL(0.0);
#define S4(j) S4_I((j))
#define S13_I(i,j,k) corr[i][j] += (data[k][i] * data[k][j]);
#define S13(i,j,k) S13_I((i),(j),(k))
#define S12_I(i,j) corr[i][j] = SCALAR_VAL(0.0);
#define S12(i,j) S12_I((i),(j))
#define S11_I(i) corr[i][i] = SCALAR_VAL(1.0);
#define S11(i) S11_I((i))
#define S10_I(i,j) data[i][j] /= SQRT_FUN(float_n) * stddev[j];
#define S10(i,j) S10_I((i),(j))
#define S15_I() corr[M-1][M-1] = SCALAR_VAL(1.0);
#define S15() S15_I()
#define S14_I(i,j) corr[j][i] = corr[i][j];
#define S14(i,j) S14_I((i),(j))

#include <isl/ctx.h>
#include <isl/space.h>
#include <isl/point.h>
#include <isl/val.h>
#include <isl/set.h>
#include <isl/map.h>

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
#include "correlation.h"

/* Variable declaration/allocation. */
DATA_TYPE float_n;
POLYBENCH_2D_ARRAY_DECL(data,DATA_TYPE,N,M,n,m);
POLYBENCH_2D_ARRAY_DECL(corr,DATA_TYPE,M,M,m,m);
POLYBENCH_1D_ARRAY_DECL(mean,DATA_TYPE,M,m);
POLYBENCH_1D_ARRAY_DECL(stddev,DATA_TYPE,M,m);
DATA_TYPE eps = SCALAR_VAL(0.1);
    
/* Array initialization. */
static
void init_array (int m,
		 int n,
		 DATA_TYPE *float_n,
		 DATA_TYPE POLYBENCH_2D(data,N,M,n,m))
{
  int i, j;

  *float_n = (DATA_TYPE)N;

  for (i = 0; i < N; i++)
    for (j = 0; j < M; j++)
      data[i][j] = (DATA_TYPE)(i*j)/M + i;
    
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int m,
		 DATA_TYPE POLYBENCH_2D(corr,M,M,m,m))

{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("corr");
  for (i = 0; i < m; i++)
    for (j = 0; j < m; j++) {
      if ((i * m + j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
      fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, corr[i][j]);
    }
  POLYBENCH_DUMP_END("corr");
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
  isl_val* p5 = isl_point_get_coordinate_val(point, isl_dim_set, 5);
  int ii5 = isl_val_get_num_si(p5);
  isl_val_free(p5);
  isl_point_free(point);
  #pragma omp task
  {
    if (ii0 == 3 && ii1 >= 0 && ii2 == 1 && ii3 >= ii1 && ii4 == 1 && ii5 >= 0) {
      for (int c1 = 128 * ii1; c1 <= min(M - 2, 128 * ii1 + 127); c1 += 1)
        for (int c3 = max(128 * ii3 + 1, c1 + 1); c3 <= min(M - 1, 128 * ii3 + 128); c3 += 1)
          for (int c5 = 128 * ii5; c5 <= min(N - 1, 128 * ii5 + 127); c5 += 1)
            S13(c1, c3, c5);
    } else if (ii0 == 3 && ii1 >= 0 && ii2 == 1 && ii3 >= ii1 && ii4 == 2 && ii5 == 0) {
      for (int c1 = 128 * ii1; c1 <= min(M - 2, 128 * ii1 + 127); c1 += 1)
        for (int c3 = max(128 * ii3 + 1, c1 + 1); c3 <= min(M - 1, 128 * ii3 + 128); c3 += 1)
          S14(c1, c3);
    } else if (ii0 == 3 && ii1 >= 0 && ii2 == 1 && ii3 >= ii1 && ii4 == 0 && ii5 == 0) {
      for (int c1 = 128 * ii1; c1 <= min(M - 2, 128 * ii1 + 127); c1 += 1)
        for (int c3 = max(128 * ii3 + 1, c1 + 1); c3 <= min(M - 1, 128 * ii3 + 128); c3 += 1)
          S12(c1, c3);
    } else if (ii0 == 1 && ii1 >= 0 && ii2 == 1 && ii3 >= 0 && ii4 == 0 && ii5 == 0) {
      for (int c1 = 128 * ii1; c1 <= min(M - 1, 128 * ii1 + 127); c1 += 1)
        for (int c3 = 128 * ii3; c3 <= min(N - 1, 128 * ii3 + 127); c3 += 1)
          S5(c1, c3);
    } else if (ii0 == 0 && ii1 >= 0 && ii2 == 1 && ii3 >= 0 && ii4 == 0 && ii5 == 0) {
      for (int c1 = 128 * ii1; c1 <= min(M - 1, 128 * ii1 + 127); c1 += 1)
        for (int c3 = 128 * ii3; c3 <= min(N - 1, 128 * ii3 + 127); c3 += 1)
          S2(c1, c3);
    } else if (ii0 == 2 && ii1 >= 0 && ii2 >= 0 && ii3 == 0 && ii4 == 0 && ii5 == 0) {
      for (int c1 = 128 * ii1; c1 <= min(N - 1, 128 * ii1 + 127); c1 += 1)
        for (int c2 = 128 * ii2; c2 <= min(M - 1, 128 * ii2 + 127); c2 += 1) {
          S9(c1, c2);
          S10(c1, c2);
        }
    } else if (ii0 == 1 && ii1 >= 0 && ii2 == 2 && ii3 == 0 && ii4 == 0 && ii5 == 0) {
      for (int c1 = 128 * ii1; c1 <= min(M - 1, 128 * ii1 + 127); c1 += 1) {
        S6(c1);
        S7(c1);
        S8(c1);
      }
    } else if (ii0 == 0 && ii1 >= 0 && ii2 == 2 && ii3 == 0 && ii4 == 0 && ii5 == 0) {
      for (int c1 = 128 * ii1; c1 <= min(M - 1, 128 * ii1 + 127); c1 += 1)
        S3(c1);
    } else if (ii0 == 3 && ii1 >= 0 && ii2 == 0 && ii3 == 0 && ii4 == 0 && ii5 == 0) {
      for (int c1 = 128 * ii1; c1 <= min(M - 2, 128 * ii1 + 127); c1 += 1)
        S11(c1);
    } else if (ii0 == 1 && ii1 >= 0 && ii2 == 0 && ii3 == 0 && ii4 == 0 && ii5 == 0) {
      for (int c1 = 128 * ii1; c1 <= min(M - 1, 128 * ii1 + 127); c1 += 1)
        S4(c1);
    } else if (ii0 == 0 && ii1 >= 0 && ii2 == 0 && ii3 == 0 && ii4 == 0 && ii5 == 0) {
      for (int c1 = 128 * ii1; c1 <= min(M - 1, 128 * ii1 + 127); c1 += 1)
        S1(c1);
    } else if (ii0 == 4 && ii1 == 0 && ii2 == 0 && ii3 == 0 && ii4 == 0 && ii5 == 0)
      S15();
  }
  return 0;
}

/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_correlation()
{
#pragma scop
  isl_ctx* ctx = isl_ctx_alloc();
  isl_map* rtile = isl_map_read_from_str(ctx, "[M, N] -> { [3, i1, 1, i3, 1, i5] -> [3, i1, 1, i3, 1, o5] : i1 >= 0 and i3 >= i1 and 128i3 <= -2 + M and i5 >= 0 and o5 >= 1 + i5 and 128o5 <= -1 + N; [2, i1, i2, 0, 0, 0] -> [3, o1, 1, o3, 1, i1] : i1 >= 0 and 128i1 <= -1 + N and 128i2 <= -1 + M and o1 >= 0 and o3 >= -1 + i2 and o3 >= o1 and o3 <= i2 and 128o3 <= -2 + M; [3, i1, 1, i3, 0, 0] -> [3, i1, 1, i3, 1, o5] : i1 >= 0 and i3 >= i1 and 128i3 <= -2 + M and o5 >= 0 and 128o5 <= -1 + N; [2, i1, i2, 0, 0, 0] -> [3, i2, 1, o3, 1, i1] : i1 >= 0 and 128i1 <= -1 + N and i2 >= 0 and o3 >= i2 and 128o3 <= -2 + M; [3, i1, 1, i3, 1, i5] -> [3, i1, 1, i3, 2, 0] : i1 >= 0 and i3 >= i1 and 128i3 <= -2 + M and i5 >= 0 and 128i5 <= -1 + N; [1, i1, 1, i3, 0, 0] -> [1, i1, 1, o3, 0, 0] : i1 >= 0 and 128i1 <= -1 + M and i3 >= 0 and o3 >= 1 + i3 and 128o3 <= -1 + N; [0, i1, 1, i3, 0, 0] -> [1, i1, 1, o3, 0, 0] : i1 >= 0 and 128i1 <= -1 + M and i3 >= 0 and 128i3 <= -1 + N and o3 >= 0 and 128o3 <= -1 + N; [0, i1, 1, i3, 0, 0] -> [0, i1, 1, o3, 0, 0] : i1 >= 0 and 128i1 <= -1 + M and i3 >= 0 and o3 >= 1 + i3 and 128o3 <= -1 + N; [1, i1, 1, i3, 0, 0] -> [2, o1, i1, 0, 0, 0] : i1 >= 0 and 128i1 <= -1 + M and i3 >= 0 and 128i3 <= -1 + N and o1 >= 0 and 128o1 <= -1 + N; [0, i1, 1, i3, 0, 0] -> [2, o1, i1, 0, 0, 0] : i1 >= 0 and 128i1 <= -1 + M and i3 >= 0 and 128i3 <= -1 + N and o1 >= 0 and 128o1 <= -1 + N; [3, i1, 1, i3, 0, 0] -> [3, i1, 1, i3, 2, 0] : i1 >= 0 and i3 >= i1 and 128i3 <= -2 + M; [0, i1, 2, 0, 0, 0] -> [1, i1, 1, o3, 0, 0] : i1 >= 0 and 128i1 <= -1 + M and o3 >= 0 and 128o3 <= -1 + N; [1, i1, 0, 0, 0, 0] -> [1, i1, 1, o3, 0, 0] : i1 >= 0 and 128i1 <= -1 + M and o3 >= 0 and 128o3 <= -1 + N; [0, i1, 0, 0, 0, 0] -> [1, i1, 1, o3, 0, 0] : i1 >= 0 and 128i1 <= -1 + M and o3 >= 0 and 128o3 <= -1 + N; [0, i1, 0, 0, 0, 0] -> [0, i1, 1, o3, 0, 0] : i1 >= 0 and 128i1 <= -1 + M and o3 >= 0 and 128o3 <= -1 + N; [1, i1, 1, i3, 0, 0] -> [1, i1, 2, 0, 0, 0] : i1 >= 0 and 128i1 <= -1 + M and i3 >= 0 and 128i3 <= -1 + N; [0, i1, 1, i3, 0, 0] -> [0, i1, 2, 0, 0, 0] : i1 >= 0 and 128i1 <= -1 + M and i3 >= 0 and 128i3 <= -1 + N; [1, i1, 2, 0, 0, 0] -> [2, o1, i1, 0, 0, 0] : i1 >= 0 and 128i1 <= -1 + M and o1 >= 0 and 128o1 <= -1 + N; [0, i1, 2, 0, 0, 0] -> [2, o1, i1, 0, 0, 0] : i1 >= 0 and 128i1 <= -1 + M and o1 >= 0 and 128o1 <= -1 + N; [1, i1, 0, 0, 0, 0] -> [2, o1, i1, 0, 0, 0] : i1 >= 0 and 128i1 <= -1 + M and o1 >= 0 and 128o1 <= -1 + N; [0, i1, 0, 0, 0, 0] -> [2, o1, i1, 0, 0, 0] : i1 >= 0 and 128i1 <= -1 + M and o1 >= 0 and 128o1 <= -1 + N; [1, i1, 1, i3, 0, 0] -> [2, i3, i1, 0, 0, 0] : i1 >= 0 and 128i1 <= -1 + M and i3 >= 0 and 128i3 <= -1 + N; [0, i1, 1, i3, 0, 0] -> [2, i3, i1, 0, 0, 0] : i1 >= 0 and 128i1 <= -1 + M and i3 >= 0 and 128i3 <= -1 + N; [1, i1, 0, 0, 0, 0] -> [1, i1, 2, 0, 0, 0] : i1 >= 0 and 128i1 <= -1 + M; [0, i1, 0, 0, 0, 0] -> [0, i1, 2, 0, 0, 0] : i1 >= 0 and 128i1 <= -1 + M }");
  isl_set* ii_set = isl_set_read_from_str(ctx, "[M, N] -> { [3, i1, 1, i3, 1, i5] : i1 >= 0 and 128i1 <= -2 + M and i3 >= 0 and 128i3 <= -2 + M and i5 >= 0 and 128i5 <= -1 + N; [3, i1, 1, i3, 2, 0] : i1 >= 0 and 128i1 <= -2 + M and i3 >= 0 and 128i3 <= -2 + M; [3, i1, 1, i3, 0, 0] : i1 >= 0 and 128i1 <= -2 + M and i3 >= 0 and 128i3 <= -2 + M; [1, i1, 1, i3, 0, 0] : i1 >= 0 and 128i1 <= -1 + M and i3 >= 0 and 128i3 <= -1 + N; [0, i1, 1, i3, 0, 0] : i1 >= 0 and 128i1 <= -1 + M and i3 >= 0 and 128i3 <= -1 + N; [2, i1, i2, 0, 0, 0] : i1 >= 0 and 128i1 <= -1 + N and i2 >= 0 and 128i2 <= -1 + M; [1, i1, 2, 0, 0, 0] : i1 >= 0 and 128i1 <= -1 + M; [0, i1, 2, 0, 0, 0] : i1 >= 0 and 128i1 <= -1 + M; [3, i1, 0, 0, 0, 0] : i1 >= 0 and 128i1 <= -2 + M; [1, i1, 0, 0, 0, 0] : i1 >= 0 and 128i1 <= -1 + M; [0, i1, 0, 0, 0, 0] : i1 >= 0 and 128i1 <= -1 + M; [4, 0, 0, 0, 0, 0] }");
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
  kernel_correlation ();

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(m, POLYBENCH_ARRAY(corr)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(data);
  POLYBENCH_FREE_ARRAY(corr);
  POLYBENCH_FREE_ARRAY(mean);
  POLYBENCH_FREE_ARRAY(stddev);

  return 0;
}
