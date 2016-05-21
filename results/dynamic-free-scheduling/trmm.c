/**
 * This version is stamped on Apr. 14, 2015
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* trmm.c: this file is part of PolyBench/C */

#include <omp.h>
#include <math.h>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

#define S2_I(i,j) B[i][j] = alpha * B[i][j];
#define S2(i,j) S2_I((i),(j))
#define S1_I(i,j,k) B[i][j] += A[k][i] * B[k][j];
#define S1(i,j,k) S1_I((i),(j),(k))

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
#include "trmm.h"

/* Variable declaration/allocation. */
DATA_TYPE alpha;
POLYBENCH_2D_ARRAY_DECL(A,DATA_TYPE,M,M,m,m);
POLYBENCH_2D_ARRAY_DECL(B,DATA_TYPE,M,N,m,n);

/* Array initialization. */
static
void init_array(int m, int n,
		DATA_TYPE *alpha,
		DATA_TYPE POLYBENCH_2D(A,M,M,m,m),
		DATA_TYPE POLYBENCH_2D(B,M,N,m,n))
{
  int i, j;

  *alpha = 1.5;
  for (i = 0; i < m; i++) {
    for (j = 0; j < i; j++) {
      A[i][j] = (DATA_TYPE)((i+j) % m)/m;
    }
    A[i][i] = 1.0;
    for (j = 0; j < n; j++) {
      B[i][j] = (DATA_TYPE)((n+(i-j)) % n)/n;
    }
 }

}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int m, int n,
		 DATA_TYPE POLYBENCH_2D(B,M,N,m,n))
{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("B");
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++) {
	if ((i * m + j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
	fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, B[i][j]);
    }
  POLYBENCH_DUMP_END("B");
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
    if (ii0 >= 0 && ii1 >= 0 && ii2 == 0) {
      for (int c0 = 256 * ii0; c0 <= min(min(M - 2, 256 * ii0 + 255), 256 * ii3 + 255); c0 += 1)
        for (int c1 = 256 * ii1; c1 <= min(N - 1, 256 * ii1 + 255); c1 += 1)
          for (int c3 = max(256 * ii3 + 1, c0 + 1); c3 <= min(M - 1, 256 * ii3 + 256); c3 += 1)
            S1(c0, c1, c3);
    } else if (ii0 >= 0 && ii1 >= 0 && ii2 == 1 && ii3 == 0)
      for (int c0 = 256 * ii0; c0 <= min(M - 1, 256 * ii0 + 255); c0 += 1)
        for (int c1 = 256 * ii1; c1 <= min(N - 1, 256 * ii1 + 255); c1 += 1)
          S2(c0, c1);
  }
  return 0;
}

/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_trmm()
{
#pragma scop
  isl_ctx* ctx = isl_ctx_alloc();
  isl_map* rtile = isl_map_read_from_str(ctx, "[M, N] -> { [i0, i1, 0, i3] -> [o0, i1, 0, o3] : i0 >= 0 and i1 >= 0 and 256i1 <= -1 + N and i3 >= i0 and o0 >= i3 and o0 <= 1 + i3 and o3 >= 1 + i3 and 256o3 <= -2 + M; [i0, i1, 0, i3] -> [o0, i1, 1, 0] : i0 >= 0 and i1 >= 0 and 256i1 <= -1 + N and i3 >= i0 and 256i3 <= -2 + M and o0 >= i3 and o0 <= 1 + i3 and 256o0 <= -1 + M; [i0, i1, 0, i3] -> [i0, i1, 0, o3] : i0 >= 0 and i1 >= 0 and 256i1 <= -1 + N and i3 >= i0 and o3 >= 1 + i3 and 256o3 <= -2 + M; [i0, i1, 0, i3] -> [i3, i1, 0, i3] : i0 >= 0 and i1 >= 0 and 256i1 <= -1 + N and i3 >= 1 + i0 and 256i3 <= -3 + M; [i0, i1, 0, i3] -> [i0, i1, 1, 0] : i0 >= 0 and i1 >= 0 and 256i1 <= -1 + N and i3 >= i0 and 256i3 <= -2 + M }");
  isl_set* ii_set = isl_set_read_from_str(ctx, "[N, M] -> { [i0, i1, 0, i3] : i0 >= 0 and 256i0 <= -2 + M and i1 >= 0 and 256i1 <= -1 + N and i3 >= 0 and 256i3 <= -2 + M; [i0, i1, 1, 0] : i0 >= 0 and 256i0 <= -1 + M and i1 >= 0 and 256i1 <= -1 + N }");
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
  init_array (m, n, &alpha, POLYBENCH_ARRAY(A), POLYBENCH_ARRAY(B));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_trmm ();

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(m, n, POLYBENCH_ARRAY(B)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(B);

  return 0;
}
