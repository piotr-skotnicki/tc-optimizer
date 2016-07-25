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

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "trmm.h"


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


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_trmm(int m, int n,
		 DATA_TYPE alpha,
		 DATA_TYPE POLYBENCH_2D(A,M,M,m,m),
		 DATA_TYPE POLYBENCH_2D(B,M,N,m,n))
{

/* ./tc ../examples/polybench/trmm.scop.c --correction-tiling --sfs-single-scheduling --omp-for-codegen --debug -b 8 */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#pragma scop
if (_PB_M >= 1) {
  #pragma omp parallel for
  for (register int ir1 = 0; ir1 < _PB_N; ir1 += 1) {
    if (_PB_M >= 2) {
      const int ii0_prim_lb = 0, ii0_prim_ub = (_PB_M - 1) / 8;
      for (register int ii0_prim = ii0_prim_lb; ii0_prim <= ii0_prim_ub; ii0_prim += 1) {
        const int ii3_prim_lb = ii0_prim, ii3_prim_ub = (_PB_M - 1) / 8;
        for (register int ii3_prim = ii3_prim_lb; ii3_prim <= ii3_prim_ub; ii3_prim += 1) {
          const int i0_lb = 8 * ii0_prim, i0_ub = min(min(_PB_M - 2, 8 * ii0_prim + 7), 8 * ii3_prim + 6);
          for (register int i0 = i0_lb; i0 <= i0_ub; i0 += 1) {
            const int i3_lb = max(8 * ii3_prim, i0 + 1), i3_ub = min(_PB_M - 1, 8 * ii3_prim + 7);
            for (register int i3 = i3_lb; i3 <= i3_ub; i3 += 1) {
              B[i0][ir1] += (A[i3][i0] * B[i3][ir1]);
            }
          }
        }
        if (ii0_prim >= 1) {
          const int i0_lb = 8 * ii0_prim, i0_ub = min(_PB_M - 1, 8 * ii0_prim + 7);
          for (register int i0 = i0_lb; i0 <= i0_ub; i0 += 1) {
            B[i0][ir1] = (alpha * B[i0][ir1]);
          }
        } else {
          const int i0_lb = 0, i0_ub = min(7, _PB_M - 1);
          for (register int i0 = i0_lb; i0 <= i0_ub; i0 += 1) {
            B[i0][ir1] = (alpha * B[i0][ir1]);
          }
        }
      }
    } else {
      B[0][ir1] = (alpha * B[0][ir1]);
    }
  }
}
#pragma endscop

}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int m = M;
  int n = N;

  /* Variable declaration/allocation. */
  DATA_TYPE alpha;
  POLYBENCH_2D_ARRAY_DECL(A,DATA_TYPE,M,M,m,m);
  POLYBENCH_2D_ARRAY_DECL(B,DATA_TYPE,M,N,m,n);

  /* Initialize array(s). */
  init_array (m, n, &alpha, POLYBENCH_ARRAY(A), POLYBENCH_ARRAY(B));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_trmm (m, n, alpha, POLYBENCH_ARRAY(A), POLYBENCH_ARRAY(B));

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
