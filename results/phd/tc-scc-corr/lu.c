/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* lu.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "lu.h"


/* Array initialization. */
static
void init_array (int n,
		 DATA_TYPE POLYBENCH_2D(A,N,N,n,n))
{
  int i, j;

  for (i = 0; i < n; i++)
    {
      for (j = 0; j <= i; j++)
	A[i][j] = (DATA_TYPE)(-j % n) / n + 1;
      for (j = i+1; j < n; j++) {
	A[i][j] = 0;
      }
      A[i][i] = 1;
    }

  /* Make the matrix positive semi-definite. */
  /* not necessary for LU, but using same code as cholesky */
  int r,s,t;
  POLYBENCH_2D_ARRAY_DECL(B, DATA_TYPE, N, N, n, n);
  for (r = 0; r < n; ++r)
    for (s = 0; s < n; ++s)
      (POLYBENCH_ARRAY(B))[r][s] = 0;
  for (t = 0; t < n; ++t)
    for (r = 0; r < n; ++r)
      for (s = 0; s < n; ++s)
	(POLYBENCH_ARRAY(B))[r][s] += A[r][t] * A[s][t];
    for (r = 0; r < n; ++r)
      for (s = 0; s < n; ++s)
	A[r][s] = (POLYBENCH_ARRAY(B))[r][s];
  POLYBENCH_FREE_ARRAY(B);

}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_2D(A,N,N,n,n))

{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("A");
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
      if ((i * n + j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
      fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, A[i][j]);
    }
  POLYBENCH_DUMP_END("A");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_lu(int n,
	       DATA_TYPE POLYBENCH_2D(A,N,N,n,n))
{
  int i, j, k;

/* TC Optimizing Compiler 0.5.1 */
/* ./tc ../examples/polybench/lu.scop.c --scc-correction-tiling --isl-wave-scheduling --omp-for-codegen --iterative-tc --debug -b 16 --align */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#pragma scop
for (register int k0 = 0; k0 <= floord(3 * _PB_N - 4, 16); k0 += 1) {
  #pragma omp parallel for
  for (register int k1 = max(0, floord(-_PB_N + 16 * k0, 32) + 1); k1 <= min(k0, (_PB_N - 1) / 16); k1 += 1) {
    if (3 * k1 >= k0 + 2 && (k0 + k1) % 2 == 0) {
      for (register int i0 = 16 * k1; i0 <= min(_PB_N - 1, 16 * k1 + 15); i0 += 1) {
        for (register int i2 = 8 * k0 - 8 * k1; i2 <= 8 * k0 - 8 * k1 + 15; i2 += 1) {
          for (register int i4 = 8 * k0 - 8 * k1; i4 < i2; i4 += 1) {
            A[i0][i2] -= (A[i0][i4] * A[i4][i2]);
          }
          A[i0][i2] /= A[i2][i2];
        }
      }
    }
    for (register int k2 = -k1 + (k0 + k1) / 2 + 1; k2 <= min(k1 - 1, k0 - k1); k2 += 1) {
      for (register int i0 = 16 * k1; i0 <= min(_PB_N - 1, 16 * k1 + 15); i0 += 1) {
        for (register int i2 = 16 * k2; i2 <= 16 * k2 + 15; i2 += 1) {
          for (register int i4 = 16 * k0 - 16 * k1 - 16 * k2; i4 <= 16 * k0 - 16 * k1 - 16 * k2 + 15; i4 += 1) {
            A[i0][i2] -= (A[i0][i4] * A[i4][i2]);
          }
        }
      }
    }
    for (register int k2 = max(k1, k0 - 2 * k1); k2 <= min(k0 - k1, (_PB_N - 1) / 16); k2 += 1) {
      {
        if (3 * k1 == k0 && 3 * k2 == k0) {
          for (register int i0 = (16 * k0 / 3) + 1; i0 <= min(_PB_N - 1, (16 * k0 / 3) + 15); i0 += 1) {
            A[i0][16 * k0 / 3] /= A[16 * k0 / 3][16 * k0 / 3];
            for (register int i2 = (16 * k0 / 3) + 1; i2 < i0; i2 += 1) {
              A[i0][i2] -= (A[i0][16 * k0 / 3] * A[16 * k0 / 3][i2]);
            }
          }
        }
        for (register int i0 = max(16 * k1, 16 * k0 - 16 * k1 - 16 * k2 + 1); i0 <= min(_PB_N - 1, 16 * k1 + 15); i0 += 1) {
          if (3 * k1 == k0 && 3 * k2 == k0) {
            for (register int i2 = (16 * k0 / 3) + 1; i2 < i0; i2 += 1) {
              for (register int i4 = (16 * k0 / 3) + 1; i4 < i2; i4 += 1) {
                A[i0][i2] -= (A[i0][i4] * A[i4][i2]);
              }
              A[i0][i2] /= A[i2][i2];
            }
          }
          for (register int i2 = max(16 * k2, i0); i2 <= min(_PB_N - 1, 16 * k2 + 15); i2 += 1) {
            for (register int i3 = 16 * k0 - 16 * k1 - 16 * k2; i3 <= min(16 * k0 - 16 * k1 - 16 * k2 + 15, i0 - 1); i3 += 1) {
              A[i0][i2] -= (A[i0][i3] * A[i3][i2]);
            }
          }
        }
      }
      if (3 * k1 >= k0 + 1 && k2 == k1) {
        for (register int i0 = 16 * k1 + 1; i0 <= min(_PB_N - 1, 16 * k1 + 15); i0 += 1) {
          for (register int i2 = 16 * k1; i2 < i0; i2 += 1) {
            for (register int i4 = 16 * k0 - 32 * k1; i4 <= 16 * k0 - 32 * k1 + 15; i4 += 1) {
              A[i0][i2] -= (A[i0][i4] * A[i4][i2]);
            }
          }
        }
      }
    }
  }
}
#pragma endscop

}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int n = N;

  /* Variable declaration/allocation. */
  POLYBENCH_2D_ARRAY_DECL(A, DATA_TYPE, N, N, n, n);

  /* Initialize array(s). */
  init_array (n, POLYBENCH_ARRAY(A));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_lu (n, POLYBENCH_ARRAY(A));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(A)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(A);

  return 0;
}
