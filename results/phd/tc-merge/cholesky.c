/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* cholesky.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "cholesky.h"


/* Array initialization. */
static
void init_array(int n,
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
    for (j = 0; j <= i; j++) {
    if ((i * n + j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
    fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, A[i][j]);
  }
  POLYBENCH_DUMP_END("A");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_cholesky(int n,
		     DATA_TYPE POLYBENCH_2D(A,N,N,n,n))
{
  int i, j, k;


/* TC Optimizing Compiler 0.5.1 */
/* ./tc ../examples/polybench/cholesky.scop.c --merge-tiling --isl-wave-scheduling --omp-for-codegen --iterative-tc --debug -b 8 --align */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#pragma scop
{
  for (register int k0 = 0; k0 < (3 * _PB_N + 1) / 8; k0 += 1) {
    #pragma omp parallel for
    for (register int k1 = max(k0 - (_PB_N + 1) / 4 + 1, (k0 + 1) / 3); k1 <= min(k0, (_PB_N - 1) / 8); k1 += 1) {
      if (3 * k1 >= k0 + 1) {
        for (register int k2 = -k1 + (k0 + k1 + 1) / 2; k2 <= min(k1 - 1, k0 - k1); k2 += 1) {
          for (register int i0 = 8 * k1; i0 <= min(_PB_N - 1, 8 * k1 + 7); i0 += 1) {
            for (register int i2 = 8 * k2; i2 <= 8 * k2 + 7; i2 += 1) {
              for (register int i4 = 8 * k0 - 8 * k1 - 8 * k2; i4 <= min(8 * k0 - 8 * k1 - 8 * k2 + 7, i2 - 1); i4 += 1) {
                A[i0][i2] -= (A[i0][i4] * A[i2][i4]);
              }
              if (k1 + 2 * k2 == k0) {
                A[i0][i2] /= A[i2][i2];
              }
            }
          }
        }
        if (k0 >= 2 * k1) {
          for (register int i0 = 8 * k1; i0 <= min(_PB_N - 1, 8 * k1 + 7); i0 += 1) {
            for (register int i2 = 8 * k1; i2 < i0; i2 += 1) {
              for (register int i4 = 8 * k0 - 16 * k1; i4 <= 8 * k0 - 16 * k1 + 7; i4 += 1) {
                A[i0][i2] -= (A[i0][i4] * A[i2][i4]);
              }
            }
            for (register int i2 = 8 * k0 - 16 * k1; i2 <= 8 * k0 - 16 * k1 + 7; i2 += 1) {
              A[i0][i0] -= (A[i0][i2] * A[i0][i2]);
            }
          }
        }
      } else if (3 * k1 == k0) {
        for (register int i0 = 8 * k0 / 3; i0 <= min(_PB_N - 1, (8 * k0 / 3) + 7); i0 += 1) {
          {
            for (register int i2 = 8 * k0 / 3; i2 < i0; i2 += 1) {
              for (register int i4 = 8 * k0 / 3; i4 < i2; i4 += 1) {
                A[i0][i2] -= (A[i0][i4] * A[i2][i4]);
              }
              A[i0][i2] /= A[i2][i2];
            }
            for (register int i2 = 8 * k0 / 3; i2 < i0; i2 += 1) {
              A[i0][i0] -= (A[i0][i2] * A[i0][i2]);
            }
          }
          A[i0][i0] = SQRT_FUN(A[i0][i0]);
        }
      }
    }
  }
  if (_PB_N % 8 <= 2) {
    for (register int i0 = -(_PB_N % 8) + _PB_N; i0 < _PB_N; i0 += 1) {
      if (i0 + 1 == _PB_N && (_PB_N - 2) % 8 == 0) {
        A[_PB_N - 1][_PB_N - 2] /= A[_PB_N - 2][_PB_N - 2];
        A[_PB_N - 1][_PB_N - 1] -= (A[_PB_N - 1][_PB_N - 2] * A[_PB_N - 1][_PB_N - 2]);
      }
      A[i0][i0] = SQRT_FUN(A[i0][i0]);
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
  kernel_cholesky (n, POLYBENCH_ARRAY(A));

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
