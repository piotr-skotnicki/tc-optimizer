/**
 * This version is stamped on Apr. 14, 2015
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

/* ./tc ../examples/polybench/cholesky.scop.c --merge-tiling --free-scheduling --omp-for-codegen -b 64 --debug */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#pragma scop
{
  if (_PB_N >= 65)
    for (int k = 0; k <= (3 * _PB_N - 6) / 64; k += 1) {
      #pragma omp parallel for
      for (int ii0 = max(k - (_PB_N + 29) / 32 + 1, (k + 1) / 3); ii0 <= min(k, (_PB_N - 1) / 64); ii0 += 1) {
        if (k >= 1) {
          for (int ii2 = -ii0 + (k + ii0 + 1) / 2; ii2 <= min(min(ii0, k - ii0), (_PB_N - 2) / 64); ii2 += 1) {
            if (3 * ii0 == k && 3 * ii2 == k) {
              A[64 * k / 3][64 * k / 3] = SQRT_FUN(A[64 * k / 3][64 * k / 3]);
              if (64 * k + 195 >= 3 * _PB_N) {
                A[(64 * k / 3) + 1][64 * k / 3] /= A[64 * k / 3][64 * k / 3];
                A[(64 * k / 3) + 1][(64 * k / 3) + 1] -= (A[(64 * k / 3) + 1][64 * k / 3] * A[(64 * k / 3) + 1][64 * k / 3]);
                A[(64 * k / 3) + 1][(64 * k / 3) + 1] = SQRT_FUN(A[(64 * k / 3) + 1][(64 * k / 3) + 1]);
              } else {
                A[(64 * k / 3) + 1][64 * k / 3] /= A[64 * k / 3][64 * k / 3];
                A[(64 * k / 3) + 1][(64 * k / 3) + 1] -= (A[(64 * k / 3) + 1][64 * k / 3] * A[(64 * k / 3) + 1][64 * k / 3]);
                A[(64 * k / 3) + 1][(64 * k / 3) + 1] = SQRT_FUN(A[(64 * k / 3) + 1][(64 * k / 3) + 1]);
              }
            }
            for (int c6 = max(max(64 * ii0, 32 * k - 32 * ii0 + 2), 64 * ii2 + 1); c6 <= min(_PB_N - 1, 64 * ii0 + 63); c6 += 1) {
              for (int c8 = 64 * ii2; c8 <= min(64 * ii2 + 63, c6 - 1); c8 += 1) {
                for (int c10 = 64 * k - 64 * ii0 - 64 * ii2; c10 <= min(64 * k - 64 * ii0 - 64 * ii2 + 63, c8 - 1); c10 += 1)
                  A[c6][c8] -= (A[c6][c10] * A[c8][c10]);
                if (ii0 + 2 * ii2 == k)
                  A[c6][c8] /= A[c8][c8];
              }
              if (_PB_N >= 64 * ii0 + 2 && ii0 + 2 * ii2 == k) {
                for (int c8 = 32 * k - 32 * ii0; c8 <= min(32 * k - 32 * ii0 + 63, c6 - 1); c8 += 1)
                  A[c6][c6] -= (A[c6][c8] * A[c6][c8]);
                if (3 * ii0 == k)
                  A[c6][c6] = SQRT_FUN(A[c6][c6]);
              } else if (3 * _PB_N >= 64 * k + 259 && 64 * ii0 + 1 == _PB_N && _PB_N + 128 * ii2 == 64 * k + 1 && c6 + 1 == _PB_N)
                for (int c8 = ((-_PB_N + 1) / 2) + 32 * k; c8 <= ((-_PB_N + 127) / 2) + 32 * k; c8 += 1)
                  A[_PB_N - 1][_PB_N - 1] -= (A[_PB_N - 1][c8] * A[_PB_N - 1][c8]);
            }
          }
          if (64 * k + 67 == 3 * _PB_N && 64 * ii0 + 1 == _PB_N)
            for (int c8 = _PB_N - 65; c8 < _PB_N - 1; c8 += 1)
              A[_PB_N - 1][_PB_N - 1] -= (A[_PB_N - 1][c8] * A[_PB_N - 1][c8]);
        } else
          for (int c6 = 0; c6 <= 63; c6 += 1) {
            for (int c8 = 0; c8 < c6; c8 += 1) {
              for (int c10 = 0; c10 < c8; c10 += 1)
                A[c6][c8] -= (A[c6][c10] * A[c8][c10]);
              A[c6][c8] /= A[c8][c8];
            }
            for (int c8 = 0; c8 < c6; c8 += 1)
              A[c6][c6] -= (A[c6][c8] * A[c6][c8]);
            A[c6][c6] = SQRT_FUN(A[c6][c6]);
          }
      }
      if (64 * k + 6 == 3 * _PB_N)
        for (int c6 = _PB_N - 2; c6 < _PB_N; c6 += 1) {
          if (c6 + 1 == _PB_N) {
            A[_PB_N - 1][_PB_N - 2] /= A[_PB_N - 2][_PB_N - 2];
            A[_PB_N - 1][_PB_N - 1] -= (A[_PB_N - 1][_PB_N - 2] * A[_PB_N - 1][_PB_N - 2]);
          }
          A[c6][c6] = SQRT_FUN(A[c6][c6]);
        }
    }
  if ((_PB_N - 1) % 64 == 0)
    A[_PB_N - 1][_PB_N - 1] = SQRT_FUN(A[_PB_N - 1][_PB_N - 1]);
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
