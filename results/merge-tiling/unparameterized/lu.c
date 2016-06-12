/**
 * This version is stamped on Apr. 14, 2015
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

/* ./tc ../examples/polybench/lu.scop.c --merge-tiling --free-scheduling --omp-for-codegen --debug -b 64 -D _PB_N=4000 --isl-union-map-tc */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#pragma scop
for (int k = 0; k <= 186; k += 1) {
  if (k >= 185) {
    for (int c6 = 3969; c6 <= 3999; c6 += 1) {
      for (int c8 = 3968; c8 < c6; c8 += 1) {
        for (int c10 = 64 * k - 7936; c10 < min(64 * k - 7872, c8); c10 += 1)
          A[c6][c8] -= (A[c6][c10] * A[c10][c8]);
        if (k == 186)
          A[c6][c8] /= A[c8][c8];
      }
      if (k == 186)
        for (int c8 = c6; c8 <= 3999; c8 += 1)
          for (int c9 = 3968; c9 < c6; c9 += 1)
            A[c6][c8] -= (A[c6][c9] * A[c9][c8]);
    }
    if (k == 185)
      for (int c6 = 3969; c6 <= 3999; c6 += 1)
        for (int c8 = c6; c8 <= 3999; c8 += 1)
          for (int c9 = 3904; c9 <= 3967; c9 += 1)
            A[c6][c8] -= (A[c6][c9] * A[c9][c8]);
  } else
    #pragma omp parallel for
    for (int ii0 = max(0, (k + 1) / 2 - 31); ii0 <= min(62, k); ii0 += 1) {
      if (k >= 1) {
        if (3 * ii0 >= k) {
          for (int ii2 = -ii0 + (k + ii0 + 1) / 2; ii2 <= min(min(61, ii0), k - ii0); ii2 += 1)
            for (int c6 = 64 * ii0 + 1; c6 <= min(3999, 64 * ii0 + 64); c6 += 1) {
              if (ii0 == 62 && 2 * ii2 + 62 == k) {
                A[c6][32 * k - 1984] /= A[32 * k - 1984][32 * k - 1984];
              } else if (ii0 + 2 * ii2 == k)
                A[c6][32 * k - 32 * ii0] /= A[32 * k - 32 * ii0][32 * k - 32 * ii0];
              for (int c8 = max(64 * ii2, 64 * k - 64 * ii0 - 64 * ii2 + 1); c8 <= min(64 * ii2 + 63, c6 - 1); c8 += 1) {
                for (int c10 = 64 * k - 64 * ii0 - 64 * ii2; c10 <= min(64 * k - 64 * ii0 - 64 * ii2 + 63, c8 - 1); c10 += 1)
                  A[c6][c8] -= (A[c6][c10] * A[c10][c8]);
                if (ii0 <= 61 && ii0 + 2 * ii2 == k) {
                  A[c6][c8] /= A[c8][c8];
                } else if (ii0 == 62 && 2 * ii2 + 62 == k)
                  A[c6][c8] /= A[c8][c8];
              }
              if (3 * ii0 == k && 3 * ii2 == k)
                for (int c8 = c6; c8 <= (64 * k / 3) + 63; c8 += 1)
                  for (int c9 = 64 * k / 3; c9 < c6; c9 += 1)
                    A[c6][c8] -= (A[c6][c9] * A[c9][c8]);
            }
          if (k >= 124 && ii0 == 62)
            for (int c6 = 3969; c6 <= 3999; c6 += 1)
              for (int c8 = 3968; c8 < c6; c8 += 1)
                for (int c10 = 64 * k - 7936; c10 < 64 * k - 7872; c10 += 1)
                  A[c6][c8] -= (A[c6][c10] * A[c10][c8]);
          if (k >= 124 && ii0 == 62)
            for (int c6 = 3969; c6 <= 3999; c6 += 1)
              for (int c8 = c6; c8 <= 3999; c8 += 1)
                for (int c9 = 64 * k - 7936; c9 < 64 * k - 7872; c9 += 1)
                  A[c6][c8] -= (A[c6][c9] * A[c9][c8]);
        }
        if (ii0 <= 61)
          for (int ii2 = max(max(ii0, k - 2 * ii0), -ii0 + (k + ii0) / 2 + 1); ii2 <= min(62, k - ii0); ii2 += 1) {
            if (2 * ii0 + ii2 >= k + 1) {
              for (int c6 = 64 * ii0 + 1; c6 <= min(64 * ii0 + 64, 64 * ii2 + 63); c6 += 1)
                for (int c8 = max(64 * ii2, c6); c8 <= min(3999, 64 * ii2 + 63); c8 += 1)
                  for (int c9 = 64 * k - 64 * ii0 - 64 * ii2; c9 <= 64 * k - 64 * ii0 - 64 * ii2 + 63; c9 += 1)
                    A[c6][c8] -= (A[c6][c9] * A[c9][c8]);
            } else
              for (int c6 = 64 * ii0 + 1; c6 <= 64 * ii0 + 64; c6 += 1)
                for (int c8 = 64 * k - 128 * ii0; c8 <= min(3999, 64 * k - 128 * ii0 + 63); c8 += 1)
                  for (int c9 = 64 * ii0; c9 < c6; c9 += 1)
                    A[c6][c8] -= (A[c6][c9] * A[c9][c8]);
          }
      } else
        for (int c6 = 1; c6 <= 64; c6 += 1) {
          for (int c8 = 0; c8 < c6; c8 += 1) {
            for (int c10 = 0; c10 < c8; c10 += 1)
              A[c6][c8] -= (A[c6][c10] * A[c10][c8]);
            A[c6][c8] /= A[c8][c8];
          }
          for (int c8 = c6; c8 <= 63; c8 += 1)
            for (int c9 = 0; c9 < c6; c9 += 1)
              A[c6][c8] -= (A[c6][c9] * A[c9][c8]);
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
