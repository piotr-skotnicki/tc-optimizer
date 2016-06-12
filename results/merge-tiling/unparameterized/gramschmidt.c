/**
 * This version is stamped on Apr. 14, 2015
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* gramschmidt.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "gramschmidt.h"


/* Array initialization. */
static
void init_array(int m, int n,
		DATA_TYPE POLYBENCH_2D(A,M,N,m,n),
		DATA_TYPE POLYBENCH_2D(R,N,N,n,n),
		DATA_TYPE POLYBENCH_2D(Q,M,N,m,n))
{
  int i, j;

  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++) {
      A[i][j] = (((DATA_TYPE) ((i*j) % m) / m )*100) + 10;
      Q[i][j] = 0.0;
    }
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      R[i][j] = 0.0;
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int m, int n,
		 DATA_TYPE POLYBENCH_2D(A,M,N,m,n),
		 DATA_TYPE POLYBENCH_2D(R,N,N,n,n),
		 DATA_TYPE POLYBENCH_2D(Q,M,N,m,n))
{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("R");
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
	if ((i*n+j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
	fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, R[i][j]);
    }
  POLYBENCH_DUMP_END("R");

  POLYBENCH_DUMP_BEGIN("Q");
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++) {
	if ((i*n+j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
	fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, Q[i][j]);
    }
  POLYBENCH_DUMP_END("Q");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
/* QR Decomposition with Modified Gram Schmidt:
 http://www.inf.ethz.ch/personal/gander/ */
static
void kernel_gramschmidt(int m, int n,
			DATA_TYPE POLYBENCH_2D(A,M,N,m,n),
			DATA_TYPE POLYBENCH_2D(R,N,N,n,n),
			DATA_TYPE POLYBENCH_2D(Q,M,N,m,n))
{
  int i, j, k;

  DATA_TYPE nrm;

/* ./tc ../examples/polybench/gramschmidt.scop.c --merge-tiling --free-scheduling --omp-for-codegen --debug -b 64 -D _PB_M=2000 -D _PB_N=2600 */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#pragma scop
{
  #pragma omp parallel for
  for (int ii0 = 0; ii0 <= 40; ii0 += 1)
    for (int ii2 = ii0; ii2 <= 40; ii2 += 1)
      for (int c6 = 64 * ii0; c6 <= min(min(2598, 64 * ii0 + 63), 64 * ii2 + 62); c6 += 1) {
        if (ii2 <= 39) {
          for (int c8 = max(64 * ii2, c6 + 1); c8 <= 64 * ii2 + 63; c8 += 1)
            R[c6][c8] = SCALAR_VAL(0.0);
        } else
          for (int c8 = max(2560, c6 + 1); c8 <= 2599; c8 += 1)
            R[c6][c8] = SCALAR_VAL(0.0);
      }
  for (int k = 1; k <= 81; k += 1) {
    if (k % 2 == 0) {
      #pragma omp parallel for
      for (int ii2 = k / 2; ii2 <= 40; ii2 += 1)
        for (int c6 = 32 * k - 64; c6 < 32 * k; c6 += 1) {
          if (ii2 <= 39) {
            for (int c8 = 64 * ii2; c8 <= 64 * ii2 + 63; c8 += 1) {
              for (int c10 = 0; c10 <= 1999; c10 += 1)
                R[c6][c8] += (Q[c10][c6] * A[c10][c8]);
              for (int c10 = 0; c10 <= 1999; c10 += 1)
                A[c10][c8] = (A[c10][c8] - (Q[c10][c6] * R[c6][c8]));
            }
          } else
            for (int c8 = 2560; c8 <= 2599; c8 += 1) {
              for (int c10 = 0; c10 <= 1999; c10 += 1)
                R[c6][c8] += (Q[c10][c6] * A[c10][c8]);
              for (int c10 = 0; c10 <= 1999; c10 += 1)
                A[c10][c8] = (A[c10][c8] - (Q[c10][c6] * R[c6][c8]));
            }
        }
    } else
      for (int c6 = 32 * k - 32; c6 <= min(2599, 32 * k + 31); c6 += 1) {
        nrm = SCALAR_VAL(0.0);
        for (int c8 = 0; c8 <= 1999; c8 += 1)
          nrm += (A[c8][c6] * A[c8][c6]);
        R[c6][c6] = SQRT_FUN(nrm);
        for (int c8 = 0; c8 <= 1999; c8 += 1)
          Q[c8][c6] = (A[c8][c6] / R[c6][c6]);
        if (k == 81) {
          for (int c8 = c6 + 1; c8 <= 2599; c8 += 1) {
            for (int c10 = 0; c10 <= 1999; c10 += 1)
              R[c6][c8] += (Q[c10][c6] * A[c10][c8]);
            for (int c10 = 0; c10 <= 1999; c10 += 1)
              A[c10][c8] = (A[c10][c8] - (Q[c10][c6] * R[c6][c8]));
          }
        } else
          for (int c8 = c6 + 1; c8 <= 32 * k + 31; c8 += 1) {
            for (int c10 = 0; c10 <= 1999; c10 += 1)
              R[c6][c8] += (Q[c10][c6] * A[c10][c8]);
            for (int c10 = 0; c10 <= 1999; c10 += 1)
              A[c10][c8] = (A[c10][c8] - (Q[c10][c6] * R[c6][c8]));
          }
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
  POLYBENCH_2D_ARRAY_DECL(A,DATA_TYPE,M,N,m,n);
  POLYBENCH_2D_ARRAY_DECL(R,DATA_TYPE,N,N,n,n);
  POLYBENCH_2D_ARRAY_DECL(Q,DATA_TYPE,M,N,m,n);

  /* Initialize array(s). */
  init_array (m, n,
	      POLYBENCH_ARRAY(A),
	      POLYBENCH_ARRAY(R),
	      POLYBENCH_ARRAY(Q));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_gramschmidt (m, n,
		      POLYBENCH_ARRAY(A),
		      POLYBENCH_ARRAY(R),
		      POLYBENCH_ARRAY(Q));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(m, n, POLYBENCH_ARRAY(A), POLYBENCH_ARRAY(R), POLYBENCH_ARRAY(Q)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(R);
  POLYBENCH_FREE_ARRAY(Q);

  return 0;
}
