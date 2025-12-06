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
/* ./tc ../examples/polybench/lu.scop.c --merge-tiling --isl-wave-scheduling --omp-for-codegen --iterative-tc --debug -b 16 --align */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#pragma scop
if (_PB_N >= 17) {
  for (register int i0 = 1; i0 <= 15; i0 += 1) {
    for (register int i2 = 0; i2 < i0; i2 += 1) {
      for (register int i4 = 0; i4 < i2; i4 += 1) {
        A[i0][i2] -= (A[i0][i4] * A[i4][i2]);
      }
      A[i0][i2] /= A[i2][i2];
    }
    for (register int i2 = i0; i2 <= 15; i2 += 1) {
      for (register int i3 = 0; i3 < i0; i3 += 1) {
        A[i0][i2] -= (A[i0][i3] * A[i3][i2]);
      }
    }
  }
  for (register int k0 = 7; k0 <= 3 * _PB_N + 1; k0 += 1) {
    if (k0 == 3 * _PB_N && (_PB_N - 2) % 16 == 0) {
      A[_PB_N - 1][_PB_N - 2] /= A[_PB_N - 2][_PB_N - 2];
      A[_PB_N - 1][_PB_N - 1] -= (A[_PB_N - 1][_PB_N - 2] * A[_PB_N - 2][_PB_N - 1]);
    }
    #pragma omp parallel for
    for (register int k1 = max(max(2, 16 * floord(-3 * _PB_N + 16 * k0 - 62, 720) + 18), 16 * floord(-_PB_N + 16 * k0 - 96, 752) + 18); k1 <= min(min(_PB_N + 1, (8 * k0 - 2) / 23), (16 * k0 - 22) / 45); k1 += 16) {
      for (register int k2 = max((15 * k1 + 2) / 16, floord(-_PB_N + 16 * k0 - 15 * k1 - 2, 32) + 1); k2 <= min(min(k1, (k0 - 1) / 3), -k1 + (k0 + k1) / 2); k2 += 1) {
        for (register int i0 = max(k1 - 2, -15 * k1 + 16 * k2 - 1); i0 <= min(_PB_N - 1, k1 + 13); i0 += 1) {
          for (register int i2 = max(16 * k0 - 15 * k1 - 32 * k2 - 2, i0); i2 <= min(_PB_N - 1, 16 * k0 - 15 * k1 - 32 * k2 + 13); i2 += 1) {
            for (register int i3 = -15 * k1 + 16 * k2 - 2; i3 <= min(-15 * k1 + 16 * k2 + 13, i0 - 1); i3 += 1) {
              A[i0][i2] -= (A[i0][i3] * A[i3][i2]);
            }
          }
        }
        if (_PB_N + 23 * k1 >= 8 * k0 + 17 && k1 + 2 * k2 == k0) {
          if (3 * k1 == k0 + 2) {
            A[(k0 - 1) / 3][(k0 - 4) / 3] /= A[(k0 - 4) / 3][(k0 - 4) / 3];
            for (register int i2 = (k0 - 1) / 3; i2 <= min(_PB_N - 1, (k0 + 41) / 3); i2 += 1) {
              A[(k0 - 1) / 3][i2] -= (A[(k0 - 1) / 3][(k0 - 4) / 3] * A[(k0 - 4) / 3][i2]);
            }
          }
          for (register int i0 = max(k1 - 2, 8 * k0 - 23 * k1 + 16); i0 <= min(_PB_N - 1, k1 + 13); i0 += 1) {
            for (register int i2 = 8 * k0 - 23 * k1 + 14; i2 <= min(8 * k0 - 23 * k1 + 29, i0 - 1); i2 += 1) {
              for (register int i4 = 8 * k0 - 23 * k1 + 14; i4 < i2; i4 += 1) {
                A[i0][i2] -= (A[i0][i4] * A[i4][i2]);
              }
              A[i0][i2] /= A[i2][i2];
            }
            if (3 * k1 == k0 + 2) {
              for (register int i2 = i0; i2 <= min(_PB_N - 1, (k0 + 41) / 3); i2 += 1) {
                for (register int i3 = (k0 - 4) / 3; i3 < i0; i3 += 1) {
                  A[i0][i2] -= (A[i0][i3] * A[i3][i2]);
                }
              }
            }
          }
        }
      }
      for (register int k2 = -k1 + (k0 + k1) / 2 + 1; k2 <= min(min(k1 - 1, ((-31 * k1 + 14) / 16) + k0), k1 + floord(_PB_N - k1, 16) - 1); k2 += 1) {
        for (register int i0 = max(k1 - 2, -15 * k1 + 16 * k2 + 15); i0 <= min(_PB_N - 1, k1 + 13); i0 += 1) {
          for (register int i2 = -15 * k1 + 16 * k2 + 14; i2 <= min(-15 * k1 + 16 * k2 + 29, i0 - 1); i2 += 1) {
            for (register int i4 = 16 * k0 - 31 * k1 - 16 * k2 + 14; i4 <= 16 * k0 - 31 * k1 - 16 * k2 + 29; i4 += 1) {
              A[i0][i2] -= (A[i0][i4] * A[i4][i2]);
            }
          }
        }
      }
    }
    if (_PB_N >= 16 * ((k0 - 4) / 46) + 1 && (k0 - 4) % 46 <= 1) {
      for (register int i0 = max(16 * ((k0 - 4) % 46) + 1, 16 * ((k0 - 4) / 46)); i0 <= min(_PB_N - 1, (8 * k0 - 9) / 23 + 14); i0 += 1) {
        if ((k0 - 4) % 46 == 0) {
          for (register int i2 = 0; i2 <= 15; i2 += 1) {
            for (register int i4 = 0; i4 < i2; i4 += 1) {
              A[i0][i2] -= (A[i0][i4] * A[i4][i2]);
            }
            A[i0][i2] /= A[i2][i2];
          }
        } else {
          for (register int i2 = 16; i2 <= min(31, i0 - 1); i2 += 1) {
            for (register int i4 = 0; i4 <= 15; i4 += 1) {
              A[i0][i2] -= (A[i0][i4] * A[i4][i2]);
            }
          }
        }
      }
    }
  }
} else if (_PB_N >= 2) {
  for (register int k0 = max(4, -2 * _PB_N + 10); k0 <= 6; k0 += 1) {
    if (k0 == 4) {
      for (register int i0 = 1; i0 < _PB_N; i0 += 1) {
        for (register int i2 = 0; i2 < i0; i2 += 1) {
          for (register int i4 = 0; i4 < i2; i4 += 1) {
            A[i0][i2] -= (A[i0][i4] * A[i4][i2]);
          }
          A[i0][i2] /= A[i2][i2];
        }
        for (register int i2 = i0; i2 < _PB_N; i2 += 1) {
          for (register int i3 = 0; i3 < i0; i3 += 1) {
            A[i0][i2] -= (A[i0][i3] * A[i3][i2]);
          }
        }
      }
    } else if (_PB_N == 2 && k0 == 6) {
      A[1][0] /= A[0][0];
      A[1][1] -= (A[1][0] * A[0][1]);
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
