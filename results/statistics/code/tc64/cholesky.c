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


/* TC Optimizing Compiler 0.2.26 */
/* ./tc ../examples/polybench/cholesky.scop.c --correction-tiling --lex-scheduling --serial-codegen -b 64 --debug --iterative-tc */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#pragma scop
for (int ii0 = 0; ii0 <= floord(_PB_N - 1, 64); ii0 += 1) {
  for (int ii2 = 0; ii2 <= min(ii0, (_PB_N - ii0 - 1) / 63 - 1); ii2 += 1) {
    for (int ii3 = max(0, -ii2 + 1); ii3 <= min(1, ii0 - ii2); ii3 += 1) {
      if (ii3 == 1) {
        if (_PB_N >= 64 * ii0 + 64) {
          for (int i0 = 64 * ii0; i0 <= 64 * ii0 + 63; i0 += 1) {
            if (ii2 == 0) {
              A[i0][0] /= A[0][0];
            }
            for (int i2 = max(1, 64 * ii2); i2 <= 64 * ii2 + 63; i2 += 1) {
              if (i2 >= 64) {
                for (int i4 = 64 * ii2; i4 < i2; i4 += 1) {
                  A[i0][i2] -= (A[i0][i4] * A[i2][i4]);
                }
              } else {
                for (int i4 = 0; i4 < i2; i4 += 1) {
                  A[i0][i2] -= (A[i0][i4] * A[i2][i4]);
                }
              }
              A[i0][i2] /= A[i2][i2];
            }
          }
        } else {
          for (int i0 = 64 * ii0; i0 < _PB_N; i0 += 1) {
            for (int i2 = 64 * ii2; i2 <= 64 * ii2 + 63; i2 += 1) {
              for (int i4 = 64 * ii2; i4 < i2; i4 += 1) {
                A[i0][i2] -= (A[i0][i4] * A[i2][i4]);
              }
              A[i0][i2] /= A[i2][i2];
            }
          }
        }
      } else {
        for (int ii4 = 0; ii4 < ii2; ii4 += 1) {
          if (_PB_N >= 64 * ii0 + 64) {
            for (int i0 = max(64 * ii0, 64 * ii2 + 1); i0 <= 64 * ii0 + 63; i0 += 1) {
              for (int i2 = 64 * ii2; i2 <= min(64 * ii2 + 63, i0 - 1); i2 += 1) {
                for (int i4 = 64 * ii4; i4 <= 64 * ii4 + 63; i4 += 1) {
                  A[i0][i2] -= (A[i0][i4] * A[i2][i4]);
                }
              }
            }
          } else {
            for (int i0 = 64 * ii0; i0 < _PB_N; i0 += 1) {
              for (int i2 = 64 * ii2; i2 <= 64 * ii2 + 63; i2 += 1) {
                for (int i4 = 64 * ii4; i4 <= 64 * ii4 + 63; i4 += 1) {
                  A[i0][i2] -= (A[i0][i4] * A[i2][i4]);
                }
              }
            }
          }
        }
      }
    }
  }
  if (64 * ii0 + 63 >= _PB_N) {
    for (int ii4 = 0; ii4 < ii0; ii4 += 1) {
      for (int i0 = 64 * ii0 + 1; i0 < _PB_N; i0 += 1) {
        for (int i2 = 64 * ii0; i2 < i0; i2 += 1) {
          for (int i4 = 64 * ii4; i4 <= 64 * ii4 + 63; i4 += 1) {
            A[i0][i2] -= (A[i0][i4] * A[i2][i4]);
          }
        }
      }
    }
  }
  for (int ii1 = max(1, -ii0 + 2); ii1 <= 2; ii1 += 1) {
    if (ii1 == 2) {
      if (_PB_N >= 64 * ii0 + 64) {
        for (int i0 = max(64, 64 * ii0); i0 <= 64 * ii0 + 63; i0 += 1) {
          for (int i2 = 64 * ii0; i2 < i0; i2 += 1) {
            for (int i4 = 64 * ii0; i4 < i2; i4 += 1) {
              A[i0][i2] -= (A[i0][i4] * A[i2][i4]);
            }
            A[i0][i2] /= A[i2][i2];
          }
          for (int i2 = 64 * ii0; i2 < i0; i2 += 1) {
            A[i0][i0] -= (A[i0][i2] * A[i0][i2]);
          }
          A[i0][i0] = SQRT_FUN(A[i0][i0]);
        }
        if (ii0 == 0) {
          A[0][0] = SQRT_FUN(A[0][0]);
          for (int i0 = 1; i0 <= 63; i0 += 1) {
            A[i0][0] /= A[0][0];
            for (int i2 = 1; i2 < i0; i2 += 1) {
              for (int i4 = 0; i4 < i2; i4 += 1) {
                A[i0][i2] -= (A[i0][i4] * A[i2][i4]);
              }
              A[i0][i2] /= A[i2][i2];
            }
            for (int i2 = 0; i2 < i0; i2 += 1) {
              A[i0][i0] -= (A[i0][i2] * A[i0][i2]);
            }
            A[i0][i0] = SQRT_FUN(A[i0][i0]);
          }
        }
      } else {
        for (int i0 = 64 * ii0; i0 < _PB_N; i0 += 1) {
          for (int i2 = 64 * ii0; i2 < i0; i2 += 1) {
            for (int i4 = 64 * ii0; i4 < i2; i4 += 1) {
              A[i0][i2] -= (A[i0][i4] * A[i2][i4]);
            }
            A[i0][i2] /= A[i2][i2];
          }
          for (int i2 = 64 * ii0; i2 < i0; i2 += 1) {
            A[i0][i0] -= (A[i0][i2] * A[i0][i2]);
          }
          A[i0][i0] = SQRT_FUN(A[i0][i0]);
        }
      }
    } else {
      for (int ii2 = 0; ii2 < ii0; ii2 += 1) {
        if (_PB_N >= 64 * ii0 + 64) {
          for (int i0 = 64 * ii0; i0 <= 64 * ii0 + 63; i0 += 1) {
            for (int i2 = 64 * ii2; i2 <= 64 * ii2 + 63; i2 += 1) {
              A[i0][i0] -= (A[i0][i2] * A[i0][i2]);
            }
          }
        } else {
          for (int i0 = 64 * ii0; i0 < _PB_N; i0 += 1) {
            for (int i2 = 64 * ii2; i2 <= 64 * ii2 + 63; i2 += 1) {
              A[i0][i0] -= (A[i0][i2] * A[i0][i2]);
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
