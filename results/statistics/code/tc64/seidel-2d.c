/**
 * This version is stamped on Apr. 14, 2015
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* seidel-2d.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "seidel-2d.h"


/* Array initialization. */
static
void init_array (int n,
		 DATA_TYPE POLYBENCH_2D(A,N,N,n,n))
{
  int i, j;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      A[i][j] = ((DATA_TYPE) i*(j+2) + 2) / n;
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
      if ((i * n + j) % 20 == 0) fprintf(POLYBENCH_DUMP_TARGET, "\n");
      fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, A[i][j]);
    }
  POLYBENCH_DUMP_END("A");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_seidel_2d(int tsteps,
		      int n,
		      DATA_TYPE POLYBENCH_2D(A,N,N,n,n))
{
  int t, i, j;

/* TC Optimizing Compiler 0.2.26 */
/* ./tc ../examples/polybench/seidel-2d.scop.c --correction-tiling --lex-scheduling --serial-codegen -b 64 --debug */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#pragma scop
for (int ii0 = 0; ii0 <= floord(_PB_TSTEPS - 1, 64); ii0 += 1) {
  for (int ii1 = 0; ii1 <= floord(_PB_N - 3, 64); ii1 += 1) {
    if (_PB_N >= 64 * ii1 + 67) {
      for (int ii2 = 0; ii2 <= (_PB_N - 3) / 64; ii2 += 1) {
        if (_PB_N >= 64 * ii2 + 67) {
          for (int i0 = 64 * ii0; i0 <= min(_PB_TSTEPS - 1, 64 * ii0 + 31); i0 += 1) {
            for (int i1 = max(1, 64 * ii0 + 64 * ii1 - i0 + 1); i1 <= min(64 * ii0 + 64 * ii1 - i0 + 64, 128 * ii0 + 64 * ii1 + 64 * ii2 - 2 * i0 + 64); i1 += 1) {
              if (i1 >= 64 * ii1 + 1) {
                for (int i2 = max(1, 128 * ii0 + 64 * ii1 + 64 * ii2 - 2 * i0 - i1 + 2); i2 <= 64 * ii1 + 64 * ii2 - i1 + 1; i2 += 1) {
                  A[i1][i2] = (((((((((A[i1 - 1][i2 - 1] + A[i1 - 1][i2]) + A[i1 - 1][i2 + 1]) + A[i1][i2 - 1]) + A[i1][i2]) + A[i1][i2 + 1]) + A[i1 + 1][i2 - 1]) + A[i1 + 1][i2]) + A[i1 + 1][i2 + 1]) / SCALAR_VAL(9.0));
                }
                for (int i2 = max(1, 64 * ii1 + 64 * ii2 - i1 + 2); i2 <= 128 * ii0 + 64 * ii1 + 64 * ii2 - 2 * i0 - i1 + 65; i2 += 1) {
                  A[i1][i2] = (((((((((A[i1 - 1][i2 - 1] + A[i1 - 1][i2]) + A[i1 - 1][i2 + 1]) + A[i1][i2 - 1]) + A[i1][i2]) + A[i1][i2 + 1]) + A[i1 + 1][i2 - 1]) + A[i1 + 1][i2]) + A[i1 + 1][i2 + 1]) / SCALAR_VAL(9.0));
                }
              } else {
                for (int i2 = max(1, 128 * ii0 + 64 * ii1 + 64 * ii2 - 2 * i0 - i1 + 2); i2 <= 128 * ii0 + 64 * ii1 + 64 * ii2 - 2 * i0 - i1 + 65; i2 += 1) {
                  A[i1][i2] = (((((((((A[i1 - 1][i2 - 1] + A[i1 - 1][i2]) + A[i1 - 1][i2 + 1]) + A[i1][i2 - 1]) + A[i1][i2]) + A[i1][i2 + 1]) + A[i1 + 1][i2 - 1]) + A[i1 + 1][i2]) + A[i1 + 1][i2 + 1]) / SCALAR_VAL(9.0));
                }
              }
            }
          }
          if (63 * ii1 + 32 * ii2 >= 32) {
            for (int i0 = 64 * ii0 + 32; i0 <= min(_PB_TSTEPS - 1, 64 * ii0 + 63); i0 += 1) {
              for (int i1 = max(1, 64 * ii0 + 64 * ii1 - i0 + 1); i1 <= min(64 * ii0 + 64 * ii1 - i0 + 64, 128 * ii0 + 64 * ii1 + 64 * ii2 - 2 * i0 + 64); i1 += 1) {
                for (int i2 = max(1, 128 * ii0 + 64 * ii1 + 64 * ii2 - 2 * i0 - i1 + 2); i2 <= 128 * ii0 + 64 * ii1 + 64 * ii2 - 2 * i0 - i1 + 65; i2 += 1) {
                  A[i1][i2] = (((((((((A[i1 - 1][i2 - 1] + A[i1 - 1][i2]) + A[i1 - 1][i2 + 1]) + A[i1][i2 - 1]) + A[i1][i2]) + A[i1][i2 + 1]) + A[i1 + 1][i2 - 1]) + A[i1 + 1][i2]) + A[i1 + 1][i2 + 1]) / SCALAR_VAL(9.0));
                }
              }
            }
          }
        } else {
          for (int i0 = 64 * ii0; i0 <= min(_PB_TSTEPS - 1, 64 * ii0 + 63); i0 += 1) {
            for (int i1 = max(1, 64 * ii0 + 64 * ii1 - i0 + 1); i1 <= 64 * ii0 + 64 * ii1 - i0 + 64; i1 += 1) {
              if (i1 >= 64 * ii1 + 1) {
                for (int i2 = max(1, 128 * ii0 + 64 * ii1 + 64 * ii2 - 2 * i0 - i1 + 2); i2 <= 64 * ii1 + 64 * ii2 - i1 + 1; i2 += 1) {
                  A[i1][i2] = (((((((((A[i1 - 1][i2 - 1] + A[i1 - 1][i2]) + A[i1 - 1][i2 + 1]) + A[i1][i2 - 1]) + A[i1][i2]) + A[i1][i2 + 1]) + A[i1 + 1][i2 - 1]) + A[i1 + 1][i2]) + A[i1 + 1][i2 + 1]) / SCALAR_VAL(9.0));
                }
                for (int i2 = 64 * ii1 + 64 * ii2 - i1 + 2; i2 < _PB_N - 1; i2 += 1) {
                  A[i1][i2] = (((((((((A[i1 - 1][i2 - 1] + A[i1 - 1][i2]) + A[i1 - 1][i2 + 1]) + A[i1][i2 - 1]) + A[i1][i2]) + A[i1][i2 + 1]) + A[i1 + 1][i2 - 1]) + A[i1 + 1][i2]) + A[i1 + 1][i2 + 1]) / SCALAR_VAL(9.0));
                }
              } else {
                for (int i2 = 128 * ii0 + 64 * ii1 + 64 * ii2 - 2 * i0 - i1 + 2; i2 < _PB_N - 1; i2 += 1) {
                  A[i1][i2] = (((((((((A[i1 - 1][i2 - 1] + A[i1 - 1][i2]) + A[i1 - 1][i2 + 1]) + A[i1][i2 - 1]) + A[i1][i2]) + A[i1][i2 + 1]) + A[i1 + 1][i2 - 1]) + A[i1 + 1][i2]) + A[i1 + 1][i2 + 1]) / SCALAR_VAL(9.0));
                }
              }
            }
          }
        }
      }
    } else if (_PB_TSTEPS >= 64 * ii0 + 2 && ii1 >= 1) {
      for (int ii2 = 0; ii2 <= ii1; ii2 += 1) {
        if (_PB_N >= 64 * ii2 + 67) {
          for (int i0 = 64 * ii0; i0 <= min(_PB_TSTEPS - 1, 64 * ii0 + 63); i0 += 1) {
            for (int i1 = 64 * ii0 + 64 * ii1 - i0 + 1; i1 <= min(_PB_N - 2, 128 * ii0 + 64 * ii1 + 64 * ii2 - 2 * i0 + 64); i1 += 1) {
              if (64 * ii1 >= i1) {
                for (int i2 = max(1, 128 * ii0 + 64 * ii1 + 64 * ii2 - 2 * i0 - i1 + 2); i2 <= 128 * ii0 + 64 * ii1 + 64 * ii2 - 2 * i0 - i1 + 65; i2 += 1) {
                  A[i1][i2] = (((((((((A[i1 - 1][i2 - 1] + A[i1 - 1][i2]) + A[i1 - 1][i2 + 1]) + A[i1][i2 - 1]) + A[i1][i2]) + A[i1][i2 + 1]) + A[i1 + 1][i2 - 1]) + A[i1 + 1][i2]) + A[i1 + 1][i2 + 1]) / SCALAR_VAL(9.0));
                }
              } else {
                for (int i2 = max(1, 128 * ii0 + 64 * ii1 + 64 * ii2 - 2 * i0 - i1 + 2); i2 <= min(64 * ii1 + 64 * ii2 - i1 + 1, 128 * ii0 + 64 * ii1 + 64 * ii2 - 2 * i0 - i1 + 65); i2 += 1) {
                  A[i1][i2] = (((((((((A[i1 - 1][i2 - 1] + A[i1 - 1][i2]) + A[i1 - 1][i2 + 1]) + A[i1][i2 - 1]) + A[i1][i2]) + A[i1][i2 + 1]) + A[i1 + 1][i2 - 1]) + A[i1 + 1][i2]) + A[i1 + 1][i2 + 1]) / SCALAR_VAL(9.0));
                }
                for (int i2 = max(1, 64 * ii1 + 64 * ii2 - i1 + 2); i2 <= 128 * ii0 + 64 * ii1 + 64 * ii2 - 2 * i0 - i1 + 65; i2 += 1) {
                  A[i1][i2] = (((((((((A[i1 - 1][i2 - 1] + A[i1 - 1][i2]) + A[i1 - 1][i2 + 1]) + A[i1][i2 - 1]) + A[i1][i2]) + A[i1][i2 + 1]) + A[i1 + 1][i2 - 1]) + A[i1 + 1][i2]) + A[i1 + 1][i2 + 1]) / SCALAR_VAL(9.0));
                }
              }
            }
          }
        } else {
          for (int i0 = 64 * ii0; i0 <= min(_PB_TSTEPS - 1, 64 * ii0 + 63); i0 += 1) {
            for (int i1 = 64 * ii0 + 64 * ii1 - i0 + 1; i1 < _PB_N - 1; i1 += 1) {
              if (64 * ii1 >= i1) {
                for (int i2 = max(1, 128 * ii0 + 128 * ii1 - 2 * i0 - i1 + 2); i2 < _PB_N - 1; i2 += 1) {
                  A[i1][i2] = (((((((((A[i1 - 1][i2 - 1] + A[i1 - 1][i2]) + A[i1 - 1][i2 + 1]) + A[i1][i2 - 1]) + A[i1][i2]) + A[i1][i2 + 1]) + A[i1 + 1][i2 - 1]) + A[i1 + 1][i2]) + A[i1 + 1][i2 + 1]) / SCALAR_VAL(9.0));
                }
              } else {
                for (int i2 = max(1, 128 * ii0 + 128 * ii1 - 2 * i0 - i1 + 2); i2 <= 128 * ii1 - i1 + 1; i2 += 1) {
                  A[i1][i2] = (((((((((A[i1 - 1][i2 - 1] + A[i1 - 1][i2]) + A[i1 - 1][i2 + 1]) + A[i1][i2 - 1]) + A[i1][i2]) + A[i1][i2 + 1]) + A[i1 + 1][i2 - 1]) + A[i1 + 1][i2]) + A[i1 + 1][i2 + 1]) / SCALAR_VAL(9.0));
                }
                for (int i2 = 128 * ii1 - i1 + 2; i2 < _PB_N - 1; i2 += 1) {
                  A[i1][i2] = (((((((((A[i1 - 1][i2 - 1] + A[i1 - 1][i2]) + A[i1 - 1][i2 + 1]) + A[i1][i2 - 1]) + A[i1][i2]) + A[i1][i2 + 1]) + A[i1 + 1][i2 - 1]) + A[i1 + 1][i2]) + A[i1 + 1][i2 + 1]) / SCALAR_VAL(9.0));
                }
              }
            }
          }
        }
      }
    } else if (64 * ii0 + 1 == _PB_TSTEPS) {
      for (int ii2 = 0; ii2 < ii1; ii2 += 1) {
        for (int i1 = 64 * ii1 + 1; i1 < _PB_N - 1; i1 += 1) {
          for (int i2 = max(1, 64 * ii1 + 64 * ii2 - i1 + 2); i2 <= 64 * ii1 + 64 * ii2 - i1 + 65; i2 += 1) {
            A[i1][i2] = (((((((((A[i1 - 1][i2 - 1] + A[i1 - 1][i2]) + A[i1 - 1][i2 + 1]) + A[i1][i2 - 1]) + A[i1][i2]) + A[i1][i2 + 1]) + A[i1 + 1][i2 - 1]) + A[i1 + 1][i2]) + A[i1 + 1][i2 + 1]) / SCALAR_VAL(9.0));
          }
        }
      }
      for (int i1 = 64 * ii1 + 1; i1 < _PB_N - 1; i1 += 1) {
        for (int i2 = max(1, 128 * ii1 - i1 + 2); i2 < _PB_N - 1; i2 += 1) {
          A[i1][i2] = (((((((((A[i1 - 1][i2 - 1] + A[i1 - 1][i2]) + A[i1 - 1][i2 + 1]) + A[i1][i2 - 1]) + A[i1][i2]) + A[i1][i2 + 1]) + A[i1 + 1][i2 - 1]) + A[i1 + 1][i2]) + A[i1 + 1][i2 + 1]) / SCALAR_VAL(9.0));
        }
      }
    } else {
      for (int i0 = 64 * ii0; i0 <= min(_PB_TSTEPS - 1, 64 * ii0 + 63); i0 += 1) {
        for (int i1 = 1; i1 < _PB_N - 1; i1 += 1) {
          for (int i2 = 1; i2 < _PB_N - 1; i2 += 1) {
            A[i1][i2] = (((((((((A[i1 - 1][i2 - 1] + A[i1 - 1][i2]) + A[i1 - 1][i2 + 1]) + A[i1][i2 - 1]) + A[i1][i2]) + A[i1][i2 + 1]) + A[i1 + 1][i2 - 1]) + A[i1 + 1][i2]) + A[i1 + 1][i2 + 1]) / SCALAR_VAL(9.0));
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
  int tsteps = TSTEPS;

  /* Variable declaration/allocation. */
  POLYBENCH_2D_ARRAY_DECL(A, DATA_TYPE, N, N, n, n);


  /* Initialize array(s). */
  init_array (n, POLYBENCH_ARRAY(A));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_seidel_2d (tsteps, n, POLYBENCH_ARRAY(A));

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
