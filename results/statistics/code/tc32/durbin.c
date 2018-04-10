/**
 * This version is stamped on Apr. 14, 2015
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* durbin.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "durbin.h"


/* Array initialization. */
static
void init_array (int n,
		 DATA_TYPE POLYBENCH_1D(r,N,n))
{
  int i, j;

  for (i = 0; i < n; i++)
    {
      r[i] = (n+1-i);
    }
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_1D(y,N,n))

{
  int i;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("y");
  for (i = 0; i < n; i++) {
    if (i % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
    fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, y[i]);
  }
  POLYBENCH_DUMP_END("y");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_durbin(int n,
		   DATA_TYPE POLYBENCH_1D(r,N,n),
		   DATA_TYPE POLYBENCH_1D(y,N,n))
{
 DATA_TYPE z[N];
 DATA_TYPE alpha;
 DATA_TYPE beta;
 DATA_TYPE sum;

 int i,k;

 y[0] = -r[0];
 beta = SCALAR_VAL(1.0);
 alpha = -r[0];

/* TC Optimizing Compiler 0.2.26 */
/* ./tc ../examples/polybench/durbin.scop.c --correction-tiling --lex-scheduling --serial-codegen -b 32 --debug --isl-union-map-tc */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#pragma scop
for (int ii0 = 0; ii0 <= floord(_PB_N - 2, 32); ii0 += 1) {
  for (int ii1 = 0; ii1 <= 6; ii1 += 1) {
    if (ii1 <= 1) {
      if (ii1 == 1) {
        sum = SCALAR_VAL(0.0);
      } else {
        beta = ((1 - (alpha * alpha)) * beta);
      }
    } else if (ii1 >= 4 && ii1 <= 5) {
      for (int ii2 = 0; ii2 <= ii0; ii2 += 1) {
        if (ii1 == 5) {
          for (int i2 = 32 * ii2; i2 <= min(32 * ii0, 32 * ii2 + 31); i2 += 1) {
            y[i2] = z[i2];
          }
          if (_PB_N >= 32 * ii0 + 3) {
            for (int i2 = 32 * ii2; i2 <= min(32 * ii0, 32 * ii2 + 31); i2 += 1) {
              sum += (r[32 * ii0 - i2 + 1] * y[i2]);
            }
          }
        } else {
          for (int i2 = 32 * ii2; i2 <= min(32 * ii0, 32 * ii2 + 31); i2 += 1) {
            z[i2] = (y[i2] + (alpha * y[32 * ii0 - i2]));
          }
        }
      }
    } else if (ii1 == 6) {
      if (_PB_N >= 32 * ii0 + 33) {
        for (int i0 = 32 * ii0 + 1; i0 <= 32 * ii0 + 32; i0 += 1) {
          if (i0 >= 32 * ii0 + 3) {
            beta = ((1 - (alpha * alpha)) * beta);
            sum = SCALAR_VAL(0.0);
          }
          if (i0 >= 32 * ii0 + 2) {
            if (i0 >= 32 * ii0 + 3) {
              sum += (r[i0 - 1] * y[0]);
              for (int i2 = 1; i2 <= 32 * ii0; i2 += 1) {
                sum += (r[i0 - i2 - 1] * y[i2]);
              }
            }
            for (int i2 = 32 * ii0 + 1; i2 < i0; i2 += 1) {
              sum += (r[i0 - i2 - 1] * y[i2]);
            }
            alpha = ((-(r[i0] + sum)) / beta);
            for (int i2 = 0; i2 < i0; i2 += 1) {
              z[i2] = (y[i2] + (alpha * y[i0 - i2 - 1]));
            }
            for (int i2 = 0; i2 < i0; i2 += 1) {
              y[i2] = z[i2];
            }
          }
          y[i0] = alpha;
        }
      } else {
        for (int i0 = 32 * ii0 + 1; i0 < _PB_N; i0 += 1) {
          if (i0 >= 32 * ii0 + 3) {
            beta = ((1 - (alpha * alpha)) * beta);
            sum = SCALAR_VAL(0.0);
          }
          if (i0 >= 32 * ii0 + 2) {
            if (i0 >= 32 * ii0 + 3) {
              for (int i2 = 0; i2 <= 32 * ii0; i2 += 1) {
                sum += (r[i0 - i2 - 1] * y[i2]);
              }
            }
            for (int i2 = 32 * ii0 + 1; i2 < i0; i2 += 1) {
              sum += (r[i0 - i2 - 1] * y[i2]);
            }
            alpha = ((-(r[i0] + sum)) / beta);
            for (int i2 = 0; i2 < i0; i2 += 1) {
              z[i2] = (y[i2] + (alpha * y[i0 - i2 - 1]));
            }
            for (int i2 = 0; i2 < i0; i2 += 1) {
              y[i2] = z[i2];
            }
          }
          y[i0] = alpha;
        }
      }
    } else if (ii1 == 3) {
      alpha = ((-(r[32 * ii0 + 1] + sum)) / beta);
      if (_PB_N >= 32 * ii0 + 3) {
        beta = ((1 - (alpha * alpha)) * beta);
        sum = SCALAR_VAL(0.0);
      }
    } else {
      for (int ii2 = 0; ii2 <= ii0; ii2 += 1) {
        if (32 * ii0 + 32 >= _PB_N) {
          for (int i2 = 32 * ii2; i2 <= min(32 * ii0, 32 * ii2 + 31); i2 += 1) {
            sum += (r[32 * ii0 - i2] * y[i2]);
          }
        } else if (ii0 >= ii2 + 1) {
          for (int i2 = 32 * ii2; i2 <= 32 * ii2 + 31; i2 += 1) {
            sum += (r[32 * ii0 - i2] * y[i2]);
          }
        } else {
          sum += (r[0] * y[32 * ii0]);
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
  POLYBENCH_1D_ARRAY_DECL(r, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(y, DATA_TYPE, N, n);


  /* Initialize array(s). */
  init_array (n, POLYBENCH_ARRAY(r));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_durbin (n,
		 POLYBENCH_ARRAY(r),
		 POLYBENCH_ARRAY(y));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(y)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(r);
  POLYBENCH_FREE_ARRAY(y);

  return 0;
}
