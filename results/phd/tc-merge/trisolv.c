/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* trisolv.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "trisolv.h"


/* Array initialization. */
static
void init_array(int n,
		DATA_TYPE POLYBENCH_2D(L,N,N,n,n),
		DATA_TYPE POLYBENCH_1D(x,N,n),
		DATA_TYPE POLYBENCH_1D(b,N,n))
{
  int i, j;

  for (i = 0; i < n; i++)
    {
      x[i] = - 999;
      b[i] =  i ;
      for (j = 0; j <= i; j++)
	L[i][j] = (DATA_TYPE) (i+n-j+1)*2/n;
    }
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_1D(x,N,n))

{
  int i;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("x");
  for (i = 0; i < n; i++) {
    fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, x[i]);
    if (i % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
  }
  POLYBENCH_DUMP_END("x");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_trisolv(int n,
		    DATA_TYPE POLYBENCH_2D(L,N,N,n,n),
		    DATA_TYPE POLYBENCH_1D(x,N,n),
		    DATA_TYPE POLYBENCH_1D(b,N,n))
{
  int i, j;

/* TC Optimizing Compiler 0.5.1 */
/* ./tc ../examples/polybench/trisolv.scop.c --merge-tiling --free-scheduling --omp-for-codegen --iterative-tc --debug -b 64 --align */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#pragma scop
{
  #pragma omp parallel for
  for (register int ii0 = 0; ii0 <= floord(_PB_N - 1, 64); ii0 += 1) {
    for (register int i0 = 64 * ii0; i0 <= min(_PB_N - 1, 64 * ii0 + 63); i0 += 1) {
      x[i0] = b[i0];
    }
  }
  if (_PB_N >= 2) {
    for (register int i0 = 0; i0 <= min(63, _PB_N - 1); i0 += 1) {
      for (register int i2 = 0; i2 < i0; i2 += 1) {
        x[i0] -= (L[i0][i2] * x[i2]);
      }
      x[i0] = (x[i0] / L[i0][i0]);
    }
  }
  for (register int k = 2; k <= min(1022, floord(_PB_N - 1, 32)); k += 1) {
    if ((k + 1) % 2 == 0) {
      for (register int i0 = 32 * k - 32; i0 <= min(_PB_N - 1, 32 * k + 31); i0 += 1) {
        for (register int i2 = 32 * k - 32; i2 < i0; i2 += 1) {
          x[i0] -= (L[i0][i2] * x[i2]);
        }
        x[i0] = (x[i0] / L[i0][i0]);
      }
    } else {
      #pragma omp parallel for
      for (register int ii0 = k / 2; ii0 <= (_PB_N - 1) / 64; ii0 += 1) {
        for (register int i0 = 64 * ii0; i0 <= min(_PB_N - 1, 64 * ii0 + 63); i0 += 1) {
          for (register int i2 = 32 * k - 64; i2 < 32 * k; i2 += 1) {
            x[i0] -= (L[i0][i2] * x[i2]);
          }
        }
      }
    }
  }
  if (_PB_N >= 66 && _PB_N <= 32672 && (_PB_N - 33) % 64 >= 33) {
    for (register int i0 = -((_PB_N - 2) % 64) + _PB_N - 2; i0 < _PB_N; i0 += 1) {
      for (register int i2 = -((_PB_N - 2) % 64) + _PB_N - 2; i2 < i0; i2 += 1) {
        x[i0] -= (L[i0][i2] * x[i2]);
      }
      x[i0] = (x[i0] / L[i0][i0]);
    }
  } else if ((_PB_N - 1) % 64 == 0) {
    for (register int k = 1023; k <= min((_PB_N - 1) / 32, (_PB_N + 32767) / 64); k += 1) {
      {
        #pragma omp parallel for
        for (register int ii2 = max(k - 513, (k + 1) / 2 - 1); ii2 < min(k - 511, (_PB_N - 1) / 64); ii2 += 1) {
          for (register int i0 = 64 * k - 32768; i0 < min(_PB_N, 64 * k - 32704); i0 += 1) {
            for (register int i2 = 64 * ii2; i2 <= min(64 * ii2 + 63, i0 - 1); i2 += 1) {
              x[i0] -= (L[i0][i2] * x[i2]);
            }
            if (ii2 + 512 == k) {
              x[i0] = (x[i0] / L[i0][i0]);
            }
          }
        }
        if (64 * k == _PB_N + 32767) {
          x[_PB_N - 1] = (x[_PB_N - 1] / L[_PB_N - 1][_PB_N - 1]);
        }
      }
      if (k >= 1024) {
        #pragma omp parallel for
        for (register int ii0 = k - 511; ii0 <= (_PB_N - 1) / 64; ii0 += 1) {
          for (register int i0 = 64 * ii0; i0 <= min(_PB_N - 1, 64 * ii0 + 63); i0 += 1) {
            for (register int i2 = 64 * k - 32832; i2 < 64 * k - 32768; i2 += 1) {
              x[i0] -= (L[i0][i2] * x[i2]);
            }
          }
        }
      }
    }
  }
  if ((_PB_N + 63) % 64 >= 1) {
    for (register int k = 1023; k <= (_PB_N + 63) / 64 + 511; k += 1) {
      #pragma omp parallel for
      for (register int ii2 = max(k - 513, (k + 1) / 2 - 1); ii2 < k - 511; ii2 += 1) {
        for (register int i0 = 64 * k - 32768; i0 < min(_PB_N, 64 * k - 32704); i0 += 1) {
          for (register int i2 = 64 * ii2; i2 <= min(64 * ii2 + 63, i0 - 1); i2 += 1) {
            x[i0] -= (L[i0][i2] * x[i2]);
          }
          if (ii2 + 512 == k) {
            x[i0] = (x[i0] / L[i0][i0]);
          }
        }
      }
      if (k >= 1024) {
        #pragma omp parallel for
        for (register int ii0 = k - 511; ii0 <= (_PB_N - 1) / 64; ii0 += 1) {
          for (register int i0 = 64 * ii0; i0 <= min(_PB_N - 1, 64 * ii0 + 63); i0 += 1) {
            for (register int i2 = 64 * k - 32832; i2 < 64 * k - 32768; i2 += 1) {
              x[i0] -= (L[i0][i2] * x[i2]);
            }
          }
        }
      }
    }
  } else if (_PB_N <= 32705) {
    x[_PB_N - 1] = (x[_PB_N - 1] / L[_PB_N - 1][_PB_N - 1]);
  }
}
#pragma endscop

}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int n = N;

  /* Variable declaration/allocation. */
  POLYBENCH_2D_ARRAY_DECL(L, DATA_TYPE, N, N, n, n);
  POLYBENCH_1D_ARRAY_DECL(x, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(b, DATA_TYPE, N, n);


  /* Initialize array(s). */
  init_array (n, POLYBENCH_ARRAY(L), POLYBENCH_ARRAY(x), POLYBENCH_ARRAY(b));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_trisolv (n, POLYBENCH_ARRAY(L), POLYBENCH_ARRAY(x), POLYBENCH_ARRAY(b));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(x)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(L);
  POLYBENCH_FREE_ARRAY(x);
  POLYBENCH_FREE_ARRAY(b);

  return 0;
}
