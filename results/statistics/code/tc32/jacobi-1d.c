/**
 * This version is stamped on Apr. 14, 2015
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* jacobi-1d.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "jacobi-1d.h"


/* Array initialization. */
static
void init_array (int n,
		 DATA_TYPE POLYBENCH_1D(A,N,n),
		 DATA_TYPE POLYBENCH_1D(B,N,n))
{
  int i;

  for (i = 0; i < n; i++)
      {
	A[i] = ((DATA_TYPE) i+ 2) / n;
	B[i] = ((DATA_TYPE) i+ 3) / n;
      }
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_1D(A,N,n))

{
  int i;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("A");
  for (i = 0; i < n; i++)
    {
      if (i % 20 == 0) fprintf(POLYBENCH_DUMP_TARGET, "\n");
      fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, A[i]);
    }
  POLYBENCH_DUMP_END("A");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_jacobi_1d(int tsteps,
			    int n,
			    DATA_TYPE POLYBENCH_1D(A,N,n),
			    DATA_TYPE POLYBENCH_1D(B,N,n))
{
  int t, i;

/* TC Optimizing Compiler 0.2.26 */
/* ./tc ../examples/polybench/jacobi-1d.scop.c --correction-tiling --lex-scheduling --serial-codegen -b 32 --debug */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#pragma scop
for (int ii0 = 0; ii0 <= floord(_PB_TSTEPS - 1, 32); ii0 += 1) {
  for (int ii1 = 0; ii1 <= min(min(1, _PB_TSTEPS - 32 * ii0 - 1), (_PB_N + 28) / 33); ii1 += 1) {
    if (ii1 == 1) {
      for (int ii2 = 0; ii2 <= (_PB_N - 3) / 32; ii2 += 1) {
        if (_PB_N >= 32 * ii2 + 35) {
          for (int i0 = 32 * ii0; i0 <= min(min(_PB_TSTEPS - 1, 32 * ii0 + 31), 32 * ii0 + 16 * ii2 + 16); i0 += 1) {
            if (i0 >= 32 * ii0 + 1) {
              for (int i2 = max(1, 64 * ii0 + 32 * ii2 - 2 * i0 + 2); i2 <= 64 * ii0 + 32 * ii2 - 2 * i0 + 33; i2 += 1) {
                B[i2] = (0.33333 * ((A[i2 - 1] + A[i2]) + A[i2 + 1]));
              }
            }
            for (int i2 = max(1, 64 * ii0 + 32 * ii2 - 2 * i0 + 1); i2 <= 64 * ii0 + 32 * ii2 - 2 * i0 + 32; i2 += 1) {
              A[i2] = (0.33333 * ((B[i2 - 1] + B[i2]) + B[i2 + 1]));
            }
          }
        } else {
          for (int i0 = 32 * ii0; i0 <= min(_PB_TSTEPS - 1, 32 * ii0 + 31); i0 += 1) {
            if (i0 >= 32 * ii0 + 1) {
              for (int i2 = max(1, 64 * ii0 + 32 * ii2 - 2 * i0 + 2); i2 < 32 * ii2; i2 += 1) {
                B[i2] = (0.33333 * ((A[i2 - 1] + A[i2]) + A[i2 + 1]));
              }
              for (int i2 = max(1, 32 * ii2); i2 < _PB_N - 1; i2 += 1) {
                B[i2] = (0.33333 * ((A[i2 - 1] + A[i2]) + A[i2 + 1]));
              }
            }
            for (int i2 = max(1, 64 * ii0 + 32 * ii2 - 2 * i0 + 1); i2 < _PB_N - 1; i2 += 1) {
              A[i2] = (0.33333 * ((B[i2 - 1] + B[i2]) + B[i2 + 1]));
            }
          }
        }
      }
    } else {
      for (int ii2 = 0; ii2 <= floord(_PB_N - 3, 32); ii2 += 1) {
        for (int i2 = 32 * ii2 + 1; i2 <= min(_PB_N - 2, 32 * ii2 + 32); i2 += 1) {
          B[i2] = (0.33333 * ((A[i2 - 1] + A[i2]) + A[i2 + 1]));
        }
      }
    }
  }
  if (_PB_N >= 5 && 32 * ii0 + 1 == _PB_TSTEPS) {
    for (int ii2 = 0; ii2 <= (_PB_N - 3) / 32; ii2 += 1) {
      for (int i2 = 32 * ii2 + 1; i2 <= min(_PB_N - 2, 32 * ii2 + 32); i2 += 1) {
        A[i2] = (0.33333 * ((B[i2 - 1] + B[i2]) + B[i2 + 1]));
      }
    }
  } else if (_PB_N <= 4) {
    for (int i0 = 32 * ii0; i0 <= min(_PB_TSTEPS - 1, 32 * ii0 + 31); i0 += 1) {
      if (i0 >= 32 * ii0 + 1) {
        for (int i2 = 1; i2 < _PB_N - 1; i2 += 1) {
          B[i2] = (0.33333 * ((A[i2 - 1] + A[i2]) + A[i2 + 1]));
        }
      }
      for (int i2 = 1; i2 < _PB_N - 1; i2 += 1) {
        A[i2] = (0.33333 * ((B[i2 - 1] + B[i2]) + B[i2 + 1]));
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
  POLYBENCH_1D_ARRAY_DECL(A, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(B, DATA_TYPE, N, n);


  /* Initialize array(s). */
  init_array (n, POLYBENCH_ARRAY(A), POLYBENCH_ARRAY(B));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_jacobi_1d(tsteps, n, POLYBENCH_ARRAY(A), POLYBENCH_ARRAY(B));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(A)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(B);

  return 0;
}
