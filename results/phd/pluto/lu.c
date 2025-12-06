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

#include <omp.h>
#include <math.h>
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

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

  int t1, t2, t3, t4, t5, t6;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
if (_PB_N >= 2) {
  for (t1=0;t1<=floord(_PB_N-1,8);t1++) {
    lbp=max(0,ceild(16*t1-_PB_N+1,16));
    ubp=min(floord(_PB_N-1,16),t1);
#pragma omp parallel for private(lbv,ubv,t3,t4,t5,t6)
    for (t2=lbp;t2<=ubp;t2++) {
      for (t3=0;t3<=min(min(floord(_PB_N-2,16),t2),t1-t2);t3++) {
        if (t2 == t3) {
          for (t4=max(16*t1-16*t2,16*t2+16);t4<=(min(_PB_N-1,16*t1-16*t2+15))-7;t4+=8) {
            for (t5=16*t2;t5<=16*t2+14;t5++) {
              A[t4][t5] /= A[t5][t5];;
              A[(t4+1)][t5] /= A[t5][t5];;
              A[(t4+2)][t5] /= A[t5][t5];;
              A[(t4+3)][t5] /= A[t5][t5];;
              A[(t4+4)][t5] /= A[t5][t5];;
              A[(t4+5)][t5] /= A[t5][t5];;
              A[(t4+6)][t5] /= A[t5][t5];;
              A[(t4+7)][t5] /= A[t5][t5];;
              for (t6=t5+1;t6<=16*t2+15;t6++) {
                A[t4][t6] -= A[t4][t5] * A[t5][t6];;
                A[(t4+1)][t6] -= A[(t4+1)][t5] * A[t5][t6];;
                A[(t4+2)][t6] -= A[(t4+2)][t5] * A[t5][t6];;
                A[(t4+3)][t6] -= A[(t4+3)][t5] * A[t5][t6];;
                A[(t4+4)][t6] -= A[(t4+4)][t5] * A[t5][t6];;
                A[(t4+5)][t6] -= A[(t4+5)][t5] * A[t5][t6];;
                A[(t4+6)][t6] -= A[(t4+6)][t5] * A[t5][t6];;
                A[(t4+7)][t6] -= A[(t4+7)][t5] * A[t5][t6];;
              }
            }
            A[t4][(16*t2+15)] /= A[(16*t2+15)][(16*t2+15)];;
            A[(t4+1)][(16*t2+15)] /= A[(16*t2+15)][(16*t2+15)];;
            A[(t4+2)][(16*t2+15)] /= A[(16*t2+15)][(16*t2+15)];;
            A[(t4+3)][(16*t2+15)] /= A[(16*t2+15)][(16*t2+15)];;
            A[(t4+4)][(16*t2+15)] /= A[(16*t2+15)][(16*t2+15)];;
            A[(t4+5)][(16*t2+15)] /= A[(16*t2+15)][(16*t2+15)];;
            A[(t4+6)][(16*t2+15)] /= A[(16*t2+15)][(16*t2+15)];;
            A[(t4+7)][(16*t2+15)] /= A[(16*t2+15)][(16*t2+15)];;
          }
          for (;t4<=min(_PB_N-1,16*t1-16*t2+15);t4++) {
            for (t5=16*t2;t5<=16*t2+14;t5++) {
              A[t4][t5] /= A[t5][t5];;
              for (t6=t5+1;t6<=16*t2+15;t6++) {
                A[t4][t6] -= A[t4][t5] * A[t5][t6];;
              }
            }
            A[t4][(16*t2+15)] /= A[(16*t2+15)][(16*t2+15)];;
          }
        }
        if (t2 >= t3+1) {
          for (t4=max(16*t1-16*t2,16*t2+16);t4<=(min(_PB_N-1,16*t1-16*t2+15))-7;t4+=8) {
            for (t5=16*t3;t5<=16*t3+15-7;t5+=8) {
              for (t6=16*t2;t6<=16*t2+15;t6++) {
                A[t4][t6] -= A[t4][t5] * A[t5][t6];;
                A[(t4+1)][t6] -= A[(t4+1)][t5] * A[t5][t6];;
                A[(t4+2)][t6] -= A[(t4+2)][t5] * A[t5][t6];;
                A[(t4+3)][t6] -= A[(t4+3)][t5] * A[t5][t6];;
                A[(t4+4)][t6] -= A[(t4+4)][t5] * A[t5][t6];;
                A[(t4+5)][t6] -= A[(t4+5)][t5] * A[t5][t6];;
                A[(t4+6)][t6] -= A[(t4+6)][t5] * A[t5][t6];;
                A[(t4+7)][t6] -= A[(t4+7)][t5] * A[t5][t6];;
                A[t4][t6] -= A[t4][(t5+1)] * A[(t5+1)][t6];;
                A[(t4+1)][t6] -= A[(t4+1)][(t5+1)] * A[(t5+1)][t6];;
                A[(t4+2)][t6] -= A[(t4+2)][(t5+1)] * A[(t5+1)][t6];;
                A[(t4+3)][t6] -= A[(t4+3)][(t5+1)] * A[(t5+1)][t6];;
                A[(t4+4)][t6] -= A[(t4+4)][(t5+1)] * A[(t5+1)][t6];;
                A[(t4+5)][t6] -= A[(t4+5)][(t5+1)] * A[(t5+1)][t6];;
                A[(t4+6)][t6] -= A[(t4+6)][(t5+1)] * A[(t5+1)][t6];;
                A[(t4+7)][t6] -= A[(t4+7)][(t5+1)] * A[(t5+1)][t6];;
                A[t4][t6] -= A[t4][(t5+2)] * A[(t5+2)][t6];;
                A[(t4+1)][t6] -= A[(t4+1)][(t5+2)] * A[(t5+2)][t6];;
                A[(t4+2)][t6] -= A[(t4+2)][(t5+2)] * A[(t5+2)][t6];;
                A[(t4+3)][t6] -= A[(t4+3)][(t5+2)] * A[(t5+2)][t6];;
                A[(t4+4)][t6] -= A[(t4+4)][(t5+2)] * A[(t5+2)][t6];;
                A[(t4+5)][t6] -= A[(t4+5)][(t5+2)] * A[(t5+2)][t6];;
                A[(t4+6)][t6] -= A[(t4+6)][(t5+2)] * A[(t5+2)][t6];;
                A[(t4+7)][t6] -= A[(t4+7)][(t5+2)] * A[(t5+2)][t6];;
                A[t4][t6] -= A[t4][(t5+3)] * A[(t5+3)][t6];;
                A[(t4+1)][t6] -= A[(t4+1)][(t5+3)] * A[(t5+3)][t6];;
                A[(t4+2)][t6] -= A[(t4+2)][(t5+3)] * A[(t5+3)][t6];;
                A[(t4+3)][t6] -= A[(t4+3)][(t5+3)] * A[(t5+3)][t6];;
                A[(t4+4)][t6] -= A[(t4+4)][(t5+3)] * A[(t5+3)][t6];;
                A[(t4+5)][t6] -= A[(t4+5)][(t5+3)] * A[(t5+3)][t6];;
                A[(t4+6)][t6] -= A[(t4+6)][(t5+3)] * A[(t5+3)][t6];;
                A[(t4+7)][t6] -= A[(t4+7)][(t5+3)] * A[(t5+3)][t6];;
                A[t4][t6] -= A[t4][(t5+4)] * A[(t5+4)][t6];;
                A[(t4+1)][t6] -= A[(t4+1)][(t5+4)] * A[(t5+4)][t6];;
                A[(t4+2)][t6] -= A[(t4+2)][(t5+4)] * A[(t5+4)][t6];;
                A[(t4+3)][t6] -= A[(t4+3)][(t5+4)] * A[(t5+4)][t6];;
                A[(t4+4)][t6] -= A[(t4+4)][(t5+4)] * A[(t5+4)][t6];;
                A[(t4+5)][t6] -= A[(t4+5)][(t5+4)] * A[(t5+4)][t6];;
                A[(t4+6)][t6] -= A[(t4+6)][(t5+4)] * A[(t5+4)][t6];;
                A[(t4+7)][t6] -= A[(t4+7)][(t5+4)] * A[(t5+4)][t6];;
                A[t4][t6] -= A[t4][(t5+5)] * A[(t5+5)][t6];;
                A[(t4+1)][t6] -= A[(t4+1)][(t5+5)] * A[(t5+5)][t6];;
                A[(t4+2)][t6] -= A[(t4+2)][(t5+5)] * A[(t5+5)][t6];;
                A[(t4+3)][t6] -= A[(t4+3)][(t5+5)] * A[(t5+5)][t6];;
                A[(t4+4)][t6] -= A[(t4+4)][(t5+5)] * A[(t5+5)][t6];;
                A[(t4+5)][t6] -= A[(t4+5)][(t5+5)] * A[(t5+5)][t6];;
                A[(t4+6)][t6] -= A[(t4+6)][(t5+5)] * A[(t5+5)][t6];;
                A[(t4+7)][t6] -= A[(t4+7)][(t5+5)] * A[(t5+5)][t6];;
                A[t4][t6] -= A[t4][(t5+6)] * A[(t5+6)][t6];;
                A[(t4+1)][t6] -= A[(t4+1)][(t5+6)] * A[(t5+6)][t6];;
                A[(t4+2)][t6] -= A[(t4+2)][(t5+6)] * A[(t5+6)][t6];;
                A[(t4+3)][t6] -= A[(t4+3)][(t5+6)] * A[(t5+6)][t6];;
                A[(t4+4)][t6] -= A[(t4+4)][(t5+6)] * A[(t5+6)][t6];;
                A[(t4+5)][t6] -= A[(t4+5)][(t5+6)] * A[(t5+6)][t6];;
                A[(t4+6)][t6] -= A[(t4+6)][(t5+6)] * A[(t5+6)][t6];;
                A[(t4+7)][t6] -= A[(t4+7)][(t5+6)] * A[(t5+6)][t6];;
                A[t4][t6] -= A[t4][(t5+7)] * A[(t5+7)][t6];;
                A[(t4+1)][t6] -= A[(t4+1)][(t5+7)] * A[(t5+7)][t6];;
                A[(t4+2)][t6] -= A[(t4+2)][(t5+7)] * A[(t5+7)][t6];;
                A[(t4+3)][t6] -= A[(t4+3)][(t5+7)] * A[(t5+7)][t6];;
                A[(t4+4)][t6] -= A[(t4+4)][(t5+7)] * A[(t5+7)][t6];;
                A[(t4+5)][t6] -= A[(t4+5)][(t5+7)] * A[(t5+7)][t6];;
                A[(t4+6)][t6] -= A[(t4+6)][(t5+7)] * A[(t5+7)][t6];;
                A[(t4+7)][t6] -= A[(t4+7)][(t5+7)] * A[(t5+7)][t6];;
              }
            }
            for (;t5<=16*t3+15;t5++) {
              for (t6=16*t2;t6<=16*t2+15;t6++) {
                A[t4][t6] -= A[t4][t5] * A[t5][t6];;
                A[(t4+1)][t6] -= A[(t4+1)][t5] * A[t5][t6];;
                A[(t4+2)][t6] -= A[(t4+2)][t5] * A[t5][t6];;
                A[(t4+3)][t6] -= A[(t4+3)][t5] * A[t5][t6];;
                A[(t4+4)][t6] -= A[(t4+4)][t5] * A[t5][t6];;
                A[(t4+5)][t6] -= A[(t4+5)][t5] * A[t5][t6];;
                A[(t4+6)][t6] -= A[(t4+6)][t5] * A[t5][t6];;
                A[(t4+7)][t6] -= A[(t4+7)][t5] * A[t5][t6];;
              }
            }
          }
          for (;t4<=min(_PB_N-1,16*t1-16*t2+15);t4++) {
            for (t5=16*t3;t5<=16*t3+15;t5++) {
              for (t6=16*t2;t6<=16*t2+15;t6++) {
                A[t4][t6] -= A[t4][t5] * A[t5][t6];;
              }
            }
          }
        }
        if ((t1 == 2*t2) && (t1 == 2*t3)) {
          if (t1%2 == 0) {
            A[(8*t1+1)][8*t1] /= A[8*t1][8*t1];;
          }
          for (t6=8*t1+1;t6<=min(_PB_N-1,8*t1+15);t6++) {
            if (t1%2 == 0) {
              A[(8*t1+1)][t6] -= A[(8*t1+1)][8*t1] * A[8*t1][t6];;
            }
          }
        }
        for (t4=max(16*t1-16*t2,16*t3+1);t4<=min(16*t2,16*t1-16*t2+15);t4++) {
          for (t5=16*t3;t5<=(min(16*t3+15,t4-1))-7;t5+=8) {
            for (t6=16*t2;t6<=min(_PB_N-1,16*t2+15);t6++) {
              A[t4][t6] -= A[t4][t5] * A[t5][t6];;
              A[t4][t6] -= A[t4][(t5+1)] * A[(t5+1)][t6];;
              A[t4][t6] -= A[t4][(t5+2)] * A[(t5+2)][t6];;
              A[t4][t6] -= A[t4][(t5+3)] * A[(t5+3)][t6];;
              A[t4][t6] -= A[t4][(t5+4)] * A[(t5+4)][t6];;
              A[t4][t6] -= A[t4][(t5+5)] * A[(t5+5)][t6];;
              A[t4][t6] -= A[t4][(t5+6)] * A[(t5+6)][t6];;
              A[t4][t6] -= A[t4][(t5+7)] * A[(t5+7)][t6];;
            }
          }
          for (;t5<=min(16*t3+15,t4-1);t5++) {
            for (t6=16*t2;t6<=min(_PB_N-1,16*t2+15);t6++) {
              A[t4][t6] -= A[t4][t5] * A[t5][t6];;
            }
          }
        }
        if ((t1 == 2*t2) && (t1 == 2*t3)) {
          for (t4=8*t1+2;t4<=min(_PB_N-1,8*t1+15);t4++) {
            for (t5=8*t1;t5<=t4-2;t5++) {
              if (t1%2 == 0) {
                A[t4][t5] /= A[t5][t5];;
              }
              for (t6=t5+1;t6<=t4-1;t6++) {
                if (t1%2 == 0) {
                  A[t4][t6] -= A[t4][t5] * A[t5][t6];;
                }
              }
              for (t6=t4;t6<=min(_PB_N-1,8*t1+15);t6++) {
                if (t1%2 == 0) {
                  A[t4][t6] -= A[t4][t5] * A[t5][t6];;
                }
              }
            }
            if (t1%2 == 0) {
              A[t4][(t4-1)] /= A[(t4-1)][(t4-1)];;
            }
            for (t6=t4;t6<=min(_PB_N-1,8*t1+15);t6++) {
              if (t1%2 == 0) {
                A[t4][t6] -= A[t4][(t4-1)] * A[(t4-1)][t6];;
              }
            }
          }
        }
        if ((t1 == 2*t2) && (t1 >= 2*t3+2)) {
          for (t4=8*t1+1;t4<=min(_PB_N-1,8*t1+15);t4++) {
            for (t5=16*t3;t5<=16*t3+15-7;t5+=8) {
              for (t6=8*t1;t6<=t4-1;t6++) {
                if (t1%2 == 0) {
                  A[t4][t6] -= A[t4][t5] * A[t5][t6];;
                  A[t4][t6] -= A[t4][(t5+1)] * A[(t5+1)][t6];;
                  A[t4][t6] -= A[t4][(t5+2)] * A[(t5+2)][t6];;
                  A[t4][t6] -= A[t4][(t5+3)] * A[(t5+3)][t6];;
                  A[t4][t6] -= A[t4][(t5+4)] * A[(t5+4)][t6];;
                  A[t4][t6] -= A[t4][(t5+5)] * A[(t5+5)][t6];;
                  A[t4][t6] -= A[t4][(t5+6)] * A[(t5+6)][t6];;
                  A[t4][t6] -= A[t4][(t5+7)] * A[(t5+7)][t6];;
                }
              }
              for (t6=t4;t6<=min(_PB_N-1,8*t1+15);t6++) {
                if (t1%2 == 0) {
                  A[t4][t6] -= A[t4][t5] * A[t5][t6];;
                  A[t4][t6] -= A[t4][(t5+1)] * A[(t5+1)][t6];;
                  A[t4][t6] -= A[t4][(t5+2)] * A[(t5+2)][t6];;
                  A[t4][t6] -= A[t4][(t5+3)] * A[(t5+3)][t6];;
                  A[t4][t6] -= A[t4][(t5+4)] * A[(t5+4)][t6];;
                  A[t4][t6] -= A[t4][(t5+5)] * A[(t5+5)][t6];;
                  A[t4][t6] -= A[t4][(t5+6)] * A[(t5+6)][t6];;
                  A[t4][t6] -= A[t4][(t5+7)] * A[(t5+7)][t6];;
                }
              }
            }
            for (;t5<=16*t3+15;t5++) {
              for (t6=8*t1;t6<=t4-1;t6++) {
                if (t1%2 == 0) {
                  A[t4][t6] -= A[t4][t5] * A[t5][t6];;
                }
              }
              for (t6=t4;t6<=min(_PB_N-1,8*t1+15);t6++) {
                if (t1%2 == 0) {
                  A[t4][t6] -= A[t4][t5] * A[t5][t6];;
                }
              }
            }
          }
        }
      }
    }
  }
}
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
