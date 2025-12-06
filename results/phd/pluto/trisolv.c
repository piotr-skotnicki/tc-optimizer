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

  int t1, t2, t3, t4, t5;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
if (_PB_N >= 1) {
  lbp=0;
  ubp=floord(_PB_N-1,64);
#pragma omp parallel for private(lbv,ubv,t3,t4,t5)
  for (t2=lbp;t2<=ubp;t2++) {
    lbv=64*t2;
    ubv=min(_PB_N-1,64*t2+63);
#pragma ivdep
#pragma vector always
    for (t3=lbv;t3<=ubv;t3++) {
      x[t3] = b[t3];;
    }
  }
  for (t2=0;t2<=floord(_PB_N-1,32);t2++) {
    lbp=max(0,ceild(64*t2-_PB_N+1,64));
    ubp=floord(t2,2);
#pragma omp parallel for private(lbv,ubv,t4,t5)
    for (t3=lbp;t3<=ubp;t3++) {
      if (t2 >= 2*t3+1) {
        for (t4=64*t2-64*t3;t4<=(min(_PB_N-1,64*t2-64*t3+63))-7;t4+=8) {
          for (t5=64*t3;t5<=64*t3+63;t5++) {
            x[t4] -= L[t4][t5] * x[t5];;
            x[(t4+1)] -= L[(t4+1)][t5] * x[t5];;
            x[(t4+2)] -= L[(t4+2)][t5] * x[t5];;
            x[(t4+3)] -= L[(t4+3)][t5] * x[t5];;
            x[(t4+4)] -= L[(t4+4)][t5] * x[t5];;
            x[(t4+5)] -= L[(t4+5)][t5] * x[t5];;
            x[(t4+6)] -= L[(t4+6)][t5] * x[t5];;
            x[(t4+7)] -= L[(t4+7)][t5] * x[t5];;
          }
        }
        for (;t4<=min(_PB_N-1,64*t2-64*t3+63);t4++) {
          for (t5=64*t3;t5<=64*t3+63;t5++) {
            x[t4] -= L[t4][t5] * x[t5];;
          }
        }
      }
      if (t2 == 2*t3) {
        if (t2%2 == 0) {
          x[32*t2] = x[32*t2] / L[32*t2][32*t2];;
        }
      }
      if (t2 == 2*t3) {
        for (t4=32*t2+1;t4<=min(_PB_N-1,32*t2+63);t4++) {
          for (t5=32*t2;t5<=t4-1;t5++) {
            if (t2%2 == 0) {
              x[t4] -= L[t4][t5] * x[t5];;
            }
          }
          if (t2%2 == 0) {
            x[t4] = x[t4] / L[t4][t4];;
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
