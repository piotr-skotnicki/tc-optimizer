/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* floyd-warshall.c: this file is part of PolyBench/C */

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
#include "floyd-warshall.h"


/* Array initialization. */
static
void init_array (int n,
		 DATA_TYPE POLYBENCH_2D(path,N,N,n,n))
{
  int i, j;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
      path[i][j] = i*j%7+1;
      if ((i+j)%13 == 0 || (i+j)%7==0 || (i+j)%11 == 0)
         path[i][j] = 999;
    }
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_2D(path,N,N,n,n))

{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("path");
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
      if ((i * n + j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
      fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, path[i][j]);
    }
  POLYBENCH_DUMP_END("path");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_floyd_warshall(int n,
			   DATA_TYPE POLYBENCH_2D(path,N,N,n,n))
{
  int i, j, k;

  int t1, t2, t3, t4, t5;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
if (_PB_N >= 1) {
  for (t1=0;t1<=_PB_N-1;t1++) {
    for (t2=0;t2<=floord(_PB_N-1,8);t2++) {
      lbp=max(0,ceild(16*t2-_PB_N+1,16));
      ubp=min(floord(_PB_N-1,16),t2);
#pragma omp parallel for private(lbv,ubv,t4,t5)
      for (t3=lbp;t3<=ubp;t3++) {
        for (t4=16*t2-16*t3;t4<=(min(_PB_N-1,16*t2-16*t3+15))-7;t4+=8) {
          for (t5=16*t3;t5<=min(_PB_N-1,16*t3+15);t5++) {
            path[t4][t5] = path[t4][t5] < path[t4][t1] + path[t1][t5] ? path[t4][t5] : path[t4][t1] + path[t1][t5];;
            path[(t4+1)][t5] = path[(t4+1)][t5] < path[(t4+1)][t1] + path[t1][t5] ? path[(t4+1)][t5] : path[(t4+1)][t1] + path[t1][t5];;
            path[(t4+2)][t5] = path[(t4+2)][t5] < path[(t4+2)][t1] + path[t1][t5] ? path[(t4+2)][t5] : path[(t4+2)][t1] + path[t1][t5];;
            path[(t4+3)][t5] = path[(t4+3)][t5] < path[(t4+3)][t1] + path[t1][t5] ? path[(t4+3)][t5] : path[(t4+3)][t1] + path[t1][t5];;
            path[(t4+4)][t5] = path[(t4+4)][t5] < path[(t4+4)][t1] + path[t1][t5] ? path[(t4+4)][t5] : path[(t4+4)][t1] + path[t1][t5];;
            path[(t4+5)][t5] = path[(t4+5)][t5] < path[(t4+5)][t1] + path[t1][t5] ? path[(t4+5)][t5] : path[(t4+5)][t1] + path[t1][t5];;
            path[(t4+6)][t5] = path[(t4+6)][t5] < path[(t4+6)][t1] + path[t1][t5] ? path[(t4+6)][t5] : path[(t4+6)][t1] + path[t1][t5];;
            path[(t4+7)][t5] = path[(t4+7)][t5] < path[(t4+7)][t1] + path[t1][t5] ? path[(t4+7)][t5] : path[(t4+7)][t1] + path[t1][t5];;
          }
        }
        for (;t4<=min(_PB_N-1,16*t2-16*t3+15);t4++) {
          for (t5=16*t3;t5<=min(_PB_N-1,16*t3+15);t5++) {
            path[t4][t5] = path[t4][t5] < path[t4][t1] + path[t1][t5] ? path[t4][t5] : path[t4][t1] + path[t1][t5];;
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
  POLYBENCH_2D_ARRAY_DECL(path, DATA_TYPE, N, N, n, n);


  /* Initialize array(s). */
  init_array (n, POLYBENCH_ARRAY(path));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_floyd_warshall (n, POLYBENCH_ARRAY(path));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(path)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(path);

  return 0;
}
