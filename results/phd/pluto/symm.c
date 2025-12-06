/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* symm.c: this file is part of PolyBench/C */

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
#include "symm.h"


/* Array initialization. */
static
void init_array(int m, int n,
		DATA_TYPE *alpha,
		DATA_TYPE *beta,
		DATA_TYPE POLYBENCH_2D(C,M,N,m,n),
		DATA_TYPE POLYBENCH_2D(A,M,M,m,m),
		DATA_TYPE POLYBENCH_2D(B,M,N,m,n))
{
  int i, j;

  *alpha = 1.5;
  *beta = 1.2;
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++) {
      C[i][j] = (DATA_TYPE) ((i+j) % 100) / m;
      B[i][j] = (DATA_TYPE) ((n+i-j) % 100) / m;
    }
  for (i = 0; i < m; i++) {
    for (j = 0; j <=i; j++)
      A[i][j] = (DATA_TYPE) ((i+j) % 100) / m;
    for (j = i+1; j < m; j++)
      A[i][j] = -999; //regions of arrays that should not be used
  }
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int m, int n,
		 DATA_TYPE POLYBENCH_2D(C,M,N,m,n))
{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("C");
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++) {
	if ((i * m + j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
	fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, C[i][j]);
    }
  POLYBENCH_DUMP_END("C");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_symm(int m, int n,
		 DATA_TYPE alpha,
		 DATA_TYPE beta,
		 DATA_TYPE POLYBENCH_2D(C,M,N,m,n),
		 DATA_TYPE POLYBENCH_2D(A,M,M,m,m),
		 DATA_TYPE POLYBENCH_2D(B,M,N,m,n))
{
  int i, j, k;
  DATA_TYPE temp2[_PB_N];

//BLAS PARAMS
//SIDE = 'L'
//UPLO = 'L'
// =>  Form  C := alpha*A*B + beta*C
// A is MxM
// B is MxN
// C is MxN
//note that due to Fortran array layout, the code below more closely resembles upper triangular case in BLAS

  int t1, t2, t3, t4, t5, t6;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
if ((_PB_M >= 1) && (_PB_N >= 1)) {
  lbp=0;
  ubp=floord(_PB_N-1,32);
#pragma omp parallel for private(lbv,ubv,t2,t3,t4,t5,t6)
  for (t1=lbp;t1<=ubp;t1++) {
    for (t2=0;t2<=floord(_PB_M-1,32);t2++) {
      if (t2 == 0) {
        lbv=32*t1;
        ubv=min(_PB_N-1,32*t1+31);
#pragma ivdep
#pragma vector always
        for (t5=lbv;t5<=ubv;t5++) {
          temp2[t5] = 0;;
          C[0][t5] = beta * C[0][t5] + alpha*B[0][t5] * A[0][0] + alpha * temp2[t5];;
        }
      }
      if ((_PB_M >= 2) && (t2 == 0)) {
        lbv=32*t1;
        ubv=min(_PB_N-1,32*t1+31);
#pragma ivdep
#pragma vector always
        for (t5=lbv;t5<=ubv;t5++) {
          C[0][t5] += alpha*B[1][t5] * A[1][0];;
          temp2[t5] = 0;;
          temp2[t5] += B[0][t5] * A[1][0];;
        }
        lbv=32*t1;
        ubv=min(_PB_N-1,32*t1+31);
#pragma ivdep
#pragma vector always
        for (t5=lbv;t5<=ubv;t5++) {
          C[1][t5] = beta * C[1][t5] + alpha*B[1][t5] * A[1][1] + alpha * temp2[t5];;
        }
      }
      for (t3=max(2,32*t2);t3<=min(_PB_M-1,32*t2+31);t3++) {
        lbv=32*t1;
        ubv=min(_PB_N-1,32*t1+31);
#pragma ivdep
#pragma vector always
        for (t5=lbv;t5<=ubv;t5++) {
          C[0][t5] += alpha*B[t3][t5] * A[t3][0];;
          temp2[t5] = 0;;
          temp2[t5] += B[0][t5] * A[t3][0];;
        }
        for (t4=1;t4<=t3-1-7;t4+=8) {
          lbv=32*t1;
          ubv=min(_PB_N-1,32*t1+31);
#pragma ivdep
#pragma vector always
          for (t5=lbv;t5<=ubv;t5++) {
            C[t4][t5] += alpha*B[t3][t5] * A[t3][t4];;
            temp2[t5] += B[t4][t5] * A[t3][t4];;
            C[(t4+1)][t5] += alpha*B[t3][t5] * A[t3][(t4+1)];;
            temp2[t5] += B[(t4+1)][t5] * A[t3][(t4+1)];;
            C[(t4+2)][t5] += alpha*B[t3][t5] * A[t3][(t4+2)];;
            temp2[t5] += B[(t4+2)][t5] * A[t3][(t4+2)];;
            C[(t4+3)][t5] += alpha*B[t3][t5] * A[t3][(t4+3)];;
            temp2[t5] += B[(t4+3)][t5] * A[t3][(t4+3)];;
            C[(t4+4)][t5] += alpha*B[t3][t5] * A[t3][(t4+4)];;
            temp2[t5] += B[(t4+4)][t5] * A[t3][(t4+4)];;
            C[(t4+5)][t5] += alpha*B[t3][t5] * A[t3][(t4+5)];;
            temp2[t5] += B[(t4+5)][t5] * A[t3][(t4+5)];;
            C[(t4+6)][t5] += alpha*B[t3][t5] * A[t3][(t4+6)];;
            temp2[t5] += B[(t4+6)][t5] * A[t3][(t4+6)];;
            C[(t4+7)][t5] += alpha*B[t3][t5] * A[t3][(t4+7)];;
            temp2[t5] += B[(t4+7)][t5] * A[t3][(t4+7)];;
          }
        }
        for (;t4<=t3-1;t4++) {
          lbv=32*t1;
          ubv=min(_PB_N-1,32*t1+31);
#pragma ivdep
#pragma vector always
          for (t5=lbv;t5<=ubv;t5++) {
            C[t4][t5] += alpha*B[t3][t5] * A[t3][t4];;
            temp2[t5] += B[t4][t5] * A[t3][t4];;
          }
        }
        lbv=32*t1;
        ubv=min(_PB_N-1,32*t1+31);
#pragma ivdep
#pragma vector always
        for (t5=lbv;t5<=ubv;t5++) {
          C[t3][t5] = beta * C[t3][t5] + alpha*B[t3][t5] * A[t3][t3] + alpha * temp2[t5];;
        }
      }
    }
  }
}

}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int m = M;
  int n = N;

  /* Variable declaration/allocation. */
  DATA_TYPE alpha;
  DATA_TYPE beta;
  POLYBENCH_2D_ARRAY_DECL(C,DATA_TYPE,M,N,m,n);
  POLYBENCH_2D_ARRAY_DECL(A,DATA_TYPE,M,M,m,m);
  POLYBENCH_2D_ARRAY_DECL(B,DATA_TYPE,M,N,m,n);

  /* Initialize array(s). */
  init_array (m, n, &alpha, &beta,
	      POLYBENCH_ARRAY(C),
	      POLYBENCH_ARRAY(A),
	      POLYBENCH_ARRAY(B));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_symm (m, n,
	       alpha, beta,
	       POLYBENCH_ARRAY(C),
	       POLYBENCH_ARRAY(A),
	       POLYBENCH_ARRAY(B));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(m, n, POLYBENCH_ARRAY(C)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(C);
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(B);

  return 0;
}
