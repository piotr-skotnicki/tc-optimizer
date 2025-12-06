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

/* TC Optimizing Compiler 0.5.1 */
/* ./tc ../examples/polybench/symm-exp.scop.c --correction-tiling --isl-wave-scheduling --omp-for-codegen --iterative-tc --debug -b 32 --align */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#pragma scop
for (register int k0 = 0; k0 <= floord(_PB_M + _PB_N - 2, 32); k0 += 1) {
  #pragma omp parallel for
  for (register int k1 = max(0, k0 - (_PB_M + 31) / 32 + 1); k1 <= min(k0, floord(_PB_N - 1, 32)); k1 += 1) {
    if (k1 == k0) {
      for (register int i1 = 32 * k0; i1 <= min(_PB_N - 1, 32 * k0 + 31); i1 += 1) {
        temp2[i1] = 0;
      }
    }
    for (register int k3 = 0; k3 < k0 - k1; k3 += 1) {
      if (k3 == 0) {
        for (register int i1 = 32 * k1; i1 <= min(_PB_N - 1, 32 * k1 + 31); i1 += 1) {
          temp2[i1] = 0;
        }
      }
      for (register int i1 = 32 * k1; i1 <= min(_PB_N - 1, 32 * k1 + 31); i1 += 1) {
        for (register int i3 = 32 * k3; i3 <= 32 * k3 + 31; i3 += 1) {
          temp2[i1] += (B[i3][i1] * A[32 * k0 - 32 * k1][i3]);
        }
      }
    }
    for (register int i0 = 32 * k0 - 32 * k1; i0 <= min(_PB_M - 1, 32 * k0 - 32 * k1 + 31); i0 += 1) {
      for (register int i1 = 32 * k1; i1 <= min(_PB_N - 1, 32 * k1 + 31); i1 += 1) {
        if (32 * k1 + i0 >= 32 * k0 + 1) {
          temp2[i1] = 0;
          for (register int i3 = 0; i3 < i0; i3 += 1) {
            if (32 * k1 + i3 >= 32 * k0) {
              C[i3][i1] += ((alpha * B[i0][i1]) * A[i0][i3]);
            }
            temp2[i1] += (B[i3][i1] * A[i0][i3]);
          }
        }
        C[i0][i1] = (((beta * C[i0][i1]) + ((alpha * B[i0][i1]) * A[i0][i0])) + (alpha * temp2[i1]));
      }
    }
    for (register int k3 = _PB_M + k0 - k1 + 1; k3 <= _PB_M + (_PB_M - 1) / 32; k3 += 1) {
      for (register int i0 = -32 * _PB_M + 32 * k3; i0 <= min(_PB_M - 1, -32 * _PB_M + 32 * k3 + 31); i0 += 1) {
        for (register int i1 = 32 * k1; i1 <= min(_PB_N - 1, 32 * k1 + 31); i1 += 1) {
          for (register int i3 = 32 * k0 - 32 * k1; i3 <= 32 * k0 - 32 * k1 + 31; i3 += 1) {
            C[i3][i1] += ((alpha * B[i0][i1]) * A[i0][i3]);
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
