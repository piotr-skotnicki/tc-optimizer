/**
 * This version is stamped on Apr. 14, 2015
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
  DATA_TYPE temp2;

//BLAS PARAMS
//SIDE = 'L'
//UPLO = 'L'
// =>  Form  C := alpha*A*B + beta*C
// A is MxM
// B is MxN
// C is MxN
//note that due to Fortran array layout, the code below more closely resembles upper triangular case in BLAS
/* ./tc ../examples/polybench/symm.scop.c --merge-tiling --free-scheduling --omp-for-codegen --debug -b 64 -D M=2000 -D N=2600 */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#pragma scop
for (int k = 0; k <= 32; k += 1) {
  #pragma omp parallel for
  for (int ii1 = 0; ii1 <= 40; ii1 += 1)
    for (int ii3 = 0; ii3 < k; ii3 += 1) {
      if (k <= 31) {
        for (int c6 = max(64 * k - 64, 64 * ii3 + 1); c6 < 64 * k; c6 += 1)
          for (int c7 = 64 * ii1; c7 <= min(2599, 64 * ii1 + 63); c7 += 1)
            for (int c9 = 64 * ii3; c9 <= min(64 * ii3 + 63, c6 - 1); c9 += 1)
              C[c9][c7] += ((alpha * B[c6][c7]) * A[c6][c9]);
      } else
        for (int c6 = max(1984, 64 * ii3 + 1); c6 <= 1999; c6 += 1) {
          if (ii1 <= 39) {
            for (int c7 = 64 * ii1; c7 <= 64 * ii1 + 63; c7 += 1)
              for (int c9 = 64 * ii3; c9 <= min(64 * ii3 + 63, c6 - 1); c9 += 1)
                C[c9][c7] += ((alpha * B[c6][c7]) * A[c6][c9]);
          } else
            for (int c7 = 2560; c7 <= 2599; c7 += 1)
              for (int c9 = 64 * ii3; c9 <= min(64 * ii3 + 63, c6 - 1); c9 += 1)
                C[c9][c7] += ((alpha * B[c6][c7]) * A[c6][c9]);
        }
    }
  if (k <= 31) {
    if (k <= 30) {
      for (int c6 = 64 * k; c6 <= 64 * k + 63; c6 += 1)
        for (int c7 = 0; c7 <= 2599; c7 += 1) {
          temp2 = 0;
          for (int c9 = 0; c9 < c6; c9 += 1)
            temp2 += (B[c9][c7] * A[c6][c9]);
          C[c6][c7] = (((beta * C[c6][c7]) + ((alpha * B[c6][c7]) * A[c6][c6])) + (alpha * temp2));
        }
    } else
      for (int c6 = 1984; c6 <= 1999; c6 += 1)
        for (int c7 = 0; c7 <= 2599; c7 += 1) {
          temp2 = 0;
          for (int c9 = 0; c9 < c6; c9 += 1)
            temp2 += (B[c9][c7] * A[c6][c9]);
          C[c6][c7] = (((beta * C[c6][c7]) + ((alpha * B[c6][c7]) * A[c6][c6])) + (alpha * temp2));
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
