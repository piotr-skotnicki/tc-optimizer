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
/* TC Optimizing Compiler 0.2.26 */
/* ./tc ../examples/polybench/symm.scop.c --correction-tiling --lex-scheduling --serial-codegen -b 64 --debug */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#pragma scop
for (int ii0 = 0; ii0 <= floord(_PB_M - 1, 64); ii0 += 1) {
  for (int ii1 = 0; ii1 <= floord(_PB_N - 1, 64); ii1 += 1) {
    temp2 = 0;
    for (int ii2 = max(1, -ii0 + 2); ii2 <= 2; ii2 += 1) {
      if (ii2 == 2) {
        if (_PB_N >= 64 * ii1 + 65) {
          for (int i1 = 64 * ii1; i1 <= 64 * ii1 + 63; i1 += 1) {
            if (i1 >= 64 * ii1 + 1) {
              temp2 = 0;
              for (int i3 = 0; i3 < 64 * ii0; i3 += 1) {
                temp2 += (B[i3][i1] * A[64 * ii0][i3]);
              }
            }
            C[64 * ii0][i1] = (((beta * C[64 * ii0][i1]) + ((alpha * B[64 * ii0][i1]) * A[64 * ii0][64 * ii0])) + (alpha * temp2));
          }
        } else {
          for (int i1 = 64 * ii1; i1 < _PB_N; i1 += 1) {
            if (i1 >= 64 * ii1 + 1) {
              temp2 = 0;
              for (int i3 = 0; i3 < 64 * ii0; i3 += 1) {
                temp2 += (B[i3][i1] * A[64 * ii0][i3]);
              }
            }
            C[64 * ii0][i1] = (((beta * C[64 * ii0][i1]) + ((alpha * B[64 * ii0][i1]) * A[64 * ii0][64 * ii0])) + (alpha * temp2));
          }
        }
        for (int i0 = 64 * ii0 + 1; i0 <= min(_PB_M - 1, 64 * ii0 + 63); i0 += 1) {
          if (64 * ii1 + 64 >= _PB_N) {
            for (int i1 = 0; i1 < 64 * ii1; i1 += 1) {
              temp2 = 0;
              for (int i3 = 0; i3 < i0; i3 += 1) {
                if (i3 >= 64 * ii0 + 1) {
                  C[i3][i1] += ((alpha * B[i0][i1]) * A[i0][i3]);
                }
                temp2 += (B[i3][i1] * A[i0][i3]);
              }
              C[i0][i1] = (((beta * C[i0][i1]) + ((alpha * B[i0][i1]) * A[i0][i0])) + (alpha * temp2));
            }
          }
          for (int i1 = 64 * ii1; i1 <= min(_PB_N - 1, 64 * ii1 + 63); i1 += 1) {
            if (64 * ii1 + 64 >= _PB_N) {
              temp2 = 0;
            }
            if (_PB_N >= 64 * ii1 + 65) {
              C[64 * ii0][i1] += ((alpha * B[i0][i1]) * A[i0][64 * ii0]);
            } else {
              for (int i3 = 0; i3 < i0; i3 += 1) {
                if (i3 >= 64 * ii0) {
                  C[i3][i1] += ((alpha * B[i0][i1]) * A[i0][i3]);
                }
                temp2 += (B[i3][i1] * A[i0][i3]);
              }
            }
            if (64 * ii1 + 64 >= _PB_N) {
              C[i0][i1] = (((beta * C[i0][i1]) + ((alpha * B[i0][i1]) * A[i0][i0])) + (alpha * temp2));
            }
          }
        }
      } else {
        for (int ii3 = 0; ii3 < ii0; ii3 += 1) {
          for (int i0 = 64 * ii0; i0 <= min(_PB_M - 1, 64 * ii0 + 63); i0 += 1) {
            for (int i1 = 64 * ii1; i1 <= min(_PB_N - 1, 64 * ii1 + 63); i1 += 1) {
              for (int i3 = 64 * ii3; i3 <= 64 * ii3 + 63; i3 += 1) {
                C[i3][i1] += ((alpha * B[i0][i1]) * A[i0][i3]);
              }
            }
          }
          for (int i3 = 64 * ii3; i3 <= 64 * ii3 + 63; i3 += 1) {
            temp2 += (B[i3][64 * ii1] * A[64 * ii0][i3]);
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
