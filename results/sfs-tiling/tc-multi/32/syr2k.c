/**
 * This version is stamped on Apr. 14, 2015
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* syr2k.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "syr2k.h"


/* Array initialization. */
static
void init_array(int n, int m,
		DATA_TYPE *alpha,
		DATA_TYPE *beta,
		DATA_TYPE POLYBENCH_2D(C,N,N,n,n),
		DATA_TYPE POLYBENCH_2D(A,N,M,n,m),
		DATA_TYPE POLYBENCH_2D(B,N,M,n,m))
{
  int i, j;

  *alpha = 1.5;
  *beta = 1.2;
  for (i = 0; i < n; i++)
    for (j = 0; j < m; j++) {
      A[i][j] = (DATA_TYPE) (i*j%n) / n;
      B[i][j] = (DATA_TYPE) (i*j%m) / m;
    }
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      C[i][j] = (DATA_TYPE) (i*j%n) / m;
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_2D(C,N,N,n,n))
{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("C");
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
	if ((i * n + j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
	fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, C[i][j]);
    }
  POLYBENCH_DUMP_END("C");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_syr2k(int n, int m,
		  DATA_TYPE alpha,
		  DATA_TYPE beta,
		  DATA_TYPE POLYBENCH_2D(C,N,N,n,n),
		  DATA_TYPE POLYBENCH_2D(A,N,M,n,m),
		  DATA_TYPE POLYBENCH_2D(B,N,M,n,m))
{

/* ./tc ../examples/polybench/syr2k.scop.c --correction-tiling --sfs-multiple-scheduling --omp-for-codegen --debug -b 32 */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#pragma scop
const int ii1_lb = 0, ii1_ub = floord(_PB_N - 1, 32);
#pragma omp parallel for
for (register int ii1 = ii1_lb; ii1 <= ii1_ub; ii1 += 1) {
  const int ii2_lb = 0, ii2_ub = (_PB_N - 1) / 32;
  for (register int ii2 = ii2_lb; ii2 <= ii2_ub; ii2 += 1) {
    const int i1_lb = 32 * ii1, i1_ub = min(_PB_N - 1, 32 * ii1 + 31);
    for (register int i1 = i1_lb; i1 <= i1_ub; i1 += 1) {
      const int i2_lb = 32 * ii2, i2_ub = min(_PB_N - 1, 32 * ii2 + 31);
      for (register int i2 = i2_lb; i2 <= i2_ub; i2 += 1) {
        C[i1][i2] *= beta;
      }
    }
    const int ii2_prim_lb = 0, ii2_prim_ub = floord(_PB_M - 1, 32);
    for (register int ii2_prim = ii2_prim_lb; ii2_prim <= ii2_prim_ub; ii2_prim += 1) {
      const int i1_lb = 32 * ii1, i1_ub = min(_PB_N - 1, 32 * ii1 + 31);
      for (register int i1 = i1_lb; i1 <= i1_ub; i1 += 1) {
        const int i2_lb = 32 * ii2_prim, i2_ub = min(_PB_M - 1, 32 * ii2_prim + 31);
        for (register int i2 = i2_lb; i2 <= i2_ub; i2 += 1) {
          const int i3_lb = 32 * ii2, i3_ub = min(_PB_N - 1, 32 * ii2 + 31);
          for (register int i3 = i3_lb; i3 <= i3_ub; i3 += 1) {
            C[i1][i3] += (((A[i3][i2] * alpha) * B[i1][i2]) + ((B[i3][i2] * alpha) * A[i1][i2]));
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
  int n = N;
  int m = M;

  /* Variable declaration/allocation. */
  DATA_TYPE alpha;
  DATA_TYPE beta;
  POLYBENCH_2D_ARRAY_DECL(C,DATA_TYPE,N,N,n,n);
  POLYBENCH_2D_ARRAY_DECL(A,DATA_TYPE,N,M,n,m);
  POLYBENCH_2D_ARRAY_DECL(B,DATA_TYPE,N,M,n,m);

  /* Initialize array(s). */
  init_array (n, m, &alpha, &beta,
	      POLYBENCH_ARRAY(C),
	      POLYBENCH_ARRAY(A),
	      POLYBENCH_ARRAY(B));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_syr2k (n, m,
		alpha, beta,
		POLYBENCH_ARRAY(C),
		POLYBENCH_ARRAY(A),
		POLYBENCH_ARRAY(B));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(C)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(C);
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(B);

  return 0;
}
