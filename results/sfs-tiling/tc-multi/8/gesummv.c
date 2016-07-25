/**
 * This version is stamped on Apr. 14, 2015
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* gesummv.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "gesummv.h"


/* Array initialization. */
static
void init_array(int n,
		DATA_TYPE *alpha,
		DATA_TYPE *beta,
		DATA_TYPE POLYBENCH_2D(A,N,N,n,n),
		DATA_TYPE POLYBENCH_2D(B,N,N,n,n),
		DATA_TYPE POLYBENCH_1D(x,N,n))
{
  int i, j;

  *alpha = 1.5;
  *beta = 1.2;
  for (i = 0; i < n; i++)
    {
      x[i] = (DATA_TYPE)( i % n) / n;
      for (j = 0; j < n; j++) {
	A[i][j] = (DATA_TYPE) (i*j % n) / n;
	B[i][j] = (DATA_TYPE) (i*j % n) / n;
      }
    }
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_1D(y,N,n))

{
  int i;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("y");
  for (i = 0; i < n; i++) {
    if (i % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
    fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, y[i]);
  }
  POLYBENCH_DUMP_END("y");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_gesummv(int n,
		    DATA_TYPE alpha,
		    DATA_TYPE beta,
		    DATA_TYPE POLYBENCH_2D(A,N,N,n,n),
		    DATA_TYPE POLYBENCH_2D(B,N,N,n,n),
		    DATA_TYPE POLYBENCH_1D(tmp,N,n),
		    DATA_TYPE POLYBENCH_1D(x,N,n),
		    DATA_TYPE POLYBENCH_1D(y,N,n))
{

/* ./tc ../examples/polybench/gesummv.scop.c --correction-tiling --sfs-multiple-scheduling --omp-for-codegen --debug -b 8 */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#pragma scop
const int ii0_lb = 0, ii0_ub = floord(_PB_N - 1, 8);
#pragma omp parallel for
for (register int ii0 = ii0_lb; ii0 <= ii0_ub; ii0 += 1) {
  {
    const int i0_lb = 8 * ii0, i0_ub = min(_PB_N - 1, 8 * ii0 + 7);
    for (register int i0 = i0_lb; i0 <= i0_ub; i0 += 1) {
      tmp[i0] = SCALAR_VAL(0.0);
    }
    const int i0_lb2 = 8 * ii0, i0_ub2 = min(_PB_N - 1, 8 * ii0 + 7);
    for (register int i0 = i0_lb2; i0 <= i0_ub2; i0 += 1) {
      y[i0] = SCALAR_VAL(0.0);
    }
    const int ii2_prim_lb = 0, ii2_prim_ub = (_PB_N - 1) / 8;
    for (register int ii2_prim = ii2_prim_lb; ii2_prim <= ii2_prim_ub; ii2_prim += 1) {
      const int i0_lb = 8 * ii0, i0_ub = min(_PB_N - 1, 8 * ii0 + 7);
      for (register int i0 = i0_lb; i0 <= i0_ub; i0 += 1) {
        const int i2_lb = 8 * ii2_prim, i2_ub = min(_PB_N - 1, 8 * ii2_prim + 7);
        for (register int i2 = i2_lb; i2 <= i2_ub; i2 += 1) {
          tmp[i0] = ((A[i0][i2] * x[i2]) + tmp[i0]);
        }
      }
      const int i0_lb2 = 8 * ii0, i0_ub2 = min(_PB_N - 1, 8 * ii0 + 7);
      for (register int i0 = i0_lb2; i0 <= i0_ub2; i0 += 1) {
        const int i2_lb = 8 * ii2_prim, i2_ub = min(_PB_N - 1, 8 * ii2_prim + 7);
        for (register int i2 = i2_lb; i2 <= i2_ub; i2 += 1) {
          y[i0] = ((B[i0][i2] * x[i2]) + y[i0]);
        }
      }
    }
  }
  const int i0_lb = 8 * ii0, i0_ub = min(_PB_N - 1, 8 * ii0 + 7);
  for (register int i0 = i0_lb; i0 <= i0_ub; i0 += 1) {
    y[i0] = ((alpha * tmp[i0]) + (beta * y[i0]));
  }
}
#pragma endscop

}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int n = N;

  /* Variable declaration/allocation. */
  DATA_TYPE alpha;
  DATA_TYPE beta;
  POLYBENCH_2D_ARRAY_DECL(A, DATA_TYPE, N, N, n, n);
  POLYBENCH_2D_ARRAY_DECL(B, DATA_TYPE, N, N, n, n);
  POLYBENCH_1D_ARRAY_DECL(tmp, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(x, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(y, DATA_TYPE, N, n);


  /* Initialize array(s). */
  init_array (n, &alpha, &beta,
	      POLYBENCH_ARRAY(A),
	      POLYBENCH_ARRAY(B),
	      POLYBENCH_ARRAY(x));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_gesummv (n, alpha, beta,
		  POLYBENCH_ARRAY(A),
		  POLYBENCH_ARRAY(B),
		  POLYBENCH_ARRAY(tmp),
		  POLYBENCH_ARRAY(x),
		  POLYBENCH_ARRAY(y));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(y)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(B);
  POLYBENCH_FREE_ARRAY(tmp);
  POLYBENCH_FREE_ARRAY(x);
  POLYBENCH_FREE_ARRAY(y);

  return 0;
}
