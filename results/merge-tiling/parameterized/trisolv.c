/**
 * This version is stamped on Apr. 14, 2015
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* trisolv.c: this file is part of PolyBench/C */

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

/* ./tc ../examples/polybench/trisolv.scop.c --merge-tiling --free-scheduling --omp-for-codegen -b 64 --debug */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#pragma scop
{
  #pragma omp parallel for
  for (int ii0 = 0; ii0 <= floord(_PB_N - 1, 64); ii0 += 1)
    for (int c4 = 64 * ii0; c4 <= min(_PB_N - 1, 64 * ii0 + 63); c4 += 1)
      x[c4] = b[c4];
  for (int k = 1; k <= (_PB_N + 30) / 32; k += 1)
    #pragma omp parallel for
    for (int ii0 = k / 2; ii0 <= min(k - 1, (_PB_N - 1) / 64); ii0 += 1) {
      if (k >= 2) {
        for (int c4 = 64 * ii0; c4 <= min(_PB_N - 1, 64 * ii0 + 63); c4 += 1) {
          for (int c6 = 64 * k - 64 * ii0 - 64; c6 < min(64 * k - 64 * ii0, c4); c6 += 1)
            x[c4] -= (L[c4][c6] * x[c6]);
          if (2 * ii0 + 1 == k)
            x[c4] = (x[c4] / L[c4][c4]);
        }
      } else
        for (int c4 = 0; c4 <= min(63, _PB_N - 1); c4 += 1) {
          for (int c6 = 0; c6 < c4; c6 += 1)
            x[c4] -= (L[c4][c6] * x[c6]);
          x[c4] = (x[c4] / L[c4][c4]);
        }
    }
  if ((_PB_N - 1) % 64 == 0)
    x[_PB_N - 1] = (x[_PB_N - 1] / L[_PB_N - 1][_PB_N - 1]);
}
#pragma endscop


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
