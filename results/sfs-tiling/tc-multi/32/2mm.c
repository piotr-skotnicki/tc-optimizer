/**
 * This version is stamped on Apr. 14, 2015
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* 2mm.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "2mm.h"


/* Array initialization. */
static
void init_array(int ni, int nj, int nk, int nl,
		DATA_TYPE *alpha,
		DATA_TYPE *beta,
		DATA_TYPE POLYBENCH_2D(A,NI,NK,ni,nk),
		DATA_TYPE POLYBENCH_2D(B,NK,NJ,nk,nj),
		DATA_TYPE POLYBENCH_2D(C,NJ,NL,nj,nl),
		DATA_TYPE POLYBENCH_2D(D,NI,NL,ni,nl))
{
  int i, j;

  *alpha = 1.5;
  *beta = 1.2;
  for (i = 0; i < ni; i++)
    for (j = 0; j < nk; j++)
      A[i][j] = (DATA_TYPE) (i*j % ni) / ni;
  for (i = 0; i < nk; i++)
    for (j = 0; j < nj; j++)
      B[i][j] = (DATA_TYPE) (i*(j+1) % nj) / nj;
  for (i = 0; i < nj; i++)
    for (j = 0; j < nl; j++)
      C[i][j] = (DATA_TYPE) (i*(j+3) % nl) / nl;
  for (i = 0; i < ni; i++)
    for (j = 0; j < nl; j++)
      D[i][j] = (DATA_TYPE) (i*(j+2) % nk) / nk;
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int ni, int nl,
		 DATA_TYPE POLYBENCH_2D(D,NI,NL,ni,nl))
{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("D");
  for (i = 0; i < ni; i++)
    for (j = 0; j < nl; j++) {
	if ((i * ni + j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
	fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, D[i][j]);
    }
  POLYBENCH_DUMP_END("D");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_2mm(int ni, int nj, int nk, int nl,
		DATA_TYPE alpha,
		DATA_TYPE beta,
		DATA_TYPE POLYBENCH_2D(tmp,NI,NJ,ni,nj),
		DATA_TYPE POLYBENCH_2D(A,NI,NK,ni,nk),
		DATA_TYPE POLYBENCH_2D(B,NK,NJ,nk,nj),
		DATA_TYPE POLYBENCH_2D(C,NJ,NL,nj,nl),
		DATA_TYPE POLYBENCH_2D(D,NI,NL,ni,nl))
{

/* ./tc ../examples/polybench/2mm.scop.c --correction-tiling --sfs-multiple-scheduling --omp-for-codegen --debug -b 32 */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#pragma scop
if (_PB_NJ >= 1) {
  const int ii1_lb = 0, ii1_ub = floord(_PB_NI - 1, 32);
  #pragma omp parallel for
  for (register int ii1 = ii1_lb; ii1 <= ii1_ub; ii1 += 1) {
    if (_PB_NL >= 1) {
      const int ii2_prim_lb = 0, ii2_prim_ub = floord(_PB_NJ - 1, 32);
      for (register int ii2_prim = ii2_prim_lb; ii2_prim <= ii2_prim_ub; ii2_prim += 1) {
        const int i1_lb = 32 * ii1, i1_ub = min(_PB_NI - 1, 32 * ii1 + 31);
        for (register int i1 = i1_lb; i1 <= i1_ub; i1 += 1) {
          const int i2_lb = 32 * ii2_prim, i2_ub = min(_PB_NJ - 1, 32 * ii2_prim + 31);
          for (register int i2 = i2_lb; i2 <= i2_ub; i2 += 1) {
            tmp[i1][i2] = SCALAR_VAL(0.0);
          }
        }
        const int ii4_prim_lb = 0, ii4_prim_ub = floord(_PB_NK - 1, 32);
        for (register int ii4_prim = ii4_prim_lb; ii4_prim <= ii4_prim_ub; ii4_prim += 1) {
          const int i1_lb = 32 * ii1, i1_ub = min(_PB_NI - 1, 32 * ii1 + 31);
          for (register int i1 = i1_lb; i1 <= i1_ub; i1 += 1) {
            const int i2_lb = 32 * ii2_prim, i2_ub = min(_PB_NJ - 1, 32 * ii2_prim + 31);
            for (register int i2 = i2_lb; i2 <= i2_ub; i2 += 1) {
              const int i4_lb = 32 * ii4_prim, i4_ub = min(_PB_NK - 1, 32 * ii4_prim + 31);
              for (register int i4 = i4_lb; i4 <= i4_ub; i4 += 1) {
                tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
              }
            }
          }
        }
      }
    } else {
      const int ii2_lb = 0, ii2_ub = floord(_PB_NJ - 1, 32);
      for (register int ii2 = ii2_lb; ii2 <= ii2_ub; ii2 += 1) {
        const int i1_lb = 32 * ii1, i1_ub = min(_PB_NI - 1, 32 * ii1 + 31);
        for (register int i1 = i1_lb; i1 <= i1_ub; i1 += 1) {
          const int i2_lb = 32 * ii2, i2_ub = min(_PB_NJ - 1, 32 * ii2 + 31);
          for (register int i2 = i2_lb; i2 <= i2_ub; i2 += 1) {
            tmp[i1][i2] = SCALAR_VAL(0.0);
          }
        }
        const int ii4_prim_lb = 0, ii4_prim_ub = floord(_PB_NK - 1, 32);
        for (register int ii4_prim = ii4_prim_lb; ii4_prim <= ii4_prim_ub; ii4_prim += 1) {
          const int i1_lb = 32 * ii1, i1_ub = min(_PB_NI - 1, 32 * ii1 + 31);
          for (register int i1 = i1_lb; i1 <= i1_ub; i1 += 1) {
            const int i2_lb = 32 * ii2, i2_ub = min(_PB_NJ - 1, 32 * ii2 + 31);
            for (register int i2 = i2_lb; i2 <= i2_ub; i2 += 1) {
              const int i4_lb = 32 * ii4_prim, i4_ub = min(_PB_NK - 1, 32 * ii4_prim + 31);
              for (register int i4 = i4_lb; i4 <= i4_ub; i4 += 1) {
                tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
              }
            }
          }
        }
      }
    }
    const int ii2_prim_lb = 0, ii2_prim_ub = floord(_PB_NL - 1, 32);
    for (register int ii2_prim = ii2_prim_lb; ii2_prim <= ii2_prim_ub; ii2_prim += 1) {
      const int i1_lb = 32 * ii1, i1_ub = min(_PB_NI - 1, 32 * ii1 + 31);
      for (register int i1 = i1_lb; i1 <= i1_ub; i1 += 1) {
        const int i2_lb = 32 * ii2_prim, i2_ub = min(_PB_NL - 1, 32 * ii2_prim + 31);
        for (register int i2 = i2_lb; i2 <= i2_ub; i2 += 1) {
          D[i1][i2] *= beta;
        }
      }
      const int ii4_prim_lb = 0, ii4_prim_ub = floord(_PB_NJ - 1, 32);
      for (register int ii4_prim = ii4_prim_lb; ii4_prim <= ii4_prim_ub; ii4_prim += 1) {
        const int i1_lb = 32 * ii1, i1_ub = min(_PB_NI - 1, 32 * ii1 + 31);
        for (register int i1 = i1_lb; i1 <= i1_ub; i1 += 1) {
          const int i2_lb = 32 * ii2_prim, i2_ub = min(_PB_NL - 1, 32 * ii2_prim + 31);
          for (register int i2 = i2_lb; i2 <= i2_ub; i2 += 1) {
            const int i4_lb = 32 * ii4_prim, i4_ub = min(_PB_NJ - 1, 32 * ii4_prim + 31);
            for (register int i4 = i4_lb; i4 <= i4_ub; i4 += 1) {
              D[i1][i2] += (tmp[i1][i4] * C[i4][i2]);
            }
          }
        }
      }
    }
  }
} else {
  const int ii1_lb = 0, ii1_ub = floord(_PB_NI - 1, 32);
  #pragma omp parallel for
  for (register int ii1 = ii1_lb; ii1 <= ii1_ub; ii1 += 1) {
    const int ii2_lb = 0, ii2_ub = floord(_PB_NL - 1, 32);
    for (register int ii2 = ii2_lb; ii2 <= ii2_ub; ii2 += 1) {
      const int i1_lb = 32 * ii1, i1_ub = min(_PB_NI - 1, 32 * ii1 + 31);
      for (register int i1 = i1_lb; i1 <= i1_ub; i1 += 1) {
        const int i2_lb = 32 * ii2, i2_ub = min(_PB_NL - 1, 32 * ii2 + 31);
        for (register int i2 = i2_lb; i2 <= i2_ub; i2 += 1) {
          D[i1][i2] *= beta;
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
  int ni = NI;
  int nj = NJ;
  int nk = NK;
  int nl = NL;

  /* Variable declaration/allocation. */
  DATA_TYPE alpha;
  DATA_TYPE beta;
  POLYBENCH_2D_ARRAY_DECL(tmp,DATA_TYPE,NI,NJ,ni,nj);
  POLYBENCH_2D_ARRAY_DECL(A,DATA_TYPE,NI,NK,ni,nk);
  POLYBENCH_2D_ARRAY_DECL(B,DATA_TYPE,NK,NJ,nk,nj);
  POLYBENCH_2D_ARRAY_DECL(C,DATA_TYPE,NJ,NL,nj,nl);
  POLYBENCH_2D_ARRAY_DECL(D,DATA_TYPE,NI,NL,ni,nl);

  /* Initialize array(s). */
  init_array (ni, nj, nk, nl, &alpha, &beta,
	      POLYBENCH_ARRAY(A),
	      POLYBENCH_ARRAY(B),
	      POLYBENCH_ARRAY(C),
	      POLYBENCH_ARRAY(D));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_2mm (ni, nj, nk, nl,
	      alpha, beta,
	      POLYBENCH_ARRAY(tmp),
	      POLYBENCH_ARRAY(A),
	      POLYBENCH_ARRAY(B),
	      POLYBENCH_ARRAY(C),
	      POLYBENCH_ARRAY(D));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(ni, nl,  POLYBENCH_ARRAY(D)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(tmp);
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(B);
  POLYBENCH_FREE_ARRAY(C);
  POLYBENCH_FREE_ARRAY(D);

  return 0;
}
