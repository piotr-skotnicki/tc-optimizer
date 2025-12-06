/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* 2mm.c: this file is part of PolyBench/C */

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
      A[i][j] = (DATA_TYPE) ((i*j+1) % ni) / ni;
  for (i = 0; i < nk; i++)
    for (j = 0; j < nj; j++)
      B[i][j] = (DATA_TYPE) (i*(j+1) % nj) / nj;
  for (i = 0; i < nj; i++)
    for (j = 0; j < nl; j++)
      C[i][j] = (DATA_TYPE) ((i*(j+3)+1) % nl) / nl;
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
  int i, j, k;

  int t1, t2, t3, t4, t5, t6, t7, t8, t9;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
if (_PB_NI >= 1) {
  lbp=0;
  ubp=floord(_PB_NI-1,32);
#pragma omp parallel for private(lbv,ubv,t3,t4,t5,t6,t7,t8,t9)
  for (t2=lbp;t2<=ubp;t2++) {
    if ((_PB_NJ >= 0) && (_PB_NL >= 0)) {
      for (t3=0;t3<=floord(_PB_NJ+_PB_NL-1,32);t3++) {
        if ((_PB_NJ >= _PB_NL+1) && (t3 <= floord(_PB_NL-1,32)) && (t3 >= ceild(_PB_NL-31,32))) {
          for (t4=32*t2;t4<=(min(_PB_NI-1,32*t2+31))-7;t4+=8) {
            lbv=32*t3;
            ubv=_PB_NL-1;
#pragma ivdep
#pragma vector always
            for (t5=lbv;t5<=ubv;t5++) {
              D[t4][t5] *= beta;;
              tmp[t4][t5] = SCALAR_VAL(0.0);;
              D[(t4+1)][t5] *= beta;;
              tmp[(t4+1)][t5] = SCALAR_VAL(0.0);;
              D[(t4+2)][t5] *= beta;;
              tmp[(t4+2)][t5] = SCALAR_VAL(0.0);;
              D[(t4+3)][t5] *= beta;;
              tmp[(t4+3)][t5] = SCALAR_VAL(0.0);;
              D[(t4+4)][t5] *= beta;;
              tmp[(t4+4)][t5] = SCALAR_VAL(0.0);;
              D[(t4+5)][t5] *= beta;;
              tmp[(t4+5)][t5] = SCALAR_VAL(0.0);;
              D[(t4+6)][t5] *= beta;;
              tmp[(t4+6)][t5] = SCALAR_VAL(0.0);;
              D[(t4+7)][t5] *= beta;;
              tmp[(t4+7)][t5] = SCALAR_VAL(0.0);;
            }
            lbv=_PB_NL;
            ubv=min(_PB_NJ-1,32*t3+31);
#pragma ivdep
#pragma vector always
            for (t5=lbv;t5<=ubv;t5++) {
              tmp[t4][t5] = SCALAR_VAL(0.0);;
              tmp[(t4+1)][t5] = SCALAR_VAL(0.0);;
              tmp[(t4+2)][t5] = SCALAR_VAL(0.0);;
              tmp[(t4+3)][t5] = SCALAR_VAL(0.0);;
              tmp[(t4+4)][t5] = SCALAR_VAL(0.0);;
              tmp[(t4+5)][t5] = SCALAR_VAL(0.0);;
              tmp[(t4+6)][t5] = SCALAR_VAL(0.0);;
              tmp[(t4+7)][t5] = SCALAR_VAL(0.0);;
            }
          }
          for (;t4<=min(_PB_NI-1,32*t2+31);t4++) {
            lbv=32*t3;
            ubv=_PB_NL-1;
#pragma ivdep
#pragma vector always
            for (t5=lbv;t5<=ubv;t5++) {
              D[t4][t5] *= beta;;
              tmp[t4][t5] = SCALAR_VAL(0.0);;
            }
            lbv=_PB_NL;
            ubv=min(_PB_NJ-1,32*t3+31);
#pragma ivdep
#pragma vector always
            for (t5=lbv;t5<=ubv;t5++) {
              tmp[t4][t5] = SCALAR_VAL(0.0);;
            }
          }
        }
        if ((_PB_NJ <= _PB_NL-1) && (t3 <= floord(_PB_NJ-1,32)) && (t3 >= ceild(_PB_NJ-31,32))) {
          for (t4=32*t2;t4<=(min(_PB_NI-1,32*t2+31))-7;t4+=8) {
            lbv=32*t3;
            ubv=_PB_NJ-1;
#pragma ivdep
#pragma vector always
            for (t5=lbv;t5<=ubv;t5++) {
              D[t4][t5] *= beta;;
              tmp[t4][t5] = SCALAR_VAL(0.0);;
              D[(t4+1)][t5] *= beta;;
              tmp[(t4+1)][t5] = SCALAR_VAL(0.0);;
              D[(t4+2)][t5] *= beta;;
              tmp[(t4+2)][t5] = SCALAR_VAL(0.0);;
              D[(t4+3)][t5] *= beta;;
              tmp[(t4+3)][t5] = SCALAR_VAL(0.0);;
              D[(t4+4)][t5] *= beta;;
              tmp[(t4+4)][t5] = SCALAR_VAL(0.0);;
              D[(t4+5)][t5] *= beta;;
              tmp[(t4+5)][t5] = SCALAR_VAL(0.0);;
              D[(t4+6)][t5] *= beta;;
              tmp[(t4+6)][t5] = SCALAR_VAL(0.0);;
              D[(t4+7)][t5] *= beta;;
              tmp[(t4+7)][t5] = SCALAR_VAL(0.0);;
            }
            lbv=_PB_NJ;
            ubv=min(_PB_NL-1,32*t3+31);
#pragma ivdep
#pragma vector always
            for (t5=lbv;t5<=ubv;t5++) {
              D[t4][t5] *= beta;;
              D[(t4+1)][t5] *= beta;;
              D[(t4+2)][t5] *= beta;;
              D[(t4+3)][t5] *= beta;;
              D[(t4+4)][t5] *= beta;;
              D[(t4+5)][t5] *= beta;;
              D[(t4+6)][t5] *= beta;;
              D[(t4+7)][t5] *= beta;;
            }
          }
          for (;t4<=min(_PB_NI-1,32*t2+31);t4++) {
            lbv=32*t3;
            ubv=_PB_NJ-1;
#pragma ivdep
#pragma vector always
            for (t5=lbv;t5<=ubv;t5++) {
              D[t4][t5] *= beta;;
              tmp[t4][t5] = SCALAR_VAL(0.0);;
            }
            lbv=_PB_NJ;
            ubv=min(_PB_NL-1,32*t3+31);
#pragma ivdep
#pragma vector always
            for (t5=lbv;t5<=ubv;t5++) {
              D[t4][t5] *= beta;;
            }
          }
        }
        if ((_PB_NJ == _PB_NL) && (t3 <= floord(_PB_NJ-1,32))) {
          for (t4=32*t2;t4<=(min(_PB_NI-1,32*t2+31))-7;t4+=8) {
            lbv=32*t3;
            ubv=min(_PB_NJ-1,32*t3+31);
#pragma ivdep
#pragma vector always
            for (t5=lbv;t5<=ubv;t5++) {
              D[t4][t5] *= beta;;
              tmp[t4][t5] = SCALAR_VAL(0.0);;
              D[(t4+1)][t5] *= beta;;
              tmp[(t4+1)][t5] = SCALAR_VAL(0.0);;
              D[(t4+2)][t5] *= beta;;
              tmp[(t4+2)][t5] = SCALAR_VAL(0.0);;
              D[(t4+3)][t5] *= beta;;
              tmp[(t4+3)][t5] = SCALAR_VAL(0.0);;
              D[(t4+4)][t5] *= beta;;
              tmp[(t4+4)][t5] = SCALAR_VAL(0.0);;
              D[(t4+5)][t5] *= beta;;
              tmp[(t4+5)][t5] = SCALAR_VAL(0.0);;
              D[(t4+6)][t5] *= beta;;
              tmp[(t4+6)][t5] = SCALAR_VAL(0.0);;
              D[(t4+7)][t5] *= beta;;
              tmp[(t4+7)][t5] = SCALAR_VAL(0.0);;
            }
          }
          for (;t4<=min(_PB_NI-1,32*t2+31);t4++) {
            lbv=32*t3;
            ubv=min(_PB_NJ-1,32*t3+31);
#pragma ivdep
#pragma vector always
            for (t5=lbv;t5<=ubv;t5++) {
              D[t4][t5] *= beta;;
              tmp[t4][t5] = SCALAR_VAL(0.0);;
            }
          }
        }
        if ((_PB_NJ <= _PB_NL-1) && (t3 <= floord(_PB_NJ-32,32))) {
          for (t4=32*t2;t4<=(min(_PB_NI-1,32*t2+31))-7;t4+=8) {
            lbv=32*t3;
            ubv=32*t3+31;
#pragma ivdep
#pragma vector always
            for (t5=lbv;t5<=ubv;t5++) {
              D[t4][t5] *= beta;;
              tmp[t4][t5] = SCALAR_VAL(0.0);;
              D[(t4+1)][t5] *= beta;;
              tmp[(t4+1)][t5] = SCALAR_VAL(0.0);;
              D[(t4+2)][t5] *= beta;;
              tmp[(t4+2)][t5] = SCALAR_VAL(0.0);;
              D[(t4+3)][t5] *= beta;;
              tmp[(t4+3)][t5] = SCALAR_VAL(0.0);;
              D[(t4+4)][t5] *= beta;;
              tmp[(t4+4)][t5] = SCALAR_VAL(0.0);;
              D[(t4+5)][t5] *= beta;;
              tmp[(t4+5)][t5] = SCALAR_VAL(0.0);;
              D[(t4+6)][t5] *= beta;;
              tmp[(t4+6)][t5] = SCALAR_VAL(0.0);;
              D[(t4+7)][t5] *= beta;;
              tmp[(t4+7)][t5] = SCALAR_VAL(0.0);;
            }
          }
          for (;t4<=min(_PB_NI-1,32*t2+31);t4++) {
            lbv=32*t3;
            ubv=32*t3+31;
#pragma ivdep
#pragma vector always
            for (t5=lbv;t5<=ubv;t5++) {
              D[t4][t5] *= beta;;
              tmp[t4][t5] = SCALAR_VAL(0.0);;
            }
          }
        }
        if ((_PB_NJ >= _PB_NL+1) && (t3 <= floord(_PB_NL-32,32))) {
          for (t4=32*t2;t4<=(min(_PB_NI-1,32*t2+31))-7;t4+=8) {
            lbv=32*t3;
            ubv=32*t3+31;
#pragma ivdep
#pragma vector always
            for (t5=lbv;t5<=ubv;t5++) {
              D[t4][t5] *= beta;;
              tmp[t4][t5] = SCALAR_VAL(0.0);;
              D[(t4+1)][t5] *= beta;;
              tmp[(t4+1)][t5] = SCALAR_VAL(0.0);;
              D[(t4+2)][t5] *= beta;;
              tmp[(t4+2)][t5] = SCALAR_VAL(0.0);;
              D[(t4+3)][t5] *= beta;;
              tmp[(t4+3)][t5] = SCALAR_VAL(0.0);;
              D[(t4+4)][t5] *= beta;;
              tmp[(t4+4)][t5] = SCALAR_VAL(0.0);;
              D[(t4+5)][t5] *= beta;;
              tmp[(t4+5)][t5] = SCALAR_VAL(0.0);;
              D[(t4+6)][t5] *= beta;;
              tmp[(t4+6)][t5] = SCALAR_VAL(0.0);;
              D[(t4+7)][t5] *= beta;;
              tmp[(t4+7)][t5] = SCALAR_VAL(0.0);;
            }
          }
          for (;t4<=min(_PB_NI-1,32*t2+31);t4++) {
            lbv=32*t3;
            ubv=32*t3+31;
#pragma ivdep
#pragma vector always
            for (t5=lbv;t5<=ubv;t5++) {
              D[t4][t5] *= beta;;
              tmp[t4][t5] = SCALAR_VAL(0.0);;
            }
          }
        }
        if ((t3 <= floord(_PB_NJ-1,32)) && (t3 >= ceild(_PB_NL,32))) {
          for (t4=32*t2;t4<=(min(_PB_NI-1,32*t2+31))-7;t4+=8) {
            lbv=32*t3;
            ubv=min(_PB_NJ-1,32*t3+31);
#pragma ivdep
#pragma vector always
            for (t5=lbv;t5<=ubv;t5++) {
              tmp[t4][t5] = SCALAR_VAL(0.0);;
              tmp[(t4+1)][t5] = SCALAR_VAL(0.0);;
              tmp[(t4+2)][t5] = SCALAR_VAL(0.0);;
              tmp[(t4+3)][t5] = SCALAR_VAL(0.0);;
              tmp[(t4+4)][t5] = SCALAR_VAL(0.0);;
              tmp[(t4+5)][t5] = SCALAR_VAL(0.0);;
              tmp[(t4+6)][t5] = SCALAR_VAL(0.0);;
              tmp[(t4+7)][t5] = SCALAR_VAL(0.0);;
            }
          }
          for (;t4<=min(_PB_NI-1,32*t2+31);t4++) {
            lbv=32*t3;
            ubv=min(_PB_NJ-1,32*t3+31);
#pragma ivdep
#pragma vector always
            for (t5=lbv;t5<=ubv;t5++) {
              tmp[t4][t5] = SCALAR_VAL(0.0);;
            }
          }
        }
        if ((t3 <= floord(_PB_NL-1,32)) && (t3 >= ceild(_PB_NJ,32))) {
          for (t4=32*t2;t4<=(min(_PB_NI-1,32*t2+31))-7;t4+=8) {
            lbv=32*t3;
            ubv=min(_PB_NL-1,32*t3+31);
#pragma ivdep
#pragma vector always
            for (t5=lbv;t5<=ubv;t5++) {
              D[t4][t5] *= beta;;
              D[(t4+1)][t5] *= beta;;
              D[(t4+2)][t5] *= beta;;
              D[(t4+3)][t5] *= beta;;
              D[(t4+4)][t5] *= beta;;
              D[(t4+5)][t5] *= beta;;
              D[(t4+6)][t5] *= beta;;
              D[(t4+7)][t5] *= beta;;
            }
          }
          for (;t4<=min(_PB_NI-1,32*t2+31);t4++) {
            lbv=32*t3;
            ubv=min(_PB_NL-1,32*t3+31);
#pragma ivdep
#pragma vector always
            for (t5=lbv;t5<=ubv;t5++) {
              D[t4][t5] *= beta;;
            }
          }
        }
      }
    }
    if (_PB_NL <= -1) {
      for (t3=0;t3<=floord(_PB_NJ-1,32);t3++) {
        for (t4=32*t2;t4<=(min(_PB_NI-1,32*t2+31))-7;t4+=8) {
          lbv=32*t3;
          ubv=min(_PB_NJ-1,32*t3+31);
#pragma ivdep
#pragma vector always
          for (t5=lbv;t5<=ubv;t5++) {
            tmp[t4][t5] = SCALAR_VAL(0.0);;
            tmp[(t4+1)][t5] = SCALAR_VAL(0.0);;
            tmp[(t4+2)][t5] = SCALAR_VAL(0.0);;
            tmp[(t4+3)][t5] = SCALAR_VAL(0.0);;
            tmp[(t4+4)][t5] = SCALAR_VAL(0.0);;
            tmp[(t4+5)][t5] = SCALAR_VAL(0.0);;
            tmp[(t4+6)][t5] = SCALAR_VAL(0.0);;
            tmp[(t4+7)][t5] = SCALAR_VAL(0.0);;
          }
        }
        for (;t4<=min(_PB_NI-1,32*t2+31);t4++) {
          lbv=32*t3;
          ubv=min(_PB_NJ-1,32*t3+31);
#pragma ivdep
#pragma vector always
          for (t5=lbv;t5<=ubv;t5++) {
            tmp[t4][t5] = SCALAR_VAL(0.0);;
          }
        }
      }
    }
    if (_PB_NJ <= -1) {
      for (t3=0;t3<=floord(_PB_NL-1,32);t3++) {
        for (t4=32*t2;t4<=(min(_PB_NI-1,32*t2+31))-7;t4+=8) {
          lbv=32*t3;
          ubv=min(_PB_NL-1,32*t3+31);
#pragma ivdep
#pragma vector always
          for (t5=lbv;t5<=ubv;t5++) {
            D[t4][t5] *= beta;;
            D[(t4+1)][t5] *= beta;;
            D[(t4+2)][t5] *= beta;;
            D[(t4+3)][t5] *= beta;;
            D[(t4+4)][t5] *= beta;;
            D[(t4+5)][t5] *= beta;;
            D[(t4+6)][t5] *= beta;;
            D[(t4+7)][t5] *= beta;;
          }
        }
        for (;t4<=min(_PB_NI-1,32*t2+31);t4++) {
          lbv=32*t3;
          ubv=min(_PB_NL-1,32*t3+31);
#pragma ivdep
#pragma vector always
          for (t5=lbv;t5<=ubv;t5++) {
            D[t4][t5] *= beta;;
          }
        }
      }
    }
  }
  if (_PB_NJ >= 1) {
    lbp=0;
    ubp=floord(_PB_NI-1,32);
#pragma omp parallel for private(lbv,ubv,t3,t4,t5,t6,t7,t8,t9)
    for (t2=lbp;t2<=ubp;t2++) {
      for (t3=0;t3<=floord(_PB_NJ-1,32);t3++) {
        if (_PB_NK >= 1) {
          for (t5=0;t5<=floord(_PB_NK-1,32);t5++) {
            for (t6=32*t2;t6<=(min(_PB_NI-1,32*t2+31))-7;t6+=8) {
              for (t8=32*t5;t8<=(min(_PB_NK-1,32*t5+31))-7;t8+=8) {
                lbv=32*t3;
                ubv=min(_PB_NJ-1,32*t3+31);
#pragma ivdep
#pragma vector always
                for (t9=lbv;t9<=ubv;t9++) {
                  tmp[t6][t9] += alpha * A[t6][t8] * B[t8][t9];;
                  tmp[(t6+1)][t9] += alpha * A[(t6+1)][t8] * B[t8][t9];;
                  tmp[(t6+2)][t9] += alpha * A[(t6+2)][t8] * B[t8][t9];;
                  tmp[(t6+3)][t9] += alpha * A[(t6+3)][t8] * B[t8][t9];;
                  tmp[(t6+4)][t9] += alpha * A[(t6+4)][t8] * B[t8][t9];;
                  tmp[(t6+5)][t9] += alpha * A[(t6+5)][t8] * B[t8][t9];;
                  tmp[(t6+6)][t9] += alpha * A[(t6+6)][t8] * B[t8][t9];;
                  tmp[(t6+7)][t9] += alpha * A[(t6+7)][t8] * B[t8][t9];;
                  tmp[t6][t9] += alpha * A[t6][(t8+1)] * B[(t8+1)][t9];;
                  tmp[(t6+1)][t9] += alpha * A[(t6+1)][(t8+1)] * B[(t8+1)][t9];;
                  tmp[(t6+2)][t9] += alpha * A[(t6+2)][(t8+1)] * B[(t8+1)][t9];;
                  tmp[(t6+3)][t9] += alpha * A[(t6+3)][(t8+1)] * B[(t8+1)][t9];;
                  tmp[(t6+4)][t9] += alpha * A[(t6+4)][(t8+1)] * B[(t8+1)][t9];;
                  tmp[(t6+5)][t9] += alpha * A[(t6+5)][(t8+1)] * B[(t8+1)][t9];;
                  tmp[(t6+6)][t9] += alpha * A[(t6+6)][(t8+1)] * B[(t8+1)][t9];;
                  tmp[(t6+7)][t9] += alpha * A[(t6+7)][(t8+1)] * B[(t8+1)][t9];;
                  tmp[t6][t9] += alpha * A[t6][(t8+2)] * B[(t8+2)][t9];;
                  tmp[(t6+1)][t9] += alpha * A[(t6+1)][(t8+2)] * B[(t8+2)][t9];;
                  tmp[(t6+2)][t9] += alpha * A[(t6+2)][(t8+2)] * B[(t8+2)][t9];;
                  tmp[(t6+3)][t9] += alpha * A[(t6+3)][(t8+2)] * B[(t8+2)][t9];;
                  tmp[(t6+4)][t9] += alpha * A[(t6+4)][(t8+2)] * B[(t8+2)][t9];;
                  tmp[(t6+5)][t9] += alpha * A[(t6+5)][(t8+2)] * B[(t8+2)][t9];;
                  tmp[(t6+6)][t9] += alpha * A[(t6+6)][(t8+2)] * B[(t8+2)][t9];;
                  tmp[(t6+7)][t9] += alpha * A[(t6+7)][(t8+2)] * B[(t8+2)][t9];;
                  tmp[t6][t9] += alpha * A[t6][(t8+3)] * B[(t8+3)][t9];;
                  tmp[(t6+1)][t9] += alpha * A[(t6+1)][(t8+3)] * B[(t8+3)][t9];;
                  tmp[(t6+2)][t9] += alpha * A[(t6+2)][(t8+3)] * B[(t8+3)][t9];;
                  tmp[(t6+3)][t9] += alpha * A[(t6+3)][(t8+3)] * B[(t8+3)][t9];;
                  tmp[(t6+4)][t9] += alpha * A[(t6+4)][(t8+3)] * B[(t8+3)][t9];;
                  tmp[(t6+5)][t9] += alpha * A[(t6+5)][(t8+3)] * B[(t8+3)][t9];;
                  tmp[(t6+6)][t9] += alpha * A[(t6+6)][(t8+3)] * B[(t8+3)][t9];;
                  tmp[(t6+7)][t9] += alpha * A[(t6+7)][(t8+3)] * B[(t8+3)][t9];;
                  tmp[t6][t9] += alpha * A[t6][(t8+4)] * B[(t8+4)][t9];;
                  tmp[(t6+1)][t9] += alpha * A[(t6+1)][(t8+4)] * B[(t8+4)][t9];;
                  tmp[(t6+2)][t9] += alpha * A[(t6+2)][(t8+4)] * B[(t8+4)][t9];;
                  tmp[(t6+3)][t9] += alpha * A[(t6+3)][(t8+4)] * B[(t8+4)][t9];;
                  tmp[(t6+4)][t9] += alpha * A[(t6+4)][(t8+4)] * B[(t8+4)][t9];;
                  tmp[(t6+5)][t9] += alpha * A[(t6+5)][(t8+4)] * B[(t8+4)][t9];;
                  tmp[(t6+6)][t9] += alpha * A[(t6+6)][(t8+4)] * B[(t8+4)][t9];;
                  tmp[(t6+7)][t9] += alpha * A[(t6+7)][(t8+4)] * B[(t8+4)][t9];;
                  tmp[t6][t9] += alpha * A[t6][(t8+5)] * B[(t8+5)][t9];;
                  tmp[(t6+1)][t9] += alpha * A[(t6+1)][(t8+5)] * B[(t8+5)][t9];;
                  tmp[(t6+2)][t9] += alpha * A[(t6+2)][(t8+5)] * B[(t8+5)][t9];;
                  tmp[(t6+3)][t9] += alpha * A[(t6+3)][(t8+5)] * B[(t8+5)][t9];;
                  tmp[(t6+4)][t9] += alpha * A[(t6+4)][(t8+5)] * B[(t8+5)][t9];;
                  tmp[(t6+5)][t9] += alpha * A[(t6+5)][(t8+5)] * B[(t8+5)][t9];;
                  tmp[(t6+6)][t9] += alpha * A[(t6+6)][(t8+5)] * B[(t8+5)][t9];;
                  tmp[(t6+7)][t9] += alpha * A[(t6+7)][(t8+5)] * B[(t8+5)][t9];;
                  tmp[t6][t9] += alpha * A[t6][(t8+6)] * B[(t8+6)][t9];;
                  tmp[(t6+1)][t9] += alpha * A[(t6+1)][(t8+6)] * B[(t8+6)][t9];;
                  tmp[(t6+2)][t9] += alpha * A[(t6+2)][(t8+6)] * B[(t8+6)][t9];;
                  tmp[(t6+3)][t9] += alpha * A[(t6+3)][(t8+6)] * B[(t8+6)][t9];;
                  tmp[(t6+4)][t9] += alpha * A[(t6+4)][(t8+6)] * B[(t8+6)][t9];;
                  tmp[(t6+5)][t9] += alpha * A[(t6+5)][(t8+6)] * B[(t8+6)][t9];;
                  tmp[(t6+6)][t9] += alpha * A[(t6+6)][(t8+6)] * B[(t8+6)][t9];;
                  tmp[(t6+7)][t9] += alpha * A[(t6+7)][(t8+6)] * B[(t8+6)][t9];;
                  tmp[t6][t9] += alpha * A[t6][(t8+7)] * B[(t8+7)][t9];;
                  tmp[(t6+1)][t9] += alpha * A[(t6+1)][(t8+7)] * B[(t8+7)][t9];;
                  tmp[(t6+2)][t9] += alpha * A[(t6+2)][(t8+7)] * B[(t8+7)][t9];;
                  tmp[(t6+3)][t9] += alpha * A[(t6+3)][(t8+7)] * B[(t8+7)][t9];;
                  tmp[(t6+4)][t9] += alpha * A[(t6+4)][(t8+7)] * B[(t8+7)][t9];;
                  tmp[(t6+5)][t9] += alpha * A[(t6+5)][(t8+7)] * B[(t8+7)][t9];;
                  tmp[(t6+6)][t9] += alpha * A[(t6+6)][(t8+7)] * B[(t8+7)][t9];;
                  tmp[(t6+7)][t9] += alpha * A[(t6+7)][(t8+7)] * B[(t8+7)][t9];;
                }
              }
              for (;t8<=min(_PB_NK-1,32*t5+31);t8++) {
                lbv=32*t3;
                ubv=min(_PB_NJ-1,32*t3+31);
#pragma ivdep
#pragma vector always
                for (t9=lbv;t9<=ubv;t9++) {
                  tmp[t6][t9] += alpha * A[t6][t8] * B[t8][t9];;
                  tmp[(t6+1)][t9] += alpha * A[(t6+1)][t8] * B[t8][t9];;
                  tmp[(t6+2)][t9] += alpha * A[(t6+2)][t8] * B[t8][t9];;
                  tmp[(t6+3)][t9] += alpha * A[(t6+3)][t8] * B[t8][t9];;
                  tmp[(t6+4)][t9] += alpha * A[(t6+4)][t8] * B[t8][t9];;
                  tmp[(t6+5)][t9] += alpha * A[(t6+5)][t8] * B[t8][t9];;
                  tmp[(t6+6)][t9] += alpha * A[(t6+6)][t8] * B[t8][t9];;
                  tmp[(t6+7)][t9] += alpha * A[(t6+7)][t8] * B[t8][t9];;
                }
              }
            }
            for (;t6<=min(_PB_NI-1,32*t2+31);t6++) {
              for (t8=32*t5;t8<=min(_PB_NK-1,32*t5+31);t8++) {
                lbv=32*t3;
                ubv=min(_PB_NJ-1,32*t3+31);
#pragma ivdep
#pragma vector always
                for (t9=lbv;t9<=ubv;t9++) {
                  tmp[t6][t9] += alpha * A[t6][t8] * B[t8][t9];;
                }
              }
            }
          }
        }
        if (_PB_NL >= 1) {
          for (t5=0;t5<=floord(_PB_NL-1,32);t5++) {
            for (t6=32*t2;t6<=(min(_PB_NI-1,32*t2+31))-7;t6+=8) {
              for (t7=32*t3;t7<=(min(_PB_NJ-1,32*t3+31))-7;t7+=8) {
                lbv=32*t5;
                ubv=min(_PB_NL-1,32*t5+31);
#pragma ivdep
#pragma vector always
                for (t9=lbv;t9<=ubv;t9++) {
                  D[t6][t9] += tmp[t6][t7] * C[t7][t9];;
                  D[(t6+1)][t9] += tmp[(t6+1)][t7] * C[t7][t9];;
                  D[(t6+2)][t9] += tmp[(t6+2)][t7] * C[t7][t9];;
                  D[(t6+3)][t9] += tmp[(t6+3)][t7] * C[t7][t9];;
                  D[(t6+4)][t9] += tmp[(t6+4)][t7] * C[t7][t9];;
                  D[(t6+5)][t9] += tmp[(t6+5)][t7] * C[t7][t9];;
                  D[(t6+6)][t9] += tmp[(t6+6)][t7] * C[t7][t9];;
                  D[(t6+7)][t9] += tmp[(t6+7)][t7] * C[t7][t9];;
                  D[t6][t9] += tmp[t6][(t7+1)] * C[(t7+1)][t9];;
                  D[(t6+1)][t9] += tmp[(t6+1)][(t7+1)] * C[(t7+1)][t9];;
                  D[(t6+2)][t9] += tmp[(t6+2)][(t7+1)] * C[(t7+1)][t9];;
                  D[(t6+3)][t9] += tmp[(t6+3)][(t7+1)] * C[(t7+1)][t9];;
                  D[(t6+4)][t9] += tmp[(t6+4)][(t7+1)] * C[(t7+1)][t9];;
                  D[(t6+5)][t9] += tmp[(t6+5)][(t7+1)] * C[(t7+1)][t9];;
                  D[(t6+6)][t9] += tmp[(t6+6)][(t7+1)] * C[(t7+1)][t9];;
                  D[(t6+7)][t9] += tmp[(t6+7)][(t7+1)] * C[(t7+1)][t9];;
                  D[t6][t9] += tmp[t6][(t7+2)] * C[(t7+2)][t9];;
                  D[(t6+1)][t9] += tmp[(t6+1)][(t7+2)] * C[(t7+2)][t9];;
                  D[(t6+2)][t9] += tmp[(t6+2)][(t7+2)] * C[(t7+2)][t9];;
                  D[(t6+3)][t9] += tmp[(t6+3)][(t7+2)] * C[(t7+2)][t9];;
                  D[(t6+4)][t9] += tmp[(t6+4)][(t7+2)] * C[(t7+2)][t9];;
                  D[(t6+5)][t9] += tmp[(t6+5)][(t7+2)] * C[(t7+2)][t9];;
                  D[(t6+6)][t9] += tmp[(t6+6)][(t7+2)] * C[(t7+2)][t9];;
                  D[(t6+7)][t9] += tmp[(t6+7)][(t7+2)] * C[(t7+2)][t9];;
                  D[t6][t9] += tmp[t6][(t7+3)] * C[(t7+3)][t9];;
                  D[(t6+1)][t9] += tmp[(t6+1)][(t7+3)] * C[(t7+3)][t9];;
                  D[(t6+2)][t9] += tmp[(t6+2)][(t7+3)] * C[(t7+3)][t9];;
                  D[(t6+3)][t9] += tmp[(t6+3)][(t7+3)] * C[(t7+3)][t9];;
                  D[(t6+4)][t9] += tmp[(t6+4)][(t7+3)] * C[(t7+3)][t9];;
                  D[(t6+5)][t9] += tmp[(t6+5)][(t7+3)] * C[(t7+3)][t9];;
                  D[(t6+6)][t9] += tmp[(t6+6)][(t7+3)] * C[(t7+3)][t9];;
                  D[(t6+7)][t9] += tmp[(t6+7)][(t7+3)] * C[(t7+3)][t9];;
                  D[t6][t9] += tmp[t6][(t7+4)] * C[(t7+4)][t9];;
                  D[(t6+1)][t9] += tmp[(t6+1)][(t7+4)] * C[(t7+4)][t9];;
                  D[(t6+2)][t9] += tmp[(t6+2)][(t7+4)] * C[(t7+4)][t9];;
                  D[(t6+3)][t9] += tmp[(t6+3)][(t7+4)] * C[(t7+4)][t9];;
                  D[(t6+4)][t9] += tmp[(t6+4)][(t7+4)] * C[(t7+4)][t9];;
                  D[(t6+5)][t9] += tmp[(t6+5)][(t7+4)] * C[(t7+4)][t9];;
                  D[(t6+6)][t9] += tmp[(t6+6)][(t7+4)] * C[(t7+4)][t9];;
                  D[(t6+7)][t9] += tmp[(t6+7)][(t7+4)] * C[(t7+4)][t9];;
                  D[t6][t9] += tmp[t6][(t7+5)] * C[(t7+5)][t9];;
                  D[(t6+1)][t9] += tmp[(t6+1)][(t7+5)] * C[(t7+5)][t9];;
                  D[(t6+2)][t9] += tmp[(t6+2)][(t7+5)] * C[(t7+5)][t9];;
                  D[(t6+3)][t9] += tmp[(t6+3)][(t7+5)] * C[(t7+5)][t9];;
                  D[(t6+4)][t9] += tmp[(t6+4)][(t7+5)] * C[(t7+5)][t9];;
                  D[(t6+5)][t9] += tmp[(t6+5)][(t7+5)] * C[(t7+5)][t9];;
                  D[(t6+6)][t9] += tmp[(t6+6)][(t7+5)] * C[(t7+5)][t9];;
                  D[(t6+7)][t9] += tmp[(t6+7)][(t7+5)] * C[(t7+5)][t9];;
                  D[t6][t9] += tmp[t6][(t7+6)] * C[(t7+6)][t9];;
                  D[(t6+1)][t9] += tmp[(t6+1)][(t7+6)] * C[(t7+6)][t9];;
                  D[(t6+2)][t9] += tmp[(t6+2)][(t7+6)] * C[(t7+6)][t9];;
                  D[(t6+3)][t9] += tmp[(t6+3)][(t7+6)] * C[(t7+6)][t9];;
                  D[(t6+4)][t9] += tmp[(t6+4)][(t7+6)] * C[(t7+6)][t9];;
                  D[(t6+5)][t9] += tmp[(t6+5)][(t7+6)] * C[(t7+6)][t9];;
                  D[(t6+6)][t9] += tmp[(t6+6)][(t7+6)] * C[(t7+6)][t9];;
                  D[(t6+7)][t9] += tmp[(t6+7)][(t7+6)] * C[(t7+6)][t9];;
                  D[t6][t9] += tmp[t6][(t7+7)] * C[(t7+7)][t9];;
                  D[(t6+1)][t9] += tmp[(t6+1)][(t7+7)] * C[(t7+7)][t9];;
                  D[(t6+2)][t9] += tmp[(t6+2)][(t7+7)] * C[(t7+7)][t9];;
                  D[(t6+3)][t9] += tmp[(t6+3)][(t7+7)] * C[(t7+7)][t9];;
                  D[(t6+4)][t9] += tmp[(t6+4)][(t7+7)] * C[(t7+7)][t9];;
                  D[(t6+5)][t9] += tmp[(t6+5)][(t7+7)] * C[(t7+7)][t9];;
                  D[(t6+6)][t9] += tmp[(t6+6)][(t7+7)] * C[(t7+7)][t9];;
                  D[(t6+7)][t9] += tmp[(t6+7)][(t7+7)] * C[(t7+7)][t9];;
                }
              }
              for (;t7<=min(_PB_NJ-1,32*t3+31);t7++) {
                lbv=32*t5;
                ubv=min(_PB_NL-1,32*t5+31);
#pragma ivdep
#pragma vector always
                for (t9=lbv;t9<=ubv;t9++) {
                  D[t6][t9] += tmp[t6][t7] * C[t7][t9];;
                  D[(t6+1)][t9] += tmp[(t6+1)][t7] * C[t7][t9];;
                  D[(t6+2)][t9] += tmp[(t6+2)][t7] * C[t7][t9];;
                  D[(t6+3)][t9] += tmp[(t6+3)][t7] * C[t7][t9];;
                  D[(t6+4)][t9] += tmp[(t6+4)][t7] * C[t7][t9];;
                  D[(t6+5)][t9] += tmp[(t6+5)][t7] * C[t7][t9];;
                  D[(t6+6)][t9] += tmp[(t6+6)][t7] * C[t7][t9];;
                  D[(t6+7)][t9] += tmp[(t6+7)][t7] * C[t7][t9];;
                }
              }
            }
            for (;t6<=min(_PB_NI-1,32*t2+31);t6++) {
              for (t7=32*t3;t7<=min(_PB_NJ-1,32*t3+31);t7++) {
                lbv=32*t5;
                ubv=min(_PB_NL-1,32*t5+31);
#pragma ivdep
#pragma vector always
                for (t9=lbv;t9<=ubv;t9++) {
                  D[t6][t9] += tmp[t6][t7] * C[t7][t9];;
                }
              }
            }
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
