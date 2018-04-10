/**
 * This version is stamped on Apr. 14, 2015
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* lu.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "lu.h"

#include <math.h>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))


/* Array initialization. */
static
void init_array (int n,
		 DATA_TYPE POLYBENCH_2D(A,N,N,n,n))
{
  int i, j;

  for (i = 0; i < n; i++)
    {
      for (j = 0; j <= i; j++)
	A[i][j] = (DATA_TYPE)(-j % n) / n + 1;
      for (j = i+1; j < n; j++) {
	A[i][j] = 0;
      }
      A[i][i] = 1;
    }

  /* Make the matrix positive semi-definite. */
  /* not necessary for LU, but using same code as cholesky */
  int r,s,t;
  POLYBENCH_2D_ARRAY_DECL(B, DATA_TYPE, N, N, n, n);
  for (r = 0; r < n; ++r)
    for (s = 0; s < n; ++s)
      (POLYBENCH_ARRAY(B))[r][s] = 0;
  for (t = 0; t < n; ++t)
    for (r = 0; r < n; ++r)
      for (s = 0; s < n; ++s)
	(POLYBENCH_ARRAY(B))[r][s] += A[r][t] * A[s][t];
    for (r = 0; r < n; ++r)
      for (s = 0; s < n; ++s)
	A[r][s] = (POLYBENCH_ARRAY(B))[r][s];
  POLYBENCH_FREE_ARRAY(B);

}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_2D(A,N,N,n,n))

{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("A");
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
      if ((i * n + j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
      fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, A[i][j]);
    }
  POLYBENCH_DUMP_END("A");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_lu(int n,
	       DATA_TYPE POLYBENCH_2D(A,N,N,n,n))
{
  int i, j, k;

/* Copyright (C) 1991-2016 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, see
   <http://www.gnu.org/licenses/>.  */
/* This header is separate from features.h so that the compiler can
   include it implicitly at the start of every compilation.  It must
   not itself include <features.h> or any other header that includes
   <features.h> because the implicit include comes before any feature
   test macros that may be defined in a source file before it first
   explicitly includes a system header.  GCC knows the name of this
   header in order to preinclude it.  */
/* glibc's intent is to support the IEC 559 math functionality, real
   and complex.  If the GCC (4.9 and later) predefined macros
   specifying compiler intent are available, use them to determine
   whether the overall intent is to support these features; otherwise,
   presume an older compiler has intent to support these features and
   define these macros by default.  */
/* wchar_t uses Unicode 8.0.0.  Version 8.0 of the Unicode Standard is
   synchronized with ISO/IEC 10646:2014, plus Amendment 1 (published
   2015-05-15).  */
/* We do not support C11 <threads.h>.  */
  int t1, t2, t3, t4, t5, t6;
 register int lbv, ubv;
/* Start of CLooG code */
if (_PB_N >= 2) {
  for (t1=0;t1<=floord(_PB_N-1,64);t1++) {
    for (t2=0;t2<=floord(_PB_N-1,64);t2++) {
      for (t3=0;t3<=min(min(floord(_PB_N-2,64),t1),t2);t3++) {
        if (t2 == t3) {
          for (t4=max(64*t1,64*t2+64);t4<=min(_PB_N-1,64*t1+63);t4++) {
            for (t5=64*t2;t5<=64*t2+62;t5++) {
              A[t4][t5] /= A[t5][t5];;
              for (t6=t5+1;t6<=64*t2+63;t6++) {
                A[t4][t6] -= A[t4][t5] * A[t5][t6];;
              }
            }
            A[t4][(64*t2+63)] /= A[(64*t2+63)][(64*t2+63)];;
          }
        }
        if (t2 >= t3+1) {
          for (t4=max(64*t1,64*t2+64);t4<=min(_PB_N-1,64*t1+63);t4++) {
            for (t5=64*t3;t5<=64*t3+63;t5++) {
              for (t6=64*t2;t6<=64*t2+63;t6++) {
                A[t4][t6] -= A[t4][t5] * A[t5][t6];;
              }
            }
          }
        }
        if ((t1 == t2) && (t1 == t3)) {
          A[(64*t1+1)][64*t1] /= A[64*t1][64*t1];;
          for (t6=64*t1+1;t6<=min(_PB_N-1,64*t1+63);t6++) {
            A[(64*t1+1)][t6] -= A[(64*t1+1)][64*t1] * A[64*t1][t6];;
          }
        }
        for (t4=max(64*t1,64*t3+1);t4<=min(64*t2,64*t1+63);t4++) {
          for (t5=64*t3;t5<=min(64*t3+63,t4-1);t5++) {
            for (t6=64*t2;t6<=min(_PB_N-1,64*t2+63);t6++) {
              A[t4][t6] -= A[t4][t5] * A[t5][t6];;
            }
          }
        }
        if ((t1 == t2) && (t1 == t3)) {
          for (t4=64*t1+2;t4<=min(_PB_N-1,64*t1+63);t4++) {
            for (t5=64*t1;t5<=t4-2;t5++) {
              A[t4][t5] /= A[t5][t5];;
              for (t6=t5+1;t6<=t4-1;t6++) {
                A[t4][t6] -= A[t4][t5] * A[t5][t6];;
              }
              for (t6=t4;t6<=min(_PB_N-1,64*t1+63);t6++) {
                A[t4][t6] -= A[t4][t5] * A[t5][t6];;
              }
            }
            A[t4][(t4-1)] /= A[(t4-1)][(t4-1)];;
            for (t6=t4;t6<=min(_PB_N-1,64*t1+63);t6++) {
              A[t4][t6] -= A[t4][(t4-1)] * A[(t4-1)][t6];;
            }
          }
        }
        if ((t1 == t2) && (t1 >= t3+1)) {
          for (t4=64*t1+1;t4<=min(_PB_N-1,64*t1+63);t4++) {
            for (t5=64*t3;t5<=64*t3+63;t5++) {
              for (t6=64*t1;t6<=t4-1;t6++) {
                A[t4][t6] -= A[t4][t5] * A[t5][t6];;
              }
              for (t6=t4;t6<=min(_PB_N-1,64*t1+63);t6++) {
                A[t4][t6] -= A[t4][t5] * A[t5][t6];;
              }
            }
          }
        }
      }
    }
  }
}
/* End of CLooG code */
}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int n = N;

  /* Variable declaration/allocation. */
  POLYBENCH_2D_ARRAY_DECL(A, DATA_TYPE, N, N, n, n);

  /* Initialize array(s). */
  init_array (n, POLYBENCH_ARRAY(A));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_lu (n, POLYBENCH_ARRAY(A));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(A)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(A);

  return 0;
}
