/**
 * This version is stamped on Apr. 14, 2015
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* adi.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "adi.h"

#include <math.h>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

/* Array initialization. */
static
void init_array (int n,
		 DATA_TYPE POLYBENCH_2D(u,N,N,n,n))
{
  int i, j;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      {
	u[i][j] =  (DATA_TYPE)(i + n-j) / n;
      }
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_2D(u,N,N,n,n))

{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("u");
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
      if ((i * n + j) % 20 == 0) fprintf(POLYBENCH_DUMP_TARGET, "\n");
      fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, u[i][j]);
    }
  POLYBENCH_DUMP_END("u");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
/* Based on a Fortran code fragment from Figure 5 of
 * "Automatic Data and Computation Decomposition on Distributed Memory Parallel Computers"
 * by Peizong Lee and Zvi Meir Kedem, TOPLAS, 2002
 */
static
void kernel_adi(int tsteps, int n,
		DATA_TYPE POLYBENCH_2D(u,N,N,n,n),
		DATA_TYPE POLYBENCH_2D(v,N,N,n,n),
		DATA_TYPE POLYBENCH_2D(p,N,N,n,n),
		DATA_TYPE POLYBENCH_2D(q,N,N,n,n))
{
  int t, i, j;
  DATA_TYPE DX, DY, DT;
  DATA_TYPE B1, B2;
  DATA_TYPE mul1, mul2;
  DATA_TYPE a, b, c, d, e, f;


  DX = SCALAR_VAL(1.0)/(DATA_TYPE)_PB_N;
  DY = SCALAR_VAL(1.0)/(DATA_TYPE)_PB_N;
  DT = SCALAR_VAL(1.0)/(DATA_TYPE)_PB_TSTEPS;
  B1 = SCALAR_VAL(2.0);
  B2 = SCALAR_VAL(1.0);
  mul1 = B1 * DT / (DX * DX);
  mul2 = B2 * DT / (DY * DY);

  a = -mul1 /  SCALAR_VAL(2.0);
  b = SCALAR_VAL(1.0)+mul1;
  c = a;
  d = -mul2 / SCALAR_VAL(2.0);
  e = SCALAR_VAL(1.0)+mul2;
  f = d;

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
  int t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
 register int lbv, ubv;
/* Start of CLooG code */
if ((_PB_N >= 3) && (_PB_TSTEPS >= 1)) {
  for (t2=1;t2<=_PB_TSTEPS;t2++) {
    for (t4=0;t4<=floord(_PB_N-2,32);t4++) {
      lbv=max(1,32*t4);
      ubv=min(_PB_N-2,32*t4+31);
#pragma ivdep
#pragma vector always
      for (t7=lbv;t7<=ubv;t7++) {
        v[0][t7] = SCALAR_VAL(1.0);;
      }
      lbv=max(1,32*t4);
      ubv=min(_PB_N-2,32*t4+31);
#pragma ivdep
#pragma vector always
      for (t7=lbv;t7<=ubv;t7++) {
        p[t7][0] = SCALAR_VAL(0.0);;
      }
      lbv=max(1,32*t4);
      ubv=min(_PB_N-2,32*t4+31);
#pragma ivdep
#pragma vector always
      for (t7=lbv;t7<=ubv;t7++) {
        q[t7][0] = v[0][t7];;
      }
      for (t6=0;t6<=floord(_PB_N-2,32);t6++) {
        for (t7=max(1,32*t4);t7<=min(_PB_N-2,32*t4+31);t7++) {
          for (t9=max(1,32*t6);t9<=min(_PB_N-2,32*t6+31);t9++) {
            p[t7][t9] = -c / (a*p[t7][t9-1]+b);;
            q[t7][t9] = (-d*u[t9][t7-1]+(SCALAR_VAL(1.0)+SCALAR_VAL(2.0)*d)*u[t9][t7] - f*u[t9][t7+1]-a*q[t7][t9-1])/(a*p[t7][t9-1]+b);;
          }
        }
      }
      lbv=max(1,32*t4);
      ubv=min(_PB_N-2,32*t4+31);
#pragma ivdep
#pragma vector always
      for (t7=lbv;t7<=ubv;t7++) {
        v[_PB_N-1][t7] = SCALAR_VAL(1.0);;
      }
      for (t6=ceild(-_PB_N-29,32);t6<=-1;t6++) {
        for (t7=max(1,32*t4);t7<=min(_PB_N-2,32*t4+31);t7++) {
          for (t9=max(32*t6,-_PB_N+2);t9<=32*t6+31;t9++) {
            v[-t9][t7] = p[t7][-t9] * v[-t9+1][t7] + q[t7][-t9];;
          }
        }
      }
    }
    for (t4=0;t4<=floord(_PB_N-2,32);t4++) {
      lbv=max(1,32*t4);
      ubv=min(_PB_N-2,32*t4+31);
#pragma ivdep
#pragma vector always
      for (t7=lbv;t7<=ubv;t7++) {
        u[t7][0] = SCALAR_VAL(1.0);;
      }
      lbv=max(1,32*t4);
      ubv=min(_PB_N-2,32*t4+31);
#pragma ivdep
#pragma vector always
      for (t7=lbv;t7<=ubv;t7++) {
        p[t7][0] = SCALAR_VAL(0.0);;
      }
      lbv=max(1,32*t4);
      ubv=min(_PB_N-2,32*t4+31);
#pragma ivdep
#pragma vector always
      for (t7=lbv;t7<=ubv;t7++) {
        q[t7][0] = u[t7][0];;
      }
      for (t6=0;t6<=floord(_PB_N-2,32);t6++) {
        for (t7=max(1,32*t4);t7<=min(_PB_N-2,32*t4+31);t7++) {
          for (t9=max(1,32*t6);t9<=min(_PB_N-2,32*t6+31);t9++) {
            p[t7][t9] = -f / (d*p[t7][t9-1]+e);;
            q[t7][t9] = (-a*v[t7-1][t9]+(SCALAR_VAL(1.0)+SCALAR_VAL(2.0)*a)*v[t7][t9] - c*v[t7+1][t9]-d*q[t7][t9-1])/(d*p[t7][t9-1]+e);;
          }
        }
      }
      lbv=max(1,32*t4);
      ubv=min(_PB_N-2,32*t4+31);
#pragma ivdep
#pragma vector always
      for (t7=lbv;t7<=ubv;t7++) {
        u[t7][_PB_N-1] = SCALAR_VAL(1.0);;
      }
      for (t6=ceild(-_PB_N-29,32);t6<=-1;t6++) {
        for (t7=max(1,32*t4);t7<=min(_PB_N-2,32*t4+31);t7++) {
          for (t9=max(32*t6,-_PB_N+2);t9<=32*t6+31;t9++) {
            u[t7][-t9] = p[t7][-t9] * u[t7][-t9+1] + q[t7][-t9];;
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
  int tsteps = TSTEPS;

  /* Variable declaration/allocation. */
  POLYBENCH_2D_ARRAY_DECL(u, DATA_TYPE, N, N, n, n);
  POLYBENCH_2D_ARRAY_DECL(v, DATA_TYPE, N, N, n, n);
  POLYBENCH_2D_ARRAY_DECL(p, DATA_TYPE, N, N, n, n);
  POLYBENCH_2D_ARRAY_DECL(q, DATA_TYPE, N, N, n, n);


  /* Initialize array(s). */
  init_array (n, POLYBENCH_ARRAY(u));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_adi (tsteps, n, POLYBENCH_ARRAY(u), POLYBENCH_ARRAY(v), POLYBENCH_ARRAY(p), POLYBENCH_ARRAY(q));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(u)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(u);
  POLYBENCH_FREE_ARRAY(v);
  POLYBENCH_FREE_ARRAY(p);
  POLYBENCH_FREE_ARRAY(q);

  return 0;
}
