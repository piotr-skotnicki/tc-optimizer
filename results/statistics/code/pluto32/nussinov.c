/**
 * This version is stamped on Apr. 14, 2015
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* nussinov.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "nussinov.h"

#include <math.h>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

/* RNA bases represented as chars, range is [0,3] */
typedef char base;

#define match(b1, b2) (((b1)+(b2)) == 3 ? 1 : 0)
#define max_score(s1, s2) ((s1 >= s2) ? s1 : s2)

/* Array initialization. */
static
void init_array (int n,
                 base POLYBENCH_1D(seq,N,n),
		 DATA_TYPE POLYBENCH_2D(table,N,N,n,n))
{
  int i, j;

  //base is AGCT/0..3
  for (i=0; i <n; i++) {
     seq[i] = (base)((i+1)%4);
  }

  for (i=0; i <n; i++)
     for (j=0; j <n; j++)
       table[i][j] = 0;
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_2D(table,N,N,n,n))

{
  int i, j;
  int t = 0;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("table");
  for (i = 0; i < n; i++) {
    for (j = i; j < n; j++) {
      if (t % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
      fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, table[i][j]);
      t++;
    }
  }
  POLYBENCH_DUMP_END("table");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
/*
  Original version by Dave Wonnacott at Haverford College <davew@cs.haverford.edu>, 
  with help from Allison Lake, Ting Zhou, and Tian Jin, 
  based on algorithm by Nussinov, described in Allison Lake's senior thesis.
*/
static
void kernel_nussinov(int n, base POLYBENCH_1D(seq,N,n),
			   DATA_TYPE POLYBENCH_2D(table,N,N,n,n))
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
  int t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
 register int lbv, ubv;
/* Start of CLooG code */
if (_PB_N >= 2) {
  for (t2=ceild(-_PB_N-29,32);t2<=0;t2++) {
    for (t4=max(0,-t2-1);t4<=floord(_PB_N-1,32);t4++) {
      if ((t2 <= floord(-_PB_N+2,32)) && (t4 >= ceild(_PB_N-31,32))) {
        table[(_PB_N-2)][(_PB_N-1)] = max_score(table[(_PB_N-2)][(_PB_N-1)], table[(_PB_N-2)][(_PB_N-1)-1]);;
        table[(_PB_N-2)][(_PB_N-1)] = max_score(table[(_PB_N-2)][(_PB_N-1)], table[(_PB_N-2)+1][(_PB_N-1)]);;
        table[(_PB_N-2)][(_PB_N-1)] = max_score(table[(_PB_N-2)][(_PB_N-1)], table[(_PB_N-2)+1][(_PB_N-1)-1]);;
      }
      if ((t2 == -t4-1) && (t2 >= ceild(-_PB_N,32))) {
        table[(-32*t2-2)][(-32*t2-1)] = max_score(table[(-32*t2-2)][(-32*t2-1)], table[(-32*t2-2)][(-32*t2-1)-1]);;
        table[(-32*t2-2)][(-32*t2-1)] = max_score(table[(-32*t2-2)][(-32*t2-1)], table[(-32*t2-2)+1][(-32*t2-1)]);;
        table[(-32*t2-2)][(-32*t2-1)] = max_score(table[(-32*t2-2)][(-32*t2-1)], table[(-32*t2-2)+1][(-32*t2-1)-1]);;
      }
      for (t5=max(max(32*t2,-_PB_N+3),-32*t4-29);t5<=min(min(0,32*t2+31),-32*t4+1);t5++) {
        table[-t5][(-t5+1)] = max_score(table[-t5][(-t5+1)], table[-t5][(-t5+1)-1]);;
        table[-t5][(-t5+1)] = max_score(table[-t5][(-t5+1)], table[-t5+1][(-t5+1)]);;
        table[-t5][(-t5+1)] = max_score(table[-t5][(-t5+1)], table[-t5+1][(-t5+1)-1]);;
        for (t7=-t5+2;t7<=min(_PB_N-1,32*t4+31);t7++) {
          table[-t5][t7] = max_score(table[-t5][t7], table[-t5][t7-1]);;
          table[-t5][t7] = max_score(table[-t5][t7], table[-t5+1][t7]);;
          table[-t5][t7] = max_score(table[-t5][t7], table[-t5+1][t7-1]+match(seq[-t5], seq[t7]));;
          for (t9=-t5+1;t9<=t7-1;t9++) {
            table[-t5][t7] = max_score(table[-t5][t7], table[-t5][t9] + table[t9+1][t7]);;
          }
        }
      }
      for (t5=max(32*t2,-32*t4+2);t5<=min(0,32*t2+31);t5++) {
        for (t7=32*t4;t7<=min(_PB_N-1,32*t4+31);t7++) {
          table[-t5][t7] = max_score(table[-t5][t7], table[-t5][t7-1]);;
          table[-t5][t7] = max_score(table[-t5][t7], table[-t5+1][t7]);;
          table[-t5][t7] = max_score(table[-t5][t7], table[-t5+1][t7-1]+match(seq[-t5], seq[t7]));;
          for (t9=-t5+1;t9<=t7-1;t9++) {
            table[-t5][t7] = max_score(table[-t5][t7], table[-t5][t9] + table[t9+1][t7]);;
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
  POLYBENCH_1D_ARRAY_DECL(seq, base, N, n);
  POLYBENCH_2D_ARRAY_DECL(table, DATA_TYPE, N, N, n, n);

  /* Initialize array(s). */
  init_array (n, POLYBENCH_ARRAY(seq), POLYBENCH_ARRAY(table));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_nussinov (n, POLYBENCH_ARRAY(seq), POLYBENCH_ARRAY(table));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(table)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(seq);
  POLYBENCH_FREE_ARRAY(table);

  return 0;
}
