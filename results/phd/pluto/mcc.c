/**
 * This version is stamped on Oct 25, 2025
 */

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
#include <stdlib.h>

#define POLYBENCH_PADDING_FACTOR 2

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "mcc.h"

static
int paired(char a, char b)
{
  return ((a == 'A' && b == 'U') || (a == 'U' && b == 'A') ||
          (a == 'G' && b == 'C') || (a == 'C' && b == 'G') ||
          (a == 'G' && b == 'U') || (a == 'U' && b == 'G'));
}

/* Array initialization. */
static
void init_array (int n,
        char POLYBENCH_1D(RNA,N,n),
        DATA_TYPE POLYBENCH_2D(Q,N,N,n,n),
        DATA_TYPE POLYBENCH_2D(Qbp,N,N,n,n))
{
  int i, j;

  srand(42);

  for (i = 0; i < _PB_N; i++) {
    RNA[i] = "ACGU"[rand() % 4];
  }

  for (i = 0; i < _PB_N + POLYBENCH_PADDING_FACTOR; i++) {
    for (j = 0; j < _PB_N + POLYBENCH_PADDING_FACTOR; j++) {
      Q[i][j] = SCALAR_VAL(1.0);
      Qbp[i][j] = SCALAR_VAL(0.0);
    }
  }
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
        DATA_TYPE POLYBENCH_2D(Q,N,N,n,n),
        DATA_TYPE POLYBENCH_2D(Qbp,N,N,n,n))
{
  int i, j;
  int t = 0;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("Q");
  for (i = 1; i <= _PB_N; i++) {
    for (j = 0; j <= _PB_N; j++) {
      if (t % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
      fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, Q[i][j]);
      t++;
    }
  }
  POLYBENCH_DUMP_END("Q");
  t = 0;
  POLYBENCH_DUMP_BEGIN("Qbp");
  for (i = 1; i <= _PB_N; i++) {
    for (j = 0; j <= _PB_N; j++) {
      if (t % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
      fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, Qbp[i][j]);
      t++;
    }
  }
  POLYBENCH_DUMP_END("Qbp");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_mcc(int n,
        int l,
        DATA_TYPE ERT,
        char POLYBENCH_1D(RNA,N,n),
        DATA_TYPE POLYBENCH_2D(Q,N,N,n,n),
        DATA_TYPE POLYBENCH_2D(Qbp,N,N,n,n))
{
  int i, j, k;

/* Copyright (C) 1991-2022 Free Software Foundation, Inc.
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
   <https://www.gnu.org/licenses/>.  */
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
/* wchar_t uses Unicode 10.0.0.  Version 10.0 of the Unicode Standard is
   synchronized with ISO/IEC 10646:2017, fifth edition, plus
   the following additions from Amendment 1 to the fifth edition:
   - 56 emoji characters
   - 285 hentaigana
   - 3 additional Zanabazar Square characters */
  int t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
if ((_PB_N >= 2) && (l >= 0) && (l <= 5)) {
  for (t2=-1;t2<=floord(_PB_N-16,16);t2++) {
    lbp=t2+1;
    ubp=min(floord(_PB_N,16),floord(16*t2+_PB_N+14,16));
#pragma omp parallel for private(lbv,ubv,t5,t6,t7,t8,t9,t10)
    for (t4=lbp;t4<=ubp;t4++) {
      if (t2 == -1) {
        for (t5=max(-_PB_N+1,-16*t4-14);t5<=-16*t4+l-15;t5++) {
          for (t7=-t5+1;t7<=min(_PB_N,16*t4+15);t7++) {
            Q[-t5][t7] = Q[-t5][t7-1];;
          }
        }
      }
      for (t5=max(max(-_PB_N+1,16*t2-16*t4),-16*t4+l-14);t5<=min(-_PB_N+l,16*t2-16*t4+15);t5++) {
        for (t7=max(16*t4,-t5+1);t7<=_PB_N;t7++) {
          Q[-t5][t7] = Q[-t5][t7-1];;
        }
      }
      if (l >= 1) {
        for (t5=max(max(16*t2-16*t4,-_PB_N+l+1),-16*t4+l-14);t5<=min(-16*t4+l,16*t2-16*t4+15);t5++) {
          for (t7=max(16*t4,-t5+1);t7<=-t5+l;t7++) {
            Q[-t5][t7] = Q[-t5][t7-1];;
          }
          for (t7=-t5+l+1;t7<=min(_PB_N,16*t4+15);t7++) {
            Q[-t5][t7] = Q[-t5][t7-1];;
            for (t9=-t5;t9<=t7-l-1;t9++) {
              Qbp[-t5][t7] = paired(RNA[-t5-1], RNA[t7-1]) * Q[-t5+1][t7-1] * ERT;;
              Q[-t5][t7] += Q[-t5][t9-1] * Qbp[t9][t7];;
            }
          }
        }
      }
      if (l >= 1) {
        for (t5=max(16*t2-16*t4,-16*t4+l+1);t5<=16*t2-16*t4+15;t5++) {
          for (t7=16*t4;t7<=min(_PB_N,16*t4+15);t7++) {
            Q[-t5][t7] = Q[-t5][t7-1];;
            for (t9=-t5;t9<=t7-l-1;t9++) {
              Qbp[-t5][t7] = paired(RNA[-t5-1], RNA[t7-1]) * Q[-t5+1][t7-1] * ERT;;
              Q[-t5][t7] += Q[-t5][t9-1] * Qbp[t9][t7];;
            }
          }
        }
      }
      if (l == 0) {
        for (t5=max(max(-_PB_N+1,16*t2-16*t4),-16*t4-14);t5<=16*t2-16*t4+15;t5++) {
          for (t7=max(16*t4,-t5+1);t7<=min(_PB_N,16*t4+15);t7++) {
            Q[-t5][t7] = Q[-t5][t7-1];;
            for (t9=-t5;t9<=t7-1;t9++) {
              Qbp[-t5][t7] = paired(RNA[-t5-1], RNA[t7-1]) * Q[-t5+1][t7-1] * ERT;;
              Q[-t5][t7] += Q[-t5][t9-1] * Qbp[t9][t7];;
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

  int l         = 1;
  DATA_TYPE Ebp = SCALAR_VAL(-1.0);
  DATA_TYPE RT  = SCALAR_VAL(1.0);
  DATA_TYPE ERT = EXP_FUN(-Ebp/RT);

  /* Variable declaration/allocation. */
  POLYBENCH_1D_ARRAY_DECL(RNA, char, N, n);
  POLYBENCH_2D_ARRAY_DECL(Q, DATA_TYPE, N, N, n, n);
  POLYBENCH_2D_ARRAY_DECL(Qbp, DATA_TYPE, N, N, n, n);

  /* Initialize array(s). */
  init_array (n, POLYBENCH_ARRAY(RNA), POLYBENCH_ARRAY(Q), POLYBENCH_ARRAY(Qbp));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_mcc (n, l, ERT, POLYBENCH_ARRAY(RNA), POLYBENCH_ARRAY(Q), POLYBENCH_ARRAY(Qbp));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(Q), POLYBENCH_ARRAY(Qbp)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(RNA);
  POLYBENCH_FREE_ARRAY(Q);
  POLYBENCH_FREE_ARRAY(Qbp);

  return 0;
}
