/**
 * This version is stamped on Apr. 14, 2015
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* deriche.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>
  
/* Include benchmark-specific header. */
#include "deriche.h"

#include <math.h>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))


/* Array initialization. */
static
void init_array (int w, int h, DATA_TYPE* alpha,
		 DATA_TYPE POLYBENCH_2D(imgIn,W,H,w,h),
		 DATA_TYPE POLYBENCH_2D(imgOut,W,H,w,h))
{
  int i, j;

  *alpha=0.25; //parameter of the filter

  //input should be between 0 and 1 (grayscale image pixel)
  for (i = 0; i < w; i++)
     for (j = 0; j < h; j++)
	imgIn[i][j] = (DATA_TYPE) ((313*i+991*j)%65536) / 65535.0f;
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int w, int h,
		 DATA_TYPE POLYBENCH_2D(imgOut,W,H,w,h))

{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("imgOut");
  for (i = 0; i < w; i++)
    for (j = 0; j < h; j++) {
      if ((i * h + j) % 20 == 0) fprintf(POLYBENCH_DUMP_TARGET, "\n");
      fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, imgOut[i][j]);
    }
  POLYBENCH_DUMP_END("imgOut");
  POLYBENCH_DUMP_FINISH;
}



/* Main computational kernel. The whole function will be timed,
   including the call and return. */
/* Original code provided by Gael Deest */
static 
void kernel_deriche(int w, int h, DATA_TYPE alpha,
       DATA_TYPE POLYBENCH_2D(imgIn, W, H, w, h),
       DATA_TYPE POLYBENCH_2D(imgOut, W, H, w, h),
       DATA_TYPE POLYBENCH_2D(y1, W, H, w, h),
       DATA_TYPE POLYBENCH_2D(y2, W, H, w, h)) {
    int i,j;
    DATA_TYPE xm1, tm1, ym1, ym2;
    DATA_TYPE xp1, xp2;
    DATA_TYPE tp1, tp2;
    DATA_TYPE yp1, yp2;
   
    DATA_TYPE k;
    DATA_TYPE a1, a2, a3, a4, a5, a6, a7, a8;
    DATA_TYPE b1, b2, c1, c2;


   k = (SCALAR_VAL(1.0)-EXP_FUN(-alpha))*(SCALAR_VAL(1.0)-EXP_FUN(-alpha))/(SCALAR_VAL(1.0)+SCALAR_VAL(2.0)*alpha*EXP_FUN(-alpha)-EXP_FUN(SCALAR_VAL(2.0)*alpha));
   a1 = a5 = k;
   a2 = a6 = k*EXP_FUN(-alpha)*(alpha-SCALAR_VAL(1.0));
   a3 = a7 = k*EXP_FUN(-alpha)*(alpha+SCALAR_VAL(1.0));
   a4 = a8 = -k*EXP_FUN(SCALAR_VAL(-2.0)*alpha);
   b1 =  POW_FUN(SCALAR_VAL(2.0),-alpha);
   b2 = -EXP_FUN(SCALAR_VAL(-2.0)*alpha);
   c1 = c2 = 1;

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
  int t1, t2, t3, t4, t5, t6, t7, t8;
 register int lbv, ubv;
/* Start of CLooG code */
if (_PB_H >= 1) {
  for (t2=0;t2<=_PB_W-1;t2++) {
    ym1 = SCALAR_VAL(0.0);;
    ym2 = SCALAR_VAL(0.0);;
    xm1 = SCALAR_VAL(0.0);;
    for (t4=0;t4<=_PB_H-1;t4++) {
      y1[t2][t4] = a1*imgIn[t2][t4] + a2*xm1 + b1*ym1 + b2*ym2;;
      xm1 = imgIn[t2][t4];;
      ym2 = ym1;;
      ym1 = y1[t2][t4];;
    }
  }
}
if (_PB_H <= 0) {
  for (t2=0;t2<=_PB_W-1;t2++) {
    ym1 = SCALAR_VAL(0.0);;
    ym2 = SCALAR_VAL(0.0);;
    xm1 = SCALAR_VAL(0.0);;
  }
}
if (_PB_H >= 1) {
  for (t2=0;t2<=_PB_W-1;t2++) {
    yp1 = SCALAR_VAL(0.0);;
    yp2 = SCALAR_VAL(0.0);;
    xp1 = SCALAR_VAL(0.0);;
    xp2 = SCALAR_VAL(0.0);;
    for (t4=-_PB_H+1;t4<=0;t4++) {
      y2[t2][-t4] = a3*xp1 + a4*xp2 + b1*yp1 + b2*yp2;;
      xp2 = xp1;;
      xp1 = imgIn[t2][-t4];;
      yp2 = yp1;;
      yp1 = y2[t2][-t4];;
    }
  }
}
if (_PB_H <= 0) {
  for (t2=0;t2<=_PB_W-1;t2++) {
    yp1 = SCALAR_VAL(0.0);;
    yp2 = SCALAR_VAL(0.0);;
    xp1 = SCALAR_VAL(0.0);;
    xp2 = SCALAR_VAL(0.0);;
  }
}
if (_PB_H >= 1) {
  for (t2=0;t2<=floord(_PB_W-1,64);t2++) {
    for (t4=0;t4<=floord(_PB_H-1,64);t4++) {
      for (t5=64*t2;t5<=min(_PB_W-1,64*t2+63);t5++) {
        lbv=64*t4;
        ubv=min(_PB_H-1,64*t4+63);
#pragma ivdep
#pragma vector always
        for (t7=lbv;t7<=ubv;t7++) {
          imgOut[t5][t7] = c1 * (y1[t5][t7] + y2[t5][t7]);;
        }
      }
    }
  }
}
if (_PB_W >= 1) {
  for (t2=0;t2<=_PB_H-1;t2++) {
    tm1 = SCALAR_VAL(0.0);;
    ym1 = SCALAR_VAL(0.0);;
    ym2 = SCALAR_VAL(0.0);;
    for (t4=0;t4<=_PB_W-1;t4++) {
      y1[t4][t2] = a5*imgOut[t4][t2] + a6*tm1 + b1*ym1 + b2*ym2;;
      tm1 = imgOut[t4][t2];;
      ym2 = ym1;;
      ym1 = y1 [t4][t2];;
    }
  }
}
if (_PB_W <= 0) {
  for (t2=0;t2<=_PB_H-1;t2++) {
    tm1 = SCALAR_VAL(0.0);;
    ym1 = SCALAR_VAL(0.0);;
    ym2 = SCALAR_VAL(0.0);;
  }
}
if (_PB_W >= 1) {
  for (t2=0;t2<=_PB_H-1;t2++) {
    tp1 = SCALAR_VAL(0.0);;
    tp2 = SCALAR_VAL(0.0);;
    yp1 = SCALAR_VAL(0.0);;
    yp2 = SCALAR_VAL(0.0);;
    for (t4=-_PB_W+1;t4<=0;t4++) {
      y2[-t4][t2] = a7*tp1 + a8*tp2 + b1*yp1 + b2*yp2;;
      tp2 = tp1;;
      tp1 = imgOut[-t4][t2];;
      yp2 = yp1;;
      yp1 = y2[-t4][t2];;
    }
  }
}
if (_PB_W <= 0) {
  for (t2=0;t2<=_PB_H-1;t2++) {
    tp1 = SCALAR_VAL(0.0);;
    tp2 = SCALAR_VAL(0.0);;
    yp1 = SCALAR_VAL(0.0);;
    yp2 = SCALAR_VAL(0.0);;
  }
}
if (_PB_H >= 1) {
  for (t2=0;t2<=floord(_PB_W-1,64);t2++) {
    for (t4=0;t4<=floord(_PB_H-1,64);t4++) {
      for (t5=64*t2;t5<=min(_PB_W-1,64*t2+63);t5++) {
        lbv=64*t4;
        ubv=min(_PB_H-1,64*t4+63);
#pragma ivdep
#pragma vector always
        for (t7=lbv;t7<=ubv;t7++) {
          imgOut[t5][t7] = c2*(y1[t5][t7] + y2[t5][t7]);;
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
  int w = W;
  int h = H;

  /* Variable declaration/allocation. */
  DATA_TYPE alpha;
  POLYBENCH_2D_ARRAY_DECL(imgIn, DATA_TYPE, W, H, w, h);
  POLYBENCH_2D_ARRAY_DECL(imgOut, DATA_TYPE, W, H, w, h);
  POLYBENCH_2D_ARRAY_DECL(y1, DATA_TYPE, W, H, w, h);
  POLYBENCH_2D_ARRAY_DECL(y2, DATA_TYPE, W, H, w, h);


  /* Initialize array(s). */
  init_array (w, h, &alpha, POLYBENCH_ARRAY(imgIn), POLYBENCH_ARRAY(imgOut));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_deriche (w, h, alpha, POLYBENCH_ARRAY(imgIn), POLYBENCH_ARRAY(imgOut), POLYBENCH_ARRAY(y1), POLYBENCH_ARRAY(y2));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(w, h, POLYBENCH_ARRAY(imgOut)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(imgIn);
  POLYBENCH_FREE_ARRAY(imgOut);
  POLYBENCH_FREE_ARRAY(y1);
  POLYBENCH_FREE_ARRAY(y2);

  return 0;
}
