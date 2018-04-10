/**
 * This version is stamped on Apr. 14, 2015
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* heat-3d.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "heat-3d.h"


/* Array initialization. */
static
void init_array (int n,
		 DATA_TYPE POLYBENCH_3D(A,N,N,N,n,n,n),
		 DATA_TYPE POLYBENCH_3D(B,N,N,N,n,n,n))
{
  int i, j, k;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      for (k = 0; k < n; k++)
        A[i][j][k] = B[i][j][k] = (DATA_TYPE) (i + j + (n-k))* 10 / (n);
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_3D(A,N,N,N,n,n,n))

{
  int i, j, k;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("A");
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      for (k = 0; k < n; k++) {
         if ((i * n * n + j * n + k) % 20 == 0) fprintf(POLYBENCH_DUMP_TARGET, "\n");
         fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, A[i][j][k]);
      }
  POLYBENCH_DUMP_END("A");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_heat_3d(int tsteps,
		      int n,
		      DATA_TYPE POLYBENCH_3D(A,N,N,N,n,n,n),
		      DATA_TYPE POLYBENCH_3D(B,N,N,N,n,n,n))
{
  int t, i, j, k;

/* TC Optimizing Compiler 0.2.26 */
/* ./tc ../examples/polybench/heat-3d.scop.c --correction-tiling --lex-scheduling --serial-codegen --debug -b 64 */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#pragma scop
for (int ii0 = 0; ii0 <= floord(_PB_TSTEPS - 1, 64); ii0 += 1) {
  for (int ii1 = 0; ii1 <= min(min(1, _PB_TSTEPS - 64 * ii0 - 1), floord(_PB_N - 3, 64)); ii1 += 1) {
    if (ii1 == 0) {
      for (int ii2 = 0; ii2 <= (_PB_N - 3) / 64; ii2 += 1) {
        for (int ii3 = 0; ii3 <= (_PB_N - 3) / 64; ii3 += 1) {
          for (int ii4 = 0; ii4 <= (_PB_N - 3) / 64; ii4 += 1) {
            for (int i2 = 64 * ii2 + 1; i2 <= min(_PB_N - 2, 64 * ii2 + 64); i2 += 1) {
              for (int i3 = 64 * ii3 + 1; i3 <= min(_PB_N - 2, 64 * ii3 + 64); i3 += 1) {
                for (int i4 = 64 * ii4 + 1; i4 <= min(_PB_N - 2, 64 * ii4 + 64); i4 += 1) {
                  B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                }
              }
            }
          }
        }
      }
    } else if (_PB_TSTEPS >= 64 * ii0 + 3) {
      for (int ii2 = 0; ii2 <= min((_PB_N - 5) / 32 - 1, (_PB_N - 3) / 64); ii2 += 1) {
        if (ii2 == 0) {
          for (int ii4 = 0; ii4 <= (_PB_N - 3) / 64; ii4 += 1) {
            if (_PB_N >= 64 * ii4 + 67) {
              for (int i0 = 64 * ii0 + 1; i0 <= min(_PB_TSTEPS, 64 * ii0 + 32); i0 += 1) {
                if (i0 >= 64 * ii0 + 2) {
                  for (int i2 = 1; i2 <= 128 * ii0 - 2 * i0 + 67; i2 += 1) {
                    if (2 * i0 + i2 >= 128 * ii0 + 7) {
                      for (int i3 = 1; i3 <= 128 * ii0 - 2 * i0 - i2 + 68; i3 += 1) {
                        if (i0 >= 64 * ii0 + 3 && 64 * ii0 + 32 * ii4 + 1 >= i0) {
                          B[i2][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 4] - (SCALAR_VAL(2.0) * A[i2][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 4])) + A[i2 - 1][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][128 * ii0 + 64 * ii4 - 2 * i0 + 4] - (SCALAR_VAL(2.0) * A[i2][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 4])) + A[i2][i3 - 1][128 * ii0 + 64 * ii4 - 2 * i0 + 4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 5] - (SCALAR_VAL(2.0) * A[i2][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 4])) + A[i2][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 3]))) + A[i2][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 4]);
                          if (i2 == 1) {
                            B[1][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 5] = ((((SCALAR_VAL(0.125) * ((A[2][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 5] - (SCALAR_VAL(2.0) * A[1][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 5])) + A[0][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 5])) + (SCALAR_VAL(0.125) * ((A[1][i3 + 1][128 * ii0 + 64 * ii4 - 2 * i0 + 5] - (SCALAR_VAL(2.0) * A[1][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 5])) + A[1][i3 - 1][128 * ii0 + 64 * ii4 - 2 * i0 + 5]))) + (SCALAR_VAL(0.125) * ((A[1][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 6] - (SCALAR_VAL(2.0) * A[1][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 5])) + A[1][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 4]))) + A[1][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 5]);
                          }
                        }
                        for (int i4 = max(max(1, 128 * ii0 + 64 * ii4 - 2 * i0 + 5), 128 * ii0 + 64 * ii4 - 2 * i0 - i2 + 7); i4 <= min(64 * ii4 - 1, 64 * ii4 - 2 * i3 + 64); i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = max(max(max(1, 128 * ii0 + 64 * ii4 - 2 * i0 + 5), 128 * ii0 + 64 * ii4 - 2 * i0 - i2 + 7), 64 * ii4 - 2 * i3 + 65); i4 < 64 * ii4; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = max(1, 64 * ii4); i4 <= min(min(128 * ii0 + 64 * ii4 - 2 * i0 + 67, 64 * ii4 - 2 * i3 + 64), 128 * ii0 + 64 * ii4 - 2 * i0 - i3 + 69); i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        if (i0 >= 64 * ii0 + 3) {
                          for (int i4 = max(max(1, 64 * ii4), 64 * ii4 - 2 * i3 + 65); i4 <= 128 * ii0 + 64 * ii4 - 2 * i0 - i2 - i3 + 69; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          for (int i4 = max(64 * ii4 - 2 * i3 + 65, 128 * ii0 + 64 * ii4 - 2 * i0 - i2 - i3 + 70); i4 <= 128 * ii0 + 64 * ii4 - 2 * i0 - i3 + 69; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        } else {
                          for (int i4 = max(max(1, 64 * ii4), 64 * ii4 - 2 * i3 + 65); i4 <= min(64 * ii4 + 63, 64 * ii4 - i3 + 65); i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                        for (int i4 = 128 * ii0 + 64 * ii4 - 2 * i0 - i3 + 70; i4 <= 128 * ii0 + 64 * ii4 - 2 * i0 + 67; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                      for (int i3 = 128 * ii0 - 2 * i0 - i2 + 69; i3 <= 128 * ii0 - 2 * i0 + 67; i3 += 1) {
                        for (int i4 = max(1, 128 * ii0 + 64 * ii4 - 2 * i0 + 4); i4 < 64 * ii4; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = max(1, 64 * ii4); i4 <= 128 * ii0 + 64 * ii4 - 2 * i0 + 67; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    } else {
                      for (int i3 = 1; i3 <= 63; i3 += 1) {
                        if (i2 + i3 <= 64) {
                          for (int i4 = max(1, 64 * ii4); i4 <= 64 * ii4 - i2 + 64; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (i2 == 2) {
                            B[2][i3][64 * ii4 + 63] = ((((SCALAR_VAL(0.125) * ((A[3][i3][64 * ii4 + 63] - (SCALAR_VAL(2.0) * A[2][i3][64 * ii4 + 63])) + A[1][i3][64 * ii4 + 63])) + (SCALAR_VAL(0.125) * ((A[2][i3 + 1][64 * ii4 + 63] - (SCALAR_VAL(2.0) * A[2][i3][64 * ii4 + 63])) + A[2][i3 - 1][64 * ii4 + 63]))) + (SCALAR_VAL(0.125) * ((A[2][i3][64 * ii4 + 64] - (SCALAR_VAL(2.0) * A[2][i3][64 * ii4 + 63])) + A[2][i3][64 * ii4 + 62]))) + A[2][i3][64 * ii4 + 63]);
                          }
                        } else {
                          for (int i4 = max(1, 64 * ii4); i4 <= 64 * ii4 + 63; i4 += 1) {
                            B[2][63][i4] = ((((SCALAR_VAL(0.125) * ((A[3][63][i4] - (SCALAR_VAL(2.0) * A[2][63][i4])) + A[1][63][i4])) + (SCALAR_VAL(0.125) * ((A[2][64][i4] - (SCALAR_VAL(2.0) * A[2][63][i4])) + A[2][62][i4]))) + (SCALAR_VAL(0.125) * ((A[2][63][i4 + 1] - (SCALAR_VAL(2.0) * A[2][63][i4])) + A[2][63][i4 - 1]))) + A[2][63][i4]);
                          }
                        }
                      }
                    }
                  }
                }
                for (int i2 = 1; i2 <= 128 * ii0 - 2 * i0 + 66; i2 += 1) {
                  for (int i3 = 1; i3 <= 128 * ii0 - 2 * i0 + 66; i3 += 1) {
                    for (int i4 = max(1, 128 * ii0 + 64 * ii4 - 2 * i0 + 3); i4 <= 128 * ii0 + 64 * ii4 - 2 * i0 + 66; i4 += 1) {
                      A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                    }
                  }
                }
              }
              if (_PB_TSTEPS >= 64 * ii0 + 33) {
                for (int i4 = max(1, 64 * ii4 - 62); i4 <= 64 * ii4 + 1; i4 += 1) {
                  B[1][1][i4] = ((((SCALAR_VAL(0.125) * ((A[2][1][i4] - (SCALAR_VAL(2.0) * A[1][1][i4])) + A[0][1][i4])) + (SCALAR_VAL(0.125) * ((A[1][2][i4] - (SCALAR_VAL(2.0) * A[1][1][i4])) + A[1][0][i4]))) + (SCALAR_VAL(0.125) * ((A[1][1][i4 + 1] - (SCALAR_VAL(2.0) * A[1][1][i4])) + A[1][1][i4 - 1]))) + A[1][1][i4]);
                }
              }
            } else {
              for (int i0 = 64 * ii0 + 1; i0 <= min(_PB_TSTEPS, 64 * ii0 + 32); i0 += 1) {
                if (i0 >= 64 * ii0 + 2) {
                  for (int i2 = 1; i2 <= 128 * ii0 - 2 * i0 + 67; i2 += 1) {
                    for (int i3 = 1; i3 <= 128 * ii0 - 2 * i0 + 67; i3 += 1) {
                      if (i0 >= 64 * ii0 + 3 && i3 == 1) {
                        for (int i4 = 128 * ii0 + 64 * ii4 - 2 * i0 + 4; i4 <= min(128 * ii0 + 64 * ii4 - 2 * i0 + 5, 128 * ii0 + 64 * ii4 - 2 * i0 + i2 + 2); i4 += 1) {
                          B[i2][1][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][1][i4] - (SCALAR_VAL(2.0) * A[i2][1][i4])) + A[i2 - 1][1][i4])) + (SCALAR_VAL(0.125) * ((A[i2][2][i4] - (SCALAR_VAL(2.0) * A[i2][1][i4])) + A[i2][0][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][1][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][1][i4])) + A[i2][1][i4 - 1]))) + A[i2][1][i4]);
                        }
                      }
                      for (int i4 = 128 * ii0 + 64 * ii4 - 2 * i0 + 4; i4 <= min(-128 * ii0 + 64 * ii4 + 2 * i0 + i2 - i3 - 5, 128 * ii0 + 64 * ii4 - 2 * i0 + floord(-i2 + i3 - 1, 5) + 4); i4 += 1) {
                        B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                      }
                      if (i0 >= 64 * ii0 + 3 && i2 >= 4 && i3 >= 2 && i2 >= i3) {
                        B[i2][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 4] - (SCALAR_VAL(2.0) * A[i2][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 4])) + A[i2 - 1][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][128 * ii0 + 64 * ii4 - 2 * i0 + 4] - (SCALAR_VAL(2.0) * A[i2][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 4])) + A[i2][i3 - 1][128 * ii0 + 64 * ii4 - 2 * i0 + 4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 5] - (SCALAR_VAL(2.0) * A[i2][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 4])) + A[i2][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 3]))) + A[i2][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 4]);
                      }
                      if (i3 >= 2) {
                        for (int i4 = 128 * ii0 + 64 * ii4 - 2 * i0 + floord(-i2 + i3 - 1, 5) + 5; i4 <= min(min(64 * ii4 - 1, 128 * ii0 + 64 * ii4 - 2 * i0 - i2 + 7), -128 * ii0 + 64 * ii4 + 2 * i0 + i2 - i3 - 5); i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                      for (int i4 = max(max(max(128 * ii0 + 64 * ii4 - 2 * i0 + 5, 128 * ii0 + 64 * ii4 - 2 * i0 - i2 + 8), 128 * ii0 + 64 * ii4 - 2 * i0 - i3 + 7), 128 * ii0 + 64 * ii4 - 2 * i0 + floord(-i2 + i3 - 1, 5) + 5); i4 < min(min(64 * ii4, 128 * ii0 + 64 * ii4 - 2 * i0 + i2 + 121 * i3 - 118), -128 * ii0 + 64 * ii4 + 2 * i0 + i2 - i3 - 4); i4 += 1) {
                        B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                      }
                      if (i3 == 1) {
                        for (int i4 = 128 * ii0 + 64 * ii4 - 2 * i0 + i2 + 3; i4 < 64 * ii4; i4 += 1) {
                          B[i2][1][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][1][i4] - (SCALAR_VAL(2.0) * A[i2][1][i4])) + A[i2 - 1][1][i4])) + (SCALAR_VAL(0.125) * ((A[i2][2][i4] - (SCALAR_VAL(2.0) * A[i2][1][i4])) + A[i2][0][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][1][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][1][i4])) + A[i2][1][i4 - 1]))) + A[i2][1][i4]);
                        }
                      }
                      for (int i4 = max(128 * ii0 + 64 * ii4 - 2 * i0 + 4, -128 * ii0 + 64 * ii4 + 2 * i0 + i2 - i3 - 4); i4 < min(64 * ii4, 128 * ii0 + 64 * ii4 - 2 * i0 + 124 * i2 - 118); i4 += 1) {
                        if ((i2 >= 2 && 2 * i0 + i4 >= 128 * ii0 + 64 * ii4 + 5) || 2 * i0 + i4 == 128 * ii0 + 64 * ii4 + 4 || i2 == 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                      if (i2 == 1) {
                        for (int i4 = max(128 * ii0 + 64 * ii4 - 2 * i0 + 6, -128 * ii0 + 64 * ii4 + 2 * i0 - i3 - 3); i4 < 64 * ii4; i4 += 1) {
                          B[1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[2][i3][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[0][i3][i4])) + (SCALAR_VAL(0.125) * ((A[1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3][i4 - 1]))) + A[1][i3][i4]);
                        }
                      }
                      for (int i4 = 64 * ii4; i4 < _PB_N - 1; i4 += 1) {
                        B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                      }
                    }
                  }
                }
                for (int i2 = 1; i2 <= 128 * ii0 - 2 * i0 + 66; i2 += 1) {
                  for (int i3 = 1; i3 <= 128 * ii0 - 2 * i0 + 66; i3 += 1) {
                    for (int i4 = 128 * ii0 + 64 * ii4 - 2 * i0 + 3; i4 < _PB_N - 1; i4 += 1) {
                      A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                    }
                  }
                }
              }
              if (_PB_TSTEPS >= 64 * ii0 + 33) {
                for (int i4 = 64 * ii4 - 62; i4 < _PB_N - 1; i4 += 1) {
                  B[1][1][i4] = ((((SCALAR_VAL(0.125) * ((A[2][1][i4] - (SCALAR_VAL(2.0) * A[1][1][i4])) + A[0][1][i4])) + (SCALAR_VAL(0.125) * ((A[1][2][i4] - (SCALAR_VAL(2.0) * A[1][1][i4])) + A[1][0][i4]))) + (SCALAR_VAL(0.125) * ((A[1][1][i4 + 1] - (SCALAR_VAL(2.0) * A[1][1][i4])) + A[1][1][i4 - 1]))) + A[1][1][i4]);
                }
              }
            }
          }
        }
        for (int ii3 = max(0, -ii2 + 1); ii3 <= (_PB_N - 3) / 64; ii3 += 1) {
          if (_PB_N >= 64 * ii3 + 67) {
            for (int ii4 = 0; ii4 <= (_PB_N - 3) / 64; ii4 += 1) {
              if (_PB_N >= 64 * ii2 + 67) {
                for (int i2 = 64 * ii2 + 1; i2 <= 64 * ii2 + 64; i2 += 1) {
                  for (int i3 = 64 * ii3 + 1; i3 <= 64 * ii3 + 64; i3 += 1) {
                    for (int i4 = 64 * ii4 + 1; i4 <= min(_PB_N - 2, 64 * ii4 + 64); i4 += 1) {
                      A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                    }
                  }
                }
                if (ii2 == 0 && 64 * ii4 + 66 >= _PB_N) {
                  for (int i2 = 1; i2 <= 63; i2 += 1) {
                    for (int i3 = 64 * ii3; i3 <= 64 * ii3 + 63; i3 += 1) {
                      for (int i4 = max(64 * ii4, 32 * ii3 + 64 * ii4 - i3 + i3 / 2 + 1); i4 < _PB_N - 1; i4 += 1) {
                        B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                      }
                    }
                  }
                  for (int i2 = 1; i2 <= 62; i2 += 1) {
                    for (int i3 = 64 * ii3 - 1; i3 <= 64 * ii3 + 62; i3 += 1) {
                      for (int i4 = max(64 * ii4 - 1, 64 * ii3 + 64 * ii4 - i3); i4 < _PB_N - 1; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                  }
                }
                if (64 * ii4 + 66 >= _PB_N) {
                  for (int i0 = max(64 * ii0 + 2, 64 * ii0 - ii2 + 3); i0 <= min(min(_PB_TSTEPS, 64 * ii0 + 64), 64 * ii0 + 32 * ii2 + 32); i0 += 1) {
                    if (i0 == 64 * ii0 + 2) {
                      for (int i2 = 64 * ii2; i2 <= 64 * ii2 + 1; i2 += 1) {
                        for (int i3 = max(1, 64 * ii2 + 64 * ii3 - i2 + 1); i3 <= 64 * ii2 + 64 * ii3 - i2 + 64; i3 += 1) {
                          for (int i4 = max(max(64 * ii4, 128 * ii2 + 64 * ii4 - 2 * i2 + 1), 32 * ii3 + 64 * ii4 - i3 + i3 / 2 + 1); i4 < _PB_N - 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                    } else {
                      for (int i2 = max(128 * ii0 + 64 * ii2 - 2 * i0 + 4, -ii2 + (59 * ii2 + 63) / 123 + 2); i2 <= 128 * ii0 + 64 * ii2 - 2 * i0 + 6; i2 += 1) {
                        if (2 * i0 + i2 == 128 * ii0 + 64 * ii2 + 6) {
                          for (int i3 = max(1, 64 * ii3 - 1); i3 <= 64 * ii3; i3 += 1) {
                            for (int i4 = 64 * ii3 + 64 * ii4 - i3; i4 < _PB_N - 1; i4 += 1) {
                              B[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 7][i3][i4] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][i4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][i4])) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][i4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][i4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][i4 - 1]))) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][i4]);
                            }
                          }
                        }
                        if (64 * ii2 >= i2 + 1) {
                          for (int i3 = max(max(1, -128 * ii0 - 64 * ii2 + 64 * ii3 + 2 * i0 + i2 - 5), 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 5); i3 <= min(-128 * ii0 - 64 * ii2 + 64 * ii3 + 2 * i0 + i2 - 3, 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 7); i3 += 1) {
                            if (2 * i0 + i2 == 128 * ii0 + 64 * ii2 + 5 && i3 == 64 * ii3 + 2) {
                              B[128 * ii0 + 64 * ii2 - 2 * i0 + 5][64 * ii3 + 2][64 * ii4] = ((((SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][64 * ii3 + 2][64 * ii4] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][64 * ii3 + 2][64 * ii4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][64 * ii3 + 2][64 * ii4])) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][64 * ii3 + 3][64 * ii4] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][64 * ii3 + 2][64 * ii4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][64 * ii3 + 1][64 * ii4]))) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][64 * ii3 + 2][64 * ii4 + 1] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][64 * ii3 + 2][64 * ii4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][64 * ii3 + 2][64 * ii4 - 1]))) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][64 * ii3 + 2][64 * ii4]);
                            }
                            for (int i4 = max(128 * ii0 + 64 * ii2 - 64 * ii3 + 64 * ii4 - 2 * i0 - i2 + i3 + 4, 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 6); i4 < _PB_N - 1; i4 += 1) {
                              if (128 * ii0 + 64 * ii2 + i3 + i4 + 4 >= 64 * ii3 + 64 * ii4 + 2 * i0 + i2 || 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                          }
                          if (2 * i0 + i2 == 128 * ii0 + 64 * ii2 + 4) {
                            for (int i3 = 64 * ii3 + 2; i3 <= 64 * ii3 + 3; i3 += 1) {
                              for (int i4 = 64 * ii4 + 1; i4 < _PB_N - 1; i4 += 1) {
                                B[128 * ii0 + 64 * ii2 - 2 * i0 + 4][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][i4] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][i3][i4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 3][i3][i4])) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][i3][i4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][i3][i4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][i3][i4 - 1]))) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][i3][i4]);
                              }
                            }
                          }
                          for (int i3 = 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 8; i3 <= 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 68; i3 += 1) {
                            if (2 * i0 + i2 == 128 * ii0 + 64 * ii2 + 6 && i3 == 64 * ii3 + 2) {
                              B[128 * ii0 + 64 * ii2 - 2 * i0 + 6][64 * ii3 + 2][64 * ii4 - 1] = ((((SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 7][64 * ii3 + 2][64 * ii4 - 1] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][64 * ii3 + 2][64 * ii4 - 1])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][64 * ii3 + 2][64 * ii4 - 1])) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][64 * ii3 + 3][64 * ii4 - 1] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][64 * ii3 + 2][64 * ii4 - 1])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][64 * ii3 + 1][64 * ii4 - 1]))) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][64 * ii3 + 2][64 * ii4] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][64 * ii3 + 2][64 * ii4 - 1])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][64 * ii3 + 2][64 * ii4 - 2]))) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][64 * ii3 + 2][64 * ii4 - 1]);
                            }
                            for (int i4 = max(128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 5, 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 8); i4 < _PB_N - 1; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                        } else {
                          for (int i3 = 64 * ii3 + 1; i3 <= 64 * ii3 + 62; i3 += 1) {
                            for (int i4 = 64 * ii4 - 1; i4 < _PB_N - 1; i4 += 1) {
                              if (i3 >= 64 * ii3 + 2 || i4 >= 64 * ii4 + 1 || 1) {
                                B[64 * ii2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[64 * ii2][i3][i4])) + A[64 * ii2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[64 * ii2][i3][i4])) + A[64 * ii2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii2][i3][i4])) + A[64 * ii2][i3][i4 - 1]))) + A[64 * ii2][i3][i4]);
                              }
                            }
                          }
                        }
                      }
                    }
                    for (int i2 = max(128 * ii0 + 64 * ii2 - 2 * i0 + 7, 2 * ii0 - (i0 + 29) / 32 + 3); i2 <= min(min(64 * ii2 + 1, 128 * ii0 + 64 * ii2 - 2 * i0 + 67), 128 * ii0 + 64 * ii2 - 32 * ii4 - 2 * i0 + _PB_N / 2 + 35); i2 += 1) {
                      if (i2 == 64 * ii2 + 1) {
                        for (int i3 = max(1, 128 * ii0 + 64 * ii3 - 2 * i0 + 4); i3 < min(64 * ii3, -_PB_N - 128 * ii0 + 64 * ii3 + 64 * ii4 + 2 * i0 - 1); i3 += 1) {
                          for (int i4 = 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 5; i4 < _PB_N - 1; i4 += 1) {
                            B[64 * ii2 + 1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][i3][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][i3][i4])) + A[64 * ii2][i3][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][i3][i4])) + A[64 * ii2 + 1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][i3][i4])) + A[64 * ii2 + 1][i3][i4 - 1]))) + A[64 * ii2 + 1][i3][i4]);
                          }
                        }
                      }
                      for (int i3 = max(1, 128 * ii0 + 64 * ii3 - 2 * i0 - i2 + 69); i3 <= min(128 * ii0 + 64 * ii3 - 2 * i0 + 70, -256 * ii0 + 64 * ii3 + 4 * i0 - 123 * i2 - 14); i3 += 1) {
                        for (int i4 = 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 70; i4 < _PB_N - 1; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                      if (_PB_N >= 131) {
                        for (int i3 = max(1, 128 * ii0 + 64 * ii3 - 2 * i0 + 71); i3 <= min(min(64 * ii3 + 1, -_PB_N - 128 * ii0 + 64 * ii3 + 64 * ii4 + 2 * i0 - 66), -256 * ii0 + 64 * ii3 + 4 * i0 - 123 * i2 - 14); i3 += 1) {
                          for (int i4 = 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 70; i4 < _PB_N - 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                        if (64 * ii2 >= i2) {
                          for (int i3 = max(max(1, 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 5), -256 * ii0 - 128 * ii2 + 64 * ii3 + 4 * i0 - 123 * i2 + 115); i3 <= min(64 * ii3, -_PB_N - 128 * ii0 - 64 * ii2 + 64 * ii3 + 64 * ii4 + 2 * i0 + i2 - 3); i3 += 1) {
                            for (int i4 = 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 6; i4 < min(_PB_N - 1, -128 * ii0 - 64 * ii2 - 64 * ii3 + 64 * ii4 + 2 * i0 + i2 + i3 - 6); i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            for (int i4 = max(-128 * ii0 - 64 * ii2 - 64 * ii3 + 64 * ii4 + 2 * i0 + i2 + i3 - 6, 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 6); i4 <= 64 * ii4 + 1; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            if (2 * i0 + i2 + i3 >= 128 * ii0 + 64 * ii2 + 64 * ii3 + 8) {
                              for (int i4 = -128 * ii0 - 64 * ii2 - 64 * ii3 + 64 * ii4 + 2 * i0 + i2 + i3 - 6; i4 < _PB_N - 1; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            } else {
                              for (int i4 = 64 * ii4 + 2; i4 < _PB_N - 1; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                          }
                        }
                        for (int i3 = max(max(1, -_PB_N - 128 * ii0 - 64 * ii2 + 64 * ii3 + 64 * ii4 + 2 * i0 + i2 - 2), 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 5); i3 <= min(64 * ii3 - 1, 128 * ii0 + 64 * ii3 - 2 * i0 + 6); i3 += 1) {
                          for (int i4 = 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 6; i4 < _PB_N - 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                        for (int i3 = max(max(max(max(1, 128 * ii0 + 64 * ii3 - 2 * i0 + 7), -_PB_N - 128 * ii0 - 64 * ii2 + 64 * ii3 + 64 * ii4 + 2 * i0 + i2 - 2), 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 5), -256 * ii0 - 128 * ii2 + 64 * ii3 + 4 * i0 - 123 * i2 + 115); i3 <= 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 7; i3 += 1) {
                          if (2 * i0 + i2 + i3 == 128 * ii0 + 64 * ii2 + 64 * ii3 + 7) {
                            for (int i4 = 64 * ii4 - 1; i4 <= 64 * ii4; i4 += 1) {
                              B[i2][128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 7][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 7][i4] - (SCALAR_VAL(2.0) * A[i2][128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 7][i4])) + A[i2 - 1][128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 7][i4])) + (SCALAR_VAL(0.125) * ((A[i2][128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 8][i4] - (SCALAR_VAL(2.0) * A[i2][128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 7][i4])) + A[i2][128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 6][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 7][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 7][i4])) + A[i2][128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 7][i4 - 1]))) + A[i2][128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 7][i4]);
                            }
                          }
                          for (int i4 = max(-128 * ii0 - 64 * ii2 - 64 * ii3 + 64 * ii4 + 2 * i0 + i2 + i3 - 6, 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 6); i4 <= 64 * ii4 + 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          for (int i4 = 64 * ii4 + 2; i4 < _PB_N - 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                        if (64 * ii4 + 2 * i0 + i2 >= _PB_N + 128 * ii0 + 64 * ii2 + 4 && 64 * ii2 >= i2 + 1 && 256 * ii0 + 128 * ii2 + 123 * i2 >= 4 * i0 + 114) {
                          for (int i4 = 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 5; i4 < _PB_N - 1; i4 += 1) {
                            B[i2][64 * ii3 + 1][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 1][i4])) + A[i2 - 1][64 * ii3 + 1][i4])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 2][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 1][i4])) + A[i2][64 * ii3][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 1][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 1][i4])) + A[i2][64 * ii3 + 1][i4 - 1]))) + A[i2][64 * ii3 + 1][i4]);
                          }
                        } else {
                          if (ii3 >= 1 && _PB_N >= 64 * ii4 + 5 && i0 == 64 * ii0 + 3 && i2 == 64 * ii2 + 1) {
                            for (int i4 = 64 * ii4 - 1; i4 < _PB_N - 1; i4 += 1) {
                              B[64 * ii2 + 1][64 * ii3][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][64 * ii3][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][64 * ii3][i4])) + A[64 * ii2][64 * ii3][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][64 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][64 * ii3][i4])) + A[64 * ii2 + 1][64 * ii3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][64 * ii3][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][64 * ii3][i4])) + A[64 * ii2 + 1][64 * ii3][i4 - 1]))) + A[64 * ii2 + 1][64 * ii3][i4]);
                            }
                          }
                          if (i2 == 64 * ii2 + 1) {
                            for (int i3 = max(1, 64 * ii3); i3 <= min(64 * ii3 + 1, -_PB_N - 128 * ii0 + 64 * ii3 + 64 * ii4 + 2 * i0 - 2); i3 += 1) {
                              for (int i4 = 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 5; i4 < min(_PB_N - 1, -128 * ii0 - 64 * ii3 + 64 * ii4 + 2 * i0 + i3 - 5); i4 += 1) {
                                B[64 * ii2 + 1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][i3][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][i3][i4])) + A[64 * ii2][i3][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][i3][i4])) + A[64 * ii2 + 1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][i3][i4])) + A[64 * ii2 + 1][i3][i4 - 1]))) + A[64 * ii2 + 1][i3][i4]);
                              }
                              if (i3 == 64 * ii3) {
                                for (int i4 = -128 * ii0 + 64 * ii4 + 2 * i0 - 5; i4 < _PB_N - 1; i4 += 1) {
                                  B[64 * ii2 + 1][64 * ii3][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][64 * ii3][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][64 * ii3][i4])) + A[64 * ii2][64 * ii3][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][64 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][64 * ii3][i4])) + A[64 * ii2 + 1][64 * ii3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][64 * ii3][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][64 * ii3][i4])) + A[64 * ii2 + 1][64 * ii3][i4 - 1]))) + A[64 * ii2 + 1][64 * ii3][i4]);
                                }
                              }
                            }
                          } else if (64 * ii4 + 2 * i0 >= _PB_N + 128 * ii0 + 4 && i2 == 64 * ii2) {
                            for (int i4 = 128 * ii0 + 64 * ii4 - 2 * i0 + 5; i4 < _PB_N - 1; i4 += 1) {
                              B[64 * ii2][64 * ii3 + 1][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][64 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[64 * ii2][64 * ii3 + 1][i4])) + A[64 * ii2 - 1][64 * ii3 + 1][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii2][64 * ii3 + 2][i4] - (SCALAR_VAL(2.0) * A[64 * ii2][64 * ii3 + 1][i4])) + A[64 * ii2][64 * ii3][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2][64 * ii3 + 1][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii2][64 * ii3 + 1][i4])) + A[64 * ii2][64 * ii3 + 1][i4 - 1]))) + A[64 * ii2][64 * ii3 + 1][i4]);
                            }
                          }
                        }
                        for (int i3 = max(64 * ii3 + 2, -256 * ii0 - 128 * ii2 + 64 * ii3 + 4 * i0 - 123 * i2 + 115); i3 <= min(-_PB_N - 128 * ii0 - 64 * ii2 + 64 * ii3 + 64 * ii4 + 2 * i0 + i2 - 3, 256 * ii2 + 64 * ii3 - 4 * i2 + 15); i3 += 1) {
                          for (int i4 = 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 5; i4 <= 128 * ii0 + 13 * ii2 - 13 * ii3 + 64 * ii4 - 2 * i0 + floord(-ii2 + ii3 - i2 + i3 - 1, 5) + 4; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          for (int i4 = max(128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 5, 128 * ii0 + 13 * ii2 - 13 * ii3 + 64 * ii4 - 2 * i0 + floord(-ii2 + ii3 - i2 + i3 - 1, 5) + 5); i4 <= 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 7; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (i2 >= 64 * ii2) {
                            for (int i4 = 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 8; i4 <= 128 * ii2 + 64 * ii4 - 2 * i2; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            if (i2 == 64 * ii2 + 1) {
                              B[64 * ii2 + 1][i3][64 * ii4 - 1] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][i3][64 * ii4 - 1] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][i3][64 * ii4 - 1])) + A[64 * ii2][i3][64 * ii4 - 1])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][i3 + 1][64 * ii4 - 1] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][i3][64 * ii4 - 1])) + A[64 * ii2 + 1][i3 - 1][64 * ii4 - 1]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][i3][64 * ii4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][i3][64 * ii4 - 1])) + A[64 * ii2 + 1][i3][64 * ii4 - 2]))) + A[64 * ii2 + 1][i3][64 * ii4 - 1]);
                            }
                            for (int i4 = max(64 * ii4, 128 * ii2 + 64 * ii4 - 2 * i2 + 1); i4 < _PB_N - 1; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          } else {
                            for (int i4 = 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 8; i4 < _PB_N - 1; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                        }
                        for (int i3 = 256 * ii2 + 64 * ii3 - 4 * i2 + 16; i3 < -_PB_N - 128 * ii0 - 64 * ii2 + 64 * ii3 + 64 * ii4 + 2 * i0 + i2 - 2; i3 += 1) {
                          for (int i4 = 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 5; i4 <= 128 * ii0 + 13 * ii2 - 13 * ii3 + 64 * ii4 - 2 * i0 + floord(-ii2 + ii3 - i2 + i3 - 1, 5) + 4; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          for (int i4 = 128 * ii0 + 13 * ii2 - 13 * ii3 + 64 * ii4 - 2 * i0 + floord(-ii2 + ii3 - i2 + i3 - 1, 5) + 5; i4 < _PB_N - 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                        if (ii2 == 1 && i2 == 1) {
                          for (int i3 = max(max(1, -_PB_N - 128 * ii0 + 64 * ii3 + 64 * ii4 + 2 * i0 - 65), 128 * ii0 + 64 * ii3 - 2 * i0 + 71); i3 <= 64 * ii3 + 1; i3 += 1) {
                            for (int i4 = 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 69; i4 < _PB_N - 1; i4 += 1) {
                              B[1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[2][i3][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[0][i3][i4])) + (SCALAR_VAL(0.125) * ((A[1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3][i4 - 1]))) + A[1][i3][i4]);
                            }
                          }
                          for (int i3 = 64 * ii3 + 2; i3 <= min(-256 * ii0 + 64 * ii3 + 4 * i0 - 137, 128 * ii0 + 64 * ii3 - 2 * i0 + 131); i3 += 1) {
                            for (int i4 = 128 * ii0 + 64 * ii4 - 2 * i0 + 68; i4 < _PB_N - 1; i4 += 1) {
                              B[1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[2][i3][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[0][i3][i4])) + (SCALAR_VAL(0.125) * ((A[1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3][i4 - 1]))) + A[1][i3][i4]);
                            }
                          }
                        } else if (i2 == 64 * ii2 + 1) {
                          for (int i3 = max(max(1, -_PB_N - 128 * ii0 + 64 * ii3 + 64 * ii4 + 2 * i0 - 1), 128 * ii0 + 64 * ii3 - 2 * i0 + 7); i3 <= min(-256 * ii0 + 64 * ii3 + 4 * i0 - 8, 128 * ii0 + 64 * ii3 - 2 * i0 + 67); i3 += 1) {
                            if (i3 >= 64 * ii3 + 2) {
                              for (int i4 = 128 * ii0 + 64 * ii4 - 2 * i0 + 4; i4 < -128 * ii0 + 64 * ii3 + 64 * ii4 + 2 * i0 - i3 - 3; i4 += 1) {
                                B[64 * ii2 + 1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][i3][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][i3][i4])) + A[64 * ii2][i3][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][i3][i4])) + A[64 * ii2 + 1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][i3][i4])) + A[64 * ii2 + 1][i3][i4 - 1]))) + A[64 * ii2 + 1][i3][i4]);
                              }
                            } else {
                              for (int i4 = 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 5; i4 < -128 * ii0 + 64 * ii3 + 64 * ii4 + 2 * i0 - i3 - 3; i4 += 1) {
                                B[64 * ii2 + 1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][i3][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][i3][i4])) + A[64 * ii2][i3][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][i3][i4])) + A[64 * ii2 + 1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][i3][i4])) + A[64 * ii2 + 1][i3][i4 - 1]))) + A[64 * ii2 + 1][i3][i4]);
                              }
                            }
                            for (int i4 = -128 * ii0 + 64 * ii3 + 64 * ii4 + 2 * i0 - i3 - 3; i4 < _PB_N - 1; i4 += 1) {
                              B[64 * ii2 + 1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][i3][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][i3][i4])) + A[64 * ii2][i3][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][i3][i4])) + A[64 * ii2 + 1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][i3][i4])) + A[64 * ii2 + 1][i3][i4 - 1]))) + A[64 * ii2 + 1][i3][i4]);
                            }
                          }
                        }
                        if (64 * ii2 >= i2) {
                          for (int i3 = max(max(max(1, -_PB_N - 128 * ii0 - 64 * ii2 + 64 * ii3 + 64 * ii4 + 2 * i0 + i2 - 2), 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 8), -256 * ii0 - 128 * ii2 + 64 * ii3 + 4 * i0 - 123 * i2 + 115); i3 <= 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 68; i3 += 1) {
                            for (int i4 = 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 5; i4 <= min(-128 * ii0 - 64 * ii2 + 64 * ii3 + 64 * ii4 + 2 * i0 + i2 - i3 - 5, 128 * ii0 + 13 * ii2 - 13 * ii3 + 64 * ii4 - 2 * i0 + floord(-ii2 + ii3 - i2 + i3 - 1, 5) + 4); i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            if (i3 >= 64 * ii3 + 2) {
                              for (int i4 = max(128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 5, 128 * ii0 + 13 * ii2 - 13 * ii3 + 64 * ii4 - 2 * i0 + floord(-ii2 + ii3 - i2 + i3 - 1, 5) + 5); i4 <= min(128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 7, -128 * ii0 - 64 * ii2 + 64 * ii3 + 64 * ii4 + 2 * i0 + i2 - i3 - 5); i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                            for (int i4 = max(128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 5, -128 * ii0 - 64 * ii2 + 64 * ii3 + 64 * ii4 + 2 * i0 + i2 - i3 - 4); i4 <= 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 6; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            if (i3 >= 64 * ii3 + 2) {
                              for (int i4 = max(128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 8, 128 * ii0 + 13 * ii2 - 13 * ii3 + 64 * ii4 - 2 * i0 + floord(-ii2 + ii3 - i2 + i3 - 1, 5) + 5); i4 < -128 * ii0 - 64 * ii2 + 64 * ii3 + 64 * ii4 + 2 * i0 + i2 - i3 - 4; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            } else {
                              for (int i4 = 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 6; i4 < -128 * ii0 - 64 * ii2 + 64 * ii3 + 64 * ii4 + 2 * i0 + i2 - i3 - 4; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                            for (int i4 = max(128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 7, -128 * ii0 - 64 * ii2 + 64 * ii3 + 64 * ii4 + 2 * i0 + i2 - i3 - 4); i4 < _PB_N - 1; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                        } else if (i2 == 64 * ii2 + 1) {
                          for (int i3 = -256 * ii0 + 64 * ii3 + 4 * i0 - 7; i3 <= 128 * ii0 + 64 * ii3 - 2 * i0 + 67; i3 += 1) {
                            for (int i4 = 128 * ii0 + 64 * ii4 - 2 * i0 + 4; i4 < _PB_N - 1; i4 += 1) {
                              B[64 * ii2 + 1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][i3][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][i3][i4])) + A[64 * ii2][i3][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][i3][i4])) + A[64 * ii2 + 1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][i3][i4])) + A[64 * ii2 + 1][i3][i4 - 1]))) + A[64 * ii2 + 1][i3][i4]);
                            }
                          }
                        }
                      }
                    }
                    for (int i2 = max(128 * ii0 + 64 * ii2 - 32 * ii4 - 2 * i0 + _PB_N / 2 + 36, -ii2 + (29 * ii2 + 31) / 61 + 2); i2 <= min(64 * ii2 + 1, 128 * ii0 + 64 * ii2 - 2 * i0 + 67); i2 += 1) {
                      for (int i3 = max(1, 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 5); i3 <= min(64 * ii3 + 1, 64 * ii2 + 64 * ii3 - i2); i3 += 1) {
                        for (int i4 = 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 6; i4 < min(_PB_N - 1, -128 * ii0 - 64 * ii2 - 64 * ii3 + 64 * ii4 + 2 * i0 + i2 + i3 - 6); i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = max(-128 * ii0 - 64 * ii2 - 64 * ii3 + 64 * ii4 + 2 * i0 + i2 + i3 - 6, 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 6); i4 <= 64 * ii4 + 1; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        if (ii2 == 1 && i2 == 1 && 2 * i0 + i3 >= 128 * ii0 + 64 * ii3 + 69 && 128 * ii0 + 64 * ii3 + 70 >= 2 * i0 + i3) {
                          for (int i4 = 64 * ii4 + 2; i4 < _PB_N - 1; i4 += 1) {
                            B[1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[2][i3][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[0][i3][i4])) + (SCALAR_VAL(0.125) * ((A[1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3][i4 - 1]))) + A[1][i3][i4]);
                          }
                        } else if (2 * i0 + i2 + i3 >= 128 * ii0 + 64 * ii2 + 64 * ii3 + 8) {
                          for (int i4 = -128 * ii0 - 64 * ii2 - 64 * ii3 + 64 * ii4 + 2 * i0 + i2 + i3 - 6; i4 < _PB_N - 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        } else {
                          for (int i4 = 64 * ii4 + 2; i4 < _PB_N - 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                      for (int i3 = max(1, 64 * ii2 + 64 * ii3 - i2 + 1); i3 <= 64 * ii3 + 1; i3 += 1) {
                        for (int i4 = 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 6; i4 < _PB_N - 1; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                      for (int i3 = 64 * ii3 + 2; i3 <= 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 68; i3 += 1) {
                        for (int i4 = 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 5; i4 <= 128 * ii0 + 13 * ii2 - 13 * ii3 + 64 * ii4 - 2 * i0 + floord(-ii2 + ii3 - i2 + i3 - 1, 5) + 4; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = max(128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 5, 128 * ii0 + 13 * ii2 - 13 * ii3 + 64 * ii4 - 2 * i0 + floord(-ii2 + ii3 - i2 + i3 - 1, 5) + 5); i4 <= 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 7; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = max(128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 8, 128 * ii0 + 13 * ii2 - 13 * ii3 + 64 * ii4 - 2 * i0 + floord(-ii2 + ii3 - i2 + i3 - 1, 5) + 5); i4 < _PB_N - 1; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    }
                    if (ii2 == 0 && i0 >= 64 * ii0 + 13) {
                      for (int i3 = 128 * ii0 + 64 * ii3 - 2 * i0 + 4; i3 <= 128 * ii0 + 64 * ii3 - 2 * i0 + 67; i3 += 1) {
                        for (int i4 = max(128 * ii0 + 64 * ii4 - 2 * i0 + 4, 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 5); i4 < min(min(_PB_N - 1, -128 * ii0 - 64 * ii3 + 64 * ii4 + 2 * i0 + i3 - 5), -128 * ii0 + 64 * ii3 + 64 * ii4 + 2 * i0 - i3 - 3); i4 += 1) {
                          B[1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[2][i3][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[0][i3][i4])) + (SCALAR_VAL(0.125) * ((A[1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3][i4 - 1]))) + A[1][i3][i4]);
                        }
                        if (2 * i0 + i3 >= 128 * ii0 + 64 * ii3 + 7) {
                          for (int i4 = -128 * ii0 - 64 * ii3 + 64 * ii4 + 2 * i0 + i3 - 5; i4 < min(_PB_N - 1, -128 * ii0 + 64 * ii3 + 64 * ii4 + 2 * i0 - i3 - 3); i4 += 1) {
                            B[1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[2][i3][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[0][i3][i4])) + (SCALAR_VAL(0.125) * ((A[1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3][i4 - 1]))) + A[1][i3][i4]);
                          }
                          for (int i4 = -128 * ii0 + 64 * ii3 + 64 * ii4 + 2 * i0 - i3 - 3; i4 < _PB_N - 1; i4 += 1) {
                            B[1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[2][i3][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[0][i3][i4])) + (SCALAR_VAL(0.125) * ((A[1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3][i4 - 1]))) + A[1][i3][i4]);
                          }
                        } else {
                          for (int i4 = max(-128 * ii0 - 64 * ii3 + 64 * ii4 + 2 * i0 + i3 - 5, 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 5); i4 < _PB_N - 1; i4 += 1) {
                            B[1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[2][i3][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[0][i3][i4])) + (SCALAR_VAL(0.125) * ((A[1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3][i4 - 1]))) + A[1][i3][i4]);
                          }
                        }
                      }
                    } else if (ii2 == 0 && 64 * ii0 + 12 >= i0) {
                      for (int i3 = 128 * ii0 + 64 * ii3 - 2 * i0 + 4; i3 <= 128 * ii0 + 64 * ii3 - 2 * i0 + 67; i3 += 1) {
                        for (int i4 = max(128 * ii0 + 64 * ii4 - 2 * i0 + 4, 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 5); i4 < min(min(_PB_N - 1, -128 * ii0 - 64 * ii3 + 64 * ii4 + 2 * i0 + i3 - 5), -128 * ii0 + 64 * ii3 + 64 * ii4 + 2 * i0 - i3 - 3); i4 += 1) {
                          B[1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[2][i3][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[0][i3][i4])) + (SCALAR_VAL(0.125) * ((A[1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3][i4 - 1]))) + A[1][i3][i4]);
                        }
                        if (2 * i0 + i3 >= 128 * ii0 + 64 * ii3 + 7) {
                          for (int i4 = -128 * ii0 - 64 * ii3 + 64 * ii4 + 2 * i0 + i3 - 5; i4 < min(_PB_N - 1, -128 * ii0 + 64 * ii3 + 64 * ii4 + 2 * i0 - i3 - 3); i4 += 1) {
                            B[1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[2][i3][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[0][i3][i4])) + (SCALAR_VAL(0.125) * ((A[1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3][i4 - 1]))) + A[1][i3][i4]);
                          }
                          for (int i4 = max(128 * ii0 + 64 * ii4 - 2 * i0 + 4, -128 * ii0 + 64 * ii3 + 64 * ii4 + 2 * i0 - i3 - 3); i4 < _PB_N - 1; i4 += 1) {
                            B[1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[2][i3][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[0][i3][i4])) + (SCALAR_VAL(0.125) * ((A[1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3][i4 - 1]))) + A[1][i3][i4]);
                          }
                        } else {
                          for (int i4 = max(-128 * ii0 - 64 * ii3 + 64 * ii4 + 2 * i0 + i3 - 5, 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 5); i4 < _PB_N - 1; i4 += 1) {
                            B[1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[2][i3][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[0][i3][i4])) + (SCALAR_VAL(0.125) * ((A[1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3][i4 - 1]))) + A[1][i3][i4]);
                          }
                        }
                      }
                    }
                    for (int i2 = 64 * ii2 + 2; i2 <= 128 * ii0 + 64 * ii2 - 2 * i0 + 67; i2 += 1) {
                      for (int i3 = max(1, 128 * ii0 + 64 * ii3 - 2 * i0 + 4); i3 < 64 * ii3; i3 += 1) {
                        if (2 * i0 + i3 >= 128 * ii0 + 64 * ii3 + 5) {
                          for (int i4 = max(128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 5, 128 * ii0 + 192 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - 3 * i2 - i3 + 13); i4 <= min(_PB_N - 2, 128 * ii0 - 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 + i2 - i3 + 3); i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (i2 == 64 * ii2 + 2) {
                            B[64 * ii2 + 2][i3][128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 5] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 3][i3][128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 5] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][i3][128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 5])) + A[64 * ii2 + 1][i3][128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 5])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][i3 + 1][128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 5] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][i3][128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 5])) + A[64 * ii2 + 2][i3 - 1][128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 5]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][i3][128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 6] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][i3][128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 5])) + A[64 * ii2 + 2][i3][128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 4]))) + A[64 * ii2 + 2][i3][128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 5]);
                          }
                        } else {
                          for (int i4 = 64 * ii4 + 1; i4 < min(_PB_N - 1, -64 * ii2 + 64 * ii4 + i2); i4 += 1) {
                            B[i2][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4] - (SCALAR_VAL(2.0) * A[i2][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4])) + A[i2 - 1][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4])) + (SCALAR_VAL(0.125) * ((A[i2][128 * ii0 + 64 * ii3 - 2 * i0 + 5][i4] - (SCALAR_VAL(2.0) * A[i2][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4])) + A[i2][128 * ii0 + 64 * ii3 - 2 * i0 + 3][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4])) + A[i2][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4 - 1]))) + A[i2][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4]);
                          }
                        }
                        for (int i4 = 128 * ii0 - 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 + i2 - i3 + 4; i4 < min(_PB_N - 1, -128 * ii0 - 64 * ii2 - 64 * ii3 + 64 * ii4 + 2 * i0 + i2 + i3 - 6); i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        if (128 * ii0 + 64 * ii3 + i2 + 2 >= 64 * ii2 + 2 * i0 + i3) {
                          for (int i4 = max(-128 * ii0 - 64 * ii2 - 64 * ii3 + 64 * ii4 + 2 * i0 + i2 + i3 - 6, 128 * ii0 - 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 + i2 - i3 + 4); i4 < min(_PB_N - 1, -128 * ii0 - 64 * ii2 + 64 * ii3 + 64 * ii4 + 2 * i0 + i2 - i3 - 4); i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        } else if (2 * i0 + i2 + i3 >= 128 * ii0 + 64 * ii2 + 64 * ii3 + 8) {
                          for (int i4 = -128 * ii0 - 64 * ii2 - 64 * ii3 + 64 * ii4 + 2 * i0 + i2 + i3 - 6; i4 < min(_PB_N - 1, -128 * ii0 - 64 * ii2 + 64 * ii3 + 64 * ii4 + 2 * i0 + i2 - i3 - 4); i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                        if (2 * i0 + i2 + i3 >= 128 * ii0 + 64 * ii2 + 64 * ii3 + 8 && 2 * i0 + i3 >= 128 * ii0 + 64 * ii3 + 5) {
                          for (int i4 = -128 * ii0 - 64 * ii2 + 64 * ii3 + 64 * ii4 + 2 * i0 + i2 - i3 - 4; i4 < _PB_N - 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        } else if (2 * i0 + i3 == 128 * ii0 + 64 * ii3 + 4) {
                          for (int i4 = -256 * ii0 - 64 * ii2 + 64 * ii4 + 4 * i0 + i2 - 8; i4 < _PB_N - 1; i4 += 1) {
                            B[i2][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4] - (SCALAR_VAL(2.0) * A[i2][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4])) + A[i2 - 1][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4])) + (SCALAR_VAL(0.125) * ((A[i2][128 * ii0 + 64 * ii3 - 2 * i0 + 5][i4] - (SCALAR_VAL(2.0) * A[i2][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4])) + A[i2][128 * ii0 + 64 * ii3 - 2 * i0 + 3][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4])) + A[i2][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4 - 1]))) + A[i2][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4]);
                          }
                        } else {
                          for (int i4 = 64 * ii4 + 1; i4 < _PB_N - 1; i4 += 1) {
                            B[64 * ii2 + 2][128 * ii0 + 64 * ii3 - 2 * i0 + 5][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 3][128 * ii0 + 64 * ii3 - 2 * i0 + 5][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][128 * ii0 + 64 * ii3 - 2 * i0 + 5][i4])) + A[64 * ii2 + 1][128 * ii0 + 64 * ii3 - 2 * i0 + 5][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][128 * ii0 + 64 * ii3 - 2 * i0 + 6][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][128 * ii0 + 64 * ii3 - 2 * i0 + 5][i4])) + A[64 * ii2 + 2][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][128 * ii0 + 64 * ii3 - 2 * i0 + 5][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][128 * ii0 + 64 * ii3 - 2 * i0 + 5][i4])) + A[64 * ii2 + 2][128 * ii0 + 64 * ii3 - 2 * i0 + 5][i4 - 1]))) + A[64 * ii2 + 2][128 * ii0 + 64 * ii3 - 2 * i0 + 5][i4]);
                          }
                        }
                      }
                      for (int i3 = max(1, 64 * ii3); i3 <= 128 * ii0 + 64 * ii3 - 2 * i0 + 67; i3 += 1) {
                        if (i0 >= 64 * ii0 + 3 && i3 == 64 * ii3 + 1) {
                          for (int i4 = 128 * ii0 + 64 * ii4 - 2 * i0 + 4; i4 <= min(128 * ii0 + 64 * ii4 - 2 * i0 + 5, 128 * ii0 - 64 * ii2 + 64 * ii4 - 2 * i0 + i2 + 2); i4 += 1) {
                            B[i2][64 * ii3 + 1][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 1][i4])) + A[i2 - 1][64 * ii3 + 1][i4])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 2][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 1][i4])) + A[i2][64 * ii3][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 1][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 1][i4])) + A[i2][64 * ii3 + 1][i4 - 1]))) + A[i2][64 * ii3 + 1][i4]);
                          }
                        }
                        if (64 * ii3 + 1 >= i3) {
                          for (int i4 = max(128 * ii0 - 64 * ii3 + 64 * ii4 - 2 * i0 + i3 + 5, 128 * ii0 + 192 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - 3 * i2 - i3 + 13); i4 <= min(128 * ii0 - 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 + i2 - i3 + 3, 64 * ii3 + 64 * ii4 - i3); i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        } else {
                          for (int i4 = 128 * ii0 + 64 * ii4 - 2 * i0 + 4; i4 < min(64 * ii4, -128 * ii0 - 64 * ii2 + 64 * ii3 + 64 * ii4 + 2 * i0 + i2 - i3 - 4); i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                        for (int i4 = max(128 * ii0 + 64 * ii4 - 2 * i0 + 4, -128 * ii0 - 64 * ii2 + 64 * ii3 + 64 * ii4 + 2 * i0 + i2 - i3 - 4); i4 < 64 * ii4; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        if (i0 >= 64 * ii0 + 3 && i2 == 64 * ii2 + 2 && i3 == 64 * ii3) {
                          B[64 * ii2 + 2][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 5] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 3][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 5] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 5])) + A[64 * ii2 + 1][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 5])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][64 * ii3 + 1][128 * ii0 + 64 * ii4 - 2 * i0 + 5] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 5])) + A[64 * ii2 + 2][64 * ii3 - 1][128 * ii0 + 64 * ii4 - 2 * i0 + 5]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 6] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 5])) + A[64 * ii2 + 2][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 4]))) + A[64 * ii2 + 2][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 5]);
                        }
                        if (64 * ii3 + 1 >= i3) {
                          for (int i4 = 128 * ii0 - 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 + i2 - i3 + 4; i4 <= 64 * ii3 + 64 * ii4 - i3; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                        for (int i4 = max(64 * ii4, 32 * ii3 + 64 * ii4 - i3 + i3 / 2 + 1); i4 < _PB_N - 1; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    }
                    for (int i2 = max(1, 128 * ii0 + 64 * ii2 - 2 * i0 + 3); i2 <= min(64 * ii2, 128 * ii0 + 64 * ii2 - 2 * i0 + 66); i2 += 1) {
                      for (int i3 = max(1, 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 4); i3 <= 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 67; i3 += 1) {
                        for (int i4 = max(128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 4, 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 5); i4 < _PB_N - 1; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                    for (int i2 = 64 * ii2 + 1; i2 <= 128 * ii0 + 64 * ii2 - 2 * i0 + 66; i2 += 1) {
                      for (int i3 = max(1, 128 * ii0 + 64 * ii3 - 2 * i0 + 3); i3 <= 128 * ii0 + 64 * ii3 - 2 * i0 + 66; i3 += 1) {
                        for (int i4 = max(128 * ii0 + 64 * ii4 - 2 * i0 + 3, 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 4); i4 < _PB_N - 1; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                  }
                  if (_PB_TSTEPS >= 64 * ii0 + 33 && ii2 == 0) {
                    for (int i3 = 64 * ii3 - 62; i3 <= 64 * ii3 + 1; i3 += 1) {
                      for (int i4 = 64 * ii3 + 64 * ii4 - i3 - 61; i4 < _PB_N - 1; i4 += 1) {
                        B[1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[2][i3][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[0][i3][i4])) + (SCALAR_VAL(0.125) * ((A[1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3][i4 - 1]))) + A[1][i3][i4]);
                      }
                    }
                  }
                } else {
                  for (int i0 = 64 * ii0 + 2; i0 <= min(min(_PB_TSTEPS, 64 * ii0 + 64), 64 * ii0 + 32 * ii2 + 33); i0 += 1) {
                    if (i0 >= 64 * ii0 + 3) {
                      for (int i2 = max(1, 128 * ii0 + 64 * ii2 - 2 * i0 + 4); i2 <= min(64 * ii0 + 64 * ii2 - i0 + 2, 128 * ii0 + 64 * ii2 - 2 * i0 + 6); i2 += 1) {
                        if (128 * ii0 + 64 * ii2 + 5 >= 2 * i0 + i2) {
                          for (int i3 = max(1, 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 5); i3 <= 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 68; i3 += 1) {
                            for (int i4 = max(max(1, 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 5), 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 6); i4 <= 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 6; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            for (int i4 = 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 7; i4 <= -128 * ii0 - 64 * ii2 + 64 * ii4 + 2 * i0 + i2 + 58; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            if (2 * i0 + i2 == 128 * ii0 + 64 * ii2 + 4) {
                              for (int i4 = 64 * ii4 + 63; i4 <= 64 * ii4 + 64; i4 += 1) {
                                B[128 * ii0 + 64 * ii2 - 2 * i0 + 4][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][i4] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][i3][i4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 3][i3][i4])) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][i3][i4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][i3][i4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][i3][i4 - 1]))) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][i3][i4]);
                              }
                            } else if (i3 == 64 * ii3) {
                              B[128 * ii0 + 64 * ii2 - 2 * i0 + 5][64 * ii3][64 * ii4 + 64] = ((((SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][64 * ii3][64 * ii4 + 64] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][64 * ii3][64 * ii4 + 64])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][64 * ii3][64 * ii4 + 64])) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][64 * ii3 + 1][64 * ii4 + 64] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][64 * ii3][64 * ii4 + 64])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][64 * ii3 - 1][64 * ii4 + 64]))) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][64 * ii3][64 * ii4 + 65] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][64 * ii3][64 * ii4 + 64])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][64 * ii3][64 * ii4 + 63]))) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][64 * ii3][64 * ii4 + 64]);
                            }
                          }
                        } else {
                          for (int i3 = max(1, 64 * ii3 - 1); i3 <= 64 * ii3 + 1; i3 += 1) {
                            for (int i4 = max(1, 64 * ii3 + 64 * ii4 - i3); i4 <= 64 * ii3 + 64 * ii4 - i3 + 63; i4 += 1) {
                              B[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 7][i3][i4] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][i4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][i4])) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][i4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][i4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][i4 - 1]))) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][i4]);
                            }
                          }
                          for (int i3 = 64 * ii3 + 2; i3 <= 64 * ii3 + 62; i3 += 1) {
                            for (int i4 = max(1, 64 * ii4 - 1); i4 <= 64 * ii4 + 62; i4 += 1) {
                              B[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 7][i3][i4] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][i4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][i4])) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][i4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][i4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][i4 - 1]))) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][i4]);
                            }
                          }
                        }
                      }
                      if (ii2 >= 1 && i0 == 64 * ii0 + 3) {
                        for (int i3 = max(1, 64 * ii3 - 1); i3 <= 64 * ii3 + 62; i3 += 1) {
                          for (int i4 = max(max(1, 64 * ii4 - 1), 64 * ii3 + 64 * ii4 - i3); i4 <= 64 * ii3 + 64 * ii4 - i3 + 63; i4 += 1) {
                            B[64 * ii2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[64 * ii2][i3][i4])) + A[64 * ii2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[64 * ii2][i3][i4])) + A[64 * ii2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii2][i3][i4])) + A[64 * ii2][i3][i4 - 1]))) + A[64 * ii2][i3][i4]);
                          }
                          for (int i4 = 64 * ii3 + 64 * ii4 - i3 + 64; i4 <= 64 * ii4 + 62; i4 += 1) {
                            B[64 * ii2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[64 * ii2][i3][i4])) + A[64 * ii2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[64 * ii2][i3][i4])) + A[64 * ii2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii2][i3][i4])) + A[64 * ii2][i3][i4 - 1]))) + A[64 * ii2][i3][i4]);
                          }
                        }
                      } else if (ii2 == 1 && i0 >= 64 * ii0 + 35) {
                        for (int i3 = max(1, 128 * ii0 + 64 * ii3 - 2 * i0 + 68); i3 <= 64 * ii3 + 1; i3 += 1) {
                          for (int i4 = max(1, 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 69); i4 <= 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 132; i4 += 1) {
                            B[1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[2][i3][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[0][i3][i4])) + (SCALAR_VAL(0.125) * ((A[1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3][i4 - 1]))) + A[1][i3][i4]);
                          }
                        }
                        for (int i3 = 64 * ii3 + 2; i3 <= 128 * ii0 + 64 * ii3 - 2 * i0 + 131; i3 += 1) {
                          for (int i4 = max(1, 128 * ii0 + 64 * ii4 - 2 * i0 + 68); i4 <= 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 132; i4 += 1) {
                            B[1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[2][i3][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[0][i3][i4])) + (SCALAR_VAL(0.125) * ((A[1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3][i4 - 1]))) + A[1][i3][i4]);
                          }
                          for (int i4 = 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 133; i4 <= 128 * ii0 + 64 * ii4 - 2 * i0 + 131; i4 += 1) {
                            B[1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[2][i3][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[0][i3][i4])) + (SCALAR_VAL(0.125) * ((A[1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3][i4 - 1]))) + A[1][i3][i4]);
                          }
                        }
                      }
                    }
                    for (int i2 = max(128 * ii0 + 64 * ii2 - 2 * i0 + 7, -ii0 + (-3 * ii0 + i0 + 28) / 61 + 1); i2 <= min(64 * ii2 + 1, 128 * ii0 + 64 * ii2 - 2 * i0 + 66); i2 += 1) {
                      if (64 * ii2 >= i2 + 1) {
                        for (int i3 = max(1, 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 5); i3 <= min(128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 + 68, 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 68); i3 += 1) {
                          for (int i4 = max(max(1, 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 5), 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 6); i4 <= 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 69; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          for (int i4 = 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 70; i4 <= 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 68; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                      for (int i3 = max(max(1, 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 + 69), 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 5); i3 <= 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 68; i3 += 1) {
                        if (ii4 == 1) {
                          for (int i4 = 128 * ii0 + 64 * ii2 - 2 * i0 - i2 + 69; i4 <= 128 * ii0 + 64 * ii2 - 2 * i0 - i2 + 70; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                        for (int i4 = max(1, 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 7); i4 <= 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 69; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 70; i4 <= 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 68; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                      if (i2 >= 64 * ii2) {
                        if (64 * ii0 + 32 * ii3 + 1 >= i0 && i2 == 64 * ii2 + 1) {
                          for (int i4 = 64 * ii4 + 1; i4 <= 64 * ii4 + 64; i4 += 1) {
                            B[64 * ii2 + 1][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4])) + A[64 * ii2][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][128 * ii0 + 64 * ii3 - 2 * i0 + 5][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4])) + A[64 * ii2 + 1][128 * ii0 + 64 * ii3 - 2 * i0 + 3][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4])) + A[64 * ii2 + 1][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4 - 1]))) + A[64 * ii2 + 1][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4]);
                          }
                          if (ii2 == 0 && ii4 == 0) {
                            for (int i3 = 128 * ii0 + 64 * ii3 - 2 * i0 + 5; i3 < 64 * ii3; i3 += 1) {
                              for (int i4 = 1; i4 <= 128 * ii0 + 64 * ii3 - 2 * i0 - i3 + 68; i4 += 1) {
                                B[1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[2][i3][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[0][i3][i4])) + (SCALAR_VAL(0.125) * ((A[1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3][i4 - 1]))) + A[1][i3][i4]);
                              }
                            }
                          }
                          if (ii2 == 0) {
                            for (int i3 = 128 * ii0 + 64 * ii3 - 2 * i0 + 5; i3 <= min(64 * ii3 - 1, 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 + 4); i3 += 1) {
                              for (int i4 = 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 5; i4 <= 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 68; i4 += 1) {
                                B[1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[2][i3][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[0][i3][i4])) + (SCALAR_VAL(0.125) * ((A[1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3][i4 - 1]))) + A[1][i3][i4]);
                              }
                            }
                          }
                        }
                        if (ii2 >= 1 && ii3 >= 1) {
                          for (int i3 = 128 * ii0 + 64 * ii3 - 2 * i0 + 5; i3 <= 64 * ii2 + 64 * ii3 - i2; i3 += 1) {
                            if (i2 == 64 * ii2 + 1 && 128 * ii0 + 64 * ii3 + 64 * ii4 + 4 >= 2 * i0 + i3) {
                              B[64 * ii2 + 1][i3][128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 5] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][i3][128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 5] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][i3][128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 5])) + A[64 * ii2][i3][128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 5])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][i3 + 1][128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 5] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][i3][128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 5])) + A[64 * ii2 + 1][i3 - 1][128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 5]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][i3][128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 6] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][i3][128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 5])) + A[64 * ii2 + 1][i3][128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 4]))) + A[64 * ii2 + 1][i3][128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 5]);
                            } else if (64 * ii0 + 32 * ii4 + 2 >= i0 && i2 == 64 * ii2 && i3 == 64 * ii3) {
                              B[64 * ii2][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 6] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 6] - (SCALAR_VAL(2.0) * A[64 * ii2][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 6])) + A[64 * ii2 - 1][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 6])) + (SCALAR_VAL(0.125) * ((A[64 * ii2][64 * ii3 + 1][128 * ii0 + 64 * ii4 - 2 * i0 + 6] - (SCALAR_VAL(2.0) * A[64 * ii2][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 6])) + A[64 * ii2][64 * ii3 - 1][128 * ii0 + 64 * ii4 - 2 * i0 + 6]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 7] - (SCALAR_VAL(2.0) * A[64 * ii2][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 6])) + A[64 * ii2][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 5]))) + A[64 * ii2][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 6]);
                            }
                            for (int i4 = max(max(1, 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 7), 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 6); i4 <= 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 69; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                        }
                        for (int i3 = max(1, 64 * ii2 + 64 * ii3 - i2 + 1); i3 <= 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 68; i3 += 1) {
                          for (int i4 = max(max(1, 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 5), 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 5); i4 <= 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 6; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (i2 >= 2) {
                            for (int i4 = max(1, 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 7); i4 <= 128 * ii3 + 64 * ii4 - 2 * i3; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                          for (int i4 = max(max(max(1, 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 7), 128 * ii0 + 64 * ii4 - 2 * i0 - 62 * i2 + 130), 128 * ii3 + 64 * ii4 - 2 * i3 + 1); i4 <= 64 * ii2 + 64 * ii4 - i2; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (ii2 == 0 && i2 == 1) {
                            for (int i4 = max(1, 128 * ii0 + 64 * ii4 - 2 * i0 + 6); i4 <= 128 * ii3 + 64 * ii4 - 2 * i3; i4 += 1) {
                              B[1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[2][i3][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[0][i3][i4])) + (SCALAR_VAL(0.125) * ((A[1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3][i4 - 1]))) + A[1][i3][i4]);
                            }
                            for (int i4 = max(max(1, 128 * ii0 + 64 * ii4 - 2 * i0 + 6), 128 * ii3 + 64 * ii4 - 2 * i3 + 1); i4 < 64 * ii4; i4 += 1) {
                              B[1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[2][i3][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[0][i3][i4])) + (SCALAR_VAL(0.125) * ((A[1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3][i4 - 1]))) + A[1][i3][i4]);
                            }
                          }
                          for (int i4 = max(max(1, 64 * ii2 + 64 * ii4 - i2 + 1), 128 * ii3 + 64 * ii4 - 2 * i3 + 1); i4 <= 128 * ii0 + 64 * ii4 - 2 * i0 + 67; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (i2 == 64 * ii2 && i3 >= 64 * ii3 + 2) {
                            B[64 * ii2][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 68] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 68] - (SCALAR_VAL(2.0) * A[64 * ii2][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 68])) + A[64 * ii2 - 1][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 68])) + (SCALAR_VAL(0.125) * ((A[64 * ii2][i3 + 1][128 * ii0 + 64 * ii4 - 2 * i0 + 68] - (SCALAR_VAL(2.0) * A[64 * ii2][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 68])) + A[64 * ii2][i3 - 1][128 * ii0 + 64 * ii4 - 2 * i0 + 68]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 69] - (SCALAR_VAL(2.0) * A[64 * ii2][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 68])) + A[64 * ii2][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 67]))) + A[64 * ii2][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 68]);
                          } else if (i2 >= 2 && i2 + i3 == 64 * ii2 + 64 * ii3 + 1) {
                            B[i2][64 * ii2 + 64 * ii3 - i2 + 1][128 * ii0 + 64 * ii4 - 2 * i0 + 68] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii2 + 64 * ii3 - i2 + 1][128 * ii0 + 64 * ii4 - 2 * i0 + 68] - (SCALAR_VAL(2.0) * A[i2][64 * ii2 + 64 * ii3 - i2 + 1][128 * ii0 + 64 * ii4 - 2 * i0 + 68])) + A[i2 - 1][64 * ii2 + 64 * ii3 - i2 + 1][128 * ii0 + 64 * ii4 - 2 * i0 + 68])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii2 + 64 * ii3 - i2 + 2][128 * ii0 + 64 * ii4 - 2 * i0 + 68] - (SCALAR_VAL(2.0) * A[i2][64 * ii2 + 64 * ii3 - i2 + 1][128 * ii0 + 64 * ii4 - 2 * i0 + 68])) + A[i2][64 * ii2 + 64 * ii3 - i2][128 * ii0 + 64 * ii4 - 2 * i0 + 68]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii2 + 64 * ii3 - i2 + 1][128 * ii0 + 64 * ii4 - 2 * i0 + 69] - (SCALAR_VAL(2.0) * A[i2][64 * ii2 + 64 * ii3 - i2 + 1][128 * ii0 + 64 * ii4 - 2 * i0 + 68])) + A[i2][64 * ii2 + 64 * ii3 - i2 + 1][128 * ii0 + 64 * ii4 - 2 * i0 + 67]))) + A[i2][64 * ii2 + 64 * ii3 - i2 + 1][128 * ii0 + 64 * ii4 - 2 * i0 + 68]);
                          } else if (ii2 == 0 && i3 == 64 * ii3) {
                            B[1][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 68] = ((((SCALAR_VAL(0.125) * ((A[2][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 68] - (SCALAR_VAL(2.0) * A[1][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 68])) + A[0][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 68])) + (SCALAR_VAL(0.125) * ((A[1][64 * ii3 + 1][128 * ii0 + 64 * ii4 - 2 * i0 + 68] - (SCALAR_VAL(2.0) * A[1][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 68])) + A[1][64 * ii3 - 1][128 * ii0 + 64 * ii4 - 2 * i0 + 68]))) + (SCALAR_VAL(0.125) * ((A[1][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 69] - (SCALAR_VAL(2.0) * A[1][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 68])) + A[1][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 67]))) + A[1][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 68]);
                          }
                        }
                      }
                    }
                    if (ii2 >= 1 && i0 >= 64 * ii0 + 33) {
                      for (int i3 = max(1, 64 * ii3 - 62); i3 <= 64 * ii3 + 1; i3 += 1) {
                        if (i0 == 64 * ii0 + 33 && 64 * ii3 + 64 * ii4 >= i3 + 62) {
                          B[64 * ii2 + 1][i3][64 * ii3 + 64 * ii4 - i3 - 61] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][i3][64 * ii3 + 64 * ii4 - i3 - 61] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][i3][64 * ii3 + 64 * ii4 - i3 - 61])) + A[64 * ii2][i3][64 * ii3 + 64 * ii4 - i3 - 61])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][i3 + 1][64 * ii3 + 64 * ii4 - i3 - 61] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][i3][64 * ii3 + 64 * ii4 - i3 - 61])) + A[64 * ii2 + 1][i3 - 1][64 * ii3 + 64 * ii4 - i3 - 61]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][i3][64 * ii3 + 64 * ii4 - i3 - 60] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][i3][64 * ii3 + 64 * ii4 - i3 - 61])) + A[64 * ii2 + 1][i3][64 * ii3 + 64 * ii4 - i3 - 62]))) + A[64 * ii2 + 1][i3][64 * ii3 + 64 * ii4 - i3 - 61]);
                          if (i3 + 62 == 64 * ii3) {
                            for (int i4 = 64 * ii4 + 2; i4 <= 64 * ii4 + 64; i4 += 1) {
                              B[64 * ii2 + 1][64 * ii3 - 62][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][64 * ii3 - 62][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][64 * ii3 - 62][i4])) + A[64 * ii2][64 * ii3 - 62][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][64 * ii3 - 61][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][64 * ii3 - 62][i4])) + A[64 * ii2 + 1][64 * ii3 - 63][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][64 * ii3 - 62][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][64 * ii3 - 62][i4])) + A[64 * ii2 + 1][64 * ii3 - 62][i4 - 1]))) + A[64 * ii2 + 1][64 * ii3 - 62][i4]);
                            }
                          }
                        }
                        if (2 * i0 + i3 >= 128 * ii0 + 64 * ii3 + 5) {
                          for (int i4 = max(max(1, 64 * ii3 + 64 * ii4 - i3 - 61), 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 6); i4 <= 64 * ii3 + 64 * ii4 - i3 + 2; i4 += 1) {
                            B[128 * ii0 + 64 * ii2 - 2 * i0 + 67][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 68][i3][i4] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 67][i3][i4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 66][i3][i4])) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 67][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 67][i3][i4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 67][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 67][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 67][i3][i4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 67][i3][i4 - 1]))) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 67][i3][i4]);
                          }
                        }
                      }
                    } else if (i0 == 64 * ii0 + 2) {
                      for (int i2 = max(1, 64 * ii2); i2 <= 64 * ii2 + 2; i2 += 1) {
                        if (64 * ii2 + 1 >= i2) {
                          for (int i3 = max(1, 64 * ii2 + 64 * ii3 - i2 + 1); i3 <= 64 * ii2 + 64 * ii3 - i2 + 64; i3 += 1) {
                            for (int i4 = max(max(1, 64 * ii2 + 64 * ii4 - i2 + 1), 128 * ii3 + 64 * ii4 - 2 * i3 + 1); i4 <= 64 * ii2 + 64 * ii4 - i2 + 64; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            if (i2 == 64 * ii2 + 1 && i3 == 64 * ii3) {
                              B[64 * ii2 + 1][64 * ii3][64 * ii4 + 64] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][64 * ii3][64 * ii4 + 64] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][64 * ii3][64 * ii4 + 64])) + A[64 * ii2][64 * ii3][64 * ii4 + 64])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][64 * ii3 + 1][64 * ii4 + 64] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][64 * ii3][64 * ii4 + 64])) + A[64 * ii2 + 1][64 * ii3 - 1][64 * ii4 + 64]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][64 * ii3][64 * ii4 + 65] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][64 * ii3][64 * ii4 + 64])) + A[64 * ii2 + 1][64 * ii3][64 * ii4 + 63]))) + A[64 * ii2 + 1][64 * ii3][64 * ii4 + 64]);
                            }
                          }
                        } else {
                          if (ii3 >= 1) {
                            for (int i4 = 64 * ii4 + 1; i4 <= 64 * ii4 + 64; i4 += 1) {
                              B[64 * ii2 + 2][64 * ii3][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 3][64 * ii3][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][64 * ii3][i4])) + A[64 * ii2 + 1][64 * ii3][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][64 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][64 * ii3][i4])) + A[64 * ii2 + 2][64 * ii3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][64 * ii3][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][64 * ii3][i4])) + A[64 * ii2 + 2][64 * ii3][i4 - 1]))) + A[64 * ii2 + 2][64 * ii3][i4]);
                            }
                          }
                          for (int i3 = 64 * ii3 + 1; i3 <= 64 * ii3 + 63; i3 += 1) {
                            for (int i4 = max(1, 64 * ii4); i4 <= 64 * ii4 + 63; i4 += 1) {
                              B[64 * ii2 + 2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 3][i3][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][i3][i4])) + A[64 * ii2 + 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][i3][i4])) + A[64 * ii2 + 2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][i3][i4])) + A[64 * ii2 + 2][i3][i4 - 1]))) + A[64 * ii2 + 2][i3][i4]);
                            }
                          }
                        }
                      }
                    }
                    for (int i2 = max(64 * ii2 + 2, 128 * ii0 + 64 * ii2 - 2 * i0 + 7); i2 <= 128 * ii0 + 64 * ii2 - 2 * i0 + 67; i2 += 1) {
                      if (i0 >= 64 * ii0 + 3 && 64 * ii0 + 32 * ii3 + 1 >= i0) {
                        for (int i4 = 64 * ii4 + 1; i4 <= 64 * ii4 + 64; i4 += 1) {
                          B[i2][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4] - (SCALAR_VAL(2.0) * A[i2][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4])) + A[i2 - 1][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4])) + (SCALAR_VAL(0.125) * ((A[i2][128 * ii0 + 64 * ii3 - 2 * i0 + 5][i4] - (SCALAR_VAL(2.0) * A[i2][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4])) + A[i2][128 * ii0 + 64 * ii3 - 2 * i0 + 3][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4])) + A[i2][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4 - 1]))) + A[i2][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4]);
                        }
                      } else if (ii3 >= 1 && i0 == 64 * ii0 + 2) {
                        for (int i4 = 64 * ii4 + 1; i4 <= 64 * ii4 + 64; i4 += 1) {
                          B[i2][64 * ii3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii3][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii3][i4])) + A[i2 - 1][64 * ii3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii3][i4])) + A[i2][64 * ii3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii3][i4])) + A[i2][64 * ii3][i4 - 1]))) + A[i2][64 * ii3][i4]);
                        }
                      }
                      for (int i3 = max(1, 128 * ii0 + 64 * ii3 - 2 * i0 + 5); i3 <= 128 * ii0 + 64 * ii3 - 2 * i0 + 67; i3 += 1) {
                        if (i0 >= 64 * ii0 + 3 && 64 * ii0 + 32 * ii4 + 1 >= i0 && i3 >= 64 * ii3 + 1) {
                          B[i2][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 4] - (SCALAR_VAL(2.0) * A[i2][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 4])) + A[i2 - 1][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][128 * ii0 + 64 * ii4 - 2 * i0 + 4] - (SCALAR_VAL(2.0) * A[i2][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 4])) + A[i2][i3 - 1][128 * ii0 + 64 * ii4 - 2 * i0 + 4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 5] - (SCALAR_VAL(2.0) * A[i2][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 4])) + A[i2][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 3]))) + A[i2][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 4]);
                        }
                        if (i3 >= 64 * ii3 && 128 * ii0 + 64 * ii2 + 64 * ii3 + 68 >= 2 * i0 + i2 + i3) {
                          for (int i4 = max(1, 128 * ii0 + 64 * ii4 - 2 * i0 + 5); i4 <= 128 * ii3 + 64 * ii4 - 2 * i3; i4 += 1) {
                            if (2 * i0 + i3 + i4 >= 128 * ii0 + 64 * ii3 + 64 * ii4 + 6 || 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                        }
                        if (128 * ii0 + 64 * ii2 + 64 * ii3 + 68 >= 2 * i0 + i2 + i3) {
                          for (int i4 = max(max(1, 128 * ii0 + 64 * ii4 - 2 * i0 + 5), 128 * ii3 + 64 * ii4 - 2 * i3 + 1); i4 < 64 * ii4; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (i3 >= 64 * ii3) {
                            for (int i4 = max(max(1, 64 * ii4), 128 * ii3 + 64 * ii4 - 2 * i3 + 1); i4 <= min(128 * ii3 + 64 * ii4 - 2 * i3 + 64, 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 69); i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            if (i0 == 64 * ii0 + 2) {
                              for (int i4 = 64 * ii2 + 64 * ii3 + 64 * ii4 - i2 - i3 + 66; i4 <= 128 * ii3 + 64 * ii4 - 2 * i3 + 64; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                              for (int i4 = max(max(1, 64 * ii4), 128 * ii3 + 64 * ii4 - 2 * i3 + 65); i4 <= min(64 * ii4 + 63, 64 * ii3 + 64 * ii4 - i3 + 65); i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                          }
                          if (i0 >= 64 * ii0 + 3) {
                            for (int i4 = max(max(1, 64 * ii4), 128 * ii3 + 64 * ii4 - 2 * i3 + 65); i4 <= 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 69; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            if (64 * ii3 >= i3 + 1) {
                              for (int i4 = max(1, 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 5); i4 <= 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 69; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                            for (int i4 = 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 70; i4 <= min(128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 68, 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 69); i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            for (int i4 = max(128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 69, 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 70); i4 <= min(128 * ii0 + 64 * ii4 - 2 * i0 + 67, 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 69); i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                        }
                        for (int i4 = max(128 * ii0 + 64 * ii4 - 2 * i0 + 68, 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 70); i4 <= 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 68; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        if (2 * i0 + i2 + i3 >= 128 * ii0 + 64 * ii2 + 64 * ii3 + 69) {
                          for (int i4 = max(1, 128 * ii0 + 64 * ii4 - 2 * i0 + 5); i4 < 64 * ii4; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          for (int i4 = max(1, 64 * ii4); i4 <= 128 * ii0 + 64 * ii4 - 2 * i0 + 67; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        } else {
                          for (int i4 = 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 70; i4 <= 128 * ii0 + 64 * ii4 - 2 * i0 + 67; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                    }
                    if (ii2 == 0 && i0 == 64 * ii0 + 33) {
                      for (int i3 = 64 * ii3 - 62; i3 <= 64 * ii3 + 1; i3 += 1) {
                        for (int i4 = max(1, 64 * ii3 + 64 * ii4 - i3 - 61); i4 <= 64 * ii3 + 64 * ii4 - i3 + 2; i4 += 1) {
                          B[1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[2][i3][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[0][i3][i4])) + (SCALAR_VAL(0.125) * ((A[1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3][i4 - 1]))) + A[1][i3][i4]);
                        }
                      }
                    }
                    if (i0 >= 64 * ii0 + 33) {
                      for (int i2 = max(1, 128 * ii0 + 64 * ii2 - 2 * i0 + 3); i2 <= 128 * ii0 + 64 * ii2 - 2 * i0 + 66; i2 += 1) {
                        for (int i3 = max(1, 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 4); i3 <= 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 67; i3 += 1) {
                          for (int i4 = max(max(1, 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 4), 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 5); i4 <= 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 68; i4 += 1) {
                            A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                          }
                          for (int i4 = 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 69; i4 <= 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 67; i4 += 1) {
                            A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                          }
                        }
                      }
                    } else {
                      for (int i2 = max(1, 128 * ii0 + 64 * ii2 - 2 * i0 + 3); i2 <= 64 * ii2; i2 += 1) {
                        for (int i3 = max(1, 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 4); i3 <= 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 67; i3 += 1) {
                          for (int i4 = max(max(1, 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 4), 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 5); i4 <= 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 68; i4 += 1) {
                            A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                          }
                          for (int i4 = 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 69; i4 <= 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 67; i4 += 1) {
                            A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                          }
                        }
                      }
                      for (int i2 = 64 * ii2 + 1; i2 <= 128 * ii0 + 64 * ii2 - 2 * i0 + 66; i2 += 1) {
                        for (int i3 = max(1, 128 * ii0 + 64 * ii3 - 2 * i0 + 3); i3 <= 128 * ii0 + 64 * ii3 - 2 * i0 + 66; i3 += 1) {
                          if (i3 >= 64 * ii3 + 1) {
                            for (int i4 = max(1, 128 * ii0 + 64 * ii4 - 2 * i0 + 3); i4 <= 128 * ii0 + 64 * ii4 - 2 * i0 + 66; i4 += 1) {
                              A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                            }
                          } else {
                            for (int i4 = max(1, 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 4); i4 <= 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 67; i4 += 1) {
                              A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                            }
                          }
                        }
                      }
                    }
                  }
                }
              } else if (ii4 >= 1) {
                if (_PB_N >= 64 * ii4 + 67) {
                  for (int i2 = 64 * ii2 + 1; i2 < _PB_N - 1; i2 += 1) {
                    for (int i3 = 64 * ii3 + 1; i3 <= 64 * ii3 + 64; i3 += 1) {
                      for (int i4 = 64 * ii4 + 1; i4 <= 64 * ii4 + 64; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                  }
                } else {
                  for (int i2 = 64 * ii2 + 1; i2 < _PB_N - 1; i2 += 1) {
                    for (int i3 = 64 * ii3 + 1; i3 <= 64 * ii3 + 64; i3 += 1) {
                      for (int i4 = 64 * ii2 + 1; i4 < _PB_N - 1; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                  }
                }
                for (int i0 = 64 * ii0 + 2; i0 <= min(_PB_TSTEPS, 64 * ii0 + 32); i0 += 1) {
                  if (_PB_N >= 64 * ii4 + 67) {
                    if (i0 >= 64 * ii0 + 3) {
                      for (int i2 = 128 * ii0 + 64 * ii2 - 2 * i0 + 4; i2 <= 128 * ii0 + 64 * ii2 - 2 * i0 + 6; i2 += 1) {
                        for (int i3 = max(1, 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 5); i3 <= 64 * ii3 + 1; i3 += 1) {
                          for (int i4 = 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 6; i4 <= 64 * ii4 + 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (i3 == 64 * ii3 + 1) {
                            for (int i4 = 64 * ii4 + 2; i4 <= 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 68; i4 += 1) {
                              B[i2][64 * ii3 + 1][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 1][i4])) + A[i2 - 1][64 * ii3 + 1][i4])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 2][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 1][i4])) + A[i2][64 * ii3][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 1][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 1][i4])) + A[i2][64 * ii3 + 1][i4 - 1]))) + A[i2][64 * ii3 + 1][i4]);
                            }
                          } else {
                            for (int i4 = 64 * ii4 + 2; i4 <= 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 69; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                        }
                        if (2 * i0 + i2 == 128 * ii0 + 64 * ii2 + 4) {
                          for (int i4 = 64 * ii4 + 1; i4 <= 64 * ii4 + 64; i4 += 1) {
                            B[128 * ii0 + 64 * ii2 - 2 * i0 + 4][64 * ii3 + 2][i4] = ((((SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][64 * ii3 + 2][i4] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][64 * ii3 + 2][i4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 3][64 * ii3 + 2][i4])) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][64 * ii3 + 3][i4] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][64 * ii3 + 2][i4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][64 * ii3 + 1][i4]))) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][64 * ii3 + 2][i4 + 1] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][64 * ii3 + 2][i4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][64 * ii3 + 2][i4 - 1]))) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][64 * ii3 + 2][i4]);
                          }
                        }
                        for (int i3 = max(64 * ii3 + 2, 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 7); i3 <= 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 68; i3 += 1) {
                          for (int i4 = 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 5; i4 <= 64 * ii4; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (64 * ii2 >= i2 + 1 && i3 == 64 * ii3 + 2) {
                            for (int i4 = 64 * ii4 + 1; i4 <= 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 68; i4 += 1) {
                              B[i2][64 * ii3 + 2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii3 + 2][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 2][i4])) + A[i2 - 1][64 * ii3 + 2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 3][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 2][i4])) + A[i2][64 * ii3 + 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 2][i4])) + A[i2][64 * ii3 + 2][i4 - 1]))) + A[i2][64 * ii3 + 2][i4]);
                            }
                          } else if (i0 == 64 * ii0 + 3 && i2 == 64 * ii2 && i3 == 64 * ii3 + 2) {
                            for (int i4 = 64 * ii4 + 1; i4 <= 64 * ii4 + 62; i4 += 1) {
                              B[64 * ii2][64 * ii3 + 2][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][64 * ii3 + 2][i4] - (SCALAR_VAL(2.0) * A[64 * ii2][64 * ii3 + 2][i4])) + A[64 * ii2 - 1][64 * ii3 + 2][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii2][64 * ii3 + 3][i4] - (SCALAR_VAL(2.0) * A[64 * ii2][64 * ii3 + 2][i4])) + A[64 * ii2][64 * ii3 + 1][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2][64 * ii3 + 2][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii2][64 * ii3 + 2][i4])) + A[64 * ii2][64 * ii3 + 2][i4 - 1]))) + A[64 * ii2][64 * ii3 + 2][i4]);
                            }
                          } else if (2 * i0 + i2 == 128 * ii0 + 64 * ii2 + 6) {
                            for (int i4 = 64 * ii4 + 1; i4 <= 64 * ii4 + 62; i4 += 1) {
                              B[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 7][i3][i4] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][i4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][i4])) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][i4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][i4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][i4 - 1]))) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][i4]);
                            }
                          } else {
                            for (int i4 = 64 * ii4 + 1; i4 <= 256 * ii0 + 128 * ii2 + 64 * ii4 - 4 * i0 - 2 * i2 + 72; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            if (2 * i0 + i2 == 128 * ii0 + 64 * ii2 + 5) {
                              B[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][64 * ii4 + 63] = ((((SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][64 * ii4 + 63] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][64 * ii4 + 63])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][i3][64 * ii4 + 63])) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3 + 1][64 * ii4 + 63] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][64 * ii4 + 63])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3 - 1][64 * ii4 + 63]))) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][64 * ii4 + 64] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][64 * ii4 + 63])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][64 * ii4 + 62]))) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][64 * ii4 + 63]);
                            }
                          }
                        }
                      }
                    }
                    for (int i2 = 128 * ii0 + 64 * ii2 - 2 * i0 + 7; i2 < 64 * ii2; i2 += 1) {
                      for (int i3 = max(1, 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 5); i3 <= 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 68; i3 += 1) {
                        if (i3 >= 64 * ii3 + 2) {
                          for (int i4 = 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 5; i4 <= 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 68; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        } else {
                          for (int i4 = 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 6; i4 <= 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 69; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                    }
                    if (i0 == 64 * ii0 + 2) {
                      for (int i2 = 64 * ii2; i2 <= min(_PB_N - 2, 64 * ii2 + 2); i2 += 1) {
                        if (64 * ii2 + 1 >= i2) {
                          for (int i3 = max(1, 64 * ii2 + 64 * ii3 - i2 + 1); i3 <= 64 * ii2 + 64 * ii3 - i2 + 64; i3 += 1) {
                            for (int i4 = max(64 * ii2 + 64 * ii4 - i2 + 1, 128 * ii3 + 64 * ii4 - 2 * i3 + 1); i4 <= 64 * ii2 + 64 * ii4 - i2 + 64; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            if (i2 == 64 * ii2 + 1 && i3 == 64 * ii3) {
                              B[64 * ii2 + 1][64 * ii3][64 * ii4 + 64] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][64 * ii3][64 * ii4 + 64] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][64 * ii3][64 * ii4 + 64])) + A[64 * ii2][64 * ii3][64 * ii4 + 64])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][64 * ii3 + 1][64 * ii4 + 64] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][64 * ii3][64 * ii4 + 64])) + A[64 * ii2 + 1][64 * ii3 - 1][64 * ii4 + 64]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][64 * ii3][64 * ii4 + 65] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][64 * ii3][64 * ii4 + 64])) + A[64 * ii2 + 1][64 * ii3][64 * ii4 + 63]))) + A[64 * ii2 + 1][64 * ii3][64 * ii4 + 64]);
                            }
                          }
                        } else {
                          if (ii3 >= 1) {
                            for (int i4 = 64 * ii4 + 1; i4 <= 64 * ii4 + 64; i4 += 1) {
                              B[64 * ii2 + 2][64 * ii3][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 3][64 * ii3][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][64 * ii3][i4])) + A[64 * ii2 + 1][64 * ii3][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][64 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][64 * ii3][i4])) + A[64 * ii2 + 2][64 * ii3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][64 * ii3][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][64 * ii3][i4])) + A[64 * ii2 + 2][64 * ii3][i4 - 1]))) + A[64 * ii2 + 2][64 * ii3][i4]);
                            }
                          }
                          for (int i3 = 64 * ii3 + 1; i3 <= 64 * ii3 + 63; i3 += 1) {
                            for (int i4 = 64 * ii4; i4 <= 64 * ii4 + 63; i4 += 1) {
                              B[64 * ii2 + 2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 3][i3][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][i3][i4])) + A[64 * ii2 + 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][i3][i4])) + A[64 * ii2 + 2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][i3][i4])) + A[64 * ii2 + 2][i3][i4 - 1]))) + A[64 * ii2 + 2][i3][i4]);
                            }
                          }
                        }
                      }
                    }
                    for (int i2 = max(64 * ii2, 128 * ii0 + 64 * ii2 - 2 * i0 + 7); i2 < _PB_N - 1; i2 += 1) {
                      if (64 * ii2 + 1 >= i2) {
                        if (i2 == 64 * ii2) {
                          for (int i3 = max(1, 128 * ii0 + 64 * ii3 - 2 * i0 + 5); i3 <= 64 * ii3; i3 += 1) {
                            for (int i4 = 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 6; i4 <= 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 69; i4 += 1) {
                              B[64 * ii2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[64 * ii2][i3][i4])) + A[64 * ii2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[64 * ii2][i3][i4])) + A[64 * ii2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii2][i3][i4])) + A[64 * ii2][i3][i4 - 1]))) + A[64 * ii2][i3][i4]);
                            }
                          }
                        } else {
                          for (int i3 = max(1, 128 * ii0 + 64 * ii3 - 2 * i0 + 4); i3 < 64 * ii3; i3 += 1) {
                            for (int i4 = 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 5; i4 <= 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 68; i4 += 1) {
                              B[64 * ii2 + 1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][i3][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][i3][i4])) + A[64 * ii2][i3][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][i3][i4])) + A[64 * ii2 + 1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][i3][i4])) + A[64 * ii2 + 1][i3][i4 - 1]))) + A[64 * ii2 + 1][i3][i4]);
                            }
                          }
                        }
                        for (int i3 = max(1, 64 * ii2 + 64 * ii3 - i2 + 1); i3 <= 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 68; i3 += 1) {
                          if (i3 >= 64 * ii3 + 2) {
                            for (int i4 = 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 5; i4 <= 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 68; i4 += 1) {
                              if (i2 + i4 >= 64 * ii2 + 64 * ii4 + 1 || 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                          } else {
                            if (i2 == 64 * ii2 + 1) {
                              B[64 * ii2 + 1][i3][128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 5] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][i3][128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 5] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][i3][128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 5])) + A[64 * ii2][i3][128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 5])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][i3 + 1][128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 5] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][i3][128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 5])) + A[64 * ii2 + 1][i3 - 1][128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 5]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][i3][128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 6] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][i3][128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 5])) + A[64 * ii2 + 1][i3][128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 4]))) + A[64 * ii2 + 1][i3][128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 5]);
                            }
                            for (int i4 = max(128 * ii0 - 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 + i2 - i3 + 5, 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 6); i4 <= 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 69; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                        }
                      } else if (i0 >= 64 * ii0 + 3) {
                        for (int i3 = max(1, 128 * ii0 + 64 * ii3 - 2 * i0 + 4); i3 <= 128 * ii0 + 64 * ii3 - 2 * i0 + 5; i3 += 1) {
                          for (int i4 = 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 5; i4 <= 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 68; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                        for (int i3 = max(1, 128 * ii0 + 64 * ii3 - 2 * i0 + 6); i3 <= min(64 * ii3 + 2, -64 * ii2 + 64 * ii3 + i2 - 1); i3 += 1) {
                          for (int i4 = max(128 * ii0 + 64 * ii4 - 2 * i0 + 4, 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 5); i4 <= min(min(64 * ii4 + 1, 128 * ii3 + 64 * ii4 - 2 * i3), 128 * ii0 - 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 + i2 - i3 + 4); i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          for (int i4 = max(128 * ii0 + 64 * ii4 - 2 * i0 + 4, 128 * ii3 + 64 * ii4 - 2 * i3 + 1); i4 <= min(64 * ii4 - 1, 128 * ii0 - 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 + i2 - i3 + 4); i4 += 1) {
                            if (2 * i0 + i3 + i4 >= 128 * ii0 + 64 * ii3 + 64 * ii4 + 7 || 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                          if (2 * i0 + i2 >= 128 * ii0 + 64 * ii2 + 67 && i3 == 64 * ii3 + 2) {
                            for (int i4 = 128 * ii0 - 64 * ii2 + 64 * ii4 - 2 * i0 + i2 + 3; i4 <= min(64 * ii4 - 1, 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 68); i4 += 1) {
                              B[i2][64 * ii3 + 2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii3 + 2][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 2][i4])) + A[i2 - 1][64 * ii3 + 2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 3][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 2][i4])) + A[i2][64 * ii3 + 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 2][i4])) + A[i2][64 * ii3 + 2][i4 - 1]))) + A[i2][64 * ii3 + 2][i4]);
                            }
                          } else if (128 * ii0 + 64 * ii2 + 64 * ii3 + 68 >= 2 * i0 + i2 + i3) {
                            for (int i4 = 128 * ii0 - 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 + i2 - i3 + 5; i4 <= min(64 * ii4 + 1, 128 * ii3 + 64 * ii4 - 2 * i3); i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            for (int i4 = max(128 * ii3 + 64 * ii4 - 2 * i3 + 1, 128 * ii0 - 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 + i2 - i3 + 5); i4 < 64 * ii4; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                          if (128 * ii0 + i2 + 2 >= 64 * ii2 + 2 * i0 && i3 == 64 * ii3 + 2) {
                            for (int i4 = 64 * ii4; i4 <= 64 * ii4 + 1; i4 += 1) {
                              B[i2][64 * ii3 + 2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii3 + 2][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 2][i4])) + A[i2 - 1][64 * ii3 + 2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 3][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 2][i4])) + A[i2][64 * ii3 + 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 2][i4])) + A[i2][64 * ii3 + 2][i4 - 1]))) + A[i2][64 * ii3 + 2][i4]);
                            }
                            for (int i4 = 64 * ii4 + 2; i4 <= min(64 * ii4 + 60, 128 * ii0 + 64 * ii4 - 2 * i0 + 67); i4 += 1) {
                              B[i2][64 * ii3 + 2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii3 + 2][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 2][i4])) + A[i2 - 1][64 * ii3 + 2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 3][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 2][i4])) + A[i2][64 * ii3 + 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 2][i4])) + A[i2][64 * ii3 + 2][i4 - 1]))) + A[i2][64 * ii3 + 2][i4]);
                            }
                            if (i0 == 64 * ii0 + 3) {
                              B[i2][64 * ii3 + 2][64 * ii4 + 61] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii3 + 2][64 * ii4 + 61] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 2][64 * ii4 + 61])) + A[i2 - 1][64 * ii3 + 2][64 * ii4 + 61])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 3][64 * ii4 + 61] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 2][64 * ii4 + 61])) + A[i2][64 * ii3 + 1][64 * ii4 + 61]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 2][64 * ii4 + 62] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 2][64 * ii4 + 61])) + A[i2][64 * ii3 + 2][64 * ii4 + 60]))) + A[i2][64 * ii3 + 2][64 * ii4 + 61]);
                            }
                          } else {
                            if (2 * i0 + i2 + i3 >= 128 * ii0 + 64 * ii2 + 64 * ii3 + 69 && 64 * ii3 + 1 >= i3) {
                              for (int i4 = 128 * ii0 - 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 + i2 - i3 + 5; i4 <= min(64 * ii4 + 1, 128 * ii3 + 64 * ii4 - 2 * i3); i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                              if (64 * ii2 + 2 * i0 >= 128 * ii0 + i2 + 5 && i3 == 64 * ii3 + 1) {
                                B[i2][64 * ii3 + 1][64 * ii4 - 1] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii3 + 1][64 * ii4 - 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 1][64 * ii4 - 1])) + A[i2 - 1][64 * ii3 + 1][64 * ii4 - 1])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 2][64 * ii4 - 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 1][64 * ii4 - 1])) + A[i2][64 * ii3][64 * ii4 - 1]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 1][64 * ii4] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 1][64 * ii4 - 1])) + A[i2][64 * ii3 + 1][64 * ii4 - 2]))) + A[i2][64 * ii3 + 1][64 * ii4 - 1]);
                              }
                            } else if (64 * ii2 + 2 * i0 >= 128 * ii0 + i2 + 3 && i3 == 64 * ii3 + 2) {
                              for (int i4 = 64 * ii4; i4 <= min(64 * ii4 + 1, 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 68); i4 += 1) {
                                B[i2][64 * ii3 + 2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii3 + 2][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 2][i4])) + A[i2 - 1][64 * ii3 + 2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 3][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 2][i4])) + A[i2][64 * ii3 + 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 2][i4])) + A[i2][64 * ii3 + 2][i4 - 1]))) + A[i2][64 * ii3 + 2][i4]);
                              }
                              for (int i4 = 64 * ii4 + 2; i4 <= 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 68; i4 += 1) {
                                B[i2][64 * ii3 + 2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii3 + 2][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 2][i4])) + A[i2 - 1][64 * ii3 + 2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 3][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 2][i4])) + A[i2][64 * ii3 + 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 2][i4])) + A[i2][64 * ii3 + 2][i4 - 1]))) + A[i2][64 * ii3 + 2][i4]);
                              }
                            }
                            if (64 * ii3 + 1 >= i3) {
                              for (int i4 = max(64 * ii4, 128 * ii3 + 64 * ii4 - 2 * i3 + 1); i4 <= 64 * ii4 + 1; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                              for (int i4 = 64 * ii4 + 2; i4 <= 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 68; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            } else {
                              for (int i4 = max(128 * ii0 - 64 * ii2 + 64 * ii4 - 2 * i0 + i2 + 3, 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 69); i4 <= 128 * ii0 + 64 * ii4 - 2 * i0 + 67; i4 += 1) {
                                B[i2][64 * ii3 + 2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii3 + 2][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 2][i4])) + A[i2 - 1][64 * ii3 + 2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 3][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 2][i4])) + A[i2][64 * ii3 + 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 2][i4])) + A[i2][64 * ii3 + 2][i4 - 1]))) + A[i2][64 * ii3 + 2][i4]);
                              }
                            }
                          }
                        }
                        if (i2 == 64 * ii2 + 2) {
                          for (int i4 = 128 * ii0 + 64 * ii4 - 2 * i0 + 4; i4 <= 128 * ii0 + 64 * ii4 - 2 * i0 + 67; i4 += 1) {
                            B[64 * ii2 + 2][64 * ii3 + 2][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 3][64 * ii3 + 2][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][64 * ii3 + 2][i4])) + A[64 * ii2 + 1][64 * ii3 + 2][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][64 * ii3 + 3][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][64 * ii3 + 2][i4])) + A[64 * ii2 + 2][64 * ii3 + 1][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][64 * ii3 + 2][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][64 * ii3 + 2][i4])) + A[64 * ii2 + 2][64 * ii3 + 2][i4 - 1]))) + A[64 * ii2 + 2][64 * ii3 + 2][i4]);
                          }
                        }
                        for (int i3 = 64 * ii3 + 3; i3 <= 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 68; i3 += 1) {
                          for (int i4 = 128 * ii0 + 64 * ii4 - 2 * i0 + 4; i4 <= min(64 * ii4, 128 * ii0 - 64 * ii2 + 64 * ii4 - 2 * i0 + i2 + 2); i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (128 * ii0 + i2 + 2 >= 64 * ii2 + 2 * i0) {
                            for (int i4 = 64 * ii4 + 1; i4 <= 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 68; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          } else {
                            for (int i4 = 128 * ii0 - 64 * ii2 + 64 * ii4 - 2 * i0 + i2 + 3; i4 <= 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 68; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                          for (int i4 = 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 69; i4 <= 128 * ii0 + 64 * ii4 - 2 * i0 + 67; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                        for (int i3 = max(64 * ii3 + 3, 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 69); i3 <= 128 * ii0 + 64 * ii3 - 2 * i0 + 67; i3 += 1) {
                          for (int i4 = 128 * ii0 + 64 * ii4 - 2 * i0 + 4; i4 <= min(64 * ii4, 128 * ii0 - 64 * ii2 + 64 * ii4 - 2 * i0 + i2 + 2); i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (128 * ii0 + i2 + 2 >= 64 * ii2 + 2 * i0) {
                            for (int i4 = 64 * ii4 + 1; i4 <= 128 * ii0 + 64 * ii4 - 2 * i0 + 67; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          } else {
                            for (int i4 = 128 * ii0 - 64 * ii2 + 64 * ii4 - 2 * i0 + i2 + 3; i4 <= 128 * ii0 + 64 * ii4 - 2 * i0 + 67; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                        }
                      } else {
                        for (int i3 = max(1, 64 * ii3); i3 <= 64 * ii2 + 64 * ii3 - i2 + 64; i3 += 1) {
                          for (int i4 = max(64 * ii4, 128 * ii3 + 64 * ii4 - 2 * i3 + 1); i4 <= 64 * ii4 + 63; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (i3 == 64 * ii3) {
                            B[i2][64 * ii3][64 * ii4 + 64] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii3][64 * ii4 + 64] - (SCALAR_VAL(2.0) * A[i2][64 * ii3][64 * ii4 + 64])) + A[i2 - 1][64 * ii3][64 * ii4 + 64])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 1][64 * ii4 + 64] - (SCALAR_VAL(2.0) * A[i2][64 * ii3][64 * ii4 + 64])) + A[i2][64 * ii3 - 1][64 * ii4 + 64]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3][64 * ii4 + 65] - (SCALAR_VAL(2.0) * A[i2][64 * ii3][64 * ii4 + 64])) + A[i2][64 * ii3][64 * ii4 + 63]))) + A[i2][64 * ii3][64 * ii4 + 64]);
                          }
                        }
                        for (int i3 = 64 * ii2 + 64 * ii3 - i2 + 65; i3 <= 64 * ii3 + 63; i3 += 1) {
                          for (int i4 = 64 * ii4; i4 <= 64 * ii4 + 63; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                    }
                  } else {
                    for (int i2 = 128 * ii0 + 64 * ii2 - 2 * i0 + 4; i2 <= 64 * ii2 + 1; i2 += 1) {
                      if (i0 >= 64 * ii0 + 3) {
                        for (int i3 = max(1, 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 5); i3 <= 64 * ii3 + 2; i3 += 1) {
                          if (i2 == 64 * ii2 + 1 && 64 * ii3 + 1 >= i3) {
                            B[64 * ii2 + 1][i3][128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i3 + 5] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][i3][128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i3 + 5] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][i3][128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i3 + 5])) + A[64 * ii2][i3][128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i3 + 5])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][i3 + 1][128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i3 + 5] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][i3][128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i3 + 5])) + A[64 * ii2 + 1][i3 - 1][128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i3 + 5]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][i3][128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i3 + 6] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][i3][128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i3 + 5])) + A[64 * ii2 + 1][i3][128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i3 + 4]))) + A[64 * ii2 + 1][i3][128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i3 + 5]);
                          }
                          if (64 * ii2 + 2 * i0 + i3 >= 128 * ii0 + 64 * ii3 + i2 + 4) {
                            for (int i4 = max(max(128 * ii0 + 128 * ii2 - 2 * i0 - i2 + 5, 128 * ii0 + 64 * ii3 - 2 * i0 + i2 - i3 + 5), 128 * ii0 + 128 * ii2 + 64 * ii3 - 2 * i0 - i2 - i3 + 6); i4 < _PB_N - 1; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          } else {
                            for (int i4 = 64 * ii2 + 2; i4 < _PB_N - 1; i4 += 1) {
                              B[64 * ii2 + 1][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4])) + A[64 * ii2][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][128 * ii0 + 64 * ii3 - 2 * i0 + 5][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4])) + A[64 * ii2 + 1][128 * ii0 + 64 * ii3 - 2 * i0 + 3][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4])) + A[64 * ii2 + 1][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4 - 1]))) + A[64 * ii2 + 1][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4]);
                            }
                          }
                        }
                        for (int i3 = 64 * ii3 + 3; i3 <= 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 68; i3 += 1) {
                          if (64 * ii3 + 63 >= i3) {
                            for (int i4 = max(-256 * ii0 - 64 * ii2 + 4 * i0 + 2 * i2 - 10, 256 * ii0 + 192 * ii2 - 4 * i0 - 2 * i2 + 9); i4 <= 64 * ii2 + 1; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            if (2 * i0 + i2 >= 128 * ii0 + 64 * ii2 + 6) {
                              for (int i4 = 128 * ii0 + 128 * ii2 - 2 * i0 - i2 + 5; i4 <= 64 * ii2 + 1; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                            for (int i4 = 64 * ii2 + 2; i4 < _PB_N - 1; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          } else {
                            for (int i4 = 64 * ii2 + 1; i4 < _PB_N - 1; i4 += 1) {
                              B[128 * ii0 + 64 * ii2 - 2 * i0 + 4][64 * ii3 + 64][i4] = ((((SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][64 * ii3 + 64][i4] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][64 * ii3 + 64][i4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 3][64 * ii3 + 64][i4])) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][64 * ii3 + 65][i4] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][64 * ii3 + 64][i4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][64 * ii3 + 63][i4]))) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][64 * ii3 + 64][i4 + 1] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][64 * ii3 + 64][i4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][64 * ii3 + 64][i4 - 1]))) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][64 * ii3 + 64][i4]);
                            }
                          }
                        }
                      } else {
                        for (int i3 = max(1, 64 * ii2 + 64 * ii3 - i2 + 1); i3 <= 64 * ii3 + 2; i3 += 1) {
                          for (int i4 = max(128 * ii2 - i2 + 1, 64 * ii2 + 32 * ii3 - i3 + i3 / 2 + 1); i4 < _PB_N - 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                        for (int i3 = 64 * ii3 + 3; i3 <= 64 * ii2 + 64 * ii3 - i2 + 64; i3 += 1) {
                          for (int i4 = 128 * ii2 - i2 + 1; i4 < _PB_N - 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                    }
                    for (int i2 = 64 * ii2 + 2; i2 < _PB_N - 1; i2 += 1) {
                      for (int i3 = max(1, 128 * ii0 + 64 * ii3 - 2 * i0 + 4); i3 < 64 * ii3; i3 += 1) {
                        for (int i4 = 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i3 + 5; i4 <= min(64 * ii2 + 1, 128 * ii0 + 64 * ii3 - 2 * i0 + i2 - i3 + 4); i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        if (64 * ii2 + 2 * i0 + i3 >= 128 * ii0 + 64 * ii3 + i2 + 4) {
                          for (int i4 = 128 * ii0 + 64 * ii3 - 2 * i0 + i2 - i3 + 5; i4 < _PB_N - 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        } else {
                          for (int i4 = 64 * ii2 + 2; i4 < _PB_N - 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                      for (int i3 = max(1, 64 * ii3); i3 <= 128 * ii0 + 64 * ii3 - 2 * i0 + 67; i3 += 1) {
                        if (64 * ii3 + 2 >= i3) {
                          for (int i4 = max(128 * ii0 + 64 * ii2 - 2 * i0 + 4, 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i3 + 5); i4 <= min(128 * ii0 + 64 * ii3 - 2 * i0 + i2 - i3 + 4, 64 * ii2 + 32 * ii3 - i3 + i3 / 2); i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        } else {
                          for (int i4 = 128 * ii0 + 64 * ii2 - 2 * i0 + 4; i4 <= min(64 * ii2 - 1, 128 * ii0 - 2 * i0 + i2 + 2); i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                        if (2 * i0 + i2 + i3 >= 128 * ii0 + 64 * ii2 + 64 * ii3 + 69) {
                          for (int i4 = max(128 * ii0 - 2 * i0 + i2 + 3, 128 * ii0 + 64 * ii3 - 2 * i0 + i2 - i3 + 5); i4 < 64 * ii2; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (64 * ii2 + 2 * i0 >= 128 * ii0 + i2 + 5 && i3 == 64 * ii3) {
                            B[i2][64 * ii3][64 * ii2] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii3][64 * ii2] - (SCALAR_VAL(2.0) * A[i2][64 * ii3][64 * ii2])) + A[i2 - 1][64 * ii3][64 * ii2])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 1][64 * ii2] - (SCALAR_VAL(2.0) * A[i2][64 * ii3][64 * ii2])) + A[i2][64 * ii3 - 1][64 * ii2]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3][64 * ii2 + 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii3][64 * ii2])) + A[i2][64 * ii3][64 * ii2 - 1]))) + A[i2][64 * ii3][64 * ii2]);
                          }
                        } else if (64 * ii3 + 2 >= i3) {
                          for (int i4 = 128 * ii0 + 64 * ii3 - 2 * i0 + i2 - i3 + 5; i4 <= 64 * ii2 + 32 * ii3 - i3 + i3 / 2; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        } else {
                          for (int i4 = 128 * ii0 - 2 * i0 + i2 + 3; i4 < 64 * ii2; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                        for (int i4 = max(64 * ii2, 64 * ii2 + 32 * ii3 - i3 + i3 / 2 + 1); i4 < _PB_N - 1; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    }
                  }
                  if (ii4 == ii2) {
                    for (int i2 = 128 * ii0 + 64 * ii2 - 2 * i0 + 3; i2 <= 64 * ii2; i2 += 1) {
                      for (int i3 = max(1, 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 4); i3 <= 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 67; i3 += 1) {
                        for (int i4 = max(128 * ii0 + 128 * ii2 - 2 * i0 - i2 + 4, 128 * ii0 + 128 * ii2 + 64 * ii3 - 2 * i0 - i2 - i3 + 5); i4 < _PB_N - 1; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                  } else {
                    for (int i2 = 128 * ii0 + 64 * ii2 - 2 * i0 + 3; i2 <= 64 * ii2; i2 += 1) {
                      for (int i3 = max(1, 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 4); i3 <= 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 67; i3 += 1) {
                        for (int i4 = max(128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 4, 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 5); i4 <= 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 68; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                        for (int i4 = 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 69; i4 <= 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 67; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                  }
                  for (int i2 = 64 * ii2 + 1; i2 < _PB_N - 1; i2 += 1) {
                    if (ii4 == ii2) {
                      for (int i3 = max(1, 128 * ii0 + 64 * ii3 - 2 * i0 + 3); i3 <= 64 * ii3; i3 += 1) {
                        for (int i4 = 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i3 + 4; i4 < _PB_N - 1; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    } else {
                      for (int i3 = max(1, 128 * ii0 + 64 * ii3 - 2 * i0 + 3); i3 <= 64 * ii3; i3 += 1) {
                        for (int i4 = 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 4; i4 <= 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 67; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                    for (int i3 = 64 * ii3 + 1; i3 <= 128 * ii0 + 64 * ii3 - 2 * i0 + 66; i3 += 1) {
                      for (int i4 = 128 * ii0 + 64 * ii4 - 2 * i0 + 3; i4 <= 64 * ii4; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                      if (_PB_N >= 64 * ii4 + 67) {
                        for (int i4 = 64 * ii4 + 1; i4 <= 128 * ii0 + 64 * ii4 - 2 * i0 + 66; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      } else {
                        for (int i4 = 64 * ii2 + 1; i4 < _PB_N - 1; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                  }
                }
                if (_PB_N >= 64 * ii4 + 67) {
                  for (int i0 = 64 * ii0 + 33; i0 <= min(_PB_TSTEPS, 64 * ii0 + 64); i0 += 1) {
                    for (int i2 = 128 * ii0 + 64 * ii2 - 2 * i0 + 4; i2 <= min(64 * ii2, 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 + 67); i2 += 1) {
                      for (int i3 = max(1, 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 5); i3 <= min(64 * ii3 + 1, 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 68); i3 += 1) {
                        for (int i4 = 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 6; i4 <= 64 * ii4 + 1; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        if (2 * i0 + i2 >= 128 * ii0 + 64 * ii2 + 7 && 2 * i0 + i2 + i3 >= 128 * ii0 + 64 * ii2 + 64 * ii3 + 7) {
                          for (int i4 = 64 * ii4 + 2; i4 <= 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 69; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        } else if (128 * ii0 + 64 * ii2 + 6 >= 2 * i0 + i2 && i3 == 64 * ii3 + 1) {
                          for (int i4 = 64 * ii4 + 2; i4 <= 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 68; i4 += 1) {
                            B[i2][64 * ii3 + 1][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 1][i4])) + A[i2 - 1][64 * ii3 + 1][i4])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 2][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 1][i4])) + A[i2][64 * ii3][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 1][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 1][i4])) + A[i2][64 * ii3 + 1][i4 - 1]))) + A[i2][64 * ii3 + 1][i4]);
                          }
                        } else {
                          for (int i4 = 64 * ii4 + 2; i4 <= 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 69; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                      if (2 * i0 + i2 == 128 * ii0 + 64 * ii2 + 4) {
                        for (int i4 = 64 * ii4 + 1; i4 <= 64 * ii4 + 64; i4 += 1) {
                          B[128 * ii0 + 64 * ii2 - 2 * i0 + 4][64 * ii3 + 2][i4] = ((((SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][64 * ii3 + 2][i4] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][64 * ii3 + 2][i4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 3][64 * ii3 + 2][i4])) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][64 * ii3 + 3][i4] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][64 * ii3 + 2][i4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][64 * ii3 + 1][i4]))) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][64 * ii3 + 2][i4 + 1] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][64 * ii3 + 2][i4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][64 * ii3 + 2][i4 - 1]))) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][64 * ii3 + 2][i4]);
                        }
                      }
                      for (int i3 = max(64 * ii3 + 2, 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 7); i3 <= 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 68; i3 += 1) {
                        for (int i4 = 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 5; i4 <= 64 * ii4; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        if (128 * ii0 + 64 * ii2 + 5 >= 2 * i0 + i2 && i3 >= 64 * ii3 + 3) {
                          for (int i4 = 64 * ii4 + 1; i4 <= 256 * ii0 + 128 * ii2 + 64 * ii4 - 4 * i0 - 2 * i2 + 72; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (2 * i0 + i2 == 128 * ii0 + 64 * ii2 + 5) {
                            B[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][64 * ii4 + 63] = ((((SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][64 * ii4 + 63] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][64 * ii4 + 63])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][i3][64 * ii4 + 63])) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3 + 1][64 * ii4 + 63] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][64 * ii4 + 63])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3 - 1][64 * ii4 + 63]))) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][64 * ii4 + 64] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][64 * ii4 + 63])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][64 * ii4 + 62]))) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][64 * ii4 + 63]);
                          }
                        } else {
                          if (128 * ii0 + 64 * ii2 + 6 >= 2 * i0 + i2 && i3 == 64 * ii3 + 2) {
                            for (int i4 = 64 * ii4 + 1; i4 <= 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 68; i4 += 1) {
                              B[i2][64 * ii3 + 2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii3 + 2][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 2][i4])) + A[i2 - 1][64 * ii3 + 2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 3][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 2][i4])) + A[i2][64 * ii3 + 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 2][i4])) + A[i2][64 * ii3 + 2][i4 - 1]))) + A[i2][64 * ii3 + 2][i4]);
                            }
                          }
                          if (2 * i0 + i2 >= 128 * ii0 + 64 * ii2 + 6 && 2 * i0 + i2 + i3 >= 128 * ii0 + 64 * ii2 + 64 * ii3 + 9) {
                            for (int i4 = 64 * ii4 + 1; i4 <= 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 68; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                        }
                      }
                    }
                    for (int i2 = 64 * ii2 + 1; i2 < _PB_N - 1; i2 += 1) {
                      for (int i3 = max(1, 128 * ii0 + 64 * ii3 - 2 * i0 + 4); i3 <= 128 * ii0 + 64 * ii3 - 2 * i0 + 67; i3 += 1) {
                        for (int i4 = 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 5; i4 <= 64 * ii4 + 1; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        if (i2 == 64 * ii2 + 1 && 2 * i0 + i3 == 128 * ii0 + 64 * ii3 + 5) {
                          for (int i4 = 64 * ii4 + 2; i4 <= 64 * ii4 + 63; i4 += 1) {
                            B[64 * ii2 + 1][128 * ii0 + 64 * ii3 - 2 * i0 + 5][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][128 * ii0 + 64 * ii3 - 2 * i0 + 5][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][128 * ii0 + 64 * ii3 - 2 * i0 + 5][i4])) + A[64 * ii2][128 * ii0 + 64 * ii3 - 2 * i0 + 5][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][128 * ii0 + 64 * ii3 - 2 * i0 + 6][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][128 * ii0 + 64 * ii3 - 2 * i0 + 5][i4])) + A[64 * ii2 + 1][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][128 * ii0 + 64 * ii3 - 2 * i0 + 5][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][128 * ii0 + 64 * ii3 - 2 * i0 + 5][i4])) + A[64 * ii2 + 1][128 * ii0 + 64 * ii3 - 2 * i0 + 5][i4 - 1]))) + A[64 * ii2 + 1][128 * ii0 + 64 * ii3 - 2 * i0 + 5][i4]);
                          }
                        } else if (128 * ii0 + 64 * ii3 + i2 + 3 >= 64 * ii2 + 2 * i0 + i3) {
                          for (int i4 = 64 * ii4 + 2; i4 <= 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 68; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        } else {
                          for (int i4 = 64 * ii4 + 2; i4 <= 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 68; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                    }
                    for (int i2 = 128 * ii0 + 64 * ii2 - 2 * i0 + 3; i2 <= min(64 * ii2, 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 + 66); i2 += 1) {
                      for (int i3 = max(1, 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 4); i3 <= 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 67; i3 += 1) {
                        for (int i4 = max(128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 4, 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 5); i4 <= 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 68; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                        for (int i4 = 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 69; i4 <= 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 67; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                    for (int i2 = 64 * ii2 + 1; i2 < _PB_N - 1; i2 += 1) {
                      for (int i3 = max(1, 128 * ii0 + 64 * ii3 - 2 * i0 + 3); i3 <= 128 * ii0 + 64 * ii3 - 2 * i0 + 66; i3 += 1) {
                        for (int i4 = 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 4; i4 <= 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 67; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                  }
                } else {
                  for (int i0 = 64 * ii0 + 33; i0 <= min(_PB_TSTEPS, 64 * ii0 + 64); i0 += 1) {
                    if (i0 >= 64 * ii0 + 34) {
                      for (int i2 = max(1, 128 * ii0 + 64 * ii2 - 2 * i0 + 4); i2 <= min(64 * ii2, 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 + 67); i2 += 1) {
                        for (int i3 = max(1, 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 5); i3 <= 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 68; i3 += 1) {
                          if (64 * ii3 + 2 >= i3) {
                            for (int i4 = max(128 * ii0 + 128 * ii2 - 2 * i0 - i2 + 5, 128 * ii0 + 128 * ii2 + 64 * ii3 - 2 * i0 - i2 - i3 + 6); i4 < _PB_N - 1; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          } else if (64 * ii3 + 63 >= i3) {
                            for (int i4 = max(-256 * ii0 - 64 * ii2 + 4 * i0 + 2 * i2 - 10, 256 * ii0 + 192 * ii2 - 4 * i0 - 2 * i2 + 9); i4 <= 64 * ii2 + 1; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            if (2 * i0 + i2 >= 128 * ii0 + 64 * ii2 + 6) {
                              for (int i4 = 128 * ii0 + 128 * ii2 - 2 * i0 - i2 + 5; i4 <= 64 * ii2 + 1; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                            for (int i4 = 64 * ii2 + 2; i4 < _PB_N - 1; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          } else {
                            for (int i4 = 64 * ii2 + 1; i4 < _PB_N - 1; i4 += 1) {
                              B[128 * ii0 + 64 * ii2 - 2 * i0 + 4][64 * ii3 + 64][i4] = ((((SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][64 * ii3 + 64][i4] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][64 * ii3 + 64][i4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 3][64 * ii3 + 64][i4])) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][64 * ii3 + 65][i4] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][64 * ii3 + 64][i4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][64 * ii3 + 63][i4]))) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][64 * ii3 + 64][i4 + 1] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][64 * ii3 + 64][i4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][64 * ii3 + 64][i4 - 1]))) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][64 * ii3 + 64][i4]);
                            }
                          }
                        }
                      }
                      for (int i3 = max(1, 128 * ii0 + 64 * ii3 - 2 * i0 + 4); i3 <= 128 * ii0 + 64 * ii3 - 2 * i0 + 67; i3 += 1) {
                        for (int i4 = 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i3 + 5; i4 < _PB_N - 1; i4 += 1) {
                          B[64 * ii2 + 1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][i3][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][i3][i4])) + A[64 * ii2][i3][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][i3][i4])) + A[64 * ii2 + 1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][i3][i4])) + A[64 * ii2 + 1][i3][i4 - 1]))) + A[64 * ii2 + 1][i3][i4]);
                        }
                      }
                      for (int i2 = 64 * ii2 + 2; i2 <= min(_PB_N - 2, 64 * ii2 + 63); i2 += 1) {
                        for (int i3 = max(1, 128 * ii0 + 64 * ii3 - 2 * i0 + 4); i3 <= 128 * ii0 + 64 * ii3 - 2 * i0 + 67; i3 += 1) {
                          for (int i4 = 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i3 + 5; i4 <= min(64 * ii2 + 1, 128 * ii0 + 64 * ii3 - 2 * i0 + i2 - i3 + 4); i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (64 * ii2 + 2 * i0 + i3 >= 128 * ii0 + 64 * ii3 + i2 + 4) {
                            for (int i4 = 128 * ii0 + 64 * ii3 - 2 * i0 + i2 - i3 + 5; i4 < _PB_N - 1; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          } else {
                            for (int i4 = 64 * ii2 + 2; i4 < _PB_N - 1; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                        }
                      }
                      if (64 * ii2 + 66 == _PB_N) {
                        for (int i3 = max(1, 128 * ii0 + 64 * ii3 - 2 * i0 + 4); i3 <= 128 * ii0 + 64 * ii3 - 2 * i0 + 67; i3 += 1) {
                          for (int i4 = _PB_N + 128 * ii0 + 64 * ii3 - 2 * i0 - i3 - 61; i4 < _PB_N - 1; i4 += 1) {
                            B[_PB_N - 2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[_PB_N - 1][i3][i4] - (SCALAR_VAL(2.0) * A[_PB_N - 2][i3][i4])) + A[_PB_N - 3][i3][i4])) + (SCALAR_VAL(0.125) * ((A[_PB_N - 2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[_PB_N - 2][i3][i4])) + A[_PB_N - 2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[_PB_N - 2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[_PB_N - 2][i3][i4])) + A[_PB_N - 2][i3][i4 - 1]))) + A[_PB_N - 2][i3][i4]);
                          }
                        }
                      }
                    } else {
                      for (int i2 = 64 * ii2 - 62; i2 <= 64 * ii2; i2 += 1) {
                        for (int i3 = max(1, 64 * ii2 + 64 * ii3 - i2 - 61); i3 <= min(64 * ii3 + 63, 64 * ii2 + 64 * ii3 - i2 + 2); i3 += 1) {
                          if (i3 >= 64 * ii3 + 3) {
                            for (int i4 = max(192 * ii2 - 2 * i2 - 123, -64 * ii2 + 2 * i2 + 122); i4 <= 64 * ii2 + 1; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            if (i2 + 60 >= 64 * ii2) {
                              for (int i4 = 128 * ii2 - i2 - 61; i4 <= 64 * ii2 + 1; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                            for (int i4 = 64 * ii2 + 2; i4 < _PB_N - 1; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          } else {
                            for (int i4 = max(128 * ii2 - i2 - 61, 128 * ii2 + 64 * ii3 - i2 - i3 - 60); i4 < _PB_N - 1; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                        }
                        if (i2 + 62 == 64 * ii2) {
                          for (int i4 = 64 * ii2 + 1; i4 < _PB_N - 1; i4 += 1) {
                            B[64 * ii2 - 62][64 * ii3 + 64][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 - 61][64 * ii3 + 64][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 - 62][64 * ii3 + 64][i4])) + A[64 * ii2 - 63][64 * ii3 + 64][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 - 62][64 * ii3 + 65][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 - 62][64 * ii3 + 64][i4])) + A[64 * ii2 - 62][64 * ii3 + 63][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 - 62][64 * ii3 + 64][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii2 - 62][64 * ii3 + 64][i4])) + A[64 * ii2 - 62][64 * ii3 + 64][i4 - 1]))) + A[64 * ii2 - 62][64 * ii3 + 64][i4]);
                          }
                        }
                      }
                      for (int i3 = max(1, 64 * ii3 - 62); i3 <= 64 * ii3 + 1; i3 += 1) {
                        for (int i4 = 64 * ii2 + 64 * ii3 - i3 - 61; i4 < _PB_N - 1; i4 += 1) {
                          B[64 * ii2 + 1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][i3][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][i3][i4])) + A[64 * ii2][i3][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][i3][i4])) + A[64 * ii2 + 1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][i3][i4])) + A[64 * ii2 + 1][i3][i4 - 1]))) + A[64 * ii2 + 1][i3][i4]);
                        }
                      }
                      for (int i2 = 64 * ii2 + 2; i2 < _PB_N - 1; i2 += 1) {
                        for (int i3 = max(1, 64 * ii3 - 62); i3 <= 64 * ii3 + 1; i3 += 1) {
                          if (i3 >= 64 * ii3) {
                            for (int i4 = 64 * ii2 + 64 * ii3 - i3 - 61; i4 <= min(64 * ii3 + i2 - i3 - 62, 64 * ii2 + 64 * ii3 - i3); i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            if (i2 + i3 >= 64 * ii2 + 64 * ii3 + 3) {
                              for (int i4 = 64 * ii3 + i2 - i3 - 61; i4 <= 64 * ii2 + 64 * ii3 - i3; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                          } else {
                            for (int i4 = 64 * ii2 + 64 * ii3 - i3 - 61; i4 <= min(64 * ii2 + 1, 64 * ii3 + i2 - i3 - 62); i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            if (i2 + i3 >= 64 * ii2 + 64 * ii3 + 3 && 64 * ii2 + i3 + 62 >= 64 * ii3 + i2) {
                              for (int i4 = 64 * ii3 + i2 - i3 - 61; i4 < _PB_N - 1; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                          }
                          if (64 * ii2 + 64 * ii3 + 2 >= i2 + i3 && 64 * ii2 + i3 + 62 >= 64 * ii3 + i2) {
                            for (int i4 = 64 * ii3 + i2 - i3 - 61; i4 < _PB_N - 1; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          } else if (64 * ii3 + i2 >= 64 * ii2 + i3 + 63 && 64 * ii3 >= i3 + 1) {
                            for (int i4 = 64 * ii2 + 2; i4 < _PB_N - 1; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          } else if (i3 >= 64 * ii3) {
                            for (int i4 = 64 * ii2 + 64 * ii3 - i3 + 1; i4 < _PB_N - 1; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                        }
                      }
                    }
                    for (int i2 = max(1, 128 * ii0 + 64 * ii2 - 2 * i0 + 3); i2 <= min(64 * ii2, 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 + 66); i2 += 1) {
                      for (int i3 = max(1, 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 4); i3 <= 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 67; i3 += 1) {
                        for (int i4 = max(128 * ii0 + 128 * ii2 - 2 * i0 - i2 + 4, 128 * ii0 + 128 * ii2 + 64 * ii3 - 2 * i0 - i2 - i3 + 5); i4 < _PB_N - 1; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                    for (int i2 = 64 * ii2 + 1; i2 < _PB_N - 1; i2 += 1) {
                      for (int i3 = max(1, 128 * ii0 + 64 * ii3 - 2 * i0 + 3); i3 <= 128 * ii0 + 64 * ii3 - 2 * i0 + 66; i3 += 1) {
                        for (int i4 = 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i3 + 4; i4 < _PB_N - 1; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                  }
                }
              } else {
                for (int i2 = 64 * ii2 + 1; i2 < _PB_N - 1; i2 += 1) {
                  for (int i3 = 64 * ii3 + 1; i3 <= 64 * ii3 + 64; i3 += 1) {
                    for (int i4 = 1; i4 <= 64; i4 += 1) {
                      A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                    }
                  }
                }
                for (int i0 = 64 * ii0 + 2; i0 <= min(_PB_TSTEPS, 64 * ii0 + 64); i0 += 1) {
                  if (i0 >= 64 * ii0 + 3) {
                    for (int i2 = max(1, 128 * ii0 + 64 * ii2 - 2 * i0 + 4); i2 <= min(64 * ii2 + 1, 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 + 67); i2 += 1) {
                      if (64 * ii0 + 32 * ii3 + 1 >= i0 && i2 == 64 * ii2 + 1) {
                        for (int i4 = 1; i4 <= 64; i4 += 1) {
                          B[64 * ii2 + 1][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4])) + A[64 * ii2][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][128 * ii0 + 64 * ii3 - 2 * i0 + 5][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4])) + A[64 * ii2 + 1][128 * ii0 + 64 * ii3 - 2 * i0 + 3][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4])) + A[64 * ii2 + 1][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4 - 1]))) + A[64 * ii2 + 1][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4]);
                        }
                      }
                      for (int i3 = max(max(1, 128 * ii0 - 64 * ii2 + 64 * ii3 - 2 * i0 + i2 + 4), 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 5); i3 <= min(64 * ii3, 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 68); i3 += 1) {
                        if (2 * i0 + i2 >= 128 * ii0 + 64 * ii2 + 7 && 64 * ii2 >= i2 && i3 == 64 * ii3) {
                          B[i2][64 * ii3][1] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii3][1] - (SCALAR_VAL(2.0) * A[i2][64 * ii3][1])) + A[i2 - 1][64 * ii3][1])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 1][1] - (SCALAR_VAL(2.0) * A[i2][64 * ii3][1])) + A[i2][64 * ii3 - 1][1]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3][2] - (SCALAR_VAL(2.0) * A[i2][64 * ii3][1])) + A[i2][64 * ii3][0]))) + A[i2][64 * ii3][1]);
                        } else if (128 * ii0 + 64 * ii2 + 6 >= 2 * i0 + i2 && i3 == 64 * ii3) {
                          B[i2][64 * ii3][1] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii3][1] - (SCALAR_VAL(2.0) * A[i2][64 * ii3][1])) + A[i2 - 1][64 * ii3][1])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 1][1] - (SCALAR_VAL(2.0) * A[i2][64 * ii3][1])) + A[i2][64 * ii3 - 1][1]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3][2] - (SCALAR_VAL(2.0) * A[i2][64 * ii3][1])) + A[i2][64 * ii3][0]))) + A[i2][64 * ii3][1]);
                        } else if (64 * ii3 >= i3 + 1) {
                          B[i2][i3][1] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][1] - (SCALAR_VAL(2.0) * A[i2][i3][1])) + A[i2 - 1][i3][1])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][1] - (SCALAR_VAL(2.0) * A[i2][i3][1])) + A[i2][i3 - 1][1]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][2] - (SCALAR_VAL(2.0) * A[i2][i3][1])) + A[i2][i3][0]))) + A[i2][i3][1]);
                        }
                        if (2 * i0 + i2 + i3 >= 128 * ii0 + 64 * ii2 + 64 * ii3 + 7) {
                          for (int i4 = max(ii2 - (58 * ii2 + i2 + 121) / 122 + 2, ii3 - (58 * ii3 + i3 + 122) / 122 + 2); i4 <= 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 - i3 + 69; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        } else {
                          for (int i4 = 2; i4 <= 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 - i3 + 69; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                      if (2 * i0 + i2 >= 128 * ii0 + 64 * ii2 + 7 && 128 * ii0 + 64 * ii2 + 67 >= 2 * i0 + i2) {
                        if (64 * ii2 >= i2 + 1) {
                          B[i2][64 * ii3 + 1][1] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii3 + 1][1] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 1][1])) + A[i2 - 1][64 * ii3 + 1][1])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 2][1] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 1][1])) + A[i2][64 * ii3][1]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 1][2] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 1][1])) + A[i2][64 * ii3 + 1][0]))) + A[i2][64 * ii3 + 1][1]);
                        }
                        for (int i4 = ii2 - (57 * ii2 + i2 + 121) / 121 + 2; i4 <= 128 * ii0 + 64 * ii2 - 2 * i0 - i2 + 68; i4 += 1) {
                          B[i2][64 * ii3 + 1][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 1][i4])) + A[i2 - 1][64 * ii3 + 1][i4])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 2][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 1][i4])) + A[i2][64 * ii3][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 1][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 1][i4])) + A[i2][64 * ii3 + 1][i4 - 1]))) + A[i2][64 * ii3 + 1][i4]);
                        }
                      } else if (128 * ii0 + 64 * ii2 + 6 >= 2 * i0 + i2) {
                        for (int i3 = 64 * ii3 + 1; i3 <= 64 * ii3 + 2; i3 += 1) {
                          if (i3 == 64 * ii3 + 1) {
                            B[i2][64 * ii3 + 1][1] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii3 + 1][1] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 1][1])) + A[i2 - 1][64 * ii3 + 1][1])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 2][1] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 1][1])) + A[i2][64 * ii3][1]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 1][2] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 1][1])) + A[i2][64 * ii3 + 1][0]))) + A[i2][64 * ii3 + 1][1]);
                          }
                          if (64 * ii2 >= i2 + 1) {
                            for (int i4 = 64 * ii3 - i3 + 3; i4 <= 128 * ii0 + 64 * ii2 - 2 * i0 - i2 + 68; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          } else {
                            if (i3 == 64 * ii3 + 2) {
                              B[64 * ii2][64 * ii3 + 2][1] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][64 * ii3 + 2][1] - (SCALAR_VAL(2.0) * A[64 * ii2][64 * ii3 + 2][1])) + A[64 * ii2 - 1][64 * ii3 + 2][1])) + (SCALAR_VAL(0.125) * ((A[64 * ii2][64 * ii3 + 3][1] - (SCALAR_VAL(2.0) * A[64 * ii2][64 * ii3 + 2][1])) + A[64 * ii2][64 * ii3 + 1][1]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2][64 * ii3 + 2][2] - (SCALAR_VAL(2.0) * A[64 * ii2][64 * ii3 + 2][1])) + A[64 * ii2][64 * ii3 + 2][0]))) + A[64 * ii2][64 * ii3 + 2][1]);
                            }
                            for (int i4 = 2; i4 <= 62; i4 += 1) {
                              B[64 * ii2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[64 * ii2][i3][i4])) + A[64 * ii2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[64 * ii2][i3][i4])) + A[64 * ii2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii2][i3][i4])) + A[64 * ii2][i3][i4 - 1]))) + A[64 * ii2][i3][i4]);
                            }
                          }
                        }
                        if (128 * ii0 + 64 * ii2 + 5 >= 2 * i0 + i2) {
                          for (int i3 = 64 * ii3 + 3; i3 <= 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 68; i3 += 1) {
                            for (int i4 = 1; i4 <= 256 * ii0 + 128 * ii2 - 4 * i0 - 2 * i2 + 72; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            if (2 * i0 + i2 == 128 * ii0 + 64 * ii2 + 5) {
                              B[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][63] = ((((SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][63] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][63])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][i3][63])) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3 + 1][63] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][63])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3 - 1][63]))) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][64] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][63])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][62]))) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][63]);
                            }
                          }
                        }
                      }
                      if (2 * i0 + i2 >= 128 * ii0 + 64 * ii2 + 6) {
                        for (int i3 = max(64 * ii3 + 2, 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 9); i3 <= 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 68; i3 += 1) {
                          for (int i4 = 1; i4 <= 128 * ii0 + 64 * ii2 - 2 * i0 - i2 + 68; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                    }
                  }
                  if (64 * ii0 + 32 * ii3 + 32 >= i0) {
                    for (int i2 = 64 * ii2 + 2; i2 <= min(min(_PB_N - 2, 64 * ii2 + 62), -128 * ii0 + 64 * ii2 + 2 * i0 - 3); i2 += 1) {
                      for (int i3 = max(1, 128 * ii0 + 64 * ii3 - 2 * i0 + 4); i3 <= 128 * ii0 + 64 * ii3 - 2 * i0 + 67; i3 += 1) {
                        if (i3 >= 64 * ii3 + 2) {
                          for (int i4 = 1; i4 <= 128 * ii0 - 2 * i0 + 67; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        } else {
                          if (128 * ii0 + 64 * ii3 + 5 >= 2 * i0 + i3) {
                            B[i2][i3][1] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][1] - (SCALAR_VAL(2.0) * A[i2][i3][1])) + A[i2 - 1][i3][1])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][1] - (SCALAR_VAL(2.0) * A[i2][i3][1])) + A[i2][i3 - 1][1]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][2] - (SCALAR_VAL(2.0) * A[i2][i3][1])) + A[i2][i3][0]))) + A[i2][i3][1]);
                          }
                          for (int i4 = max(1, 64 * ii0 + 32 * ii3 - i0 - i3 + (i3 + 1) / 2 + 4); i4 <= 128 * ii0 + 64 * ii3 - 2 * i0 - i3 + 68; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                    }
                  }
                  if (_PB_N >= 64 * ii2 + 65 && i0 >= 64 * ii0 + 34) {
                    for (int i3 = max(1, 128 * ii0 + 64 * ii3 - 2 * i0 + 4); i3 <= 128 * ii0 + 64 * ii3 - 2 * i0 + 67; i3 += 1) {
                      if (2 * i0 + i3 >= 128 * ii0 + 64 * ii3 + 6) {
                        B[64 * ii2 + 63][i3][1] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 64][i3][1] - (SCALAR_VAL(2.0) * A[64 * ii2 + 63][i3][1])) + A[64 * ii2 + 62][i3][1])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 63][i3 + 1][1] - (SCALAR_VAL(2.0) * A[64 * ii2 + 63][i3][1])) + A[64 * ii2 + 63][i3 - 1][1]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 63][i3][2] - (SCALAR_VAL(2.0) * A[64 * ii2 + 63][i3][1])) + A[64 * ii2 + 63][i3][0]))) + A[64 * ii2 + 63][i3][1]);
                      }
                      for (int i4 = -2 * ii0 - ii3 + (-6 * ii0 - 3 * ii3 + 2 * i0 + i3 - 6) / 61 + 2; i4 <= 128 * ii0 + 64 * ii3 - 2 * i0 - i3 + 68; i4 += 1) {
                        B[64 * ii2 + 63][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 64][i3][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 63][i3][i4])) + A[64 * ii2 + 62][i3][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 63][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 63][i3][i4])) + A[64 * ii2 + 63][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 63][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii2 + 63][i3][i4])) + A[64 * ii2 + 63][i3][i4 - 1]))) + A[64 * ii2 + 63][i3][i4]);
                      }
                    }
                    if (64 * ii2 + 66 == _PB_N) {
                      for (int i3 = max(1, 128 * ii0 + 64 * ii3 - 2 * i0 + 4); i3 <= 128 * ii0 + 64 * ii3 - 2 * i0 + 67; i3 += 1) {
                        if (2 * i0 + i3 >= 128 * ii0 + 64 * ii3 + 6) {
                          B[_PB_N - 2][i3][1] = ((((SCALAR_VAL(0.125) * ((A[_PB_N - 1][i3][1] - (SCALAR_VAL(2.0) * A[_PB_N - 2][i3][1])) + A[_PB_N - 3][i3][1])) + (SCALAR_VAL(0.125) * ((A[_PB_N - 2][i3 + 1][1] - (SCALAR_VAL(2.0) * A[_PB_N - 2][i3][1])) + A[_PB_N - 2][i3 - 1][1]))) + (SCALAR_VAL(0.125) * ((A[_PB_N - 2][i3][2] - (SCALAR_VAL(2.0) * A[_PB_N - 2][i3][1])) + A[_PB_N - 2][i3][0]))) + A[_PB_N - 2][i3][1]);
                        }
                        for (int i4 = -2 * ii0 - ii3 + (-6 * ii0 - 3 * ii3 + 2 * i0 + i3 - 6) / 61 + 2; i4 <= 128 * ii0 + 64 * ii3 - 2 * i0 - i3 + 68; i4 += 1) {
                          B[_PB_N - 2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[_PB_N - 1][i3][i4] - (SCALAR_VAL(2.0) * A[_PB_N - 2][i3][i4])) + A[_PB_N - 3][i3][i4])) + (SCALAR_VAL(0.125) * ((A[_PB_N - 2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[_PB_N - 2][i3][i4])) + A[_PB_N - 2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[_PB_N - 2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[_PB_N - 2][i3][i4])) + A[_PB_N - 2][i3][i4 - 1]))) + A[_PB_N - 2][i3][i4]);
                        }
                      }
                    }
                  } else if (i0 == 64 * ii0 + 2) {
                    for (int i2 = 64 * ii2; i2 <= min(_PB_N - 2, 64 * ii2 + 2); i2 += 1) {
                      if (64 * ii2 + 1 >= i2) {
                        for (int i3 = max(1, 64 * ii2 + 64 * ii3 - i2 + 1); i3 <= 64 * ii2 + 64 * ii3 - i2 + 64; i3 += 1) {
                          for (int i4 = 1; i4 <= 64 * ii2 - i2 + 64; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (i2 == 64 * ii2 + 1 && i3 == 64 * ii3) {
                            B[64 * ii2 + 1][64 * ii3][64] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][64 * ii3][64] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][64 * ii3][64])) + A[64 * ii2][64 * ii3][64])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][64 * ii3 + 1][64] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][64 * ii3][64])) + A[64 * ii2 + 1][64 * ii3 - 1][64]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][64 * ii3][65] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][64 * ii3][64])) + A[64 * ii2 + 1][64 * ii3][63]))) + A[64 * ii2 + 1][64 * ii3][64]);
                          }
                        }
                      } else {
                        if (ii3 >= 1) {
                          for (int i4 = 1; i4 <= 64; i4 += 1) {
                            B[64 * ii2 + 2][64 * ii3][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 3][64 * ii3][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][64 * ii3][i4])) + A[64 * ii2 + 1][64 * ii3][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][64 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][64 * ii3][i4])) + A[64 * ii2 + 2][64 * ii3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][64 * ii3][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][64 * ii3][i4])) + A[64 * ii2 + 2][64 * ii3][i4 - 1]))) + A[64 * ii2 + 2][64 * ii3][i4]);
                          }
                        }
                        for (int i3 = 64 * ii3 + 1; i3 <= 64 * ii3 + 63; i3 += 1) {
                          for (int i4 = 1; i4 <= 63; i4 += 1) {
                            B[64 * ii2 + 2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 3][i3][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][i3][i4])) + A[64 * ii2 + 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][i3][i4])) + A[64 * ii2 + 2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][i3][i4])) + A[64 * ii2 + 2][i3][i4 - 1]))) + A[64 * ii2 + 2][i3][i4]);
                          }
                        }
                      }
                    }
                  }
                  for (int i2 = max(-128 * ii0 + 64 * ii2 + 2 * i0 - 2, 128 * ii0 + 64 * ii2 - 2 * i0 + 7); i2 <= min(_PB_N - 2, 64 * ii2 + 62); i2 += 1) {
                    if (ii3 >= 1 && i0 == 64 * ii0 + 2) {
                      for (int i4 = 1; i4 <= 64; i4 += 1) {
                        B[i2][64 * ii3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii3][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii3][i4])) + A[i2 - 1][64 * ii3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii3][i4])) + A[i2][64 * ii3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii3][i4])) + A[i2][64 * ii3][i4 - 1]))) + A[i2][64 * ii3][i4]);
                      }
                    } else if (i0 >= 64 * ii0 + 3) {
                      for (int i3 = max(1, 128 * ii0 + 64 * ii3 - 2 * i0 + 4); i3 <= min(64 * ii3 + 1, 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 68); i3 += 1) {
                        B[i2][i3][1] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][1] - (SCALAR_VAL(2.0) * A[i2][i3][1])) + A[i2 - 1][i3][1])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][1] - (SCALAR_VAL(2.0) * A[i2][i3][1])) + A[i2][i3 - 1][1]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][2] - (SCALAR_VAL(2.0) * A[i2][i3][1])) + A[i2][i3][0]))) + A[i2][i3][1]);
                        for (int i4 = 2; i4 <= 128 * ii0 + 64 * ii3 - 2 * i0 - i3 + 68; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    }
                    for (int i3 = 64 * ii3 + 2; i3 <= min(-128 * ii0 + 64 * ii3 + 2 * i0 - 5, 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 68); i3 += 1) {
                      for (int i4 = 1; i4 <= 128 * ii0 - 2 * i0 + 67; i4 += 1) {
                        B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                      }
                    }
                    for (int i3 = max(1, 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 69); i3 < 64 * ii3; i3 += 1) {
                      for (int i4 = 1; i4 <= 128 * ii0 + 64 * ii3 - 2 * i0 - i3 + 68; i4 += 1) {
                        B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                      }
                    }
                    for (int i3 = max(max(1, 64 * ii3), 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 69); i3 <= min(64 * ii3 + 2, -128 * ii0 + 64 * ii3 + 2 * i0 - 5); i3 += 1) {
                      if (64 * ii3 + 1 >= i3) {
                        for (int i4 = 1; i4 <= 128 * ii0 + 64 * ii3 - 2 * i0 - i3 + 68; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      } else {
                        for (int i4 = 1; i4 <= 128 * ii0 - 2 * i0 + 67; i4 += 1) {
                          B[i2][64 * ii3 + 2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii3 + 2][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 2][i4])) + A[i2 - 1][64 * ii3 + 2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 3][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 2][i4])) + A[i2][64 * ii3 + 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 2][i4])) + A[i2][64 * ii3 + 2][i4 - 1]))) + A[i2][64 * ii3 + 2][i4]);
                        }
                      }
                    }
                    if (i0 == 64 * ii0 + 3 && i2 >= 64 * ii2 + 61) {
                      for (int i4 = 1; i4 <= 61; i4 += 1) {
                        B[i2][64 * ii3 + 2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii3 + 2][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 2][i4])) + A[i2 - 1][64 * ii3 + 2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 3][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 2][i4])) + A[i2][64 * ii3 + 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 2][i4])) + A[i2][64 * ii3 + 2][i4 - 1]))) + A[i2][64 * ii3 + 2][i4]);
                      }
                    }
                    for (int i3 = max(64 * ii3 + 1, -128 * ii0 + 64 * ii3 + 2 * i0 - 4); i3 <= 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 68; i3 += 1) {
                      if (i0 >= 64 * ii0 + 3) {
                        for (int i4 = 1; i4 <= 128 * ii0 - 2 * i0 + 67; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      } else {
                        for (int i4 = 1; i4 <= min(63, 64 * ii3 - i3 + 65); i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = 64 * ii3 - i3 + 66; i4 <= 63; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    }
                    for (int i3 = max(64 * ii3 + 3, 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 69); i3 <= 128 * ii0 + 64 * ii3 - 2 * i0 + 67; i3 += 1) {
                      for (int i4 = 1; i4 <= 128 * ii0 - 2 * i0 + 67; i4 += 1) {
                        B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                      }
                    }
                  }
                  if (64 * ii0 + 32 * ii3 + 32 >= i0 && 64 * ii0 + 33 >= i0) {
                    for (int i2 = 64 * ii2 + 63; i2 < _PB_N - 1; i2 += 1) {
                      for (int i3 = max(1, 128 * ii0 + 64 * ii3 - 2 * i0 + 4); i3 < 64 * ii3; i3 += 1) {
                        if (2 * i0 + i3 >= 128 * ii0 + 64 * ii3 + 6) {
                          B[i2][i3][1] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][1] - (SCALAR_VAL(2.0) * A[i2][i3][1])) + A[i2 - 1][i3][1])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][1] - (SCALAR_VAL(2.0) * A[i2][i3][1])) + A[i2][i3 - 1][1]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][2] - (SCALAR_VAL(2.0) * A[i2][i3][1])) + A[i2][i3][0]))) + A[i2][i3][1]);
                        }
                        for (int i4 = -2 * ii0 - ii3 + (-8 * ii0 - 4 * ii3 + 2 * i0 + i3 - 6) / 60 + 2; i4 <= 128 * ii0 + 64 * ii3 - 2 * i0 - i3 + 68; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                      for (int i3 = max(1, 64 * ii3); i3 <= 128 * ii0 + 64 * ii3 - 2 * i0 + 67; i3 += 1) {
                        if (i0 == 64 * ii0 + 2 && 64 * ii3 + 1 >= i3) {
                          for (int i4 = 1; i4 <= min(63, 128 * ii3 - 2 * i3 + 64); i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (i2 == 64 * ii2 + 63 && i3 == 64 * ii3 + 1) {
                            B[64 * ii2 + 63][64 * ii3 + 1][63] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 64][64 * ii3 + 1][63] - (SCALAR_VAL(2.0) * A[64 * ii2 + 63][64 * ii3 + 1][63])) + A[64 * ii2 + 62][64 * ii3 + 1][63])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 63][64 * ii3 + 2][63] - (SCALAR_VAL(2.0) * A[64 * ii2 + 63][64 * ii3 + 1][63])) + A[64 * ii2 + 63][64 * ii3][63]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 63][64 * ii3 + 1][64] - (SCALAR_VAL(2.0) * A[64 * ii2 + 63][64 * ii3 + 1][63])) + A[64 * ii2 + 63][64 * ii3 + 1][62]))) + A[64 * ii2 + 63][64 * ii3 + 1][63]);
                          }
                        } else if (i3 >= 64 * ii3 + 2 && 128 * ii0 + 64 * ii3 + i2 + 3 >= 64 * ii2 + 2 * i0 + i3) {
                          for (int i4 = 1; i4 <= min(128 * ii3 - 2 * i3 + 64, 128 * ii0 + 64 * ii3 - 2 * i0 - i3 + 69); i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        } else if (i2 == 64 * ii2 + 63 && 2 * i0 + i3 == 128 * ii0 + 64 * ii3 + 67) {
                          for (int i4 = 1; i4 <= min(min(2, -256 * ii0 + 4 * i0 - 70), 128 * ii0 - 2 * i0 + 67); i4 += 1) {
                            B[64 * ii2 + 63][128 * ii0 + 64 * ii3 - 2 * i0 + 67][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 64][128 * ii0 + 64 * ii3 - 2 * i0 + 67][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 63][128 * ii0 + 64 * ii3 - 2 * i0 + 67][i4])) + A[64 * ii2 + 62][128 * ii0 + 64 * ii3 - 2 * i0 + 67][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 63][128 * ii0 + 64 * ii3 - 2 * i0 + 68][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 63][128 * ii0 + 64 * ii3 - 2 * i0 + 67][i4])) + A[64 * ii2 + 63][128 * ii0 + 64 * ii3 - 2 * i0 + 66][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 63][128 * ii0 + 64 * ii3 - 2 * i0 + 67][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii2 + 63][128 * ii0 + 64 * ii3 - 2 * i0 + 67][i4])) + A[64 * ii2 + 63][128 * ii0 + 64 * ii3 - 2 * i0 + 67][i4 - 1]))) + A[64 * ii2 + 63][128 * ii0 + 64 * ii3 - 2 * i0 + 67][i4]);
                          }
                        }
                        if (2 * i0 + i2 + i3 >= 128 * ii0 + 64 * ii2 + 64 * ii3 + 69) {
                          for (int i4 = max(1, 128 * ii3 - 2 * i3 + 65); i4 <= min(63, 128 * ii0 + 64 * ii3 - 2 * i0 - i3 + 69); i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (i0 >= 64 * ii0 + 3 && 64 * ii3 + 1 >= i3 && 128 * ii0 + 64 * ii3 + i2 + 3 >= 64 * ii2 + 2 * i0 + i3) {
                            for (int i4 = 1; i4 <= 128 * ii0 - 2 * i0 + 67; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                        }
                        if (i3 == 64 * ii3) {
                          B[i2][64 * ii3][128 * ii0 - 2 * i0 + 68] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii3][128 * ii0 - 2 * i0 + 68] - (SCALAR_VAL(2.0) * A[i2][64 * ii3][128 * ii0 - 2 * i0 + 68])) + A[i2 - 1][64 * ii3][128 * ii0 - 2 * i0 + 68])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 1][128 * ii0 - 2 * i0 + 68] - (SCALAR_VAL(2.0) * A[i2][64 * ii3][128 * ii0 - 2 * i0 + 68])) + A[i2][64 * ii3 - 1][128 * ii0 - 2 * i0 + 68]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3][128 * ii0 - 2 * i0 + 69] - (SCALAR_VAL(2.0) * A[i2][64 * ii3][128 * ii0 - 2 * i0 + 68])) + A[i2][64 * ii3][128 * ii0 - 2 * i0 + 67]))) + A[i2][64 * ii3][128 * ii0 - 2 * i0 + 68]);
                        }
                        for (int i4 = 128 * ii0 + 64 * ii3 - 2 * i0 - i3 + 70; i4 <= 128 * ii0 - 2 * i0 + 67; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    }
                  } else if (ii3 == 0 && i0 == 64 * ii0 + 33) {
                    for (int i2 = 64 * ii2 + 2; i2 < _PB_N - 1; i2 += 1) {
                      B[i2][1][1] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][1][1] - (SCALAR_VAL(2.0) * A[i2][1][1])) + A[i2 - 1][1][1])) + (SCALAR_VAL(0.125) * ((A[i2][2][1] - (SCALAR_VAL(2.0) * A[i2][1][1])) + A[i2][0][1]))) + (SCALAR_VAL(0.125) * ((A[i2][1][2] - (SCALAR_VAL(2.0) * A[i2][1][1])) + A[i2][1][0]))) + A[i2][1][1]);
                    }
                  }
                  for (int i2 = max(1, 128 * ii0 + 64 * ii2 - 2 * i0 + 3); i2 <= min(64 * ii2, 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 + 66); i2 += 1) {
                    for (int i3 = max(1, 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 4); i3 <= 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 67; i3 += 1) {
                      for (int i4 = 1; i4 <= 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 - i3 + 68; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                      for (int i4 = 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 - i3 + 69; i4 <= 128 * ii0 + 64 * ii2 - 2 * i0 - i2 + 67; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                  }
                  for (int i2 = 64 * ii2 + 1; i2 < _PB_N - 1; i2 += 1) {
                    for (int i3 = max(1, 128 * ii0 + 64 * ii3 - 2 * i0 + 3); i3 <= 128 * ii0 + 64 * ii3 - 2 * i0 + 66; i3 += 1) {
                      if (i3 >= 64 * ii3 + 1) {
                        for (int i4 = 1; i4 <= 128 * ii0 - 2 * i0 + 66; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      } else {
                        for (int i4 = 1; i4 <= 128 * ii0 + 64 * ii3 - 2 * i0 - i3 + 67; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                  }
                }
              }
            }
          } else {
            for (int ii4 = 0; ii4 < ii3; ii4 += 1) {
              if (_PB_N >= 64 * ii2 + 67) {
                for (int i0 = 64 * ii0 + 1; i0 <= min(_PB_TSTEPS, 64 * ii0 + 32); i0 += 1) {
                  if (i0 >= 64 * ii0 + 2) {
                    if (i0 >= 64 * ii0 + 3 && 64 * ii0 + 32 * ii2 + 2 >= i0) {
                      for (int i2 = 128 * ii0 + 64 * ii2 - 2 * i0 + 4; i2 <= 128 * ii0 + 64 * ii2 - 2 * i0 + 6; i2 += 1) {
                        for (int i3 = 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 5; i3 <= min(64 * ii3 + 1, 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 6); i3 += 1) {
                          for (int i4 = max(1, 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 6); i4 <= min(-128 * ii0 - 64 * ii2 + 64 * ii4 + 2 * i0 + i2 + 58, 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 69); i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (2 * i0 + i2 + i3 == 128 * ii0 + 64 * ii2 + 64 * ii3 + 5) {
                            for (int i4 = -128 * ii0 - 64 * ii2 + 64 * ii4 + 2 * i0 + i2 + 59; i4 <= 85 * ii0 + 43 * ii2 + 64 * ii4 - i0 - i2 + floord(ii0 - ii2 - i0 + i2 + 1, 3) + 67; i4 += 1) {
                              B[i2][128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 5][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 5][i4] - (SCALAR_VAL(2.0) * A[i2][128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 5][i4])) + A[i2 - 1][128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 5][i4])) + (SCALAR_VAL(0.125) * ((A[i2][128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 6][i4] - (SCALAR_VAL(2.0) * A[i2][128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 5][i4])) + A[i2][128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 4][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 5][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 5][i4])) + A[i2][128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 5][i4 - 1]))) + A[i2][128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 5][i4]);
                            }
                          }
                        }
                        for (int i3 = 64 * ii0 + 32 * ii2 + 64 * ii3 - i0 - i2 + (i2 + 1) / 2 + 4; i3 < _PB_N - 1; i3 += 1) {
                          if (64 * ii3 + 2 >= i3) {
                            for (int i4 = max(1, 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 5); i4 <= 64 * ii4; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            if (i0 >= 64 * ii0 + 4 && 2 * i0 + i2 == 128 * ii0 + 64 * ii2 + 6 && i3 == 64 * ii3 + 1) {
                              B[128 * ii0 + 64 * ii2 - 2 * i0 + 6][64 * ii3 + 1][64 * ii4 + 1] = ((((SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 7][64 * ii3 + 1][64 * ii4 + 1] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][64 * ii3 + 1][64 * ii4 + 1])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][64 * ii3 + 1][64 * ii4 + 1])) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][64 * ii3 + 2][64 * ii4 + 1] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][64 * ii3 + 1][64 * ii4 + 1])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][64 * ii3][64 * ii4 + 1]))) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][64 * ii3 + 1][64 * ii4 + 2] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][64 * ii3 + 1][64 * ii4 + 1])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][64 * ii3 + 1][64 * ii4]))) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][64 * ii3 + 1][64 * ii4 + 1]);
                            } else if (i0 == 64 * ii0 + 3 && i2 == 64 * ii2 && i3 == 64 * ii3 + 1) {
                              B[64 * ii2][64 * ii3 + 1][64 * ii4 + 1] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][64 * ii3 + 1][64 * ii4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii2][64 * ii3 + 1][64 * ii4 + 1])) + A[64 * ii2 - 1][64 * ii3 + 1][64 * ii4 + 1])) + (SCALAR_VAL(0.125) * ((A[64 * ii2][64 * ii3 + 2][64 * ii4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii2][64 * ii3 + 1][64 * ii4 + 1])) + A[64 * ii2][64 * ii3][64 * ii4 + 1]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2][64 * ii3 + 1][64 * ii4 + 2] - (SCALAR_VAL(2.0) * A[64 * ii2][64 * ii3 + 1][64 * ii4 + 1])) + A[64 * ii2][64 * ii3 + 1][64 * ii4]))) + A[64 * ii2][64 * ii3 + 1][64 * ii4 + 1]);
                            }
                          } else if (ii4 >= 1 && 2 * i0 + i2 == 128 * ii0 + 64 * ii2 + 6) {
                            B[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][64 * ii4 - 1] = ((((SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 7][i3][64 * ii4 - 1] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][64 * ii4 - 1])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][64 * ii4 - 1])) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3 + 1][64 * ii4 - 1] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][64 * ii4 - 1])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3 - 1][64 * ii4 - 1]))) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][64 * ii4] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][64 * ii4 - 1])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][64 * ii4 - 2]))) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][64 * ii4 - 1]);
                          }
                          for (int i4 = max(max(max(64 * ii4, 63 * ii4 + 1), 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 5), 43 * ii3 + 64 * ii4 - i3 + (-ii3 + i3 - 1) / 3 + 3); i4 <= 64 * ii4 + 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (128 * ii0 + 64 * ii2 + 5 >= 2 * i0 + i2) {
                            for (int i4 = 64 * ii4 + 2; i4 <= 256 * ii0 + 128 * ii2 + 64 * ii4 - 4 * i0 - 2 * i2 + 72; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            if (2 * i0 + i2 == 128 * ii0 + 64 * ii2 + 5) {
                              B[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][64 * ii4 + 63] = ((((SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][64 * ii4 + 63] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][64 * ii4 + 63])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][i3][64 * ii4 + 63])) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3 + 1][64 * ii4 + 63] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][64 * ii4 + 63])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3 - 1][64 * ii4 + 63]))) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][64 * ii4 + 64] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][64 * ii4 + 63])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][64 * ii4 + 62]))) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][64 * ii4 + 63]);
                            }
                          } else {
                            for (int i4 = 64 * ii4 + 2; i4 <= 64 * ii4 + 62; i4 += 1) {
                              B[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 7][i3][i4] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][i4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][i4])) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][i4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][i4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][i4 - 1]))) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][i4]);
                            }
                          }
                        }
                      }
                    }
                    for (int i2 = max(1, 128 * ii0 + 64 * ii2 - 2 * i0 + 7); i2 < 64 * ii2; i2 += 1) {
                      for (int i3 = 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 5; i3 < _PB_N - 1; i3 += 1) {
                        if (64 * ii3 + 2 >= i3) {
                          for (int i4 = max(max(1, 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 5), 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 6); i4 <= min(64 * ii4 + 1, 64 * ii3 + 64 * ii4 - i3 + 2); i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (2 * i0 + i2 + i3 >= 128 * ii0 + 64 * ii2 + 64 * ii3 + 7 && 64 * ii3 + 1 >= i3) {
                            for (int i4 = 64 * ii4 + 2; i4 <= 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 69; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          } else if (128 * ii0 + 64 * ii2 + 64 * ii3 + 6 >= 2 * i0 + i2 + i3) {
                            for (int i4 = 64 * ii4 + 2; i4 <= 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 69; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                        } else {
                          for (int i4 = max(1, 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 5); i4 < 64 * ii4; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                        if (i3 >= 64 * ii3 + 2) {
                          for (int i4 = max(max(64 * ii4, 63 * ii4 + 1), 43 * ii3 + 64 * ii4 - i3 + (-ii3 + i3 - 1) / 3 + 3); i4 <= 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 68; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                    }
                    if (64 * ii3 + 4 >= _PB_N && 64 * ii3 + 2 * i0 >= _PB_N + 128 * ii0 + 3) {
                      for (int i2 = max(max(1, 64 * ii2), 128 * ii0 + 64 * ii2 - 2 * i0 + 7); i2 <= 128 * ii0 + 64 * ii2 - 2 * i0 + 67; i2 += 1) {
                        for (int i3 = max(128 * ii0 + 64 * ii3 - 2 * i0 + 4, 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 5); i3 < _PB_N - 1; i3 += 1) {
                          if (i2 + i3 >= 64 * ii2 + 64 * ii3 + 1 && 64 * ii3 >= i3 + 1) {
                            for (int i4 = max(1, 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 5); i4 <= 64 * ii4 + 1; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            if (64 * ii2 + 2 * i0 + i3 >= 128 * ii0 + 64 * ii3 + i2 + 4) {
                              for (int i4 = 64 * ii4 + 2; i4 <= min(128 * ii0 + 64 * ii4 - 2 * i0 + 67, 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 69); i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                          } else {
                            if (64 * ii2 + 64 * ii3 >= i2 + i3) {
                              for (int i4 = max(max(1, 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 5), 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 6); i4 <= 64 * ii4 + 1; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            } else {
                              for (int i4 = max(max(max(1, 128 * ii0 + 64 * ii4 - 2 * i0 + 4), 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 5), 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 5); i4 <= 64 * ii2 + 64 * ii4 - i2; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                            if (i3 >= 64 * ii3) {
                              for (int i4 = max(max(max(1, 128 * ii0 + 64 * ii4 - 2 * i0 + 4), 64 * ii2 + 64 * ii4 - i2 + 1), 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 5); i4 <= 32 * ii3 + 64 * ii4 - i3 + i3 / 2; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                            if (2 * i0 + i2 + i3 >= 128 * ii0 + 64 * ii2 + 64 * ii3 + 7 && 64 * ii2 + 2 * i0 + i3 >= 128 * ii0 + 64 * ii3 + i2 + 4 && 64 * ii2 + 64 * ii3 >= i2 + i3) {
                              for (int i4 = 64 * ii4 + 2; i4 <= 128 * ii0 + 64 * ii4 - 2 * i0 + 67; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            } else if (i2 + i3 >= 64 * ii2 + 64 * ii3 + 1 && i3 >= 64 * ii3) {
                              for (int i4 = max(max(1, 64 * ii2 + 64 * ii4 - i2 + 1), 32 * ii3 + 64 * ii4 - i3 + i3 / 2 + 1); i4 <= 128 * ii0 + 64 * ii4 - 2 * i0 + 67; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                          }
                          if (2 * i0 + i2 + i3 >= 128 * ii0 + 64 * ii2 + 64 * ii3 + 7 && 64 * ii2 + 2 * i0 + i3 >= 128 * ii0 + 64 * ii3 + i2 + 4 && 64 * ii3 >= i3) {
                            for (int i4 = 128 * ii0 + 64 * ii4 - 2 * i0 + 68; i4 <= 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 69; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                          if (64 * ii2 + 2 * i0 + i3 >= 128 * ii0 + 64 * ii3 + i2 + 4 && 64 * ii3 >= i3 + 1) {
                            for (int i4 = 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 70; i4 <= 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 68; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            if (2 * i0 + i3 >= 128 * ii0 + 64 * ii3 + 5 && 128 * ii0 + 64 * ii2 + 64 * ii3 + 6 >= 2 * i0 + i2 + i3) {
                              for (int i4 = 64 * ii4 + 2; i4 <= 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 69; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                          } else if (64 * ii2 + 2 * i0 >= 128 * ii0 + i2 + 4 && i2 >= 64 * ii2 + 2 && i3 == 64 * ii3) {
                            B[i2][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 68] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 68] - (SCALAR_VAL(2.0) * A[i2][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 68])) + A[i2 - 1][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 68])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 1][128 * ii0 + 64 * ii4 - 2 * i0 + 68] - (SCALAR_VAL(2.0) * A[i2][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 68])) + A[i2][64 * ii3 - 1][128 * ii0 + 64 * ii4 - 2 * i0 + 68]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 69] - (SCALAR_VAL(2.0) * A[i2][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 68])) + A[i2][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 67]))) + A[i2][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 68]);
                          } else if (128 * ii0 + i2 + 3 >= 64 * ii2 + 2 * i0 && i3 == 64 * ii3) {
                            B[i2][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 68] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 68] - (SCALAR_VAL(2.0) * A[i2][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 68])) + A[i2 - 1][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 68])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 1][128 * ii0 + 64 * ii4 - 2 * i0 + 68] - (SCALAR_VAL(2.0) * A[i2][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 68])) + A[i2][64 * ii3 - 1][128 * ii0 + 64 * ii4 - 2 * i0 + 68]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 69] - (SCALAR_VAL(2.0) * A[i2][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 68])) + A[i2][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 67]))) + A[i2][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 68]);
                          } else if (64 * ii3 >= i3 + 1 && 128 * ii0 + 64 * ii3 + i2 + 3 >= 64 * ii2 + 2 * i0 + i3) {
                            for (int i4 = 64 * ii4 + 2; i4 <= 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 68; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          } else if (i2 == 64 * ii2 && i3 >= 64 * ii3 + 1) {
                            B[64 * ii2][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 68] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 68] - (SCALAR_VAL(2.0) * A[64 * ii2][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 68])) + A[64 * ii2 - 1][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 68])) + (SCALAR_VAL(0.125) * ((A[64 * ii2][i3 + 1][128 * ii0 + 64 * ii4 - 2 * i0 + 68] - (SCALAR_VAL(2.0) * A[64 * ii2][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 68])) + A[64 * ii2][i3 - 1][128 * ii0 + 64 * ii4 - 2 * i0 + 68]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 69] - (SCALAR_VAL(2.0) * A[64 * ii2][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 68])) + A[64 * ii2][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 67]))) + A[64 * ii2][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 68]);
                          }
                        }
                      }
                    } else {
                      if (i0 == 64 * ii0 + 2) {
                        for (int i2 = max(1, 64 * ii2); i2 <= 64 * ii2 + 2; i2 += 1) {
                          if (i2 >= 64 * ii2 + 1) {
                            for (int i4 = 64 * ii4 + 1; i4 <= 64 * ii4 + 64; i4 += 1) {
                              B[i2][64 * ii3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii3][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii3][i4])) + A[i2 - 1][64 * ii3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii3][i4])) + A[i2][64 * ii3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii3][i4])) + A[i2][64 * ii3][i4 - 1]))) + A[i2][64 * ii3][i4]);
                            }
                          }
                          for (int i3 = 64 * ii3 + 1; i3 < _PB_N - 1; i3 += 1) {
                            for (int i4 = max(max(1, 64 * ii4), 64 * ii2 + 64 * ii4 - i2 + 1); i4 <= 64 * ii2 + 64 * ii4 - i2 + 64; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            if (i2 == 64 * ii2 + 2) {
                              B[64 * ii2 + 2][i3][64 * ii4 + 63] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 3][i3][64 * ii4 + 63] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][i3][64 * ii4 + 63])) + A[64 * ii2 + 1][i3][64 * ii4 + 63])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][i3 + 1][64 * ii4 + 63] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][i3][64 * ii4 + 63])) + A[64 * ii2 + 2][i3 - 1][64 * ii4 + 63]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][i3][64 * ii4 + 64] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][i3][64 * ii4 + 63])) + A[64 * ii2 + 2][i3][64 * ii4 + 62]))) + A[64 * ii2 + 2][i3][64 * ii4 + 63]);
                            }
                          }
                        }
                      }
                      if (64 * ii3 + 4 >= _PB_N && _PB_N + 128 * ii0 + 2 >= 64 * ii3 + 2 * i0) {
                        for (int i2 = 128 * ii0 + 64 * ii2 - 2 * i0 + 7; i2 <= 128 * ii0 + 64 * ii2 - 2 * i0 + 67; i2 += 1) {
                          if (i0 == 64 * ii0 + 2) {
                            for (int i3 = 64 * ii3; i3 < _PB_N - 1; i3 += 1) {
                              for (int i4 = max(1, 32 * ii3 + 64 * ii4 - i3 + i3 / 2 + 1); i4 <= 64 * ii4 + 63; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                              if (i3 == 64 * ii3) {
                                B[i2][64 * ii3][64 * ii4 + 64] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii3][64 * ii4 + 64] - (SCALAR_VAL(2.0) * A[i2][64 * ii3][64 * ii4 + 64])) + A[i2 - 1][64 * ii3][64 * ii4 + 64])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 1][64 * ii4 + 64] - (SCALAR_VAL(2.0) * A[i2][64 * ii3][64 * ii4 + 64])) + A[i2][64 * ii3 - 1][64 * ii4 + 64]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3][64 * ii4 + 65] - (SCALAR_VAL(2.0) * A[i2][64 * ii3][64 * ii4 + 64])) + A[i2][64 * ii3][64 * ii4 + 63]))) + A[i2][64 * ii3][64 * ii4 + 64]);
                              }
                            }
                          } else {
                            for (int i3 = _PB_N - 6; i3 < min(_PB_N - 1, _PB_N + 64 * ii4 - 2); i3 += 1) {
                              if (_PB_N >= i3 + 3) {
                                if (i3 + 4 >= _PB_N) {
                                  for (int i4 = max(1, _PB_N + 64 * ii4 - i3 - 5); i4 < _PB_N + 64 * ii4 - i3 - 3; i4 += 1) {
                                    B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                                  }
                                  for (int i4 = max(1, _PB_N + 64 * ii4 - i3 - 3); i4 <= 64 * ii4 + 1; i4 += 1) {
                                    B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                                  }
                                  if (64 * ii2 + i3 + 6 >= _PB_N + i2) {
                                    for (int i4 = 64 * ii4 + 2; i4 <= 64 * ii4 + 61; i4 += 1) {
                                      B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                                    }
                                  }
                                } else {
                                  for (int i4 = max(1, _PB_N + 64 * ii4 - i3 - 5); i4 <= 64 * ii4 + 1; i4 += 1) {
                                    B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                                  }
                                }
                                if (64 * ii2 + 2 >= i2 && 64 * ii2 + 60 * i3 + 301 >= 60 * _PB_N + i2) {
                                  for (int i4 = -60 * _PB_N + 64 * ii4 + 60 * i3 + 302; i4 <= _PB_N + 64 * ii4 - i3 + 58; i4 += 1) {
                                    B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                                  }
                                } else if (_PB_N + i2 >= 64 * ii2 + i3 + 7) {
                                  for (int i4 = 64 * ii4 + 2; i4 <= _PB_N + 64 * ii4 - i3 + 58; i4 += 1) {
                                    B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                                  }
                                }
                              } else {
                                for (int i4 = 64 * ii4 - 2; i4 < 64 * ii4; i4 += 1) {
                                  B[i2][_PB_N - 2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][_PB_N - 2][i4] - (SCALAR_VAL(2.0) * A[i2][_PB_N - 2][i4])) + A[i2 - 1][_PB_N - 2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][_PB_N - 1][i4] - (SCALAR_VAL(2.0) * A[i2][_PB_N - 2][i4])) + A[i2][_PB_N - 3][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][_PB_N - 2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][_PB_N - 2][i4])) + A[i2][_PB_N - 2][i4 - 1]))) + A[i2][_PB_N - 2][i4]);
                                }
                                if (64 * ii2 + 3 >= i2) {
                                  for (int i4 = 64 * ii4; i4 <= 64 * ii4 + 1; i4 += 1) {
                                    B[i2][_PB_N - 2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][_PB_N - 2][i4] - (SCALAR_VAL(2.0) * A[i2][_PB_N - 2][i4])) + A[i2 - 1][_PB_N - 2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][_PB_N - 1][i4] - (SCALAR_VAL(2.0) * A[i2][_PB_N - 2][i4])) + A[i2][_PB_N - 3][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][_PB_N - 2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][_PB_N - 2][i4])) + A[i2][_PB_N - 2][i4 - 1]))) + A[i2][_PB_N - 2][i4]);
                                  }
                                  for (int i4 = 64 * ii4 + 2; i4 <= 64 * ii4 + 61; i4 += 1) {
                                    B[i2][_PB_N - 2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][_PB_N - 2][i4] - (SCALAR_VAL(2.0) * A[i2][_PB_N - 2][i4])) + A[i2 - 1][_PB_N - 2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][_PB_N - 1][i4] - (SCALAR_VAL(2.0) * A[i2][_PB_N - 2][i4])) + A[i2][_PB_N - 3][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][_PB_N - 2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][_PB_N - 2][i4])) + A[i2][_PB_N - 2][i4 - 1]))) + A[i2][_PB_N - 2][i4]);
                                  }
                                } else {
                                  B[i2][_PB_N - 2][64 * ii4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][_PB_N - 2][64 * ii4] - (SCALAR_VAL(2.0) * A[i2][_PB_N - 2][64 * ii4])) + A[i2 - 1][_PB_N - 2][64 * ii4])) + (SCALAR_VAL(0.125) * ((A[i2][_PB_N - 1][64 * ii4] - (SCALAR_VAL(2.0) * A[i2][_PB_N - 2][64 * ii4])) + A[i2][_PB_N - 3][64 * ii4]))) + (SCALAR_VAL(0.125) * ((A[i2][_PB_N - 2][64 * ii4 + 1] - (SCALAR_VAL(2.0) * A[i2][_PB_N - 2][64 * ii4])) + A[i2][_PB_N - 2][64 * ii4 - 1]))) + A[i2][_PB_N - 2][64 * ii4]);
                                  for (int i4 = 64 * ii4 + 1; i4 <= 64 * ii4 + 60; i4 += 1) {
                                    B[i2][_PB_N - 2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][_PB_N - 2][i4] - (SCALAR_VAL(2.0) * A[i2][_PB_N - 2][i4])) + A[i2 - 1][_PB_N - 2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][_PB_N - 1][i4] - (SCALAR_VAL(2.0) * A[i2][_PB_N - 2][i4])) + A[i2][_PB_N - 3][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][_PB_N - 2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][_PB_N - 2][i4])) + A[i2][_PB_N - 2][i4 - 1]))) + A[i2][_PB_N - 2][i4]);
                                  }
                                  B[i2][_PB_N - 2][64 * ii4 + 61] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][_PB_N - 2][64 * ii4 + 61] - (SCALAR_VAL(2.0) * A[i2][_PB_N - 2][64 * ii4 + 61])) + A[i2 - 1][_PB_N - 2][64 * ii4 + 61])) + (SCALAR_VAL(0.125) * ((A[i2][_PB_N - 1][64 * ii4 + 61] - (SCALAR_VAL(2.0) * A[i2][_PB_N - 2][64 * ii4 + 61])) + A[i2][_PB_N - 3][64 * ii4 + 61]))) + (SCALAR_VAL(0.125) * ((A[i2][_PB_N - 2][64 * ii4 + 62] - (SCALAR_VAL(2.0) * A[i2][_PB_N - 2][64 * ii4 + 61])) + A[i2][_PB_N - 2][64 * ii4 + 60]))) + A[i2][_PB_N - 2][64 * ii4 + 61]);
                                }
                              }
                            }
                            if (ii4 == 0) {
                              for (int i4 = 1; i4 <= 61; i4 += 1) {
                                B[i2][_PB_N - 2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][_PB_N - 2][i4] - (SCALAR_VAL(2.0) * A[i2][_PB_N - 2][i4])) + A[i2 - 1][_PB_N - 2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][_PB_N - 1][i4] - (SCALAR_VAL(2.0) * A[i2][_PB_N - 2][i4])) + A[i2][_PB_N - 3][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][_PB_N - 2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][_PB_N - 2][i4])) + A[i2][_PB_N - 2][i4 - 1]))) + A[i2][_PB_N - 2][i4]);
                              }
                            }
                          }
                        }
                      } else {
                        for (int i2 = max(max(1, 64 * ii2), 128 * ii0 + 64 * ii2 - 2 * i0 + 7); i2 <= 128 * ii0 + 64 * ii2 - 2 * i0 + 67; i2 += 1) {
                          if (64 * ii2 + 2 * i0 >= 128 * ii0 + i2 + 3) {
                            for (int i3 = max(128 * ii0 + 64 * ii3 - 2 * i0 + 4, 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 5); i3 <= 64 * ii3 + 1; i3 += 1) {
                              if (i2 + i3 >= 64 * ii2 + 64 * ii3 + 1 && i3 >= 64 * ii3) {
                                for (int i4 = max(max(1, 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 5), 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 6); i4 <= 64 * ii2 + 64 * ii4 - i2; i4 += 1) {
                                  B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                                }
                              }
                              if (i3 >= 64 * ii3) {
                                for (int i4 = max(max(1, 64 * ii2 + 64 * ii4 - i2 + 1), 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 5); i4 <= 64 * ii3 + 64 * ii4 - i3; i4 += 1) {
                                  B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                                }
                                if (i2 == 64 * ii2 && i3 == 64 * ii3) {
                                  for (int i4 = max(1, 128 * ii0 + 64 * ii4 - 2 * i0 + 6); i4 <= 64 * ii4 + 1; i4 += 1) {
                                    B[64 * ii2][64 * ii3][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][64 * ii3][i4] - (SCALAR_VAL(2.0) * A[64 * ii2][64 * ii3][i4])) + A[64 * ii2 - 1][64 * ii3][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii2][64 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[64 * ii2][64 * ii3][i4])) + A[64 * ii2][64 * ii3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2][64 * ii3][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii2][64 * ii3][i4])) + A[64 * ii2][64 * ii3][i4 - 1]))) + A[64 * ii2][64 * ii3][i4]);
                                  }
                                }
                                if (i2 == 64 * ii2 && i3 == 64 * ii3) {
                                  for (int i4 = 64 * ii4 + 2; i4 <= 128 * ii0 + 64 * ii4 - 2 * i0 + 68; i4 += 1) {
                                    B[64 * ii2][64 * ii3][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][64 * ii3][i4] - (SCALAR_VAL(2.0) * A[64 * ii2][64 * ii3][i4])) + A[64 * ii2 - 1][64 * ii3][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii2][64 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[64 * ii2][64 * ii3][i4])) + A[64 * ii2][64 * ii3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2][64 * ii3][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii2][64 * ii3][i4])) + A[64 * ii2][64 * ii3][i4 - 1]))) + A[64 * ii2][64 * ii3][i4]);
                                  }
                                } else if (i2 + i3 >= 64 * ii2 + 64 * ii3 + 1) {
                                  for (int i4 = max(max(1, 64 * ii2 + 64 * ii4 - i2 + 1), 64 * ii3 + 64 * ii4 - i3 + 1); i4 <= 64 * ii4 + 1; i4 += 1) {
                                    B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                                  }
                                  if (64 * ii2 + 2 * i0 + i3 >= 128 * ii0 + 64 * ii3 + i2 + 4) {
                                    for (int i4 = 64 * ii4 + 2; i4 <= 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 68; i4 += 1) {
                                      B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                                    }
                                  }
                                }
                              } else {
                                for (int i4 = max(max(1, 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 5), 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 6); i4 <= 64 * ii4 + 1; i4 += 1) {
                                  B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                                }
                                if (2 * i0 + i3 >= 128 * ii0 + 64 * ii3 + 5 && 128 * ii0 + 64 * ii2 + 64 * ii3 + 6 >= 2 * i0 + i2 + i3) {
                                  for (int i4 = 64 * ii4 + 2; i4 <= 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 69; i4 += 1) {
                                    B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                                  }
                                } else if (2 * i0 + i2 + i3 >= 128 * ii0 + 64 * ii2 + 64 * ii3 + 7 && 64 * ii2 + 2 * i0 + i3 >= 128 * ii0 + 64 * ii3 + i2 + 4) {
                                  for (int i4 = 64 * ii4 + 2; i4 <= min(128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 68, 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 69); i4 += 1) {
                                    B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                                  }
                                }
                                if (64 * ii2 + 2 * i0 + i3 >= 128 * ii0 + 64 * ii3 + i2 + 4) {
                                  for (int i4 = 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 70; i4 <= 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 68; i4 += 1) {
                                    B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                                  }
                                }
                              }
                              if (i2 == 64 * ii2 && 2 * i0 + i3 >= 128 * ii0 + 64 * ii3 + 7) {
                                B[64 * ii2][i3][128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 69] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][i3][128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 69] - (SCALAR_VAL(2.0) * A[64 * ii2][i3][128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 69])) + A[64 * ii2 - 1][i3][128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 69])) + (SCALAR_VAL(0.125) * ((A[64 * ii2][i3 + 1][128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 69] - (SCALAR_VAL(2.0) * A[64 * ii2][i3][128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 69])) + A[64 * ii2][i3 - 1][128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 69]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2][i3][128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 70] - (SCALAR_VAL(2.0) * A[64 * ii2][i3][128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 69])) + A[64 * ii2][i3][128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 68]))) + A[64 * ii2][i3][128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 69]);
                              } else if (128 * ii0 + 64 * ii3 + i2 + 3 >= 64 * ii2 + 2 * i0 + i3) {
                                for (int i4 = 64 * ii4 + 2; i4 <= 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 68; i4 += 1) {
                                  B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                                }
                              }
                            }
                            if (ii4 == 0 && 64 * ii2 + 1 >= i2) {
                              for (int i3 = 64 * ii3 + 2; i3 <= min(min(_PB_N - 2, 64 * ii3 + 31), 128 * ii0 + 64 * ii3 - 2 * i0 + 68); i3 += 1) {
                                for (int i4 = 1; i4 <= 128 * ii0 + 64 * ii2 - 2 * i0 - i2 + 68; i4 += 1) {
                                  B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                                }
                              }
                              for (int i3 = 128 * ii0 + 64 * ii3 - 2 * i0 + 69; i3 <= min(_PB_N - 2, 64 * ii3 + 31); i3 += 1) {
                                for (int i4 = 1; i4 <= 128 * ii0 + 64 * ii2 - 2 * i0 - i2 + 68; i4 += 1) {
                                  B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                                }
                              }
                              if (_PB_N >= 64 * ii3 + 34 && i2 == 64 * ii2 + 1) {
                                for (int i4 = 1; i4 <= 128 * ii0 - 2 * i0 + 67; i4 += 1) {
                                  B[64 * ii2 + 1][64 * ii3 + 32][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][64 * ii3 + 32][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][64 * ii3 + 32][i4])) + A[64 * ii2][64 * ii3 + 32][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][64 * ii3 + 33][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][64 * ii3 + 32][i4])) + A[64 * ii2 + 1][64 * ii3 + 31][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][64 * ii3 + 32][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][64 * ii3 + 32][i4])) + A[64 * ii2 + 1][64 * ii3 + 32][i4 - 1]))) + A[64 * ii2 + 1][64 * ii3 + 32][i4]);
                                }
                              }
                              for (int i3 = -64 * ii2 + 64 * ii3 + i2 + 32; i3 < _PB_N - 1; i3 += 1) {
                                for (int i4 = 1; i4 <= 128 * ii0 + 64 * ii2 - 2 * i0 - i2 + 68; i4 += 1) {
                                  B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                                }
                              }
                            } else if (ii4 >= 1 && 64 * ii2 + 1 >= i2) {
                              for (int i3 = 64 * ii3 + 2; i3 < _PB_N - 1; i3 += 1) {
                                if (i3 == 64 * ii3 + 2) {
                                  for (int i4 = 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 5; i4 <= 64 * ii4; i4 += 1) {
                                    B[i2][64 * ii3 + 2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii3 + 2][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 2][i4])) + A[i2 - 1][64 * ii3 + 2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 3][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 2][i4])) + A[i2][64 * ii3 + 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 2][i4])) + A[i2][64 * ii3 + 2][i4 - 1]))) + A[i2][64 * ii3 + 2][i4]);
                                  }
                                } else {
                                  for (int i4 = 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 5; i4 < 64 * ii4; i4 += 1) {
                                    B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                                  }
                                }
                                for (int i4 = max(64 * ii4, 43 * ii3 + 64 * ii4 - i3 + (-ii3 + i3 - 1) / 3 + 3); i4 <= 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 68; i4 += 1) {
                                  B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                                }
                              }
                            } else {
                              for (int i3 = 64 * ii3 + 2; i3 < _PB_N - 1; i3 += 1) {
                                if (i3 == 64 * ii3 + 2) {
                                  for (int i4 = max(1, 128 * ii0 + 64 * ii4 - 2 * i0 + 4); i4 < 64 * ii4; i4 += 1) {
                                    B[i2][64 * ii3 + 2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii3 + 2][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 2][i4])) + A[i2 - 1][64 * ii3 + 2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 3][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 2][i4])) + A[i2][64 * ii3 + 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 2][i4])) + A[i2][64 * ii3 + 2][i4 - 1]))) + A[i2][64 * ii3 + 2][i4]);
                                  }
                                } else {
                                  for (int i4 = max(1, 128 * ii0 + 64 * ii4 - 2 * i0 + 4); i4 < 64 * ii4; i4 += 1) {
                                    B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                                  }
                                }
                                for (int i4 = max(64 * ii4, 63 * ii4 + 1); i4 <= 128 * ii0 + 64 * ii4 - 2 * i0 + 67; i4 += 1) {
                                  B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                                }
                              }
                            }
                          } else if (i0 >= 64 * ii0 + 3) {
                            for (int i3 = 128 * ii0 + 64 * ii3 - 2 * i0 + 4; i3 <= min(_PB_N - 2, 64 * ii3 + 64 * ii4 + 1); i3 += 1) {
                              if (64 * ii3 >= i3 + 1) {
                                for (int i4 = max(1, 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 5); i4 <= 64 * ii4 + 1; i4 += 1) {
                                  B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                                }
                                for (int i4 = 64 * ii4 + 2; i4 <= 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 68; i4 += 1) {
                                  B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                                }
                              } else {
                                if (64 * ii3 + 2 >= i3) {
                                  for (int i4 = max(max(1, 128 * ii0 + 64 * ii4 - 2 * i0 + 4), 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 5); i4 <= 32 * ii3 + 64 * ii4 - i3 + i3 / 2; i4 += 1) {
                                    B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                                  }
                                } else {
                                  for (int i4 = 128 * ii0 + 64 * ii4 - 2 * i0 + 4; i4 < 64 * ii4; i4 += 1) {
                                    B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                                  }
                                }
                                for (int i4 = max(max(1, 64 * ii4), 32 * ii3 + 64 * ii4 - i3 + i3 / 2 + 1); i4 <= 128 * ii0 + 64 * ii4 - 2 * i0 + 67; i4 += 1) {
                                  B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                                }
                                if (i3 == 64 * ii3) {
                                  B[i2][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 68] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 68] - (SCALAR_VAL(2.0) * A[i2][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 68])) + A[i2 - 1][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 68])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 1][128 * ii0 + 64 * ii4 - 2 * i0 + 68] - (SCALAR_VAL(2.0) * A[i2][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 68])) + A[i2][64 * ii3 - 1][128 * ii0 + 64 * ii4 - 2 * i0 + 68]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 69] - (SCALAR_VAL(2.0) * A[i2][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 68])) + A[i2][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 67]))) + A[i2][64 * ii3][128 * ii0 + 64 * ii4 - 2 * i0 + 68]);
                                }
                              }
                            }
                            if (ii4 == 0) {
                              for (int i3 = 64 * ii3 + 2; i3 <= min(_PB_N - 2, 64 * ii3 + 31); i3 += 1) {
                                for (int i4 = 1; i4 <= 128 * ii0 - 2 * i0 + 67; i4 += 1) {
                                  B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                                }
                              }
                              if (_PB_N >= 64 * ii3 + 34) {
                                for (int i4 = 1; i4 <= 128 * ii0 - 2 * i0 + 67; i4 += 1) {
                                  B[i2][64 * ii3 + 32][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii3 + 32][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 32][i4])) + A[i2 - 1][64 * ii3 + 32][i4])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 33][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 32][i4])) + A[i2][64 * ii3 + 31][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 32][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 32][i4])) + A[i2][64 * ii3 + 32][i4 - 1]))) + A[i2][64 * ii3 + 32][i4]);
                                }
                              }
                              for (int i3 = 64 * ii3 + 33; i3 < _PB_N - 1; i3 += 1) {
                                for (int i4 = 1; i4 <= 128 * ii0 - 2 * i0 + 67; i4 += 1) {
                                  B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                                }
                              }
                            }
                          } else {
                            for (int i3 = 64 * ii3; i3 <= min(min(_PB_N - 2, 64 * ii3 + 32), 64 * ii3 + 32 * ii4 + 31); i3 += 1) {
                              for (int i4 = max(max(1, 64 * ii4), 32 * ii3 + 64 * ii4 - i3 + i3 / 2 + 1); i4 <= 64 * ii4 + 63; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                              if (i3 == 64 * ii3) {
                                B[i2][64 * ii3][64 * ii4 + 64] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii3][64 * ii4 + 64] - (SCALAR_VAL(2.0) * A[i2][64 * ii3][64 * ii4 + 64])) + A[i2 - 1][64 * ii3][64 * ii4 + 64])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 1][64 * ii4 + 64] - (SCALAR_VAL(2.0) * A[i2][64 * ii3][64 * ii4 + 64])) + A[i2][64 * ii3 - 1][64 * ii4 + 64]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3][64 * ii4 + 65] - (SCALAR_VAL(2.0) * A[i2][64 * ii3][64 * ii4 + 64])) + A[i2][64 * ii3][64 * ii4 + 63]))) + A[i2][64 * ii3][64 * ii4 + 64]);
                              }
                            }
                            if (_PB_N >= 64 * ii3 + 34 && ii4 == 0) {
                              for (int i4 = 1; i4 <= 33; i4 += 1) {
                                B[i2][64 * ii3 + 32][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii3 + 32][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 32][i4])) + A[i2 - 1][64 * ii3 + 32][i4])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 33][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 32][i4])) + A[i2][64 * ii3 + 31][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 32][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 32][i4])) + A[i2][64 * ii3 + 32][i4 - 1]))) + A[i2][64 * ii3 + 32][i4]);
                              }
                              for (int i4 = 34; i4 <= 63; i4 += 1) {
                                B[i2][64 * ii3 + 32][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii3 + 32][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 32][i4])) + A[i2 - 1][64 * ii3 + 32][i4])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 33][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 32][i4])) + A[i2][64 * ii3 + 31][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 32][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 32][i4])) + A[i2][64 * ii3 + 32][i4 - 1]))) + A[i2][64 * ii3 + 32][i4]);
                              }
                            }
                            for (int i3 = 64 * ii3 + 33; i3 < _PB_N - 1; i3 += 1) {
                              for (int i4 = max(1, 64 * ii4); i4 <= 64 * ii4 + 63; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                  for (int i2 = max(1, 128 * ii0 + 64 * ii2 - 2 * i0 + 3); i2 <= 64 * ii2; i2 += 1) {
                    for (int i3 = 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 4; i3 < _PB_N - 1; i3 += 1) {
                      for (int i4 = max(max(1, 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 4), 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 5); i4 <= 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 68; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                      for (int i4 = max(1, 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 69); i4 <= 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 67; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                  }
                  for (int i2 = 64 * ii2 + 1; i2 <= 128 * ii0 + 64 * ii2 - 2 * i0 + 66; i2 += 1) {
                    for (int i3 = 128 * ii0 + 64 * ii3 - 2 * i0 + 3; i3 < _PB_N - 1; i3 += 1) {
                      if (i3 >= 64 * ii3 + 1) {
                        for (int i4 = max(1, 128 * ii0 + 64 * ii4 - 2 * i0 + 3); i4 <= 128 * ii0 + 64 * ii4 - 2 * i0 + 66; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      } else {
                        for (int i4 = max(1, 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 4); i4 <= 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 67; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                  }
                }
                for (int i0 = 64 * ii0 + 33; i0 <= min(min(_PB_TSTEPS, 64 * ii0 + 64), 64 * ii0 + 32 * ii2 + 33); i0 += 1) {
                  for (int i2 = max(1, 128 * ii0 + 64 * ii2 - 2 * i0 + 4); i2 <= 128 * ii0 + 64 * ii2 - 2 * i0 + 67; i2 += 1) {
                    if (i0 == 64 * ii0 + 33 && i2 == 64 * ii2 + 1) {
                      for (int i4 = 64 * ii4 + 1; i4 <= 64 * ii4 + 64; i4 += 1) {
                        B[64 * ii2 + 1][64 * ii3 - 62][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][64 * ii3 - 62][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][64 * ii3 - 62][i4])) + A[64 * ii2][64 * ii3 - 62][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][64 * ii3 - 61][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][64 * ii3 - 62][i4])) + A[64 * ii2 + 1][64 * ii3 - 63][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][64 * ii3 - 62][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][64 * ii3 - 62][i4])) + A[64 * ii2 + 1][64 * ii3 - 62][i4 - 1]))) + A[64 * ii2 + 1][64 * ii3 - 62][i4]);
                      }
                    }
                    for (int i3 = max(128 * ii0 + 64 * ii3 - 2 * i0 + 5, 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 5); i3 <= min(min(64 * ii3 + 1, 128 * ii0 + 64 * ii3 - 2 * i0 + 124 * i2 - 58), 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 6); i3 += 1) {
                      for (int i4 = max(1, 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 6); i4 <= 64 * ii4 + 1; i4 += 1) {
                        B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                      }
                      for (int i4 = 64 * ii4 + 2; i4 <= min(-128 * ii0 - 64 * ii2 + 64 * ii4 + 2 * i0 + i2 + 58, 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 69); i4 += 1) {
                        B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                      }
                      if (2 * i0 + i2 + i3 == 128 * ii0 + 64 * ii2 + 64 * ii3 + 5) {
                        for (int i4 = -128 * ii0 - 64 * ii2 + 64 * ii4 + 2 * i0 + i2 + 59; i4 <= 85 * ii0 + 43 * ii2 + 64 * ii4 - i0 - i2 + floord(ii0 - ii2 - i0 + i2 + 1, 3) + 67; i4 += 1) {
                          B[i2][128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 5][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 5][i4] - (SCALAR_VAL(2.0) * A[i2][128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 5][i4])) + A[i2 - 1][128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 5][i4])) + (SCALAR_VAL(0.125) * ((A[i2][128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 6][i4] - (SCALAR_VAL(2.0) * A[i2][128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 5][i4])) + A[i2][128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 4][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 5][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 5][i4])) + A[i2][128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 5][i4 - 1]))) + A[i2][128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 5][i4]);
                        }
                      }
                    }
                    for (int i3 = 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 7; i3 <= min(64 * ii2 + 64 * ii3 - i2, -64 * ii2 + 64 * ii3 + i2 + 2); i3 += 1) {
                      for (int i4 = max(1, 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 6); i4 <= 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 69; i4 += 1) {
                        B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                      }
                    }
                    if (ii2 == 1 && i0 >= 64 * ii0 + 35 && i2 == 1) {
                      for (int i3 = 128 * ii0 + 64 * ii3 - 2 * i0 + 68; i3 <= 128 * ii0 + 64 * ii3 - 2 * i0 + 69; i3 += 1) {
                        for (int i4 = max(1, 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 69); i4 <= 128 * ii0 + 64 * ii3 + 64 * ii4 - 2 * i0 - i3 + 132; i4 += 1) {
                          B[1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[2][i3][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[0][i3][i4])) + (SCALAR_VAL(0.125) * ((A[1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3][i4 - 1]))) + A[1][i3][i4]);
                        }
                      }
                    }
                    if (2 * i0 + i2 >= 128 * ii0 + 64 * ii2 + 7) {
                      for (int i3 = max(-64 * ii2 + 64 * ii3 + i2 + 3, 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 7); i3 <= min(min(_PB_N - 2, 64 * ii3 + 2), 64 * ii3 + 64 * ii4 + 1); i3 += 1) {
                        for (int i4 = max(max(1, 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 5), 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 6); i4 <= min(64 * ii4 + 1, 64 * ii3 + 64 * ii4 - i3 + 2); i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        if (64 * ii3 + 1 >= i3) {
                          for (int i4 = 64 * ii4 + 2; i4 <= 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 69; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        } else {
                          for (int i4 = 64 * ii4 + 1; i4 <= 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 68; i4 += 1) {
                            B[i2][64 * ii3 + 2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii3 + 2][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 2][i4])) + A[i2 - 1][64 * ii3 + 2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 3][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 2][i4])) + A[i2][64 * ii3 + 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 2][i4])) + A[i2][64 * ii3 + 2][i4 - 1]))) + A[i2][64 * ii3 + 2][i4]);
                          }
                        }
                      }
                      if (_PB_N >= 64 * ii3 + 4 && ii4 == 0 && 64 * ii2 >= i2 + 1) {
                        for (int i4 = 1; i4 <= 128 * ii0 + 64 * ii2 - 2 * i0 - i2 + 68; i4 += 1) {
                          B[i2][64 * ii3 + 2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii3 + 2][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 2][i4])) + A[i2 - 1][64 * ii3 + 2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 3][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 2][i4])) + A[i2][64 * ii3 + 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii3 + 2][i4])) + A[i2][64 * ii3 + 2][i4 - 1]))) + A[i2][64 * ii3 + 2][i4]);
                        }
                      }
                    }
                    for (int i3 = 64 * ii2 + 64 * ii3 - i2 + 1; i3 <= min(min(_PB_N - 2, 64 * ii3 + 64 * ii4 + 2), -64 * ii2 + 64 * ii3 + i2 + 2); i3 += 1) {
                      for (int i4 = max(max(1, 64 * ii2 + 64 * ii4 - i2 - 61), 64 * ii2 + 64 * ii3 + 64 * ii4 - i2 - i3 - 60); i4 <= 64 * ii3 + 64 * ii4 - i3 + 2; i4 += 1) {
                        B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                      }
                      if (i3 >= 64 * ii3 + 2) {
                        for (int i4 = 64 * ii3 + 64 * ii4 - i3 + 3; i4 <= 64 * ii2 + 64 * ii4 - i2 + 2; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      } else if (i2 == 64 * ii2 && i3 == 64 * ii3 + 1) {
                        B[64 * ii2][64 * ii3 + 1][64 * ii4 + 2] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][64 * ii3 + 1][64 * ii4 + 2] - (SCALAR_VAL(2.0) * A[64 * ii2][64 * ii3 + 1][64 * ii4 + 2])) + A[64 * ii2 - 1][64 * ii3 + 1][64 * ii4 + 2])) + (SCALAR_VAL(0.125) * ((A[64 * ii2][64 * ii3 + 2][64 * ii4 + 2] - (SCALAR_VAL(2.0) * A[64 * ii2][64 * ii3 + 1][64 * ii4 + 2])) + A[64 * ii2][64 * ii3][64 * ii4 + 2]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2][64 * ii3 + 1][64 * ii4 + 3] - (SCALAR_VAL(2.0) * A[64 * ii2][64 * ii3 + 1][64 * ii4 + 2])) + A[64 * ii2][64 * ii3 + 1][64 * ii4 + 1]))) + A[64 * ii2][64 * ii3 + 1][64 * ii4 + 2]);
                      }
                    }
                    if (ii4 >= 1 && 2 * i0 + i2 >= 128 * ii0 + 64 * ii2 + 7) {
                      for (int i3 = max(64 * ii3 + 3, -64 * ii2 + 64 * ii3 + i2 + 3); i3 < _PB_N - 1; i3 += 1) {
                        for (int i4 = 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 5; i4 < 64 * ii4; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = 64 * ii4; i4 <= 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 68; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    } else {
                      if (ii2 == 1 && i0 == 64 * ii0 + 34 && i2 == 1) {
                        for (int i3 = 64 * ii3; i3 <= 64 * ii3 + 1; i3 += 1) {
                          for (int i4 = max(1, 64 * ii3 + 64 * ii4 - i3 + 1); i4 <= 64 * ii3 + 64 * ii4 - i3 + 64; i4 += 1) {
                            B[1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[2][i3][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[0][i3][i4])) + (SCALAR_VAL(0.125) * ((A[1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3][i4 - 1]))) + A[1][i3][i4]);
                          }
                        }
                      }
                      if (_PB_N >= 64 * ii3 + 4 && ii4 == 0 && 2 * i0 + i2 == 128 * ii0 + 64 * ii2 + 5) {
                        for (int i4 = 1; i4 <= 63; i4 += 1) {
                          B[128 * ii0 + 64 * ii2 - 2 * i0 + 5][64 * ii3 + 2][i4] = ((((SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][64 * ii3 + 2][i4] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][64 * ii3 + 2][i4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][64 * ii3 + 2][i4])) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][64 * ii3 + 3][i4] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][64 * ii3 + 2][i4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][64 * ii3 + 1][i4]))) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][64 * ii3 + 2][i4 + 1] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][64 * ii3 + 2][i4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][64 * ii3 + 2][i4 - 1]))) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][64 * ii3 + 2][i4]);
                        }
                      } else if (_PB_N >= 64 * ii3 + 4 && ii4 >= 1 && 2 * i0 + i2 == 128 * ii0 + 64 * ii2 + 5) {
                        for (int i4 = 64 * ii4; i4 <= 64 * ii4 + 63; i4 += 1) {
                          B[128 * ii0 + 64 * ii2 - 2 * i0 + 5][64 * ii3 + 2][i4] = ((((SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][64 * ii3 + 2][i4] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][64 * ii3 + 2][i4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][64 * ii3 + 2][i4])) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][64 * ii3 + 3][i4] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][64 * ii3 + 2][i4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][64 * ii3 + 1][i4]))) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][64 * ii3 + 2][i4 + 1] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][64 * ii3 + 2][i4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][64 * ii3 + 2][i4 - 1]))) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][64 * ii3 + 2][i4]);
                        }
                      }
                      if (128 * ii0 + 64 * ii2 + 5 >= 2 * i0 + i2) {
                        for (int i3 = -128 * ii0 - 64 * ii2 + 64 * ii3 + 2 * i0 + i2 - 2; i3 <= min(_PB_N - 2, 64 * ii3 + 64 * ii4 + 2); i3 += 1) {
                          for (int i4 = 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 5; i4 <= 256 * ii0 + 128 * ii2 + 64 * ii4 - 4 * i0 - 2 * i2 + 72; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (2 * i0 + i2 == 128 * ii0 + 64 * ii2 + 5) {
                            B[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][64 * ii4 + 63] = ((((SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][64 * ii4 + 63] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][64 * ii4 + 63])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][i3][64 * ii4 + 63])) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3 + 1][64 * ii4 + 63] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][64 * ii4 + 63])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3 - 1][64 * ii4 + 63]))) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][64 * ii4 + 64] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][64 * ii4 + 63])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][64 * ii4 + 62]))) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][64 * ii4 + 63]);
                          }
                        }
                        if (ii4 == 0) {
                          for (int i3 = 64 * ii3 + 3; i3 < _PB_N - 1; i3 += 1) {
                            for (int i4 = 1; i4 <= 256 * ii0 + 128 * ii2 - 4 * i0 - 2 * i2 + 72; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            if (2 * i0 + i2 == 128 * ii0 + 64 * ii2 + 5) {
                              B[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][63] = ((((SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][63] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][63])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][i3][63])) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3 + 1][63] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][63])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3 - 1][63]))) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][64] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][63])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][62]))) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][63]);
                            }
                          }
                        }
                      } else {
                        if (ii4 == 0 && 2 * i0 + i2 >= 128 * ii0 + 64 * ii2 + 7) {
                          for (int i3 = 64 * ii3 + 3; i3 < _PB_N - 1; i3 += 1) {
                            for (int i4 = 1; i4 <= 128 * ii0 + 64 * ii2 - 2 * i0 - i2 + 68; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                        }
                        if (2 * i0 + i2 == 128 * ii0 + 64 * ii2 + 6) {
                          for (int i3 = 64 * ii3 + 1; i3 < _PB_N - 1; i3 += 1) {
                            if (64 * ii3 + 2 >= i3) {
                              for (int i4 = max(1, 64 * ii4 - 1); i4 <= 64 * ii3 + 64 * ii4 - i3 + 2; i4 += 1) {
                                B[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 7][i3][i4] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][i4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][i4])) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][i4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][i4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][i4 - 1]))) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][i4]);
                              }
                            } else if (ii4 >= 1) {
                              B[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][64 * ii4 - 1] = ((((SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 7][i3][64 * ii4 - 1] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][64 * ii4 - 1])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][64 * ii4 - 1])) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3 + 1][64 * ii4 - 1] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][64 * ii4 - 1])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3 - 1][64 * ii4 - 1]))) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][64 * ii4] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][64 * ii4 - 1])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][64 * ii4 - 2]))) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][64 * ii4 - 1]);
                            }
                            for (int i4 = max(max(64 * ii4, 63 * ii4 + 1), 43 * ii3 + 64 * ii4 - i3 + (-ii3 + i3 - 1) / 3 + 3); i4 <= 64 * ii4 + 62; i4 += 1) {
                              B[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 7][i3][i4] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][i4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][i3][i4])) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][i4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][i4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][i4 - 1]))) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 6][i3][i4]);
                            }
                          }
                        }
                      }
                    }
                  }
                  for (int i2 = max(1, 128 * ii0 + 64 * ii2 - 2 * i0 + 3); i2 <= 128 * ii0 + 64 * ii2 - 2 * i0 + 66; i2 += 1) {
                    for (int i3 = 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 4; i3 <= 64 * ii3 + 1; i3 += 1) {
                      for (int i4 = max(1, 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 5); i4 <= 128 * ii0 + 64 * ii2 + 64 * ii3 + 64 * ii4 - 2 * i0 - i2 - i3 + 68; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                    for (int i3 = 64 * ii3 + 2; i3 < _PB_N - 1; i3 += 1) {
                      for (int i4 = max(1, 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 4); i4 <= 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 67; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                  }
                }
              } else {
                for (int i2 = 64 * ii2 + 1; i2 < _PB_N - 1; i2 += 1) {
                  for (int i3 = 64 * ii2 + 1; i3 < _PB_N - 1; i3 += 1) {
                    for (int i4 = 64 * ii4 + 1; i4 <= 64 * ii4 + 64; i4 += 1) {
                      A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                    }
                  }
                }
                for (int i0 = 64 * ii0 + 2; i0 <= min(_PB_TSTEPS, 64 * ii0 + 64); i0 += 1) {
                  if (i0 >= 64 * ii0 + 3) {
                    for (int i2 = max(1, 128 * ii0 + 64 * ii2 - 2 * i0 + 4); i2 <= min(64 * ii2 - 1, 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 + 67); i2 += 1) {
                      for (int i3 = 128 * ii0 + 128 * ii2 - 2 * i0 - i2 + 5; i3 < _PB_N - 1; i3 += 1) {
                        if (i3 >= 64 * ii2 + 1) {
                          for (int i4 = max(1, 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 5); i4 <= 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 68; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        } else {
                          for (int i4 = max(1, 128 * ii0 + 128 * ii2 + 64 * ii4 - 2 * i0 - i2 - i3 + 6); i4 <= 128 * ii0 + 128 * ii2 + 64 * ii4 - 2 * i0 - i2 - i3 + 69; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                    }
                    if (ii4 == 0) {
                      for (int i2 = 128 * ii0 + 64 * ii2 - 2 * i0 + 68; i2 < 64 * ii2; i2 += 1) {
                        for (int i3 = max(1, 128 * ii0 + 128 * ii2 - 2 * i0 - i2 + 5); i3 <= 128 * ii0 + 128 * ii2 - 2 * i0 - i2 + 68; i3 += 1) {
                          for (int i4 = 1; i4 <= 128 * ii0 + 128 * ii2 - 2 * i0 - i2 - i3 + 69; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                      if (ii2 == 1 && i0 >= 64 * ii0 + 35) {
                        for (int i2 = 64; i2 <= 65; i2 += 1) {
                          for (int i3 = 1; i3 <= 128 * ii0 - 2 * i0 - i2 + 196; i3 += 1) {
                            for (int i4 = 1; i4 <= 128 * ii0 - 2 * i0 - i2 - i3 + 197; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                        }
                        for (int i2 = 66; i2 < _PB_N - 1; i2 += 1) {
                          for (int i3 = 1; i3 <= 128 * ii0 - 2 * i0 + 131; i3 += 1) {
                            for (int i4 = 1; i4 <= 128 * ii0 - 2 * i0 - i3 + 132; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                        }
                      }
                    }
                    if (64 * ii0 + 32 * ii2 + 32 * ii4 + 2 >= i0) {
                      for (int i2 = 64 * ii2; i2 < _PB_N - 1; i2 += 1) {
                        if (64 * ii0 + 32 * ii2 + 1 >= i0 && i2 >= 64 * ii2 + 1 && 64 * ii2 + 62 >= i2) {
                          for (int i4 = 64 * ii4 + 1; i4 <= 64 * ii4 + 64; i4 += 1) {
                            B[i2][128 * ii0 + 64 * ii2 - 2 * i0 + 4][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][128 * ii0 + 64 * ii2 - 2 * i0 + 4][i4] - (SCALAR_VAL(2.0) * A[i2][128 * ii0 + 64 * ii2 - 2 * i0 + 4][i4])) + A[i2 - 1][128 * ii0 + 64 * ii2 - 2 * i0 + 4][i4])) + (SCALAR_VAL(0.125) * ((A[i2][128 * ii0 + 64 * ii2 - 2 * i0 + 5][i4] - (SCALAR_VAL(2.0) * A[i2][128 * ii0 + 64 * ii2 - 2 * i0 + 4][i4])) + A[i2][128 * ii0 + 64 * ii2 - 2 * i0 + 3][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][128 * ii0 + 64 * ii2 - 2 * i0 + 4][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][128 * ii0 + 64 * ii2 - 2 * i0 + 4][i4])) + A[i2][128 * ii0 + 64 * ii2 - 2 * i0 + 4][i4 - 1]))) + A[i2][128 * ii0 + 64 * ii2 - 2 * i0 + 4][i4]);
                          }
                        } else if (64 * ii0 + 32 * ii2 + 1 >= i0 && i2 >= 64 * ii2 + 63) {
                          for (int i4 = 64 * ii4 + 1; i4 <= 64 * ii4 + 64; i4 += 1) {
                            B[i2][128 * ii0 + 64 * ii2 - 2 * i0 + 4][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][128 * ii0 + 64 * ii2 - 2 * i0 + 4][i4] - (SCALAR_VAL(2.0) * A[i2][128 * ii0 + 64 * ii2 - 2 * i0 + 4][i4])) + A[i2 - 1][128 * ii0 + 64 * ii2 - 2 * i0 + 4][i4])) + (SCALAR_VAL(0.125) * ((A[i2][128 * ii0 + 64 * ii2 - 2 * i0 + 5][i4] - (SCALAR_VAL(2.0) * A[i2][128 * ii0 + 64 * ii2 - 2 * i0 + 4][i4])) + A[i2][128 * ii0 + 64 * ii2 - 2 * i0 + 3][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][128 * ii0 + 64 * ii2 - 2 * i0 + 4][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][128 * ii0 + 64 * ii2 - 2 * i0 + 4][i4])) + A[i2][128 * ii0 + 64 * ii2 - 2 * i0 + 4][i4 - 1]))) + A[i2][128 * ii0 + 64 * ii2 - 2 * i0 + 4][i4]);
                          }
                        }
                        if (64 * ii2 + 62 >= i2) {
                          for (int i3 = 128 * ii0 + 64 * ii2 - 2 * i0 + 5; i3 <= min(64 * ii2, 128 * ii0 + 128 * ii2 + 64 * ii4 - 2 * i0 - i2 + 68); i3 += 1) {
                            if (64 * ii2 >= i3 + 1) {
                              for (int i4 = max(max(1, 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i3 + 5), 128 * ii0 + 128 * ii2 + 64 * ii4 - 2 * i0 - i2 - i3 + 6); i4 <= 128 * ii0 + 128 * ii2 + 64 * ii4 - 2 * i0 - i2 - i3 + 69; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                              for (int i4 = 128 * ii0 + 128 * ii2 + 64 * ii4 - 2 * i0 - i2 - i3 + 70; i4 <= 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i3 + 68; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            } else {
                              for (int i4 = max(max(1, 128 * ii0 + 64 * ii4 - 2 * i0 + 5), 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 6); i4 <= min(64 * ii4, 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 69); i4 += 1) {
                                B[i2][64 * ii2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii2][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii2][i4])) + A[i2 - 1][64 * ii2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii2 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii2][i4])) + A[i2][64 * ii2 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii2][i4])) + A[i2][64 * ii2][i4 - 1]))) + A[i2][64 * ii2][i4]);
                              }
                              if (i2 == 64 * ii2) {
                                for (int i4 = 64 * ii4 + 1; i4 <= 128 * ii0 + 64 * ii4 - 2 * i0 + 69; i4 += 1) {
                                  B[64 * ii2][64 * ii2][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][64 * ii2][i4] - (SCALAR_VAL(2.0) * A[64 * ii2][64 * ii2][i4])) + A[64 * ii2 - 1][64 * ii2][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii2][64 * ii2 + 1][i4] - (SCALAR_VAL(2.0) * A[64 * ii2][64 * ii2][i4])) + A[64 * ii2][64 * ii2 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2][64 * ii2][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii2][64 * ii2][i4])) + A[64 * ii2][64 * ii2][i4 - 1]))) + A[64 * ii2][64 * ii2][i4]);
                                }
                              }
                              for (int i4 = 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 70; i4 <= min(64 * ii4, 128 * ii0 + 64 * ii4 - 2 * i0 + 68); i4 += 1) {
                                B[i2][64 * ii2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii2][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii2][i4])) + A[i2 - 1][64 * ii2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii2 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii2][i4])) + A[i2][64 * ii2 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii2][i4])) + A[i2][64 * ii2][i4 - 1]))) + A[i2][64 * ii2][i4]);
                              }
                              if (i2 >= 64 * ii2 + 1) {
                                for (int i4 = 64 * ii4 + 1; i4 <= 128 * ii0 + 64 * ii4 - 2 * i0 + 68; i4 += 1) {
                                  B[i2][64 * ii2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii2][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii2][i4])) + A[i2 - 1][64 * ii2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii2 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii2][i4])) + A[i2][64 * ii2 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii2][i4])) + A[i2][64 * ii2][i4 - 1]))) + A[i2][64 * ii2][i4]);
                                }
                              }
                            }
                          }
                          for (int i3 = 128 * ii0 + 128 * ii2 + 64 * ii4 - 2 * i0 - i2 + 69; i3 <= min(64 * ii2, 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 + 67); i3 += 1) {
                            for (int i4 = 1; i4 <= 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i3 + 68; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                          if (ii4 == 0 && 2 * i0 + i2 >= 128 * ii0 + 64 * ii2 + 68) {
                            for (int i3 = 64 * ii2 + 1; i3 <= min(min(64 * ii2 + 2, -128 * ii0 + 64 * ii2 + 2 * i0 - 5), 128 * ii0 + 64 * ii2 - 2 * i0 + 68); i3 += 1) {
                              for (int i4 = 1; i4 <= 128 * ii0 - 2 * i0 + 67; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                          }
                          if (ii4 == 1 && 2 * i0 + i2 >= 128 * ii0 + 64 * ii2 + 132) {
                            for (int i3 = 64 * ii2 + 1; i3 < _PB_N - 1; i3 += 1) {
                              for (int i4 = 1; i4 <= 128 * ii0 - 2 * i0 + 131; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                          }
                          if (ii4 == 0 && 2 * i0 + i2 >= 128 * ii0 + 64 * ii2 + 68) {
                            for (int i3 = 64 * ii2 + 3; i3 <= min(-128 * ii0 + 64 * ii2 + 2 * i0 - 5, 128 * ii0 + 64 * ii2 - 2 * i0 + 68); i3 += 1) {
                              for (int i4 = 1; i4 <= 128 * ii0 - 2 * i0 + 67; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                          }
                          if (2 * i0 + i2 >= 128 * ii0 + 64 * ii2 + 7 && 128 * ii0 + 64 * ii2 + 64 * ii4 + 67 >= 2 * i0 + i2) {
                            for (int i3 = 64 * ii2 + 1; i3 <= min(min(min(min(_PB_N - 2, -128 * ii0 + 64 * ii2 + 2 * i0 - 5), 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 + 68), 128 * ii0 + 64 * ii2 - 2 * i0 + 69), 128 * ii0 - 2 * i0 + i2 + 68); i3 += 1) {
                              for (int i4 = max(max(1, 128 * ii0 + 64 * ii4 - 2 * i0 + 4), 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 5); i4 <= 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 68; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                              for (int i4 = 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 69; i4 < 64 * ii4; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                              for (int i4 = max(64 * ii4, 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 69); i4 <= 128 * ii0 + 64 * ii4 - 2 * i0 + 67; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                          }
                        } else {
                          for (int i3 = 128 * ii0 + 64 * ii2 - 2 * i0 + 5; i3 <= min(64 * ii2, 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 + 67); i3 += 1) {
                            if (64 * ii2 >= i3 + 1) {
                              for (int i4 = max(1, 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i3 + 5); i4 <= 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i3 + 68; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            } else {
                              for (int i4 = max(1, 128 * ii0 + 64 * ii4 - 2 * i0 + 5); i4 <= min(64 * ii4, 128 * ii0 + 64 * ii4 - 2 * i0 + 68); i4 += 1) {
                                B[i2][64 * ii2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii2][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii2][i4])) + A[i2 - 1][64 * ii2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii2 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii2][i4])) + A[i2][64 * ii2 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii2][i4])) + A[i2][64 * ii2][i4 - 1]))) + A[i2][64 * ii2][i4]);
                              }
                              for (int i4 = 64 * ii4 + 1; i4 <= 128 * ii0 + 64 * ii4 - 2 * i0 + 68; i4 += 1) {
                                B[i2][64 * ii2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii2][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii2][i4])) + A[i2 - 1][64 * ii2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii2 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii2][i4])) + A[i2][64 * ii2 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii2][i4])) + A[i2][64 * ii2][i4 - 1]))) + A[i2][64 * ii2][i4]);
                              }
                            }
                          }
                          for (int i3 = 64 * ii2 + 1; i3 <= min(min(64 * ii2 + 2, -128 * ii0 + 64 * ii2 + 2 * i0 - 5), 128 * ii0 + 64 * ii2 - 2 * i0 + 69); i3 += 1) {
                            for (int i4 = max(1, 128 * ii0 + 64 * ii4 - 2 * i0 + 4); i4 <= 128 * ii0 + 64 * ii4 - 2 * i0 + 67; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                          if (i0 >= 64 * ii0 + 34) {
                            for (int i3 = max(64 * ii2 + 1, 128 * ii0 + 64 * ii2 - 2 * i0 + 70); i3 < _PB_N - 1; i3 += 1) {
                              for (int i4 = max(1, 128 * ii0 + 64 * ii4 - 2 * i0 + 4); i4 <= 128 * ii0 + 64 * ii4 - 2 * i0 + 67; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                          }
                          for (int i3 = 64 * ii2 + 3; i3 <= min(min(-128 * ii0 + 64 * ii2 + 2 * i0 - 5, 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 + 68), 128 * ii0 + 64 * ii2 - 2 * i0 + 69); i3 += 1) {
                            for (int i4 = max(1, 128 * ii0 + 64 * ii4 - 2 * i0 + 4); i4 <= 128 * ii0 + 64 * ii4 - 2 * i0 + 67; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                        }
                        if (2 * i0 + i2 >= 128 * ii0 + 64 * ii2 + 7) {
                          for (int i3 = -128 * ii0 + 64 * ii2 + 2 * i0 - 4; i3 <= min(min(min(_PB_N - 2, 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 + 68), 128 * ii0 + 64 * ii2 - 2 * i0 + 69), 128 * ii0 - 2 * i0 + i2 + 68); i3 += 1) {
                            for (int i4 = max(max(1, 128 * ii0 + 64 * ii4 - 2 * i0 + 4), 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 5); i4 <= 128 * ii0 + 64 * ii4 - 2 * i0 + 5; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            for (int i4 = max(1, 128 * ii0 + 64 * ii4 - 2 * i0 + 6); i4 <= 64 * ii2 + 64 * ii4 - i2; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            for (int i4 = max(max(1, 128 * ii0 + 64 * ii4 - 2 * i0 + 6), 64 * ii2 + 64 * ii4 - i2 + 1); i4 <= min(64 * ii4 - 1, 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 68); i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            for (int i4 = max(max(1, 128 * ii0 + 64 * ii4 - 2 * i0 + 6), 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 69); i4 < 64 * ii4; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            for (int i4 = max(max(1, 64 * ii4), 64 * ii2 + 64 * ii4 - i2 + 1); i4 <= 128 * ii0 + 64 * ii4 - 2 * i0 + 67; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            if (i2 == 64 * ii2) {
                              B[64 * ii2][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 68] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 68] - (SCALAR_VAL(2.0) * A[64 * ii2][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 68])) + A[64 * ii2 - 1][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 68])) + (SCALAR_VAL(0.125) * ((A[64 * ii2][i3 + 1][128 * ii0 + 64 * ii4 - 2 * i0 + 68] - (SCALAR_VAL(2.0) * A[64 * ii2][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 68])) + A[64 * ii2][i3 - 1][128 * ii0 + 64 * ii4 - 2 * i0 + 68]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 69] - (SCALAR_VAL(2.0) * A[64 * ii2][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 68])) + A[64 * ii2][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 67]))) + A[64 * ii2][i3][128 * ii0 + 64 * ii4 - 2 * i0 + 68]);
                            }
                          }
                          if (ii4 == 0 && _PB_N + 2 * i0 >= 128 * ii0 + 64 * ii2 + 71 && 2 * i0 + i2 >= 128 * ii0 + 64 * ii2 + 68) {
                            for (int i4 = 1; i4 <= 128 * ii0 - 2 * i0 + 67; i4 += 1) {
                              B[i2][128 * ii0 + 64 * ii2 - 2 * i0 + 69][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][128 * ii0 + 64 * ii2 - 2 * i0 + 69][i4] - (SCALAR_VAL(2.0) * A[i2][128 * ii0 + 64 * ii2 - 2 * i0 + 69][i4])) + A[i2 - 1][128 * ii0 + 64 * ii2 - 2 * i0 + 69][i4])) + (SCALAR_VAL(0.125) * ((A[i2][128 * ii0 + 64 * ii2 - 2 * i0 + 70][i4] - (SCALAR_VAL(2.0) * A[i2][128 * ii0 + 64 * ii2 - 2 * i0 + 69][i4])) + A[i2][128 * ii0 + 64 * ii2 - 2 * i0 + 68][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][128 * ii0 + 64 * ii2 - 2 * i0 + 69][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][128 * ii0 + 64 * ii2 - 2 * i0 + 69][i4])) + A[i2][128 * ii0 + 64 * ii2 - 2 * i0 + 69][i4 - 1]))) + A[i2][128 * ii0 + 64 * ii2 - 2 * i0 + 69][i4]);
                            }
                          }
                          if (ii4 == 0 && 2 * i0 + i2 >= 128 * ii0 + 64 * ii2 + 68 && 64 * ii2 + 62 >= i2) {
                            for (int i3 = 128 * ii0 + 64 * ii2 - 2 * i0 + 70; i3 <= min(_PB_N - 2, 128 * ii0 - 2 * i0 + i2 + 68); i3 += 1) {
                              for (int i4 = 1; i4 <= 128 * ii0 - 2 * i0 + 67; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                          } else if (64 * ii0 + 33 >= i0 && i2 >= 64 * ii2 + 63) {
                            for (int i3 = 128 * ii0 + 64 * ii2 - 2 * i0 + 70; i3 < _PB_N - 1; i3 += 1) {
                              for (int i4 = max(1, 128 * ii0 + 64 * ii4 - 2 * i0 + 4); i4 <= 128 * ii0 + 64 * ii4 - 2 * i0 + 67; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                          }
                          if (ii4 == 0 && 2 * i0 + i2 >= 128 * ii0 + 64 * ii2 + 68) {
                            for (int i3 = 128 * ii0 - 2 * i0 + i2 + 69; i3 < _PB_N - 1; i3 += 1) {
                              for (int i4 = 1; i4 <= 128 * ii0 - 2 * i0 + 67; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                          }
                          if (ii4 == 0 && _PB_N + 2 * i0 >= 128 * ii0 + 64 * ii2 + 71 && i2 >= 64 * ii2 + 1 && 128 * ii0 + 64 * ii2 + 67 >= 2 * i0 + i2) {
                            for (int i4 = 1; i4 <= 128 * ii0 - 2 * i0 + 67; i4 += 1) {
                              B[i2][128 * ii0 + 64 * ii2 - 2 * i0 + 69][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][128 * ii0 + 64 * ii2 - 2 * i0 + 69][i4] - (SCALAR_VAL(2.0) * A[i2][128 * ii0 + 64 * ii2 - 2 * i0 + 69][i4])) + A[i2 - 1][128 * ii0 + 64 * ii2 - 2 * i0 + 69][i4])) + (SCALAR_VAL(0.125) * ((A[i2][128 * ii0 + 64 * ii2 - 2 * i0 + 70][i4] - (SCALAR_VAL(2.0) * A[i2][128 * ii0 + 64 * ii2 - 2 * i0 + 69][i4])) + A[i2][128 * ii0 + 64 * ii2 - 2 * i0 + 68][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][128 * ii0 + 64 * ii2 - 2 * i0 + 69][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][128 * ii0 + 64 * ii2 - 2 * i0 + 69][i4])) + A[i2][128 * ii0 + 64 * ii2 - 2 * i0 + 69][i4 - 1]))) + A[i2][128 * ii0 + 64 * ii2 - 2 * i0 + 69][i4]);
                            }
                          }
                        } else {
                          for (int i3 = 64 * ii2 + 1; i3 <= min(_PB_N - 2, 64 * ii2 + 62); i3 += 1) {
                            for (int i4 = max(1, 64 * ii4 - 1); i4 <= 64 * ii4 + 62; i4 += 1) {
                              B[64 * ii2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[64 * ii2][i3][i4])) + A[64 * ii2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[64 * ii2][i3][i4])) + A[64 * ii2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii2][i3][i4])) + A[64 * ii2][i3][i4 - 1]))) + A[64 * ii2][i3][i4]);
                            }
                          }
                        }
                        if (64 * ii2 + 62 >= i2 && 128 * ii0 + 64 * ii2 + 64 * ii4 + 67 >= 2 * i0 + i2) {
                          for (int i3 = max(64 * ii2 + 1, 128 * ii0 + 64 * ii2 - 2 * i0 + 70); i3 <= min(_PB_N - 2, 128 * ii0 - 2 * i0 + i2 + 68); i3 += 1) {
                            for (int i4 = max(1, 128 * ii0 + 64 * ii4 - 2 * i0 + 4); i4 <= 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 68; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            for (int i4 = 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 69; i4 <= min(64 * ii4 - 1, 128 * ii0 + 64 * ii4 - 2 * i0 + 67); i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            for (int i4 = max(64 * ii4, 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 69); i4 <= 128 * ii0 + 64 * ii4 - 2 * i0 + 67; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                          for (int i3 = max(64 * ii2 + 1, 128 * ii0 - 2 * i0 + i2 + 69); i3 < _PB_N - 1; i3 += 1) {
                            for (int i4 = max(max(1, 128 * ii0 + 64 * ii4 - 2 * i0 + 4), 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 5); i4 <= 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 68; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            for (int i4 = 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 69; i4 <= min(64 * ii4 - 1, 128 * ii0 + 64 * ii4 - 2 * i0 + 67); i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            for (int i4 = max(64 * ii4, 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 69); i4 <= 128 * ii0 + 64 * ii4 - 2 * i0 + 67; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                        }
                      }
                    }
                  } else {
                    for (int i2 = 64 * ii2; i2 <= min(_PB_N - 2, 64 * ii2 + 2); i2 += 1) {
                      for (int i3 = max(64 * ii2, 128 * ii2 - i2 + 1); i3 < _PB_N - 1; i3 += 1) {
                        for (int i4 = max(max(max(1, 64 * ii4), 64 * ii2 + 64 * ii4 - i2 + 1), 128 * ii2 + 64 * ii4 - 2 * i3 + 1); i4 <= 64 * ii2 + 64 * ii4 - i2 + 64; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        if (i2 == 64 * ii2 + 2 && i3 >= 64 * ii2 + 1) {
                          B[64 * ii2 + 2][i3][64 * ii4 + 63] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 3][i3][64 * ii4 + 63] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][i3][64 * ii4 + 63])) + A[64 * ii2 + 1][i3][64 * ii4 + 63])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][i3 + 1][64 * ii4 + 63] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][i3][64 * ii4 + 63])) + A[64 * ii2 + 2][i3 - 1][64 * ii4 + 63]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][i3][64 * ii4 + 64] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][i3][64 * ii4 + 63])) + A[64 * ii2 + 2][i3][64 * ii4 + 62]))) + A[64 * ii2 + 2][i3][64 * ii4 + 63]);
                        } else if (i3 == 64 * ii2) {
                          for (int i4 = 64 * ii2 + 64 * ii4 - i2 + 65; i4 <= 64 * ii4 + 64; i4 += 1) {
                            B[i2][64 * ii2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii2][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii2][i4])) + A[i2 - 1][64 * ii2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii2 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii2][i4])) + A[i2][64 * ii2 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii2][i4])) + A[i2][64 * ii2][i4 - 1]))) + A[i2][64 * ii2][i4]);
                          }
                        }
                      }
                    }
                    for (int i2 = 64 * ii2 + 3; i2 < _PB_N - 1; i2 += 1) {
                      for (int i3 = 64 * ii2; i3 < _PB_N - 1; i3 += 1) {
                        for (int i4 = max(max(1, 64 * ii4), 128 * ii2 + 64 * ii4 - 2 * i3 + 1); i4 <= 128 * ii2 + 64 * ii4 - 2 * i3 + 64; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = max(max(1, 64 * ii4), 128 * ii2 + 64 * ii4 - 2 * i3 + 65); i4 <= 64 * ii4 + 63; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    }
                  }
                  for (int i2 = max(1, 128 * ii0 + 64 * ii2 - 2 * i0 + 3); i2 <= 64 * ii2; i2 += 1) {
                    for (int i3 = max(1, 128 * ii0 + 128 * ii2 - 2 * i0 - i2 + 4); i3 <= min(64 * ii2 + 1, 128 * ii0 + 128 * ii2 + 64 * ii4 - 2 * i0 - i2 + 67); i3 += 1) {
                      for (int i4 = max(1, 128 * ii0 + 128 * ii2 + 64 * ii4 - 2 * i0 - i2 - i3 + 5); i4 <= 128 * ii0 + 128 * ii2 + 64 * ii4 - 2 * i0 - i2 - i3 + 68; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                    for (int i3 = 64 * ii2 + 2; i3 < _PB_N - 1; i3 += 1) {
                      for (int i4 = max(1, 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 4); i4 <= 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i2 + 67; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                  }
                  for (int i2 = 64 * ii2 + 1; i2 < _PB_N - 1; i2 += 1) {
                    for (int i3 = max(1, 128 * ii0 + 64 * ii2 - 2 * i0 + 3); i3 <= min(64 * ii2, 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 + 66); i3 += 1) {
                      for (int i4 = max(1, 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i3 + 4); i4 <= 128 * ii0 + 64 * ii2 + 64 * ii4 - 2 * i0 - i3 + 67; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                    for (int i3 = 64 * ii2 + 1; i3 < _PB_N - 1; i3 += 1) {
                      for (int i4 = max(1, 128 * ii0 + 64 * ii4 - 2 * i0 + 3); i4 <= 128 * ii0 + 64 * ii4 - 2 * i0 + 66; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                  }
                }
              }
            }
            if (_PB_N >= 64 * ii2 + 67) {
              for (int i0 = 64 * ii0 + 1; i0 <= min(min(_PB_TSTEPS, 64 * ii0 + 64), 64 * ii0 + 32 * ii2 + 33); i0 += 1) {
                if (i0 >= 64 * ii0 + 2) {
                  if (i0 >= 64 * ii0 + 3) {
                    for (int i2 = max(1, 128 * ii0 + 64 * ii2 - 2 * i0 + 4); i2 <= 128 * ii0 + 64 * ii2 - 2 * i0 + 67; i2 += 1) {
                      if (_PB_N + 2 * i0 + i2 >= 128 * ii0 + 64 * ii2 + 64 * ii3 + 8) {
                        if (i2 == 64 * ii2 + 1) {
                          for (int i4 = 64 * ii3 + 1; i4 < _PB_N - 1; i4 += 1) {
                            B[64 * ii2 + 1][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4])) + A[64 * ii2][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][128 * ii0 + 64 * ii3 - 2 * i0 + 5][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4])) + A[64 * ii2 + 1][128 * ii0 + 64 * ii3 - 2 * i0 + 3][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4])) + A[64 * ii2 + 1][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4 - 1]))) + A[64 * ii2 + 1][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4]);
                          }
                        }
                        for (int i3 = max(128 * ii0 + 64 * ii3 - 2 * i0 + 5, 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 5); i3 <= -_PB_N + 128 * ii0 + 64 * ii2 + 128 * ii3 - 2 * i0 - i2 + 9; i3 += 1) {
                          if (_PB_N + 2 * i0 + i2 == 128 * ii0 + 64 * ii2 + 64 * ii3 + 8 && i3 == 64 * ii3 + 1) {
                            B[-_PB_N + 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 + 8][64 * ii3 + 1][_PB_N - 3] = ((((SCALAR_VAL(0.125) * ((A[-_PB_N + 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 + 9][64 * ii3 + 1][_PB_N - 3] - (SCALAR_VAL(2.0) * A[-_PB_N + 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 + 8][64 * ii3 + 1][_PB_N - 3])) + A[-_PB_N + 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 + 7][64 * ii3 + 1][_PB_N - 3])) + (SCALAR_VAL(0.125) * ((A[-_PB_N + 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 + 8][64 * ii3 + 2][_PB_N - 3] - (SCALAR_VAL(2.0) * A[-_PB_N + 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 + 8][64 * ii3 + 1][_PB_N - 3])) + A[-_PB_N + 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 + 8][64 * ii3][_PB_N - 3]))) + (SCALAR_VAL(0.125) * ((A[-_PB_N + 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 + 8][64 * ii3 + 1][_PB_N - 2] - (SCALAR_VAL(2.0) * A[-_PB_N + 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 + 8][64 * ii3 + 1][_PB_N - 3])) + A[-_PB_N + 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 + 8][64 * ii3 + 1][_PB_N - 4]))) + A[-_PB_N + 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 + 8][64 * ii3 + 1][_PB_N - 3]);
                          }
                          for (int i4 = max(128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 6, 128 * ii0 + 64 * ii2 + 128 * ii3 - 2 * i0 - i2 - i3 + 6); i4 < _PB_N - 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                        if (64 * ii3 + 4 == _PB_N && 2 * i0 + i2 == 128 * ii0 + 64 * ii2 + 4) {
                          for (int i4 = _PB_N - 3; i4 < _PB_N - 1; i4 += 1) {
                            B[128 * ii0 + 64 * ii2 - 2 * i0 + 4][_PB_N - 2][i4] = ((((SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][_PB_N - 2][i4] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][_PB_N - 2][i4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 3][_PB_N - 2][i4])) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][_PB_N - 1][i4] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][_PB_N - 2][i4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][_PB_N - 3][i4]))) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][_PB_N - 2][i4 + 1] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][_PB_N - 2][i4])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][_PB_N - 2][i4 - 1]))) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][_PB_N - 2][i4]);
                          }
                        } else {
                          if (i2 >= 64 * ii2 + 2) {
                            for (int i4 = 64 * ii3 + 1; i4 < _PB_N - 1; i4 += 1) {
                              B[i2][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4] - (SCALAR_VAL(2.0) * A[i2][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4])) + A[i2 - 1][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4])) + (SCALAR_VAL(0.125) * ((A[i2][128 * ii0 + 64 * ii3 - 2 * i0 + 5][i4] - (SCALAR_VAL(2.0) * A[i2][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4])) + A[i2][128 * ii0 + 64 * ii3 - 2 * i0 + 3][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4])) + A[i2][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4 - 1]))) + A[i2][128 * ii0 + 64 * ii3 - 2 * i0 + 4][i4]);
                            }
                          }
                          if (_PB_N + 2 * i0 + i2 >= 128 * ii0 + 64 * ii2 + 64 * ii3 + 9) {
                            for (int i3 = max(max(128 * ii0 + 64 * ii3 - 2 * i0 + 5, 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 5), -_PB_N + 128 * ii0 + 64 * ii2 + 128 * ii3 - 2 * i0 - i2 + 10); i3 < _PB_N - 1; i3 += 1) {
                              if (64 * ii3 + 2 * i2 >= 128 * ii2 + i3 + 2 && i3 >= 64 * ii3 + 2 && 64 * ii3 + 4 >= i3) {
                                B[i2][i3][128 * ii0 + 64 * ii3 - 2 * i0 + 4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][128 * ii0 + 64 * ii3 - 2 * i0 + 4] - (SCALAR_VAL(2.0) * A[i2][i3][128 * ii0 + 64 * ii3 - 2 * i0 + 4])) + A[i2 - 1][i3][128 * ii0 + 64 * ii3 - 2 * i0 + 4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][128 * ii0 + 64 * ii3 - 2 * i0 + 4] - (SCALAR_VAL(2.0) * A[i2][i3][128 * ii0 + 64 * ii3 - 2 * i0 + 4])) + A[i2][i3 - 1][128 * ii0 + 64 * ii3 - 2 * i0 + 4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][128 * ii0 + 64 * ii3 - 2 * i0 + 5] - (SCALAR_VAL(2.0) * A[i2][i3][128 * ii0 + 64 * ii3 - 2 * i0 + 4])) + A[i2][i3][128 * ii0 + 64 * ii3 - 2 * i0 + 3]))) + A[i2][i3][128 * ii0 + 64 * ii3 - 2 * i0 + 4]);
                              } else {
                                if (64 * ii3 + 1 >= i3) {
                                  for (int i4 = 128 * ii0 + 128 * ii3 - 2 * i0 - i3 + 5; i4 <= min(128 * ii0 + 64 * ii3 - 2 * i0 + 5, 128 * ii0 - 128 * ii2 + 64 * ii3 - 2 * i0 + 2 * i2 + 1); i4 += 1) {
                                    B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                                  }
                                }
                                if (128 * ii2 + i3 + 1 >= 64 * ii3 + 2 * i2) {
                                  for (int i4 = max(128 * ii0 + 64 * ii3 - 2 * i0 + 4, 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 5); i4 <= min(128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 6, 128 * ii0 + 64 * ii2 - 128 * ii3 - 2 * i0 - i2 + 3 * i3 + 2); i4 += 1) {
                                    B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                                  }
                                } else if (i2 >= 64 * ii2 + 2 && 64 * ii3 >= i3 + 1) {
                                  B[i2][i3][128 * ii0 + 128 * ii3 - 2 * i0 - i3 + 5] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][128 * ii0 + 128 * ii3 - 2 * i0 - i3 + 5] - (SCALAR_VAL(2.0) * A[i2][i3][128 * ii0 + 128 * ii3 - 2 * i0 - i3 + 5])) + A[i2 - 1][i3][128 * ii0 + 128 * ii3 - 2 * i0 - i3 + 5])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][128 * ii0 + 128 * ii3 - 2 * i0 - i3 + 5] - (SCALAR_VAL(2.0) * A[i2][i3][128 * ii0 + 128 * ii3 - 2 * i0 - i3 + 5])) + A[i2][i3 - 1][128 * ii0 + 128 * ii3 - 2 * i0 - i3 + 5]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][128 * ii0 + 128 * ii3 - 2 * i0 - i3 + 6] - (SCALAR_VAL(2.0) * A[i2][i3][128 * ii0 + 128 * ii3 - 2 * i0 - i3 + 5])) + A[i2][i3][128 * ii0 + 128 * ii3 - 2 * i0 - i3 + 4]))) + A[i2][i3][128 * ii0 + 128 * ii3 - 2 * i0 - i3 + 5]);
                                }
                              }
                              for (int i4 = max(max(128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 6, 128 * ii0 - 128 * ii2 + 128 * ii3 - 2 * i0 + 2 * i2 - i3 + 3), 128 * ii0 + 64 * ii2 + 128 * ii3 - 2 * i0 - i2 - i3 + 6); i4 <= 128 * ii0 + 64 * ii2 + 128 * ii3 - 2 * i0 - i2 - i3 + 7; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                              for (int i4 = max(max(max(max(128 * ii0 + 64 * ii3 - 2 * i0 + 4, 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 7), 128 * ii0 + 128 * ii3 - 2 * i0 - i3 + 6), 128 * ii0 + 64 * ii2 + 128 * ii3 - 2 * i0 - i2 - i3 + 8), 128 * ii0 + 85 * ii3 - 2 * i0 - (-ii3 + i3 + 1) / 3 + 6); i4 < _PB_N - 1; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                          }
                        }
                      } else {
                        B[128 * ii0 + 64 * ii2 - 2 * i0 + 4][_PB_N - 2][_PB_N - 2] = ((((SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 5][_PB_N - 2][_PB_N - 2] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][_PB_N - 2][_PB_N - 2])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 3][_PB_N - 2][_PB_N - 2])) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][_PB_N - 1][_PB_N - 2] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][_PB_N - 2][_PB_N - 2])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][_PB_N - 3][_PB_N - 2]))) + (SCALAR_VAL(0.125) * ((A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][_PB_N - 2][_PB_N - 1] - (SCALAR_VAL(2.0) * A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][_PB_N - 2][_PB_N - 2])) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][_PB_N - 2][_PB_N - 3]))) + A[128 * ii0 + 64 * ii2 - 2 * i0 + 4][_PB_N - 2][_PB_N - 2]);
                      }
                    }
                  } else {
                    for (int i2 = max(1, 64 * ii2); i2 <= 64 * ii2 + 63; i2 += 1) {
                      for (int i3 = max(64 * ii3, 64 * ii2 + 64 * ii3 - i2 + 1); i3 < _PB_N - 1; i3 += 1) {
                        for (int i4 = max(max(64 * ii3, 64 * ii2 + 64 * ii3 - i2 + 1), 192 * ii3 - 2 * i3 + 1); i4 < _PB_N - 1; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    }
                  }
                }
                if (64 * ii0 + 32 * ii2 + 32 >= i0) {
                  for (int i2 = max(1, 128 * ii0 + 64 * ii2 - 2 * i0 + 3); i2 <= min(64 * ii2, 128 * ii0 + 64 * ii2 - 2 * i0 + 66); i2 += 1) {
                    for (int i3 = 128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 4; i3 < _PB_N - 1; i3 += 1) {
                      for (int i4 = max(128 * ii0 + 64 * ii2 + 64 * ii3 - 2 * i0 - i2 + 4, 128 * ii0 + 64 * ii2 + 128 * ii3 - 2 * i0 - i2 - i3 + 5); i4 < _PB_N - 1; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                  }
                  for (int i2 = 64 * ii2 + 1; i2 <= 128 * ii0 + 64 * ii2 - 2 * i0 + 66; i2 += 1) {
                    for (int i3 = 128 * ii0 + 64 * ii3 - 2 * i0 + 3; i3 <= 64 * ii3; i3 += 1) {
                      for (int i4 = 128 * ii0 + 128 * ii3 - 2 * i0 - i3 + 4; i4 < _PB_N - 1; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                    for (int i3 = 64 * ii3 + 1; i3 < _PB_N - 1; i3 += 1) {
                      for (int i4 = 128 * ii0 + 64 * ii3 - 2 * i0 + 3; i4 < _PB_N - 1; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                  }
                }
              }
            } else {
              for (int i0 = 64 * ii0 + 1; i0 <= min(_PB_TSTEPS, 64 * ii0 + 64); i0 += 1) {
                if (i0 >= 64 * ii0 + 2) {
                  for (int i2 = max(1, 128 * ii0 + 64 * ii2 - 2 * i0 + 4); i2 < _PB_N - 1; i2 += 1) {
                    if (64 * ii2 + 3 == _PB_N && i2 + 2 == _PB_N) {
                      for (int i3 = _PB_N + 128 * ii0 - 2 * i0 + 1; i3 <= min(_PB_N - 4, _PB_N + 128 * ii0 - 2 * i0 + 3); i3 += 1) {
                        for (int i4 = 2 * _PB_N + 128 * ii0 - 2 * i0 - i3 - 1; i4 < _PB_N - 1; i4 += 1) {
                          B[_PB_N - 2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[_PB_N - 1][i3][i4] - (SCALAR_VAL(2.0) * A[_PB_N - 2][i3][i4])) + A[_PB_N - 3][i3][i4])) + (SCALAR_VAL(0.125) * ((A[_PB_N - 2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[_PB_N - 2][i3][i4])) + A[_PB_N - 2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[_PB_N - 2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[_PB_N - 2][i3][i4])) + A[_PB_N - 2][i3][i4 - 1]))) + A[_PB_N - 2][i3][i4]);
                        }
                      }
                    } else if (64 * ii2 + 1 >= i2) {
                      for (int i3 = max(1, 128 * ii0 + 128 * ii2 - 2 * i0 - i2 + 5); i3 <= min(min(64 * ii2 - 1, 128 * ii0 + 171 * ii2 - 2 * i0 - 2 * i2 + (-ii2 + i2 + 1) / 3 + 8), -43 * ii2 + 2 * i2 - (-ii2 + i2 + 1) / 3 - 1); i3 += 1) {
                        for (int i4 = 128 * ii0 + 192 * ii2 - 2 * i0 - i2 - i3 + 6; i4 < _PB_N - 1; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    }
                    if (i0 == 64 * ii0 + 3 && i2 == 64 * ii2 + 1) {
                      for (int i4 = 64 * ii2 - 1; i4 < _PB_N - 1; i4 += 1) {
                        B[64 * ii2 + 1][64 * ii2][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][64 * ii2][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][64 * ii2][i4])) + A[64 * ii2][64 * ii2][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][64 * ii2 + 1][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][64 * ii2][i4])) + A[64 * ii2 + 1][64 * ii2 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][64 * ii2][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][64 * ii2][i4])) + A[64 * ii2 + 1][64 * ii2][i4 - 1]))) + A[64 * ii2 + 1][64 * ii2][i4]);
                      }
                    }
                    if (ii2 == 1) {
                      for (int i3 = max(max(1, 128 * ii0 - 2 * i0 - i2 + 133), 2 * i2 - (i2 + 3) / 3 - 42); i3 <= min(-128 * ii0 + 2 * i0 + 2 * i2 - (i2 + 3) / 3 - 110, 128 * ii0 - 2 * i0 - 2 * i2 + i2 / 3 + 179); i3 += 1) {
                        for (int i4 = 128 * ii0 - 2 * i0 - i2 - i3 + 198; i4 < _PB_N - 1; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    }
                    for (int i3 = max(max(128 * ii0 + 128 * ii2 - 2 * i0 - i2 + 5, -128 * ii0 - 107 * ii2 + 2 * i0 + 2 * i2 - (-ii2 + i2 + 1) / 3 - 3), -43 * ii2 + 2 * i2 - (-ii2 + i2 + 1) / 3); i3 <= min(43 * ii2 + (-ii2 + i2 + 2) / 3, 128 * ii0 + 171 * ii2 - 2 * i0 - 2 * i2 + (-ii2 + i2 + 1) / 3 + 8); i3 += 1) {
                      for (int i4 = 128 * ii0 + 192 * ii2 - 2 * i0 - i2 - i3 + 6; i4 <= 21 * ii2 + i2 - (-ii2 + i2 + 1) / 3; i4 += 1) {
                        B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                      }
                      for (int i4 = max(128 * ii0 + 192 * ii2 - 2 * i0 - i2 - i3 + 6, 21 * ii2 + i2 - (-ii2 + i2 + 1) / 3 + 1); i4 < _PB_N - 1; i4 += 1) {
                        B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                      }
                    }
                    if (i0 >= 64 * ii0 + 3 && i2 >= 64 * ii2 + 2) {
                      for (int i3 = max(1, 128 * ii0 + 64 * ii2 - 2 * i0 + 4); i3 <= min(57 * ii2 + (-ii2 + i2 - 2) / 9 + 1, 128 * ii0 + 121 * ii2 - 2 * i0 + (-ii2 + i2 - 2) / 9 + 4); i3 += 1) {
                        for (int i4 = max(max(1, 128 * ii0 + 64 * ii2 - 2 * i0 + 4), 128 * ii0 + 128 * ii2 - 2 * i0 - i3 + 5); i4 <= 128 * ii0 + 121 * ii2 - 2 * i0 - i3 + (-ii2 + i2 + 7) / 9 + 4; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        if (i2 == 64 * ii2 + 2 && 2 * i0 + i3 >= 128 * ii0 + 64 * ii2 + 5) {
                          B[64 * ii2 + 2][i3][128 * ii0 + 128 * ii2 - 2 * i0 - i3 + 6] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 3][i3][128 * ii0 + 128 * ii2 - 2 * i0 - i3 + 6] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][i3][128 * ii0 + 128 * ii2 - 2 * i0 - i3 + 6])) + A[64 * ii2 + 1][i3][128 * ii0 + 128 * ii2 - 2 * i0 - i3 + 6])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][i3 + 1][128 * ii0 + 128 * ii2 - 2 * i0 - i3 + 6] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][i3][128 * ii0 + 128 * ii2 - 2 * i0 - i3 + 6])) + A[64 * ii2 + 2][i3 - 1][128 * ii0 + 128 * ii2 - 2 * i0 - i3 + 6]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][i3][128 * ii0 + 128 * ii2 - 2 * i0 - i3 + 7] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][i3][128 * ii0 + 128 * ii2 - 2 * i0 - i3 + 6])) + A[64 * ii2 + 2][i3][128 * ii0 + 128 * ii2 - 2 * i0 - i3 + 5]))) + A[64 * ii2 + 2][i3][128 * ii0 + 128 * ii2 - 2 * i0 - i3 + 6]);
                        }
                        for (int i4 = 128 * ii0 + 121 * ii2 - 2 * i0 - i3 + (-ii2 + i2 + 7) / 9 + 5; i4 <= min(192 * ii2 - 2 * i3, 128 * ii0 - 21 * ii2 - 2 * i0 + i3 + (-ii2 + i2 + 1) / 3 + 2); i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = max(192 * ii2 - 2 * i3 + 1, 128 * ii0 + 121 * ii2 - 2 * i0 - i3 + (-ii2 + i2 + 7) / 9 + 5); i4 <= min(64 * ii2 - 1, 128 * ii0 - 21 * ii2 - 2 * i0 + i3 + (-ii2 + i2 + 1) / 3 + 2); i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = max(max(128 * ii0 + 192 * ii2 - 2 * i0 - i2 - i3 + 9, 128 * ii0 - 21 * ii2 - 2 * i0 + i3 + (-ii2 + i2 + 1) / 3 + 3), 128 * ii0 + 121 * ii2 - 2 * i0 - i3 + (-ii2 + i2 + 7) / 9 + 5); i4 <= min(128 * ii0 + 128 * ii2 - 2 * i0 - i3 + 7, 21 * ii2 + i2 - (-ii2 + i2 + 1) / 3); i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = max(max(128 * ii0 + 128 * ii2 - 2 * i0 - i3 + 8, 128 * ii0 - 21 * ii2 - 2 * i0 + i3 + (-ii2 + i2 + 1) / 3 + 3), 128 * ii0 + 121 * ii2 - 2 * i0 - i3 + (-ii2 + i2 - 2) / 9 + 6); i4 <= min(min(_PB_N - 2, 192 * ii2 - 2 * i3), -64 * ii0 - 64 * ii2 + i0 + i2 + (i2 + i3) / 2 - 4); i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        if (64 * ii2 >= i3 + 1) {
                          for (int i4 = max(192 * ii2 - 2 * i3 + 1, 128 * ii0 + 121 * ii2 - 2 * i0 - i3 + (-ii2 + i2 + 7) / 9 + 5); i4 <= 128 * ii0 - 21 * ii2 - 2 * i0 + i3 + (-ii2 + i2 + 1) / 3 + 2; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          for (int i4 = max(192 * ii2 - 2 * i3 + 1, 128 * ii0 - 21 * ii2 - 2 * i0 + i3 + (-ii2 + i2 + 1) / 3 + 3); i4 < min(_PB_N - 1, -64 * ii0 - 64 * ii2 + i0 + i2 + (i2 + i3) / 2 - 3); i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                        for (int i4 = max(max(192 * ii2 - 2 * i3 + 1, 128 * ii0 + 128 * ii2 - 2 * i0 - i3 + 8), 128 * ii0 - 21 * ii2 - 2 * i0 + i3 + (-ii2 + i2 + 1) / 3 + 3); i4 < 64 * ii2; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        if (64 * ii2 >= i3 + 1) {
                          for (int i4 = -64 * ii0 - 64 * ii2 + i0 + i2 + (i2 + i3) / 2 - 3; i4 < _PB_N - 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        } else {
                          if (i0 == 64 * ii0 + 3 && 64 * ii2 + 19 >= i2 && 64 * ii2 + 1 >= i3 && i2 + 6 * i3 >= 448 * ii2 + 14) {
                            B[i2][i3][128 * ii2 - i3 + 1] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][128 * ii2 - i3 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][128 * ii2 - i3 + 1])) + A[i2 - 1][i3][128 * ii2 - i3 + 1])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][128 * ii2 - i3 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][128 * ii2 - i3 + 1])) + A[i2][i3 - 1][128 * ii2 - i3 + 1]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][128 * ii2 - i3 + 2] - (SCALAR_VAL(2.0) * A[i2][i3][128 * ii2 - i3 + 1])) + A[i2][i3][128 * ii2 - i3]))) + A[i2][i3][128 * ii2 - i3 + 1]);
                          }
                          for (int i4 = max(max(max(64 * ii2, 192 * ii2 - 2 * i3 + 1), 128 * ii0 + 128 * ii2 - 2 * i0 - i3 + 8), 128 * ii0 + 121 * ii2 - 2 * i0 - i3 + (-ii2 + i2 - 2) / 9 + 6); i4 < _PB_N - 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                    }
                    if (ii2 == 1 && i2 >= 66) {
                      for (int i3 = 128 * ii0 - 2 * i0 + (i2 - 3) / 9 + 126; i3 <= (i2 - 3) / 9 + 58; i3 += 1) {
                        if (i2 == 66 && 2 * i0 + i3 == 128 * ii0 + 133) {
                          B[66][128 * ii0 - 2 * i0 + 133][1] = ((((SCALAR_VAL(0.125) * ((A[67][128 * ii0 - 2 * i0 + 133][1] - (SCALAR_VAL(2.0) * A[66][128 * ii0 - 2 * i0 + 133][1])) + A[65][128 * ii0 - 2 * i0 + 133][1])) + (SCALAR_VAL(0.125) * ((A[66][128 * ii0 - 2 * i0 + 134][1] - (SCALAR_VAL(2.0) * A[66][128 * ii0 - 2 * i0 + 133][1])) + A[66][128 * ii0 - 2 * i0 + 132][1]))) + (SCALAR_VAL(0.125) * ((A[66][128 * ii0 - 2 * i0 + 133][2] - (SCALAR_VAL(2.0) * A[66][128 * ii0 - 2 * i0 + 133][1])) + A[66][128 * ii0 - 2 * i0 + 133][0]))) + A[66][128 * ii0 - 2 * i0 + 133][1]);
                        }
                        for (int i4 = 1; i4 < 128 * ii0 - 2 * i0 + i3 + i2 / 3 - 18; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = max(max(1, 128 * ii0 - 2 * i0 - i2 - i3 + 201), 128 * ii0 - 2 * i0 + i3 + i2 / 3 - 18); i4 <= 128 * ii0 - 2 * i0 - i3 + 135; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = max(max(1, 128 * ii0 - 2 * i0 - i3 + 136), 128 * ii0 - 2 * i0 + i3 + i2 / 3 - 18); i4 <= min(_PB_N - 2, -2 * i3 + 192); i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        if (i3 <= 63) {
                          for (int i4 = -2 * i3 + 193; i4 < _PB_N - 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                        for (int i4 = -2 * i3 + 193; i4 <= 63; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        if (i3 >= 64) {
                          for (int i4 = max(64, -2 * i3 + 193); i4 < _PB_N - 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                    }
                    if (i0 >= 64 * ii0 + 3 && i2 >= 64 * ii2 + 2) {
                      for (int i3 = 57 * ii2 + (-ii2 + i2 - 2) / 9 + 2; i3 < min(min(min(_PB_N - 1, -128 * ii0 + 2 * i0 + i2 - 2), _PB_N - 128 * ii0 + 21 * ii2 + 2 * i0 - (-ii2 + i2 + 1) / 3 - 4), -43 * ii2 + 2 * i2 - (-ii2 + i2 + 1) / 3); i3 += 1) {
                        for (int i4 = max(1, 128 * ii0 + 64 * ii2 - 2 * i0 + 4); i4 <= min(64 * ii2 - 1, 128 * ii0 - 21 * ii2 - 2 * i0 + i3 + (-ii2 + i2 + 1) / 3 + 2); i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = max(1, 128 * ii0 - 21 * ii2 - 2 * i0 + i3 + (-ii2 + i2 + 1) / 3 + 3); i4 < 64 * ii2; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = 64 * ii2; i4 < _PB_N - 1; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    } else if (64 * ii2 + 1 >= i2) {
                      for (int i3 = max(1, 128 * ii0 + 171 * ii2 - 2 * i0 - 2 * i2 + (-ii2 + i2 + 1) / 3 + 9); i3 < -43 * ii2 + 2 * i2 - (-ii2 + i2 + 1) / 3; i3 += 1) {
                        for (int i4 = max(1, 128 * ii0 + 192 * ii2 - 2 * i0 - i2 - i3 + 6); i4 < _PB_N - 1; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    }
                    if (ii2 == 1) {
                      for (int i3 = max(2 * i2 - (i2 + 3) / 3 - 42, 128 * ii0 - 2 * i0 - 2 * i2 + i2 / 3 + 180); i3 < min(_PB_N - 1, -128 * ii0 + 2 * i0 + 2 * i2 - (i2 + 3) / 3 - 109); i3 += 1) {
                        for (int i4 = max(1, 128 * ii0 - 2 * i0 - i2 - i3 + 198); i4 < _PB_N - 1; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    }
                    for (int i3 = max(max(128 * ii0 + 171 * ii2 - 2 * i0 - 2 * i2 + (-ii2 + i2 + 1) / 3 + 9, -128 * ii0 - 107 * ii2 + 2 * i0 + 2 * i2 - (-ii2 + i2 + 1) / 3 - 3), -43 * ii2 + 2 * i2 - (-ii2 + i2 + 1) / 3); i3 <= 43 * ii2 + (-ii2 + i2 + 2) / 3; i3 += 1) {
                      for (int i4 = 128 * ii0 + 192 * ii2 - 2 * i0 - i2 - i3 + 6; i4 <= 128 * ii0 + 192 * ii2 - 2 * i0 - i2 - i3 + 8; i4 += 1) {
                        B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                      }
                      for (int i4 = 128 * ii0 + 192 * ii2 - 2 * i0 - i2 - i3 + 9; i4 < _PB_N - 1; i4 += 1) {
                        B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                      }
                    }
                    if (i0 >= 64 * ii0 + 3) {
                      for (int i3 = max(max(max(128 * ii0 + 128 * ii2 - 2 * i0 - i2 + 5, 43 * ii2 + (-ii2 + i2 + 2) / 3 + 1), -128 * ii0 - 107 * ii2 + 2 * i0 + 2 * i2 - (-ii2 + i2 + 1) / 3 - 3), -43 * ii2 + 2 * i2 - (-ii2 + i2 + 1) / 3); i3 < min(min(_PB_N - 1, -128 * ii0 + 2 * i0 + i2 - 2), _PB_N - 128 * ii0 - 107 * ii2 + 2 * i0 + 2 * i2 - (-ii2 + i2 + 1) / 3 - 6); i3 += 1) {
                        for (int i4 = max(max(max(1, 128 * ii0 + 64 * ii2 - 2 * i0 + 4), 128 * ii0 + 128 * ii2 - 2 * i0 - i2 + 5), 128 * ii0 + 192 * ii2 - 2 * i0 - i2 - i3 + 6); i4 <= 128 * ii0 + 107 * ii2 - 2 * i0 - 2 * i2 + i3 + (-ii2 + i2 + 1) / 3 + 4; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = 128 * ii0 + 107 * ii2 - 2 * i0 - 2 * i2 + i3 + (-ii2 + i2 + 1) / 3 + 5; i4 <= min(64 * ii2 - 1, 128 * ii0 - 21 * ii2 - 2 * i0 + i3 + (-ii2 + i2 + 1) / 3 + 2); i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = max(64 * ii2, 128 * ii0 + 107 * ii2 - 2 * i0 - 2 * i2 + i3 + (-ii2 + i2 + 1) / 3 + 5); i4 <= 128 * ii0 - 21 * ii2 - 2 * i0 + i3 + (-ii2 + i2 + 1) / 3 + 2; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = 128 * ii0 + 107 * ii2 - 2 * i0 - 2 * i2 + i3 + (-ii2 + i2 + 1) / 3 + 5; i4 <= min(128 * ii0 + 192 * ii2 - 2 * i0 - i2 - i3 + 8, 21 * ii2 + i2 - (-ii2 + i2 + 1) / 3); i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = max(max(128 * ii0 + 192 * ii2 - 2 * i0 - i2 - i3 + 9, 128 * ii0 - 21 * ii2 - 2 * i0 + i3 + (-ii2 + i2 + 1) / 3 + 3), 128 * ii0 + 107 * ii2 - 2 * i0 - 2 * i2 + i3 + (-ii2 + i2 + 1) / 3 + 5); i4 < min(_PB_N - 1, -64 * ii0 - 64 * ii2 + i0 + i2 + (i2 + i3) / 2 - 3); i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        if (64 * ii2 >= i2 + 1) {
                          for (int i4 = max(max(128 * ii0 + 107 * ii2 - 2 * i0 - 2 * i2 + i3 + (-ii2 + i2 + 1) / 3 + 5, 21 * ii2 + i2 - (-ii2 + i2 + 1) / 3 + 1), -64 * ii0 - 64 * ii2 + i0 + i2 + (i2 + i3) / 2 - 3); i4 < _PB_N - 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        } else {
                          for (int i4 = max(max(128 * ii2 - i2 + 1, 128 * ii0 + 107 * ii2 - 2 * i0 - 2 * i2 + i3 + (-ii2 + i2 + 1) / 3 + 5), -64 * ii0 - 64 * ii2 + i0 + i2 + (i2 + i3) / 2 - 3); i4 < _PB_N - 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                      for (int i3 = -128 * ii0 + 2 * i0 + i2 - 2; i3 < min(min(_PB_N - 1, _PB_N - 128 * ii0 - 107 * ii2 + 2 * i0 + 2 * i2 - (-ii2 + i2 + 1) / 3 - 6), _PB_N - 128 * ii0 + 21 * ii2 + 2 * i0 - (-ii2 + i2 + 1) / 3 - 4); i3 += 1) {
                        for (int i4 = max(128 * ii0 + 64 * ii2 - 2 * i0 + 4, 128 * ii0 + 128 * ii2 - 2 * i0 - i2 + 5); i4 < _PB_N - 1; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                      for (int i3 = _PB_N - 128 * ii0 + 21 * ii2 + 2 * i0 - (-ii2 + i2 + 1) / 3 - 4; i3 < min(_PB_N - 1, -43 * ii2 + 2 * i2 - (-ii2 + i2 + 1) / 3); i3 += 1) {
                        for (int i4 = 128 * ii0 + 64 * ii2 - 2 * i0 + 4; i4 < _PB_N - 1; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                      for (int i3 = max(128 * ii0 + 128 * ii2 - 2 * i0 - i2 + 5, _PB_N - 128 * ii0 - 107 * ii2 + 2 * i0 + 2 * i2 - (-ii2 + i2 + 1) / 3 - 6); i3 < _PB_N - 1; i3 += 1) {
                        for (int i4 = max(128 * ii0 + 128 * ii2 - 2 * i0 - i2 + 5, 128 * ii0 + 192 * ii2 - 2 * i0 - i2 - i3 + 6); i4 < _PB_N - 1; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                      for (int i3 = max(_PB_N - 128 * ii0 + 21 * ii2 + 2 * i0 - (-ii2 + i2 + 1) / 3 - 4, -43 * ii2 + 2 * i2 - (-ii2 + i2 + 1) / 3); i3 < _PB_N - 1; i3 += 1) {
                        for (int i4 = 128 * ii0 + 64 * ii2 - 2 * i0 + 4; i4 < _PB_N - 1; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    }
                    if (i0 == 64 * ii0 + 2) {
                      for (int i3 = max(64 * ii2, 128 * ii2 - i2 + 1); i3 <= min(min(_PB_N - 2, i2 + 1), _PB_N + 21 * ii2 - (-ii2 + i2 + 1) / 3 - 1); i3 += 1) {
                        for (int i4 = max(max(64 * ii2, 128 * ii2 - i2 + 1), 192 * ii2 - 2 * i3 + 1); i4 < _PB_N - 1; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                      for (int i3 = i2 + 2; i3 < min(_PB_N - 1, _PB_N + 21 * ii2 - (-ii2 + i2 + 1) / 3); i3 += 1) {
                        for (int i4 = max(64 * ii2, 128 * ii2 - i2 + 1); i4 < _PB_N - 1; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                      for (int i3 = _PB_N + 21 * ii2 - (-ii2 + i2 + 1) / 3; i3 < _PB_N - 1; i3 += 1) {
                        for (int i4 = 64 * ii2; i4 < _PB_N - 1; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    }
                  }
                }
                for (int i2 = max(1, 128 * ii0 + 64 * ii2 - 2 * i0 + 3); i2 <= 64 * ii2; i2 += 1) {
                  for (int i3 = max(1, 128 * ii0 + 128 * ii2 - 2 * i0 - i2 + 4); i3 < _PB_N - 1; i3 += 1) {
                    for (int i4 = max(max(1, 128 * ii0 + 128 * ii2 - 2 * i0 - i2 + 4), 128 * ii0 + 192 * ii2 - 2 * i0 - i2 - i3 + 5); i4 < _PB_N - 1; i4 += 1) {
                      A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                    }
                  }
                }
                for (int i2 = 64 * ii2 + 1; i2 < _PB_N - 1; i2 += 1) {
                  for (int i3 = max(1, 128 * ii0 + 64 * ii2 - 2 * i0 + 3); i3 <= 64 * ii2; i3 += 1) {
                    for (int i4 = max(1, 128 * ii0 + 128 * ii2 - 2 * i0 - i3 + 4); i4 < _PB_N - 1; i4 += 1) {
                      A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                    }
                  }
                  for (int i3 = 64 * ii2 + 1; i3 < _PB_N - 1; i3 += 1) {
                    for (int i4 = max(1, 128 * ii0 + 64 * ii2 - 2 * i0 + 3); i4 < _PB_N - 1; i4 += 1) {
                      A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                    }
                  }
                }
              }
            }
          }
        }
      }
      if (_PB_N <= 68) {
        for (int ii3 = 0; ii3 <= 1; ii3 += 1) {
          for (int ii4 = 0; ii4 <= 1; ii4 += 1) {
            if (_PB_N == 67) {
              if (ii3 == 0 && ii4 == 0) {
                for (int i0 = 64 * ii0 + 1; i0 <= 64 * ii0 + 3; i0 += 1) {
                  if (i0 >= 64 * ii0 + 2) {
                    if (i0 == 64 * ii0 + 3) {
                      for (int i2 = 62; i2 <= 65; i2 += 1) {
                        if (i2 >= 64) {
                          for (int i3 = 1; i3 <= -i2 + 126; i3 += 1) {
                            if (2 * i2 + i3 >= 131) {
                              for (int i4 = 1; i4 <= -i2 + 126; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            } else {
                              if (i3 == 2) {
                                B[64][2][1] = ((((SCALAR_VAL(0.125) * ((A[65][2][1] - (SCALAR_VAL(2.0) * A[64][2][1])) + A[63][2][1])) + (SCALAR_VAL(0.125) * ((A[64][3][1] - (SCALAR_VAL(2.0) * A[64][2][1])) + A[64][1][1]))) + (SCALAR_VAL(0.125) * ((A[64][2][2] - (SCALAR_VAL(2.0) * A[64][2][1])) + A[64][2][0]))) + A[64][2][1]);
                              } else {
                                B[64][1][1] = ((((SCALAR_VAL(0.125) * ((A[65][1][1] - (SCALAR_VAL(2.0) * A[64][1][1])) + A[63][1][1])) + (SCALAR_VAL(0.125) * ((A[64][2][1] - (SCALAR_VAL(2.0) * A[64][1][1])) + A[64][0][1]))) + (SCALAR_VAL(0.125) * ((A[64][1][2] - (SCALAR_VAL(2.0) * A[64][1][1])) + A[64][1][0]))) + A[64][1][1]);
                              }
                              for (int i4 = 2; i4 <= 62; i4 += 1) {
                                B[64][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[65][i3][i4] - (SCALAR_VAL(2.0) * A[64][i3][i4])) + A[63][i3][i4])) + (SCALAR_VAL(0.125) * ((A[64][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[64][i3][i4])) + A[64][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[64][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[64][i3][i4])) + A[64][i3][i4 - 1]))) + A[64][i3][i4]);
                              }
                            }
                          }
                        } else {
                          for (int i3 = 1; i3 <= 2; i3 += 1) {
                            if (i3 == 1) {
                              B[i2][1][1] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][1][1] - (SCALAR_VAL(2.0) * A[i2][1][1])) + A[i2 - 1][1][1])) + (SCALAR_VAL(0.125) * ((A[i2][2][1] - (SCALAR_VAL(2.0) * A[i2][1][1])) + A[i2][0][1]))) + (SCALAR_VAL(0.125) * ((A[i2][1][2] - (SCALAR_VAL(2.0) * A[i2][1][1])) + A[i2][1][0]))) + A[i2][1][1]);
                            }
                            for (int i4 = -i3 + 3; i4 <= -i2 + 126; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                          for (int i3 = 3; i3 <= -i2 + 126; i3 += 1) {
                            for (int i4 = 1; i4 <= -2 * i2 + 188; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            if (i2 == 63) {
                              B[63][i3][63] = ((((SCALAR_VAL(0.125) * ((A[64][i3][63] - (SCALAR_VAL(2.0) * A[63][i3][63])) + A[62][i3][63])) + (SCALAR_VAL(0.125) * ((A[63][i3 + 1][63] - (SCALAR_VAL(2.0) * A[63][i3][63])) + A[63][i3 - 1][63]))) + (SCALAR_VAL(0.125) * ((A[63][i3][64] - (SCALAR_VAL(2.0) * A[63][i3][63])) + A[63][i3][62]))) + A[63][i3][63]);
                            }
                          }
                        }
                      }
                    } else {
                      for (int i2 = 64; i2 <= 65; i2 += 1) {
                        for (int i3 = 1; i3 <= -i2 + 128; i3 += 1) {
                          for (int i4 = 1; i4 <= -i2 + 128; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                    }
                  }
                  for (int i2 = 128 * ii0 - 2 * i0 + 67; i2 <= 65; i2 += 1) {
                    for (int i3 = 1; i3 <= 128 * ii0 - 2 * i0 - i2 + 131; i3 += 1) {
                      for (int i4 = 1; i4 <= 128 * ii0 - 2 * i0 - i2 + 131; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                  }
                }
              }
              for (int i0 = max(64 * ii0 + 1, 64 * ii0 - 3 * ii3 - 3 * ii4 + 4); i0 <= min(min(_PB_TSTEPS, 64 * ii0 + 64), 64 * ii0 + 32 * ii3 + 31 * ii4 + 33); i0 += 1) {
                if (2 * ii4 + i0 >= 64 * ii0 + 4 && ii3 + i0 >= 64 * ii0 + 3 && 64 * ii0 + 32 * ii4 + 33 >= i0 && 64 * ii0 + 31 * ii3 + 33 >= i0) {
                  if (ii3 + ii4 <= 1) {
                    if (ii3 == 0 && ii4 == 1) {
                      for (int i2 = 128 * ii0 - 2 * i0 + 68; i2 <= 128 * ii0 - 2 * i0 + 69; i2 += 1) {
                        for (int i3 = 1; i3 <= 128 * ii0 - 2 * i0 - i2 + 132; i3 += 1) {
                          if (i3 <= 63) {
                            for (int i4 = 128 * ii0 - 2 * i0 - i2 + 133; i4 <= 65; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          } else {
                            B[128 * ii0 - 2 * i0 + 68][64][65] = ((((SCALAR_VAL(0.125) * ((A[128 * ii0 - 2 * i0 + 69][64][65] - (SCALAR_VAL(2.0) * A[128 * ii0 - 2 * i0 + 68][64][65])) + A[128 * ii0 - 2 * i0 + 67][64][65])) + (SCALAR_VAL(0.125) * ((A[128 * ii0 - 2 * i0 + 68][65][65] - (SCALAR_VAL(2.0) * A[128 * ii0 - 2 * i0 + 68][64][65])) + A[128 * ii0 - 2 * i0 + 68][63][65]))) + (SCALAR_VAL(0.125) * ((A[128 * ii0 - 2 * i0 + 68][64][66] - (SCALAR_VAL(2.0) * A[128 * ii0 - 2 * i0 + 68][64][65])) + A[128 * ii0 - 2 * i0 + 68][64][64]))) + A[128 * ii0 - 2 * i0 + 68][64][65]);
                          }
                        }
                      }
                    } else if (ii3 == 0) {
                      for (int i2 = 128 * ii0 - 2 * i0 + 68; i2 <= 128 * ii0 - 2 * i0 + 69; i2 += 1) {
                        for (int i3 = 1; i3 <= 128 * ii0 - 2 * i0 - i2 + 132; i3 += 1) {
                          if (i3 >= 3) {
                            for (int i4 = 1; i4 <= 256 * ii0 - 4 * i0 - 2 * i2 + 200; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            if (2 * i0 + i2 == 128 * ii0 + 69) {
                              B[128 * ii0 - 2 * i0 + 69][i3][63] = ((((SCALAR_VAL(0.125) * ((A[128 * ii0 - 2 * i0 + 70][i3][63] - (SCALAR_VAL(2.0) * A[128 * ii0 - 2 * i0 + 69][i3][63])) + A[128 * ii0 - 2 * i0 + 68][i3][63])) + (SCALAR_VAL(0.125) * ((A[128 * ii0 - 2 * i0 + 69][i3 + 1][63] - (SCALAR_VAL(2.0) * A[128 * ii0 - 2 * i0 + 69][i3][63])) + A[128 * ii0 - 2 * i0 + 69][i3 - 1][63]))) + (SCALAR_VAL(0.125) * ((A[128 * ii0 - 2 * i0 + 69][i3][64] - (SCALAR_VAL(2.0) * A[128 * ii0 - 2 * i0 + 69][i3][63])) + A[128 * ii0 - 2 * i0 + 69][i3][62]))) + A[128 * ii0 - 2 * i0 + 69][i3][63]);
                            }
                          } else {
                            if (i3 == 1) {
                              B[i2][1][1] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][1][1] - (SCALAR_VAL(2.0) * A[i2][1][1])) + A[i2 - 1][1][1])) + (SCALAR_VAL(0.125) * ((A[i2][2][1] - (SCALAR_VAL(2.0) * A[i2][1][1])) + A[i2][0][1]))) + (SCALAR_VAL(0.125) * ((A[i2][1][2] - (SCALAR_VAL(2.0) * A[i2][1][1])) + A[i2][1][0]))) + A[i2][1][1]);
                            }
                            for (int i4 = -i3 + 3; i4 <= 128 * ii0 - 2 * i0 - i2 + 132; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                        }
                      }
                      for (int i3 = 1; i3 <= 62; i3 += 1) {
                        if (i3 == 1) {
                          B[128 * ii0 - 2 * i0 + 70][1][1] = ((((SCALAR_VAL(0.125) * ((A[128 * ii0 - 2 * i0 + 71][1][1] - (SCALAR_VAL(2.0) * A[128 * ii0 - 2 * i0 + 70][1][1])) + A[128 * ii0 - 2 * i0 + 69][1][1])) + (SCALAR_VAL(0.125) * ((A[128 * ii0 - 2 * i0 + 70][2][1] - (SCALAR_VAL(2.0) * A[128 * ii0 - 2 * i0 + 70][1][1])) + A[128 * ii0 - 2 * i0 + 70][0][1]))) + (SCALAR_VAL(0.125) * ((A[128 * ii0 - 2 * i0 + 70][1][2] - (SCALAR_VAL(2.0) * A[128 * ii0 - 2 * i0 + 70][1][1])) + A[128 * ii0 - 2 * i0 + 70][1][0]))) + A[128 * ii0 - 2 * i0 + 70][1][1]);
                        }
                        for (int i4 = max(1, -i3 + 3); i4 <= 62; i4 += 1) {
                          B[128 * ii0 - 2 * i0 + 70][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[128 * ii0 - 2 * i0 + 71][i3][i4] - (SCALAR_VAL(2.0) * A[128 * ii0 - 2 * i0 + 70][i3][i4])) + A[128 * ii0 - 2 * i0 + 69][i3][i4])) + (SCALAR_VAL(0.125) * ((A[128 * ii0 - 2 * i0 + 70][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[128 * ii0 - 2 * i0 + 70][i3][i4])) + A[128 * ii0 - 2 * i0 + 70][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[128 * ii0 - 2 * i0 + 70][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[128 * ii0 - 2 * i0 + 70][i3][i4])) + A[128 * ii0 - 2 * i0 + 70][i3][i4 - 1]))) + A[128 * ii0 - 2 * i0 + 70][i3][i4]);
                        }
                      }
                    }
                    for (int i2 = max(128 * ii0 - 2 * ii3 - 2 * i0 + 70, 128 * ii0 - 3 * ii3 - 3 * ii4 - 2 * i0 + 71); i2 <= 63; i2 += 1) {
                      if (ii3 == 0) {
                        for (int i3 = 1; i3 <= ii4 + 1; i3 += 1) {
                          if (ii4 == 0 && i3 == 1) {
                            B[i2][1][1] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][1][1] - (SCALAR_VAL(2.0) * A[i2][1][1])) + A[i2 - 1][1][1])) + (SCALAR_VAL(0.125) * ((A[i2][2][1] - (SCALAR_VAL(2.0) * A[i2][1][1])) + A[i2][0][1]))) + (SCALAR_VAL(0.125) * ((A[i2][1][2] - (SCALAR_VAL(2.0) * A[i2][1][1])) + A[i2][1][0]))) + A[i2][1][1]);
                          }
                          for (int i4 = max(2, 128 * ii0 + 60 * ii4 - 2 * i0 - i2 + 73); i4 <= min(65, 128 * ii0 + 62 * ii4 - 2 * i0 - i2 + 132); i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                      if (i0 >= 64 * ii0 + 4) {
                        for (int i3 = max(ii4 + 2, 128 * ii0 + 60 * ii3 - 2 * i0 - i2 + 73); i3 <= min(64, 128 * ii0 + 61 * ii3 - 2 * i0 - i2 + 132); i3 += 1) {
                          if (ii4 == 0) {
                            for (int i4 = 1; i4 <= min(128 * ii0 + 61 * ii3 - 2 * i0 - i2 + 132, 128 * ii0 - 2 * i0 - i2 - i3 + 197); i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          } else {
                            for (int i4 = 128 * ii0 - 2 * i0 - i2 + 133; i4 <= 65; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                        }
                        if (ii3 == 1 && ii4 == 0) {
                          for (int i4 = 1; i4 <= 128 * ii0 - 2 * i0 - i2 + 132; i4 += 1) {
                            B[i2][65][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][65][i4] - (SCALAR_VAL(2.0) * A[i2][65][i4])) + A[i2 - 1][65][i4])) + (SCALAR_VAL(0.125) * ((A[i2][66][i4] - (SCALAR_VAL(2.0) * A[i2][65][i4])) + A[i2][64][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][65][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][65][i4])) + A[i2][65][i4 - 1]))) + A[i2][65][i4]);
                          }
                        }
                      }
                    }
                    for (int i2 = 64; i2 <= min(65, 43 * ii0 + ii3 + ii4 - i0 + (-ii0 - ii3 - ii4 + i0 + 1) / 3 + 86); i2 += 1) {
                      for (int i3 = max(max(1, 128 * ii0 + 59 * ii3 - 2 * ii4 - 2 * i0 + 9), 128 * ii0 + 60 * ii3 - 2 * ii4 - 2 * i0 - 61 * i2 + 3913); i3 <= min(128 * ii0 - 2 * i0 + 68, 128 * ii0 + 2 * ii3 - 2 * i0 - i2 + 132); i3 += 1) {
                        if (ii4 == 0 || 1 == 0 || (ii3 == 0 && i3 >= 3)) {
                          for (int i4 = max(1, 128 * ii0 + 60 * ii4 - 2 * i0 - i2 + 73); i4 <= min(min(-59 * ii3 + 5 * ii4 + 60, 128 * ii0 + 62 * ii4 - 2 * i0 - i2 + 132), 128 * ii0 + 61 * ii4 - 2 * i0 + i3 + 66); i4 += 1) {
                            if ((ii4 == 0 && 57 * ii3 + 64 >= 2 * i3 + i4 && 128 * ii0 + 69 >= 2 * i0 + i3 + i4) || ii4 == 1 || i3 >= 2 || 1 == 0) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                          if (ii3 == 1 && ii4 == 0 && i0 >= 64 * ii0 + 4 && i2 == 65 && 2 * i0 + i3 == 128 * ii0 + 68) {
                            for (int i4 = 2; i4 <= 64; i4 += 1) {
                              B[65][128 * ii0 - 2 * i0 + 68][i4] = ((((SCALAR_VAL(0.125) * ((A[66][128 * ii0 - 2 * i0 + 68][i4] - (SCALAR_VAL(2.0) * A[65][128 * ii0 - 2 * i0 + 68][i4])) + A[64][128 * ii0 - 2 * i0 + 68][i4])) + (SCALAR_VAL(0.125) * ((A[65][128 * ii0 - 2 * i0 + 69][i4] - (SCALAR_VAL(2.0) * A[65][128 * ii0 - 2 * i0 + 68][i4])) + A[65][128 * ii0 - 2 * i0 + 67][i4]))) + (SCALAR_VAL(0.125) * ((A[65][128 * ii0 - 2 * i0 + 68][i4 + 1] - (SCALAR_VAL(2.0) * A[65][128 * ii0 - 2 * i0 + 68][i4])) + A[65][128 * ii0 - 2 * i0 + 68][i4 - 1]))) + A[65][128 * ii0 - 2 * i0 + 68][i4]);
                            }
                          } else if (ii3 == 0 && ii4 == 0 && i2 == 64 && i3 == 1) {
                            B[64][1][128 * ii0 - 2 * i0 + 68] = ((((SCALAR_VAL(0.125) * ((A[65][1][128 * ii0 - 2 * i0 + 68] - (SCALAR_VAL(2.0) * A[64][1][128 * ii0 - 2 * i0 + 68])) + A[63][1][128 * ii0 - 2 * i0 + 68])) + (SCALAR_VAL(0.125) * ((A[64][2][128 * ii0 - 2 * i0 + 68] - (SCALAR_VAL(2.0) * A[64][1][128 * ii0 - 2 * i0 + 68])) + A[64][0][128 * ii0 - 2 * i0 + 68]))) + (SCALAR_VAL(0.125) * ((A[64][1][128 * ii0 - 2 * i0 + 69] - (SCALAR_VAL(2.0) * A[64][1][128 * ii0 - 2 * i0 + 68])) + A[64][1][128 * ii0 - 2 * i0 + 67]))) + A[64][1][128 * ii0 - 2 * i0 + 68]);
                          }
                        }
                        if (ii3 == 0 && ii4 == 1 && i3 <= 2) {
                          if (i2 == 65 && i3 == 1) {
                            B[65][1][128 * ii0 - 2 * i0 + 68] = ((((SCALAR_VAL(0.125) * ((A[66][1][128 * ii0 - 2 * i0 + 68] - (SCALAR_VAL(2.0) * A[65][1][128 * ii0 - 2 * i0 + 68])) + A[64][1][128 * ii0 - 2 * i0 + 68])) + (SCALAR_VAL(0.125) * ((A[65][2][128 * ii0 - 2 * i0 + 68] - (SCALAR_VAL(2.0) * A[65][1][128 * ii0 - 2 * i0 + 68])) + A[65][0][128 * ii0 - 2 * i0 + 68]))) + (SCALAR_VAL(0.125) * ((A[65][1][128 * ii0 - 2 * i0 + 69] - (SCALAR_VAL(2.0) * A[65][1][128 * ii0 - 2 * i0 + 68])) + A[65][1][128 * ii0 - 2 * i0 + 67]))) + A[65][1][128 * ii0 - 2 * i0 + 68]);
                          }
                          if (128 * ii0 + 259 >= 2 * i0 + 3 * i2) {
                            for (int i4 = max(128 * ii0 - 2 * i0 - i2 + 133, 128 * ii0 - 2 * i0 - i3 + 70); i4 <= min(65, 128 * ii0 - 2 * i0 - i2 + 194); i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            if (i0 == 64 * ii0 + 33 && i2 == 64) {
                              B[64][i3][65] = ((((SCALAR_VAL(0.125) * ((A[65][i3][65] - (SCALAR_VAL(2.0) * A[64][i3][65])) + A[63][i3][65])) + (SCALAR_VAL(0.125) * ((A[64][i3 + 1][65] - (SCALAR_VAL(2.0) * A[64][i3][65])) + A[64][i3 - 1][65]))) + (SCALAR_VAL(0.125) * ((A[64][i3][66] - (SCALAR_VAL(2.0) * A[64][i3][65])) + A[64][i3][64]))) + A[64][i3][65]);
                            }
                          } else {
                            for (int i4 = 3; i4 <= 65; i4 += 1) {
                              B[65][1][i4] = ((((SCALAR_VAL(0.125) * ((A[66][1][i4] - (SCALAR_VAL(2.0) * A[65][1][i4])) + A[64][1][i4])) + (SCALAR_VAL(0.125) * ((A[65][2][i4] - (SCALAR_VAL(2.0) * A[65][1][i4])) + A[65][0][i4]))) + (SCALAR_VAL(0.125) * ((A[65][1][i4 + 1] - (SCALAR_VAL(2.0) * A[65][1][i4])) + A[65][1][i4 - 1]))) + A[65][1][i4]);
                            }
                          }
                        }
                      }
                      if (ii3 == 1 && ii4 == 0) {
                        for (int i3 = 128 * ii0 - 2 * i0 + 69; i3 <= 65; i3 += 1) {
                          if (i3 <= 64) {
                            if (2 * i0 + i3 == 128 * ii0 + 69) {
                              B[i2][128 * ii0 - 2 * i0 + 69][1] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][128 * ii0 - 2 * i0 + 69][1] - (SCALAR_VAL(2.0) * A[i2][128 * ii0 - 2 * i0 + 69][1])) + A[i2 - 1][128 * ii0 - 2 * i0 + 69][1])) + (SCALAR_VAL(0.125) * ((A[i2][128 * ii0 - 2 * i0 + 70][1] - (SCALAR_VAL(2.0) * A[i2][128 * ii0 - 2 * i0 + 69][1])) + A[i2][128 * ii0 - 2 * i0 + 68][1]))) + (SCALAR_VAL(0.125) * ((A[i2][128 * ii0 - 2 * i0 + 69][2] - (SCALAR_VAL(2.0) * A[i2][128 * ii0 - 2 * i0 + 69][1])) + A[i2][128 * ii0 - 2 * i0 + 69][0]))) + A[i2][128 * ii0 - 2 * i0 + 69][1]);
                            }
                            for (int i4 = max(1, 128 * ii0 - 2 * i0 - i3 + 71); i4 <= 128 * ii0 - 2 * i0 - i2 - i3 + 197; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          } else {
                            for (int i4 = 1; i4 <= 128 * ii0 - 2 * i0 - i2 + 132; i4 += 1) {
                              B[i2][65][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][65][i4] - (SCALAR_VAL(2.0) * A[i2][65][i4])) + A[i2 - 1][65][i4])) + (SCALAR_VAL(0.125) * ((A[i2][66][i4] - (SCALAR_VAL(2.0) * A[i2][65][i4])) + A[i2][64][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][65][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][65][i4])) + A[i2][65][i4 - 1]))) + A[i2][65][i4]);
                            }
                          }
                        }
                      }
                    }
                    if (ii3 == 0 && ii4 == 0 && i0 == 64 * ii0 + 33) {
                      B[65][1][1] = ((((SCALAR_VAL(0.125) * ((A[66][1][1] - (SCALAR_VAL(2.0) * A[65][1][1])) + A[64][1][1])) + (SCALAR_VAL(0.125) * ((A[65][2][1] - (SCALAR_VAL(2.0) * A[65][1][1])) + A[65][0][1]))) + (SCALAR_VAL(0.125) * ((A[65][1][2] - (SCALAR_VAL(2.0) * A[65][1][1])) + A[65][1][0]))) + A[65][1][1]);
                    }
                  } else if (i0 >= 64 * ii0 + 3) {
                    for (int i2 = max(1, 128 * ii0 - 2 * i0 + 68); i2 <= min(64, 77 * ii0 - i0 - (ii0 + i0 + 4) / 5 + 105); i2 += 1) {
                      for (int i3 = max(1, 128 * ii0 - 2 * i0 - i2 + 133); i3 <= (i2 + 1) / 3 + 43; i3 += 1) {
                        for (int i4 = max(1, 128 * ii0 - 2 * i0 - i2 - i3 + 198); i4 <= min(128 * ii0 - 2 * i0 - i2 - i3 + 200, i2 - (i2 + 3) / 3 + 22); i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = max(1, 128 * ii0 - 2 * i0 - i2 - i3 + 201); i4 <= min(65, -64 * ii0 + i0 + i2 + (i2 + i3) / 2 - 68); i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = max(max(128 * ii0 - 2 * i0 - i2 - i3 + 198, i2 - (i2 + 3) / 3 + 23), -64 * ii0 + i0 + i2 + (i2 + i3) / 2 - 67); i4 <= 65; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                      for (int i3 = (i2 + 1) / 3 + 44; i3 < -128 * ii0 + 2 * i0 + 2 * i2 - (i2 + 3) / 3 - 109; i3 += 1) {
                        if (2 * i0 + i2 + i3 == 128 * ii0 + 199) {
                          B[i2][128 * ii0 - 2 * i0 - i2 + 199][1] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][128 * ii0 - 2 * i0 - i2 + 199][1] - (SCALAR_VAL(2.0) * A[i2][128 * ii0 - 2 * i0 - i2 + 199][1])) + A[i2 - 1][128 * ii0 - 2 * i0 - i2 + 199][1])) + (SCALAR_VAL(0.125) * ((A[i2][128 * ii0 - 2 * i0 - i2 + 200][1] - (SCALAR_VAL(2.0) * A[i2][128 * ii0 - 2 * i0 - i2 + 199][1])) + A[i2][128 * ii0 - 2 * i0 - i2 + 198][1]))) + (SCALAR_VAL(0.125) * ((A[i2][128 * ii0 - 2 * i0 - i2 + 199][2] - (SCALAR_VAL(2.0) * A[i2][128 * ii0 - 2 * i0 - i2 + 199][1])) + A[i2][128 * ii0 - 2 * i0 - i2 + 199][0]))) + A[i2][128 * ii0 - 2 * i0 - i2 + 199][1]);
                        }
                        for (int i4 = max(1, 128 * ii0 - 2 * i0 - i2 - i3 + 201); i4 <= 65; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                      for (int i3 = max(max(128 * ii0 - 2 * i0 - i2 + 133, -128 * ii0 + 2 * i0 + 2 * i2 - (i2 + 3) / 3 - 109), (i2 + 1) / 3 + 44); i3 <= 65; i3 += 1) {
                        for (int i4 = max(1, 128 * ii0 - 2 * i0 - i2 - i3 + 198); i4 <= min(65, 128 * ii0 - 2 * i0 - 2 * i2 + i3 + i2 / 3 + 111); i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = 128 * ii0 - 2 * i0 - 2 * i2 + i3 + i2 / 3 + 112; i4 <= min(128 * ii0 - 2 * i0 - i2 - i3 + 200, i2 - (i2 + 3) / 3 + 22); i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = max(128 * ii0 - 2 * i0 - i2 - i3 + 201, 128 * ii0 - 2 * i0 - 2 * i2 + i3 + i2 / 3 + 112); i4 <= min(65, -64 * ii0 + i0 + i2 + (i2 + i3) / 2 - 68); i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        if (i2 <= 63) {
                          for (int i4 = max(max(i2 - (i2 + 3) / 3 + 23, 128 * ii0 - 2 * i0 - 2 * i2 + i3 + i2 / 3 + 112), -64 * ii0 + i0 + i2 + (i2 + i3) / 2 - 67); i4 <= 65; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        } else if (64 * ii0 + 4 >= i0) {
                          B[64][65][65] = ((((SCALAR_VAL(0.125) * ((A[65][65][65] - (SCALAR_VAL(2.0) * A[64][65][65])) + A[63][65][65])) + (SCALAR_VAL(0.125) * ((A[64][66][65] - (SCALAR_VAL(2.0) * A[64][65][65])) + A[64][64][65]))) + (SCALAR_VAL(0.125) * ((A[64][65][66] - (SCALAR_VAL(2.0) * A[64][65][65])) + A[64][65][64]))) + A[64][65][65]);
                        }
                      }
                    }
                    for (int i2 = 77 * ii0 - i0 - (ii0 + i0 + 4) / 5 + 106; i2 <= min(64, 43 * ii0 - i0 + (-ii0 + i0) / 3 + 88); i2 += 1) {
                      for (int i3 = 1; i3 <= 65; i3 += 1) {
                        for (int i4 = max(1, 128 * ii0 - 2 * i0 - i2 - i3 + 198); i4 <= min(128 * ii0 - 2 * i0 - i2 - i3 + 200, i2 - (i2 + 3) / 3 + 22); i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = max(1, 128 * ii0 - 2 * i0 - i2 - i3 + 201); i4 <= min(65, -64 * ii0 + i0 + i2 + (i2 + i3) / 2 - 68); i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = max(i2 - (i2 + 3) / 3 + 23, -64 * ii0 + i0 + i2 + (i2 + i3) / 2 - 67); i4 <= 65; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    }
                    for (int i2 = 43 * ii0 - i0 + (-ii0 + i0) / 3 + 89; i2 <= 63; i2 += 1) {
                      for (int i3 = 1; i3 <= 65; i3 += 1) {
                        for (int i4 = max(1, 128 * ii0 - 2 * i0 - i2 - i3 + 198); i4 <= 65; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    }
                    if (i0 >= 64 * ii0 + 37) {
                      for (int i3 = 1; i3 <= 65; i3 += 1) {
                        for (int i4 = max(1, 128 * ii0 - 2 * i0 - i3 + 134); i4 <= 65; i4 += 1) {
                          B[64][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[65][i3][i4] - (SCALAR_VAL(2.0) * A[64][i3][i4])) + A[63][i3][i4])) + (SCALAR_VAL(0.125) * ((A[64][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[64][i3][i4])) + A[64][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[64][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[64][i3][i4])) + A[64][i3][i4 - 1]))) + A[64][i3][i4]);
                        }
                      }
                    }
                    for (int i3 = max(1, 128 * ii0 - 2 * i0 + 68); i3 <= 65; i3 += 1) {
                      for (int i4 = max(1, 128 * ii0 - 2 * i0 - i3 + 133); i4 <= 65; i4 += 1) {
                        B[65][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[66][i3][i4] - (SCALAR_VAL(2.0) * A[65][i3][i4])) + A[64][i3][i4])) + (SCALAR_VAL(0.125) * ((A[65][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[65][i3][i4])) + A[65][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[65][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[65][i3][i4])) + A[65][i3][i4 - 1]))) + A[65][i3][i4]);
                      }
                    }
                  } else {
                    for (int i2 = 64; i2 <= 65; i2 += 1) {
                      for (int i3 = -i2 + 129; i3 <= 65; i3 += 1) {
                        for (int i4 = max(-i2 + 129, -2 * i3 + 193); i4 <= 65; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    }
                  }
                } else if (ii3 == 0 && ii4 == 1 && i0 >= 64 * ii0 + 34) {
                  for (int i2 = 1; i2 <= 128 * ii0 - 2 * i0 + 131; i2 += 1) {
                    for (int i3 = 1; i3 <= min(2, 128 * ii0 - 2 * i0 - i2 + 132); i3 += 1) {
                      for (int i4 = 128 * ii0 - 2 * i0 - i2 + 133; i4 <= 65; i4 += 1) {
                        B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                      }
                    }
                    if (i0 == 64 * ii0 + 34 && i2 == 1) {
                      for (int i3 = 3; i3 <= 63; i3 += 1) {
                        for (int i4 = 64; i4 <= 65; i4 += 1) {
                          B[1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[2][i3][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[0][i3][i4])) + (SCALAR_VAL(0.125) * ((A[1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3][i4 - 1]))) + A[1][i3][i4]);
                        }
                      }
                    } else {
                      for (int i3 = 3; i3 <= 128 * ii0 - 2 * i0 - i2 + 132; i3 += 1) {
                        for (int i4 = 128 * ii0 - 2 * i0 - i2 + 133; i4 <= 65; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    }
                  }
                } else if (ii3 == 1 && ii4 == 0 && i0 >= 64 * ii0 + 34) {
                  for (int i2 = 1; i2 <= 65; i2 += 1) {
                    for (int i3 = max(1, 128 * ii0 - 2 * i0 - i2 + 133); i3 <= min(65, 128 * ii0 - 2 * i0 - i2 + 196); i3 += 1) {
                      if (i0 == 64 * ii0 + 34 && i3 == 1) {
                        B[i2][1][1] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][1][1] - (SCALAR_VAL(2.0) * A[i2][1][1])) + A[i2 - 1][1][1])) + (SCALAR_VAL(0.125) * ((A[i2][2][1] - (SCALAR_VAL(2.0) * A[i2][1][1])) + A[i2][0][1]))) + (SCALAR_VAL(0.125) * ((A[i2][1][2] - (SCALAR_VAL(2.0) * A[i2][1][1])) + A[i2][1][0]))) + A[i2][1][1]);
                      }
                      for (int i4 = max(1, 128 * ii0 - 2 * i0 - i3 + 71); i4 <= 128 * ii0 - 2 * i0 - i2 - i3 + 197; i4 += 1) {
                        B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                      }
                    }
                  }
                } else if (ii3 == 1 && ii4 == 0 && i0 >= 64 * ii0 + 2) {
                  if (i0 == 64 * ii0 + 3) {
                    for (int i2 = 62; i2 <= 65; i2 += 1) {
                      for (int i3 = -i2 + 127; i3 <= 64; i3 += 1) {
                        if (i3 <= 63) {
                          B[i2][i3][1] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][1] - (SCALAR_VAL(2.0) * A[i2][i3][1])) + A[i2 - 1][i3][1])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][1] - (SCALAR_VAL(2.0) * A[i2][i3][1])) + A[i2][i3 - 1][1]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][2] - (SCALAR_VAL(2.0) * A[i2][i3][1])) + A[i2][i3][0]))) + A[i2][i3][1]);
                          if (i2 == 65 && i3 == 62) {
                            for (int i4 = 2; i4 <= 64; i4 += 1) {
                              B[65][62][i4] = ((((SCALAR_VAL(0.125) * ((A[66][62][i4] - (SCALAR_VAL(2.0) * A[65][62][i4])) + A[64][62][i4])) + (SCALAR_VAL(0.125) * ((A[65][63][i4] - (SCALAR_VAL(2.0) * A[65][62][i4])) + A[65][61][i4]))) + (SCALAR_VAL(0.125) * ((A[65][62][i4 + 1] - (SCALAR_VAL(2.0) * A[65][62][i4])) + A[65][62][i4 - 1]))) + A[65][62][i4]);
                            }
                          }
                        }
                        if (i3 >= 63) {
                          for (int i4 = -i3 + 65; i4 <= -i2 - i3 + 191; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                      for (int i4 = 1; i4 <= -i2 + 126; i4 += 1) {
                        B[i2][65][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][65][i4] - (SCALAR_VAL(2.0) * A[i2][65][i4])) + A[i2 - 1][65][i4])) + (SCALAR_VAL(0.125) * ((A[i2][66][i4] - (SCALAR_VAL(2.0) * A[i2][65][i4])) + A[i2][64][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][65][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][65][i4])) + A[i2][65][i4 - 1]))) + A[i2][65][i4]);
                      }
                    }
                  } else {
                    for (int i2 = 64; i2 <= 65; i2 += 1) {
                      for (int i3 = -i2 + 129; i3 <= 65; i3 += 1) {
                        for (int i4 = 1; i4 <= -i2 + 128; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        if (i2 == 65 && i3 == 64) {
                          B[65][64][64] = ((((SCALAR_VAL(0.125) * ((A[66][64][64] - (SCALAR_VAL(2.0) * A[65][64][64])) + A[64][64][64])) + (SCALAR_VAL(0.125) * ((A[65][65][64] - (SCALAR_VAL(2.0) * A[65][64][64])) + A[65][63][64]))) + (SCALAR_VAL(0.125) * ((A[65][64][65] - (SCALAR_VAL(2.0) * A[65][64][64])) + A[65][64][63]))) + A[65][64][64]);
                        }
                      }
                    }
                  }
                } else if (ii3 == 0 && ii4 == 1 && i0 == 64 * ii0 + 2) {
                  for (int i2 = 64; i2 <= 65; i2 += 1) {
                    for (int i3 = 1; i3 <= -i2 + 128; i3 += 1) {
                      for (int i4 = -i2 + 129; i4 <= 65; i4 += 1) {
                        B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                      }
                    }
                  }
                }
                if (ii3 == 1) {
                  for (int i2 = max(1, 128 * ii0 - 2 * i0 + 67); i2 <= 64; i2 += 1) {
                    for (int i3 = max(1, 128 * ii0 - 2 * i0 - i2 + 132); i3 <= min(65, 128 * ii0 + 62 * ii4 - 2 * i0 - i2 + 195); i3 += 1) {
                      for (int i4 = max(1, 128 * ii0 + 64 * ii4 - 2 * i0 - i2 - i3 + 133); i4 <= min(65, 128 * ii0 + 126 * ii4 - 2 * i0 - i2 - i3 + 196); i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                  }
                  if (ii4 == 0 && i0 >= 64 * ii0 + 33) {
                    for (int i3 = 1; i3 <= 128 * ii0 - 2 * i0 + 130; i3 += 1) {
                      for (int i4 = 1; i4 <= 128 * ii0 - 2 * i0 - i3 + 131; i4 += 1) {
                        A[65][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[66][i3][i4] - (SCALAR_VAL(2.0) * B[65][i3][i4])) + B[64][i3][i4])) + (SCALAR_VAL(0.125) * ((B[65][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[65][i3][i4])) + B[65][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[65][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[65][i3][i4])) + B[65][i3][i4 - 1]))) + B[65][i3][i4]);
                      }
                    }
                  }
                } else if (ii4 == 1) {
                  for (int i2 = max(1, 128 * ii0 - 2 * i0 + 67); i2 <= min(64, 128 * ii0 - 2 * i0 + 130); i2 += 1) {
                    for (int i3 = 1; i3 <= 128 * ii0 - 2 * i0 - i2 + 131; i3 += 1) {
                      for (int i4 = 128 * ii0 - 2 * i0 - i2 + 132; i4 <= 65; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                  }
                }
                if (ii4 == ii3) {
                  for (int i2 = max(64 * ii3 + 1, 128 * ii0 - 2 * i0 + 67); i2 <= ii3 + 64; i2 += 1) {
                    for (int i3 = max(1, 128 * ii0 + 58 * ii3 - 2 * i0 + 9); i3 <= min(64, 128 * ii0 + 127 * ii3 - 2 * i0 - i2 + 131); i3 += 1) {
                      if (ii3 == 0) {
                        for (int i4 = 1; i4 <= 128 * ii0 - 2 * i0 - i2 + 131; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      } else {
                        for (int i4 = max(1, 128 * ii0 - 2 * i0 - i3 + 132); i4 <= 65; i4 += 1) {
                          A[65][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[66][i3][i4] - (SCALAR_VAL(2.0) * B[65][i3][i4])) + B[64][i3][i4])) + (SCALAR_VAL(0.125) * ((B[65][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[65][i3][i4])) + B[65][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[65][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[65][i3][i4])) + B[65][i3][i4 - 1]))) + B[65][i3][i4]);
                        }
                      }
                    }
                    if (ii3 == 1 && i2 == 65) {
                      for (int i4 = max(1, 128 * ii0 - 2 * i0 + 67); i4 <= 65; i4 += 1) {
                        A[65][65][i4] = ((((SCALAR_VAL(0.125) * ((B[66][65][i4] - (SCALAR_VAL(2.0) * B[65][65][i4])) + B[64][65][i4])) + (SCALAR_VAL(0.125) * ((B[65][66][i4] - (SCALAR_VAL(2.0) * B[65][65][i4])) + B[65][64][i4]))) + (SCALAR_VAL(0.125) * ((B[65][65][i4 + 1] - (SCALAR_VAL(2.0) * B[65][65][i4])) + B[65][65][i4 - 1]))) + B[65][65][i4]);
                      }
                    }
                  }
                }
                if (ii3 + ii4 <= 1 && 64 * ii0 + 32 >= i0) {
                  if (ii4 == 0) {
                    for (int i3 = max(1, 128 * ii0 + 62 * ii3 - 2 * i0 + 5); i3 <= min(64, 128 * ii0 + 62 * ii3 - 2 * i0 + 66); i3 += 1) {
                      for (int i4 = 1; i4 <= min(128 * ii0 + 62 * ii3 - 2 * i0 + 66, 128 * ii0 - 2 * i0 - i3 + 131); i4 += 1) {
                        A[65][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[66][i3][i4] - (SCALAR_VAL(2.0) * B[65][i3][i4])) + B[64][i3][i4])) + (SCALAR_VAL(0.125) * ((B[65][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[65][i3][i4])) + B[65][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[65][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[65][i3][i4])) + B[65][i3][i4 - 1]))) + B[65][i3][i4]);
                      }
                    }
                    if (ii3 == 1) {
                      for (int i4 = 1; i4 <= 128 * ii0 - 2 * i0 + 66; i4 += 1) {
                        A[65][65][i4] = ((((SCALAR_VAL(0.125) * ((B[66][65][i4] - (SCALAR_VAL(2.0) * B[65][65][i4])) + B[64][65][i4])) + (SCALAR_VAL(0.125) * ((B[65][66][i4] - (SCALAR_VAL(2.0) * B[65][65][i4])) + B[65][64][i4]))) + (SCALAR_VAL(0.125) * ((B[65][65][i4 + 1] - (SCALAR_VAL(2.0) * B[65][65][i4])) + B[65][65][i4 - 1]))) + B[65][65][i4]);
                      }
                    }
                  } else {
                    for (int i3 = 1; i3 <= 128 * ii0 - 2 * i0 + 66; i3 += 1) {
                      for (int i4 = 128 * ii0 - 2 * i0 + 67; i4 <= 65; i4 += 1) {
                        A[65][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[66][i3][i4] - (SCALAR_VAL(2.0) * B[65][i3][i4])) + B[64][i3][i4])) + (SCALAR_VAL(0.125) * ((B[65][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[65][i3][i4])) + B[65][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[65][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[65][i3][i4])) + B[65][i3][i4 - 1]))) + B[65][i3][i4]);
                      }
                    }
                  }
                }
              }
              if (ii3 == 0 && ii4 == 0) {
                for (int i0 = 64 * ii0 + 34; i0 <= min(_PB_TSTEPS, 64 * ii0 + 64); i0 += 1) {
                  for (int i2 = 1; i2 <= 128 * ii0 - 2 * i0 + 131; i2 += 1) {
                    for (int i4 = 1; i4 <= 128 * ii0 - 2 * i0 - i2 + 132; i4 += 1) {
                      B[i2][1][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][1][i4] - (SCALAR_VAL(2.0) * A[i2][1][i4])) + A[i2 - 1][1][i4])) + (SCALAR_VAL(0.125) * ((A[i2][2][i4] - (SCALAR_VAL(2.0) * A[i2][1][i4])) + A[i2][0][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][1][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][1][i4])) + A[i2][1][i4 - 1]))) + A[i2][1][i4]);
                    }
                    if (i0 == 64 * ii0 + 34 && i2 <= 2) {
                      for (int i4 = 1; i4 <= -i2 + 64; i4 += 1) {
                        B[i2][2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][2][i4] - (SCALAR_VAL(2.0) * A[i2][2][i4])) + A[i2 - 1][2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][3][i4] - (SCALAR_VAL(2.0) * A[i2][2][i4])) + A[i2][1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][2][i4])) + A[i2][2][i4 - 1]))) + A[i2][2][i4]);
                      }
                    }
                    if (2 * i0 + i2 >= 128 * ii0 + 70) {
                      for (int i3 = max(2, 128 * ii0 - 2 * i0 - i2 + 73); i3 <= 128 * ii0 - 2 * i0 - i2 + 132; i3 += 1) {
                        for (int i4 = 1; i4 <= 128 * ii0 - 2 * i0 - i2 + 132; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    } else {
                      for (int i3 = 3; i3 <= 63; i3 += 1) {
                        for (int i4 = 1; i4 <= 63; i4 += 1) {
                          B[1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[2][i3][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[0][i3][i4])) + (SCALAR_VAL(0.125) * ((A[1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3][i4 - 1]))) + A[1][i3][i4]);
                        }
                      }
                    }
                  }
                  for (int i2 = 1; i2 <= 128 * ii0 - 2 * i0 + 130; i2 += 1) {
                    for (int i3 = 1; i3 <= 128 * ii0 - 2 * i0 - i2 + 131; i3 += 1) {
                      for (int i4 = 1; i4 <= 128 * ii0 - 2 * i0 - i2 + 131; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                  }
                }
              }
            } else if (ii3 >= ii4) {
              if (ii3 == 0 && ii4 == 0) {
                for (int i0 = 64 * ii0 + 1; i0 <= 64 * ii0 + 2; i0 += 1) {
                  if (i0 == 64 * ii0 + 2) {
                    for (int i2 = 64; i2 <= 66; i2 += 1) {
                      for (int i3 = 1; i3 <= -i2 + 128; i3 += 1) {
                        for (int i4 = 1; i4 <= -i2 + 128; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        if (i2 == 66) {
                          B[66][i3][63] = ((((SCALAR_VAL(0.125) * ((A[67][i3][63] - (SCALAR_VAL(2.0) * A[66][i3][63])) + A[65][i3][63])) + (SCALAR_VAL(0.125) * ((A[66][i3 + 1][63] - (SCALAR_VAL(2.0) * A[66][i3][63])) + A[66][i3 - 1][63]))) + (SCALAR_VAL(0.125) * ((A[66][i3][64] - (SCALAR_VAL(2.0) * A[66][i3][63])) + A[66][i3][62]))) + A[66][i3][63]);
                        }
                      }
                      if (i2 == 66) {
                        for (int i4 = 1; i4 <= 63; i4 += 1) {
                          B[66][63][i4] = ((((SCALAR_VAL(0.125) * ((A[67][63][i4] - (SCALAR_VAL(2.0) * A[66][63][i4])) + A[65][63][i4])) + (SCALAR_VAL(0.125) * ((A[66][64][i4] - (SCALAR_VAL(2.0) * A[66][63][i4])) + A[66][62][i4]))) + (SCALAR_VAL(0.125) * ((A[66][63][i4 + 1] - (SCALAR_VAL(2.0) * A[66][63][i4])) + A[66][63][i4 - 1]))) + A[66][63][i4]);
                        }
                      }
                    }
                  }
                  if (i0 == 64 * ii0 + 2) {
                    for (int i2 = 63; i2 <= 64; i2 += 1) {
                      for (int i3 = 1; i3 <= -i2 + 127; i3 += 1) {
                        for (int i4 = 1; i4 <= -i2 + 127; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                  }
                  for (int i2 = 65; i2 <= 66; i2 += 1) {
                    for (int i3 = 1; i3 <= 128 * ii0 - 2 * i0 + 66; i3 += 1) {
                      for (int i4 = 1; i4 <= 128 * ii0 - 2 * i0 + 66; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                  }
                }
              } else {
                for (int i0 = 64 * ii0 + 1; i0 <= 64 * ii0 - ii4 + 3; i0 += 1) {
                  if (i0 >= 64 * ii0 + 2) {
                    if (ii4 == 0 && i0 == 64 * ii0 + 3) {
                      for (int i3 = 65; i3 <= 66; i3 += 1) {
                        for (int i4 = 1; i4 <= 64; i4 += 1) {
                          B[62][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[63][i3][i4] - (SCALAR_VAL(2.0) * A[62][i3][i4])) + A[61][i3][i4])) + (SCALAR_VAL(0.125) * ((A[62][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[62][i3][i4])) + A[62][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[62][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[62][i3][i4])) + A[62][i3][i4 - 1]))) + A[62][i3][i4]);
                        }
                      }
                    }
                    for (int i2 = max(ii4 + 63, 64 * ii0 - ii4 - i0 + 66); i2 <= 66; i2 += 1) {
                      if (ii4 == 0 && i0 == 64 * ii0 + 3) {
                        for (int i3 = max(62, -i2 + 127); i3 <= -i2 + i2 / 2 + 96; i3 += 1) {
                          if (i3 <= 63) {
                            B[i2][i3][1] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][1] - (SCALAR_VAL(2.0) * A[i2][i3][1])) + A[i2 - 1][i3][1])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][1] - (SCALAR_VAL(2.0) * A[i2][i3][1])) + A[i2][i3 - 1][1]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][2] - (SCALAR_VAL(2.0) * A[i2][i3][1])) + A[i2][i3][0]))) + A[i2][i3][1]);
                          }
                          if (i3 >= 63) {
                            for (int i4 = -i3 + 65; i4 <= -i2 - i3 + 191; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            if (i2 == 66 && i3 == 63) {
                              B[66][63][63] = ((((SCALAR_VAL(0.125) * ((A[67][63][63] - (SCALAR_VAL(2.0) * A[66][63][63])) + A[65][63][63])) + (SCALAR_VAL(0.125) * ((A[66][64][63] - (SCALAR_VAL(2.0) * A[66][63][63])) + A[66][62][63]))) + (SCALAR_VAL(0.125) * ((A[66][63][64] - (SCALAR_VAL(2.0) * A[66][63][63])) + A[66][63][62]))) + A[66][63][63]);
                            }
                          } else {
                            for (int i4 = 2; i4 <= 64; i4 += 1) {
                              B[i2][62][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][62][i4] - (SCALAR_VAL(2.0) * A[i2][62][i4])) + A[i2 - 1][62][i4])) + (SCALAR_VAL(0.125) * ((A[i2][63][i4] - (SCALAR_VAL(2.0) * A[i2][62][i4])) + A[i2][61][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][62][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][62][i4])) + A[i2][62][i4 - 1]))) + A[i2][62][i4]);
                            }
                          }
                        }
                      }
                      if (ii4 == 0 && i2 >= 65) {
                        if (i0 == 64 * ii0 + 3) {
                          for (int i4 = 1; i4 <= 62; i4 += 1) {
                            B[i2][64][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64][i4] - (SCALAR_VAL(2.0) * A[i2][64][i4])) + A[i2 - 1][64][i4])) + (SCALAR_VAL(0.125) * ((A[i2][65][i4] - (SCALAR_VAL(2.0) * A[i2][64][i4])) + A[i2][63][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][64][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][64][i4])) + A[i2][64][i4 - 1]))) + A[i2][64][i4]);
                          }
                        } else {
                          for (int i4 = 1; i4 <= 64; i4 += 1) {
                            B[i2][64][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64][i4] - (SCALAR_VAL(2.0) * A[i2][64][i4])) + A[i2 - 1][64][i4])) + (SCALAR_VAL(0.125) * ((A[i2][65][i4] - (SCALAR_VAL(2.0) * A[i2][64][i4])) + A[i2][63][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][64][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][64][i4])) + A[i2][64][i4 - 1]))) + A[i2][64][i4]);
                          }
                        }
                      }
                      if (ii4 == 0 && i2 <= 65) {
                        for (int i3 = 65; i3 <= 66; i3 += 1) {
                          if (i0 == 64 * ii0 + 3) {
                            for (int i4 = 1; i4 <= -i2 + 126; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          } else {
                            for (int i4 = 1; i4 <= -i2 + 128; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                        }
                      } else {
                        for (int i3 = max(-ii4 + 65, -i2 + 129); i3 <= 66; i3 += 1) {
                          if (ii4 == 1 && i0 == 64 * ii0 + 2) {
                            for (int i4 = max(max(64, -i2 + 129), -2 * i3 + 193); i4 <= 66; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          } else if (i0 == 64 * ii0 + 3) {
                            for (int i4 = 1; i4 <= 61; i4 += 1) {
                              B[66][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[67][i3][i4] - (SCALAR_VAL(2.0) * A[66][i3][i4])) + A[65][i3][i4])) + (SCALAR_VAL(0.125) * ((A[66][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[66][i3][i4])) + A[66][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[66][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[66][i3][i4])) + A[66][i3][i4 - 1]))) + A[66][i3][i4]);
                            }
                          } else {
                            for (int i4 = 1; i4 <= 63; i4 += 1) {
                              B[66][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[67][i3][i4] - (SCALAR_VAL(2.0) * A[66][i3][i4])) + A[65][i3][i4])) + (SCALAR_VAL(0.125) * ((A[66][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[66][i3][i4])) + A[66][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[66][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[66][i3][i4])) + A[66][i3][i4 - 1]))) + A[66][i3][i4]);
                            }
                          }
                        }
                      }
                    }
                  }
                  for (int i2 = 128 * ii0 - 2 * i0 + 67; i2 <= 64; i2 += 1) {
                    for (int i3 = 128 * ii0 - 2 * i0 - i2 + 132; i3 <= 66; i3 += 1) {
                      for (int i4 = max(max(1, 66 * ii4 - i2 + 62), 66 * ii4 - i2 - i3 + 127); i4 <= min(66, 128 * ii0 + 4 * ii4 - 2 * i0 - i2 - i3 + 196); i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                      if (ii4 == 0 && i3 == 66) {
                        A[i2][66][128 * ii0 - 2 * i0 - i2 + 131] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][66][128 * ii0 - 2 * i0 - i2 + 131] - (SCALAR_VAL(2.0) * B[i2][66][128 * ii0 - 2 * i0 - i2 + 131])) + B[i2 - 1][66][128 * ii0 - 2 * i0 - i2 + 131])) + (SCALAR_VAL(0.125) * ((B[i2][67][128 * ii0 - 2 * i0 - i2 + 131] - (SCALAR_VAL(2.0) * B[i2][66][128 * ii0 - 2 * i0 - i2 + 131])) + B[i2][65][128 * ii0 - 2 * i0 - i2 + 131]))) + (SCALAR_VAL(0.125) * ((B[i2][66][128 * ii0 - 2 * i0 - i2 + 132] - (SCALAR_VAL(2.0) * B[i2][66][128 * ii0 - 2 * i0 - i2 + 131])) + B[i2][66][128 * ii0 - 2 * i0 - i2 + 130]))) + B[i2][66][128 * ii0 - 2 * i0 - i2 + 131]);
                      }
                    }
                  }
                  for (int i2 = 65; i2 <= 66; i2 += 1) {
                    if (ii4 == 1 && i0 == 64 * ii0 + 2) {
                      for (int i3 = 63; i3 <= 64; i3 += 1) {
                        for (int i4 = -i3 + 128; i4 <= 66; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    } else if (ii4 == 0) {
                      for (int i3 = 128 * ii0 - 2 * i0 + 67; i3 <= 64; i3 += 1) {
                        for (int i4 = 1; i4 <= 128 * ii0 - 2 * i0 - i3 + 131; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                    for (int i3 = 65; i3 <= 66; i3 += 1) {
                      if (i0 >= 64 * ii0 + ii4 + 1) {
                        for (int i4 = 62 * ii4 + 1; i4 <= 128 * ii0 + 2 * ii4 - 2 * i0 + 66; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                      if (ii4 == 1) {
                        for (int i4 = 65; i4 <= 66; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                  }
                }
              }
              for (int i0 = 64 * ii0 + ii3 - ii4 + 3; i0 <= min(_PB_TSTEPS, 64 * ii0 + 31 * ii4 + 33); i0 += 1) {
                if (ii3 == 1 && ii4 == 0) {
                  for (int i3 = 65; i3 <= 66; i3 += 1) {
                    for (int i4 = 1; i4 <= 64; i4 += 1) {
                      B[128 * ii0 - 2 * i0 + 68][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[128 * ii0 - 2 * i0 + 69][i3][i4] - (SCALAR_VAL(2.0) * A[128 * ii0 - 2 * i0 + 68][i3][i4])) + A[128 * ii0 - 2 * i0 + 67][i3][i4])) + (SCALAR_VAL(0.125) * ((A[128 * ii0 - 2 * i0 + 68][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[128 * ii0 - 2 * i0 + 68][i3][i4])) + A[128 * ii0 - 2 * i0 + 68][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[128 * ii0 - 2 * i0 + 68][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[128 * ii0 - 2 * i0 + 68][i3][i4])) + A[128 * ii0 - 2 * i0 + 68][i3][i4 - 1]))) + A[128 * ii0 - 2 * i0 + 68][i3][i4]);
                    }
                  }
                  for (int i2 = 128 * ii0 - 2 * i0 + 69; i2 <= 64; i2 += 1) {
                    for (int i3 = 128 * ii0 - 2 * i0 - i2 + 133; i3 <= 66; i3 += 1) {
                      if (i3 >= 65) {
                        for (int i4 = 1; i4 <= 128 * ii0 - 2 * i0 - i2 + 132; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      } else {
                        if (i2 == 64 && 2 * i0 + i3 == 128 * ii0 + 69) {
                          B[64][128 * ii0 - 2 * i0 + 69][1] = ((((SCALAR_VAL(0.125) * ((A[65][128 * ii0 - 2 * i0 + 69][1] - (SCALAR_VAL(2.0) * A[64][128 * ii0 - 2 * i0 + 69][1])) + A[63][128 * ii0 - 2 * i0 + 69][1])) + (SCALAR_VAL(0.125) * ((A[64][128 * ii0 - 2 * i0 + 70][1] - (SCALAR_VAL(2.0) * A[64][128 * ii0 - 2 * i0 + 69][1])) + A[64][128 * ii0 - 2 * i0 + 68][1]))) + (SCALAR_VAL(0.125) * ((A[64][128 * ii0 - 2 * i0 + 69][2] - (SCALAR_VAL(2.0) * A[64][128 * ii0 - 2 * i0 + 69][1])) + A[64][128 * ii0 - 2 * i0 + 69][0]))) + A[64][128 * ii0 - 2 * i0 + 69][1]);
                        }
                        for (int i4 = max(1, 128 * ii0 - 2 * i0 - i3 + 71); i4 <= 128 * ii0 - 2 * i0 - i2 - i3 + 197; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    }
                  }
                } else if (ii3 == 1 && ii4 == 1) {
                  for (int i2 = max(1, 128 * ii0 - 2 * i0 + 68); i2 <= 96 * ii0 - 2 * i0 + (i0 + 1) / 2 + 66; i2 += 1) {
                    for (int i3 = 128 * ii0 - 2 * i0 - i2 + 133; i3 <= 66; i3 += 1) {
                      for (int i4 = max(128 * ii0 - 2 * i0 - i2 + 133, 128 * ii0 - 2 * i0 - i2 - i3 + 198); i4 <= 66; i4 += 1) {
                        B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                      }
                    }
                  }
                }
                for (int i2 = max(max(63 * ii3 - 64 * ii4 + 2, 128 * ii0 - 2 * i0 + 68), 96 * ii0 + 15 * ii4 - 2 * i0 + (ii4 + i0) / 2 + 52); i2 <= min(min(66, 128 * ii0 + ii3 + 42 * ii4 - 2 * i0 + 131), 43 * ii0 - i0 + (-ii0 + i0 - 1) / 3 + 89); i2 += 1) {
                  if (ii3 == 1) {
                    if (ii4 == 0) {
                      for (int i4 = 1; i4 <= 64; i4 += 1) {
                        B[i2][128 * ii0 - 2 * i0 + 68][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][128 * ii0 - 2 * i0 + 68][i4] - (SCALAR_VAL(2.0) * A[i2][128 * ii0 - 2 * i0 + 68][i4])) + A[i2 - 1][128 * ii0 - 2 * i0 + 68][i4])) + (SCALAR_VAL(0.125) * ((A[i2][128 * ii0 - 2 * i0 + 69][i4] - (SCALAR_VAL(2.0) * A[i2][128 * ii0 - 2 * i0 + 68][i4])) + A[i2][128 * ii0 - 2 * i0 + 67][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][128 * ii0 - 2 * i0 + 68][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][128 * ii0 - 2 * i0 + 68][i4])) + A[i2][128 * ii0 - 2 * i0 + 68][i4 - 1]))) + A[i2][128 * ii0 - 2 * i0 + 68][i4]);
                      }
                      for (int i3 = 128 * ii0 - 2 * i0 + 69; i3 <= 63; i3 += 1) {
                        if (2 * i0 + i3 == 128 * ii0 + 69) {
                          B[i2][128 * ii0 - 2 * i0 + 69][1] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][128 * ii0 - 2 * i0 + 69][1] - (SCALAR_VAL(2.0) * A[i2][128 * ii0 - 2 * i0 + 69][1])) + A[i2 - 1][128 * ii0 - 2 * i0 + 69][1])) + (SCALAR_VAL(0.125) * ((A[i2][128 * ii0 - 2 * i0 + 70][1] - (SCALAR_VAL(2.0) * A[i2][128 * ii0 - 2 * i0 + 69][1])) + A[i2][128 * ii0 - 2 * i0 + 68][1]))) + (SCALAR_VAL(0.125) * ((A[i2][128 * ii0 - 2 * i0 + 69][2] - (SCALAR_VAL(2.0) * A[i2][128 * ii0 - 2 * i0 + 69][1])) + A[i2][128 * ii0 - 2 * i0 + 69][0]))) + A[i2][128 * ii0 - 2 * i0 + 69][1]);
                        }
                        for (int i4 = max(1, 128 * ii0 - 2 * i0 - i3 + 71); i4 <= 128 * ii0 - 2 * i0 - i2 - i3 + 197; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        if (i2 == 66) {
                          B[66][i3][128 * ii0 - 2 * i0 - i3 + 132] = ((((SCALAR_VAL(0.125) * ((A[67][i3][128 * ii0 - 2 * i0 - i3 + 132] - (SCALAR_VAL(2.0) * A[66][i3][128 * ii0 - 2 * i0 - i3 + 132])) + A[65][i3][128 * ii0 - 2 * i0 - i3 + 132])) + (SCALAR_VAL(0.125) * ((A[66][i3 + 1][128 * ii0 - 2 * i0 - i3 + 132] - (SCALAR_VAL(2.0) * A[66][i3][128 * ii0 - 2 * i0 - i3 + 132])) + A[66][i3 - 1][128 * ii0 - 2 * i0 - i3 + 132]))) + (SCALAR_VAL(0.125) * ((A[66][i3][128 * ii0 - 2 * i0 - i3 + 133] - (SCALAR_VAL(2.0) * A[66][i3][128 * ii0 - 2 * i0 - i3 + 132])) + A[66][i3][128 * ii0 - 2 * i0 - i3 + 131]))) + A[66][i3][128 * ii0 - 2 * i0 - i3 + 132]);
                        }
                      }
                    }
                    if (i2 >= ii4 + 65) {
                      for (int i3 = max(-63 * ii4 + 64, 128 * ii0 - 2 * i0 + 68); i3 <= min(128 * ii0 - 2 * i0 + 132, ii4 + 64); i3 += 1) {
                        for (int i4 = max(1, 128 * ii0 + 60 * ii4 - 2 * i0 - i3 + 73); i4 <= min(min(65, 128 * ii0 + 65 * ii4 - 2 * i0 + 68), 128 * ii0 - 2 * i0 - i3 + 134); i4 += 1) {
                          if (128 * ii0 + 133 >= 2 * i0 + i3 + i4 || 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                        if (ii4 == 1 && i2 == 66 && 128 * ii0 + 69 >= 2 * i0 + i3) {
                          B[66][i3][66] = ((((SCALAR_VAL(0.125) * ((A[67][i3][66] - (SCALAR_VAL(2.0) * A[66][i3][66])) + A[65][i3][66])) + (SCALAR_VAL(0.125) * ((A[66][i3 + 1][66] - (SCALAR_VAL(2.0) * A[66][i3][66])) + A[66][i3 - 1][66]))) + (SCALAR_VAL(0.125) * ((A[66][i3][67] - (SCALAR_VAL(2.0) * A[66][i3][66])) + A[66][i3][65]))) + A[66][i3][66]);
                        } else if (ii4 == 1 && i2 == 66) {
                          for (int i4 = 128 * ii0 - 2 * i0 - i3 + 135; i4 <= 66; i4 += 1) {
                            B[66][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[67][i3][i4] - (SCALAR_VAL(2.0) * A[66][i3][i4])) + A[65][i3][i4])) + (SCALAR_VAL(0.125) * ((A[66][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[66][i3][i4])) + A[66][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[66][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[66][i3][i4])) + A[66][i3][i4 - 1]))) + A[66][i3][i4]);
                          }
                        }
                      }
                      if (ii4 == 0 && i2 == 65) {
                        for (int i3 = 65; i3 <= 66; i3 += 1) {
                          for (int i4 = 1; i4 <= 128 * ii0 - 2 * i0 + 67; i4 += 1) {
                            B[65][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[66][i3][i4] - (SCALAR_VAL(2.0) * A[65][i3][i4])) + A[64][i3][i4])) + (SCALAR_VAL(0.125) * ((A[65][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[65][i3][i4])) + A[65][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[65][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[65][i3][i4])) + A[65][i3][i4 - 1]))) + A[65][i3][i4]);
                          }
                        }
                      }
                    } else {
                      for (int i3 = max(1, 128 * ii0 - 2 * i0 - i2 + 133); i3 <= min(min(min(128 * ii0 - 2 * i0 + 132, -i2 + 128), 128 * ii0 - 2 * i0 - 2 * i2 + i2 / 3 + 179), (i2 + 1) / 3 + 43); i3 += 1) {
                        for (int i4 = 128 * ii0 - 2 * i0 - i2 - i3 + 198; i4 <= 66; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                      for (int i3 = max(1, 128 * ii0 - 2 * i0 - 2 * i2 + i2 / 3 + 180); i3 <= min(min(128 * ii0 - 2 * i0 + 132, -i2 + 128), (i2 + 1) / 3 + 43); i3 += 1) {
                        for (int i4 = 128 * ii0 - 2 * i0 - i2 - i3 + 198; i4 <= 66; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                      if (i2 == 65) {
                        for (int i3 = 64; i3 <= min(65, 128 * ii0 - 2 * i0 + 132); i3 += 1) {
                          for (int i4 = 128 * ii0 - 2 * i0 - i3 + 133; i4 <= 66; i4 += 1) {
                            B[65][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[66][i3][i4] - (SCALAR_VAL(2.0) * A[65][i3][i4])) + A[64][i3][i4])) + (SCALAR_VAL(0.125) * ((A[65][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[65][i3][i4])) + A[65][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[65][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[65][i3][i4])) + A[65][i3][i4 - 1]))) + A[65][i3][i4]);
                          }
                        }
                      }
                    }
                    if (ii4 == 1 && i0 == 64 * ii0 + 34 && i2 >= 65) {
                      for (int i4 = 1; i4 <= 66; i4 += 1) {
                        B[i2][65][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][65][i4] - (SCALAR_VAL(2.0) * A[i2][65][i4])) + A[i2 - 1][65][i4])) + (SCALAR_VAL(0.125) * ((A[i2][66][i4] - (SCALAR_VAL(2.0) * A[i2][65][i4])) + A[i2][64][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][65][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][65][i4])) + A[i2][65][i4 - 1]))) + A[i2][65][i4]);
                      }
                    }
                    if (ii4 == 1 && i2 == 66) {
                      for (int i4 = max(1, 128 * ii0 - 2 * i0 + 68); i4 <= 66; i4 += 1) {
                        B[66][66][i4] = ((((SCALAR_VAL(0.125) * ((A[67][66][i4] - (SCALAR_VAL(2.0) * A[66][66][i4])) + A[65][66][i4])) + (SCALAR_VAL(0.125) * ((A[66][67][i4] - (SCALAR_VAL(2.0) * A[66][66][i4])) + A[66][65][i4]))) + (SCALAR_VAL(0.125) * ((A[66][66][i4 + 1] - (SCALAR_VAL(2.0) * A[66][66][i4])) + A[66][66][i4 - 1]))) + A[66][66][i4]);
                      }
                    }
                    for (int i3 = 128 * ii0 - 2 * i0 + 133; i3 <= min(128 * ii0 - 2 * i0 - 2 * i2 + i2 / 3 + 179, (i2 + 1) / 3 + 43); i3 += 1) {
                      for (int i4 = 128 * ii0 - 2 * i0 - i2 - i3 + 198; i4 <= 66; i4 += 1) {
                        B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                      }
                    }
                    for (int i3 = max(128 * ii0 - 2 * i0 + 133, 128 * ii0 - 2 * i0 - 2 * i2 + i2 / 3 + 180); i3 <= min(-i2 + 128, (i2 + 1) / 3 + 43); i3 += 1) {
                      for (int i4 = max(1, 128 * ii0 - 2 * i0 - i2 - i3 + 198); i4 <= 66; i4 += 1) {
                        B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                      }
                    }
                    if (i2 <= 63) {
                      for (int i3 = (i2 + 1) / 3 + 44; i3 <= min(66, -128 * ii0 + 2 * i0 + 2 * i2 - (i2 + 3) / 3 - 110); i3 += 1) {
                        if (2 * i0 + i2 + i3 == 128 * ii0 + 199) {
                          B[i2][128 * ii0 - 2 * i0 - i2 + 199][1] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][128 * ii0 - 2 * i0 - i2 + 199][1] - (SCALAR_VAL(2.0) * A[i2][128 * ii0 - 2 * i0 - i2 + 199][1])) + A[i2 - 1][128 * ii0 - 2 * i0 - i2 + 199][1])) + (SCALAR_VAL(0.125) * ((A[i2][128 * ii0 - 2 * i0 - i2 + 200][1] - (SCALAR_VAL(2.0) * A[i2][128 * ii0 - 2 * i0 - i2 + 199][1])) + A[i2][128 * ii0 - 2 * i0 - i2 + 198][1]))) + (SCALAR_VAL(0.125) * ((A[i2][128 * ii0 - 2 * i0 - i2 + 199][2] - (SCALAR_VAL(2.0) * A[i2][128 * ii0 - 2 * i0 - i2 + 199][1])) + A[i2][128 * ii0 - 2 * i0 - i2 + 199][0]))) + A[i2][128 * ii0 - 2 * i0 - i2 + 199][1]);
                        }
                        for (int i4 = max(1, 128 * ii0 - 2 * i0 - i2 - i3 + 201); i4 <= 66; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    } else {
                      if (ii4 == 1 && i0 == 64 * ii0 + 35 && i2 == 65) {
                        for (int i3 = 64; i3 <= 65; i3 += 1) {
                          if (i3 == 64) {
                            B[65][64][1] = ((((SCALAR_VAL(0.125) * ((A[66][64][1] - (SCALAR_VAL(2.0) * A[65][64][1])) + A[64][64][1])) + (SCALAR_VAL(0.125) * ((A[65][65][1] - (SCALAR_VAL(2.0) * A[65][64][1])) + A[65][63][1]))) + (SCALAR_VAL(0.125) * ((A[65][64][2] - (SCALAR_VAL(2.0) * A[65][64][1])) + A[65][64][0]))) + A[65][64][1]);
                          }
                          for (int i4 = -i3 + 66; i4 <= 66; i4 += 1) {
                            B[65][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[66][i3][i4] - (SCALAR_VAL(2.0) * A[65][i3][i4])) + A[64][i3][i4])) + (SCALAR_VAL(0.125) * ((A[65][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[65][i3][i4])) + A[65][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[65][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[65][i3][i4])) + A[65][i3][i4 - 1]))) + A[65][i3][i4]);
                          }
                        }
                      }
                      if (ii4 == 1 && i0 + i2 >= 64 * ii0 + 99) {
                        for (int i3 = max(2 * i2 - (i2 + 3) / 3 - 42, (i2 + 1) / 3 + 44); i3 <= 66; i3 += 1) {
                          if (i0 + i2 == 64 * ii0 + 99 && i0 + i3 == 64 * ii0 + 100) {
                            B[64 * ii0 - i0 + 99][64 * ii0 - i0 + 100][1] = ((((SCALAR_VAL(0.125) * ((A[64 * ii0 - i0 + 100][64 * ii0 - i0 + 100][1] - (SCALAR_VAL(2.0) * A[64 * ii0 - i0 + 99][64 * ii0 - i0 + 100][1])) + A[64 * ii0 - i0 + 98][64 * ii0 - i0 + 100][1])) + (SCALAR_VAL(0.125) * ((A[64 * ii0 - i0 + 99][64 * ii0 - i0 + 101][1] - (SCALAR_VAL(2.0) * A[64 * ii0 - i0 + 99][64 * ii0 - i0 + 100][1])) + A[64 * ii0 - i0 + 99][64 * ii0 - i0 + 99][1]))) + (SCALAR_VAL(0.125) * ((A[64 * ii0 - i0 + 99][64 * ii0 - i0 + 100][2] - (SCALAR_VAL(2.0) * A[64 * ii0 - i0 + 99][64 * ii0 - i0 + 100][1])) + A[64 * ii0 - i0 + 99][64 * ii0 - i0 + 100][0]))) + A[64 * ii0 - i0 + 99][64 * ii0 - i0 + 100][1]);
                          }
                          for (int i4 = max(1, 128 * ii0 - 2 * i0 - i2 - i3 + 201); i4 <= 66; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                    }
                    if (65 * ii4 + i2 >= 66) {
                      for (int i3 = max(max(ii4 + (-ii4 + i2 + 2) / 3 + 43, -128 * ii0 + ii4 + 2 * i0 + 2 * i2 - (i2 + 3) / 3 - 110), 2 * ii4 + 2 * i2 - (i2 + 3) / 3 - 44); i3 <= 66; i3 += 1) {
                        for (int i4 = max(max(1, 128 * ii0 + 58 * ii4 - 2 * i0 - i2 + 75), 128 * ii0 + 58 * ii4 - 2 * i0 - i2 - i3 + 140); i4 <= min(min(66, 128 * ii0 + 108 * ii4 - 2 * i0 + 67), 128 * ii0 - ii4 - 2 * i0 - 2 * i2 + i3 + i2 / 3 + 112); i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        if (ii4 == 1) {
                          for (int i4 = 128 * ii0 - 2 * i0 - 2 * i2 + i3 + i2 / 3 + 112; i4 <= min(128 * ii0 - 2 * i0 - i2 - i3 + 200, i2 - (i2 + 3) / 3 + 22); i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          for (int i4 = max(128 * ii0 - 2 * i0 - i2 - i3 + 201, 128 * ii0 - 2 * i0 - 2 * i2 + i3 + i2 / 3 + 112); i4 <= min(min(66, -i2 + 128), -64 * ii0 + i0 + i2 + (i2 + i3) / 2 - 68); i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (i2 == 63 && 2 * i0 + i3 >= 128 * ii0 + 79) {
                            B[63][i3][66] = ((((SCALAR_VAL(0.125) * ((A[64][i3][66] - (SCALAR_VAL(2.0) * A[63][i3][66])) + A[62][i3][66])) + (SCALAR_VAL(0.125) * ((A[63][i3 + 1][66] - (SCALAR_VAL(2.0) * A[63][i3][66])) + A[63][i3 - 1][66]))) + (SCALAR_VAL(0.125) * ((A[63][i3][67] - (SCALAR_VAL(2.0) * A[63][i3][66])) + A[63][i3][65]))) + A[63][i3][66]);
                          }
                          if (i2 <= 63) {
                            for (int i4 = max(max(i2 - (i2 + 3) / 3 + 23, 128 * ii0 - 2 * i0 - 2 * i2 + i3 + i2 / 3 + 112), -64 * ii0 + i0 + i2 + (i2 + i3) / 2 - 67); i4 <= 66; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          } else {
                            for (int i4 = -i2 + 129; i4 <= 66; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                        }
                      }
                    }
                  } else if (i2 <= 65) {
                    if (128 * ii0 + 70 >= 2 * i0 + i2) {
                      for (int i3 = 1; i3 <= 2; i3 += 1) {
                        if (i3 == 2) {
                          B[i2][2][1] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][2][1] - (SCALAR_VAL(2.0) * A[i2][2][1])) + A[i2 - 1][2][1])) + (SCALAR_VAL(0.125) * ((A[i2][3][1] - (SCALAR_VAL(2.0) * A[i2][2][1])) + A[i2][1][1]))) + (SCALAR_VAL(0.125) * ((A[i2][2][2] - (SCALAR_VAL(2.0) * A[i2][2][1])) + A[i2][2][0]))) + A[i2][2][1]);
                        } else {
                          B[i2][1][1] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][1][1] - (SCALAR_VAL(2.0) * A[i2][1][1])) + A[i2 - 1][1][1])) + (SCALAR_VAL(0.125) * ((A[i2][2][1] - (SCALAR_VAL(2.0) * A[i2][1][1])) + A[i2][0][1]))) + (SCALAR_VAL(0.125) * ((A[i2][1][2] - (SCALAR_VAL(2.0) * A[i2][1][1])) + A[i2][1][0]))) + A[i2][1][1]);
                        }
                        for (int i4 = 2; i4 <= 128 * ii0 - 2 * i0 - i2 + 132; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                      if (128 * ii0 + 69 >= 2 * i0 + i2) {
                        for (int i3 = 3; i3 <= (i2 - 1) / 3 + 43; i3 += 1) {
                          for (int i4 = 1; i4 <= 256 * ii0 - 4 * i0 - 2 * i2 + 200; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (2 * i0 + i2 == 128 * ii0 + 69) {
                            B[128 * ii0 - 2 * i0 + 69][i3][63] = ((((SCALAR_VAL(0.125) * ((A[128 * ii0 - 2 * i0 + 70][i3][63] - (SCALAR_VAL(2.0) * A[128 * ii0 - 2 * i0 + 69][i3][63])) + A[128 * ii0 - 2 * i0 + 68][i3][63])) + (SCALAR_VAL(0.125) * ((A[128 * ii0 - 2 * i0 + 69][i3 + 1][63] - (SCALAR_VAL(2.0) * A[128 * ii0 - 2 * i0 + 69][i3][63])) + A[128 * ii0 - 2 * i0 + 69][i3 - 1][63]))) + (SCALAR_VAL(0.125) * ((A[128 * ii0 - 2 * i0 + 69][i3][64] - (SCALAR_VAL(2.0) * A[128 * ii0 - 2 * i0 + 69][i3][63])) + A[128 * ii0 - 2 * i0 + 69][i3][62]))) + A[128 * ii0 - 2 * i0 + 69][i3][63]);
                          }
                        }
                        for (int i3 = (i2 - 1) / 3 + 44; i3 <= 128 * ii0 - 2 * i0 - i2 + 132; i3 += 1) {
                          for (int i4 = 1; i4 <= 256 * ii0 - 4 * i0 - 2 * i2 + 200; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (2 * i0 + i2 == 128 * ii0 + 69) {
                            B[128 * ii0 - 2 * i0 + 69][i3][63] = ((((SCALAR_VAL(0.125) * ((A[128 * ii0 - 2 * i0 + 70][i3][63] - (SCALAR_VAL(2.0) * A[128 * ii0 - 2 * i0 + 69][i3][63])) + A[128 * ii0 - 2 * i0 + 68][i3][63])) + (SCALAR_VAL(0.125) * ((A[128 * ii0 - 2 * i0 + 69][i3 + 1][63] - (SCALAR_VAL(2.0) * A[128 * ii0 - 2 * i0 + 69][i3][63])) + A[128 * ii0 - 2 * i0 + 69][i3 - 1][63]))) + (SCALAR_VAL(0.125) * ((A[128 * ii0 - 2 * i0 + 69][i3][64] - (SCALAR_VAL(2.0) * A[128 * ii0 - 2 * i0 + 69][i3][63])) + A[128 * ii0 - 2 * i0 + 69][i3][62]))) + A[128 * ii0 - 2 * i0 + 69][i3][63]);
                          }
                        }
                      }
                    } else {
                      if (i2 <= 63) {
                        B[i2][1][1] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][1][1] - (SCALAR_VAL(2.0) * A[i2][1][1])) + A[i2 - 1][1][1])) + (SCALAR_VAL(0.125) * ((A[i2][2][1] - (SCALAR_VAL(2.0) * A[i2][1][1])) + A[i2][0][1]))) + (SCALAR_VAL(0.125) * ((A[i2][1][2] - (SCALAR_VAL(2.0) * A[i2][1][1])) + A[i2][1][0]))) + A[i2][1][1]);
                      }
                      for (int i4 = -((i2 + 54) / 59) + 3; i4 <= 128 * ii0 - 2 * i0 - i2 + 132; i4 += 1) {
                        B[i2][1][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][1][i4] - (SCALAR_VAL(2.0) * A[i2][1][i4])) + A[i2 - 1][1][i4])) + (SCALAR_VAL(0.125) * ((A[i2][2][i4] - (SCALAR_VAL(2.0) * A[i2][1][i4])) + A[i2][0][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][1][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][1][i4])) + A[i2][1][i4 - 1]))) + A[i2][1][i4]);
                      }
                    }
                    if (2 * i0 + i2 >= 128 * ii0 + 70) {
                      for (int i3 = max(2, 128 * ii0 - 2 * i0 - i2 + 73); i3 <= 128 * ii0 - 2 * i0 - i2 + 132; i3 += 1) {
                        for (int i4 = 1; i4 <= 128 * ii0 - 2 * i0 - i2 + 132; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    }
                  } else {
                    for (int i3 = 1; i3 <= 128 * ii0 - 2 * i0 + 67; i3 += 1) {
                      for (int i4 = 1; i4 <= min(128 * ii0 - 2 * i0 + 66, 128 * ii0 - 2 * i0 - i3 + 69); i4 += 1) {
                        B[66][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[67][i3][i4] - (SCALAR_VAL(2.0) * A[66][i3][i4])) + A[65][i3][i4])) + (SCALAR_VAL(0.125) * ((A[66][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[66][i3][i4])) + A[66][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[66][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[66][i3][i4])) + A[66][i3][i4 - 1]))) + A[66][i3][i4]);
                      }
                      if (128 * ii0 + 66 >= 2 * i0 + i3) {
                        for (int i4 = 128 * ii0 - 2 * i0 - i3 + 70; i4 <= 128 * ii0 - 2 * i0 + 66; i4 += 1) {
                          B[66][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[67][i3][i4] - (SCALAR_VAL(2.0) * A[66][i3][i4])) + A[65][i3][i4])) + (SCALAR_VAL(0.125) * ((A[66][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[66][i3][i4])) + A[66][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[66][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[66][i3][i4])) + A[66][i3][i4 - 1]))) + A[66][i3][i4]);
                        }
                      } else {
                        for (int i4 = 3; i4 <= 128 * ii0 - 2 * i0 + 66; i4 += 1) {
                          B[66][128 * ii0 - 2 * i0 + 67][i4] = ((((SCALAR_VAL(0.125) * ((A[67][128 * ii0 - 2 * i0 + 67][i4] - (SCALAR_VAL(2.0) * A[66][128 * ii0 - 2 * i0 + 67][i4])) + A[65][128 * ii0 - 2 * i0 + 67][i4])) + (SCALAR_VAL(0.125) * ((A[66][128 * ii0 - 2 * i0 + 68][i4] - (SCALAR_VAL(2.0) * A[66][128 * ii0 - 2 * i0 + 67][i4])) + A[66][128 * ii0 - 2 * i0 + 66][i4]))) + (SCALAR_VAL(0.125) * ((A[66][128 * ii0 - 2 * i0 + 67][i4 + 1] - (SCALAR_VAL(2.0) * A[66][128 * ii0 - 2 * i0 + 67][i4])) + A[66][128 * ii0 - 2 * i0 + 67][i4 - 1]))) + A[66][128 * ii0 - 2 * i0 + 67][i4]);
                        }
                      }
                      B[66][i3][128 * ii0 - 2 * i0 + 67] = ((((SCALAR_VAL(0.125) * ((A[67][i3][128 * ii0 - 2 * i0 + 67] - (SCALAR_VAL(2.0) * A[66][i3][128 * ii0 - 2 * i0 + 67])) + A[65][i3][128 * ii0 - 2 * i0 + 67])) + (SCALAR_VAL(0.125) * ((A[66][i3 + 1][128 * ii0 - 2 * i0 + 67] - (SCALAR_VAL(2.0) * A[66][i3][128 * ii0 - 2 * i0 + 67])) + A[66][i3 - 1][128 * ii0 - 2 * i0 + 67]))) + (SCALAR_VAL(0.125) * ((A[66][i3][128 * ii0 - 2 * i0 + 68] - (SCALAR_VAL(2.0) * A[66][i3][128 * ii0 - 2 * i0 + 67])) + A[66][i3][128 * ii0 - 2 * i0 + 66]))) + A[66][i3][128 * ii0 - 2 * i0 + 67]);
                    }
                  }
                }
                if (ii3 == 0 && ii4 == 0 && i0 == 64 * ii0 + 33) {
                  B[66][1][1] = ((((SCALAR_VAL(0.125) * ((A[67][1][1] - (SCALAR_VAL(2.0) * A[66][1][1])) + A[65][1][1])) + (SCALAR_VAL(0.125) * ((A[66][2][1] - (SCALAR_VAL(2.0) * A[66][1][1])) + A[66][0][1]))) + (SCALAR_VAL(0.125) * ((A[66][1][2] - (SCALAR_VAL(2.0) * A[66][1][1])) + A[66][1][0]))) + A[66][1][1]);
                } else if (ii3 == 1) {
                  for (int i2 = 43 * ii0 - i0 + (-ii0 + i0 - 1) / 3 + 90; i2 <= 66; i2 += 1) {
                    for (int i3 = 1; i3 <= 66; i3 += 1) {
                      if (i2 == 66 && 128 * ii0 + 132 >= 2 * i0 + i3) {
                        B[66][i3][128 * ii0 - 2 * i0 - i3 + 133] = ((((SCALAR_VAL(0.125) * ((A[67][i3][128 * ii0 - 2 * i0 - i3 + 133] - (SCALAR_VAL(2.0) * A[66][i3][128 * ii0 - 2 * i0 - i3 + 133])) + A[65][i3][128 * ii0 - 2 * i0 - i3 + 133])) + (SCALAR_VAL(0.125) * ((A[66][i3 + 1][128 * ii0 - 2 * i0 - i3 + 133] - (SCALAR_VAL(2.0) * A[66][i3][128 * ii0 - 2 * i0 - i3 + 133])) + A[66][i3 - 1][128 * ii0 - 2 * i0 - i3 + 133]))) + (SCALAR_VAL(0.125) * ((A[66][i3][128 * ii0 - 2 * i0 - i3 + 134] - (SCALAR_VAL(2.0) * A[66][i3][128 * ii0 - 2 * i0 - i3 + 133])) + A[66][i3][128 * ii0 - 2 * i0 - i3 + 132]))) + A[66][i3][128 * ii0 - 2 * i0 - i3 + 133]);
                      }
                      for (int i4 = max(max(1, 128 * ii0 - 2 * i0 - i2 - i3 + 198), 128 * ii0 - 2 * i0 - i3 + (i2 - 3) / 9 + 127); i4 <= 66; i4 += 1) {
                        B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                      }
                    }
                  }
                }
                if (ii4 == 0) {
                  for (int i2 = 128 * ii0 - 2 * i0 + 67; i2 <= 64; i2 += 1) {
                    for (int i3 = max(1, 128 * ii0 + 64 * ii3 - 2 * i0 - i2 + 68); i3 <= min(66, 128 * ii0 + 65 * ii3 - 2 * i0 - i2 + 131); i3 += 1) {
                      for (int i4 = 1; i4 <= min(128 * ii0 + 63 * ii3 - 2 * i0 - i2 + 131, 128 * ii0 - 2 * i0 - i2 - i3 + 196); i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                      if (ii3 == 1 && i3 == 66) {
                        A[i2][66][128 * ii0 - 2 * i0 - i2 + 131] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][66][128 * ii0 - 2 * i0 - i2 + 131] - (SCALAR_VAL(2.0) * B[i2][66][128 * ii0 - 2 * i0 - i2 + 131])) + B[i2 - 1][66][128 * ii0 - 2 * i0 - i2 + 131])) + (SCALAR_VAL(0.125) * ((B[i2][67][128 * ii0 - 2 * i0 - i2 + 131] - (SCALAR_VAL(2.0) * B[i2][66][128 * ii0 - 2 * i0 - i2 + 131])) + B[i2][65][128 * ii0 - 2 * i0 - i2 + 131]))) + (SCALAR_VAL(0.125) * ((B[i2][66][128 * ii0 - 2 * i0 - i2 + 132] - (SCALAR_VAL(2.0) * B[i2][66][128 * ii0 - 2 * i0 - i2 + 131])) + B[i2][66][128 * ii0 - 2 * i0 - i2 + 130]))) + B[i2][66][128 * ii0 - 2 * i0 - i2 + 131]);
                      }
                    }
                  }
                  if (ii3 == 1 && i0 == 64 * ii0 + 33) {
                    for (int i2 = 65; i2 <= 66; i2 += 1) {
                      for (int i3 = 1; i3 <= 64; i3 += 1) {
                        for (int i4 = 1; i4 <= -i3 + 65; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                  } else if (64 * ii0 + 32 >= i0) {
                    for (int i2 = 65; i2 <= 66; i2 += 1) {
                      for (int i3 = max(1, 128 * ii0 + 60 * ii3 - 2 * i0 + 7); i3 <= min(64, 128 * ii0 + 62 * ii3 - 2 * i0 + 66); i3 += 1) {
                        for (int i4 = 1; i4 <= min(128 * ii0 + 62 * ii3 - 2 * i0 + 66, 128 * ii0 - 2 * i0 - i3 + 131); i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                      if (ii3 == 1) {
                        for (int i3 = 65; i3 <= 66; i3 += 1) {
                          for (int i4 = 1; i4 <= 128 * ii0 - 2 * i0 + 66; i4 += 1) {
                            A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                          }
                        }
                      }
                    }
                  }
                } else {
                  for (int i2 = max(1, 128 * ii0 - 2 * i0 + 67); i2 <= 66; i2 += 1) {
                    for (int i3 = max(max(1, 128 * ii0 - 2 * i0 + 67), 128 * ii0 - 2 * i0 - i2 + 132); i3 <= 66; i3 += 1) {
                      for (int i4 = max(max(max(max(1, 128 * ii0 - 2 * i0 + 67), 128 * ii0 - 2 * i0 - i2 + 132), 128 * ii0 - 2 * i0 - i3 + 132), 128 * ii0 - 2 * i0 - i2 - i3 + 197); i4 <= 66; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                  }
                }
              }
              if (ii3 == 0 && ii4 == 0) {
                for (int i0 = 64 * ii0 + 34; i0 <= min(_PB_TSTEPS, 64 * ii0 + 64); i0 += 1) {
                  for (int i2 = 1; i2 <= 128 * ii0 - 2 * i0 + 131; i2 += 1) {
                    for (int i4 = 1; i4 <= 128 * ii0 - 2 * i0 - i2 + 132; i4 += 1) {
                      B[i2][1][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][1][i4] - (SCALAR_VAL(2.0) * A[i2][1][i4])) + A[i2 - 1][1][i4])) + (SCALAR_VAL(0.125) * ((A[i2][2][i4] - (SCALAR_VAL(2.0) * A[i2][1][i4])) + A[i2][0][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][1][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][1][i4])) + A[i2][1][i4 - 1]))) + A[i2][1][i4]);
                    }
                    if (i0 == 64 * ii0 + 34 && i2 <= 2) {
                      for (int i4 = 1; i4 <= -i2 + 64; i4 += 1) {
                        B[i2][2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][2][i4] - (SCALAR_VAL(2.0) * A[i2][2][i4])) + A[i2 - 1][2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][3][i4] - (SCALAR_VAL(2.0) * A[i2][2][i4])) + A[i2][1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][2][i4])) + A[i2][2][i4 - 1]))) + A[i2][2][i4]);
                      }
                    }
                    if (2 * i0 + i2 >= 128 * ii0 + 70) {
                      for (int i3 = max(2, 128 * ii0 - 2 * i0 - i2 + 73); i3 <= 128 * ii0 - 2 * i0 - i2 + 132; i3 += 1) {
                        for (int i4 = 1; i4 <= 128 * ii0 - 2 * i0 - i2 + 132; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    } else {
                      for (int i3 = 3; i3 <= 63; i3 += 1) {
                        for (int i4 = 1; i4 <= 63; i4 += 1) {
                          B[1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[2][i3][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[0][i3][i4])) + (SCALAR_VAL(0.125) * ((A[1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3][i4 - 1]))) + A[1][i3][i4]);
                        }
                      }
                    }
                  }
                  for (int i2 = 1; i2 <= 128 * ii0 - 2 * i0 + 130; i2 += 1) {
                    for (int i3 = 1; i3 <= 128 * ii0 - 2 * i0 - i2 + 131; i3 += 1) {
                      for (int i4 = 1; i4 <= 128 * ii0 - 2 * i0 - i2 + 131; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                  }
                }
              } else if (ii4 == 0) {
                for (int i0 = 64 * ii0 + 34; i0 <= min(_PB_TSTEPS, 64 * ii0 + 64); i0 += 1) {
                  for (int i2 = 1; i2 <= 66; i2 += 1) {
                    for (int i3 = max(1, 128 * ii0 - 2 * i0 - i2 + 133); i3 <= min(64, 128 * ii0 - 2 * i0 - i2 + 196); i3 += 1) {
                      if (i0 == 64 * ii0 + 34 && i3 == 1) {
                        B[i2][1][1] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][1][1] - (SCALAR_VAL(2.0) * A[i2][1][1])) + A[i2 - 1][1][1])) + (SCALAR_VAL(0.125) * ((A[i2][2][1] - (SCALAR_VAL(2.0) * A[i2][1][1])) + A[i2][0][1]))) + (SCALAR_VAL(0.125) * ((A[i2][1][2] - (SCALAR_VAL(2.0) * A[i2][1][1])) + A[i2][1][0]))) + A[i2][1][1]);
                      }
                      for (int i4 = max(1, 128 * ii0 - 2 * i0 - i3 + 71); i4 <= 128 * ii0 - 2 * i0 - i2 - i3 + 197; i4 += 1) {
                        B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                      }
                      if (i2 == 66) {
                        B[66][i3][128 * ii0 - 2 * i0 - i3 + 132] = ((((SCALAR_VAL(0.125) * ((A[67][i3][128 * ii0 - 2 * i0 - i3 + 132] - (SCALAR_VAL(2.0) * A[66][i3][128 * ii0 - 2 * i0 - i3 + 132])) + A[65][i3][128 * ii0 - 2 * i0 - i3 + 132])) + (SCALAR_VAL(0.125) * ((A[66][i3 + 1][128 * ii0 - 2 * i0 - i3 + 132] - (SCALAR_VAL(2.0) * A[66][i3][128 * ii0 - 2 * i0 - i3 + 132])) + A[66][i3 - 1][128 * ii0 - 2 * i0 - i3 + 132]))) + (SCALAR_VAL(0.125) * ((A[66][i3][128 * ii0 - 2 * i0 - i3 + 133] - (SCALAR_VAL(2.0) * A[66][i3][128 * ii0 - 2 * i0 - i3 + 132])) + A[66][i3][128 * ii0 - 2 * i0 - i3 + 131]))) + A[66][i3][128 * ii0 - 2 * i0 - i3 + 132]);
                      }
                    }
                    if (i2 == 66) {
                      B[66][128 * ii0 - 2 * i0 + 131][1] = ((((SCALAR_VAL(0.125) * ((A[67][128 * ii0 - 2 * i0 + 131][1] - (SCALAR_VAL(2.0) * A[66][128 * ii0 - 2 * i0 + 131][1])) + A[65][128 * ii0 - 2 * i0 + 131][1])) + (SCALAR_VAL(0.125) * ((A[66][128 * ii0 - 2 * i0 + 132][1] - (SCALAR_VAL(2.0) * A[66][128 * ii0 - 2 * i0 + 131][1])) + A[66][128 * ii0 - 2 * i0 + 130][1]))) + (SCALAR_VAL(0.125) * ((A[66][128 * ii0 - 2 * i0 + 131][2] - (SCALAR_VAL(2.0) * A[66][128 * ii0 - 2 * i0 + 131][1])) + A[66][128 * ii0 - 2 * i0 + 131][0]))) + A[66][128 * ii0 - 2 * i0 + 131][1]);
                    }
                    for (int i3 = 65; i3 <= 66; i3 += 1) {
                      for (int i4 = 1; i4 <= 128 * ii0 - 2 * i0 - i2 + 132; i4 += 1) {
                        B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                      }
                    }
                  }
                  for (int i2 = 1; i2 <= 66; i2 += 1) {
                    if (i2 >= 65) {
                      for (int i3 = 1; i3 <= 128 * ii0 - 2 * i0 + 130; i3 += 1) {
                        for (int i4 = 1; i4 <= 128 * ii0 - 2 * i0 - i3 + 131; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    } else {
                      for (int i3 = max(1, 128 * ii0 - 2 * i0 - i2 + 132); i3 <= min(65, 128 * ii0 - 2 * i0 - i2 + 195); i3 += 1) {
                        for (int i4 = 1; i4 <= 128 * ii0 - 2 * i0 - i2 - i3 + 196; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                      for (int i4 = 1; i4 <= 128 * ii0 - 2 * i0 - i2 + 131; i4 += 1) {
                        A[i2][66][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][66][i4] - (SCALAR_VAL(2.0) * B[i2][66][i4])) + B[i2 - 1][66][i4])) + (SCALAR_VAL(0.125) * ((B[i2][67][i4] - (SCALAR_VAL(2.0) * B[i2][66][i4])) + B[i2][65][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][66][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][66][i4])) + B[i2][66][i4 - 1]))) + B[i2][66][i4]);
                      }
                    }
                  }
                }
              }
            } else {
              for (int i0 = 64 * ii0 + 1; i0 <= min(_PB_TSTEPS, 64 * ii0 + 64); i0 += 1) {
                if (i0 >= 64 * ii0 + 2) {
                  if (i0 >= 64 * ii0 + 3) {
                    for (int i2 = max(1, 128 * ii0 - 2 * i0 + 68); i2 <= min(66, 128 * ii0 - 2 * i0 + 131); i2 += 1) {
                      for (int i3 = 1; i3 < i2 - 63; i3 += 1) {
                        for (int i4 = 128 * ii0 - 2 * i0 + 68; i4 <= 128 * ii0 - 2 * i0 + i2 - i3 + 4; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = 128 * ii0 - 2 * i0 + i2 - i3 + 5; i4 <= 66; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                      if (i2 <= 65) {
                        for (int i3 = max(1, i2 - 63); i3 <= 128 * ii0 - 2 * i0 - i2 + 132; i3 += 1) {
                          if (i3 <= 2) {
                            for (int i4 = 128 * ii0 - 2 * i0 - i2 + 133; i4 <= min(66, -i2 + 128); i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            if (i2 == 63) {
                              B[63][i3][66] = ((((SCALAR_VAL(0.125) * ((A[64][i3][66] - (SCALAR_VAL(2.0) * A[63][i3][66])) + A[62][i3][66])) + (SCALAR_VAL(0.125) * ((A[63][i3 + 1][66] - (SCALAR_VAL(2.0) * A[63][i3][66])) + A[63][i3 - 1][66]))) + (SCALAR_VAL(0.125) * ((A[63][i3][67] - (SCALAR_VAL(2.0) * A[63][i3][66])) + A[63][i3][65]))) + A[63][i3][66]);
                            } else if (i2 >= 64) {
                              for (int i4 = -i2 + 129; i4 <= 66; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                          } else if (i3 <= 63) {
                            for (int i4 = max(-256 * ii0 + 4 * i0 + 2 * i2 - 74, 256 * ii0 - 4 * i0 - 2 * i2 + 201); i4 <= 65; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            if (2 * i0 + i2 >= 128 * ii0 + 70) {
                              for (int i4 = 128 * ii0 - 2 * i0 - i2 + 133; i4 <= 65; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                            if (i2 <= 63) {
                              B[i2][i3][66] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][66] - (SCALAR_VAL(2.0) * A[i2][i3][66])) + A[i2 - 1][i3][66])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][66] - (SCALAR_VAL(2.0) * A[i2][i3][66])) + A[i2][i3 - 1][66]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][67] - (SCALAR_VAL(2.0) * A[i2][i3][66])) + A[i2][i3][65]))) + A[i2][i3][66]);
                            } else {
                              B[i2][i3][66] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][66] - (SCALAR_VAL(2.0) * A[i2][i3][66])) + A[i2 - 1][i3][66])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][66] - (SCALAR_VAL(2.0) * A[i2][i3][66])) + A[i2][i3 - 1][66]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][67] - (SCALAR_VAL(2.0) * A[i2][i3][66])) + A[i2][i3][65]))) + A[i2][i3][66]);
                            }
                          } else {
                            for (int i4 = 65; i4 <= 66; i4 += 1) {
                              B[128 * ii0 - 2 * i0 + 68][64][i4] = ((((SCALAR_VAL(0.125) * ((A[128 * ii0 - 2 * i0 + 69][64][i4] - (SCALAR_VAL(2.0) * A[128 * ii0 - 2 * i0 + 68][64][i4])) + A[128 * ii0 - 2 * i0 + 67][64][i4])) + (SCALAR_VAL(0.125) * ((A[128 * ii0 - 2 * i0 + 68][65][i4] - (SCALAR_VAL(2.0) * A[128 * ii0 - 2 * i0 + 68][64][i4])) + A[128 * ii0 - 2 * i0 + 68][63][i4]))) + (SCALAR_VAL(0.125) * ((A[128 * ii0 - 2 * i0 + 68][64][i4 + 1] - (SCALAR_VAL(2.0) * A[128 * ii0 - 2 * i0 + 68][64][i4])) + A[128 * ii0 - 2 * i0 + 68][64][i4 - 1]))) + A[128 * ii0 - 2 * i0 + 68][64][i4]);
                            }
                          }
                        }
                      } else {
                        for (int i3 = 3; i3 <= 128 * ii0 - 2 * i0 + 67; i3 += 1) {
                          for (int i4 = 128 * ii0 - 2 * i0 + 68; i4 <= 66; i4 += 1) {
                            B[66][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[67][i3][i4] - (SCALAR_VAL(2.0) * A[66][i3][i4])) + A[65][i3][i4])) + (SCALAR_VAL(0.125) * ((A[66][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[66][i3][i4])) + A[66][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[66][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[66][i3][i4])) + A[66][i3][i4 - 1]))) + A[66][i3][i4]);
                          }
                        }
                      }
                    }
                    if (i0 == 64 * ii0 + 33) {
                      for (int i4 = 2; i4 <= 66; i4 += 1) {
                        B[66][1][i4] = ((((SCALAR_VAL(0.125) * ((A[67][1][i4] - (SCALAR_VAL(2.0) * A[66][1][i4])) + A[65][1][i4])) + (SCALAR_VAL(0.125) * ((A[66][2][i4] - (SCALAR_VAL(2.0) * A[66][1][i4])) + A[66][0][i4]))) + (SCALAR_VAL(0.125) * ((A[66][1][i4 + 1] - (SCALAR_VAL(2.0) * A[66][1][i4])) + A[66][1][i4 - 1]))) + A[66][1][i4]);
                      }
                    }
                  } else {
                    for (int i2 = 64; i2 <= 66; i2 += 1) {
                      for (int i3 = 1; i3 <= -i2 + 128; i3 += 1) {
                        for (int i4 = max(64, -i2 + 129); i4 <= 66; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                      if (i2 == 66) {
                        for (int i4 = 64; i4 <= 66; i4 += 1) {
                          B[66][63][i4] = ((((SCALAR_VAL(0.125) * ((A[67][63][i4] - (SCALAR_VAL(2.0) * A[66][63][i4])) + A[65][63][i4])) + (SCALAR_VAL(0.125) * ((A[66][64][i4] - (SCALAR_VAL(2.0) * A[66][63][i4])) + A[66][62][i4]))) + (SCALAR_VAL(0.125) * ((A[66][63][i4 + 1] - (SCALAR_VAL(2.0) * A[66][63][i4])) + A[66][63][i4 - 1]))) + A[66][63][i4]);
                        }
                      }
                    }
                  }
                }
                for (int i2 = max(1, 128 * ii0 - 2 * i0 + 67); i2 <= min(64, 128 * ii0 - 2 * i0 + 130); i2 += 1) {
                  for (int i3 = 1; i3 <= 128 * ii0 - 2 * i0 - i2 + 131; i3 += 1) {
                    for (int i4 = 128 * ii0 - 2 * i0 - i2 + 132; i4 <= 66; i4 += 1) {
                      A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                    }
                  }
                }
                for (int i2 = 65; i2 <= 66; i2 += 1) {
                  for (int i3 = 1; i3 <= 128 * ii0 - 2 * i0 + 66; i3 += 1) {
                    for (int i4 = 128 * ii0 - 2 * i0 + 67; i4 <= 66; i4 += 1) {
                      A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                    }
                  }
                }
              }
            }
          }
        }
      }
    } else {
      for (int ii2 = 0; ii2 <= (_PB_N - 3) / 64; ii2 += 1) {
        for (int ii3 = 0; ii3 <= (_PB_N - 3) / 64; ii3 += 1) {
          if (_PB_N >= 64 * ii3 + 67) {
            for (int ii4 = 0; ii4 <= (_PB_N - 3) / 64; ii4 += 1) {
              for (int i0 = _PB_TSTEPS - 1; i0 <= _PB_TSTEPS; i0 += 1) {
                if (_PB_N >= 64 * ii4 + 67 && i0 == _PB_TSTEPS) {
                  if (_PB_N >= 64 * ii2 + 67) {
                    for (int i2 = max(1, 64 * ii2); i2 <= 64 * ii2 + 63; i2 += 1) {
                      if (i2 >= 64 * ii2 + 2) {
                        if (ii3 >= 1 && i2 == 64 * ii2 + 2) {
                          for (int i4 = 64 * ii4 + 1; i4 <= 64 * ii4 + 64; i4 += 1) {
                            B[64 * ii2 + 2][64 * ii3][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 3][64 * ii3][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][64 * ii3][i4])) + A[64 * ii2 + 1][64 * ii3][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][64 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][64 * ii3][i4])) + A[64 * ii2 + 2][64 * ii3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][64 * ii3][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][64 * ii3][i4])) + A[64 * ii2 + 2][64 * ii3][i4 - 1]))) + A[64 * ii2 + 2][64 * ii3][i4]);
                          }
                        }
                        for (int i3 = max(max(64 * ii3, 63 * ii3 + 1), 64 * ii2 + 64 * ii3 - i2 + 3); i3 <= min(64 * ii3 + 63, -64 * ii2 + 64 * ii3 + i2 + 60); i3 += 1) {
                          if (i2 + i3 >= 64 * ii2 + 64 * ii3 + 65) {
                            for (int i4 = max(1, 64 * ii4); i4 <= 64 * ii4 + 63; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          } else if (i2 >= 64 * ii2 + 3) {
                            for (int i4 = max(max(1, 64 * ii4), 128 * ii3 + 64 * ii4 - 2 * i3 + 1); i4 <= 128 * ii3 + 64 * ii4 - 2 * i3 + 64; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            for (int i4 = max(max(1, 64 * ii4), 128 * ii3 + 64 * ii4 - 2 * i3 + 65); i4 <= 64 * ii4 + 63; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          } else {
                            for (int i4 = max(1, 64 * ii4); i4 <= 64 * ii4 + 63; i4 += 1) {
                              B[64 * ii2 + 2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 3][i3][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][i3][i4])) + A[64 * ii2 + 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][i3][i4])) + A[64 * ii2 + 2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][i3][i4])) + A[64 * ii2 + 2][i3][i4 - 1]))) + A[64 * ii2 + 2][i3][i4]);
                            }
                          }
                        }
                        if (i2 == 64 * ii2 + 2) {
                          for (int i4 = max(1, 64 * ii4); i4 <= 64 * ii4 + 63; i4 += 1) {
                            B[64 * ii2 + 2][64 * ii3 + 63][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 3][64 * ii3 + 63][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][64 * ii3 + 63][i4])) + A[64 * ii2 + 1][64 * ii3 + 63][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][64 * ii3 + 64][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][64 * ii3 + 63][i4])) + A[64 * ii2 + 2][64 * ii3 + 62][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][64 * ii3 + 63][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][64 * ii3 + 63][i4])) + A[64 * ii2 + 2][64 * ii3 + 63][i4 - 1]))) + A[64 * ii2 + 2][64 * ii3 + 63][i4]);
                          }
                        }
                      } else {
                        for (int i3 = max(1, 64 * ii2 + 64 * ii3 - i2 + 1); i3 <= 64 * ii2 + 64 * ii3 - i2 + 64; i3 += 1) {
                          for (int i4 = max(max(1, 64 * ii2 + 64 * ii4 - i2 + 1), 128 * ii3 + 64 * ii4 - 2 * i3 + 1); i4 <= 64 * ii2 + 64 * ii4 - i2 + 64; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (i2 == 64 * ii2 + 1 && i3 == 64 * ii3) {
                            B[64 * ii2 + 1][64 * ii3][64 * ii4 + 64] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][64 * ii3][64 * ii4 + 64] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][64 * ii3][64 * ii4 + 64])) + A[64 * ii2][64 * ii3][64 * ii4 + 64])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][64 * ii3 + 1][64 * ii4 + 64] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][64 * ii3][64 * ii4 + 64])) + A[64 * ii2 + 1][64 * ii3 - 1][64 * ii4 + 64]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][64 * ii3][64 * ii4 + 65] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][64 * ii3][64 * ii4 + 64])) + A[64 * ii2 + 1][64 * ii3][64 * ii4 + 63]))) + A[64 * ii2 + 1][64 * ii3][64 * ii4 + 64]);
                          }
                        }
                      }
                    }
                  } else {
                    for (int i2 = 64 * ii2; i2 <= min(_PB_N - 2, 64 * ii2 + 2); i2 += 1) {
                      if (64 * ii2 + 1 >= i2) {
                        for (int i3 = max(1, 64 * ii2 + 64 * ii3 - i2 + 1); i3 <= 64 * ii2 + 64 * ii3 - i2 + 64; i3 += 1) {
                          for (int i4 = max(max(1, 64 * ii2 + 64 * ii4 - i2 + 1), 128 * ii3 + 64 * ii4 - 2 * i3 + 1); i4 <= 64 * ii2 + 64 * ii4 - i2 + 64; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (i2 == 64 * ii2 + 1 && i3 == 64 * ii3) {
                            B[64 * ii2 + 1][64 * ii3][64 * ii4 + 64] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][64 * ii3][64 * ii4 + 64] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][64 * ii3][64 * ii4 + 64])) + A[64 * ii2][64 * ii3][64 * ii4 + 64])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][64 * ii3 + 1][64 * ii4 + 64] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][64 * ii3][64 * ii4 + 64])) + A[64 * ii2 + 1][64 * ii3 - 1][64 * ii4 + 64]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][64 * ii3][64 * ii4 + 65] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][64 * ii3][64 * ii4 + 64])) + A[64 * ii2 + 1][64 * ii3][64 * ii4 + 63]))) + A[64 * ii2 + 1][64 * ii3][64 * ii4 + 64]);
                          }
                        }
                      } else {
                        if (ii3 >= 1) {
                          for (int i4 = 64 * ii4 + 1; i4 <= 64 * ii4 + 64; i4 += 1) {
                            B[64 * ii2 + 2][64 * ii3][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 3][64 * ii3][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][64 * ii3][i4])) + A[64 * ii2 + 1][64 * ii3][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][64 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][64 * ii3][i4])) + A[64 * ii2 + 2][64 * ii3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][64 * ii3][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][64 * ii3][i4])) + A[64 * ii2 + 2][64 * ii3][i4 - 1]))) + A[64 * ii2 + 2][64 * ii3][i4]);
                          }
                        }
                        for (int i3 = 64 * ii3 + 1; i3 <= 64 * ii3 + 63; i3 += 1) {
                          for (int i4 = max(1, 64 * ii4); i4 <= 64 * ii4 + 63; i4 += 1) {
                            B[64 * ii2 + 2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 3][i3][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][i3][i4])) + A[64 * ii2 + 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][i3][i4])) + A[64 * ii2 + 2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][i3][i4])) + A[64 * ii2 + 2][i3][i4 - 1]))) + A[64 * ii2 + 2][i3][i4]);
                          }
                        }
                      }
                    }
                    for (int i2 = 64 * ii2 + 3; i2 < _PB_N - 1; i2 += 1) {
                      for (int i3 = max(1, 64 * ii3); i3 <= 64 * ii3 + 63; i3 += 1) {
                        for (int i4 = max(max(1, 64 * ii4), 128 * ii3 + 64 * ii4 - 2 * i3 + 1); i4 <= 128 * ii3 + 64 * ii4 - 2 * i3 + 64; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = max(max(1, 64 * ii4), 128 * ii3 + 64 * ii4 - 2 * i3 + 65); i4 <= 64 * ii4 + 63; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    }
                  }
                }
                if (_PB_N >= 64 * ii4 + 67) {
                  if (i0 == _PB_TSTEPS) {
                    for (int i2 = max(1, 64 * ii2 - 1); i2 <= 64 * ii2; i2 += 1) {
                      for (int i3 = max(1, 64 * ii2 + 64 * ii3 - i2); i3 <= 64 * ii2 + 64 * ii3 - i2 + 63; i3 += 1) {
                        for (int i4 = max(max(1, 64 * ii2 + 64 * ii4 - i2), 64 * ii2 + 64 * ii3 + 64 * ii4 - i2 - i3 + 1); i4 <= 64 * ii2 + 64 * ii3 + 64 * ii4 - i2 - i3 + 64; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                        for (int i4 = 64 * ii2 + 64 * ii3 + 64 * ii4 - i2 - i3 + 65; i4 <= 64 * ii2 + 64 * ii4 - i2 + 63; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                  }
                  if (_PB_N >= 64 * ii2 + 67) {
                    for (int i2 = 64 * ii2 + 1; i2 <= 2 * _PB_TSTEPS + 64 * ii2 - 2 * i0 + 62; i2 += 1) {
                      if (i0 == _PB_TSTEPS) {
                        for (int i3 = max(1, 64 * ii3 - 1); i3 <= 64 * ii3; i3 += 1) {
                          for (int i4 = max(1, 64 * ii3 + 64 * ii4 - i3); i4 <= 64 * ii3 + 64 * ii4 - i3 + 63; i4 += 1) {
                            A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                          }
                        }
                      }
                      for (int i3 = 64 * ii3 + 1; i3 <= 2 * _PB_TSTEPS + 64 * ii3 - 2 * i0 + 62; i3 += 1) {
                        for (int i4 = max(1, 2 * _PB_TSTEPS + 64 * ii4 - 2 * i0 - 1); i4 <= 2 * _PB_TSTEPS + 64 * ii4 - 2 * i0 + 62; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                  } else {
                    for (int i2 = 64 * ii2 + 1; i2 < _PB_N - 1; i2 += 1) {
                      if (i0 == _PB_TSTEPS) {
                        for (int i3 = max(1, 64 * ii3 - 1); i3 <= 64 * ii3; i3 += 1) {
                          for (int i4 = max(1, 64 * ii3 + 64 * ii4 - i3); i4 <= 64 * ii3 + 64 * ii4 - i3 + 63; i4 += 1) {
                            A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                          }
                        }
                      }
                      for (int i3 = 64 * ii3 + 1; i3 <= 2 * _PB_TSTEPS + 64 * ii3 - 2 * i0 + 62; i3 += 1) {
                        if (i0 == _PB_TSTEPS) {
                          for (int i4 = max(1, 64 * ii4 - 1); i4 <= 64 * ii4; i4 += 1) {
                            A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                          }
                        }
                        for (int i4 = 64 * ii4 + 1; i4 <= 2 * _PB_TSTEPS + 64 * ii4 - 2 * i0 + 62; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                  }
                } else {
                  if (_PB_N >= 64 * ii2 + 67 && i0 == _PB_TSTEPS) {
                    for (int i2 = max(1, 64 * ii2); i2 <= 64 * ii2 + 63; i2 += 1) {
                      for (int i3 = max(max(1, 64 * ii3), 64 * ii2 + 64 * ii3 - i2 + 1); i3 <= 64 * ii2 + 64 * ii3 - i2 + 64; i3 += 1) {
                        for (int i4 = max(max(64 * ii4, 128 * ii2 + 64 * ii4 - 2 * i2 + 1), 32 * ii3 + 64 * ii4 - i3 + i3 / 2 + 1); i4 < _PB_N - 1; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                      for (int i3 = 64 * ii2 + 64 * ii3 - i2 + 65; i3 <= 64 * ii3 + 63; i3 += 1) {
                        for (int i4 = 64 * ii4; i4 < _PB_N - 1; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    }
                  }
                  if (_PB_N >= 64 * ii2 + 67) {
                    if (i0 == _PB_TSTEPS) {
                      for (int i2 = max(1, 64 * ii2 - 1); i2 <= 64 * ii2; i2 += 1) {
                        for (int i3 = max(1, 64 * ii2 + 64 * ii3 - i2); i3 <= 64 * ii2 + 64 * ii3 - i2 + 63; i3 += 1) {
                          for (int i4 = max(64 * ii2 + 64 * ii4 - i2, 64 * ii2 + 64 * ii3 + 64 * ii4 - i2 - i3 + 1); i4 < _PB_N - 1; i4 += 1) {
                            A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                          }
                        }
                      }
                    }
                    for (int i2 = 64 * ii2 + 1; i2 <= 2 * _PB_TSTEPS + 64 * ii2 - 2 * i0 + 62; i2 += 1) {
                      if (i0 == _PB_TSTEPS) {
                        for (int i3 = max(1, 64 * ii3 - 1); i3 <= 64 * ii3; i3 += 1) {
                          for (int i4 = 64 * ii3 + 64 * ii4 - i3; i4 < _PB_N - 1; i4 += 1) {
                            A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                          }
                        }
                      }
                      for (int i3 = 64 * ii3 + 1; i3 <= 2 * _PB_TSTEPS + 64 * ii3 - 2 * i0 + 62; i3 += 1) {
                        for (int i4 = 2 * _PB_TSTEPS + 64 * ii4 - 2 * i0 - 1; i4 < _PB_N - 1; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                  } else {
                    if (i0 == _PB_TSTEPS) {
                      for (int i2 = 64 * ii2; i2 < _PB_N - 1; i2 += 1) {
                        for (int i3 = max(max(1, 64 * ii3), 64 * ii2 + 64 * ii3 - i2 + 1); i3 <= 64 * ii2 + 64 * ii3 - i2 + 64; i3 += 1) {
                          for (int i4 = max(max(64 * ii2, 128 * ii2 - i2 + 1), 64 * ii2 + 32 * ii3 - i3 + i3 / 2 + 1); i4 < _PB_N - 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                        for (int i3 = 64 * ii2 + 64 * ii3 - i2 + 65; i3 <= 64 * ii3 + 63; i3 += 1) {
                          for (int i4 = 64 * ii2; i4 < _PB_N - 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                    }
                    if (i0 == _PB_TSTEPS) {
                      for (int i2 = 64 * ii2 - 1; i2 <= 64 * ii2; i2 += 1) {
                        for (int i3 = max(1, 64 * ii2 + 64 * ii3 - i2); i3 <= 64 * ii2 + 64 * ii3 - i2 + 63; i3 += 1) {
                          for (int i4 = max(128 * ii2 - i2, 128 * ii2 + 64 * ii3 - i2 - i3 + 1); i4 < _PB_N - 1; i4 += 1) {
                            A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                          }
                        }
                      }
                    }
                    for (int i2 = 64 * ii2 + 1; i2 < _PB_N - 1; i2 += 1) {
                      if (i0 == _PB_TSTEPS) {
                        for (int i3 = max(1, 64 * ii3 - 1); i3 <= 64 * ii3; i3 += 1) {
                          for (int i4 = 64 * ii2 + 64 * ii3 - i3; i4 < _PB_N - 1; i4 += 1) {
                            A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                          }
                        }
                      }
                      for (int i3 = 64 * ii3 + 1; i3 <= 2 * _PB_TSTEPS + 64 * ii3 - 2 * i0 + 62; i3 += 1) {
                        if (i0 == _PB_TSTEPS) {
                          for (int i4 = 64 * ii2 - 1; i4 <= 64 * ii2; i4 += 1) {
                            A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                          }
                        }
                        for (int i4 = 64 * ii2 + 1; i4 < _PB_N - 1; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                  }
                }
              }
            }
          } else if (_PB_N >= 64 * ii2 + 67) {
            for (int ii4 = 0; ii4 <= ii3; ii4 += 1) {
              for (int i0 = _PB_TSTEPS - 1; i0 <= _PB_TSTEPS; i0 += 1) {
                if (i0 == _PB_TSTEPS) {
                  for (int i2 = max(1, 64 * ii2); i2 <= 64 * ii2 + 63; i2 += 1) {
                    if (_PB_N >= 64 * ii4 + 67 && i2 >= 64 * ii2 + 2) {
                      if (i2 == 64 * ii2 + 2) {
                        for (int i4 = 64 * ii4 + 1; i4 <= 64 * ii4 + 64; i4 += 1) {
                          B[64 * ii2 + 2][64 * ii3][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 3][64 * ii3][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][64 * ii3][i4])) + A[64 * ii2 + 1][64 * ii3][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][64 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][64 * ii3][i4])) + A[64 * ii2 + 2][64 * ii3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][64 * ii3][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][64 * ii3][i4])) + A[64 * ii2 + 2][64 * ii3][i4 - 1]))) + A[64 * ii2 + 2][64 * ii3][i4]);
                        }
                      }
                      for (int i3 = max(64 * ii3, 64 * ii2 + 64 * ii3 - i2 + 3); i3 < _PB_N - 1; i3 += 1) {
                        for (int i4 = max(max(1, 64 * ii4), 32 * ii3 + 64 * ii4 - i3 + i3 / 2 + 1); i4 <= 64 * ii4 + 63; i4 += 1) {
                          if (64 * ii4 + i2 + 60 >= 64 * ii2 + i4 || 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                        if (i3 == 64 * ii3) {
                          B[i2][64 * ii3][64 * ii4 + 64] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii3][64 * ii4 + 64] - (SCALAR_VAL(2.0) * A[i2][64 * ii3][64 * ii4 + 64])) + A[i2 - 1][64 * ii3][64 * ii4 + 64])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3 + 1][64 * ii4 + 64] - (SCALAR_VAL(2.0) * A[i2][64 * ii3][64 * ii4 + 64])) + A[i2][64 * ii3 - 1][64 * ii4 + 64]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii3][64 * ii4 + 65] - (SCALAR_VAL(2.0) * A[i2][64 * ii3][64 * ii4 + 64])) + A[i2][64 * ii3][64 * ii4 + 63]))) + A[i2][64 * ii3][64 * ii4 + 64]);
                        }
                      }
                    } else if (_PB_N >= 64 * ii4 + 67 && 64 * ii2 + 1 >= i2) {
                      for (int i3 = 64 * ii2 + 64 * ii3 - i2 + 1; i3 < _PB_N - 1; i3 += 1) {
                        for (int i4 = max(max(1, 64 * ii2 + 64 * ii4 - i2 + 1), 32 * ii3 + 64 * ii4 - i3 + i3 / 2 + 1); i4 <= 64 * ii2 + 64 * ii4 - i2 + 64; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        if (i2 == 64 * ii2 + 1 && i3 == 64 * ii3) {
                          B[64 * ii2 + 1][64 * ii3][64 * ii4 + 64] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][64 * ii3][64 * ii4 + 64] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][64 * ii3][64 * ii4 + 64])) + A[64 * ii2][64 * ii3][64 * ii4 + 64])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][64 * ii3 + 1][64 * ii4 + 64] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][64 * ii3][64 * ii4 + 64])) + A[64 * ii2 + 1][64 * ii3 - 1][64 * ii4 + 64]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 1][64 * ii3][64 * ii4 + 65] - (SCALAR_VAL(2.0) * A[64 * ii2 + 1][64 * ii3][64 * ii4 + 64])) + A[64 * ii2 + 1][64 * ii3][64 * ii4 + 63]))) + A[64 * ii2 + 1][64 * ii3][64 * ii4 + 64]);
                        }
                      }
                    } else {
                      for (int i3 = max(64 * ii3, 64 * ii2 + 64 * ii3 - i2 + 1); i3 < _PB_N - 1; i3 += 1) {
                        for (int i4 = max(max(64 * ii3, 64 * ii2 + 64 * ii3 - i2 + 1), 192 * ii3 - 2 * i3 + 1); i4 < _PB_N - 1; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    }
                  }
                }
                if (_PB_N >= 64 * ii4 + 67 && i0 == _PB_TSTEPS) {
                  for (int i2 = max(1, 64 * ii2 - 1); i2 <= 64 * ii2; i2 += 1) {
                    for (int i3 = 64 * ii2 + 64 * ii3 - i2; i3 <= 64 * ii3 + 1; i3 += 1) {
                      for (int i4 = max(1, 64 * ii2 + 64 * ii3 + 64 * ii4 - i2 - i3 + 1); i4 <= 64 * ii2 + 64 * ii3 + 64 * ii4 - i2 - i3 + 64; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                    for (int i3 = 64 * ii3 + 2; i3 < _PB_N - 1; i3 += 1) {
                      for (int i4 = max(1, 64 * ii2 + 64 * ii4 - i2); i4 <= 64 * ii2 + 64 * ii4 - i2 + 63; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                  }
                } else if (ii4 == ii3 && i0 == _PB_TSTEPS) {
                  for (int i2 = max(1, 64 * ii2 - 1); i2 <= 64 * ii2; i2 += 1) {
                    for (int i3 = 64 * ii2 + 64 * ii3 - i2; i3 < _PB_N - 1; i3 += 1) {
                      for (int i4 = max(64 * ii2 + 64 * ii3 - i2, 64 * ii2 + 128 * ii3 - i2 - i3 + 1); i4 < _PB_N - 1; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                  }
                }
                for (int i2 = 64 * ii2 + 1; i2 <= 2 * _PB_TSTEPS + 64 * ii2 - 2 * i0 + 62; i2 += 1) {
                  if (_PB_N >= 64 * ii4 + 67 && i0 == _PB_TSTEPS) {
                    for (int i3 = 64 * ii3 - 1; i3 <= 64 * ii3; i3 += 1) {
                      for (int i4 = max(1, 64 * ii3 + 64 * ii4 - i3); i4 <= 64 * ii3 + 64 * ii4 - i3 + 63; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                  } else if (ii4 == ii3 && i0 == _PB_TSTEPS) {
                    for (int i3 = 64 * ii3 - 1; i3 <= 64 * ii3; i3 += 1) {
                      for (int i4 = 128 * ii3 - i3; i4 < _PB_N - 1; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                  }
                  for (int i3 = 64 * ii3 + 1; i3 < _PB_N - 1; i3 += 1) {
                    if (_PB_N >= 64 * ii4 + 67) {
                      for (int i4 = max(1, 2 * _PB_TSTEPS + 64 * ii4 - 2 * i0 - 1); i4 <= 2 * _PB_TSTEPS + 64 * ii4 - 2 * i0 + 62; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    } else {
                      if (i0 == _PB_TSTEPS) {
                        for (int i4 = 64 * ii3 - 1; i4 <= 64 * ii3; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                      for (int i4 = 64 * ii3 + 1; i4 < _PB_N - 1; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                  }
                }
              }
            }
          } else {
            for (int ii4 = 0; ii4 <= ii2; ii4 += 1) {
              if (_PB_N >= 64 * ii4 + 67) {
                for (int i0 = _PB_TSTEPS - 1; i0 <= _PB_TSTEPS; i0 += 1) {
                  if (i0 == _PB_TSTEPS) {
                    for (int i2 = 64 * ii2; i2 < _PB_N - 1; i2 += 1) {
                      for (int i3 = max(64 * ii2, 128 * ii2 - i2 + 1); i3 < _PB_N - 1; i3 += 1) {
                        if (i2 >= 64 * ii2 + 3) {
                          for (int i4 = max(max(1, 64 * ii4), 128 * ii2 + 64 * ii4 - 2 * i3 + 1); i4 <= 128 * ii2 + 64 * ii4 - 2 * i3 + 64; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          for (int i4 = max(max(1, 64 * ii4), 128 * ii2 + 64 * ii4 - 2 * i3 + 65); i4 <= 64 * ii4 + 63; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        } else {
                          for (int i4 = max(max(max(1, 64 * ii4), 64 * ii2 + 64 * ii4 - i2 + 1), 128 * ii2 + 64 * ii4 - 2 * i3 + 1); i4 <= 64 * ii2 + 64 * ii4 - i2 + 64; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (i2 == 64 * ii2 + 2 && i3 >= 64 * ii2 + 1) {
                            B[64 * ii2 + 2][i3][64 * ii4 + 63] = ((((SCALAR_VAL(0.125) * ((A[64 * ii2 + 3][i3][64 * ii4 + 63] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][i3][64 * ii4 + 63])) + A[64 * ii2 + 1][i3][64 * ii4 + 63])) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][i3 + 1][64 * ii4 + 63] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][i3][64 * ii4 + 63])) + A[64 * ii2 + 2][i3 - 1][64 * ii4 + 63]))) + (SCALAR_VAL(0.125) * ((A[64 * ii2 + 2][i3][64 * ii4 + 64] - (SCALAR_VAL(2.0) * A[64 * ii2 + 2][i3][64 * ii4 + 63])) + A[64 * ii2 + 2][i3][64 * ii4 + 62]))) + A[64 * ii2 + 2][i3][64 * ii4 + 63]);
                          } else if (i3 == 64 * ii2) {
                            for (int i4 = 64 * ii2 + 64 * ii4 - i2 + 65; i4 <= 64 * ii4 + 64; i4 += 1) {
                              B[i2][64 * ii2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii2][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii2][i4])) + A[i2 - 1][64 * ii2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii2 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii2][i4])) + A[i2][64 * ii2 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii2][i4])) + A[i2][64 * ii2][i4 - 1]))) + A[i2][64 * ii2][i4]);
                            }
                          }
                        }
                      }
                    }
                  }
                  if (i0 == _PB_TSTEPS) {
                    for (int i2 = 64 * ii2 - 1; i2 <= 64 * ii2; i2 += 1) {
                      for (int i3 = 128 * ii2 - i2; i3 < _PB_N - 1; i3 += 1) {
                        for (int i4 = max(max(1, 64 * ii2 + 64 * ii4 - i2), 128 * ii2 + 64 * ii4 - i2 - i3 + 1); i4 <= 128 * ii2 + 64 * ii4 - i2 - i3 + 64; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                        for (int i4 = 128 * ii2 + 64 * ii4 - i2 - i3 + 65; i4 <= 64 * ii2 + 64 * ii4 - i2 + 63; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                  }
                  for (int i2 = 64 * ii2 + 1; i2 < _PB_N - 1; i2 += 1) {
                    if (i0 == _PB_TSTEPS) {
                      for (int i3 = 64 * ii2 - 1; i3 <= 64 * ii2; i3 += 1) {
                        for (int i4 = max(1, 64 * ii2 + 64 * ii4 - i3); i4 <= 64 * ii2 + 64 * ii4 - i3 + 63; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                    for (int i3 = 64 * ii2 + 1; i3 < _PB_N - 1; i3 += 1) {
                      for (int i4 = max(1, 2 * _PB_TSTEPS + 64 * ii4 - 2 * i0 - 1); i4 <= 2 * _PB_TSTEPS + 64 * ii4 - 2 * i0 + 62; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                  }
                }
              } else {
                for (int i0 = _PB_TSTEPS - 1; i0 <= _PB_TSTEPS; i0 += 1) {
                  if (i0 == _PB_TSTEPS) {
                    for (int i2 = 64 * ii2; i2 < _PB_N - 1; i2 += 1) {
                      for (int i3 = max(64 * ii2, 128 * ii2 - i2 + 1); i3 < _PB_N - 1; i3 += 1) {
                        for (int i4 = max(max(64 * ii2, 128 * ii2 - i2 + 1), 192 * ii2 - 2 * i3 + 1); i4 < _PB_N - 1; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    }
                  }
                  if (i0 == _PB_TSTEPS) {
                    for (int i2 = 64 * ii2 - 1; i2 <= 64 * ii2; i2 += 1) {
                      for (int i3 = 128 * ii2 - i2; i3 < _PB_N - 1; i3 += 1) {
                        for (int i4 = max(128 * ii2 - i2, 192 * ii2 - i2 - i3 + 1); i4 < _PB_N - 1; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                  }
                  for (int i2 = 64 * ii2 + 1; i2 < _PB_N - 1; i2 += 1) {
                    if (i0 == _PB_TSTEPS) {
                      for (int i3 = 64 * ii2 - 1; i3 <= 64 * ii2; i3 += 1) {
                        for (int i4 = 128 * ii2 - i3; i4 < _PB_N - 1; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                    for (int i3 = 64 * ii2 + 1; i3 < _PB_N - 1; i3 += 1) {
                      if (i0 == _PB_TSTEPS) {
                        for (int i4 = 64 * ii2 - 1; i4 <= 64 * ii2; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                      for (int i4 = 64 * ii2 + 1; i4 < _PB_N - 1; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
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
  }
  if (_PB_N >= 67 && 64 * ii0 + 1 == _PB_TSTEPS) {
    for (int ii2 = 0; ii2 <= (_PB_N - 3) / 64; ii2 += 1) {
      if (_PB_N >= 64 * ii2 + 67) {
        for (int ii3 = 0; ii3 <= (_PB_N - 3) / 64; ii3 += 1) {
          for (int ii4 = 0; ii4 < (_PB_N - 3) / 64; ii4 += 1) {
            for (int i2 = 64 * ii2 + 1; i2 <= 64 * ii2 + 64; i2 += 1) {
              for (int i3 = 64 * ii3 + 1; i3 <= min(_PB_N - 2, 64 * ii3 + 64); i3 += 1) {
                for (int i4 = 64 * ii4 + 1; i4 <= 64 * ii4 + 64; i4 += 1) {
                  A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                }
              }
            }
          }
          if (_PB_N >= 64 * ii3 + 67) {
            for (int i2 = 64 * ii2 + 1; i2 <= 64 * ii2 + 64; i2 += 1) {
              for (int i3 = 64 * ii3 + 1; i3 <= 64 * ii3 + 64; i3 += 1) {
                for (int i4 = -((_PB_N - 3) % 64) + _PB_N - 2; i4 < _PB_N - 1; i4 += 1) {
                  A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                }
              }
            }
          } else {
            for (int i2 = 64 * ii2 + 1; i2 <= 64 * ii2 + 64; i2 += 1) {
              for (int i3 = 64 * ii3 + 1; i3 < _PB_N - 1; i3 += 1) {
                for (int i4 = 64 * ii3 + 1; i4 < _PB_N - 1; i4 += 1) {
                  A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                }
              }
            }
          }
        }
      } else {
        for (int ii3 = 0; ii3 <= ii2; ii3 += 1) {
          if (_PB_N >= 64 * ii3 + 67 || 1) {
            for (int ii4 = 0; ii4 < ii2; ii4 += 1) {
              for (int i2 = 64 * ii2 + 1; i2 < _PB_N - 1; i2 += 1) {
                for (int i3 = 64 * ii3 + 1; i3 <= min(_PB_N - 2, 64 * ii3 + 64); i3 += 1) {
                  for (int i4 = 64 * ii4 + 1; i4 <= 64 * ii4 + 64; i4 += 1) {
                    A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                  }
                }
              }
            }
            if (ii3 == ii2) {
              for (int i2 = 64 * ii2 + 1; i2 < _PB_N - 1; i2 += 1) {
                for (int i3 = 64 * ii2 + 1; i3 < _PB_N - 1; i3 += 1) {
                  for (int i4 = 64 * ii2 + 1; i4 < _PB_N - 1; i4 += 1) {
                    A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                  }
                }
              }
            } else {
              for (int i2 = 64 * ii2 + 1; i2 < _PB_N - 1; i2 += 1) {
                for (int i3 = 64 * ii3 + 1; i3 <= 64 * ii3 + 64; i3 += 1) {
                  for (int i4 = 64 * ii2 + 1; i4 < _PB_N - 1; i4 += 1) {
                    A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                  }
                }
              }
            }
          }
        }
      }
    }
  } else if (_PB_N <= 66) {
    for (int i0 = 64 * ii0 + 1; i0 <= min(_PB_TSTEPS, 64 * ii0 + 64); i0 += 1) {
      if (i0 >= 64 * ii0 + 2) {
        for (int i2 = 1; i2 < _PB_N - 1; i2 += 1) {
          for (int i3 = 1; i3 < _PB_N - 1; i3 += 1) {
            for (int i4 = 1; i4 < _PB_N - 1; i4 += 1) {
              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
            }
          }
        }
      }
      for (int i2 = 1; i2 < _PB_N - 1; i2 += 1) {
        for (int i3 = 1; i3 < _PB_N - 1; i3 += 1) {
          for (int i4 = 1; i4 < _PB_N - 1; i4 += 1) {
            A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
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
  int tsteps = TSTEPS;

  /* Variable declaration/allocation. */
  POLYBENCH_3D_ARRAY_DECL(A, DATA_TYPE, N, N, N, n, n, n);
  POLYBENCH_3D_ARRAY_DECL(B, DATA_TYPE, N, N, N, n, n, n);


  /* Initialize array(s). */
  init_array (n, POLYBENCH_ARRAY(A), POLYBENCH_ARRAY(B));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_heat_3d (tsteps, n, POLYBENCH_ARRAY(A), POLYBENCH_ARRAY(B));

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
