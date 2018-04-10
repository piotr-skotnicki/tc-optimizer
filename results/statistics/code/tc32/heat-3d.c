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
/* ./tc ../examples/polybench/heat-3d.scop.c --correction-tiling --lex-scheduling --serial-codegen --debug -b 32 */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#pragma scop
for (int ii0 = 0; ii0 <= floord(_PB_TSTEPS - 1, 32); ii0 += 1) {
  for (int ii1 = 0; ii1 <= min(min(1, _PB_TSTEPS - 32 * ii0 - 1), floord(_PB_N - 3, 32)); ii1 += 1) {
    if (ii1 == 0) {
      for (int ii2 = 0; ii2 <= (_PB_N - 3) / 32; ii2 += 1) {
        for (int ii3 = 0; ii3 <= (_PB_N - 3) / 32; ii3 += 1) {
          for (int ii4 = 0; ii4 <= (_PB_N - 3) / 32; ii4 += 1) {
            for (int i2 = 32 * ii2 + 1; i2 <= min(_PB_N - 2, 32 * ii2 + 32); i2 += 1) {
              for (int i3 = 32 * ii3 + 1; i3 <= min(_PB_N - 2, 32 * ii3 + 32); i3 += 1) {
                for (int i4 = 32 * ii4 + 1; i4 <= min(_PB_N - 2, 32 * ii4 + 32); i4 += 1) {
                  B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                }
              }
            }
          }
        }
      }
    } else if (_PB_TSTEPS >= 32 * ii0 + 3) {
      for (int ii2 = 0; ii2 <= min((_PB_N - 5) / 16 - 1, (_PB_N - 3) / 32); ii2 += 1) {
        for (int ii3 = 0; ii3 <= (_PB_N - 3) / 32; ii3 += 1) {
          if (_PB_N >= 32 * ii2 + 35) {
            for (int ii4 = 0; ii4 <= (_PB_N - 3) / 32; ii4 += 1) {
              if (_PB_N >= 32 * ii4 + 35) {
                if (_PB_N >= 32 * ii3 + 35) {
                  for (int i2 = 32 * ii2 + 1; i2 <= 32 * ii2 + 32; i2 += 1) {
                    for (int i3 = 32 * ii3 + 1; i3 <= 32 * ii3 + 32; i3 += 1) {
                      for (int i4 = 32 * ii4 + 1; i4 <= 32 * ii4 + 32; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                  }
                } else {
                  for (int i2 = 32 * ii2 + 1; i2 <= 32 * ii2 + 32; i2 += 1) {
                    for (int i3 = 32 * ii3 + 1; i3 < _PB_N - 1; i3 += 1) {
                      for (int i4 = 32 * ii4 + 1; i4 <= 32 * ii4 + 32; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                  }
                }
                for (int i0 = 32 * ii0 + 2; i0 <= min(_PB_TSTEPS, 32 * ii0 + 17); i0 += 1) {
                  if (ii2 == 0 && ii3 == 0 && ii4 == 0 && i0 >= 32 * ii0 + 10 && 32 * ii0 + 16 >= i0) {
                    for (int i3 = 1; i3 <= 64 * ii0 - 2 * i0 + 35; i3 += 1) {
                      for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 + 35; i4 += 1) {
                        B[1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[2][i3][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[0][i3][i4])) + (SCALAR_VAL(0.125) * ((A[1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3][i4 - 1]))) + A[1][i3][i4]);
                      }
                    }
                  }
                  if (_PB_N >= 32 * ii3 + 35 && i0 >= 32 * ii0 + 3) {
                    for (int i2 = max(64 * ii0 + 32 * ii2 - 2 * i0 + 4, -5 * ii0 - ii2 + floord(3 * ii0 - ii2 + i0 - 3, 7) + 1); i2 <= 64 * ii0 + 32 * ii2 - 2 * i0 + 5; i2 += 1) {
                      for (int i3 = max(1, 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 5); i3 <= 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 36; i3 += 1) {
                        for (int i4 = max(max(1, 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 5), 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 6); i4 <= -64 * ii0 - 32 * ii2 + 32 * ii4 + 2 * i0 + i2 + 26; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        if (2 * i0 + i2 == 64 * ii0 + 32 * ii2 + 5 && i3 == 32 * ii3) {
                          B[64 * ii0 + 32 * ii2 - 2 * i0 + 5][32 * ii3][32 * ii4 + 32] = ((((SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][32 * ii3][32 * ii4 + 32] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][32 * ii3][32 * ii4 + 32])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][32 * ii3][32 * ii4 + 32])) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][32 * ii3 + 1][32 * ii4 + 32] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][32 * ii3][32 * ii4 + 32])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][32 * ii3 - 1][32 * ii4 + 32]))) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][32 * ii3][32 * ii4 + 33] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][32 * ii3][32 * ii4 + 32])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][32 * ii3][32 * ii4 + 31]))) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][32 * ii3][32 * ii4 + 32]);
                        } else if (2 * i0 + i2 == 64 * ii0 + 32 * ii2 + 4) {
                          for (int i4 = 32 * ii4 + 31; i4 <= 32 * ii4 + 32; i4 += 1) {
                            B[64 * ii0 + 32 * ii2 - 2 * i0 + 4][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3][i4] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][i3][i4])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 3][i3][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][i3][i4])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][i3][i4])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][i3][i4 - 1]))) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][i3][i4]);
                          }
                        }
                      }
                    }
                  }
                  if (ii2 >= 1 && _PB_N >= 32 * ii3 + 35) {
                    for (int i2 = 64 * ii0 + 32 * ii2 - 2 * i0 + 6; i2 <= 32 * ii2; i2 += 1) {
                      if (32 * ii2 >= i2 + 1) {
                        for (int i3 = max(1, 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 5); i3 <= 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 36; i3 += 1) {
                          for (int i4 = max(max(1, 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 5), 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 6); i4 <= 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 37; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          for (int i4 = 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 38; i4 <= 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 36; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      } else {
                        if (ii3 >= 1) {
                          for (int i3 = 64 * ii0 + 32 * ii3 - 2 * i0 + 5; i3 <= 32 * ii3; i3 += 1) {
                            for (int i4 = max(1, 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 6); i4 <= 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 37; i4 += 1) {
                              B[32 * ii2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3][i4 - 1]))) + A[32 * ii2][i3][i4]);
                            }
                          }
                        }
                        for (int i3 = 32 * ii3 + 1; i3 <= 64 * ii0 + 32 * ii3 - 2 * i0 + 36; i3 += 1) {
                          if (i0 >= 32 * ii0 + 4) {
                            for (int i4 = max(1, 64 * ii0 + 32 * ii4 - 2 * i0 + 5); i4 <= 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 37; i4 += 1) {
                              B[32 * ii2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3][i4 - 1]))) + A[32 * ii2][i3][i4]);
                            }
                            for (int i4 = 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 38; i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 + 36; i4 += 1) {
                              B[32 * ii2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3][i4 - 1]))) + A[32 * ii2][i3][i4]);
                            }
                          } else {
                            for (int i4 = max(1, 32 * ii4 - 1); i4 <= 32 * ii4; i4 += 1) {
                              B[32 * ii2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3][i4 - 1]))) + A[32 * ii2][i3][i4]);
                            }
                            for (int i4 = 32 * ii4 + 1; i4 <= 32 * ii4 + 30; i4 += 1) {
                              B[32 * ii2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3][i4 - 1]))) + A[32 * ii2][i3][i4]);
                            }
                          }
                        }
                      }
                    }
                  }
                  if (ii3 == 0 && ii4 == 0 && i0 >= 32 * ii0 + 3 && 32 * ii0 + 232 * ii2 + 9 >= i0) {
                    for (int i3 = 1; i3 <= 64 * ii0 - 2 * i0 + 35; i3 += 1) {
                      for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 + 35; i4 += 1) {
                        B[32 * ii2 + 1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3][i4 - 1]))) + A[32 * ii2 + 1][i3][i4]);
                      }
                    }
                  } else if (ii3 == 0 && i0 >= 32 * ii0 + 3 && 32 * ii0 + 16 * ii4 + 1 >= i0) {
                    for (int i3 = 1; i3 <= 64 * ii0 - 2 * i0 + 35; i3 += 1) {
                      for (int i4 = 64 * ii0 + 32 * ii4 - 2 * i0 + 4; i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 + 5; i4 += 1) {
                        B[32 * ii2 + 1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3][i4 - 1]))) + A[32 * ii2 + 1][i3][i4]);
                      }
                      if (ii2 == 0) {
                        for (int i4 = 64 * ii0 + 32 * ii4 - 2 * i0 + 6; i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 - i3 + 36; i4 += 1) {
                          B[1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[2][i3][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[0][i3][i4])) + (SCALAR_VAL(0.125) * ((A[1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3][i4 - 1]))) + A[1][i3][i4]);
                        }
                        for (int i4 = 64 * ii0 + 32 * ii4 - 2 * i0 - i3 + 37; i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 + 35; i4 += 1) {
                          B[1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[2][i3][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[0][i3][i4])) + (SCALAR_VAL(0.125) * ((A[1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3][i4 - 1]))) + A[1][i3][i4]);
                        }
                      } else {
                        for (int i4 = 64 * ii0 + 32 * ii4 - 2 * i0 + 6; i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 + 35; i4 += 1) {
                          B[32 * ii2 + 1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3][i4 - 1]))) + A[32 * ii2 + 1][i3][i4]);
                        }
                      }
                    }
                  } else if (_PB_N >= 32 * ii3 + 35 && i0 >= 32 * ii0 + 3 && 32 * ii0 + 16 * ii3 + 1 >= i0) {
                    for (int i3 = 64 * ii0 + 32 * ii3 - 2 * i0 + 4; i3 <= min(min(32 * ii3, -64 * ii0 - 59 * ii2 + 32 * ii3 + 2 * i0 - 5), 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 + 4); i3 += 1) {
                      B[32 * ii2 + 1][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5])) + A[32 * ii2][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3 + 1][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5])) + A[32 * ii2 + 1][i3 - 1][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 6] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5])) + A[32 * ii2 + 1][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 4]))) + A[32 * ii2 + 1][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5]);
                      if (ii2 == 0 && 2 * i0 + i3 >= 64 * ii0 + 32 * ii3 + 5) {
                        for (int i4 = 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 6; i4 <= 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 36; i4 += 1) {
                          B[1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[2][i3][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[0][i3][i4])) + (SCALAR_VAL(0.125) * ((A[1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3][i4 - 1]))) + A[1][i3][i4]);
                        }
                      } else {
                        for (int i4 = 32 * ii4 + 2; i4 <= 32 * ii4 + 32; i4 += 1) {
                          B[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4])) + A[32 * ii2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 5][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4])) + A[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 3][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4])) + A[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4 - 1]))) + A[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4]);
                        }
                      }
                    }
                    if (ii2 == 0 && ii4 == 0) {
                      for (int i3 = 64 * ii0 + 32 * ii3 - 2 * i0 + 5; i3 < 32 * ii3; i3 += 1) {
                        for (int i4 = 1; i4 <= 64 * ii0 + 32 * ii3 - 2 * i0 - i3 + 36; i4 += 1) {
                          B[1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[2][i3][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[0][i3][i4])) + (SCALAR_VAL(0.125) * ((A[1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3][i4 - 1]))) + A[1][i3][i4]);
                        }
                      }
                      for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 + 36; i4 += 1) {
                        B[1][32 * ii3][i4] = ((((SCALAR_VAL(0.125) * ((A[2][32 * ii3][i4] - (SCALAR_VAL(2.0) * A[1][32 * ii3][i4])) + A[0][32 * ii3][i4])) + (SCALAR_VAL(0.125) * ((A[1][32 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[1][32 * ii3][i4])) + A[1][32 * ii3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[1][32 * ii3][i4 + 1] - (SCALAR_VAL(2.0) * A[1][32 * ii3][i4])) + A[1][32 * ii3][i4 - 1]))) + A[1][32 * ii3][i4]);
                      }
                      for (int i3 = 32 * ii3 + 1; i3 <= min(-64 * ii0 + 32 * ii3 + 2 * i0 - 5, 64 * ii0 + 32 * ii3 - 2 * i0 + 35); i3 += 1) {
                        for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 + 35; i4 += 1) {
                          B[1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[2][i3][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[0][i3][i4])) + (SCALAR_VAL(0.125) * ((A[1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3][i4 - 1]))) + A[1][i3][i4]);
                        }
                      }
                    } else if (128 * ii0 + 59 * ii2 + 8 >= 4 * i0) {
                      for (int i4 = 32 * ii4 + 1; i4 <= 32 * ii4 + 32; i4 += 1) {
                        B[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4])) + A[32 * ii2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 5][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4])) + A[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 3][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4])) + A[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4 - 1]))) + A[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4]);
                      }
                    } else if (ii2 == 0 && 32 * ii0 + 16 * ii4 + 2 >= i0) {
                      for (int i3 = 32 * ii3 + 1; i3 <= min(-64 * ii0 + 32 * ii3 + 2 * i0 - 5, 64 * ii0 + 32 * ii3 - 2 * i0 + 35); i3 += 1) {
                        for (int i4 = 64 * ii0 + 32 * ii4 - 2 * i0 + 4; i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 + 35; i4 += 1) {
                          B[1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[2][i3][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[0][i3][i4])) + (SCALAR_VAL(0.125) * ((A[1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3][i4 - 1]))) + A[1][i3][i4]);
                        }
                      }
                    }
                    for (int i3 = max(-64 * ii0 - 59 * ii2 + 32 * ii3 + 2 * i0 - 4, 64 * ii0 + 32 * ii3 - 2 * i0 + 5); i3 <= 64 * ii0 + 32 * ii3 - 2 * i0 + 35; i3 += 1) {
                      for (int i4 = max(max(1, 64 * ii0 + 32 * ii4 - 2 * i0 + 4), 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5); i4 < 32 * ii4; i4 += 1) {
                        if (2 * i0 + i4 >= 64 * ii0 + 32 * ii4 + 4) {
                          B[32 * ii2 + 1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3][i4 - 1]))) + A[32 * ii2 + 1][i3][i4]);
                        }
                      }
                      for (int i4 = max(1, 32 * ii4); i4 <= min(64 * ii3 + 32 * ii4 - 2 * i3, 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 36); i4 += 1) {
                        if (2 * i0 + i3 + i4 >= 64 * ii0 + 32 * ii3 + 32 * ii4 + 6 || 1) {
                          B[32 * ii2 + 1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3][i4 - 1]))) + A[32 * ii2 + 1][i3][i4]);
                        }
                      }
                      if (32 * ii3 >= i3 + 1) {
                        for (int i4 = 64 * ii3 + 32 * ii4 - 2 * i3 + 1; i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 + 35; i4 += 1) {
                          B[32 * ii2 + 1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3][i4 - 1]))) + A[32 * ii2 + 1][i3][i4]);
                        }
                      } else {
                        for (int i4 = max(max(1, 32 * ii4), 64 * ii3 + 32 * ii4 - 2 * i3 + 1); i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 + 35; i4 += 1) {
                          B[32 * ii2 + 1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3][i4 - 1]))) + A[32 * ii2 + 1][i3][i4]);
                        }
                      }
                      for (int i4 = max(64 * ii0 + 32 * ii4 - 2 * i0 + 36, 64 * ii3 + 32 * ii4 - 2 * i3 + 1); i4 <= 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 36; i4 += 1) {
                        B[32 * ii2 + 1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3][i4 - 1]))) + A[32 * ii2 + 1][i3][i4]);
                      }
                    }
                  } else {
                    if (ii2 >= 1 && i0 == 32 * ii0 + 2) {
                      for (int i3 = 32 * ii3 + 1; i3 <= min(_PB_N - 2, 32 * ii3 + 32); i3 += 1) {
                        for (int i4 = 32 * ii4 + 1; i4 <= 32 * ii4 + 32; i4 += 1) {
                          B[32 * ii2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3][i4 - 1]))) + A[32 * ii2][i3][i4]);
                        }
                      }
                    }
                    if (ii3 == 0 && i0 == 32 * ii0 + 2) {
                      for (int i3 = 1; i3 <= 31; i3 += 1) {
                        for (int i4 = max(1, 32 * ii4); i4 <= 32 * ii4 + 31; i4 += 1) {
                          B[32 * ii2 + 1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3][i4 - 1]))) + A[32 * ii2 + 1][i3][i4]);
                        }
                      }
                    } else if (ii3 >= 1 && _PB_N >= 32 * ii3 + 35 && i0 == 32 * ii0 + 2) {
                      for (int i3 = 32 * ii3; i3 <= 32 * ii3 + 31; i3 += 1) {
                        for (int i4 = max(max(1, 32 * ii4), 64 * ii3 + 32 * ii4 - 2 * i3 + 1); i4 <= 32 * ii4 + 31; i4 += 1) {
                          B[32 * ii2 + 1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3][i4 - 1]))) + A[32 * ii2 + 1][i3][i4]);
                        }
                        if (i3 == 32 * ii3) {
                          B[32 * ii2 + 1][32 * ii3][32 * ii4 + 32] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][32 * ii3][32 * ii4 + 32] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][32 * ii3][32 * ii4 + 32])) + A[32 * ii2][32 * ii3][32 * ii4 + 32])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][32 * ii3 + 1][32 * ii4 + 32] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][32 * ii3][32 * ii4 + 32])) + A[32 * ii2 + 1][32 * ii3 - 1][32 * ii4 + 32]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][32 * ii3][32 * ii4 + 33] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][32 * ii3][32 * ii4 + 32])) + A[32 * ii2 + 1][32 * ii3][32 * ii4 + 31]))) + A[32 * ii2 + 1][32 * ii3][32 * ii4 + 32]);
                        }
                      }
                    }
                    if (_PB_N >= 32 * ii3 + 35 && i0 == 32 * ii0 + 2) {
                      for (int i3 = max(32 * ii3, 31 * ii3 + 1); i3 <= 32 * ii3 + 31; i3 += 1) {
                        for (int i4 = max(max(1, 32 * ii4), 64 * ii3 + 32 * ii4 - 2 * i3 + 1); i4 <= min(32 * ii4 + 31, -32 * ii3 + 32 * ii4 + i3 + 30); i4 += 1) {
                          B[32 * ii2 + 2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 3][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][i3][i4])) + A[32 * ii2 + 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][i3][i4])) + A[32 * ii2 + 2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][i3][i4])) + A[32 * ii2 + 2][i3][i4 - 1]))) + A[32 * ii2 + 2][i3][i4]);
                        }
                        if (i3 == 32 * ii3) {
                          for (int i4 = 32 * ii4 + 31; i4 <= 32 * ii4 + 32; i4 += 1) {
                            B[32 * ii2 + 2][32 * ii3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 3][32 * ii3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][32 * ii3][i4])) + A[32 * ii2 + 1][32 * ii3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][32 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][32 * ii3][i4])) + A[32 * ii2 + 2][32 * ii3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][32 * ii3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][32 * ii3][i4])) + A[32 * ii2 + 2][32 * ii3][i4 - 1]))) + A[32 * ii2 + 2][32 * ii3][i4]);
                          }
                        }
                      }
                    }
                  }
                  if (_PB_N >= 32 * ii3 + 35) {
                    for (int i2 = max(32 * ii2 + 2, 64 * ii0 + 32 * ii2 - 2 * i0 + 7); i2 <= 64 * ii0 + 32 * ii2 - 2 * i0 + 35; i2 += 1) {
                      if (i0 >= 32 * ii0 + 3 && 32 * ii0 + 16 * ii3 + 1 >= i0 && 64 * ii0 + i2 + 1 >= 32 * ii2 + 2 * i0) {
                        for (int i4 = 32 * ii4 + 1; i4 <= 32 * ii4 + 32; i4 += 1) {
                          B[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4])) + A[i2 - 1][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 5][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4])) + A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 3][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4])) + A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4 - 1]))) + A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4]);
                        }
                      } else if (32 * ii0 + 16 * ii3 + 1 >= i0 && 32 * ii2 + 2 * i0 >= 64 * ii0 + i2 + 2) {
                        for (int i4 = 32 * ii4 + 1; i4 <= 32 * ii4 + 32; i4 += 1) {
                          B[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4])) + A[i2 - 1][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 5][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4])) + A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 3][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4])) + A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4 - 1]))) + A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4]);
                        }
                      }
                      for (int i3 = max(1, 64 * ii0 + 32 * ii3 - 2 * i0 + 5); i3 <= min(32 * ii3 - 1, 32 * ii2 + 32 * ii3 - i2 + 2); i3 += 1) {
                        if (64 * ii0 + 32 * ii3 + 32 * ii4 + 4 >= 2 * i0 + i3) {
                          B[i2][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5] - (SCALAR_VAL(2.0) * A[i2][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5])) + A[i2 - 1][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5] - (SCALAR_VAL(2.0) * A[i2][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5])) + A[i2][i3 - 1][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 6] - (SCALAR_VAL(2.0) * A[i2][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5])) + A[i2][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 4]))) + A[i2][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5]);
                        }
                        for (int i4 = max(1, 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 6); i4 <= 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 36; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                      for (int i3 = max(max(1, 64 * ii0 + 32 * ii3 - 2 * i0 + 5), 32 * ii2 + 32 * ii3 - i2 + 3); i3 < 32 * ii3; i3 += 1) {
                        for (int i4 = max(1, 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5); i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 + 35; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = 64 * ii0 + 32 * ii4 - 2 * i0 + 36; i4 <= 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 36; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                      for (int i3 = max(1, 32 * ii3); i3 <= 64 * ii0 + 32 * ii3 - 2 * i0 + 35; i3 += 1) {
                        if (i0 >= 32 * ii0 + 3 && 32 * ii0 + 16 * ii4 + 1 >= i0 && i3 >= 32 * ii3 + 1) {
                          B[i2][i3][64 * ii0 + 32 * ii4 - 2 * i0 + 4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][64 * ii0 + 32 * ii4 - 2 * i0 + 4] - (SCALAR_VAL(2.0) * A[i2][i3][64 * ii0 + 32 * ii4 - 2 * i0 + 4])) + A[i2 - 1][i3][64 * ii0 + 32 * ii4 - 2 * i0 + 4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][64 * ii0 + 32 * ii4 - 2 * i0 + 4] - (SCALAR_VAL(2.0) * A[i2][i3][64 * ii0 + 32 * ii4 - 2 * i0 + 4])) + A[i2][i3 - 1][64 * ii0 + 32 * ii4 - 2 * i0 + 4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][64 * ii0 + 32 * ii4 - 2 * i0 + 5] - (SCALAR_VAL(2.0) * A[i2][i3][64 * ii0 + 32 * ii4 - 2 * i0 + 4])) + A[i2][i3][64 * ii0 + 32 * ii4 - 2 * i0 + 3]))) + A[i2][i3][64 * ii0 + 32 * ii4 - 2 * i0 + 4]);
                        } else if (i0 >= 32 * ii0 + 3 && 32 * ii0 + 16 * ii4 + 2 >= i0 && i3 == 32 * ii3) {
                          B[i2][32 * ii3][64 * ii0 + 32 * ii4 - 2 * i0 + 5] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii3][64 * ii0 + 32 * ii4 - 2 * i0 + 5] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][64 * ii0 + 32 * ii4 - 2 * i0 + 5])) + A[i2 - 1][32 * ii3][64 * ii0 + 32 * ii4 - 2 * i0 + 5])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 1][64 * ii0 + 32 * ii4 - 2 * i0 + 5] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][64 * ii0 + 32 * ii4 - 2 * i0 + 5])) + A[i2][32 * ii3 - 1][64 * ii0 + 32 * ii4 - 2 * i0 + 5]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3][64 * ii0 + 32 * ii4 - 2 * i0 + 6] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][64 * ii0 + 32 * ii4 - 2 * i0 + 5])) + A[i2][32 * ii3][64 * ii0 + 32 * ii4 - 2 * i0 + 4]))) + A[i2][32 * ii3][64 * ii0 + 32 * ii4 - 2 * i0 + 5]);
                        }
                        if (64 * ii0 + 32 * ii2 + 32 * ii3 + 36 >= 2 * i0 + i2 + i3) {
                          for (int i4 = max(max(1, 64 * ii0 + 32 * ii4 - 2 * i0 + 5), 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 6); i4 <= 64 * ii3 + 32 * ii4 - 2 * i3; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          for (int i4 = max(max(1, 64 * ii0 + 32 * ii4 - 2 * i0 + 5), 64 * ii3 + 32 * ii4 - 2 * i3 + 1); i4 < 32 * ii4; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          for (int i4 = max(max(1, 32 * ii4), 64 * ii3 + 32 * ii4 - 2 * i3 + 1); i4 <= min(min(64 * ii0 + 32 * ii4 - 2 * i0 + 35, 64 * ii3 + 32 * ii4 - 2 * i3 + 32), 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 37); i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (i0 == 32 * ii0 + 2) {
                            for (int i4 = max(max(1, 32 * ii4), 64 * ii3 + 32 * ii4 - 2 * i3 + 33); i4 <= min(32 * ii4 + 31, 32 * ii3 + 32 * ii4 - i3 + 33); i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          } else if (i0 >= 32 * ii0 + 3) {
                            for (int i4 = max(max(1, 32 * ii4), 64 * ii3 + 32 * ii4 - 2 * i3 + 33); i4 <= 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 37; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            for (int i4 = max(64 * ii3 + 32 * ii4 - 2 * i3 + 33, 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 38); i4 <= min(64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 36, 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 37); i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            for (int i4 = max(64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 37, 64 * ii3 + 32 * ii4 - 2 * i3 + 33); i4 <= 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 37; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                          if (i3 == 32 * ii3) {
                            B[i2][32 * ii3][64 * ii0 + 32 * ii4 - 2 * i0 + 36] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii3][64 * ii0 + 32 * ii4 - 2 * i0 + 36] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][64 * ii0 + 32 * ii4 - 2 * i0 + 36])) + A[i2 - 1][32 * ii3][64 * ii0 + 32 * ii4 - 2 * i0 + 36])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 1][64 * ii0 + 32 * ii4 - 2 * i0 + 36] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][64 * ii0 + 32 * ii4 - 2 * i0 + 36])) + A[i2][32 * ii3 - 1][64 * ii0 + 32 * ii4 - 2 * i0 + 36]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3][64 * ii0 + 32 * ii4 - 2 * i0 + 37] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][64 * ii0 + 32 * ii4 - 2 * i0 + 36])) + A[i2][32 * ii3][64 * ii0 + 32 * ii4 - 2 * i0 + 35]))) + A[i2][32 * ii3][64 * ii0 + 32 * ii4 - 2 * i0 + 36]);
                          }
                          for (int i4 = 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 38; i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 + 35; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        } else {
                          for (int i4 = max(1, 64 * ii0 + 32 * ii4 - 2 * i0 + 5); i4 < 32 * ii4; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          for (int i4 = max(1, 32 * ii4); i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 + 35; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                    }
                  } else {
                    if (_PB_N >= 32 * ii3 + 4 && i0 == 32 * ii0 + 2) {
                      for (int i3 = 32 * ii3; i3 < _PB_N - 1; i3 += 1) {
                        for (int i4 = max(max(1, 32 * ii4), 16 * ii3 + 32 * ii4 - i3 + i3 / 2 + 1); i4 <= 32 * ii4 + 31; i4 += 1) {
                          B[32 * ii2 + 1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3][i4 - 1]))) + A[32 * ii2 + 1][i3][i4]);
                        }
                        if (i3 == 32 * ii3) {
                          B[32 * ii2 + 1][32 * ii3][32 * ii4 + 32] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][32 * ii3][32 * ii4 + 32] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][32 * ii3][32 * ii4 + 32])) + A[32 * ii2][32 * ii3][32 * ii4 + 32])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][32 * ii3 + 1][32 * ii4 + 32] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][32 * ii3][32 * ii4 + 32])) + A[32 * ii2 + 1][32 * ii3 - 1][32 * ii4 + 32]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][32 * ii3][32 * ii4 + 33] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][32 * ii3][32 * ii4 + 32])) + A[32 * ii2 + 1][32 * ii3][32 * ii4 + 31]))) + A[32 * ii2 + 1][32 * ii3][32 * ii4 + 32]);
                        }
                      }
                    } else if (32 * ii3 + 3 == _PB_N && i0 == 32 * ii0 + 2) {
                      for (int i3 = _PB_N - 3; i3 < _PB_N - 1; i3 += 1) {
                        for (int i4 = max(1, _PB_N + 32 * ii4 - i3 - 2); i4 <= _PB_N + 32 * ii4 - i3 + 29; i4 += 1) {
                          B[32 * ii2 + 1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3][i4 - 1]))) + A[32 * ii2 + 1][i3][i4]);
                        }
                      }
                    }
                    if (i0 == 32 * ii0 + 2) {
                      for (int i3 = 32 * ii3; i3 < _PB_N - 1; i3 += 1) {
                        for (int i4 = max(max(1, 32 * ii4), 16 * ii3 + 32 * ii4 - i3 + i3 / 2 + 1); i4 <= min(32 * ii4 + 31, -32 * ii3 + 32 * ii4 + i3 + 30); i4 += 1) {
                          B[32 * ii2 + 2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 3][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][i3][i4])) + A[32 * ii2 + 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][i3][i4])) + A[32 * ii2 + 2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][i3][i4])) + A[32 * ii2 + 2][i3][i4 - 1]))) + A[32 * ii2 + 2][i3][i4]);
                        }
                        if (i3 == 32 * ii3) {
                          for (int i4 = 32 * ii4 + 31; i4 <= 32 * ii4 + 32; i4 += 1) {
                            B[32 * ii2 + 2][32 * ii3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 3][32 * ii3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][32 * ii3][i4])) + A[32 * ii2 + 1][32 * ii3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][32 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][32 * ii3][i4])) + A[32 * ii2 + 2][32 * ii3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][32 * ii3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][32 * ii3][i4])) + A[32 * ii2 + 2][32 * ii3][i4 - 1]))) + A[32 * ii2 + 2][32 * ii3][i4]);
                          }
                        }
                      }
                      for (int i2 = 32 * ii2 + 3; i2 <= 32 * ii2 + 31; i2 += 1) {
                        for (int i3 = 32 * ii3; i3 < _PB_N - 1; i3 += 1) {
                          for (int i4 = max(max(1, 32 * ii4), 16 * ii3 + 32 * ii4 - i3 + i3 / 2 + 1); i4 <= 32 * ii4 + 31; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (i3 == 32 * ii3) {
                            B[i2][32 * ii3][32 * ii4 + 32] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii3][32 * ii4 + 32] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][32 * ii4 + 32])) + A[i2 - 1][32 * ii3][32 * ii4 + 32])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 1][32 * ii4 + 32] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][32 * ii4 + 32])) + A[i2][32 * ii3 - 1][32 * ii4 + 32]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3][32 * ii4 + 33] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][32 * ii4 + 32])) + A[i2][32 * ii3][32 * ii4 + 31]))) + A[i2][32 * ii3][32 * ii4 + 32]);
                          }
                        }
                      }
                    }
                  }
                  if (32 * ii3 + 34 >= _PB_N && i0 >= 32 * ii0 + 3) {
                    if (_PB_N >= 32 * ii3 + 4) {
                      for (int i2 = max(1, 64 * ii0 + 32 * ii2 - 2 * i0 + 4); i2 < 32 * ii2; i2 += 1) {
                        for (int i3 = 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 5; i3 <= min(32 * ii3 + 1, 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 6); i3 += 1) {
                          for (int i4 = max(1, 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 6); i4 <= min(-64 * ii0 - 32 * ii2 + 32 * ii4 + 2 * i0 + i2 + 26, 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 37); i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (2 * i0 + i2 + i3 == 64 * ii0 + 32 * ii2 + 32 * ii3 + 5) {
                            for (int i4 = -64 * ii0 - 32 * ii2 + 32 * ii4 + 2 * i0 + i2 + 27; i4 <= 43 * ii0 + 21 * ii2 + 32 * ii4 - i0 - i2 + floord(-ii0 + ii2 - i0 + i2 + 1, 3) + 35; i4 += 1) {
                              B[i2][64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 5][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 5][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 5][i4])) + A[i2 - 1][64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 5][i4])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 6][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 5][i4])) + A[i2][64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 4][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 5][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 5][i4])) + A[i2][64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 5][i4 - 1]))) + A[i2][64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 5][i4]);
                            }
                          }
                        }
                        if (2 * i0 + i2 >= 64 * ii0 + 32 * ii2 + 7) {
                          for (int i3 = 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 7; i3 <= 32 * ii3 + 1; i3 += 1) {
                            for (int i4 = max(1, 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 6); i4 <= 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 37; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                        } else if (2 * i0 + i2 == 64 * ii0 + 32 * ii2 + 6) {
                          for (int i4 = max(1, 32 * ii4 - 1); i4 <= 32 * ii4 + 30; i4 += 1) {
                            B[64 * ii0 + 32 * ii2 - 2 * i0 + 6][32 * ii3 + 1][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 7][32 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][32 * ii3 + 1][i4])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][32 * ii3 + 1][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][32 * ii3 + 2][i4] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][32 * ii3 + 1][i4])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][32 * ii3][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][32 * ii3 + 1][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][32 * ii3 + 1][i4])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][32 * ii3 + 1][i4 - 1]))) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][32 * ii3 + 1][i4]);
                          }
                        }
                        for (int i3 = 32 * ii3 + 2; i3 < _PB_N - 1; i3 += 1) {
                          if (i3 == 32 * ii3 + 2) {
                            for (int i4 = max(1, 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 5); i4 <= 32 * ii4; i4 += 1) {
                              B[i2][32 * ii3 + 2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii3 + 2][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 2][i4])) + A[i2 - 1][32 * ii3 + 2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 3][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 2][i4])) + A[i2][32 * ii3 + 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 2][i4])) + A[i2][32 * ii3 + 2][i4 - 1]))) + A[i2][32 * ii3 + 2][i4]);
                            }
                          } else {
                            for (int i4 = max(1, 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 5); i4 < 32 * ii4; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                          for (int i4 = max(max(max(32 * ii4, 31 * ii4 + 1), 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 5), 21 * ii3 + 32 * ii4 - i3 + (ii3 + i3 - 1) / 3 + 3); i4 <= 32 * ii4 + 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (64 * ii0 + 32 * ii2 + 5 >= 2 * i0 + i2) {
                            for (int i4 = 32 * ii4 + 2; i4 <= 128 * ii0 + 64 * ii2 + 32 * ii4 - 4 * i0 - 2 * i2 + 40; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            if (2 * i0 + i2 == 64 * ii0 + 32 * ii2 + 5) {
                              B[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3][32 * ii4 + 31] = ((((SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3][32 * ii4 + 31] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3][32 * ii4 + 31])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][i3][32 * ii4 + 31])) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3 + 1][32 * ii4 + 31] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3][32 * ii4 + 31])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3 - 1][32 * ii4 + 31]))) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3][32 * ii4 + 32] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3][32 * ii4 + 31])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3][32 * ii4 + 30]))) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3][32 * ii4 + 31]);
                            }
                          } else {
                            for (int i4 = 32 * ii4 + 2; i4 <= 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 36; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                        }
                      }
                      if (ii2 >= 1 && 32 * ii0 + 10 >= i0 && 32 * ii3 + 2 * i0 >= _PB_N + 64 * ii0 + 3) {
                        for (int i3 = 64 * ii0 + 32 * ii3 - 2 * i0 + 5; i3 <= 32 * ii3; i3 += 1) {
                          for (int i4 = max(1, 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 6); i4 <= 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 37; i4 += 1) {
                            B[32 * ii2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3][i4 - 1]))) + A[32 * ii2][i3][i4]);
                          }
                        }
                        for (int i3 = 32 * ii3 + 1; i3 < _PB_N - 1; i3 += 1) {
                          if (32 * ii3 + 2 >= i3) {
                            for (int i4 = max(1, 64 * ii0 + 32 * ii4 - 2 * i0 + 5); i4 <= 32 * ii3 + 32 * ii4 - i3 + 2; i4 += 1) {
                              B[32 * ii2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3][i4 - 1]))) + A[32 * ii2][i3][i4]);
                            }
                          } else {
                            for (int i4 = max(1, 64 * ii0 + 32 * ii4 - 2 * i0 + 5); i4 < 32 * ii4; i4 += 1) {
                              B[32 * ii2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3][i4 - 1]))) + A[32 * ii2][i3][i4]);
                            }
                          }
                          for (int i4 = max(max(32 * ii4, 31 * ii4 + 1), 32 * ii3 + 32 * ii4 - i3 + 3); i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 + 36; i4 += 1) {
                            B[32 * ii2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3][i4 - 1]))) + A[32 * ii2][i3][i4]);
                          }
                        }
                      } else if (ii2 >= 1 && ii4 >= 1 && i0 >= 32 * ii0 + 11) {
                        for (int i3 = 64 * ii0 + 32 * ii3 - 2 * i0 + 5; i3 < _PB_N - 1; i3 += 1) {
                          if (32 * ii3 + 2 >= i3) {
                            for (int i4 = max(64 * ii0 + 32 * ii4 - 2 * i0 + 5, 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 6); i4 <= min(32 * ii4 + 1, 32 * ii3 + 32 * ii4 - i3 + 2); i4 += 1) {
                              B[32 * ii2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3][i4 - 1]))) + A[32 * ii2][i3][i4]);
                            }
                          } else {
                            for (int i4 = 64 * ii0 + 32 * ii4 - 2 * i0 + 5; i4 < 32 * ii4; i4 += 1) {
                              B[32 * ii2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3][i4 - 1]))) + A[32 * ii2][i3][i4]);
                            }
                          }
                          if (i3 >= 32 * ii3 + 2) {
                            for (int i4 = max(32 * ii4, 21 * ii3 + 32 * ii4 - i3 + (ii3 + i3 - 1) / 3 + 3); i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 + 36; i4 += 1) {
                              B[32 * ii2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3][i4 - 1]))) + A[32 * ii2][i3][i4]);
                            }
                          } else {
                            for (int i4 = 32 * ii4 + 2; i4 <= 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 37; i4 += 1) {
                              B[32 * ii2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3][i4 - 1]))) + A[32 * ii2][i3][i4]);
                            }
                          }
                        }
                      } else if (ii2 >= 1 && ii4 == 0 && i0 >= 32 * ii0 + 11) {
                        for (int i3 = 64 * ii0 + 32 * ii3 - 2 * i0 + 5; i3 < _PB_N - 1; i3 += 1) {
                          if (i3 >= 32 * ii3 + 2) {
                            for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 + 36; i4 += 1) {
                              B[32 * ii2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3][i4 - 1]))) + A[32 * ii2][i3][i4]);
                            }
                          } else {
                            if (32 * ii3 >= i3) {
                              B[32 * ii2][i3][1] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][1] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][1])) + A[32 * ii2 - 1][i3][1])) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3 + 1][1] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][1])) + A[32 * ii2][i3 - 1][1]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3][2] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][1])) + A[32 * ii2][i3][0]))) + A[32 * ii2][i3][1]);
                            } else {
                              B[32 * ii2][32 * ii3 + 1][1] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][32 * ii3 + 1][1] - (SCALAR_VAL(2.0) * A[32 * ii2][32 * ii3 + 1][1])) + A[32 * ii2 - 1][32 * ii3 + 1][1])) + (SCALAR_VAL(0.125) * ((A[32 * ii2][32 * ii3 + 2][1] - (SCALAR_VAL(2.0) * A[32 * ii2][32 * ii3 + 1][1])) + A[32 * ii2][32 * ii3][1]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2][32 * ii3 + 1][2] - (SCALAR_VAL(2.0) * A[32 * ii2][32 * ii3 + 1][1])) + A[32 * ii2][32 * ii3 + 1][0]))) + A[32 * ii2][32 * ii3 + 1][1]);
                            }
                            for (int i4 = 2; i4 <= 64 * ii0 + 32 * ii3 - 2 * i0 - i3 + 37; i4 += 1) {
                              B[32 * ii2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3][i4 - 1]))) + A[32 * ii2][i3][i4]);
                            }
                          }
                        }
                      } else if (ii2 >= 1 && i0 == 32 * ii0 + 3) {
                        for (int i3 = 32 * ii3 - 1; i3 < _PB_N - 1; i3 += 1) {
                          if (32 * ii3 + 2 >= i3) {
                            for (int i4 = max(max(1, 32 * ii4 - 1), 32 * ii3 + 32 * ii4 - i3); i4 <= min(32 * ii4 + 1, 32 * ii3 + 32 * ii4 - i3 + 2); i4 += 1) {
                              B[32 * ii2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3][i4 - 1]))) + A[32 * ii2][i3][i4]);
                            }
                          } else if (ii4 >= 1) {
                            B[32 * ii2][i3][32 * ii4 - 1] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][32 * ii4 - 1] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][32 * ii4 - 1])) + A[32 * ii2 - 1][i3][32 * ii4 - 1])) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3 + 1][32 * ii4 - 1] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][32 * ii4 - 1])) + A[32 * ii2][i3 - 1][32 * ii4 - 1]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3][32 * ii4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][32 * ii4 - 1])) + A[32 * ii2][i3][32 * ii4 - 2]))) + A[32 * ii2][i3][32 * ii4 - 1]);
                          }
                          if (i3 >= 32 * ii3 + 1) {
                            for (int i4 = max(max(32 * ii4, 31 * ii4 + 1), 32 * ii3 + 32 * ii4 - i3 + 3); i4 <= 32 * ii4 + 30; i4 += 1) {
                              B[32 * ii2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3][i4 - 1]))) + A[32 * ii2][i3][i4]);
                            }
                          } else {
                            for (int i4 = 32 * ii4 + 2; i4 <= 32 * ii3 + 32 * ii4 - i3 + 31; i4 += 1) {
                              B[32 * ii2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3][i4 - 1]))) + A[32 * ii2][i3][i4]);
                            }
                          }
                        }
                      } else if (ii2 >= 1 && i0 >= 32 * ii0 + 4 && 32 * ii0 + 10 >= i0 && _PB_N + 64 * ii0 + 2 >= 32 * ii3 + 2 * i0) {
                        for (int i3 = 64 * ii0 + 32 * ii3 - 2 * i0 + 5; i3 < _PB_N - 1; i3 += 1) {
                          if (i3 >= 32 * ii3 + 1 && 32 * ii3 + 2 >= i3) {
                            for (int i4 = max(1, 64 * ii0 + 32 * ii4 - 2 * i0 + 5); i4 <= 32 * ii4; i4 += 1) {
                              B[32 * ii2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3][i4 - 1]))) + A[32 * ii2][i3][i4]);
                            }
                          } else if (32 * ii3 >= i3) {
                            for (int i4 = max(1, 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 6); i4 <= 32 * ii4 + 1; i4 += 1) {
                              B[32 * ii2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3][i4 - 1]))) + A[32 * ii2][i3][i4]);
                            }
                            for (int i4 = 32 * ii4 + 2; i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 + 35; i4 += 1) {
                              B[32 * ii2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3][i4 - 1]))) + A[32 * ii2][i3][i4]);
                            }
                          } else {
                            for (int i4 = max(1, 64 * ii0 + 32 * ii4 - 2 * i0 + 5); i4 < 32 * ii4; i4 += 1) {
                              B[32 * ii2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3][i4 - 1]))) + A[32 * ii2][i3][i4]);
                            }
                            if (ii4 >= 1 && 64 * ii0 + 32 * ii3 + 37 >= 2 * i0 + i3) {
                              B[32 * ii2][i3][32 * ii4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][32 * ii4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][32 * ii4])) + A[32 * ii2 - 1][i3][32 * ii4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3 + 1][32 * ii4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][32 * ii4])) + A[32 * ii2][i3 - 1][32 * ii4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3][32 * ii4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][32 * ii4])) + A[32 * ii2][i3][32 * ii4 - 1]))) + A[32 * ii2][i3][32 * ii4]);
                            }
                          }
                          if (i3 >= 32 * ii3 + 1) {
                            for (int i4 = 32 * ii4 + 1; i4 <= min(min(64 * ii0 + 32 * ii4 - 2 * i0 + 35, 64 * ii3 + 32 * ii4 - 2 * i3 + 32), 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 37); i4 += 1) {
                              B[32 * ii2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3][i4 - 1]))) + A[32 * ii2][i3][i4]);
                            }
                          }
                          for (int i4 = 64 * ii0 + 32 * ii4 - 2 * i0 + 36; i4 <= 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 37; i4 += 1) {
                            B[32 * ii2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3][i4 - 1]))) + A[32 * ii2][i3][i4]);
                          }
                          for (int i4 = max(32 * ii4 + 1, 64 * ii3 + 32 * ii4 - 2 * i3 + 33); i4 <= 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 37; i4 += 1) {
                            B[32 * ii2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3][i4 - 1]))) + A[32 * ii2][i3][i4]);
                          }
                          for (int i4 = max(max(32 * ii4, 31 * ii4 + 1), 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 38); i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 + 36; i4 += 1) {
                            B[32 * ii2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3][i4 - 1]))) + A[32 * ii2][i3][i4]);
                          }
                        }
                      }
                      for (int i3 = 64 * ii0 + 32 * ii3 - 2 * i0 + 4; i3 <= min(_PB_N - 2, 32 * ii3 + 32 * ii4 + 1); i3 += 1) {
                        if (32 * ii3 + 2 >= i3) {
                          for (int i4 = max(max(1, 64 * ii0 + 32 * ii4 - 2 * i0 + 4), 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5); i4 <= min(32 * ii4 + 1, 16 * ii3 + 32 * ii4 - i3 + i3 / 2); i4 += 1) {
                            B[32 * ii2 + 1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3][i4 - 1]))) + A[32 * ii2 + 1][i3][i4]);
                          }
                          if (i3 + 1 == 32 * ii3) {
                            B[32 * ii2 + 1][32 * ii3 - 1][32 * ii4 + 1] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][32 * ii3 - 1][32 * ii4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][32 * ii3 - 1][32 * ii4 + 1])) + A[32 * ii2][32 * ii3 - 1][32 * ii4 + 1])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][32 * ii3][32 * ii4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][32 * ii3 - 1][32 * ii4 + 1])) + A[32 * ii2 + 1][32 * ii3 - 2][32 * ii4 + 1]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][32 * ii3 - 1][32 * ii4 + 2] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][32 * ii3 - 1][32 * ii4 + 1])) + A[32 * ii2 + 1][32 * ii3 - 1][32 * ii4]))) + A[32 * ii2 + 1][32 * ii3 - 1][32 * ii4 + 1]);
                          }
                        } else {
                          for (int i4 = 64 * ii0 + 32 * ii4 - 2 * i0 + 4; i4 < 32 * ii4; i4 += 1) {
                            B[32 * ii2 + 1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3][i4 - 1]))) + A[32 * ii2 + 1][i3][i4]);
                          }
                        }
                        if (i3 >= 32 * ii3) {
                          for (int i4 = max(max(1, 32 * ii4), 16 * ii3 + 32 * ii4 - i3 + i3 / 2 + 1); i4 <= min(min(64 * ii0 + 32 * ii4 - 2 * i0 + 35, 64 * ii3 + 32 * ii4 - 2 * i3 + 32), 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 37); i4 += 1) {
                            B[32 * ii2 + 1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3][i4 - 1]))) + A[32 * ii2 + 1][i3][i4]);
                          }
                        }
                        for (int i4 = max(32 * ii4, 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 38); i4 <= min(64 * ii0 + 32 * ii4 - 2 * i0 + 35, 64 * ii3 + 32 * ii4 - 2 * i3 + 32); i4 += 1) {
                          B[32 * ii2 + 1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3][i4 - 1]))) + A[32 * ii2 + 1][i3][i4]);
                        }
                        for (int i4 = max(32 * ii4, 64 * ii3 + 32 * ii4 - 2 * i3 + 33); i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 + 35; i4 += 1) {
                          B[32 * ii2 + 1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3][i4 - 1]))) + A[32 * ii2 + 1][i3][i4]);
                        }
                        if (i3 == 32 * ii3) {
                          B[32 * ii2 + 1][32 * ii3][64 * ii0 + 32 * ii4 - 2 * i0 + 36] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][32 * ii3][64 * ii0 + 32 * ii4 - 2 * i0 + 36] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][32 * ii3][64 * ii0 + 32 * ii4 - 2 * i0 + 36])) + A[32 * ii2][32 * ii3][64 * ii0 + 32 * ii4 - 2 * i0 + 36])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][32 * ii3 + 1][64 * ii0 + 32 * ii4 - 2 * i0 + 36] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][32 * ii3][64 * ii0 + 32 * ii4 - 2 * i0 + 36])) + A[32 * ii2 + 1][32 * ii3 - 1][64 * ii0 + 32 * ii4 - 2 * i0 + 36]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][32 * ii3][64 * ii0 + 32 * ii4 - 2 * i0 + 37] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][32 * ii3][64 * ii0 + 32 * ii4 - 2 * i0 + 36])) + A[32 * ii2 + 1][32 * ii3][64 * ii0 + 32 * ii4 - 2 * i0 + 35]))) + A[32 * ii2 + 1][32 * ii3][64 * ii0 + 32 * ii4 - 2 * i0 + 36]);
                        } else if (32 * ii3 >= i3 + 1) {
                          for (int i4 = 32 * ii4 + 2; i4 <= 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 36; i4 += 1) {
                            B[32 * ii2 + 1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3][i4 - 1]))) + A[32 * ii2 + 1][i3][i4]);
                          }
                        }
                      }
                      if (ii4 == 0) {
                        for (int i3 = 32 * ii3 + 2; i3 <= min(min(_PB_N - 2, 32 * ii3 + 15), 64 * ii0 + 32 * ii3 - 2 * i0 + 36); i3 += 1) {
                          B[32 * ii2 + 1][i3][1] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3][1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][1])) + A[32 * ii2][i3][1])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3 + 1][1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][1])) + A[32 * ii2 + 1][i3 - 1][1]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][2] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][1])) + A[32 * ii2 + 1][i3][0]))) + A[32 * ii2 + 1][i3][1]);
                          for (int i4 = 2; i4 <= 64 * ii0 - 2 * i0 + 35; i4 += 1) {
                            B[32 * ii2 + 1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3][i4 - 1]))) + A[32 * ii2 + 1][i3][i4]);
                          }
                        }
                        if (_PB_N >= 32 * ii3 + 18 && 32 * ii0 + 10 >= i0) {
                          for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 + 35; i4 += 1) {
                            B[32 * ii2 + 1][32 * ii3 + 16][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][32 * ii3 + 16][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][32 * ii3 + 16][i4])) + A[32 * ii2][32 * ii3 + 16][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][32 * ii3 + 17][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][32 * ii3 + 16][i4])) + A[32 * ii2 + 1][32 * ii3 + 15][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][32 * ii3 + 16][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][32 * ii3 + 16][i4])) + A[32 * ii2 + 1][32 * ii3 + 16][i4 - 1]))) + A[32 * ii2 + 1][32 * ii3 + 16][i4]);
                          }
                        }
                        for (int i3 = 32 * ii3 + 17; i3 <= min(_PB_N - 2, 64 * ii0 + 32 * ii3 - 2 * i0 + 36); i3 += 1) {
                          for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 + 35; i4 += 1) {
                            B[32 * ii2 + 1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3][i4 - 1]))) + A[32 * ii2 + 1][i3][i4]);
                          }
                        }
                        for (int i3 = 64 * ii0 + 32 * ii3 - 2 * i0 + 37; i3 < _PB_N - 1; i3 += 1) {
                          for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 + 35; i4 += 1) {
                            B[32 * ii2 + 1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3][i4 - 1]))) + A[32 * ii2 + 1][i3][i4]);
                          }
                        }
                      }
                      for (int i2 = 32 * ii2 + 2; i2 <= min(-64 * ii0 + 32 * ii2 + 2 * i0 - 3, 64 * ii0 + 32 * ii2 - 2 * i0 + 35); i2 += 1) {
                        for (int i3 = 64 * ii0 + 32 * ii3 - 2 * i0 + 4; i3 < _PB_N - 1; i3 += 1) {
                          if (32 * ii3 + 2 >= i3) {
                            for (int i4 = max(max(1, 64 * ii0 + 32 * ii4 - 2 * i0 + 4), 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5); i4 <= min(32 * ii4 + 1, 32 * ii3 + 32 * ii4 - i3 + 2); i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            if (32 * ii3 + 1 >= i3) {
                              for (int i4 = 32 * ii4 + 2; i4 <= 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 36; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                          } else {
                            for (int i4 = max(1, 64 * ii0 + 32 * ii4 - 2 * i0 + 4); i4 < 32 * ii4; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                          if (i3 >= 32 * ii3 + 2) {
                            for (int i4 = max(max(32 * ii4, 31 * ii4 + 1), 21 * ii3 + 32 * ii4 - i3 + (ii3 + i3 - 1) / 3 + 3); i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 + 35; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                        }
                      }
                    } else {
                      if (ii2 >= 1) {
                        for (int i2 = 64 * ii0 + 32 * ii2 - 2 * i0 + 4; i2 <= 32 * ii2; i2 += 1) {
                          for (int i3 = _PB_N + 64 * ii0 + 32 * ii2 - 2 * i0 - i2 + 2; i3 < _PB_N - 1; i3 += 1) {
                            for (int i4 = max(1, _PB_N + 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 - i3 + 3); i4 <= 32 * ii4; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            if (32 * ii2 >= i2 + 1 && i3 + 2 == _PB_N) {
                              B[i2][_PB_N - 2][32 * ii4 + 1] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][_PB_N - 2][32 * ii4 + 1] - (SCALAR_VAL(2.0) * A[i2][_PB_N - 2][32 * ii4 + 1])) + A[i2 - 1][_PB_N - 2][32 * ii4 + 1])) + (SCALAR_VAL(0.125) * ((A[i2][_PB_N - 1][32 * ii4 + 1] - (SCALAR_VAL(2.0) * A[i2][_PB_N - 2][32 * ii4 + 1])) + A[i2][_PB_N - 3][32 * ii4 + 1]))) + (SCALAR_VAL(0.125) * ((A[i2][_PB_N - 2][32 * ii4 + 2] - (SCALAR_VAL(2.0) * A[i2][_PB_N - 2][32 * ii4 + 1])) + A[i2][_PB_N - 2][32 * ii4]))) + A[i2][_PB_N - 2][32 * ii4 + 1]);
                            } else if (_PB_N >= i3 + 3) {
                              B[i2][i3][32 * ii4 + 1] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][32 * ii4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][32 * ii4 + 1])) + A[i2 - 1][i3][32 * ii4 + 1])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][32 * ii4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][32 * ii4 + 1])) + A[i2][i3 - 1][32 * ii4 + 1]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][32 * ii4 + 2] - (SCALAR_VAL(2.0) * A[i2][i3][32 * ii4 + 1])) + A[i2][i3][32 * ii4]))) + A[i2][i3][32 * ii4 + 1]);
                            }
                            if (_PB_N + 64 * ii0 + 32 * ii2 + 3 >= 2 * i0 + i2 + i3) {
                              for (int i4 = 32 * ii4 + 2; i4 <= min(-64 * ii0 - 32 * ii2 + 32 * ii4 + 2 * i0 + i2 + 26, _PB_N + 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 - i3 + 34); i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                              if (2 * i0 + i2 + i3 == _PB_N + 64 * ii0 + 32 * ii2 + 2) {
                                for (int i4 = -64 * ii0 - 32 * ii2 + 32 * ii4 + 2 * i0 + i2 + 27; i4 <= 43 * ii0 + 21 * ii2 + 32 * ii4 - i0 - i2 + floord(-ii0 + ii2 - i0 + i2 + 1, 3) + 35; i4 += 1) {
                                  B[i2][_PB_N + 64 * ii0 + 32 * ii2 - 2 * i0 - i2 + 2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][_PB_N + 64 * ii0 + 32 * ii2 - 2 * i0 - i2 + 2][i4] - (SCALAR_VAL(2.0) * A[i2][_PB_N + 64 * ii0 + 32 * ii2 - 2 * i0 - i2 + 2][i4])) + A[i2 - 1][_PB_N + 64 * ii0 + 32 * ii2 - 2 * i0 - i2 + 2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][_PB_N + 64 * ii0 + 32 * ii2 - 2 * i0 - i2 + 3][i4] - (SCALAR_VAL(2.0) * A[i2][_PB_N + 64 * ii0 + 32 * ii2 - 2 * i0 - i2 + 2][i4])) + A[i2][_PB_N + 64 * ii0 + 32 * ii2 - 2 * i0 - i2 + 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][_PB_N + 64 * ii0 + 32 * ii2 - 2 * i0 - i2 + 2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][_PB_N + 64 * ii0 + 32 * ii2 - 2 * i0 - i2 + 2][i4])) + A[i2][_PB_N + 64 * ii0 + 32 * ii2 - 2 * i0 - i2 + 2][i4 - 1]))) + A[i2][_PB_N + 64 * ii0 + 32 * ii2 - 2 * i0 - i2 + 2][i4]);
                                }
                              }
                            } else {
                              if (i0 == 32 * ii0 + 3 && i2 == 32 * ii2 && i3 + 2 == _PB_N) {
                                B[32 * ii2][_PB_N - 2][32 * ii4 + 1] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][_PB_N - 2][32 * ii4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2][_PB_N - 2][32 * ii4 + 1])) + A[32 * ii2 - 1][_PB_N - 2][32 * ii4 + 1])) + (SCALAR_VAL(0.125) * ((A[32 * ii2][_PB_N - 1][32 * ii4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2][_PB_N - 2][32 * ii4 + 1])) + A[32 * ii2][_PB_N - 3][32 * ii4 + 1]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2][_PB_N - 2][32 * ii4 + 2] - (SCALAR_VAL(2.0) * A[32 * ii2][_PB_N - 2][32 * ii4 + 1])) + A[32 * ii2][_PB_N - 2][32 * ii4]))) + A[32 * ii2][_PB_N - 2][32 * ii4 + 1]);
                              } else if (i2 == 32 * ii2 && i3 + 2 == _PB_N) {
                                B[32 * ii2][_PB_N - 2][32 * ii4 + 1] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][_PB_N - 2][32 * ii4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2][_PB_N - 2][32 * ii4 + 1])) + A[32 * ii2 - 1][_PB_N - 2][32 * ii4 + 1])) + (SCALAR_VAL(0.125) * ((A[32 * ii2][_PB_N - 1][32 * ii4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2][_PB_N - 2][32 * ii4 + 1])) + A[32 * ii2][_PB_N - 3][32 * ii4 + 1]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2][_PB_N - 2][32 * ii4 + 2] - (SCALAR_VAL(2.0) * A[32 * ii2][_PB_N - 2][32 * ii4 + 1])) + A[32 * ii2][_PB_N - 2][32 * ii4]))) + A[32 * ii2][_PB_N - 2][32 * ii4 + 1]);
                              }
                              for (int i4 = 32 * ii4 + 2; i4 <= _PB_N + 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 - i3 + 34; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                          }
                        }
                      }
                      for (int i2 = 32 * ii2 + 1; i2 <= min(-64 * ii0 + 32 * ii2 + 2 * i0 - 3, 64 * ii0 + 32 * ii2 - 2 * i0 + 35); i2 += 1) {
                        for (int i3 = _PB_N + 64 * ii0 - 2 * i0 + 1; i3 < _PB_N - 1; i3 += 1) {
                          if (_PB_N >= i3 + 4) {
                            for (int i4 = max(1, _PB_N + 64 * ii0 + 32 * ii4 - 2 * i0 - i3 + 2); i4 <= 32 * ii4 + 1; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            if (i2 == 32 * ii2 + 1 && 2 * i0 + i3 == _PB_N + 64 * ii0 + 2) {
                              for (int i4 = 32 * ii4 + 2; i4 <= 32 * ii4 + 31; i4 += 1) {
                                B[32 * ii2 + 1][_PB_N + 64 * ii0 - 2 * i0 + 2][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][_PB_N + 64 * ii0 - 2 * i0 + 2][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][_PB_N + 64 * ii0 - 2 * i0 + 2][i4])) + A[32 * ii2][_PB_N + 64 * ii0 - 2 * i0 + 2][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][_PB_N + 64 * ii0 - 2 * i0 + 3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][_PB_N + 64 * ii0 - 2 * i0 + 2][i4])) + A[32 * ii2 + 1][_PB_N + 64 * ii0 - 2 * i0 + 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][_PB_N + 64 * ii0 - 2 * i0 + 2][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][_PB_N + 64 * ii0 - 2 * i0 + 2][i4])) + A[32 * ii2 + 1][_PB_N + 64 * ii0 - 2 * i0 + 2][i4 - 1]))) + A[32 * ii2 + 1][_PB_N + 64 * ii0 - 2 * i0 + 2][i4]);
                              }
                            }
                            if (2 * i0 + i2 + i3 >= _PB_N + 64 * ii0 + 32 * ii2 + 4 && 32 * ii2 + 2 * i0 + i3 >= _PB_N + 64 * ii0 + i2 + 1) {
                              for (int i4 = 32 * ii4 + 2; i4 <= _PB_N + 64 * ii0 + 32 * ii4 - 2 * i0 - i3 + 33; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                          } else {
                            for (int i4 = max(1, _PB_N + 64 * ii0 + 32 * ii4 - 2 * i0 - i3 + 2); i4 < _PB_N + 32 * ii4 - i3 - 2; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            if (64 * ii0 + i2 + 3 == 32 * ii2 + 2 * i0 && i3 + 3 == _PB_N) {
                              B[-64 * ii0 + 32 * ii2 + 2 * i0 - 3][_PB_N - 3][32 * ii4 + 1] = ((((SCALAR_VAL(0.125) * ((A[-64 * ii0 + 32 * ii2 + 2 * i0 - 2][_PB_N - 3][32 * ii4 + 1] - (SCALAR_VAL(2.0) * A[-64 * ii0 + 32 * ii2 + 2 * i0 - 3][_PB_N - 3][32 * ii4 + 1])) + A[-64 * ii0 + 32 * ii2 + 2 * i0 - 4][_PB_N - 3][32 * ii4 + 1])) + (SCALAR_VAL(0.125) * ((A[-64 * ii0 + 32 * ii2 + 2 * i0 - 3][_PB_N - 2][32 * ii4 + 1] - (SCALAR_VAL(2.0) * A[-64 * ii0 + 32 * ii2 + 2 * i0 - 3][_PB_N - 3][32 * ii4 + 1])) + A[-64 * ii0 + 32 * ii2 + 2 * i0 - 3][_PB_N - 4][32 * ii4 + 1]))) + (SCALAR_VAL(0.125) * ((A[-64 * ii0 + 32 * ii2 + 2 * i0 - 3][_PB_N - 3][32 * ii4 + 2] - (SCALAR_VAL(2.0) * A[-64 * ii0 + 32 * ii2 + 2 * i0 - 3][_PB_N - 3][32 * ii4 + 1])) + A[-64 * ii0 + 32 * ii2 + 2 * i0 - 3][_PB_N - 3][32 * ii4]))) + A[-64 * ii0 + 32 * ii2 + 2 * i0 - 3][_PB_N - 3][32 * ii4 + 1]);
                            } else {
                              for (int i4 = max(1, _PB_N + 32 * ii4 - i3 - 2); i4 <= _PB_N + 64 * ii0 + 32 * ii4 - 2 * i0 - i3 + 33; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                          }
                          if (_PB_N + 64 * ii0 + i2 >= 32 * ii2 + 2 * i0 + i3) {
                            for (int i4 = 32 * ii4 + 2; i4 <= _PB_N + 64 * ii0 + 32 * ii4 - 2 * i0 - i3 + 33; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                        }
                      }
                    }
                    for (int i2 = -64 * ii0 + 32 * ii2 + 2 * i0 - 2; i2 <= 64 * ii0 + 32 * ii2 - 2 * i0 + 35; i2 += 1) {
                      for (int i3 = 64 * ii0 + 32 * ii3 - 2 * i0 + 4; i3 < 32 * ii3; i3 += 1) {
                        for (int i4 = max(1, 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5); i4 <= 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 36; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                      for (int i4 = max(1, 64 * ii0 + 32 * ii4 - 2 * i0 + 5); i4 <= 32 * ii4 + 1; i4 += 1) {
                        B[i2][32 * ii3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii3][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][i4])) + A[i2 - 1][32 * ii3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][i4])) + A[i2][32 * ii3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][i4])) + A[i2][32 * ii3][i4 - 1]))) + A[i2][32 * ii3][i4]);
                      }
                      for (int i4 = 32 * ii4 + 2; i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 + 36; i4 += 1) {
                        B[i2][32 * ii3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii3][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][i4])) + A[i2 - 1][32 * ii3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][i4])) + A[i2][32 * ii3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][i4])) + A[i2][32 * ii3][i4 - 1]))) + A[i2][32 * ii3][i4]);
                      }
                      for (int i3 = 32 * ii3 + 1; i3 <= min(min(_PB_N - 2, 32 * ii3 + 16), 32 * ii3 + 16 * ii4 + 15); i3 += 1) {
                        if (32 * ii3 + 2 >= i3) {
                          for (int i4 = max(1, 64 * ii0 + 32 * ii4 - 2 * i0 + 4); i4 < 32 * ii4; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        } else {
                          for (int i4 = max(1, 64 * ii0 + 32 * ii4 - 2 * i0 + 4); i4 < 32 * ii4; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                        for (int i4 = max(1, 32 * ii4); i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 + 35; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                      for (int i3 = 32 * ii3 + 17; i3 <= min(_PB_N - 2, 32 * ii3 + 16 * ii4 + 15); i3 += 1) {
                        for (int i4 = 64 * ii0 + 32 * ii4 - 2 * i0 + 4; i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 + 35; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                      for (int i3 = 32 * ii3 + 16 * ii4 + 16; i3 < _PB_N - 1; i3 += 1) {
                        if (32 * ii3 + 34 == _PB_N && ii4 == 1 && i3 + 2 == _PB_N) {
                          for (int i4 = 64 * ii0 - 2 * i0 + 36; i4 <= 31; i4 += 1) {
                            B[i2][_PB_N - 2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][_PB_N - 2][i4] - (SCALAR_VAL(2.0) * A[i2][_PB_N - 2][i4])) + A[i2 - 1][_PB_N - 2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][_PB_N - 1][i4] - (SCALAR_VAL(2.0) * A[i2][_PB_N - 2][i4])) + A[i2][_PB_N - 3][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][_PB_N - 2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][_PB_N - 2][i4])) + A[i2][_PB_N - 2][i4 - 1]))) + A[i2][_PB_N - 2][i4]);
                          }
                        }
                        for (int i4 = max(1, 32 * ii4); i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 + 35; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    }
                  } else if (ii2 == 0 && ii3 == 0 && ii4 == 0 && i0 == 32 * ii0 + 17) {
                    B[1][1][1] = ((((SCALAR_VAL(0.125) * ((A[2][1][1] - (SCALAR_VAL(2.0) * A[1][1][1])) + A[0][1][1])) + (SCALAR_VAL(0.125) * ((A[1][2][1] - (SCALAR_VAL(2.0) * A[1][1][1])) + A[1][0][1]))) + (SCALAR_VAL(0.125) * ((A[1][1][2] - (SCALAR_VAL(2.0) * A[1][1][1])) + A[1][1][0]))) + A[1][1][1]);
                  }
                  for (int i2 = max(1, 64 * ii0 + 32 * ii2 - 2 * i0 + 3); i2 <= 64 * ii0 + 32 * ii2 - 2 * i0 + 34; i2 += 1) {
                    if (i2 >= 32 * ii2 + 1) {
                      for (int i3 = max(1, 64 * ii0 + 32 * ii3 - 2 * i0 + 3); i3 <= 32 * ii3; i3 += 1) {
                        for (int i4 = max(1, 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 4); i4 <= 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 35; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                      if (_PB_N >= 32 * ii3 + 35) {
                        for (int i3 = 32 * ii3 + 1; i3 <= 64 * ii0 + 32 * ii3 - 2 * i0 + 34; i3 += 1) {
                          for (int i4 = max(1, 64 * ii0 + 32 * ii4 - 2 * i0 + 3); i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 + 34; i4 += 1) {
                            A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                          }
                        }
                      } else {
                        for (int i3 = 32 * ii3 + 1; i3 < _PB_N - 1; i3 += 1) {
                          for (int i4 = max(1, 64 * ii0 + 32 * ii4 - 2 * i0 + 3); i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 + 34; i4 += 1) {
                            A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                          }
                        }
                      }
                    } else if (_PB_N >= 32 * ii3 + 35) {
                      for (int i3 = max(1, 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 4); i3 <= 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 35; i3 += 1) {
                        for (int i4 = max(max(1, 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 4), 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 5); i4 <= 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 36; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                        for (int i4 = 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 37; i4 <= 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 35; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    } else {
                      for (int i3 = 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 4; i3 <= 32 * ii3 + 1; i3 += 1) {
                        for (int i4 = max(1, 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 5); i4 <= 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 36; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                      for (int i3 = 32 * ii3 + 2; i3 < _PB_N - 1; i3 += 1) {
                        for (int i4 = max(1, 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 4); i4 <= 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 35; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                  }
                }
                if (_PB_TSTEPS >= 32 * ii0 + 18 && ii2 == 1 && ii3 == 0) {
                  for (int i2 = 1; i2 <= 31; i2 += 1) {
                    for (int i3 = 1; i3 <= -i2 + 32; i3 += 1) {
                      for (int i4 = max(1, 32 * ii4 - i2 + 1); i4 <= 32 * ii4 - i2 + 32; i4 += 1) {
                        B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                      }
                    }
                  }
                  for (int i2 = 1; i2 <= 30; i2 += 1) {
                    for (int i3 = 1; i3 <= -i2 + 31; i3 += 1) {
                      for (int i4 = max(1, 32 * ii4 - i2); i4 <= 32 * ii4 - i2 + 31; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                  }
                }
                if (_PB_N >= 32 * ii3 + 35) {
                  for (int i0 = 32 * ii0 + 18; i0 <= min(min(min(_PB_TSTEPS, 32 * ii0 + 32), 32 * ii0 + 16 * ii2 + 2), 32 * ii0 + 16 * ii2 + ii3 + 1); i0 += 1) {
                    for (int i2 = max(1, 64 * ii0 + 32 * ii2 - 2 * i0 + 4); i2 <= 64 * ii0 + 32 * ii2 - 2 * i0 + 35; i2 += 1) {
                      for (int i3 = max(1, 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 5); i3 <= 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 36; i3 += 1) {
                        for (int i4 = max(max(1, 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 5), 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 6); i4 <= 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 6; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        if (64 * ii0 + 32 * ii2 + 5 >= 2 * i0 + i2) {
                          for (int i4 = 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 7; i4 <= -64 * ii0 - 32 * ii2 + 32 * ii4 + 2 * i0 + i2 + 26; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (2 * i0 + i2 == 64 * ii0 + 32 * ii2 + 5 && i3 == 32 * ii3) {
                            B[64 * ii0 + 32 * ii2 - 2 * i0 + 5][32 * ii3][32 * ii4 + 32] = ((((SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][32 * ii3][32 * ii4 + 32] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][32 * ii3][32 * ii4 + 32])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][32 * ii3][32 * ii4 + 32])) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][32 * ii3 + 1][32 * ii4 + 32] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][32 * ii3][32 * ii4 + 32])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][32 * ii3 - 1][32 * ii4 + 32]))) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][32 * ii3][32 * ii4 + 33] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][32 * ii3][32 * ii4 + 32])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][32 * ii3][32 * ii4 + 31]))) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][32 * ii3][32 * ii4 + 32]);
                          }
                          if (2 * i0 + i2 == 64 * ii0 + 32 * ii2 + 4) {
                            for (int i4 = 32 * ii4 + 31; i4 <= 32 * ii4 + 32; i4 += 1) {
                              B[64 * ii0 + 32 * ii2 - 2 * i0 + 4][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3][i4] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][i3][i4])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 3][i3][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][i3][i4])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][i3][i4])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][i3][i4 - 1]))) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][i3][i4]);
                            }
                          }
                        } else {
                          for (int i4 = max(max(1, 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 7), 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 6); i4 <= 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 37; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          for (int i4 = 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 38; i4 <= 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 36; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                    }
                    for (int i2 = max(1, 64 * ii0 + 32 * ii2 - 2 * i0 + 3); i2 <= 64 * ii0 + 32 * ii2 - 2 * i0 + 34; i2 += 1) {
                      for (int i3 = max(1, 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 4); i3 <= 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 35; i3 += 1) {
                        for (int i4 = max(max(1, 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 4), 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 5); i4 <= 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 36; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                        for (int i4 = 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 37; i4 <= 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 35; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                  }
                  if (ii2 == 1) {
                    for (int i0 = 32 * ii0 + 19; i0 <= min(_PB_TSTEPS, 32 * ii0 + 32); i0 += 1) {
                      for (int i2 = 1; i2 <= 64 * ii0 - 2 * i0 + 67; i2 += 1) {
                        for (int i3 = max(1, 64 * ii0 + 32 * ii3 - 2 * i0 - i2 + 37); i3 <= 64 * ii0 + 32 * ii3 - 2 * i0 - i2 + 68; i3 += 1) {
                          for (int i4 = max(max(1, 64 * ii0 + 32 * ii4 - 2 * i0 - i2 + 37), 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 38); i4 <= 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 69; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          for (int i4 = 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 70; i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 - i2 + 68; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                      for (int i2 = 1; i2 <= 64 * ii0 - 2 * i0 + 66; i2 += 1) {
                        for (int i3 = max(1, 64 * ii0 + 32 * ii3 - 2 * i0 - i2 + 36); i3 <= 64 * ii0 + 32 * ii3 - 2 * i0 - i2 + 67; i3 += 1) {
                          for (int i4 = max(max(1, 64 * ii0 + 32 * ii4 - 2 * i0 - i2 + 36), 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 37); i4 <= 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 68; i4 += 1) {
                            A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                          }
                          for (int i4 = 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 69; i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 - i2 + 67; i4 += 1) {
                            A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                          }
                        }
                      }
                    }
                  }
                } else {
                  for (int i0 = 32 * ii0 + 18; i0 <= min(min(_PB_TSTEPS, 32 * ii0 + 32), 32 * ii0 + 16 * ii2 + 17); i0 += 1) {
                    for (int i2 = max(1, 64 * ii0 + 32 * ii2 - 2 * i0 + 4); i2 <= 64 * ii0 + 32 * ii2 - 2 * i0 + 35; i2 += 1) {
                      if (2 * i0 + i2 >= 64 * ii0 + 32 * ii2 + 7) {
                        for (int i3 = 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 5; i3 <= 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 6; i3 += 1) {
                          for (int i4 = max(1, 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 6); i4 <= 32 * ii4 + 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          for (int i4 = 32 * ii4 + 2; i4 <= min(64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 + 60 * i2 - i3 - 24, 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 37); i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (ii2 == 1 && i2 == 1) {
                            for (int i4 = 32 * ii4 + 2; i4 <= 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 68; i4 += 1) {
                              B[1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[2][i3][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[0][i3][i4])) + (SCALAR_VAL(0.125) * ((A[1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3][i4 - 1]))) + A[1][i3][i4]);
                            }
                          }
                        }
                        for (int i3 = 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 7; i3 <= min(32 * ii3 + 1, 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 35); i3 += 1) {
                          for (int i4 = max(1, 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 6); i4 <= 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 37; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                        if (2 * i0 + i2 == 64 * ii0 + 32 * ii2 + 35) {
                          for (int i4 = max(1, 32 * ii4 - 30); i4 <= 32 * ii4 + 1; i4 += 1) {
                            B[64 * ii0 + 32 * ii2 - 2 * i0 + 35][32 * ii3 + 1][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 36][32 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 35][32 * ii3 + 1][i4])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 34][32 * ii3 + 1][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 35][32 * ii3 + 2][i4] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 35][32 * ii3 + 1][i4])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 35][32 * ii3][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 35][32 * ii3 + 1][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 35][32 * ii3 + 1][i4])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 35][32 * ii3 + 1][i4 - 1]))) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 35][32 * ii3 + 1][i4]);
                          }
                        }
                        if (_PB_N >= 32 * ii3 + 4 && ii4 == 0) {
                          for (int i4 = 1; i4 <= 64 * ii0 + 32 * ii2 - 2 * i0 - i2 + 36; i4 += 1) {
                            B[i2][32 * ii3 + 2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii3 + 2][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 2][i4])) + A[i2 - 1][32 * ii3 + 2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 3][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 2][i4])) + A[i2][32 * ii3 + 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 2][i4])) + A[i2][32 * ii3 + 2][i4 - 1]))) + A[i2][32 * ii3 + 2][i4]);
                          }
                        } else if (_PB_N >= 32 * ii3 + 4 && ii4 >= 1) {
                          for (int i4 = 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 5; i4 <= 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 36; i4 += 1) {
                            B[i2][32 * ii3 + 2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii3 + 2][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 2][i4])) + A[i2 - 1][32 * ii3 + 2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 3][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 2][i4])) + A[i2][32 * ii3 + 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 2][i4])) + A[i2][32 * ii3 + 2][i4 - 1]))) + A[i2][32 * ii3 + 2][i4]);
                          }
                        }
                        for (int i3 = 32 * ii3 + 3; i3 < _PB_N - 1; i3 += 1) {
                          for (int i4 = max(1, 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 5); i4 <= 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 36; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      } else if (64 * ii0 + 32 * ii2 + 5 >= 2 * i0 + i2) {
                        for (int i3 = 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 5; i3 <= min(min(_PB_N - 2, 32 * ii3 + 32 * ii4 + 1), -64 * ii0 - 32 * ii2 + 32 * ii3 + 2 * i0 + i2 - 3); i3 += 1) {
                          for (int i4 = max(max(1, 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 5), 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 6); i4 <= min(32 * ii4 + 1, 32 * ii3 + 32 * ii4 - i3 + 2); i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (ii2 == 1 && i0 == 32 * ii0 + 18 && i2 == 1 && 32 * ii3 + 1 >= i3) {
                            for (int i4 = 32 * ii4 + 2; i4 <= 32 * ii4 + 31; i4 += 1) {
                              B[1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[2][i3][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[0][i3][i4])) + (SCALAR_VAL(0.125) * ((A[1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3][i4 - 1]))) + A[1][i3][i4]);
                            }
                          }
                          if (32 * ii3 + 1 >= i3) {
                            for (int i4 = 32 * ii4 + 2; i4 <= min(-64 * ii0 - 32 * ii2 + 32 * ii4 + 2 * i0 + i2 + 26, 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 + 60 * i2 - i3 - 24); i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            if (2 * i0 + i2 + i3 == 64 * ii0 + 32 * ii2 + 32 * ii3 + 5) {
                              for (int i4 = -64 * ii0 - 32 * ii2 + 32 * ii4 + 2 * i0 + i2 + 27; i4 <= 43 * ii0 + 21 * ii2 + 32 * ii4 - i0 - i2 + floord(-ii0 + ii2 - i0 + i2 + 1, 3) + 35; i4 += 1) {
                                B[i2][64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 5][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 5][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 5][i4])) + A[i2 - 1][64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 5][i4])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 6][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 5][i4])) + A[i2][64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 4][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 5][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 5][i4])) + A[i2][64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 5][i4 - 1]))) + A[i2][64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 5][i4]);
                              }
                            }
                          } else {
                            for (int i4 = 32 * ii4 + 1; i4 <= 32 * ii4 + 31; i4 += 1) {
                              B[64 * ii0 + 32 * ii2 - 2 * i0 + 5][32 * ii3 + 2][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][32 * ii3 + 2][i4] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][32 * ii3 + 2][i4])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][32 * ii3 + 2][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][32 * ii3 + 3][i4] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][32 * ii3 + 2][i4])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][32 * ii3 + 1][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][32 * ii3 + 2][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][32 * ii3 + 2][i4])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][32 * ii3 + 2][i4 - 1]))) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][32 * ii3 + 2][i4]);
                            }
                          }
                        }
                        if (_PB_N >= 32 * ii3 + 4 && ii4 == 0 && 2 * i0 + i2 == 64 * ii0 + 32 * ii2 + 5) {
                          for (int i4 = 1; i4 <= 31; i4 += 1) {
                            B[64 * ii0 + 32 * ii2 - 2 * i0 + 5][32 * ii3 + 2][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][32 * ii3 + 2][i4] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][32 * ii3 + 2][i4])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][32 * ii3 + 2][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][32 * ii3 + 3][i4] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][32 * ii3 + 2][i4])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][32 * ii3 + 1][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][32 * ii3 + 2][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][32 * ii3 + 2][i4])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][32 * ii3 + 2][i4 - 1]))) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][32 * ii3 + 2][i4]);
                          }
                        }
                        for (int i3 = -64 * ii0 - 32 * ii2 + 32 * ii3 + 2 * i0 + i2 - 2; i3 < _PB_N - 1; i3 += 1) {
                          for (int i4 = max(31 * ii4 + 1, 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 5); i4 <= 128 * ii0 + 64 * ii2 + 32 * ii4 - 4 * i0 - 2 * i2 + 40; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (2 * i0 + i2 == 64 * ii0 + 32 * ii2 + 5) {
                            B[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3][32 * ii4 + 31] = ((((SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3][32 * ii4 + 31] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3][32 * ii4 + 31])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][i3][32 * ii4 + 31])) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3 + 1][32 * ii4 + 31] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3][32 * ii4 + 31])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3 - 1][32 * ii4 + 31]))) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3][32 * ii4 + 32] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3][32 * ii4 + 31])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3][32 * ii4 + 30]))) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3][32 * ii4 + 31]);
                          }
                        }
                      } else {
                        for (int i3 = 32 * ii3 - 1; i3 <= 32 * ii3; i3 += 1) {
                          for (int i4 = max(1, 32 * ii3 + 32 * ii4 - i3); i4 <= 32 * ii4 + 1; i4 += 1) {
                            B[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 7][i3][i4] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3][i4])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3][i4])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3][i4])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3][i4 - 1]))) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3][i4]);
                          }
                          for (int i4 = 32 * ii4 + 2; i4 <= 32 * ii3 + 32 * ii4 - i3 + 31; i4 += 1) {
                            B[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 7][i3][i4] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3][i4])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3][i4])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3][i4])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3][i4 - 1]))) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3][i4]);
                          }
                        }
                        for (int i3 = 32 * ii3 + 1; i3 < _PB_N - 1; i3 += 1) {
                          if (32 * ii3 + 2 >= i3) {
                            for (int i4 = max(1, 32 * ii4 - 1); i4 <= 32 * ii3 + 32 * ii4 - i3 + 2; i4 += 1) {
                              B[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 7][i3][i4] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3][i4])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3][i4])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3][i4])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3][i4 - 1]))) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3][i4]);
                            }
                          } else if (ii4 >= 1) {
                            B[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3][32 * ii4 - 1] = ((((SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 7][i3][32 * ii4 - 1] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3][32 * ii4 - 1])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3][32 * ii4 - 1])) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3 + 1][32 * ii4 - 1] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3][32 * ii4 - 1])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3 - 1][32 * ii4 - 1]))) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3][32 * ii4] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3][32 * ii4 - 1])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3][32 * ii4 - 2]))) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3][32 * ii4 - 1]);
                          }
                          for (int i4 = max(max(32 * ii4, 31 * ii4 + 1), 21 * ii3 + 32 * ii4 - i3 + (ii3 + i3 - 1) / 3 + 3); i4 <= 32 * ii4 + 30; i4 += 1) {
                            B[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 7][i3][i4] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3][i4])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3][i4])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3][i4])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3][i4 - 1]))) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3][i4]);
                          }
                        }
                      }
                    }
                    for (int i2 = max(1, 64 * ii0 + 32 * ii2 - 2 * i0 + 3); i2 <= 64 * ii0 + 32 * ii2 - 2 * i0 + 34; i2 += 1) {
                      for (int i3 = 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 4; i3 <= 32 * ii3 + 1; i3 += 1) {
                        for (int i4 = max(1, 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 5); i4 <= 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 36; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                      for (int i3 = 32 * ii3 + 2; i3 < _PB_N - 1; i3 += 1) {
                        for (int i4 = max(1, 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 4); i4 <= 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 36; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                        for (int i4 = max(1, 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 37); i4 <= 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 35; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                  }
                }
              } else if (_PB_N >= 32 * ii3 + 35) {
                for (int i2 = 32 * ii2 + 1; i2 <= 32 * ii2 + 32; i2 += 1) {
                  for (int i3 = 32 * ii3 + 1; i3 <= 32 * ii3 + 32; i3 += 1) {
                    for (int i4 = 32 * ii4 + 1; i4 < _PB_N - 1; i4 += 1) {
                      A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                    }
                  }
                }
                for (int i0 = 32 * ii0 + 2; i0 <= min(min(_PB_TSTEPS, 32 * ii0 + 32), 32 * ii0 + 16 * ii2 + 17); i0 += 1) {
                  if (ii2 >= 1) {
                    for (int i2 = max(1, 64 * ii0 + 32 * ii2 - 2 * i0 + 4); i2 <= min(32 * ii2 - 1, 64 * ii0 + 32 * ii2 - 2 * i0 + 35); i2 += 1) {
                      if (2 * i0 + i2 >= 64 * ii0 + 32 * ii2 + 6) {
                        for (int i3 = max(1, 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 5); i3 <= min(32 * ii3, 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 7); i3 += 1) {
                          if (2 * i0 + i2 + i3 == 64 * ii0 + 32 * ii2 + 32 * ii3 + 7) {
                            for (int i4 = 32 * ii4 - 1; i4 <= 32 * ii4; i4 += 1) {
                              B[i2][64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 7][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 7][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 7][i4])) + A[i2 - 1][64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 7][i4])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 8][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 7][i4])) + A[i2][64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 6][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 7][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 7][i4])) + A[i2][64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 7][i4 - 1]))) + A[i2][64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 7][i4]);
                            }
                          }
                          for (int i4 = max(-64 * ii0 - 32 * ii2 - 32 * ii3 + 32 * ii4 + 2 * i0 + i2 + i3 - 6, 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 6); i4 <= 32 * ii4 + 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (ii2 == 1 && i2 == 1 && 2 * i0 + i3 >= 64 * ii0 + 32 * ii3 + 37) {
                            for (int i4 = 32 * ii4 + 2; i4 < _PB_N - 1; i4 += 1) {
                              B[1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[2][i3][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[0][i3][i4])) + (SCALAR_VAL(0.125) * ((A[1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3][i4 - 1]))) + A[1][i3][i4]);
                            }
                          } else {
                            for (int i4 = 32 * ii4 + 2; i4 < _PB_N - 1; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                        }
                      }
                      for (int i3 = max(1, 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 8); i3 <= min(min(-_PB_N - 64 * ii0 - 32 * ii2 + 32 * ii3 + 32 * ii4 + 2 * i0 + i2 - 3, 128 * ii2 + 32 * ii3 - 4 * i2 + 5), 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 36); i3 += 1) {
                        if (i3 >= 32 * ii3 + 2) {
                          for (int i4 = 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 5; i4 <= 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 7; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          for (int i4 = 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 8; i4 < _PB_N - 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        } else {
                          for (int i4 = 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 6; i4 < _PB_N - 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                      for (int i3 = 128 * ii2 + 32 * ii3 - 4 * i2 + 6; i3 <= min(-_PB_N - 64 * ii0 - 32 * ii2 + 32 * ii3 + 32 * ii4 + 2 * i0 + i2 - 3, 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 36); i3 += 1) {
                        for (int i4 = 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 5; i4 < _PB_N - 1; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                      for (int i3 = max(max(1, -64 * ii0 - 32 * ii2 + 32 * ii3 + 2 * i0 + i2 - 5), 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 5); i3 <= min(-64 * ii0 - 32 * ii2 + 32 * ii3 + 2 * i0 + i2 - 3, 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 7); i3 += 1) {
                        if (2 * i0 + i2 == 64 * ii0 + 32 * ii2 + 5 && i3 == 32 * ii3 + 2) {
                          B[64 * ii0 + 32 * ii2 - 2 * i0 + 5][32 * ii3 + 2][32 * ii4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][32 * ii3 + 2][32 * ii4] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][32 * ii3 + 2][32 * ii4])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][32 * ii3 + 2][32 * ii4])) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][32 * ii3 + 3][32 * ii4] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][32 * ii3 + 2][32 * ii4])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][32 * ii3 + 1][32 * ii4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][32 * ii3 + 2][32 * ii4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][32 * ii3 + 2][32 * ii4])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][32 * ii3 + 2][32 * ii4 - 1]))) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][32 * ii3 + 2][32 * ii4]);
                        }
                        for (int i4 = max(64 * ii0 + 32 * ii2 - 32 * ii3 + 32 * ii4 - 2 * i0 - i2 + i3 + 4, 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 6); i4 < _PB_N - 1; i4 += 1) {
                          if (64 * ii0 + 32 * ii2 + i3 + i4 + 4 >= 32 * ii3 + 32 * ii4 + 2 * i0 + i2 || 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                      if (2 * i0 + i2 == 64 * ii0 + 32 * ii2 + 4) {
                        for (int i3 = 32 * ii3 + 2; i3 <= 32 * ii3 + 3; i3 += 1) {
                          for (int i4 = 32 * ii4 + 1; i4 < _PB_N - 1; i4 += 1) {
                            B[64 * ii0 + 32 * ii2 - 2 * i0 + 4][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3][i4] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][i3][i4])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 3][i3][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][i3][i4])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][i3][i4])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][i3][i4 - 1]))) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][i3][i4]);
                          }
                        }
                      }
                      if (64 * ii0 + 32 * ii2 + 6 >= 2 * i0 + i2) {
                        for (int i3 = 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 8; i3 <= 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 36; i3 += 1) {
                          if (2 * i0 + i2 == 64 * ii0 + 32 * ii2 + 6 && i3 == 32 * ii3 + 2) {
                            B[64 * ii0 + 32 * ii2 - 2 * i0 + 6][32 * ii3 + 2][32 * ii4 - 1] = ((((SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 7][32 * ii3 + 2][32 * ii4 - 1] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][32 * ii3 + 2][32 * ii4 - 1])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][32 * ii3 + 2][32 * ii4 - 1])) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][32 * ii3 + 3][32 * ii4 - 1] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][32 * ii3 + 2][32 * ii4 - 1])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][32 * ii3 + 1][32 * ii4 - 1]))) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][32 * ii3 + 2][32 * ii4] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][32 * ii3 + 2][32 * ii4 - 1])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][32 * ii3 + 2][32 * ii4 - 2]))) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][32 * ii3 + 2][32 * ii4 - 1]);
                          }
                          for (int i4 = max(64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 5, 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 8); i4 < _PB_N - 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      } else {
                        for (int i3 = max(max(1, -_PB_N - 64 * ii0 - 32 * ii2 + 32 * ii3 + 32 * ii4 + 2 * i0 + i2 - 2), 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 8); i3 <= 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 36; i3 += 1) {
                          if (i3 >= 32 * ii3 + 2) {
                            for (int i4 = 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 5; i4 <= min(64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 7, -64 * ii0 - 32 * ii2 + 32 * ii3 + 32 * ii4 + 2 * i0 + i2 - i3 - 5); i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            if (ii2 == 1 && i2 == 1) {
                              for (int i4 = 64 * ii0 + 32 * ii4 - 2 * i0 + 39; i4 < -64 * ii0 + 32 * ii3 + 32 * ii4 + 2 * i0 - i3 - 35; i4 += 1) {
                                B[1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[2][i3][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[0][i3][i4])) + (SCALAR_VAL(0.125) * ((A[1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3][i4 - 1]))) + A[1][i3][i4]);
                              }
                            }
                            if (i2 >= 2) {
                              for (int i4 = 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 8; i4 < -64 * ii0 - 32 * ii2 + 32 * ii3 + 32 * ii4 + 2 * i0 + i2 - i3 - 4; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                          } else {
                            for (int i4 = 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 6; i4 < -64 * ii0 - 32 * ii2 - 32 * ii3 + 32 * ii4 + 2 * i0 + i2 + i3 - 6; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                          for (int i4 = -64 * ii0 - 32 * ii2 - 32 * ii3 + 32 * ii4 + 2 * i0 + i2 + i3 - 6; i4 < -64 * ii0 - 32 * ii2 + 32 * ii3 + 32 * ii4 + 2 * i0 + i2 - i3 - 4; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          for (int i4 = max(64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 5, -64 * ii0 - 32 * ii2 + 32 * ii3 + 32 * ii4 + 2 * i0 + i2 - i3 - 4); i4 < _PB_N - 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                    }
                    for (int i2 = 32 * ii2; i2 <= 64 * ii0 + 32 * ii2 - 2 * i0 + 35; i2 += 1) {
                      for (int i3 = max(1, 64 * ii0 + 32 * ii3 - 2 * i0 + 4); i3 <= min(32 * ii3 - 1, 64 * ii0 - 32 * ii2 + 32 * ii3 - 2 * i0 + i2 + 2); i3 += 1) {
                        if (2 * i0 + i3 == 64 * ii0 + 32 * ii3 + 4) {
                          for (int i4 = 32 * ii4 + 1; i4 < min(_PB_N - 1, -32 * ii2 + 32 * ii4 + i2); i4 += 1) {
                            B[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4])) + A[i2 - 1][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 5][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4])) + A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 3][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4])) + A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4 - 1]))) + A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4]);
                          }
                        } else {
                          for (int i4 = 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5; i4 < min(_PB_N - 1, -64 * ii0 - 32 * ii2 - 32 * ii3 + 32 * ii4 + 2 * i0 + i2 + i3 - 6); i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                        for (int i4 = max(-64 * ii0 - 32 * ii2 - 32 * ii3 + 32 * ii4 + 2 * i0 + i2 + i3 - 6, 64 * ii0 - 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 + i2 - i3 + 4); i4 < _PB_N - 1; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                      if (i0 >= 32 * ii0 + 3) {
                        for (int i3 = max(max(1, 64 * ii0 - 32 * ii2 + 32 * ii3 - 2 * i0 + i2 + 3), 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 5); i3 <= min(32 * ii3, _PB_N + 64 * ii0 + 32 * ii2 + 32 * ii3 - 32 * ii4 - 2 * i0 - i2 + 4); i3 += 1) {
                          if (i2 == 32 * ii2 + 2) {
                            B[32 * ii2 + 2][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 3][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5])) + A[32 * ii2 + 1][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3 + 1][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5])) + A[32 * ii2 + 2][i3 - 1][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 6] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5])) + A[32 * ii2 + 2][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 4]))) + A[32 * ii2 + 2][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5]);
                          }
                          for (int i4 = max(64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5, 64 * ii0 + 96 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - 3 * i2 - i3 + 13); i4 <= 64 * ii0 - 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 + i2 - i3 + 3; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          for (int i4 = max(64 * ii0 - 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 + i2 - i3 + 4, 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 6); i4 <= min(32 * ii4, -64 * ii0 - 32 * ii2 - 32 * ii3 + 32 * ii4 + 2 * i0 + i2 + i3 - 7); i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (i0 >= 32 * ii0 + 4 && i2 == 32 * ii2 && i3 == 32 * ii3) {
                            B[32 * ii2][32 * ii3][32 * ii4 + 1] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][32 * ii3][32 * ii4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2][32 * ii3][32 * ii4 + 1])) + A[32 * ii2 - 1][32 * ii3][32 * ii4 + 1])) + (SCALAR_VAL(0.125) * ((A[32 * ii2][32 * ii3 + 1][32 * ii4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2][32 * ii3][32 * ii4 + 1])) + A[32 * ii2][32 * ii3 - 1][32 * ii4 + 1]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2][32 * ii3][32 * ii4 + 2] - (SCALAR_VAL(2.0) * A[32 * ii2][32 * ii3][32 * ii4 + 1])) + A[32 * ii2][32 * ii3][32 * ii4]))) + A[32 * ii2][32 * ii3][32 * ii4 + 1]);
                          } else if (2 * i0 + i2 + i3 >= 64 * ii0 + 32 * ii2 + 32 * ii3 + 8 && 32 * ii3 >= i3 + 1) {
                            B[i2][i3][32 * ii4 + 1] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][32 * ii4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][32 * ii4 + 1])) + A[i2 - 1][i3][32 * ii4 + 1])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][32 * ii4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][32 * ii4 + 1])) + A[i2][i3 - 1][32 * ii4 + 1]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][32 * ii4 + 2] - (SCALAR_VAL(2.0) * A[i2][i3][32 * ii4 + 1])) + A[i2][i3][32 * ii4]))) + A[i2][i3][32 * ii4 + 1]);
                          }
                          if (i2 == 32 * ii2 && i3 == 32 * ii3) {
                            for (int i4 = 32 * ii4 + 2; i4 < -64 * ii0 + 32 * ii4 + 2 * i0 - 6; i4 += 1) {
                              B[32 * ii2][32 * ii3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][32 * ii3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][32 * ii3][i4])) + A[32 * ii2 - 1][32 * ii3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2][32 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][32 * ii3][i4])) + A[32 * ii2][32 * ii3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2][32 * ii3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2][32 * ii3][i4])) + A[32 * ii2][32 * ii3][i4 - 1]))) + A[32 * ii2][32 * ii3][i4]);
                            }
                            if (i0 == 32 * ii0 + 3) {
                              for (int i4 = 32 * ii4; i4 <= 32 * ii4 + 1; i4 += 1) {
                                B[32 * ii2][32 * ii3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][32 * ii3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][32 * ii3][i4])) + A[32 * ii2 - 1][32 * ii3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2][32 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][32 * ii3][i4])) + A[32 * ii2][32 * ii3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2][32 * ii3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2][32 * ii3][i4])) + A[32 * ii2][32 * ii3][i4 - 1]))) + A[32 * ii2][32 * ii3][i4]);
                              }
                            }
                          } else {
                            if (32 * ii3 >= i3 + 1) {
                              for (int i4 = 32 * ii4 + 2; i4 < -64 * ii0 - 32 * ii2 - 32 * ii3 + 32 * ii4 + 2 * i0 + i2 + i3 - 6; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                            if (32 * ii3 >= i3 + 1) {
                              for (int i4 = max(-64 * ii0 - 32 * ii2 - 32 * ii3 + 32 * ii4 + 2 * i0 + i2 + i3 - 6, 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 6); i4 <= 32 * ii4 + 1; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            } else {
                              B[i2][32 * ii3][32 * ii4 + 1] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii3][32 * ii4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][32 * ii4 + 1])) + A[i2 - 1][32 * ii3][32 * ii4 + 1])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 1][32 * ii4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][32 * ii4 + 1])) + A[i2][32 * ii3 - 1][32 * ii4 + 1]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3][32 * ii4 + 2] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][32 * ii4 + 1])) + A[i2][32 * ii3][32 * ii4]))) + A[i2][32 * ii3][32 * ii4 + 1]);
                              for (int i4 = 32 * ii4 + 2; i4 < -64 * ii0 - 32 * ii2 + 32 * ii4 + 2 * i0 + i2 - 6; i4 += 1) {
                                B[i2][32 * ii3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii3][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][i4])) + A[i2 - 1][32 * ii3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][i4])) + A[i2][32 * ii3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][i4])) + A[i2][32 * ii3][i4 - 1]))) + A[i2][32 * ii3][i4]);
                              }
                            }
                          }
                          for (int i4 = max(32 * ii4 + 2, -64 * ii0 - 32 * ii2 - 32 * ii3 + 32 * ii4 + 2 * i0 + i2 + i3 - 6); i4 < _PB_N - 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                      if (32 * ii2 + 1 >= i2) {
                        for (int i3 = max(1, _PB_N + 64 * ii0 + 32 * ii2 + 32 * ii3 - 32 * ii4 - 2 * i0 - i2 + 5); i3 <= 32 * ii2 + 32 * ii3 - i2; i3 += 1) {
                          for (int i4 = 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 6; i4 < _PB_N - 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                        if (ii3 >= 1 && i0 == 32 * ii0 + 2 && i2 == 32 * ii2 + 1) {
                          for (int i4 = 32 * ii4 + 1; i4 < _PB_N - 1; i4 += 1) {
                            B[32 * ii2 + 1][32 * ii3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][32 * ii3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][32 * ii3][i4])) + A[32 * ii2][32 * ii3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][32 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][32 * ii3][i4])) + A[32 * ii2 + 1][32 * ii3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][32 * ii3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][32 * ii3][i4])) + A[32 * ii2 + 1][32 * ii3][i4 - 1]))) + A[32 * ii2 + 1][32 * ii3][i4]);
                          }
                        } else if (ii3 >= 1 && 32 * ii4 + 2 * i0 >= _PB_N + 64 * ii0 + 4 && i2 == 32 * ii2 + 1) {
                          for (int i4 = 64 * ii0 + 32 * ii4 - 2 * i0 + 5; i4 < _PB_N - 1; i4 += 1) {
                            B[32 * ii2 + 1][32 * ii3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][32 * ii3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][32 * ii3][i4])) + A[32 * ii2][32 * ii3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][32 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][32 * ii3][i4])) + A[32 * ii2 + 1][32 * ii3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][32 * ii3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][32 * ii3][i4])) + A[32 * ii2 + 1][32 * ii3][i4 - 1]))) + A[32 * ii2 + 1][32 * ii3][i4]);
                          }
                        }
                        for (int i3 = 32 * ii3 + 1; i3 <= 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 36; i3 += 1) {
                          for (int i4 = 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 5; i4 <= min(-64 * ii0 - 32 * ii2 + 32 * ii3 + 32 * ii4 + 2 * i0 + i2 - i3 - 5, 64 * ii0 + 6 * ii2 - 6 * ii3 + 32 * ii4 - 2 * i0 + floord(2 * ii2 - 2 * ii3 - i2 + i3 - 1, 5) + 4); i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          for (int i4 = 64 * ii0 + 6 * ii2 - 6 * ii3 + 32 * ii4 - 2 * i0 + floord(2 * ii2 - 2 * ii3 - i2 + i3 - 1, 5) + 5; i4 < min(_PB_N - 1, -64 * ii0 - 32 * ii2 + 32 * ii3 + 32 * ii4 + 2 * i0 + i2 - i3 - 4); i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (i0 == 32 * ii0 + 3 && i2 == 32 * ii2) {
                            for (int i4 = max(32 * ii4 - 1, 32 * ii3 + 32 * ii4 - i3 + 2); i4 <= 32 * ii4; i4 += 1) {
                              B[32 * ii2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3][i4 - 1]))) + A[32 * ii2][i3][i4]);
                            }
                          }
                          if (64 * ii0 + 32 * ii2 + 6 >= 2 * i0 + i2) {
                            for (int i4 = max(32 * ii4, 64 * ii2 + 32 * ii4 - 2 * i2 + 1); i4 < _PB_N - 1; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          } else {
                            for (int i4 = max(64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 5, -64 * ii0 - 32 * ii2 + 32 * ii3 + 32 * ii4 + 2 * i0 + i2 - i3 - 4); i4 < _PB_N - 1; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                        }
                      } else {
                        for (int i3 = max(max(1, 64 * ii0 - 32 * ii2 + 32 * ii3 - 2 * i0 + i2 + 3), _PB_N + 64 * ii0 + 32 * ii2 + 32 * ii3 - 32 * ii4 - 2 * i0 - i2 + 5); i3 < 32 * ii3; i3 += 1) {
                          for (int i4 = max(64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5, 64 * ii0 + 96 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - 3 * i2 - i3 + 13); i4 <= 64 * ii0 - 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 + i2 - i3 + 3; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (i2 == 32 * ii2 + 2) {
                            B[32 * ii2 + 2][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 3][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5])) + A[32 * ii2 + 1][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3 + 1][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5])) + A[32 * ii2 + 2][i3 - 1][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 6] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5])) + A[32 * ii2 + 2][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 4]))) + A[32 * ii2 + 2][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5]);
                          }
                          for (int i4 = 64 * ii0 - 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 + i2 - i3 + 4; i4 < _PB_N - 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                        if (ii3 >= 1 && 32 * ii4 + 2 * i0 + i2 >= _PB_N + 64 * ii0 + 32 * ii2 + 5 && 32 * ii2 + 2 * i0 >= 64 * ii0 + i2 + 3) {
                          for (int i4 = max(64 * ii0 + 32 * ii4 - 2 * i0 + 5, 64 * ii0 + 96 * ii2 + 32 * ii4 - 2 * i0 - 3 * i2 + 13); i4 <= 64 * ii0 - 32 * ii2 + 32 * ii4 - 2 * i0 + i2 + 3; i4 += 1) {
                            B[i2][32 * ii3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii3][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][i4])) + A[i2 - 1][32 * ii3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][i4])) + A[i2][32 * ii3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][i4])) + A[i2][32 * ii3][i4 - 1]))) + A[i2][32 * ii3][i4]);
                          }
                          if (i2 == 32 * ii2 + 2) {
                            B[32 * ii2 + 2][32 * ii3][64 * ii0 + 32 * ii4 - 2 * i0 + 5] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 3][32 * ii3][64 * ii0 + 32 * ii4 - 2 * i0 + 5] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][32 * ii3][64 * ii0 + 32 * ii4 - 2 * i0 + 5])) + A[32 * ii2 + 1][32 * ii3][64 * ii0 + 32 * ii4 - 2 * i0 + 5])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][32 * ii3 + 1][64 * ii0 + 32 * ii4 - 2 * i0 + 5] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][32 * ii3][64 * ii0 + 32 * ii4 - 2 * i0 + 5])) + A[32 * ii2 + 2][32 * ii3 - 1][64 * ii0 + 32 * ii4 - 2 * i0 + 5]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][32 * ii3][64 * ii0 + 32 * ii4 - 2 * i0 + 6] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][32 * ii3][64 * ii0 + 32 * ii4 - 2 * i0 + 5])) + A[32 * ii2 + 2][32 * ii3][64 * ii0 + 32 * ii4 - 2 * i0 + 4]))) + A[32 * ii2 + 2][32 * ii3][64 * ii0 + 32 * ii4 - 2 * i0 + 5]);
                          }
                          for (int i4 = 64 * ii0 - 32 * ii2 + 32 * ii4 - 2 * i0 + i2 + 4; i4 < _PB_N - 1; i4 += 1) {
                            B[i2][32 * ii3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii3][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][i4])) + A[i2 - 1][32 * ii3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][i4])) + A[i2][32 * ii3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][i4])) + A[i2][32 * ii3][i4 - 1]))) + A[i2][32 * ii3][i4]);
                          }
                        } else if (ii3 >= 1 && 64 * ii0 + i2 + 2 >= 32 * ii2 + 2 * i0) {
                          if (i0 == 32 * ii0 + 2) {
                            for (int i4 = 32 * ii4 + 1; i4 < min(_PB_N - 1, -32 * ii2 + 32 * ii4 + i2); i4 += 1) {
                              B[i2][32 * ii3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii3][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][i4])) + A[i2 - 1][32 * ii3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][i4])) + A[i2][32 * ii3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][i4])) + A[i2][32 * ii3][i4 - 1]))) + A[i2][32 * ii3][i4]);
                            }
                          } else {
                            for (int i4 = 64 * ii0 + 32 * ii4 - 2 * i0 + 5; i4 <= min(_PB_N - 2, 64 * ii0 - 32 * ii2 + 32 * ii4 - 2 * i0 + i2 + 3); i4 += 1) {
                              B[i2][32 * ii3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii3][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][i4])) + A[i2 - 1][32 * ii3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][i4])) + A[i2][32 * ii3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][i4])) + A[i2][32 * ii3][i4 - 1]))) + A[i2][32 * ii3][i4]);
                            }
                          }
                          for (int i4 = 64 * ii0 - 32 * ii2 + 32 * ii4 - 2 * i0 + i2 + 4; i4 < _PB_N - 1; i4 += 1) {
                            B[i2][32 * ii3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii3][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][i4])) + A[i2 - 1][32 * ii3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][i4])) + A[i2][32 * ii3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][i4])) + A[i2][32 * ii3][i4 - 1]))) + A[i2][32 * ii3][i4]);
                          }
                        }
                        for (int i3 = 32 * ii3 + 1; i3 <= min(min(-128 * ii0 - 32 * ii2 + 32 * ii3 + 4 * i0 + i2 - 9, -32 * ii2 + 32 * ii3 + i2), 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 36); i3 += 1) {
                          if (i0 >= 32 * ii0 + 3) {
                            for (int i4 = 64 * ii0 + 32 * ii4 - 2 * i0 + 4; i4 < min(_PB_N - 1, -64 * ii0 - 32 * ii2 + 32 * ii3 + 32 * ii4 + 2 * i0 + i2 - i3 - 4); i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            for (int i4 = -64 * ii0 - 32 * ii2 + 32 * ii3 + 32 * ii4 + 2 * i0 + i2 - i3 - 4; i4 < _PB_N - 1; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          } else {
                            for (int i4 = 32 * ii4; i4 < _PB_N - 1; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                        }
                        for (int i3 = -32 * ii2 + 32 * ii3 + i2 + 1; i3 <= min(-128 * ii0 - 32 * ii2 + 32 * ii3 + 4 * i0 + i2 - 9, 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 36); i3 += 1) {
                          for (int i4 = 64 * ii0 + 32 * ii4 - 2 * i0 + 4; i4 <= min(-64 * ii0 - 32 * ii2 + 32 * ii3 + 32 * ii4 + 2 * i0 + i2 - i3 - 5, 64 * ii0 + 6 * ii2 - 6 * ii3 + 32 * ii4 - 2 * i0 + floord(2 * ii2 - 2 * ii3 - i2 + i3 - 1, 5) + 4); i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (i2 == 32 * ii2 + 2 && 32 * ii3 + 4 * i0 >= 128 * ii0 + i3 + 8 && 32 * ii3 + 7 >= i3) {
                            B[32 * ii2 + 2][i3][64 * ii0 + 32 * ii4 - 2 * i0 + 5] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 3][i3][64 * ii0 + 32 * ii4 - 2 * i0 + 5] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][i3][64 * ii0 + 32 * ii4 - 2 * i0 + 5])) + A[32 * ii2 + 1][i3][64 * ii0 + 32 * ii4 - 2 * i0 + 5])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3 + 1][64 * ii0 + 32 * ii4 - 2 * i0 + 5] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][i3][64 * ii0 + 32 * ii4 - 2 * i0 + 5])) + A[32 * ii2 + 2][i3 - 1][64 * ii0 + 32 * ii4 - 2 * i0 + 5]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3][64 * ii0 + 32 * ii4 - 2 * i0 + 6] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][i3][64 * ii0 + 32 * ii4 - 2 * i0 + 5])) + A[32 * ii2 + 2][i3][64 * ii0 + 32 * ii4 - 2 * i0 + 4]))) + A[32 * ii2 + 2][i3][64 * ii0 + 32 * ii4 - 2 * i0 + 5]);
                          }
                          for (int i4 = max(64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 8, 64 * ii0 + 6 * ii2 - 6 * ii3 + 32 * ii4 - 2 * i0 + floord(2 * ii2 - 2 * ii3 - i2 + i3 - 1, 5) + 5); i4 < min(32 * ii4, -64 * ii0 - 32 * ii2 + 32 * ii3 + 32 * ii4 + 2 * i0 + i2 - i3 - 4); i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          for (int i4 = 32 * ii4; i4 < min(_PB_N - 1, -64 * ii0 - 32 * ii2 + 32 * ii3 + 32 * ii4 + 2 * i0 + i2 - i3 - 4); i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          for (int i4 = -64 * ii0 - 32 * ii2 + 32 * ii3 + 32 * ii4 + 2 * i0 + i2 - i3 - 4; i4 < _PB_N - 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                        for (int i3 = -128 * ii0 - 32 * ii2 + 32 * ii3 + 4 * i0 + i2 - 8; i3 <= 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 36; i3 += 1) {
                          if (i0 >= 32 * ii0 + 3) {
                            for (int i4 = 64 * ii0 + 32 * ii4 - 2 * i0 + 4; i4 < _PB_N - 1; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          } else {
                            for (int i4 = 32 * ii4; i4 < _PB_N - 1; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                        }
                      }
                      for (int i3 = 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 37; i3 <= 64 * ii0 + 32 * ii3 - 2 * i0 + 35; i3 += 1) {
                        for (int i4 = 64 * ii0 + 32 * ii4 - 2 * i0 + 4; i4 <= min(-64 * ii0 - 32 * ii2 + 32 * ii3 + 32 * ii4 + 2 * i0 + i2 - i3 - 5, 64 * ii0 + 6 * ii2 - 6 * ii3 + 32 * ii4 - 2 * i0 + floord(2 * ii2 - 2 * ii3 - i2 + i3 - 1, 5) + 4); i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = max(64 * ii0 + 32 * ii4 - 2 * i0 + 4, 64 * ii0 + 6 * ii2 - 6 * ii3 + 32 * ii4 - 2 * i0 + floord(2 * ii2 - 2 * ii3 - i2 + i3 - 1, 5) + 5); i4 < min(32 * ii4, -64 * ii0 - 32 * ii2 + 32 * ii3 + 32 * ii4 + 2 * i0 + i2 - i3 - 4); i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = max(64 * ii0 + 32 * ii4 - 2 * i0 + 4, -64 * ii0 - 32 * ii2 + 32 * ii3 + 32 * ii4 + 2 * i0 + i2 - i3 - 4); i4 < 32 * ii4; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = 32 * ii4; i4 < _PB_N - 1; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    }
                  } else if (32 * ii4 + 4 * i0 >= _PB_N + 128 * ii0 + 37) {
                    for (int i2 = 1; i2 <= 64 * ii0 - 2 * i0 + 35; i2 += 1) {
                      for (int i3 = max(1, 64 * ii0 + 32 * ii3 - 2 * i0 + 4); i3 <= 64 * ii0 + 32 * ii3 - 2 * i0 + 35; i3 += 1) {
                        if (i3 >= 32 * ii3 + 2) {
                          for (int i4 = 64 * ii0 + 32 * ii4 - 2 * i0 + 4; i4 <= 64 * ii0 - 6 * ii3 + 32 * ii4 - 2 * i0 + (-2 * ii3 - i2 + i3 + 14) / 5 + 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          for (int i4 = max(64 * ii0 + 32 * ii4 - 2 * i0 + 4, 64 * ii0 - 6 * ii3 + 32 * ii4 - 2 * i0 + (-2 * ii3 - i2 + i3 + 14) / 5 + 2); i4 < _PB_N - 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        } else {
                          if (i3 == 32 * ii3 + 1) {
                            for (int i4 = 64 * ii0 + 32 * ii4 - 2 * i0 + 4; i4 <= min(64 * ii0 + 32 * ii4 - 2 * i0 + 5, 64 * ii0 + 32 * ii4 - 2 * i0 + i2 + 2); i4 += 1) {
                              B[i2][32 * ii3 + 1][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 1][i4])) + A[i2 - 1][32 * ii3 + 1][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 2][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 1][i4])) + A[i2][32 * ii3][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 1][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 1][i4])) + A[i2][32 * ii3 + 1][i4 - 1]))) + A[i2][32 * ii3 + 1][i4]);
                            }
                          }
                          if (2 * i0 + i3 >= 64 * ii0 + 32 * ii3 + 5) {
                            for (int i4 = max(max(64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5, 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - 3 * i2 - i3 + 13), 64 * ii0 - 6 * ii3 + 32 * ii4 - 2 * i0 + (-2 * ii3 + i3 - 1) / 5 + 6); i4 <= min(_PB_N - 2, 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 + i2 - i3 + 3); i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            if (i2 == 2 && 32 * ii3 >= i3) {
                              B[2][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5] = ((((SCALAR_VAL(0.125) * ((A[3][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5] - (SCALAR_VAL(2.0) * A[2][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5])) + A[1][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5])) + (SCALAR_VAL(0.125) * ((A[2][i3 + 1][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5] - (SCALAR_VAL(2.0) * A[2][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5])) + A[2][i3 - 1][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5]))) + (SCALAR_VAL(0.125) * ((A[2][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 6] - (SCALAR_VAL(2.0) * A[2][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5])) + A[2][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 4]))) + A[2][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5]);
                            }
                          }
                          for (int i4 = 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 + i2 - i3 + 4; i4 <= min(min(_PB_N - 2, -64 * ii0 - 32 * ii3 + 32 * ii4 + 2 * i0 + i2 + i3 - 7), 16 * ii3 + 32 * ii4 - i3 + i3 / 2); i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (32 * ii3 >= i3 + 1) {
                            for (int i4 = 16 * ii3 + 32 * ii4 - i3 + i3 / 2 + 1; i4 < min(_PB_N - 1, -64 * ii0 - 32 * ii3 + 32 * ii4 + 2 * i0 + i2 + i3 - 6); i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                          for (int i4 = max(-64 * ii0 - 32 * ii3 + 32 * ii4 + 2 * i0 + i2 + i3 - 6, 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 + i2 - i3 + 4); i4 <= 32 * ii4 + 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (i2 == 1 && 2 * i0 + i3 >= 64 * ii0 + 32 * ii3 + 5 && 64 * ii0 + 32 * ii3 + 6 >= 2 * i0 + i3) {
                            for (int i4 = 32 * ii4 + 2; i4 < _PB_N - 1; i4 += 1) {
                              B[1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[2][i3][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[0][i3][i4])) + (SCALAR_VAL(0.125) * ((A[1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3][i4 - 1]))) + A[1][i3][i4]);
                            }
                          } else {
                            if (i2 <= 2 && 2 * i0 + i3 == 64 * ii0 + 32 * ii3 + i2 + 3) {
                              for (int i4 = 32 * ii4 + 2; i4 < _PB_N - 1; i4 += 1) {
                                B[i2][64 * ii0 + 32 * ii3 - 2 * i0 + i2 + 3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + i2 + 3][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + i2 + 3][i4])) + A[i2 - 1][64 * ii0 + 32 * ii3 - 2 * i0 + i2 + 3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + i2 + 4][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + i2 + 3][i4])) + A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + i2 + 2][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + i2 + 3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + i2 + 3][i4])) + A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + i2 + 3][i4 - 1]))) + A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + i2 + 3][i4]);
                              }
                            } else {
                              if (i3 >= 32 * ii3) {
                                for (int i4 = 32 * ii3 + 32 * ii4 - i3 + 1; i4 < min(_PB_N - 1, -64 * ii0 - 32 * ii3 + 32 * ii4 + 2 * i0 + i2 + i3 - 6); i4 += 1) {
                                  B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                                }
                              }
                              if (2 * i0 + i2 + i3 >= 64 * ii0 + 32 * ii3 + 8 && 2 * i0 + i3 >= 64 * ii0 + 32 * ii3 + i2 + 3) {
                                for (int i4 = -64 * ii0 - 32 * ii3 + 32 * ii4 + 2 * i0 + i2 + i3 - 6; i4 < _PB_N - 1; i4 += 1) {
                                  B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                                }
                              }
                            }
                            if (2 * i0 + i3 == 64 * ii0 + 32 * ii3 + 4) {
                              for (int i4 = 32 * ii4 + 1; i4 < min(_PB_N - 1, 32 * ii4 + i2); i4 += 1) {
                                B[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4])) + A[i2 - 1][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 5][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4])) + A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 3][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4])) + A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4 - 1]))) + A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4]);
                              }
                            }
                            if (64 * ii0 + 32 * ii3 + i2 + 2 >= 2 * i0 + i3) {
                              for (int i4 = max(-64 * ii0 - 32 * ii3 + 32 * ii4 + 2 * i0 + i2 + i3 - 6, 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 + i2 - i3 + 4); i4 < _PB_N - 1; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                          }
                        }
                      }
                    }
                  } else if (i0 >= 32 * ii0 + 3) {
                    for (int i3 = max(1, 64 * ii0 + 32 * ii3 - 2 * i0 + 4); i3 < 32 * ii3; i3 += 1) {
                      for (int i4 = 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5; i4 < _PB_N - 1; i4 += 1) {
                        B[1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[2][i3][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[0][i3][i4])) + (SCALAR_VAL(0.125) * ((A[1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3][i4 - 1]))) + A[1][i3][i4]);
                      }
                    }
                    for (int i3 = max(1, 32 * ii3); i3 <= min(-128 * ii0 + 32 * ii3 + 4 * i0 - 9, 64 * ii0 + 32 * ii3 - 2 * i0 + 35); i3 += 1) {
                      if (32 * ii3 + 1 >= i3) {
                        for (int i4 = 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5; i4 < min(_PB_N - 1, -64 * ii0 - 32 * ii3 + 32 * ii4 + 2 * i0 + i3 - 5); i4 += 1) {
                          B[1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[2][i3][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[0][i3][i4])) + (SCALAR_VAL(0.125) * ((A[1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3][i4 - 1]))) + A[1][i3][i4]);
                        }
                        if (i0 >= 32 * ii0 + 4 && i3 == 32 * ii3) {
                          for (int i4 = -64 * ii0 + 32 * ii4 + 2 * i0 - 5; i4 < min(_PB_N - 1, -64 * ii0 + 32 * ii4 + 2 * i0 - 3); i4 += 1) {
                            B[1][32 * ii3][i4] = ((((SCALAR_VAL(0.125) * ((A[2][32 * ii3][i4] - (SCALAR_VAL(2.0) * A[1][32 * ii3][i4])) + A[0][32 * ii3][i4])) + (SCALAR_VAL(0.125) * ((A[1][32 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[1][32 * ii3][i4])) + A[1][32 * ii3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[1][32 * ii3][i4 + 1] - (SCALAR_VAL(2.0) * A[1][32 * ii3][i4])) + A[1][32 * ii3][i4 - 1]))) + A[1][32 * ii3][i4]);
                          }
                        } else if (i0 == 32 * ii0 + 3 && i3 == 32 * ii3) {
                          for (int i4 = 32 * ii4 + 1; i4 < _PB_N - 1; i4 += 1) {
                            B[1][32 * ii3][i4] = ((((SCALAR_VAL(0.125) * ((A[2][32 * ii3][i4] - (SCALAR_VAL(2.0) * A[1][32 * ii3][i4])) + A[0][32 * ii3][i4])) + (SCALAR_VAL(0.125) * ((A[1][32 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[1][32 * ii3][i4])) + A[1][32 * ii3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[1][32 * ii3][i4 + 1] - (SCALAR_VAL(2.0) * A[1][32 * ii3][i4])) + A[1][32 * ii3][i4 - 1]))) + A[1][32 * ii3][i4]);
                          }
                        }
                      } else {
                        for (int i4 = 64 * ii0 + 32 * ii4 - 2 * i0 + 4; i4 < min(_PB_N - 1, -64 * ii0 + 32 * ii3 + 32 * ii4 + 2 * i0 - i3 - 3); i4 += 1) {
                          B[1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[2][i3][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[0][i3][i4])) + (SCALAR_VAL(0.125) * ((A[1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3][i4 - 1]))) + A[1][i3][i4]);
                        }
                      }
                      if (2 * i0 + i3 >= 64 * ii0 + 32 * ii3 + 7) {
                        for (int i4 = -64 * ii0 + 32 * ii3 + 32 * ii4 + 2 * i0 - i3 - 3; i4 < _PB_N - 1; i4 += 1) {
                          B[1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[2][i3][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[0][i3][i4])) + (SCALAR_VAL(0.125) * ((A[1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3][i4 - 1]))) + A[1][i3][i4]);
                        }
                      }
                    }
                    for (int i3 = -128 * ii0 + 32 * ii3 + 4 * i0 - 8; i3 <= 64 * ii0 + 32 * ii3 - 2 * i0 + 35; i3 += 1) {
                      if (128 * ii0 + i3 + 8 == 32 * ii3 + 4 * i0) {
                        B[1][-128 * ii0 + 32 * ii3 + 4 * i0 - 8][64 * ii0 + 32 * ii4 - 2 * i0 + 4] = ((((SCALAR_VAL(0.125) * ((A[2][-128 * ii0 + 32 * ii3 + 4 * i0 - 8][64 * ii0 + 32 * ii4 - 2 * i0 + 4] - (SCALAR_VAL(2.0) * A[1][-128 * ii0 + 32 * ii3 + 4 * i0 - 8][64 * ii0 + 32 * ii4 - 2 * i0 + 4])) + A[0][-128 * ii0 + 32 * ii3 + 4 * i0 - 8][64 * ii0 + 32 * ii4 - 2 * i0 + 4])) + (SCALAR_VAL(0.125) * ((A[1][-128 * ii0 + 32 * ii3 + 4 * i0 - 7][64 * ii0 + 32 * ii4 - 2 * i0 + 4] - (SCALAR_VAL(2.0) * A[1][-128 * ii0 + 32 * ii3 + 4 * i0 - 8][64 * ii0 + 32 * ii4 - 2 * i0 + 4])) + A[1][-128 * ii0 + 32 * ii3 + 4 * i0 - 9][64 * ii0 + 32 * ii4 - 2 * i0 + 4]))) + (SCALAR_VAL(0.125) * ((A[1][-128 * ii0 + 32 * ii3 + 4 * i0 - 8][64 * ii0 + 32 * ii4 - 2 * i0 + 5] - (SCALAR_VAL(2.0) * A[1][-128 * ii0 + 32 * ii3 + 4 * i0 - 8][64 * ii0 + 32 * ii4 - 2 * i0 + 4])) + A[1][-128 * ii0 + 32 * ii3 + 4 * i0 - 8][64 * ii0 + 32 * ii4 - 2 * i0 + 3]))) + A[1][-128 * ii0 + 32 * ii3 + 4 * i0 - 8][64 * ii0 + 32 * ii4 - 2 * i0 + 4]);
                      }
                      for (int i4 = max(64 * ii0 + 32 * ii4 - 2 * i0 + 4, -64 * ii0 + 32 * ii3 + 32 * ii4 + 2 * i0 - i3 - 3); i4 < _PB_N - 1; i4 += 1) {
                        B[1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[2][i3][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[0][i3][i4])) + (SCALAR_VAL(0.125) * ((A[1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3][i4 - 1]))) + A[1][i3][i4]);
                      }
                    }
                    for (int i2 = 2; i2 <= 64 * ii0 - 2 * i0 + 35; i2 += 1) {
                      for (int i3 = max(1, 64 * ii0 + 32 * ii3 - 2 * i0 + 4); i3 <= 64 * ii0 + 32 * ii3 - 2 * i0 + 35; i3 += 1) {
                        if (i3 == 32 * ii3 + 1) {
                          for (int i4 = 64 * ii0 + 32 * ii4 - 2 * i0 + 4; i4 <= min(64 * ii0 + 32 * ii4 - 2 * i0 + 5, 64 * ii0 + 32 * ii4 - 2 * i0 + i2 + 2); i4 += 1) {
                            B[i2][32 * ii3 + 1][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 1][i4])) + A[i2 - 1][32 * ii3 + 1][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 2][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 1][i4])) + A[i2][32 * ii3][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 1][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 1][i4])) + A[i2][32 * ii3 + 1][i4 - 1]))) + A[i2][32 * ii3 + 1][i4]);
                          }
                        }
                        if (2 * i0 + i3 >= 64 * ii0 + 32 * ii3 + 5 && 32 * ii3 + 1 >= i3) {
                          for (int i4 = max(max(64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5, 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - 3 * i2 - i3 + 13), 64 * ii0 - 6 * ii3 + 32 * ii4 - 2 * i0 + (-2 * ii3 + i3 - 1) / 5 + 6); i4 <= min(min(_PB_N - 2, 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 + i2 - i3 + 3), 16 * ii3 + 32 * ii4 - i3 + i3 / 2); i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (32 * ii3 >= i3 + 1) {
                            for (int i4 = 16 * ii3 + 32 * ii4 - i3 + i3 / 2 + 1; i4 <= min(_PB_N - 2, 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 + i2 - i3 + 3); i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                          if (i2 == 2 && 32 * ii3 >= i3) {
                            B[2][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5] = ((((SCALAR_VAL(0.125) * ((A[3][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5] - (SCALAR_VAL(2.0) * A[2][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5])) + A[1][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5])) + (SCALAR_VAL(0.125) * ((A[2][i3 + 1][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5] - (SCALAR_VAL(2.0) * A[2][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5])) + A[2][i3 - 1][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5]))) + (SCALAR_VAL(0.125) * ((A[2][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 6] - (SCALAR_VAL(2.0) * A[2][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5])) + A[2][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 4]))) + A[2][i3][64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5]);
                          }
                        }
                        if (32 * ii3 + 1 >= i3) {
                          for (int i4 = 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 + i2 - i3 + 4; i4 <= min(min(_PB_N - 2, -64 * ii0 - 32 * ii3 + 32 * ii4 + 2 * i0 + i2 + i3 - 7), 16 * ii3 + 32 * ii4 - i3 + i3 / 2); i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (32 * ii3 >= i3 + 1) {
                            for (int i4 = max(64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 + i2 - i3 + 4, 16 * ii3 + 32 * ii4 - i3 + i3 / 2 + 1); i4 < min(_PB_N - 1, -64 * ii0 - 32 * ii3 + 32 * ii4 + 2 * i0 + i2 + i3 - 6); i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            if (2 * i0 + i3 == 64 * ii0 + 32 * ii3 + 4) {
                              for (int i4 = 32 * ii4 + 1; i4 < min(_PB_N - 1, 32 * ii4 + i2); i4 += 1) {
                                B[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4])) + A[i2 - 1][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 5][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4])) + A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 3][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4])) + A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4 - 1]))) + A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4]);
                              }
                            }
                            if (64 * ii0 + 32 * ii3 + i2 + 2 >= 2 * i0 + i3) {
                              for (int i4 = max(-64 * ii0 - 32 * ii3 + 32 * ii4 + 2 * i0 + i2 + i3 - 6, 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 + i2 - i3 + 4); i4 < min(_PB_N - 1, -64 * ii0 + 32 * ii3 + 32 * ii4 + 2 * i0 + i2 - i3 - 4); i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                              if (2 * i0 + i3 == 64 * ii0 + 32 * ii3 + 4) {
                                for (int i4 = -128 * ii0 + 32 * ii4 + 4 * i0 + i2 - 8; i4 < _PB_N - 1; i4 += 1) {
                                  B[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4])) + A[i2 - 1][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 5][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4])) + A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 3][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4])) + A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4 - 1]))) + A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4]);
                                }
                              }
                            } else if (2 * i0 + i2 + i3 >= 64 * ii0 + 32 * ii3 + 8) {
                              for (int i4 = -64 * ii0 - 32 * ii3 + 32 * ii4 + 2 * i0 + i2 + i3 - 6; i4 < min(_PB_N - 1, -64 * ii0 + 32 * ii3 + 32 * ii4 + 2 * i0 + i2 - i3 - 4); i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                            if (2 * i0 + i2 + i3 >= 64 * ii0 + 32 * ii3 + 8 && 2 * i0 + i3 >= 64 * ii0 + 32 * ii3 + 5) {
                              for (int i4 = -64 * ii0 + 32 * ii3 + 32 * ii4 + 2 * i0 + i2 - i3 - 4; i4 < _PB_N - 1; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                          }
                        } else {
                          for (int i4 = 64 * ii0 + 32 * ii4 - 2 * i0 + 4; i4 < min(32 * ii4, -64 * ii0 + 32 * ii3 + 32 * ii4 + 2 * i0 + i2 - i3 - 4); i4 += 1) {
                            if (320 * ii0 + 160 * ii4 + i3 + 19 >= 32 * ii3 + 10 * i0 + i2 + 5 * i4 || 2 * i0 + i2 + i4 >= 64 * ii0 + 32 * ii4 + 8 || 64 * ii0 + 32 * ii3 + 36 >= 2 * i0 + i2 + i3 || 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                        }
                        for (int i4 = max(64 * ii0 + 32 * ii4 - 2 * i0 + 4, -64 * ii0 + 32 * ii3 + 32 * ii4 + 2 * i0 + i2 - i3 - 4); i4 < 32 * ii4; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        if (i2 == 2 && 2 * i0 + i3 == 64 * ii0 + 32 * ii3 + 5) {
                          for (int i4 = 32 * ii4 + 1; i4 < _PB_N - 1; i4 += 1) {
                            B[2][64 * ii0 + 32 * ii3 - 2 * i0 + 5][i4] = ((((SCALAR_VAL(0.125) * ((A[3][64 * ii0 + 32 * ii3 - 2 * i0 + 5][i4] - (SCALAR_VAL(2.0) * A[2][64 * ii0 + 32 * ii3 - 2 * i0 + 5][i4])) + A[1][64 * ii0 + 32 * ii3 - 2 * i0 + 5][i4])) + (SCALAR_VAL(0.125) * ((A[2][64 * ii0 + 32 * ii3 - 2 * i0 + 6][i4] - (SCALAR_VAL(2.0) * A[2][64 * ii0 + 32 * ii3 - 2 * i0 + 5][i4])) + A[2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4]))) + (SCALAR_VAL(0.125) * ((A[2][64 * ii0 + 32 * ii3 - 2 * i0 + 5][i4 + 1] - (SCALAR_VAL(2.0) * A[2][64 * ii0 + 32 * ii3 - 2 * i0 + 5][i4])) + A[2][64 * ii0 + 32 * ii3 - 2 * i0 + 5][i4 - 1]))) + A[2][64 * ii0 + 32 * ii3 - 2 * i0 + 5][i4]);
                          }
                        } else if (i3 >= 32 * ii3) {
                          for (int i4 = max(32 * ii4, 16 * ii3 + 32 * ii4 - i3 + i3 / 2 + 1); i4 < _PB_N - 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                    }
                  } else {
                    for (int i2 = 1; i2 <= 31; i2 += 1) {
                      for (int i3 = max(1, 32 * ii3); i3 <= 32 * ii3 + 31; i3 += 1) {
                        for (int i4 = max(32 * ii4, 16 * ii3 + 32 * ii4 - i3 + i3 / 2 + 1); i4 < _PB_N - 1; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    }
                  }
                  for (int i2 = max(1, 64 * ii0 + 32 * ii2 - 2 * i0 + 3); i2 <= 64 * ii0 + 32 * ii2 - 2 * i0 + 34; i2 += 1) {
                    if (i2 >= 32 * ii2 + 1) {
                      for (int i3 = max(1, 64 * ii0 + 32 * ii3 - 2 * i0 + 3); i3 <= 64 * ii0 + 32 * ii3 - 2 * i0 + 34; i3 += 1) {
                        for (int i4 = max(64 * ii0 + 32 * ii4 - 2 * i0 + 3, 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 4); i4 < _PB_N - 1; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    } else {
                      for (int i3 = max(1, 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 4); i3 <= 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 35; i3 += 1) {
                        for (int i4 = max(64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 4, 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 5); i4 < _PB_N - 1; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                  }
                }
              } else {
                for (int i0 = 32 * ii0 + 1; i0 <= min(min(_PB_TSTEPS, 32 * ii0 + 32), 32 * ii0 + 16 * ii2 + 17); i0 += 1) {
                  if (i0 >= 32 * ii0 + 2) {
                    if (i0 >= 32 * ii0 + 3) {
                      for (int i2 = max(1, 64 * ii0 + 32 * ii2 - 2 * i0 + 4); i2 <= 64 * ii0 + 32 * ii2 - 2 * i0 + 35; i2 += 1) {
                        if (_PB_N + 2 * i0 + i2 >= 64 * ii0 + 32 * ii2 + 32 * ii3 + 8) {
                          if (i2 == 32 * ii2 + 1) {
                            for (int i4 = 32 * ii3 + 1; i4 < _PB_N - 1; i4 += 1) {
                              B[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4])) + A[32 * ii2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 5][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4])) + A[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 3][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4])) + A[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4 - 1]))) + A[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4]);
                            }
                          }
                          for (int i3 = max(64 * ii0 + 32 * ii3 - 2 * i0 + 5, 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 5); i3 <= -_PB_N + 64 * ii0 + 32 * ii2 + 64 * ii3 - 2 * i0 - i2 + 9; i3 += 1) {
                            if (_PB_N + 2 * i0 + i2 == 64 * ii0 + 32 * ii2 + 32 * ii3 + 8 && i3 == 32 * ii3 + 1) {
                              B[-_PB_N + 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 + 8][32 * ii3 + 1][_PB_N - 3] = ((((SCALAR_VAL(0.125) * ((A[-_PB_N + 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 + 9][32 * ii3 + 1][_PB_N - 3] - (SCALAR_VAL(2.0) * A[-_PB_N + 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 + 8][32 * ii3 + 1][_PB_N - 3])) + A[-_PB_N + 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 + 7][32 * ii3 + 1][_PB_N - 3])) + (SCALAR_VAL(0.125) * ((A[-_PB_N + 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 + 8][32 * ii3 + 2][_PB_N - 3] - (SCALAR_VAL(2.0) * A[-_PB_N + 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 + 8][32 * ii3 + 1][_PB_N - 3])) + A[-_PB_N + 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 + 8][32 * ii3][_PB_N - 3]))) + (SCALAR_VAL(0.125) * ((A[-_PB_N + 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 + 8][32 * ii3 + 1][_PB_N - 2] - (SCALAR_VAL(2.0) * A[-_PB_N + 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 + 8][32 * ii3 + 1][_PB_N - 3])) + A[-_PB_N + 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 + 8][32 * ii3 + 1][_PB_N - 4]))) + A[-_PB_N + 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 + 8][32 * ii3 + 1][_PB_N - 3]);
                            }
                            for (int i4 = max(64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 6, 64 * ii0 + 32 * ii2 + 64 * ii3 - 2 * i0 - i2 - i3 + 6); i4 < _PB_N - 1; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                          if (32 * ii3 + 4 == _PB_N && 2 * i0 + i2 == 64 * ii0 + 32 * ii2 + 4) {
                            for (int i4 = _PB_N - 3; i4 < _PB_N - 1; i4 += 1) {
                              B[64 * ii0 + 32 * ii2 - 2 * i0 + 4][_PB_N - 2][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][_PB_N - 2][i4] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][_PB_N - 2][i4])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 3][_PB_N - 2][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][_PB_N - 1][i4] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][_PB_N - 2][i4])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][_PB_N - 3][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][_PB_N - 2][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][_PB_N - 2][i4])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][_PB_N - 2][i4 - 1]))) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][_PB_N - 2][i4]);
                            }
                          } else {
                            if (i2 >= 32 * ii2 + 2) {
                              for (int i4 = 32 * ii3 + 1; i4 < _PB_N - 1; i4 += 1) {
                                B[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4])) + A[i2 - 1][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 5][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4])) + A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 3][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4])) + A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4 - 1]))) + A[i2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4]);
                              }
                            }
                            if (_PB_N + 2 * i0 + i2 >= 64 * ii0 + 32 * ii2 + 32 * ii3 + 9) {
                              for (int i3 = max(max(64 * ii0 + 32 * ii3 - 2 * i0 + 5, 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 5), -_PB_N + 64 * ii0 + 32 * ii2 + 64 * ii3 - 2 * i0 - i2 + 10); i3 < _PB_N - 1; i3 += 1) {
                                if (64 * ii2 + i3 + 1 >= 32 * ii3 + 2 * i2) {
                                  for (int i4 = max(64 * ii0 + 32 * ii3 - 2 * i0 + 4, 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 5); i4 <= min(64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 6, 64 * ii0 + 32 * ii2 - 64 * ii3 - 2 * i0 - i2 + 3 * i3 + 2); i4 += 1) {
                                    B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                                  }
                                }
                                if (32 * ii3 + 1 >= i3) {
                                  for (int i4 = 64 * ii0 + 64 * ii3 - 2 * i0 - i3 + 5; i4 <= min(64 * ii0 + 32 * ii3 - 2 * i0 + 5, 64 * ii0 - 64 * ii2 + 32 * ii3 - 2 * i0 + 2 * i2 + 1); i4 += 1) {
                                    B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                                  }
                                  if (i2 >= 32 * ii2 + 2 && 32 * ii3 >= i3 + 1) {
                                    B[i2][i3][64 * ii0 + 64 * ii3 - 2 * i0 - i3 + 5] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][64 * ii0 + 64 * ii3 - 2 * i0 - i3 + 5] - (SCALAR_VAL(2.0) * A[i2][i3][64 * ii0 + 64 * ii3 - 2 * i0 - i3 + 5])) + A[i2 - 1][i3][64 * ii0 + 64 * ii3 - 2 * i0 - i3 + 5])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][64 * ii0 + 64 * ii3 - 2 * i0 - i3 + 5] - (SCALAR_VAL(2.0) * A[i2][i3][64 * ii0 + 64 * ii3 - 2 * i0 - i3 + 5])) + A[i2][i3 - 1][64 * ii0 + 64 * ii3 - 2 * i0 - i3 + 5]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][64 * ii0 + 64 * ii3 - 2 * i0 - i3 + 6] - (SCALAR_VAL(2.0) * A[i2][i3][64 * ii0 + 64 * ii3 - 2 * i0 - i3 + 5])) + A[i2][i3][64 * ii0 + 64 * ii3 - 2 * i0 - i3 + 4]))) + A[i2][i3][64 * ii0 + 64 * ii3 - 2 * i0 - i3 + 5]);
                                  }
                                } else if (32 * ii3 + 2 * i2 >= 64 * ii2 + i3 + 2 && 32 * ii3 + 4 >= i3) {
                                  B[i2][i3][64 * ii0 + 32 * ii3 - 2 * i0 + 4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][64 * ii0 + 32 * ii3 - 2 * i0 + 4] - (SCALAR_VAL(2.0) * A[i2][i3][64 * ii0 + 32 * ii3 - 2 * i0 + 4])) + A[i2 - 1][i3][64 * ii0 + 32 * ii3 - 2 * i0 + 4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 4] - (SCALAR_VAL(2.0) * A[i2][i3][64 * ii0 + 32 * ii3 - 2 * i0 + 4])) + A[i2][i3 - 1][64 * ii0 + 32 * ii3 - 2 * i0 + 4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][64 * ii0 + 32 * ii3 - 2 * i0 + 5] - (SCALAR_VAL(2.0) * A[i2][i3][64 * ii0 + 32 * ii3 - 2 * i0 + 4])) + A[i2][i3][64 * ii0 + 32 * ii3 - 2 * i0 + 3]))) + A[i2][i3][64 * ii0 + 32 * ii3 - 2 * i0 + 4]);
                                }
                                for (int i4 = max(max(64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 6, 64 * ii0 - 64 * ii2 + 64 * ii3 - 2 * i0 + 2 * i2 - i3 + 3), 64 * ii0 + 32 * ii2 + 64 * ii3 - 2 * i0 - i2 - i3 + 6); i4 <= 64 * ii0 + 32 * ii2 + 64 * ii3 - 2 * i0 - i2 - i3 + 7; i4 += 1) {
                                  B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                                }
                                for (int i4 = max(max(max(max(64 * ii0 + 32 * ii3 - 2 * i0 + 4, 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 7), 64 * ii0 + 64 * ii3 - 2 * i0 - i3 + 6), 64 * ii0 + 32 * ii2 + 64 * ii3 - 2 * i0 - i2 - i3 + 8), 64 * ii0 + 43 * ii3 - 2 * i0 - (ii3 + i3 + 1) / 3 + 6); i4 < _PB_N - 1; i4 += 1) {
                                  B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                                }
                              }
                            }
                          }
                        } else {
                          B[64 * ii0 + 32 * ii2 - 2 * i0 + 4][_PB_N - 2][_PB_N - 2] = ((((SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][_PB_N - 2][_PB_N - 2] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][_PB_N - 2][_PB_N - 2])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 3][_PB_N - 2][_PB_N - 2])) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][_PB_N - 1][_PB_N - 2] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][_PB_N - 2][_PB_N - 2])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][_PB_N - 3][_PB_N - 2]))) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][_PB_N - 2][_PB_N - 1] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][_PB_N - 2][_PB_N - 2])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][_PB_N - 2][_PB_N - 3]))) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][_PB_N - 2][_PB_N - 2]);
                        }
                      }
                    } else {
                      for (int i2 = max(1, 32 * ii2); i2 <= 32 * ii2 + 31; i2 += 1) {
                        for (int i3 = max(32 * ii3, 32 * ii2 + 32 * ii3 - i2 + 1); i3 < _PB_N - 1; i3 += 1) {
                          for (int i4 = max(max(32 * ii3, 32 * ii2 + 32 * ii3 - i2 + 1), 96 * ii3 - 2 * i3 + 1); i4 < _PB_N - 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                    }
                  }
                  if (32 * ii0 + 16 * ii2 + 16 >= i0) {
                    for (int i2 = max(1, 64 * ii0 + 32 * ii2 - 2 * i0 + 3); i2 <= min(32 * ii2, 64 * ii0 + 32 * ii2 - 2 * i0 + 34); i2 += 1) {
                      for (int i3 = 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 4; i3 < _PB_N - 1; i3 += 1) {
                        for (int i4 = max(64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 4, 64 * ii0 + 32 * ii2 + 64 * ii3 - 2 * i0 - i2 - i3 + 5); i4 < _PB_N - 1; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                    for (int i2 = 32 * ii2 + 1; i2 <= 64 * ii0 + 32 * ii2 - 2 * i0 + 34; i2 += 1) {
                      for (int i3 = 64 * ii0 + 32 * ii3 - 2 * i0 + 3; i3 <= 32 * ii3; i3 += 1) {
                        for (int i4 = 64 * ii0 + 64 * ii3 - 2 * i0 - i3 + 4; i4 < _PB_N - 1; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                      for (int i3 = 32 * ii3 + 1; i3 < _PB_N - 1; i3 += 1) {
                        for (int i4 = 64 * ii0 + 32 * ii3 - 2 * i0 + 3; i4 < _PB_N - 1; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                  }
                }
              }
            }
          } else if (_PB_N >= 32 * ii3 + 35) {
            for (int ii4 = 0; ii4 <= ii2; ii4 += 1) {
              if (ii4 >= 1) {
                if (_PB_N >= 32 * ii4 + 35) {
                  for (int i2 = 32 * ii2 + 1; i2 < _PB_N - 1; i2 += 1) {
                    for (int i3 = 32 * ii3 + 1; i3 <= 32 * ii3 + 32; i3 += 1) {
                      for (int i4 = 32 * ii4 + 1; i4 <= 32 * ii4 + 32; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                  }
                } else {
                  for (int i2 = 32 * ii2 + 1; i2 < _PB_N - 1; i2 += 1) {
                    for (int i3 = 32 * ii3 + 1; i3 <= 32 * ii3 + 32; i3 += 1) {
                      for (int i4 = 32 * ii2 + 1; i4 < _PB_N - 1; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                  }
                }
                for (int i0 = 32 * ii0 + 2; i0 <= min(_PB_TSTEPS, 32 * ii0 + 16); i0 += 1) {
                  if (_PB_N >= 32 * ii4 + 35) {
                    if (i0 >= 32 * ii0 + 3) {
                      for (int i2 = 64 * ii0 + 32 * ii2 - 2 * i0 + 4; i2 < _PB_N - 1; i2 += 1) {
                        if (i2 >= 32 * ii2 + 1) {
                          for (int i3 = max(1, 64 * ii0 + 32 * ii3 - 2 * i0 + 4); i3 <= min(32 * ii3 + 2, -16 * ii2 + 32 * ii3 + (i2 + 1) / 2); i3 += 1) {
                            for (int i4 = max(64 * ii0 + 32 * ii4 - 2 * i0 + 4, 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5); i4 <= min(min(32 * ii4 + 1, 32 * ii3 + 32 * ii4 - i3 + 2), 64 * ii0 - 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 + i2 - i3 + 4); i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            if (64 * ii0 + 32 * ii2 + 32 * ii3 + 36 >= 2 * i0 + i2 + i3) {
                              for (int i4 = 64 * ii0 - 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 + i2 - i3 + 5; i4 <= min(32 * ii4 + 1, 32 * ii3 + 32 * ii4 - i3 + 2); i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                              if (32 * ii2 + 2 * i0 >= 64 * ii0 + i2 + 3 && i3 == 32 * ii3 + 2) {
                                B[i2][32 * ii3 + 2][32 * ii4 + 1] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii3 + 2][32 * ii4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 2][32 * ii4 + 1])) + A[i2 - 1][32 * ii3 + 2][32 * ii4 + 1])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 3][32 * ii4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 2][32 * ii4 + 1])) + A[i2][32 * ii3 + 1][32 * ii4 + 1]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 2][32 * ii4 + 2] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 2][32 * ii4 + 1])) + A[i2][32 * ii3 + 2][32 * ii4]))) + A[i2][32 * ii3 + 2][32 * ii4 + 1]);
                              }
                            }
                            if (32 * ii2 + 2 * i0 >= 64 * ii0 + i2 + 3 && i3 == 32 * ii3 + 2) {
                              for (int i4 = 32 * ii4 + 2; i4 <= 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 36; i4 += 1) {
                                B[i2][32 * ii3 + 2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii3 + 2][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 2][i4])) + A[i2 - 1][32 * ii3 + 2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 3][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 2][i4])) + A[i2][32 * ii3 + 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 2][i4])) + A[i2][32 * ii3 + 2][i4 - 1]))) + A[i2][32 * ii3 + 2][i4]);
                              }
                            }
                            if (2 * i0 + i2 >= 64 * ii0 + 32 * ii2 + 35 && i3 == 32 * ii3 + 2) {
                              for (int i4 = 32 * ii4 + 1; i4 <= 64 * ii0 - 32 * ii2 + 32 * ii4 - 2 * i0 + i2 + 2; i4 += 1) {
                                B[i2][32 * ii3 + 2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii3 + 2][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 2][i4])) + A[i2 - 1][32 * ii3 + 2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 3][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 2][i4])) + A[i2][32 * ii3 + 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 2][i4])) + A[i2][32 * ii3 + 2][i4 - 1]))) + A[i2][32 * ii3 + 2][i4]);
                              }
                              if (64 * ii0 + i2 + 2 >= 32 * ii2 + 2 * i0) {
                                for (int i4 = 64 * ii0 - 32 * ii2 + 32 * ii4 - 2 * i0 + i2 + 3; i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 + 35; i4 += 1) {
                                  B[i2][32 * ii3 + 2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii3 + 2][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 2][i4])) + A[i2 - 1][32 * ii3 + 2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 3][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 2][i4])) + A[i2][32 * ii3 + 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 2][i4])) + A[i2][32 * ii3 + 2][i4 - 1]))) + A[i2][32 * ii3 + 2][i4]);
                                }
                              }
                              for (int i4 = 64 * ii0 - 32 * ii2 + 32 * ii4 - 2 * i0 + i2 + 3; i4 <= 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 36; i4 += 1) {
                                B[i2][32 * ii3 + 2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii3 + 2][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 2][i4])) + A[i2 - 1][32 * ii3 + 2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 3][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 2][i4])) + A[i2][32 * ii3 + 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 2][i4])) + A[i2][32 * ii3 + 2][i4 - 1]))) + A[i2][32 * ii3 + 2][i4]);
                              }
                            } else {
                              if (64 * ii0 + i2 + 2 >= 32 * ii2 + 2 * i0 && 64 * ii0 + 32 * ii2 + 34 >= 2 * i0 + i2 && i3 == 32 * ii3 + 2) {
                                for (int i4 = 32 * ii4 + 1; i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 + 35; i4 += 1) {
                                  B[i2][32 * ii3 + 2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii3 + 2][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 2][i4])) + A[i2 - 1][32 * ii3 + 2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 3][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 2][i4])) + A[i2][32 * ii3 + 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 2][i4])) + A[i2][32 * ii3 + 2][i4 - 1]))) + A[i2][32 * ii3 + 2][i4]);
                                }
                              } else {
                                if (32 * ii3 + 1 >= i3 && 64 * ii0 + 32 * ii3 + i2 + 3 >= 32 * ii2 + 2 * i0 + i3) {
                                  for (int i4 = 32 * ii4 + 2; i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 + 35; i4 += 1) {
                                    B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                                  }
                                  if (32 * ii3 >= i3 + 1) {
                                    for (int i4 = 64 * ii0 + 32 * ii4 - 2 * i0 + 36; i4 <= 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 36; i4 += 1) {
                                      B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                                    }
                                  }
                                } else if (i2 == 32 * ii2 + 1 && 2 * i0 + i3 == 64 * ii0 + 32 * ii3 + 5) {
                                  for (int i4 = 32 * ii4 + 2; i4 <= 32 * ii4 + 31; i4 += 1) {
                                    B[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 5][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][64 * ii0 + 32 * ii3 - 2 * i0 + 5][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 5][i4])) + A[32 * ii2][64 * ii0 + 32 * ii3 - 2 * i0 + 5][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 6][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 5][i4])) + A[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 5][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 5][i4])) + A[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 5][i4 - 1]))) + A[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 5][i4]);
                                  }
                                }
                                if (2 * i0 + i2 + i3 >= 64 * ii0 + 32 * ii2 + 32 * ii3 + 37 && 32 * ii3 + 1 >= i3) {
                                  for (int i4 = 64 * ii0 - 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 + i2 - i3 + 5; i4 <= 32 * ii4 + 1; i4 += 1) {
                                    B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                                  }
                                }
                              }
                              if (64 * ii0 + i2 + 3 >= 32 * ii2 + 2 * i0 && i3 == 32 * ii3) {
                                B[i2][32 * ii3][64 * ii0 + 32 * ii4 - 2 * i0 + 36] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii3][64 * ii0 + 32 * ii4 - 2 * i0 + 36] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][64 * ii0 + 32 * ii4 - 2 * i0 + 36])) + A[i2 - 1][32 * ii3][64 * ii0 + 32 * ii4 - 2 * i0 + 36])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 1][64 * ii0 + 32 * ii4 - 2 * i0 + 36] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][64 * ii0 + 32 * ii4 - 2 * i0 + 36])) + A[i2][32 * ii3 - 1][64 * ii0 + 32 * ii4 - 2 * i0 + 36]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3][64 * ii0 + 32 * ii4 - 2 * i0 + 37] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][64 * ii0 + 32 * ii4 - 2 * i0 + 36])) + A[i2][32 * ii3][64 * ii0 + 32 * ii4 - 2 * i0 + 35]))) + A[i2][32 * ii3][64 * ii0 + 32 * ii4 - 2 * i0 + 36]);
                              }
                            }
                            if (32 * ii2 + 2 * i0 >= 64 * ii0 + i2 + 3 && i3 == 32 * ii3 + 2) {
                              for (int i4 = max(64 * ii0 - 32 * ii2 + 32 * ii4 - 2 * i0 + i2 + 3, 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 37); i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 + 35; i4 += 1) {
                                B[i2][32 * ii3 + 2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii3 + 2][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 2][i4])) + A[i2 - 1][32 * ii3 + 2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 3][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 2][i4])) + A[i2][32 * ii3 + 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 2][i4])) + A[i2][32 * ii3 + 2][i4 - 1]))) + A[i2][32 * ii3 + 2][i4]);
                              }
                            } else if (2 * i0 + i2 + i3 >= 64 * ii0 + 32 * ii2 + 32 * ii3 + 7 && 32 * ii2 + 2 * i0 + i3 >= 64 * ii0 + 32 * ii3 + i2 + 4 && 32 * ii3 + 1 >= i3) {
                              for (int i4 = 32 * ii4 + 2; i4 <= 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 36; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                          }
                          if (32 * ii2 + 2 >= i2) {
                            for (int i4 = 64 * ii0 + 32 * ii4 - 2 * i0 + 4; i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 + 35; i4 += 1) {
                              B[i2][32 * ii3 + 2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii3 + 2][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 2][i4])) + A[i2 - 1][32 * ii3 + 2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 3][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 2][i4])) + A[i2][32 * ii3 + 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 2][i4])) + A[i2][32 * ii3 + 2][i4 - 1]))) + A[i2][32 * ii3 + 2][i4]);
                            }
                          }
                          for (int i3 = 32 * ii3 + 3; i3 <= 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 36; i3 += 1) {
                            B[i2][i3][64 * ii0 + 32 * ii4 - 2 * i0 + 4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][64 * ii0 + 32 * ii4 - 2 * i0 + 4] - (SCALAR_VAL(2.0) * A[i2][i3][64 * ii0 + 32 * ii4 - 2 * i0 + 4])) + A[i2 - 1][i3][64 * ii0 + 32 * ii4 - 2 * i0 + 4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][64 * ii0 + 32 * ii4 - 2 * i0 + 4] - (SCALAR_VAL(2.0) * A[i2][i3][64 * ii0 + 32 * ii4 - 2 * i0 + 4])) + A[i2][i3 - 1][64 * ii0 + 32 * ii4 - 2 * i0 + 4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][64 * ii0 + 32 * ii4 - 2 * i0 + 5] - (SCALAR_VAL(2.0) * A[i2][i3][64 * ii0 + 32 * ii4 - 2 * i0 + 4])) + A[i2][i3][64 * ii0 + 32 * ii4 - 2 * i0 + 3]))) + A[i2][i3][64 * ii0 + 32 * ii4 - 2 * i0 + 4]);
                            for (int i4 = 64 * ii0 + 32 * ii4 - 2 * i0 + 5; i4 <= min(32 * ii4 - 1, 64 * ii0 - 32 * ii2 + 32 * ii4 - 2 * i0 + i2 + 2); i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            for (int i4 = max(64 * ii0 + 32 * ii4 - 2 * i0 + 5, 64 * ii0 - 32 * ii2 + 32 * ii4 - 2 * i0 + i2 + 3); i4 < 32 * ii4; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            for (int i4 = 32 * ii4; i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 + 35; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                          for (int i3 = max(32 * ii3 + 3, 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 37); i3 <= 64 * ii0 + 32 * ii3 - 2 * i0 + 35; i3 += 1) {
                            for (int i4 = 64 * ii0 + 32 * ii4 - 2 * i0 + 4; i4 <= min(32 * ii4, 64 * ii0 - 32 * ii2 + 32 * ii4 - 2 * i0 + i2 + 2); i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            if (64 * ii0 + i2 + 2 >= 32 * ii2 + 2 * i0) {
                              for (int i4 = 32 * ii4 + 1; i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 + 35; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            } else {
                              for (int i4 = 64 * ii0 - 32 * ii2 + 32 * ii4 - 2 * i0 + i2 + 3; i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 + 35; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                          }
                        } else {
                          for (int i3 = max(1, 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 5); i3 <= 32 * ii3 + 1; i3 += 1) {
                            for (int i4 = 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 6; i4 <= 32 * ii4 + 1; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            if (2 * i0 + i2 >= 64 * ii0 + 32 * ii2 + 7 && 2 * i0 + i2 + i3 >= 64 * ii0 + 32 * ii2 + 32 * ii3 + 7) {
                              for (int i4 = 32 * ii4 + 2; i4 <= 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 37; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            } else if (32 * ii3 + 2 * i0 + i2 >= 64 * ii0 + 32 * ii2 + i3 + 5 && 64 * ii0 + 32 * ii2 + 32 * ii3 + 6 >= 2 * i0 + i2 + i3) {
                              for (int i4 = 32 * ii4 + 2; i4 <= 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 37; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            } else {
                              for (int i4 = 32 * ii4 + 2; i4 <= 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 36; i4 += 1) {
                                B[i2][32 * ii3 + 1][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 1][i4])) + A[i2 - 1][32 * ii3 + 1][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 2][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 1][i4])) + A[i2][32 * ii3][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 1][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 1][i4])) + A[i2][32 * ii3 + 1][i4 - 1]))) + A[i2][32 * ii3 + 1][i4]);
                              }
                            }
                          }
                          if (2 * i0 + i2 >= 64 * ii0 + 32 * ii2 + 7) {
                            for (int i4 = 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 5; i4 <= 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 36; i4 += 1) {
                              B[i2][32 * ii3 + 2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii3 + 2][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 2][i4])) + A[i2 - 1][32 * ii3 + 2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 3][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 2][i4])) + A[i2][32 * ii3 + 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 2][i4])) + A[i2][32 * ii3 + 2][i4 - 1]))) + A[i2][32 * ii3 + 2][i4]);
                            }
                          } else {
                            for (int i4 = 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 5; i4 <= 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 36; i4 += 1) {
                              B[i2][32 * ii3 + 2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii3 + 2][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 2][i4])) + A[i2 - 1][32 * ii3 + 2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 3][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 2][i4])) + A[i2][32 * ii3 + 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 2][i4])) + A[i2][32 * ii3 + 2][i4 - 1]))) + A[i2][32 * ii3 + 2][i4]);
                            }
                          }
                          for (int i3 = 32 * ii3 + 3; i3 <= 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 36; i3 += 1) {
                            for (int i4 = 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 5; i4 <= 32 * ii4; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            if (2 * i0 + i2 >= 64 * ii0 + 32 * ii2 + 6) {
                              for (int i4 = 32 * ii4 + 1; i4 <= 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 36; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            } else {
                              for (int i4 = 32 * ii4 + 1; i4 <= 128 * ii0 + 64 * ii2 + 32 * ii4 - 4 * i0 - 2 * i2 + 40; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                              if (2 * i0 + i2 == 64 * ii0 + 32 * ii2 + 5) {
                                B[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3][32 * ii4 + 31] = ((((SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3][32 * ii4 + 31] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3][32 * ii4 + 31])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][i3][32 * ii4 + 31])) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3 + 1][32 * ii4 + 31] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3][32 * ii4 + 31])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3 - 1][32 * ii4 + 31]))) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3][32 * ii4 + 32] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3][32 * ii4 + 31])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3][32 * ii4 + 30]))) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3][32 * ii4 + 31]);
                              }
                            }
                          }
                        }
                      }
                    } else {
                      for (int i2 = 32 * ii2; i2 <= min(_PB_N - 2, 32 * ii2 + 2); i2 += 1) {
                        for (int i3 = max(max(1, 32 * ii3), 32 * ii2 + 32 * ii3 - i2 + 1); i3 <= 32 * ii2 + 32 * ii3 - i2 + 32; i3 += 1) {
                          for (int i4 = max(max(32 * ii4, 32 * ii2 + 32 * ii4 - i2 + 1), 64 * ii3 + 32 * ii4 - 2 * i3 + 1); i4 <= 32 * ii2 + 32 * ii4 - i2 + 32; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (i3 == 32 * ii3) {
                            for (int i4 = 32 * ii2 + 32 * ii4 - i2 + 33; i4 <= 32 * ii4 + 32; i4 += 1) {
                              B[i2][32 * ii3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii3][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][i4])) + A[i2 - 1][32 * ii3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][i4])) + A[i2][32 * ii3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][i4])) + A[i2][32 * ii3][i4 - 1]))) + A[i2][32 * ii3][i4]);
                            }
                          } else if (i2 == 32 * ii2 + 2) {
                            B[32 * ii2 + 2][i3][32 * ii4 + 31] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 3][i3][32 * ii4 + 31] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][i3][32 * ii4 + 31])) + A[32 * ii2 + 1][i3][32 * ii4 + 31])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3 + 1][32 * ii4 + 31] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][i3][32 * ii4 + 31])) + A[32 * ii2 + 2][i3 - 1][32 * ii4 + 31]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3][32 * ii4 + 32] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][i3][32 * ii4 + 31])) + A[32 * ii2 + 2][i3][32 * ii4 + 30]))) + A[32 * ii2 + 2][i3][32 * ii4 + 31]);
                          }
                        }
                        if (i2 == 32 * ii2 + 2) {
                          for (int i4 = 32 * ii4; i4 <= 32 * ii4 + 31; i4 += 1) {
                            B[32 * ii2 + 2][32 * ii3 + 31][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 3][32 * ii3 + 31][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][32 * ii3 + 31][i4])) + A[32 * ii2 + 1][32 * ii3 + 31][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][32 * ii3 + 32][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][32 * ii3 + 31][i4])) + A[32 * ii2 + 2][32 * ii3 + 30][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][32 * ii3 + 31][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][32 * ii3 + 31][i4])) + A[32 * ii2 + 2][32 * ii3 + 31][i4 - 1]))) + A[32 * ii2 + 2][32 * ii3 + 31][i4]);
                          }
                        }
                      }
                      for (int i2 = 32 * ii2 + 3; i2 <= min(_PB_N - 2, 32 * ii2 + 31); i2 += 1) {
                        for (int i3 = max(1, 32 * ii3); i3 <= 32 * ii3 + 31; i3 += 1) {
                          for (int i4 = max(32 * ii4, 64 * ii3 + 32 * ii4 - 2 * i3 + 1); i4 <= 32 * ii4 + 31; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (i3 == 32 * ii3) {
                            B[i2][32 * ii3][32 * ii4 + 32] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii3][32 * ii4 + 32] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][32 * ii4 + 32])) + A[i2 - 1][32 * ii3][32 * ii4 + 32])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 1][32 * ii4 + 32] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][32 * ii4 + 32])) + A[i2][32 * ii3 - 1][32 * ii4 + 32]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3][32 * ii4 + 33] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][32 * ii4 + 32])) + A[i2][32 * ii3][32 * ii4 + 31]))) + A[i2][32 * ii3][32 * ii4 + 32]);
                          }
                        }
                      }
                      if (32 * ii2 + 34 == _PB_N) {
                        for (int i3 = max(32 * ii3, 31 * ii3 + 1); i3 <= 32 * ii3 + 31; i3 += 1) {
                          for (int i4 = max(32 * ii4, 32 * ii3 + 32 * ii4 - i3 + 1); i4 <= 32 * ii4 + 31; i4 += 1) {
                            B[_PB_N - 2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[_PB_N - 1][i3][i4] - (SCALAR_VAL(2.0) * A[_PB_N - 2][i3][i4])) + A[_PB_N - 3][i3][i4])) + (SCALAR_VAL(0.125) * ((A[_PB_N - 2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[_PB_N - 2][i3][i4])) + A[_PB_N - 2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[_PB_N - 2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[_PB_N - 2][i3][i4])) + A[_PB_N - 2][i3][i4 - 1]))) + A[_PB_N - 2][i3][i4]);
                          }
                          if (i3 == 32 * ii3) {
                            B[_PB_N - 2][32 * ii3][32 * ii4 + 32] = ((((SCALAR_VAL(0.125) * ((A[_PB_N - 1][32 * ii3][32 * ii4 + 32] - (SCALAR_VAL(2.0) * A[_PB_N - 2][32 * ii3][32 * ii4 + 32])) + A[_PB_N - 3][32 * ii3][32 * ii4 + 32])) + (SCALAR_VAL(0.125) * ((A[_PB_N - 2][32 * ii3 + 1][32 * ii4 + 32] - (SCALAR_VAL(2.0) * A[_PB_N - 2][32 * ii3][32 * ii4 + 32])) + A[_PB_N - 2][32 * ii3 - 1][32 * ii4 + 32]))) + (SCALAR_VAL(0.125) * ((A[_PB_N - 2][32 * ii3][32 * ii4 + 33] - (SCALAR_VAL(2.0) * A[_PB_N - 2][32 * ii3][32 * ii4 + 32])) + A[_PB_N - 2][32 * ii3][32 * ii4 + 31]))) + A[_PB_N - 2][32 * ii3][32 * ii4 + 32]);
                          }
                        }
                      }
                    }
                  } else {
                    for (int i2 = 64 * ii0 + 32 * ii2 - 2 * i0 + 4; i2 < _PB_N - 1; i2 += 1) {
                      if (32 * ii2 >= i2 + 1) {
                        for (int i3 = max(1, 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 5); i3 <= 32 * ii3 + 2; i3 += 1) {
                          for (int i4 = max(64 * ii0 + 64 * ii2 - 2 * i0 - i2 + 5, 64 * ii0 + 64 * ii2 + 32 * ii3 - 2 * i0 - i2 - i3 + 6); i4 < _PB_N - 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                        for (int i3 = 32 * ii3 + 3; i3 <= 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 36; i3 += 1) {
                          if (32 * ii3 + 31 >= i3) {
                            if (2 * i0 + i2 >= 64 * ii0 + 32 * ii2 + 6) {
                              for (int i4 = 64 * ii0 + 64 * ii2 - 2 * i0 - i2 + 5; i4 <= 32 * ii2 + 1; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                            for (int i4 = max(-128 * ii0 - 32 * ii2 + 4 * i0 + 2 * i2 - 10, 128 * ii0 + 96 * ii2 - 4 * i0 - 2 * i2 + 9); i4 <= 32 * ii2 + 1; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            for (int i4 = 32 * ii2 + 2; i4 < _PB_N - 1; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          } else {
                            for (int i4 = 32 * ii2 + 1; i4 < _PB_N - 1; i4 += 1) {
                              B[64 * ii0 + 32 * ii2 - 2 * i0 + 4][32 * ii3 + 32][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][32 * ii3 + 32][i4] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][32 * ii3 + 32][i4])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 3][32 * ii3 + 32][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][32 * ii3 + 33][i4] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][32 * ii3 + 32][i4])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][32 * ii3 + 31][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][32 * ii3 + 32][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][32 * ii3 + 32][i4])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][32 * ii3 + 32][i4 - 1]))) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][32 * ii3 + 32][i4]);
                            }
                          }
                        }
                      } else {
                        if (i2 == 32 * ii2) {
                          for (int i3 = max(1, 64 * ii0 + 32 * ii3 - 2 * i0 + 5); i3 <= 32 * ii3; i3 += 1) {
                            for (int i4 = 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i3 + 6; i4 < _PB_N - 1; i4 += 1) {
                              B[32 * ii2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3][i4 - 1]))) + A[32 * ii2][i3][i4]);
                            }
                          }
                        } else {
                          for (int i3 = max(1, 64 * ii0 + 32 * ii3 - 2 * i0 + 4); i3 <= min(32 * ii3 - 1, 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 36); i3 += 1) {
                            for (int i4 = 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i3 + 5; i4 <= min(32 * ii2 + 1, 64 * ii0 + 32 * ii3 - 2 * i0 + i2 - i3 + 4); i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            if (32 * ii2 + 2 * i0 + i3 >= 64 * ii0 + 32 * ii3 + i2 + 4) {
                              for (int i4 = 64 * ii0 + 32 * ii3 - 2 * i0 + i2 - i3 + 5; i4 < _PB_N - 1; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            } else {
                              for (int i4 = 32 * ii2 + 2; i4 < _PB_N - 1; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                          }
                        }
                        for (int i3 = max(max(1, 32 * ii3), 32 * ii2 + 32 * ii3 - i2 + 1); i3 <= min(32 * ii3 + 2, 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 36); i3 += 1) {
                          for (int i4 = max(64 * ii0 + 32 * ii2 - 2 * i0 + 4, 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i3 + 5); i4 <= min(64 * ii0 + 32 * ii3 - 2 * i0 + i2 - i3 + 4, 32 * ii2 + 16 * ii3 - i3 + i3 / 2); i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (i0 >= 32 * ii0 + 3 && 64 * ii0 + i2 + 3 == 32 * ii2 + 2 * i0 && 32 * ii3 + 1 >= i3) {
                            B[-64 * ii0 + 32 * ii2 + 2 * i0 - 3][i3][32 * ii2 + 32 * ii3 - i3 + 1] = ((((SCALAR_VAL(0.125) * ((A[-64 * ii0 + 32 * ii2 + 2 * i0 - 2][i3][32 * ii2 + 32 * ii3 - i3 + 1] - (SCALAR_VAL(2.0) * A[-64 * ii0 + 32 * ii2 + 2 * i0 - 3][i3][32 * ii2 + 32 * ii3 - i3 + 1])) + A[-64 * ii0 + 32 * ii2 + 2 * i0 - 4][i3][32 * ii2 + 32 * ii3 - i3 + 1])) + (SCALAR_VAL(0.125) * ((A[-64 * ii0 + 32 * ii2 + 2 * i0 - 3][i3 + 1][32 * ii2 + 32 * ii3 - i3 + 1] - (SCALAR_VAL(2.0) * A[-64 * ii0 + 32 * ii2 + 2 * i0 - 3][i3][32 * ii2 + 32 * ii3 - i3 + 1])) + A[-64 * ii0 + 32 * ii2 + 2 * i0 - 3][i3 - 1][32 * ii2 + 32 * ii3 - i3 + 1]))) + (SCALAR_VAL(0.125) * ((A[-64 * ii0 + 32 * ii2 + 2 * i0 - 3][i3][32 * ii2 + 32 * ii3 - i3 + 2] - (SCALAR_VAL(2.0) * A[-64 * ii0 + 32 * ii2 + 2 * i0 - 3][i3][32 * ii2 + 32 * ii3 - i3 + 1])) + A[-64 * ii0 + 32 * ii2 + 2 * i0 - 3][i3][32 * ii2 + 32 * ii3 - i3]))) + A[-64 * ii0 + 32 * ii2 + 2 * i0 - 3][i3][32 * ii2 + 32 * ii3 - i3 + 1]);
                            if (i3 == 32 * ii3) {
                              for (int i4 = 32 * ii2 + 2; i4 < _PB_N - 1; i4 += 1) {
                                B[-64 * ii0 + 32 * ii2 + 2 * i0 - 3][32 * ii3][i4] = ((((SCALAR_VAL(0.125) * ((A[-64 * ii0 + 32 * ii2 + 2 * i0 - 2][32 * ii3][i4] - (SCALAR_VAL(2.0) * A[-64 * ii0 + 32 * ii2 + 2 * i0 - 3][32 * ii3][i4])) + A[-64 * ii0 + 32 * ii2 + 2 * i0 - 4][32 * ii3][i4])) + (SCALAR_VAL(0.125) * ((A[-64 * ii0 + 32 * ii2 + 2 * i0 - 3][32 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[-64 * ii0 + 32 * ii2 + 2 * i0 - 3][32 * ii3][i4])) + A[-64 * ii0 + 32 * ii2 + 2 * i0 - 3][32 * ii3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[-64 * ii0 + 32 * ii2 + 2 * i0 - 3][32 * ii3][i4 + 1] - (SCALAR_VAL(2.0) * A[-64 * ii0 + 32 * ii2 + 2 * i0 - 3][32 * ii3][i4])) + A[-64 * ii0 + 32 * ii2 + 2 * i0 - 3][32 * ii3][i4 - 1]))) + A[-64 * ii0 + 32 * ii2 + 2 * i0 - 3][32 * ii3][i4]);
                              }
                            }
                          }
                          if (i0 >= 32 * ii0 + 3 && 32 * ii2 + 2 * i0 >= 64 * ii0 + i2 + 3 && 32 * ii2 + 2 * i0 + i3 >= 64 * ii0 + 32 * ii3 + i2 + 4) {
                            for (int i4 = max(64 * ii0 + 64 * ii2 - 2 * i0 - i2 + 5, 64 * ii0 + 32 * ii3 - 2 * i0 + i2 - i3 + 5); i4 < _PB_N - 1; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          } else if (64 * ii0 + i2 + 2 >= 32 * ii2 + 2 * i0) {
                            for (int i4 = 32 * ii2 + 16 * ii3 - i3 + i3 / 2 + 1; i4 < _PB_N - 1; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                          if (i0 == 32 * ii0 + 2 && 32 * ii2 + 1 >= i2) {
                            for (int i4 = max(64 * ii2 - i2 + 1, 32 * ii2 + 16 * ii3 - i3 + i3 / 2 + 1); i4 < _PB_N - 1; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                        }
                        if (32 * ii2 + 1 >= i2) {
                          for (int i3 = 32 * ii3 + 3; i3 <= 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 36; i3 += 1) {
                            if (2 * i0 + i2 >= 64 * ii0 + 32 * ii2 + 6) {
                              for (int i4 = 64 * ii0 + 64 * ii2 - 2 * i0 - i2 + 5; i4 <= 32 * ii2 + 1; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                              for (int i4 = 32 * ii2 + 2; i4 < _PB_N - 1; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            } else {
                              for (int i4 = 64 * ii2 - i2 + 1; i4 < _PB_N - 1; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                          }
                        }
                        for (int i3 = max(1, 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 37); i3 <= min(32 * ii3 - 1, 64 * ii0 - 32 * ii2 + 32 * ii3 - 2 * i0 + i2 + 3); i3 += 1) {
                          for (int i4 = 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i3 + 5; i4 < _PB_N - 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                        for (int i3 = max(max(1, 64 * ii0 - 32 * ii2 + 32 * ii3 - 2 * i0 + i2 + 4), 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 37); i3 < 32 * ii3; i3 += 1) {
                          for (int i4 = 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i3 + 5; i4 < _PB_N - 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                        if (i0 == 32 * ii0 + 2 && i2 >= 32 * ii2 + 2) {
                          for (int i3 = 32 * ii3 + 3; i3 <= 32 * ii2 + 32 * ii3 - i2 + 32; i3 += 1) {
                            for (int i4 = 32 * ii2; i4 < _PB_N - 1; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                        } else if (i0 >= 32 * ii0 + 3 && i2 >= 32 * ii2 + 2) {
                          for (int i3 = 32 * ii3 + 3; i3 <= 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 36; i3 += 1) {
                            for (int i4 = 64 * ii0 + 32 * ii2 - 2 * i0 + 4; i4 <= min(32 * ii2 - 1, 64 * ii0 - 2 * i0 + i2 + 2); i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            for (int i4 = 64 * ii0 - 2 * i0 + i2 + 3; i4 < 32 * ii2; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            for (int i4 = 32 * ii2; i4 < _PB_N - 1; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                        }
                        for (int i3 = max(max(1, 32 * ii3), 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 37); i3 <= 64 * ii0 + 32 * ii3 - 2 * i0 + 35; i3 += 1) {
                          if (32 * ii3 + 2 >= i3) {
                            for (int i4 = max(64 * ii0 + 32 * ii2 - 2 * i0 + 4, 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i3 + 5); i4 <= min(64 * ii0 + 32 * ii3 - 2 * i0 + i2 - i3 + 4, 32 * ii2 + 16 * ii3 - i3 + i3 / 2); i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          } else {
                            for (int i4 = 64 * ii0 + 32 * ii2 - 2 * i0 + 4; i4 <= min(32 * ii2 - 1, 64 * ii0 - 2 * i0 + i2 + 2); i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                          for (int i4 = max(64 * ii0 - 2 * i0 + i2 + 3, 64 * ii0 + 32 * ii3 - 2 * i0 + i2 - i3 + 5); i4 < 32 * ii2; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (32 * ii2 + 2 * i0 >= 64 * ii0 + i2 + 5 && i3 == 32 * ii3) {
                            B[i2][32 * ii3][32 * ii2] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii3][32 * ii2] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][32 * ii2])) + A[i2 - 1][32 * ii3][32 * ii2])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 1][32 * ii2] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][32 * ii2])) + A[i2][32 * ii3 - 1][32 * ii2]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3][32 * ii2 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][32 * ii2])) + A[i2][32 * ii3][32 * ii2 - 1]))) + A[i2][32 * ii3][32 * ii2]);
                          }
                          for (int i4 = max(32 * ii2, 32 * ii2 + 16 * ii3 - i3 + i3 / 2 + 1); i4 < _PB_N - 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                    }
                  }
                  if (ii4 == ii2) {
                    for (int i2 = 64 * ii0 + 32 * ii2 - 2 * i0 + 3; i2 <= 32 * ii2; i2 += 1) {
                      for (int i3 = max(1, 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 4); i3 <= 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 35; i3 += 1) {
                        for (int i4 = max(64 * ii0 + 64 * ii2 - 2 * i0 - i2 + 4, 64 * ii0 + 64 * ii2 + 32 * ii3 - 2 * i0 - i2 - i3 + 5); i4 < _PB_N - 1; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                  } else {
                    for (int i2 = 64 * ii0 + 32 * ii2 - 2 * i0 + 3; i2 <= 32 * ii2; i2 += 1) {
                      for (int i3 = max(1, 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 4); i3 <= 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 35; i3 += 1) {
                        for (int i4 = max(64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 4, 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 5); i4 <= 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 36; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                        for (int i4 = 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 37; i4 <= 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 35; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                  }
                  for (int i2 = 32 * ii2 + 1; i2 < _PB_N - 1; i2 += 1) {
                    if (ii4 == ii2) {
                      for (int i3 = max(1, 64 * ii0 + 32 * ii3 - 2 * i0 + 3); i3 <= 32 * ii3; i3 += 1) {
                        for (int i4 = 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i3 + 4; i4 < _PB_N - 1; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    } else {
                      for (int i3 = max(1, 64 * ii0 + 32 * ii3 - 2 * i0 + 3); i3 <= 32 * ii3; i3 += 1) {
                        for (int i4 = 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 4; i4 <= 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 35; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                    for (int i3 = 32 * ii3 + 1; i3 <= 64 * ii0 + 32 * ii3 - 2 * i0 + 34; i3 += 1) {
                      for (int i4 = 64 * ii0 + 32 * ii4 - 2 * i0 + 3; i4 <= 32 * ii4; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                      if (_PB_N >= 32 * ii4 + 35) {
                        for (int i4 = 32 * ii4 + 1; i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 + 34; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      } else {
                        for (int i4 = 32 * ii2 + 1; i4 < _PB_N - 1; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                  }
                }
                if (_PB_N >= 32 * ii4 + 35) {
                  for (int i0 = 32 * ii0 + 17; i0 <= min(_PB_TSTEPS, 32 * ii0 + 32); i0 += 1) {
                    for (int i2 = 64 * ii0 + 32 * ii2 - 2 * i0 + 4; i2 <= min(32 * ii2, 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 + 35); i2 += 1) {
                      for (int i3 = max(1, 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 5); i3 <= min(32 * ii3 + 1, 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 36); i3 += 1) {
                        for (int i4 = 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 6; i4 <= 32 * ii4 + 1; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        if (2 * i0 + i2 >= 64 * ii0 + 32 * ii2 + 7 && 2 * i0 + i2 + i3 >= 64 * ii0 + 32 * ii2 + 32 * ii3 + 7) {
                          for (int i4 = 32 * ii4 + 2; i4 <= 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 37; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        } else if (32 * ii3 + 2 * i0 + i2 >= 64 * ii0 + 32 * ii2 + i3 + 5 && 64 * ii0 + 32 * ii2 + 32 * ii3 + 6 >= 2 * i0 + i2 + i3) {
                          for (int i4 = 32 * ii4 + 2; i4 <= 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 37; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        } else {
                          for (int i4 = 32 * ii4 + 2; i4 <= 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 36; i4 += 1) {
                            B[i2][32 * ii3 + 1][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 1][i4])) + A[i2 - 1][32 * ii3 + 1][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 2][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 1][i4])) + A[i2][32 * ii3][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 1][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 1][i4])) + A[i2][32 * ii3 + 1][i4 - 1]))) + A[i2][32 * ii3 + 1][i4]);
                          }
                        }
                      }
                      if (2 * i0 + i2 >= 64 * ii0 + 32 * ii2 + 7 && 32 * ii2 >= i2 + 1 && 64 * ii0 + 32 * ii2 + 34 >= 2 * i0 + i2) {
                        for (int i4 = 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 5; i4 <= 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 36; i4 += 1) {
                          B[i2][32 * ii3 + 2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii3 + 2][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 2][i4])) + A[i2 - 1][32 * ii3 + 2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 3][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 2][i4])) + A[i2][32 * ii3 + 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 2][i4])) + A[i2][32 * ii3 + 2][i4 - 1]))) + A[i2][32 * ii3 + 2][i4]);
                        }
                      } else if (64 * ii0 + 32 * ii2 + 6 >= 2 * i0 + i2) {
                        for (int i4 = 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 5; i4 <= 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 36; i4 += 1) {
                          B[i2][32 * ii3 + 2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii3 + 2][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 2][i4])) + A[i2 - 1][32 * ii3 + 2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 3][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 2][i4])) + A[i2][32 * ii3 + 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 2][i4])) + A[i2][32 * ii3 + 2][i4 - 1]))) + A[i2][32 * ii3 + 2][i4]);
                        }
                      }
                      for (int i3 = 32 * ii3 + 3; i3 <= 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 36; i3 += 1) {
                        for (int i4 = 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 5; i4 <= 32 * ii4; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        if (64 * ii0 + 32 * ii2 + 5 >= 2 * i0 + i2) {
                          for (int i4 = 32 * ii4 + 1; i4 <= 128 * ii0 + 64 * ii2 + 32 * ii4 - 4 * i0 - 2 * i2 + 40; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (2 * i0 + i2 == 64 * ii0 + 32 * ii2 + 5) {
                            B[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3][32 * ii4 + 31] = ((((SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3][32 * ii4 + 31] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3][32 * ii4 + 31])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][i3][32 * ii4 + 31])) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3 + 1][32 * ii4 + 31] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3][32 * ii4 + 31])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3 - 1][32 * ii4 + 31]))) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3][32 * ii4 + 32] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3][32 * ii4 + 31])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3][32 * ii4 + 30]))) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3][32 * ii4 + 31]);
                          }
                        } else {
                          for (int i4 = 32 * ii4 + 1; i4 <= 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 36; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                      if (i0 == 32 * ii0 + 17 && i2 == 32 * ii2) {
                        for (int i4 = 32 * ii4 - 29; i4 <= 32 * ii4 + 2; i4 += 1) {
                          B[32 * ii2][32 * ii3 + 2][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][32 * ii3 + 2][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][32 * ii3 + 2][i4])) + A[32 * ii2 - 1][32 * ii3 + 2][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2][32 * ii3 + 3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][32 * ii3 + 2][i4])) + A[32 * ii2][32 * ii3 + 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2][32 * ii3 + 2][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2][32 * ii3 + 2][i4])) + A[32 * ii2][32 * ii3 + 2][i4 - 1]))) + A[32 * ii2][32 * ii3 + 2][i4]);
                        }
                      }
                    }
                    for (int i2 = 32 * ii2 + 1; i2 < _PB_N - 1; i2 += 1) {
                      for (int i3 = max(1, 64 * ii0 + 32 * ii3 - 2 * i0 + 4); i3 <= 64 * ii0 + 32 * ii3 - 2 * i0 + 35; i3 += 1) {
                        for (int i4 = 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 5; i4 <= 32 * ii4 + 1; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        if (i2 == 32 * ii2 + 1 && 2 * i0 + i3 == 64 * ii0 + 32 * ii3 + 5) {
                          for (int i4 = 32 * ii4 + 2; i4 <= 32 * ii4 + 31; i4 += 1) {
                            B[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 5][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][64 * ii0 + 32 * ii3 - 2 * i0 + 5][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 5][i4])) + A[32 * ii2][64 * ii0 + 32 * ii3 - 2 * i0 + 5][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 6][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 5][i4])) + A[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 5][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 5][i4])) + A[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 5][i4 - 1]))) + A[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 5][i4]);
                          }
                        } else if (32 * ii3 >= i3 + 1 && 64 * ii0 + 32 * ii3 + i2 + 3 >= 32 * ii2 + 2 * i0 + i3) {
                          for (int i4 = 32 * ii4 + 2; i4 <= 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 36; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        } else if (i0 == 32 * ii0 + 17 && i2 >= 32 * ii2 + 31 && i3 == 32 * ii3) {
                          B[i2][32 * ii3][32 * ii4 + 2] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii3][32 * ii4 + 2] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][32 * ii4 + 2])) + A[i2 - 1][32 * ii3][32 * ii4 + 2])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 1][32 * ii4 + 2] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][32 * ii4 + 2])) + A[i2][32 * ii3 - 1][32 * ii4 + 2]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3][32 * ii4 + 3] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][32 * ii4 + 2])) + A[i2][32 * ii3][32 * ii4 + 1]))) + A[i2][32 * ii3][32 * ii4 + 2]);
                        } else if (32 * ii2 + 2 * i0 + i3 >= 64 * ii0 + 32 * ii3 + i2 + 4) {
                          for (int i4 = 32 * ii4 + 2; i4 <= 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 36; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                    }
                    for (int i2 = 64 * ii0 + 32 * ii2 - 2 * i0 + 3; i2 <= min(32 * ii2, 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 + 34); i2 += 1) {
                      for (int i3 = max(1, 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 4); i3 <= 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 35; i3 += 1) {
                        for (int i4 = max(64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 4, 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 5); i4 <= 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 36; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                        for (int i4 = 64 * ii0 + 32 * ii2 + 32 * ii3 + 32 * ii4 - 2 * i0 - i2 - i3 + 37; i4 <= 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 35; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                    for (int i2 = 32 * ii2 + 1; i2 < _PB_N - 1; i2 += 1) {
                      for (int i3 = max(1, 64 * ii0 + 32 * ii3 - 2 * i0 + 3); i3 <= 64 * ii0 + 32 * ii3 - 2 * i0 + 34; i3 += 1) {
                        for (int i4 = 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 4; i4 <= 64 * ii0 + 32 * ii3 + 32 * ii4 - 2 * i0 - i3 + 35; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                  }
                } else {
                  for (int i0 = 32 * ii0 + 17; i0 <= min(_PB_TSTEPS, 32 * ii0 + 32); i0 += 1) {
                    if (i0 >= 32 * ii0 + 18) {
                      for (int i2 = max(64 * ii0 + 32 * ii2 - 2 * i0 + 4, ii2 - (29 * ii2 + 61) / 61 + 1); i2 <= min(32 * ii2, 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 + 35); i2 += 1) {
                        for (int i3 = max(1, 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 5); i3 <= 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 36; i3 += 1) {
                          if (32 * ii3 + 2 >= i3) {
                            for (int i4 = max(64 * ii0 + 64 * ii2 - 2 * i0 - i2 + 5, 64 * ii0 + 64 * ii2 + 32 * ii3 - 2 * i0 - i2 - i3 + 6); i4 < _PB_N - 1; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          } else if (32 * ii3 + 31 >= i3) {
                            if (2 * i0 + i2 >= 64 * ii0 + 32 * ii2 + 6) {
                              for (int i4 = 64 * ii0 + 64 * ii2 - 2 * i0 - i2 + 5; i4 <= 32 * ii2 + 1; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                            for (int i4 = max(-128 * ii0 - 32 * ii2 + 4 * i0 + 2 * i2 - 10, 128 * ii0 + 96 * ii2 - 4 * i0 - 2 * i2 + 9); i4 <= 32 * ii2 + 1; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            for (int i4 = 32 * ii2 + 2; i4 < _PB_N - 1; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          } else {
                            for (int i4 = 32 * ii2 + 1; i4 < _PB_N - 1; i4 += 1) {
                              B[64 * ii0 + 32 * ii2 - 2 * i0 + 4][32 * ii3 + 32][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][32 * ii3 + 32][i4] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][32 * ii3 + 32][i4])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 3][32 * ii3 + 32][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][32 * ii3 + 33][i4] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][32 * ii3 + 32][i4])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][32 * ii3 + 31][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][32 * ii3 + 32][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][32 * ii3 + 32][i4])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][32 * ii3 + 32][i4 - 1]))) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][32 * ii3 + 32][i4]);
                            }
                          }
                        }
                      }
                      for (int i3 = max(1, 64 * ii0 + 32 * ii3 - 2 * i0 + 4); i3 <= 64 * ii0 + 32 * ii3 - 2 * i0 + 35; i3 += 1) {
                        for (int i4 = 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i3 + 5; i4 < _PB_N - 1; i4 += 1) {
                          B[32 * ii2 + 1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3][i4 - 1]))) + A[32 * ii2 + 1][i3][i4]);
                        }
                      }
                      for (int i2 = 32 * ii2 + 2; i2 <= min(_PB_N - 2, 32 * ii2 + 31); i2 += 1) {
                        for (int i3 = max(1, 64 * ii0 + 32 * ii3 - 2 * i0 + 4); i3 <= 64 * ii0 + 32 * ii3 - 2 * i0 + 35; i3 += 1) {
                          for (int i4 = 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i3 + 5; i4 <= min(32 * ii2 + 1, 64 * ii0 + 32 * ii3 - 2 * i0 + i2 - i3 + 4); i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (32 * ii2 + 2 * i0 + i3 >= 64 * ii0 + 32 * ii3 + i2 + 4) {
                            for (int i4 = 64 * ii0 + 32 * ii3 - 2 * i0 + i2 - i3 + 5; i4 < _PB_N - 1; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          } else {
                            for (int i4 = 32 * ii2 + 2; i4 < _PB_N - 1; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                        }
                      }
                      if (32 * ii2 + 34 == _PB_N) {
                        for (int i3 = max(1, 64 * ii0 + 32 * ii3 - 2 * i0 + 4); i3 <= 64 * ii0 + 32 * ii3 - 2 * i0 + 35; i3 += 1) {
                          for (int i4 = _PB_N + 64 * ii0 + 32 * ii3 - 2 * i0 - i3 - 29; i4 < _PB_N - 1; i4 += 1) {
                            B[_PB_N - 2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[_PB_N - 1][i3][i4] - (SCALAR_VAL(2.0) * A[_PB_N - 2][i3][i4])) + A[_PB_N - 3][i3][i4])) + (SCALAR_VAL(0.125) * ((A[_PB_N - 2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[_PB_N - 2][i3][i4])) + A[_PB_N - 2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[_PB_N - 2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[_PB_N - 2][i3][i4])) + A[_PB_N - 2][i3][i4 - 1]))) + A[_PB_N - 2][i3][i4]);
                          }
                        }
                      }
                    } else {
                      for (int i2 = 32 * ii2 - 30; i2 <= 32 * ii2; i2 += 1) {
                        for (int i3 = max(1, 32 * ii2 + 32 * ii3 - i2 - 29); i3 <= 32 * ii3 + 2; i3 += 1) {
                          for (int i4 = max(64 * ii2 - i2 - 29, 64 * ii2 + 32 * ii3 - i2 - i3 - 28); i4 < _PB_N - 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                        for (int i3 = 32 * ii3 + 3; i3 <= min(32 * ii3 + 31, 32 * ii2 + 32 * ii3 - i2 + 2); i3 += 1) {
                          if (i2 + 28 >= 32 * ii2) {
                            for (int i4 = 64 * ii2 - i2 - 29; i4 <= 32 * ii2 + 1; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                          for (int i4 = max(96 * ii2 - 2 * i2 - 59, -32 * ii2 + 2 * i2 + 58); i4 <= 32 * ii2 + 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          for (int i4 = 32 * ii2 + 2; i4 < _PB_N - 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                        if (i2 + 30 == 32 * ii2) {
                          for (int i4 = 32 * ii2 + 1; i4 < _PB_N - 1; i4 += 1) {
                            B[32 * ii2 - 30][32 * ii3 + 32][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 - 29][32 * ii3 + 32][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 - 30][32 * ii3 + 32][i4])) + A[32 * ii2 - 31][32 * ii3 + 32][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 - 30][32 * ii3 + 33][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 - 30][32 * ii3 + 32][i4])) + A[32 * ii2 - 30][32 * ii3 + 31][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 - 30][32 * ii3 + 32][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 - 30][32 * ii3 + 32][i4])) + A[32 * ii2 - 30][32 * ii3 + 32][i4 - 1]))) + A[32 * ii2 - 30][32 * ii3 + 32][i4]);
                          }
                        }
                      }
                      for (int i3 = max(1, 32 * ii3 - 30); i3 <= 32 * ii3 + 1; i3 += 1) {
                        for (int i4 = 32 * ii2 + 32 * ii3 - i3 - 29; i4 < _PB_N - 1; i4 += 1) {
                          B[32 * ii2 + 1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3][i4 - 1]))) + A[32 * ii2 + 1][i3][i4]);
                        }
                      }
                      for (int i2 = 32 * ii2 + 2; i2 < _PB_N - 1; i2 += 1) {
                        for (int i3 = max(1, 32 * ii3 - 30); i3 <= 32 * ii3 + 1; i3 += 1) {
                          if (i3 >= 32 * ii3) {
                            for (int i4 = 32 * ii2 + 32 * ii3 - i3 - 29; i4 <= min(32 * ii3 + i2 - i3 - 30, 32 * ii2 + 32 * ii3 - i3); i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            if (i2 + i3 >= 32 * ii2 + 32 * ii3 + 3) {
                              for (int i4 = 32 * ii3 + i2 - i3 - 29; i4 <= 32 * ii2 + 32 * ii3 - i3; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                          } else {
                            for (int i4 = 32 * ii2 + 32 * ii3 - i3 - 29; i4 <= min(32 * ii2 + 1, 32 * ii3 + i2 - i3 - 30); i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            if (i2 + i3 >= 32 * ii2 + 32 * ii3 + 3 && 32 * ii2 + i3 + 30 >= 32 * ii3 + i2) {
                              for (int i4 = 32 * ii3 + i2 - i3 - 29; i4 < _PB_N - 1; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                          }
                          if (32 * ii2 + 32 * ii3 + 2 >= i2 + i3 && 32 * ii2 + i3 + 30 >= 32 * ii3 + i2) {
                            for (int i4 = 32 * ii3 + i2 - i3 - 29; i4 < _PB_N - 1; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          } else if (32 * ii3 + i2 >= 32 * ii2 + i3 + 31 && 32 * ii3 >= i3 + 1) {
                            for (int i4 = 32 * ii2 + 2; i4 < _PB_N - 1; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          } else if (i3 >= 32 * ii3) {
                            for (int i4 = 32 * ii2 + 32 * ii3 - i3 + 1; i4 < _PB_N - 1; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                        }
                      }
                    }
                    for (int i2 = max(1, 64 * ii0 + 32 * ii2 - 2 * i0 + 3); i2 <= min(32 * ii2, 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 + 34); i2 += 1) {
                      for (int i3 = max(1, 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 4); i3 <= 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 35; i3 += 1) {
                        for (int i4 = max(64 * ii0 + 64 * ii2 - 2 * i0 - i2 + 4, 64 * ii0 + 64 * ii2 + 32 * ii3 - 2 * i0 - i2 - i3 + 5); i4 < _PB_N - 1; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                    for (int i2 = 32 * ii2 + 1; i2 < _PB_N - 1; i2 += 1) {
                      for (int i3 = max(1, 64 * ii0 + 32 * ii3 - 2 * i0 + 3); i3 <= 64 * ii0 + 32 * ii3 - 2 * i0 + 34; i3 += 1) {
                        for (int i4 = 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i3 + 4; i4 < _PB_N - 1; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                  }
                }
              } else {
                for (int i0 = 32 * ii0 + 1; i0 <= min(_PB_TSTEPS, 32 * ii0 + 32); i0 += 1) {
                  if (_PB_N >= 32 * ii2 + 4 && i0 >= 32 * ii0 + 2 && 32 * ii0 + 16 * ii3 + 16 >= i0 && _PB_N + 2 * i0 >= 64 * ii0 + 32 * ii2 + 9) {
                    if (i0 >= 32 * ii0 + 3) {
                      for (int i2 = 64 * ii0 + 32 * ii2 - 2 * i0 + 4; i2 <= min(32 * ii2 + 1, 64 * ii0 + 32 * ii2 - 2 * i0 + 34); i2 += 1) {
                        if (2 * i0 + i2 >= 64 * ii0 + 32 * ii2 + 7) {
                          if (32 * ii0 + 16 * ii3 + 1 >= i0 && i2 == 32 * ii2 + 1) {
                            for (int i4 = 1; i4 <= 32; i4 += 1) {
                              B[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4])) + A[32 * ii2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 5][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4])) + A[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 3][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4])) + A[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4 - 1]))) + A[32 * ii2 + 1][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4]);
                            }
                          }
                          for (int i3 = max(max(1, 64 * ii0 - 32 * ii2 + 32 * ii3 - 2 * i0 + i2 + 4), 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 5); i3 <= min(32 * ii3 + 1, 32 * ii2 + 32 * ii3 - i2); i3 += 1) {
                            B[i2][i3][1] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][1] - (SCALAR_VAL(2.0) * A[i2][i3][1])) + A[i2 - 1][i3][1])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][1] - (SCALAR_VAL(2.0) * A[i2][i3][1])) + A[i2][i3 - 1][1]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][2] - (SCALAR_VAL(2.0) * A[i2][i3][1])) + A[i2][i3][0]))) + A[i2][i3][1]);
                            if (2 * i0 + i2 + i3 >= 64 * ii0 + 32 * ii2 + 32 * ii3 + 7) {
                              for (int i4 = 2; i4 <= 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 - i3 + 37; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            } else {
                              for (int i4 = 2; i4 <= 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 - i3 + 37; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                          }
                          if (i2 >= 32 * ii2) {
                            for (int i3 = max(1, 32 * ii2 + 32 * ii3 - i2 + 1); i3 <= 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 36; i3 += 1) {
                              if (32 * ii3 + 1 >= i3) {
                                B[i2][i3][1] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][1] - (SCALAR_VAL(2.0) * A[i2][i3][1])) + A[i2 - 1][i3][1])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][1] - (SCALAR_VAL(2.0) * A[i2][i3][1])) + A[i2][i3 - 1][1]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][2] - (SCALAR_VAL(2.0) * A[i2][i3][1])) + A[i2][i3][0]))) + A[i2][i3][1]);
                                for (int i4 = 2; i4 <= 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 - i3 + 37; i4 += 1) {
                                  B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                                }
                              } else {
                                for (int i4 = 1; i4 <= 64 * ii0 + 32 * ii2 - 2 * i0 - i2 + 36; i4 += 1) {
                                  B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                                }
                              }
                            }
                          } else {
                            for (int i3 = 32 * ii3 + 2; i3 <= 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 36; i3 += 1) {
                              for (int i4 = 1; i4 <= 64 * ii0 + 32 * ii2 - 2 * i0 - i2 + 36; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                          }
                        } else {
                          if (i0 == 32 * ii0 + 3 && i2 == 32 * ii2) {
                            for (int i3 = max(1, 32 * ii3 - 1); i3 <= 32 * ii3; i3 += 1) {
                              for (int i4 = 1; i4 <= 32 * ii3 - i3 + 31; i4 += 1) {
                                B[32 * ii2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3][i4 - 1]))) + A[32 * ii2][i3][i4]);
                              }
                            }
                            for (int i3 = 32 * ii3 + 1; i3 <= 32 * ii3 + 2; i3 += 1) {
                              for (int i4 = 1; i4 <= 30; i4 += 1) {
                                B[32 * ii2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3][i4 - 1]))) + A[32 * ii2][i3][i4]);
                              }
                            }
                          } else {
                            for (int i3 = max(1, 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 5); i3 <= 32 * ii3; i3 += 1) {
                              B[i2][i3][1] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][1] - (SCALAR_VAL(2.0) * A[i2][i3][1])) + A[i2 - 1][i3][1])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][1] - (SCALAR_VAL(2.0) * A[i2][i3][1])) + A[i2][i3 - 1][1]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][2] - (SCALAR_VAL(2.0) * A[i2][i3][1])) + A[i2][i3][0]))) + A[i2][i3][1]);
                              for (int i4 = 2; i4 <= 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 - i3 + 37; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                            for (int i3 = 32 * ii3 + 1; i3 <= 32 * ii3 + 2; i3 += 1) {
                              if (i3 == 32 * ii3 + 1) {
                                B[i2][32 * ii3 + 1][1] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii3 + 1][1] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 1][1])) + A[i2 - 1][32 * ii3 + 1][1])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 2][1] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 1][1])) + A[i2][32 * ii3][1]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 1][2] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 1][1])) + A[i2][32 * ii3 + 1][0]))) + A[i2][32 * ii3 + 1][1]);
                              }
                              for (int i4 = 32 * ii3 - i3 + 3; i4 <= 64 * ii0 + 32 * ii2 - 2 * i0 - i2 + 36; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                          }
                          if (2 * i0 + i2 == 64 * ii0 + 32 * ii2 + 6) {
                            for (int i3 = 32 * ii3 + 3; i3 <= 32 * ii3 + 30; i3 += 1) {
                              for (int i4 = 1; i4 <= 30; i4 += 1) {
                                B[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 7][i3][i4] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3][i4])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3][i4])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3][i4])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3][i4 - 1]))) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3][i4]);
                              }
                            }
                          } else {
                            for (int i3 = 32 * ii3 + 3; i3 <= 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 36; i3 += 1) {
                              for (int i4 = 1; i4 <= 128 * ii0 + 64 * ii2 - 4 * i0 - 2 * i2 + 40; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                              if (2 * i0 + i2 == 64 * ii0 + 32 * ii2 + 5) {
                                B[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3][31] = ((((SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3][31] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3][31])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][i3][31])) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3 + 1][31] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3][31])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3 - 1][31]))) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3][32] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3][31])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3][30]))) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3][31]);
                              }
                            }
                          }
                        }
                      }
                    }
                    for (int i2 = 64 * ii0 + 32 * ii2 - 2 * i0 + 35; i2 <= 32 * ii2; i2 += 1) {
                      for (int i3 = max(1, 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 5); i3 <= 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 36; i3 += 1) {
                        for (int i4 = 1; i4 <= 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 - i3 + 37; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    }
                    if (ii3 == 1 && i0 >= 32 * ii0 + 18) {
                      for (int i3 = 1; i3 <= 64 * ii0 - 2 * i0 + 67; i3 += 1) {
                        B[32 * ii2 + 1][i3][1] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3][1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][1])) + A[32 * ii2][i3][1])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3 + 1][1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][1])) + A[32 * ii2 + 1][i3 - 1][1]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][2] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][1])) + A[32 * ii2 + 1][i3][0]))) + A[32 * ii2 + 1][i3][1]);
                        if (2 * i0 + i3 >= 64 * ii0 + 38) {
                          for (int i4 = 2; i4 <= 64 * ii0 - 2 * i0 - i3 + 68; i4 += 1) {
                            B[32 * ii2 + 1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3][i4 - 1]))) + A[32 * ii2 + 1][i3][i4]);
                          }
                        } else {
                          for (int i4 = 2; i4 <= 31; i4 += 1) {
                            B[32 * ii2 + 1][1][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][1][i4])) + A[32 * ii2][1][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][2][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][1][i4])) + A[32 * ii2 + 1][0][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][1][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][1][i4])) + A[32 * ii2 + 1][1][i4 - 1]))) + A[32 * ii2 + 1][1][i4]);
                          }
                        }
                      }
                    } else if (i0 >= 32 * ii0 + 18 && 32 * ii0 + 16 * ii3 + 1 >= i0) {
                      for (int i3 = 64 * ii0 + 32 * ii3 - 2 * i0 + 4; i3 <= 64 * ii0 + 32 * ii3 - 2 * i0 + 35; i3 += 1) {
                        for (int i4 = 1; i4 <= 64 * ii0 + 32 * ii3 - 2 * i0 - i3 + 36; i4 += 1) {
                          B[32 * ii2 + 1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3][i4 - 1]))) + A[32 * ii2 + 1][i3][i4]);
                        }
                      }
                    } else if (i0 == 32 * ii0 + 17) {
                      for (int i3 = 32 * ii3 - 30; i3 <= 32 * ii3 + 1; i3 += 1) {
                        for (int i4 = 1; i4 <= 32 * ii3 - i3 + 2; i4 += 1) {
                          if (i3 + 60 >= 32 * ii3 + 30 * i4 || 32 * ii3 >= i3 + 29 || 1) {
                            B[32 * ii2 + 1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3][i4 - 1]))) + A[32 * ii2 + 1][i3][i4]);
                          }
                        }
                      }
                    }
                    for (int i2 = 32 * ii2 + 2; i2 <= min(min(_PB_N - 2, 32 * ii2 + 31), -64 * ii0 + 32 * ii2 + 2 * i0 - 3); i2 += 1) {
                      for (int i3 = max(1, 64 * ii0 + 32 * ii3 - 2 * i0 + 4); i3 <= min(32 * ii3 - 1, 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 36); i3 += 1) {
                        for (int i4 = 1; i4 <= 64 * ii0 + 32 * ii3 - 2 * i0 - i3 + 36; i4 += 1) {
                          if (2 * i0 + i3 + 2 * i4 >= 64 * ii0 + 32 * ii3 + 8 || 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                      for (int i3 = max(1, 32 * ii3); i3 <= 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 36; i3 += 1) {
                        if (i3 >= 32 * ii3 + 2) {
                          for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 + 35; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        } else if (32 * ii2 + 2 * i0 + i3 >= 64 * ii0 + 32 * ii3 + i2 + 4) {
                          for (int i4 = 1; i4 <= 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 - i3 + 37; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          for (int i4 = 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 - i3 + 38; i4 <= 64 * ii0 + 32 * ii3 - 2 * i0 - i3 + 36; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        } else {
                          for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 + 35; i4 += 1) {
                            B[-64 * ii0 + 32 * ii2 + 2 * i0 - 3][32 * ii3][i4] = ((((SCALAR_VAL(0.125) * ((A[-64 * ii0 + 32 * ii2 + 2 * i0 - 2][32 * ii3][i4] - (SCALAR_VAL(2.0) * A[-64 * ii0 + 32 * ii2 + 2 * i0 - 3][32 * ii3][i4])) + A[-64 * ii0 + 32 * ii2 + 2 * i0 - 4][32 * ii3][i4])) + (SCALAR_VAL(0.125) * ((A[-64 * ii0 + 32 * ii2 + 2 * i0 - 3][32 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[-64 * ii0 + 32 * ii2 + 2 * i0 - 3][32 * ii3][i4])) + A[-64 * ii0 + 32 * ii2 + 2 * i0 - 3][32 * ii3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[-64 * ii0 + 32 * ii2 + 2 * i0 - 3][32 * ii3][i4 + 1] - (SCALAR_VAL(2.0) * A[-64 * ii0 + 32 * ii2 + 2 * i0 - 3][32 * ii3][i4])) + A[-64 * ii0 + 32 * ii2 + 2 * i0 - 3][32 * ii3][i4 - 1]))) + A[-64 * ii0 + 32 * ii2 + 2 * i0 - 3][32 * ii3][i4]);
                          }
                          B[-64 * ii0 + 32 * ii2 + 2 * i0 - 3][32 * ii3][64 * ii0 - 2 * i0 + 36] = ((((SCALAR_VAL(0.125) * ((A[-64 * ii0 + 32 * ii2 + 2 * i0 - 2][32 * ii3][64 * ii0 - 2 * i0 + 36] - (SCALAR_VAL(2.0) * A[-64 * ii0 + 32 * ii2 + 2 * i0 - 3][32 * ii3][64 * ii0 - 2 * i0 + 36])) + A[-64 * ii0 + 32 * ii2 + 2 * i0 - 4][32 * ii3][64 * ii0 - 2 * i0 + 36])) + (SCALAR_VAL(0.125) * ((A[-64 * ii0 + 32 * ii2 + 2 * i0 - 3][32 * ii3 + 1][64 * ii0 - 2 * i0 + 36] - (SCALAR_VAL(2.0) * A[-64 * ii0 + 32 * ii2 + 2 * i0 - 3][32 * ii3][64 * ii0 - 2 * i0 + 36])) + A[-64 * ii0 + 32 * ii2 + 2 * i0 - 3][32 * ii3 - 1][64 * ii0 - 2 * i0 + 36]))) + (SCALAR_VAL(0.125) * ((A[-64 * ii0 + 32 * ii2 + 2 * i0 - 3][32 * ii3][64 * ii0 - 2 * i0 + 37] - (SCALAR_VAL(2.0) * A[-64 * ii0 + 32 * ii2 + 2 * i0 - 3][32 * ii3][64 * ii0 - 2 * i0 + 36])) + A[-64 * ii0 + 32 * ii2 + 2 * i0 - 3][32 * ii3][64 * ii0 - 2 * i0 + 35]))) + A[-64 * ii0 + 32 * ii2 + 2 * i0 - 3][32 * ii3][64 * ii0 - 2 * i0 + 36]);
                        }
                      }
                      for (int i3 = max(1, 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 37); i3 <= min(32 * ii3 + 1, 64 * ii0 + 32 * ii3 - 2 * i0 + 35); i3 += 1) {
                        for (int i4 = 1; i4 <= 64 * ii0 + 32 * ii3 - 2 * i0 - i3 + 36; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                      for (int i3 = max(32 * ii3 + 2, 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 37); i3 <= 64 * ii0 + 32 * ii3 - 2 * i0 + 35; i3 += 1) {
                        for (int i4 = 1; i4 <= 64 * ii0 + 32 * ii2 - 2 * i0 - i2 + 36; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = max(1, 64 * ii0 + 32 * ii2 - 2 * i0 - i2 + 37); i4 <= 64 * ii0 - 2 * i0 + 35; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    }
                    if (32 * ii2 + 34 == _PB_N && i0 >= 32 * ii0 + 18) {
                      for (int i3 = max(1, 64 * ii0 + 32 * ii3 - 2 * i0 + 4); i3 <= 64 * ii0 + 32 * ii3 - 2 * i0 + 35; i3 += 1) {
                        if (2 * i0 + i3 >= 64 * ii0 + 32 * ii3 + 6) {
                          B[_PB_N - 2][i3][1] = ((((SCALAR_VAL(0.125) * ((A[_PB_N - 1][i3][1] - (SCALAR_VAL(2.0) * A[_PB_N - 2][i3][1])) + A[_PB_N - 3][i3][1])) + (SCALAR_VAL(0.125) * ((A[_PB_N - 2][i3 + 1][1] - (SCALAR_VAL(2.0) * A[_PB_N - 2][i3][1])) + A[_PB_N - 2][i3 - 1][1]))) + (SCALAR_VAL(0.125) * ((A[_PB_N - 2][i3][2] - (SCALAR_VAL(2.0) * A[_PB_N - 2][i3][1])) + A[_PB_N - 2][i3][0]))) + A[_PB_N - 2][i3][1]);
                        }
                        for (int i4 = -2 * ii0 - ii3 + (-6 * ii0 - 3 * ii3 + 2 * i0 + i3 - 6) / 29 + 2; i4 <= 64 * ii0 + 32 * ii3 - 2 * i0 - i3 + 36; i4 += 1) {
                          B[_PB_N - 2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[_PB_N - 1][i3][i4] - (SCALAR_VAL(2.0) * A[_PB_N - 2][i3][i4])) + A[_PB_N - 3][i3][i4])) + (SCALAR_VAL(0.125) * ((A[_PB_N - 2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[_PB_N - 2][i3][i4])) + A[_PB_N - 2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[_PB_N - 2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[_PB_N - 2][i3][i4])) + A[_PB_N - 2][i3][i4 - 1]))) + A[_PB_N - 2][i3][i4]);
                        }
                      }
                    } else if (i0 == 32 * ii0 + 2) {
                      for (int i2 = 32 * ii2; i2 <= 32 * ii2 + 2; i2 += 1) {
                        for (int i3 = max(max(1, 32 * ii3), 32 * ii2 + 32 * ii3 - i2 + 1); i3 <= 32 * ii2 + 32 * ii3 - i2 + 32; i3 += 1) {
                          for (int i4 = 1; i4 <= 32 * ii2 - i2 + 32; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (i3 == 32 * ii3) {
                            for (int i4 = 32 * ii2 - i2 + 33; i4 <= 32; i4 += 1) {
                              B[i2][32 * ii3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii3][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][i4])) + A[i2 - 1][32 * ii3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][i4])) + A[i2][32 * ii3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][i4])) + A[i2][32 * ii3][i4 - 1]))) + A[i2][32 * ii3][i4]);
                            }
                          } else if (i2 == 32 * ii2 + 2) {
                            B[32 * ii2 + 2][i3][31] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 3][i3][31] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][i3][31])) + A[32 * ii2 + 1][i3][31])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3 + 1][31] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][i3][31])) + A[32 * ii2 + 2][i3 - 1][31]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3][32] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][i3][31])) + A[32 * ii2 + 2][i3][30]))) + A[32 * ii2 + 2][i3][31]);
                          }
                        }
                        if (i2 == 32 * ii2 + 2) {
                          for (int i4 = 1; i4 <= 31; i4 += 1) {
                            B[32 * ii2 + 2][32 * ii3 + 31][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 3][32 * ii3 + 31][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][32 * ii3 + 31][i4])) + A[32 * ii2 + 1][32 * ii3 + 31][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][32 * ii3 + 32][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][32 * ii3 + 31][i4])) + A[32 * ii2 + 2][32 * ii3 + 30][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][32 * ii3 + 31][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][32 * ii3 + 31][i4])) + A[32 * ii2 + 2][32 * ii3 + 31][i4 - 1]))) + A[32 * ii2 + 2][32 * ii3 + 31][i4]);
                          }
                        }
                      }
                    }
                    for (int i2 = max(-64 * ii0 + 32 * ii2 + 2 * i0 - 2, 64 * ii0 + 32 * ii2 - 2 * i0 + 7); i2 < _PB_N - 1; i2 += 1) {
                      for (int i3 = max(1, 64 * ii0 + 32 * ii3 - 2 * i0 + 4); i3 <= min(32 * ii3 - 1, 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 36); i3 += 1) {
                        for (int i4 = 1; i4 <= 64 * ii0 + 32 * ii3 - 2 * i0 - i3 + 36; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                      for (int i3 = max(1, 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 37); i3 < 32 * ii3; i3 += 1) {
                        if (2 * i0 + i3 >= 64 * ii0 + 32 * ii3 + 6) {
                          B[i2][i3][1] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][1] - (SCALAR_VAL(2.0) * A[i2][i3][1])) + A[i2 - 1][i3][1])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][1] - (SCALAR_VAL(2.0) * A[i2][i3][1])) + A[i2][i3 - 1][1]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][2] - (SCALAR_VAL(2.0) * A[i2][i3][1])) + A[i2][i3][0]))) + A[i2][i3][1]);
                        }
                        for (int i4 = -2 * ii0 - ii3 + (-8 * ii0 - 4 * ii3 + 2 * i0 + i3 - 6) / 28 + 2; i4 <= 64 * ii0 + 32 * ii3 - 2 * i0 - i3 + 36; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                      if (ii3 >= 1 && 64 * ii0 + 32 * ii2 + 36 >= 2 * i0 + i2) {
                        for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 + 36; i4 += 1) {
                          B[i2][32 * ii3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii3][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][i4])) + A[i2 - 1][32 * ii3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][i4])) + A[i2][32 * ii3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][i4])) + A[i2][32 * ii3][i4 - 1]))) + A[i2][32 * ii3][i4]);
                        }
                      }
                      for (int i3 = 32 * ii3 + 1; i3 <= min(-64 * ii0 + 32 * ii3 + 2 * i0 - 5, 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 36); i3 += 1) {
                        for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 + 35; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                      for (int i3 = max(32 * ii3 + 1, -64 * ii0 + 32 * ii3 + 2 * i0 - 4); i3 <= 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 36; i3 += 1) {
                        if (i0 >= 32 * ii0 + 3) {
                          for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 + 35; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        } else {
                          for (int i4 = 1; i4 <= min(31, 32 * ii3 - i3 + 33); i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          for (int i4 = 32 * ii3 - i3 + 34; i4 <= 31; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                      for (int i3 = max(max(1, 32 * ii3), 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 37); i3 <= 64 * ii0 + 32 * ii3 - 2 * i0 + 35; i3 += 1) {
                        if (i0 == 32 * ii0 + 2) {
                          for (int i4 = 1; i4 <= 64 * ii3 - 2 * i3 + 32; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          for (int i4 = max(1, 64 * ii3 - 2 * i3 + 33); i4 <= min(31, 32 * ii3 - i3 + 33); i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        } else if (i3 >= 32 * ii3 + 2) {
                          for (int i4 = 1; i4 <= 64 * ii0 + 32 * ii3 - 2 * i0 - i3 + 37; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                        for (int i4 = 64 * ii0 + 32 * ii3 - 2 * i0 - i3 + 38; i4 <= 64 * ii0 - 2 * i0 + 35; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        if (i0 >= 32 * ii0 + 3 && 32 * ii3 + 1 >= i3) {
                          for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 + 35; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (i3 == 32 * ii3) {
                            B[i2][32 * ii3][64 * ii0 - 2 * i0 + 36] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii3][64 * ii0 - 2 * i0 + 36] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][64 * ii0 - 2 * i0 + 36])) + A[i2 - 1][32 * ii3][64 * ii0 - 2 * i0 + 36])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 1][64 * ii0 - 2 * i0 + 36] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][64 * ii0 - 2 * i0 + 36])) + A[i2][32 * ii3 - 1][64 * ii0 - 2 * i0 + 36]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3][64 * ii0 - 2 * i0 + 37] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][64 * ii0 - 2 * i0 + 36])) + A[i2][32 * ii3][64 * ii0 - 2 * i0 + 35]))) + A[i2][32 * ii3][64 * ii0 - 2 * i0 + 36]);
                          }
                        }
                      }
                    }
                  } else if (ii3 == 0 && i0 >= 32 * ii0 + 17) {
                    for (int i2 = max(1, 64 * ii0 + 32 * ii2 - 2 * i0 + 4); i2 <= 64 * ii0 + 32 * ii2 - 2 * i0 + 35; i2 += 1) {
                      for (int i4 = 1; i4 <= 64 * ii0 + 32 * ii2 - 2 * i0 - i2 + 36; i4 += 1) {
                        B[i2][1][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][1][i4] - (SCALAR_VAL(2.0) * A[i2][1][i4])) + A[i2 - 1][1][i4])) + (SCALAR_VAL(0.125) * ((A[i2][2][i4] - (SCALAR_VAL(2.0) * A[i2][1][i4])) + A[i2][0][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][1][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][1][i4])) + A[i2][1][i4 - 1]))) + A[i2][1][i4]);
                      }
                      if (64 * ii0 + 32 * ii2 + 6 >= 2 * i0 + i2) {
                        for (int i4 = 1; i4 <= 64 * ii0 + 32 * ii2 - 2 * i0 - i2 + 36; i4 += 1) {
                          B[i2][2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][2][i4] - (SCALAR_VAL(2.0) * A[i2][2][i4])) + A[i2 - 1][2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][3][i4] - (SCALAR_VAL(2.0) * A[i2][2][i4])) + A[i2][1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][2][i4])) + A[i2][2][i4 - 1]))) + A[i2][2][i4]);
                        }
                        if (64 * ii0 + 32 * ii2 + 5 >= 2 * i0 + i2) {
                          for (int i3 = 3; i3 <= 64 * ii0 + 32 * ii2 - 2 * i0 - i2 + 36; i3 += 1) {
                            for (int i4 = 1; i4 <= 128 * ii0 + 64 * ii2 - 4 * i0 - 2 * i2 + 40; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            if (2 * i0 + i2 == 64 * ii0 + 32 * ii2 + 5) {
                              B[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3][31] = ((((SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 6][i3][31] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3][31])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 4][i3][31])) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3 + 1][31] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3][31])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3 - 1][31]))) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3][32] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3][31])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3][30]))) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 5][i3][31]);
                            }
                          }
                        }
                      }
                      if (2 * i0 + i2 >= 64 * ii0 + 32 * ii2 + 6 && 32 * ii2 >= i2 + 1) {
                        for (int i3 = max(2, 64 * ii0 + 32 * ii2 - 2 * i0 - i2 + 9); i3 <= 64 * ii0 + 32 * ii2 - 2 * i0 - i2 + 36; i3 += 1) {
                          for (int i4 = 1; i4 <= 64 * ii0 + 32 * ii2 - 2 * i0 - i2 + 36; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      } else if (i0 == 32 * ii0 + 17 && i2 == 32 * ii2) {
                        for (int i4 = 1; i4 <= 2; i4 += 1) {
                          B[32 * ii2][2][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][2][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][2][i4])) + A[32 * ii2 - 1][2][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2][3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][2][i4])) + A[32 * ii2][1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2][2][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2][2][i4])) + A[32 * ii2][2][i4 - 1]))) + A[32 * ii2][2][i4]);
                        }
                      }
                    }
                    if (i0 == 32 * ii0 + 17) {
                      for (int i2 = 32 * ii2 + 2; i2 < _PB_N - 1; i2 += 1) {
                        B[i2][1][1] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][1][1] - (SCALAR_VAL(2.0) * A[i2][1][1])) + A[i2 - 1][1][1])) + (SCALAR_VAL(0.125) * ((A[i2][2][1] - (SCALAR_VAL(2.0) * A[i2][1][1])) + A[i2][0][1]))) + (SCALAR_VAL(0.125) * ((A[i2][1][2] - (SCALAR_VAL(2.0) * A[i2][1][1])) + A[i2][1][0]))) + A[i2][1][1]);
                      }
                    }
                  } else if (32 * ii2 + 4 >= _PB_N && i0 == 32 * ii0 + 2) {
                    for (int i2 = 32 * ii2; i2 < _PB_N - 1; i2 += 1) {
                      for (int i3 = max(max(1, 32 * ii3), 32 * ii2 + 32 * ii3 - i2 + 1); i3 <= 32 * ii2 + 32 * ii3 - i2 + 32; i3 += 1) {
                        for (int i4 = 1; i4 <= 32 * ii2 - i2 + 32; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        if (i3 == 32 * ii3) {
                          for (int i4 = 32 * ii2 - i2 + 33; i4 <= 32; i4 += 1) {
                            B[i2][32 * ii3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii3][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][i4])) + A[i2 - 1][32 * ii3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][i4])) + A[i2][32 * ii3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][i4])) + A[i2][32 * ii3][i4 - 1]))) + A[i2][32 * ii3][i4]);
                          }
                        } else if (32 * ii2 + 4 == _PB_N && i2 + 2 == _PB_N) {
                          B[_PB_N - 2][i3][31] = ((((SCALAR_VAL(0.125) * ((A[_PB_N - 1][i3][31] - (SCALAR_VAL(2.0) * A[_PB_N - 2][i3][31])) + A[_PB_N - 3][i3][31])) + (SCALAR_VAL(0.125) * ((A[_PB_N - 2][i3 + 1][31] - (SCALAR_VAL(2.0) * A[_PB_N - 2][i3][31])) + A[_PB_N - 2][i3 - 1][31]))) + (SCALAR_VAL(0.125) * ((A[_PB_N - 2][i3][32] - (SCALAR_VAL(2.0) * A[_PB_N - 2][i3][31])) + A[_PB_N - 2][i3][30]))) + A[_PB_N - 2][i3][31]);
                        }
                      }
                      if (32 * ii2 + 4 == _PB_N && i2 + 2 == _PB_N) {
                        for (int i4 = 1; i4 <= 31; i4 += 1) {
                          B[_PB_N - 2][32 * ii3 + 31][i4] = ((((SCALAR_VAL(0.125) * ((A[_PB_N - 1][32 * ii3 + 31][i4] - (SCALAR_VAL(2.0) * A[_PB_N - 2][32 * ii3 + 31][i4])) + A[_PB_N - 3][32 * ii3 + 31][i4])) + (SCALAR_VAL(0.125) * ((A[_PB_N - 2][32 * ii3 + 32][i4] - (SCALAR_VAL(2.0) * A[_PB_N - 2][32 * ii3 + 31][i4])) + A[_PB_N - 2][32 * ii3 + 30][i4]))) + (SCALAR_VAL(0.125) * ((A[_PB_N - 2][32 * ii3 + 31][i4 + 1] - (SCALAR_VAL(2.0) * A[_PB_N - 2][32 * ii3 + 31][i4])) + A[_PB_N - 2][32 * ii3 + 31][i4 - 1]))) + A[_PB_N - 2][32 * ii3 + 31][i4]);
                        }
                      }
                    }
                  } else if (32 * ii2 + 3 == _PB_N && i0 >= 32 * ii0 + 3) {
                    for (int i2 = _PB_N + 64 * ii0 - 2 * i0 + 1; i2 < _PB_N - 1; i2 += 1) {
                      if (_PB_N >= i2 + 4) {
                        for (int i3 = max(1, _PB_N + 64 * ii0 + 32 * ii3 - 2 * i0 - i2 + 2); i3 <= min(32 * ii3 + 1, _PB_N + 64 * ii0 + 32 * ii3 - 2 * i0 - i2 + 33); i3 += 1) {
                          B[i2][i3][1] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][1] - (SCALAR_VAL(2.0) * A[i2][i3][1])) + A[i2 - 1][i3][1])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][1] - (SCALAR_VAL(2.0) * A[i2][i3][1])) + A[i2][i3 - 1][1]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][2] - (SCALAR_VAL(2.0) * A[i2][i3][1])) + A[i2][i3][0]))) + A[i2][i3][1]);
                          if (_PB_N + 64 * ii0 + 3 >= 2 * i0 + i2 && i3 == 32 * ii3 + 1) {
                            for (int i4 = 2; i4 <= _PB_N + 64 * ii0 - 2 * i0 - i2 + 33; i4 += 1) {
                              B[i2][32 * ii3 + 1][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 1][i4])) + A[i2 - 1][32 * ii3 + 1][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 2][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 1][i4])) + A[i2][32 * ii3][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 1][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 1][i4])) + A[i2][32 * ii3 + 1][i4 - 1]))) + A[i2][32 * ii3 + 1][i4]);
                            }
                          } else if (2 * i0 + i2 >= _PB_N + 64 * ii0 + 4 && 2 * i0 + i2 + i3 >= _PB_N + 64 * ii0 + 32 * ii3 + 4) {
                            for (int i4 = 2; i4 <= _PB_N + 64 * ii0 + 32 * ii3 - 2 * i0 - i2 - i3 + 34; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          } else {
                            for (int i4 = 2; i4 <= _PB_N + 64 * ii0 + 32 * ii3 - 2 * i0 - i2 - i3 + 34; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                        }
                        if (_PB_N + 64 * ii0 + 3 >= 2 * i0 + i2) {
                          for (int i4 = 1; i4 <= _PB_N + 64 * ii0 - 2 * i0 - i2 + 33; i4 += 1) {
                            B[i2][32 * ii3 + 2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii3 + 2][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 2][i4])) + A[i2 - 1][32 * ii3 + 2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 3][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 2][i4])) + A[i2][32 * ii3 + 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii3 + 2][i4])) + A[i2][32 * ii3 + 2][i4 - 1]))) + A[i2][32 * ii3 + 2][i4]);
                          }
                        }
                        if (2 * i0 + i2 >= _PB_N + 64 * ii0 + 3) {
                          for (int i3 = max(32 * ii3 + 2, _PB_N + 64 * ii0 + 32 * ii3 - 2 * i0 - i2 + 6); i3 <= _PB_N + 64 * ii0 + 32 * ii3 - 2 * i0 - i2 + 33; i3 += 1) {
                            for (int i4 = 1; i4 <= _PB_N + 64 * ii0 - 2 * i0 - i2 + 33; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                        } else {
                          for (int i3 = 32 * ii3 + 3; i3 <= _PB_N + 64 * ii0 + 32 * ii3 - 2 * i0 - i2 + 33; i3 += 1) {
                            for (int i4 = 1; i4 <= 2 * _PB_N + 128 * ii0 - 4 * i0 - 2 * i2 + 34; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            if (2 * i0 + i2 == _PB_N + 64 * ii0 + 2) {
                              B[_PB_N + 64 * ii0 - 2 * i0 + 2][i3][31] = ((((SCALAR_VAL(0.125) * ((A[_PB_N + 64 * ii0 - 2 * i0 + 3][i3][31] - (SCALAR_VAL(2.0) * A[_PB_N + 64 * ii0 - 2 * i0 + 2][i3][31])) + A[_PB_N + 64 * ii0 - 2 * i0 + 1][i3][31])) + (SCALAR_VAL(0.125) * ((A[_PB_N + 64 * ii0 - 2 * i0 + 2][i3 + 1][31] - (SCALAR_VAL(2.0) * A[_PB_N + 64 * ii0 - 2 * i0 + 2][i3][31])) + A[_PB_N + 64 * ii0 - 2 * i0 + 2][i3 - 1][31]))) + (SCALAR_VAL(0.125) * ((A[_PB_N + 64 * ii0 - 2 * i0 + 2][i3][32] - (SCALAR_VAL(2.0) * A[_PB_N + 64 * ii0 - 2 * i0 + 2][i3][31])) + A[_PB_N + 64 * ii0 - 2 * i0 + 2][i3][30]))) + A[_PB_N + 64 * ii0 - 2 * i0 + 2][i3][31]);
                            }
                          }
                        }
                      } else if (2 * i0 + i2 >= _PB_N + 64 * ii0 + 4) {
                        if (32 * ii0 + 17 >= i0 && i2 + 3 == _PB_N) {
                          for (int i3 = max(1, 64 * ii0 + 32 * ii3 - 2 * i0 + 5); i3 <= 32 * ii3; i3 += 1) {
                            for (int i4 = 1; i4 <= 64 * ii0 + 32 * ii3 - 2 * i0 - i3 + 37; i4 += 1) {
                              B[_PB_N - 3][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[_PB_N - 2][i3][i4] - (SCALAR_VAL(2.0) * A[_PB_N - 3][i3][i4])) + A[_PB_N - 4][i3][i4])) + (SCALAR_VAL(0.125) * ((A[_PB_N - 3][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[_PB_N - 3][i3][i4])) + A[_PB_N - 3][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[_PB_N - 3][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[_PB_N - 3][i3][i4])) + A[_PB_N - 3][i3][i4 - 1]))) + A[_PB_N - 3][i3][i4]);
                            }
                          }
                        } else {
                          if (32 * ii0 + 16 * ii3 + 1 >= i0 && i2 + 2 == _PB_N) {
                            for (int i4 = 1; i4 <= 32; i4 += 1) {
                              B[_PB_N - 2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4] = ((((SCALAR_VAL(0.125) * ((A[_PB_N - 1][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4] - (SCALAR_VAL(2.0) * A[_PB_N - 2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4])) + A[_PB_N - 3][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4])) + (SCALAR_VAL(0.125) * ((A[_PB_N - 2][64 * ii0 + 32 * ii3 - 2 * i0 + 5][i4] - (SCALAR_VAL(2.0) * A[_PB_N - 2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4])) + A[_PB_N - 2][64 * ii0 + 32 * ii3 - 2 * i0 + 3][i4]))) + (SCALAR_VAL(0.125) * ((A[_PB_N - 2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4 + 1] - (SCALAR_VAL(2.0) * A[_PB_N - 2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4])) + A[_PB_N - 2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4 - 1]))) + A[_PB_N - 2][64 * ii0 + 32 * ii3 - 2 * i0 + 4][i4]);
                            }
                          }
                          if (i2 + 2 == _PB_N) {
                            for (int i3 = max(1, 64 * ii0 + 32 * ii3 - 2 * i0 + 5); i3 <= min(32 * ii3 - 1, 64 * ii0 + 32 * ii3 - 2 * i0 + 35); i3 += 1) {
                              for (int i4 = 1; i4 <= 64 * ii0 + 32 * ii3 - 2 * i0 - i3 + 36; i4 += 1) {
                                B[_PB_N - 2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[_PB_N - 1][i3][i4] - (SCALAR_VAL(2.0) * A[_PB_N - 2][i3][i4])) + A[_PB_N - 3][i3][i4])) + (SCALAR_VAL(0.125) * ((A[_PB_N - 2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[_PB_N - 2][i3][i4])) + A[_PB_N - 2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[_PB_N - 2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[_PB_N - 2][i3][i4])) + A[_PB_N - 2][i3][i4 - 1]))) + A[_PB_N - 2][i3][i4]);
                              }
                            }
                          }
                        }
                        for (int i3 = max(1, _PB_N + 32 * ii3 - i2 - 2); i3 <= _PB_N + 64 * ii0 + 32 * ii3 - 2 * i0 - i2 + 33; i3 += 1) {
                          if (32 * ii3 + 1 >= i3) {
                            for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 + 35; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            if (i2 + i3 + 2 == _PB_N + 32 * ii3) {
                              B[i2][_PB_N + 32 * ii3 - i2 - 2][64 * ii0 - 2 * i0 + 36] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][_PB_N + 32 * ii3 - i2 - 2][64 * ii0 - 2 * i0 + 36] - (SCALAR_VAL(2.0) * A[i2][_PB_N + 32 * ii3 - i2 - 2][64 * ii0 - 2 * i0 + 36])) + A[i2 - 1][_PB_N + 32 * ii3 - i2 - 2][64 * ii0 - 2 * i0 + 36])) + (SCALAR_VAL(0.125) * ((A[i2][_PB_N + 32 * ii3 - i2 - 1][64 * ii0 - 2 * i0 + 36] - (SCALAR_VAL(2.0) * A[i2][_PB_N + 32 * ii3 - i2 - 2][64 * ii0 - 2 * i0 + 36])) + A[i2][_PB_N + 32 * ii3 - i2 - 3][64 * ii0 - 2 * i0 + 36]))) + (SCALAR_VAL(0.125) * ((A[i2][_PB_N + 32 * ii3 - i2 - 2][64 * ii0 - 2 * i0 + 37] - (SCALAR_VAL(2.0) * A[i2][_PB_N + 32 * ii3 - i2 - 2][64 * ii0 - 2 * i0 + 36])) + A[i2][_PB_N + 32 * ii3 - i2 - 2][64 * ii0 - 2 * i0 + 35]))) + A[i2][_PB_N + 32 * ii3 - i2 - 2][64 * ii0 - 2 * i0 + 36]);
                            }
                          } else {
                            for (int i4 = 1; i4 <= _PB_N + 64 * ii0 - 2 * i0 - i2 + 33; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                        }
                        if (i0 >= 32 * ii0 + 18 && i2 + 3 == _PB_N) {
                          for (int i3 = max(1, 64 * ii0 + 32 * ii3 - 2 * i0 + 5); i3 <= 64 * ii0 + 32 * ii3 - 2 * i0 + 36; i3 += 1) {
                            for (int i4 = 1; i4 <= 64 * ii0 + 32 * ii3 - 2 * i0 - i3 + 37; i4 += 1) {
                              B[_PB_N - 3][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[_PB_N - 2][i3][i4] - (SCALAR_VAL(2.0) * A[_PB_N - 3][i3][i4])) + A[_PB_N - 4][i3][i4])) + (SCALAR_VAL(0.125) * ((A[_PB_N - 3][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[_PB_N - 3][i3][i4])) + A[_PB_N - 3][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[_PB_N - 3][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[_PB_N - 3][i3][i4])) + A[_PB_N - 3][i3][i4 - 1]))) + A[_PB_N - 3][i3][i4]);
                            }
                          }
                        }
                      } else {
                        for (int i3 = max(1, 32 * ii3 - 1); i3 <= 32 * ii3; i3 += 1) {
                          for (int i4 = 1; i4 <= 32 * ii3 - i3 + 31; i4 += 1) {
                            B[_PB_N - 3][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[_PB_N - 2][i3][i4] - (SCALAR_VAL(2.0) * A[_PB_N - 3][i3][i4])) + A[_PB_N - 4][i3][i4])) + (SCALAR_VAL(0.125) * ((A[_PB_N - 3][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[_PB_N - 3][i3][i4])) + A[_PB_N - 3][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[_PB_N - 3][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[_PB_N - 3][i3][i4])) + A[_PB_N - 3][i3][i4 - 1]))) + A[_PB_N - 3][i3][i4]);
                          }
                        }
                        for (int i3 = 32 * ii3 + 1; i3 <= 32 * ii3 + 30; i3 += 1) {
                          for (int i4 = 1; i4 <= 30; i4 += 1) {
                            B[_PB_N - 3][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[_PB_N - 2][i3][i4] - (SCALAR_VAL(2.0) * A[_PB_N - 3][i3][i4])) + A[_PB_N - 4][i3][i4])) + (SCALAR_VAL(0.125) * ((A[_PB_N - 3][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[_PB_N - 3][i3][i4])) + A[_PB_N - 3][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[_PB_N - 3][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[_PB_N - 3][i3][i4])) + A[_PB_N - 3][i3][i4 - 1]))) + A[_PB_N - 3][i3][i4]);
                          }
                        }
                      }
                    }
                  }
                  for (int i2 = max(1, 64 * ii0 + 32 * ii2 - 2 * i0 + 3); i2 <= min(32 * ii2, 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 + 34); i2 += 1) {
                    for (int i3 = max(1, 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 4); i3 <= 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 35; i3 += 1) {
                      for (int i4 = 1; i4 <= 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 - i3 + 36; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                      for (int i4 = 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 - i3 + 37; i4 <= 64 * ii0 + 32 * ii2 - 2 * i0 - i2 + 35; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                  }
                  for (int i2 = 32 * ii2 + 1; i2 < _PB_N - 1; i2 += 1) {
                    for (int i3 = max(1, 64 * ii0 + 32 * ii3 - 2 * i0 + 3); i3 <= 64 * ii0 + 32 * ii3 - 2 * i0 + 34; i3 += 1) {
                      if (i3 >= 32 * ii3 + 1) {
                        for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 + 34; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      } else {
                        for (int i4 = 1; i4 <= 64 * ii0 + 32 * ii3 - 2 * i0 - i3 + 35; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                  }
                }
              }
            }
          } else {
            for (int ii4 = 0; ii4 < ii2; ii4 += 1) {
              for (int i0 = 32 * ii0 + 1; i0 <= min(min(_PB_TSTEPS, 32 * ii0 + 32), 32 * ii0 + 16 * ii4 + 16); i0 += 1) {
                if (32 * ii0 + 10 >= i0 && 32 * ii2 + 2 * i0 >= _PB_N + 64 * ii0 + 3) {
                  for (int i2 = 64 * ii0 + 32 * ii2 - 2 * i0 + 4; i2 < _PB_N - 1; i2 += 1) {
                    for (int i3 = max(64 * ii0 + 32 * ii2 - 2 * i0 + 4, 64 * ii0 + 64 * ii2 - 2 * i0 - i2 + 5); i3 < _PB_N - 1; i3 += 1) {
                      if (i3 >= 32 * ii2 + 1) {
                        for (int i4 = max(max(1, 64 * ii0 + 32 * ii4 - 2 * i0 + 4), 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 5); i4 <= 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 36; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 37; i4 <= min(64 * ii0 + 32 * ii4 - 2 * i0 + 35, 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i3 + 37); i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = max(64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 37, 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i3 + 38); i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 + 35; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      } else {
                        for (int i4 = max(max(1, 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i3 + 5), 64 * ii0 + 64 * ii2 + 32 * ii4 - 2 * i0 - i2 - i3 + 6); i4 <= min(32 * ii4 + 1, 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i3 + 6); i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        if (2 * i0 + i3 == 64 * ii0 + 32 * ii2 + 4) {
                          for (int i4 = 32 * ii4 + 2; i4 <= 32 * ii4 + 32; i4 += 1) {
                            B[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 4][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii0 + 32 * ii2 - 2 * i0 + 4][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 4][i4])) + A[i2 - 1][64 * ii0 + 32 * ii2 - 2 * i0 + 4][i4])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 5][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 4][i4])) + A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 3][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 4][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 4][i4])) + A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 4][i4 - 1]))) + A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 4][i4]);
                          }
                        } else {
                          for (int i4 = max(max(1, 64 * ii0 + 64 * ii2 + 32 * ii4 - 2 * i0 - i2 - i3 + 6), 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i3 + 7); i4 <= min(64 * ii0 + 32 * ii4 - 2 * i0 + 35, 64 * ii0 + 64 * ii2 + 32 * ii4 - 2 * i0 - i2 - i3 + 37); i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                        for (int i4 = 64 * ii0 + 64 * ii2 + 32 * ii4 - 2 * i0 - i2 - i3 + 38; i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 + 35; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        if (2 * i0 + i3 >= 64 * ii0 + 32 * ii2 + 5 && 32 * ii2 >= i3 + 1) {
                          for (int i4 = 64 * ii0 + 32 * ii4 - 2 * i0 + 36; i4 <= 64 * ii0 + 64 * ii2 + 32 * ii4 - 2 * i0 - i2 - i3 + 37; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          for (int i4 = max(64 * ii0 + 32 * ii4 - 2 * i0 + 36, 64 * ii0 + 64 * ii2 + 32 * ii4 - 2 * i0 - i2 - i3 + 38); i4 <= 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i3 + 36; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        } else if (32 * ii2 >= i2 && i3 == 32 * ii2) {
                          for (int i4 = 64 * ii0 + 32 * ii4 - 2 * i0 + 36; i4 <= 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 37; i4 += 1) {
                            B[i2][32 * ii2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii2][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2 - 1][32 * ii2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii2 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2][32 * ii2 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2][32 * ii2][i4 - 1]))) + A[i2][32 * ii2][i4]);
                          }
                        } else if (i3 == 32 * ii2) {
                          B[i2][32 * ii2][64 * ii0 + 32 * ii4 - 2 * i0 + 36] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii2][64 * ii0 + 32 * ii4 - 2 * i0 + 36] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][64 * ii0 + 32 * ii4 - 2 * i0 + 36])) + A[i2 - 1][32 * ii2][64 * ii0 + 32 * ii4 - 2 * i0 + 36])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii2 + 1][64 * ii0 + 32 * ii4 - 2 * i0 + 36] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][64 * ii0 + 32 * ii4 - 2 * i0 + 36])) + A[i2][32 * ii2 - 1][64 * ii0 + 32 * ii4 - 2 * i0 + 36]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii2][64 * ii0 + 32 * ii4 - 2 * i0 + 37] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][64 * ii0 + 32 * ii4 - 2 * i0 + 36])) + A[i2][32 * ii2][64 * ii0 + 32 * ii4 - 2 * i0 + 35]))) + A[i2][32 * ii2][64 * ii0 + 32 * ii4 - 2 * i0 + 36]);
                        }
                      }
                    }
                  }
                } else if (i0 >= 32 * ii0 + 11) {
                  for (int i2 = 64 * ii0 + 32 * ii2 - 2 * i0 + 4; i2 < _PB_N - 1; i2 += 1) {
                    if (i2 >= 32 * ii2) {
                      for (int i3 = max(64 * ii0 + 32 * ii2 - 2 * i0 + 4, 64 * ii0 + 64 * ii2 - 2 * i0 - i2 + 5); i3 <= min(32 * ii2, 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 + 5); i3 += 1) {
                        for (int i4 = max(max(1, 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i3 + 5), 64 * ii0 + 64 * ii2 + 32 * ii4 - 2 * i0 - i2 - i3 + 6); i4 <= min(32 * ii4 + 1, 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i3 + 6); i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        if (i2 >= 32 * ii2 + 1 && i3 == 32 * ii2) {
                          for (int i4 = 64 * ii0 + 32 * ii4 - 2 * i0 + 7; i4 <= min(32 * ii4, 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 37); i4 += 1) {
                            B[i2][32 * ii2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii2][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2 - 1][32 * ii2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii2 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2][32 * ii2 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2][32 * ii2][i4 - 1]))) + A[i2][32 * ii2][i4]);
                          }
                        } else if (i2 == 32 * ii2 && i3 == 32 * ii2) {
                          for (int i4 = 64 * ii0 + 32 * ii4 - 2 * i0 + 7; i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 + 37; i4 += 1) {
                            B[32 * ii2][32 * ii2][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][32 * ii2][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][32 * ii2][i4])) + A[32 * ii2 - 1][32 * ii2][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2][32 * ii2 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][32 * ii2][i4])) + A[32 * ii2][32 * ii2 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2][32 * ii2][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2][32 * ii2][i4])) + A[32 * ii2][32 * ii2][i4 - 1]))) + A[32 * ii2][32 * ii2][i4]);
                          }
                        } else if (2 * i0 + i3 >= 64 * ii0 + 32 * ii2 + 5) {
                          for (int i4 = 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i3 + 7; i4 <= 64 * ii0 + 64 * ii2 + 32 * ii4 - 2 * i0 - i2 - i3 + 37; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        } else {
                          for (int i4 = 32 * ii4 + 2; i4 <= 32 * ii4 + 32; i4 += 1) {
                            B[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 4][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii0 + 32 * ii2 - 2 * i0 + 4][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 4][i4])) + A[i2 - 1][64 * ii0 + 32 * ii2 - 2 * i0 + 4][i4])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 5][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 4][i4])) + A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 3][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 4][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 4][i4])) + A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 4][i4 - 1]))) + A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 4][i4]);
                          }
                        }
                        if (i3 == 32 * ii2) {
                          for (int i4 = max(64 * ii0 + 32 * ii4 - 2 * i0 + 7, 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 38); i4 <= min(32 * ii4, 64 * ii0 + 32 * ii4 - 2 * i0 + 36); i4 += 1) {
                            B[i2][32 * ii2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii2][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2 - 1][32 * ii2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii2 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2][32 * ii2 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2][32 * ii2][i4 - 1]))) + A[i2][32 * ii2][i4]);
                          }
                          if (i2 >= 32 * ii2 + 1) {
                            for (int i4 = 32 * ii4 + 1; i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 + 36; i4 += 1) {
                              B[i2][32 * ii2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii2][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2 - 1][32 * ii2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii2 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2][32 * ii2 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2][32 * ii2][i4 - 1]))) + A[i2][32 * ii2][i4]);
                            }
                          }
                        } else if (2 * i0 + i3 >= 64 * ii0 + 32 * ii2 + 5) {
                          for (int i4 = max(64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i3 + 7, 64 * ii0 + 64 * ii2 + 32 * ii4 - 2 * i0 - i2 - i3 + 38); i4 <= 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i3 + 36; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                      for (int i3 = 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 + 6; i3 <= 32 * ii2; i3 += 1) {
                        if (32 * ii2 >= i3 + 1) {
                          for (int i4 = 1; i4 <= 64 * ii0 + 64 * ii2 + 32 * ii4 - 2 * i0 - i2 - i3 + 37; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          for (int i4 = max(1, 64 * ii0 + 64 * ii2 + 32 * ii4 - 2 * i0 - i2 - i3 + 38); i4 <= 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i3 + 36; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        } else {
                          if (i2 >= 32 * ii2 + 1) {
                            for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 + 36; i4 += 1) {
                              B[i2][32 * ii2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii2][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2 - 1][32 * ii2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii2 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2][32 * ii2 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2][32 * ii2][i4 - 1]))) + A[i2][32 * ii2][i4]);
                            }
                          }
                          if (i2 == 32 * ii2) {
                            for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 + 36; i4 += 1) {
                              B[32 * ii2][32 * ii2][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][32 * ii2][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][32 * ii2][i4])) + A[32 * ii2 - 1][32 * ii2][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2][32 * ii2 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][32 * ii2][i4])) + A[32 * ii2][32 * ii2 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2][32 * ii2][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2][32 * ii2][i4])) + A[32 * ii2][32 * ii2][i4 - 1]))) + A[32 * ii2][32 * ii2][i4]);
                            }
                          }
                          for (int i4 = max(1, 64 * ii0 - 2 * i0 + 37); i4 <= 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 37; i4 += 1) {
                            B[i2][32 * ii2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii2][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2 - 1][32 * ii2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii2 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2][32 * ii2 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2][32 * ii2][i4 - 1]))) + A[i2][32 * ii2][i4]);
                          }
                          if (ii4 == 1) {
                            for (int i4 = max(1, 64 * ii0 + 32 * ii2 - 2 * i0 - i2 + 70); i4 <= 64 * ii0 - 2 * i0 + 68; i4 += 1) {
                              B[i2][32 * ii2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii2][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2 - 1][32 * ii2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii2 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2][32 * ii2 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2][32 * ii2][i4 - 1]))) + A[i2][32 * ii2][i4]);
                            }
                          }
                        }
                      }
                      if (ii4 == 0 && 2 * i0 + i2 >= 64 * ii0 + 32 * ii2 + 36 && 32 * ii2 + 30 >= i2) {
                        for (int i3 = 32 * ii2 + 1; i3 <= 32 * ii2 + 2; i3 += 1) {
                          for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 + 35; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                        for (int i3 = 32 * ii2 + 3; i3 < _PB_N - 1; i3 += 1) {
                          for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 + 35; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      } else {
                        if (ii4 == 1 && 2 * i0 + i2 >= 64 * ii0 + 32 * ii2 + 68 && 32 * ii2 + 30 >= i2) {
                          for (int i3 = 32 * ii2 + 1; i3 < _PB_N - 1; i3 += 1) {
                            for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 + 67; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                        } else if (32 * ii0 + 17 >= i0 && i2 >= 32 * ii2 + 31) {
                          for (int i3 = 32 * ii2 + 1; i3 <= 32 * ii2 + 2; i3 += 1) {
                            for (int i4 = max(1, 64 * ii0 + 32 * ii4 - 2 * i0 + 4); i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 + 35; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                        } else if (i0 >= 32 * ii0 + 18 && i2 >= 32 * ii2 + 31) {
                          for (int i3 = 32 * ii2 + 1; i3 < _PB_N - 1; i3 += 1) {
                            for (int i4 = max(1, 64 * ii0 + 32 * ii4 - 2 * i0 + 4); i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 + 35; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                        }
                        if (32 * ii0 + 17 >= i0 && i2 >= 32 * ii2 + 31) {
                          for (int i3 = 32 * ii2 + 3; i3 < _PB_N - 1; i3 += 1) {
                            for (int i4 = max(1, 64 * ii0 + 32 * ii4 - 2 * i0 + 4); i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 + 35; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                        }
                      }
                    } else {
                      for (int i3 = 64 * ii0 + 64 * ii2 - 2 * i0 - i2 + 5; i3 <= 32 * ii2; i3 += 1) {
                        for (int i4 = max(1, 64 * ii0 + 64 * ii2 + 32 * ii4 - 2 * i0 - i2 - i3 + 6); i4 <= 64 * ii0 + 64 * ii2 + 32 * ii4 - 2 * i0 - i2 - i3 + 37; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    }
                    if (32 * ii2 + 30 >= i2 && 64 * ii0 + 32 * ii2 + 32 * ii4 + 35 >= 2 * i0 + i2) {
                      for (int i3 = 32 * ii2 + 1; i3 < _PB_N - 1; i3 += 1) {
                        for (int i4 = max(max(1, 64 * ii0 + 32 * ii4 - 2 * i0 + 4), 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 5); i4 <= 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 36; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 37; i4 <= min(32 * ii4 - 1, 64 * ii0 + 32 * ii4 - 2 * i0 + 35); i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = max(32 * ii4, 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 37); i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 + 35; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    }
                  }
                } else if (i0 >= 32 * ii0 + 2 && _PB_N + 2 * i0 >= 64 * ii0 + 32 * ii2 + 9) {
                  if (i0 == 32 * ii0 + 2) {
                    for (int i2 = 32 * ii2; i2 <= 32 * ii2 + 2; i2 += 1) {
                      for (int i3 = max(32 * ii2, 64 * ii2 - i2 + 1); i3 < _PB_N - 1; i3 += 1) {
                        for (int i4 = max(max(max(1, 32 * ii4), 32 * ii2 + 32 * ii4 - i2 + 1), 64 * ii2 + 32 * ii4 - 2 * i3 + 1); i4 <= 32 * ii2 + 32 * ii4 - i2 + 32; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        if (i3 == 32 * ii2) {
                          for (int i4 = 32 * ii2 + 32 * ii4 - i2 + 33; i4 <= 32 * ii4 + 32; i4 += 1) {
                            B[i2][32 * ii2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii2][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2 - 1][32 * ii2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii2 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2][32 * ii2 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2][32 * ii2][i4 - 1]))) + A[i2][32 * ii2][i4]);
                          }
                        } else if (i2 == 32 * ii2 + 2) {
                          B[32 * ii2 + 2][i3][32 * ii4 + 31] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 3][i3][32 * ii4 + 31] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][i3][32 * ii4 + 31])) + A[32 * ii2 + 1][i3][32 * ii4 + 31])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3 + 1][32 * ii4 + 31] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][i3][32 * ii4 + 31])) + A[32 * ii2 + 2][i3 - 1][32 * ii4 + 31]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3][32 * ii4 + 32] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][i3][32 * ii4 + 31])) + A[32 * ii2 + 2][i3][32 * ii4 + 30]))) + A[32 * ii2 + 2][i3][32 * ii4 + 31]);
                        }
                      }
                    }
                  }
                  for (int i2 = 64 * ii0 + 32 * ii2 - 2 * i0 + 4; i2 < 32 * ii2; i2 += 1) {
                    for (int i3 = 64 * ii0 + 64 * ii2 - 2 * i0 - i2 + 5; i3 < _PB_N - 1; i3 += 1) {
                      if (i3 >= 32 * ii2 + 1) {
                        for (int i4 = max(1, 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 5); i4 <= 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 36; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      } else {
                        for (int i4 = max(1, 64 * ii0 + 64 * ii2 + 32 * ii4 - 2 * i0 - i2 - i3 + 6); i4 <= 64 * ii0 + 64 * ii2 + 32 * ii4 - 2 * i0 - i2 - i3 + 37; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    }
                  }
                  if (i0 == 32 * ii0 + 3) {
                    for (int i3 = 32 * ii2 - 1; i3 < _PB_N - 1; i3 += 1) {
                      if (i3 >= 32 * ii2 + 1) {
                        for (int i4 = max(1, 32 * ii4 - 1); i4 <= 32 * ii4 + 30; i4 += 1) {
                          B[32 * ii2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3][i4 - 1]))) + A[32 * ii2][i3][i4]);
                        }
                      } else {
                        for (int i4 = max(1, 32 * ii2 + 32 * ii4 - i3); i4 <= 32 * ii2 + 32 * ii4 - i3 + 31; i4 += 1) {
                          B[32 * ii2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2][i3][i4])) + A[32 * ii2][i3][i4 - 1]))) + A[32 * ii2][i3][i4]);
                        }
                      }
                    }
                  }
                  for (int i2 = max(32 * ii2, 64 * ii0 + 32 * ii2 - 2 * i0 + 7); i2 < _PB_N - 1; i2 += 1) {
                    for (int i3 = max(64 * ii0 + 32 * ii2 - 2 * i0 + 4, 64 * ii0 + 64 * ii2 - 2 * i0 - i2 + 5); i3 <= min(32 * ii2 - 1, 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 + 5); i3 += 1) {
                      for (int i4 = max(max(1, 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i3 + 5), 64 * ii0 + 64 * ii2 + 32 * ii4 - 2 * i0 - i2 - i3 + 6); i4 <= min(32 * ii4 + 1, 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i3 + 6); i4 += 1) {
                        B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                      }
                      if (2 * i0 + i3 >= 64 * ii0 + 32 * ii2 + 5) {
                        for (int i4 = 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i3 + 7; i4 <= 64 * ii0 + 64 * ii2 + 32 * ii4 - 2 * i0 - i2 - i3 + 37; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = max(64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i3 + 7, 64 * ii0 + 64 * ii2 + 32 * ii4 - 2 * i0 - i2 - i3 + 38); i4 <= 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i3 + 36; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      } else {
                        for (int i4 = 32 * ii4 + 2; i4 <= 32 * ii4 + 32; i4 += 1) {
                          B[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 4][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii0 + 32 * ii2 - 2 * i0 + 4][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 4][i4])) + A[i2 - 1][64 * ii0 + 32 * ii2 - 2 * i0 + 4][i4])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 5][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 4][i4])) + A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 3][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 4][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 4][i4])) + A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 4][i4 - 1]))) + A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 4][i4]);
                        }
                      }
                    }
                    if (ii4 == 0) {
                      for (int i3 = 64 * ii0 + 32 * ii2 - 2 * i0 + 6; i3 <= min(32 * ii2 - 1, 64 * ii0 + 64 * ii2 - 2 * i0 - i2 + 36); i3 += 1) {
                        for (int i4 = 1; i4 <= 64 * ii0 + 64 * ii2 - 2 * i0 - i2 - i3 + 37; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = 64 * ii0 + 64 * ii2 - 2 * i0 - i2 - i3 + 38; i4 <= 64 * ii0 + 32 * ii2 - 2 * i0 - i3 + 36; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                      for (int i3 = max(64 * ii0 + 32 * ii2 - 2 * i0 + 6, 64 * ii0 + 64 * ii2 - 2 * i0 - i2 + 37); i3 < 32 * ii2; i3 += 1) {
                        for (int i4 = 1; i4 <= 64 * ii0 + 32 * ii2 - 2 * i0 - i3 + 36; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                      if (i2 == 32 * ii2) {
                        for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 + 37; i4 += 1) {
                          B[32 * ii2][32 * ii2][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][32 * ii2][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][32 * ii2][i4])) + A[32 * ii2 - 1][32 * ii2][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2][32 * ii2 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][32 * ii2][i4])) + A[32 * ii2][32 * ii2 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2][32 * ii2][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2][32 * ii2][i4])) + A[32 * ii2][32 * ii2][i4 - 1]))) + A[32 * ii2][32 * ii2][i4]);
                        }
                      }
                    } else if (i2 == 32 * ii2) {
                      for (int i4 = 64 * ii0 + 32 * ii4 - 2 * i0 + 6; i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 + 37; i4 += 1) {
                        B[32 * ii2][32 * ii2][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][32 * ii2][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][32 * ii2][i4])) + A[32 * ii2 - 1][32 * ii2][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2][32 * ii2 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2][32 * ii2][i4])) + A[32 * ii2][32 * ii2 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2][32 * ii2][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2][32 * ii2][i4])) + A[32 * ii2][32 * ii2][i4 - 1]))) + A[32 * ii2][32 * ii2][i4]);
                      }
                    }
                    if (i0 == 32 * ii0 + 2) {
                      for (int i4 = 32 * ii4 + 1; i4 <= 32 * ii4 + 32; i4 += 1) {
                        B[i2][32 * ii2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii2][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2 - 1][32 * ii2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii2 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2][32 * ii2 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2][32 * ii2][i4 - 1]))) + A[i2][32 * ii2][i4]);
                      }
                    } else if (i2 >= 32 * ii2 + 1 && 32 * ii2 + 30 >= i2 && 64 * ii0 + 32 * ii2 + 32 * ii4 + 35 >= 2 * i0 + i2) {
                      for (int i4 = max(1, 64 * ii0 + 32 * ii4 - 2 * i0 + 5); i4 <= min(64 * ii0 + 32 * ii4 - 2 * i0 + 35, 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 37); i4 += 1) {
                        B[i2][32 * ii2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii2][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2 - 1][32 * ii2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii2 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2][32 * ii2 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2][32 * ii2][i4 - 1]))) + A[i2][32 * ii2][i4]);
                      }
                      for (int i4 = 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 38; i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 + 35; i4 += 1) {
                        B[i2][32 * ii2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii2][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2 - 1][32 * ii2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii2 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2][32 * ii2 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2][32 * ii2][i4 - 1]))) + A[i2][32 * ii2][i4]);
                      }
                      B[i2][32 * ii2][64 * ii0 + 32 * ii4 - 2 * i0 + 36] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii2][64 * ii0 + 32 * ii4 - 2 * i0 + 36] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][64 * ii0 + 32 * ii4 - 2 * i0 + 36])) + A[i2 - 1][32 * ii2][64 * ii0 + 32 * ii4 - 2 * i0 + 36])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii2 + 1][64 * ii0 + 32 * ii4 - 2 * i0 + 36] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][64 * ii0 + 32 * ii4 - 2 * i0 + 36])) + A[i2][32 * ii2 - 1][64 * ii0 + 32 * ii4 - 2 * i0 + 36]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii2][64 * ii0 + 32 * ii4 - 2 * i0 + 37] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][64 * ii0 + 32 * ii4 - 2 * i0 + 36])) + A[i2][32 * ii2][64 * ii0 + 32 * ii4 - 2 * i0 + 35]))) + A[i2][32 * ii2][64 * ii0 + 32 * ii4 - 2 * i0 + 36]);
                    } else if (ii4 == 0 && 2 * i0 + i2 >= 64 * ii0 + 32 * ii2 + 36 && 32 * ii2 + 30 >= i2) {
                      for (int i3 = 32 * ii2; i3 < -64 * ii0 + 32 * ii2 + 2 * i0 - 4; i3 += 1) {
                        if (i3 >= 32 * ii2 + 1) {
                          for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 + 35; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        } else {
                          if (2 * i0 + i2 == 64 * ii0 + 32 * ii2 + 36) {
                            B[64 * ii0 + 32 * ii2 - 2 * i0 + 36][32 * ii2][1] = ((((SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 37][32 * ii2][1] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 36][32 * ii2][1])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 35][32 * ii2][1])) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 36][32 * ii2 + 1][1] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 36][32 * ii2][1])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 36][32 * ii2 - 1][1]))) + (SCALAR_VAL(0.125) * ((A[64 * ii0 + 32 * ii2 - 2 * i0 + 36][32 * ii2][2] - (SCALAR_VAL(2.0) * A[64 * ii0 + 32 * ii2 - 2 * i0 + 36][32 * ii2][1])) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 36][32 * ii2][0]))) + A[64 * ii0 + 32 * ii2 - 2 * i0 + 36][32 * ii2][1]);
                          }
                          for (int i4 = max(1, 64 * ii0 + 32 * ii2 - 2 * i0 - i2 + 38); i4 <= 64 * ii0 - 2 * i0 + 36; i4 += 1) {
                            B[i2][32 * ii2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii2][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2 - 1][32 * ii2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii2 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2][32 * ii2 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2][32 * ii2][i4 - 1]))) + A[i2][32 * ii2][i4]);
                          }
                        }
                      }
                    }
                    if (i2 >= 32 * ii2 + 31) {
                      for (int i3 = 32 * ii2; i3 < -64 * ii0 + 32 * ii2 + 2 * i0 - 4; i3 += 1) {
                        if (i3 >= 32 * ii2 + 1) {
                          for (int i4 = max(1, 64 * ii0 + 32 * ii4 - 2 * i0 + 4); i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 + 35; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        } else {
                          for (int i4 = max(1, 64 * ii0 + 32 * ii4 - 2 * i0 + 5); i4 <= 32 * ii4; i4 += 1) {
                            B[i2][32 * ii2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii2][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2 - 1][32 * ii2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii2 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2][32 * ii2 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2][32 * ii2][i4 - 1]))) + A[i2][32 * ii2][i4]);
                          }
                          for (int i4 = 32 * ii4 + 1; i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 + 35; i4 += 1) {
                            B[i2][32 * ii2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii2][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2 - 1][32 * ii2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii2 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2][32 * ii2 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2][32 * ii2][i4 - 1]))) + A[i2][32 * ii2][i4]);
                          }
                          B[i2][32 * ii2][64 * ii0 + 32 * ii4 - 2 * i0 + 36] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii2][64 * ii0 + 32 * ii4 - 2 * i0 + 36] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][64 * ii0 + 32 * ii4 - 2 * i0 + 36])) + A[i2 - 1][32 * ii2][64 * ii0 + 32 * ii4 - 2 * i0 + 36])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii2 + 1][64 * ii0 + 32 * ii4 - 2 * i0 + 36] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][64 * ii0 + 32 * ii4 - 2 * i0 + 36])) + A[i2][32 * ii2 - 1][64 * ii0 + 32 * ii4 - 2 * i0 + 36]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii2][64 * ii0 + 32 * ii4 - 2 * i0 + 37] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][64 * ii0 + 32 * ii4 - 2 * i0 + 36])) + A[i2][32 * ii2][64 * ii0 + 32 * ii4 - 2 * i0 + 35]))) + A[i2][32 * ii2][64 * ii0 + 32 * ii4 - 2 * i0 + 36]);
                        }
                      }
                    } else if (64 * ii0 + 32 * ii2 + 32 * ii4 + 35 >= 2 * i0 + i2) {
                      for (int i3 = 32 * ii2 + 1; i3 < -64 * ii0 + 32 * ii2 + 2 * i0 - 4; i3 += 1) {
                        for (int i4 = max(max(1, 64 * ii0 + 32 * ii4 - 2 * i0 + 4), 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 5); i4 <= 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 36; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 37; i4 < 32 * ii4; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = max(32 * ii4, 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 37); i4 <= min(64 * ii0 + 32 * ii4 - 2 * i0 + 35, 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i3 + 37); i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = max(64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 37, 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i3 + 38); i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 + 35; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    }
                    for (int i3 = max(32 * ii2 + 1, -64 * ii0 + 32 * ii2 + 2 * i0 - 4); i3 <= min(min(min(_PB_N - 2, 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 + 36), 64 * ii0 + 32 * ii2 - 2 * i0 + 37), 64 * ii0 - 2 * i0 + i2 + 36); i3 += 1) {
                      if (i0 >= 32 * ii0 + 3) {
                        for (int i4 = max(max(1, 64 * ii0 + 32 * ii4 - 2 * i0 + 4), 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 5); i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 + 5; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                      for (int i4 = max(1, 64 * ii0 + 32 * ii4 - 2 * i0 + 6); i4 <= 32 * ii2 + 32 * ii4 - i2; i4 += 1) {
                        B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                      }
                      for (int i4 = max(max(1, 64 * ii0 + 32 * ii4 - 2 * i0 + 6), 32 * ii2 + 32 * ii4 - i2 + 1); i4 <= min(32 * ii4 - 1, 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 36); i4 += 1) {
                        B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                      }
                      for (int i4 = max(max(1, 64 * ii0 + 32 * ii4 - 2 * i0 + 6), 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 37); i4 < 32 * ii4; i4 += 1) {
                        B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                      }
                      if (i0 == 32 * ii0 + 2) {
                        for (int i4 = max(1, 32 * ii4); i4 <= min(32 * ii2 + 32 * ii4 - i2 + 32, 64 * ii2 + 32 * ii4 - 2 * i3 + 32); i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      } else {
                        for (int i4 = max(max(1, 32 * ii4), 32 * ii2 + 32 * ii4 - i2 + 1); i4 <= min(64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 36, 64 * ii2 + 32 * ii4 - 2 * i3 + 32); i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                      for (int i4 = max(max(1, 32 * ii4), 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 37); i4 <= 64 * ii2 + 32 * ii4 - 2 * i3 + 32; i4 += 1) {
                        B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                      }
                      for (int i4 = max(max(max(1, 32 * ii4), 32 * ii2 + 32 * ii4 - i2 + 1), 64 * ii2 + 32 * ii4 - 2 * i3 + 33); i4 <= min(32 * ii4 + 31, 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i3 + 37); i4 += 1) {
                        B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                      }
                      if (i0 == 32 * ii0 + 2) {
                        for (int i4 = 32 * ii2 + 32 * ii4 - i3 + 34; i4 <= 32 * ii2 + 32 * ii4 - i2 + 32; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      } else {
                        for (int i4 = 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i3 + 38; i4 <= 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 36; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                      for (int i4 = max(64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 37, 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i3 + 38); i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 + 35; i4 += 1) {
                        B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                      }
                    }
                    if (ii4 == 0 && _PB_N + 2 * i0 >= 64 * ii0 + 32 * ii2 + 39 && 2 * i0 + i2 >= 64 * ii0 + 32 * ii2 + 36) {
                      for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 + 35; i4 += 1) {
                        B[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 37][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii0 + 32 * ii2 - 2 * i0 + 37][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 37][i4])) + A[i2 - 1][64 * ii0 + 32 * ii2 - 2 * i0 + 37][i4])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 38][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 37][i4])) + A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 36][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 37][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 37][i4])) + A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 37][i4 - 1]))) + A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 37][i4]);
                      }
                    }
                    if (ii4 == 0 && 2 * i0 + i2 >= 64 * ii0 + 32 * ii2 + 36 && 32 * ii2 + 30 >= i2) {
                      for (int i3 = 64 * ii0 + 32 * ii2 - 2 * i0 + 38; i3 < _PB_N - 1; i3 += 1) {
                        for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 + 35; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    } else if (i2 >= 32 * ii2 + 31) {
                      for (int i3 = 64 * ii0 + 32 * ii2 - 2 * i0 + 38; i3 < _PB_N - 1; i3 += 1) {
                        for (int i4 = max(1, 64 * ii0 + 32 * ii4 - 2 * i0 + 4); i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 + 35; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    } else {
                      if (ii4 == 0 && _PB_N + 2 * i0 >= 64 * ii0 + 32 * ii2 + 39 && i2 >= 32 * ii2 + 1 && 64 * ii0 + 32 * ii2 + 35 >= 2 * i0 + i2) {
                        for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 + 35; i4 += 1) {
                          B[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 37][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii0 + 32 * ii2 - 2 * i0 + 37][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 37][i4])) + A[i2 - 1][64 * ii0 + 32 * ii2 - 2 * i0 + 37][i4])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 38][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 37][i4])) + A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 36][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 37][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 37][i4])) + A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 37][i4 - 1]))) + A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 37][i4]);
                        }
                      }
                      for (int i3 = 64 * ii0 + 32 * ii2 - 2 * i0 + 38; i3 <= min(_PB_N - 2, 64 * ii0 - 2 * i0 + i2 + 36); i3 += 1) {
                        for (int i4 = max(1, 64 * ii0 + 32 * ii4 - 2 * i0 + 4); i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 + 35; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    }
                    for (int i3 = 64 * ii0 - 2 * i0 + i2 + 37; i3 < _PB_N - 1; i3 += 1) {
                      for (int i4 = max(max(1, 64 * ii0 + 32 * ii4 - 2 * i0 + 4), 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 5); i4 <= 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 36; i4 += 1) {
                        B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                      }
                      for (int i4 = 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 37; i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 + 35; i4 += 1) {
                        B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                      }
                    }
                  }
                } else if (32 * ii2 + 4 >= _PB_N && i0 == 32 * ii0 + 2) {
                  for (int i2 = 32 * ii2; i2 < _PB_N - 1; i2 += 1) {
                    for (int i3 = max(32 * ii2, 64 * ii2 - i2 + 1); i3 < _PB_N - 1; i3 += 1) {
                      for (int i4 = max(max(max(1, 32 * ii4), 32 * ii2 + 32 * ii4 - i2 + 1), 64 * ii2 + 32 * ii4 - 2 * i3 + 1); i4 <= 32 * ii2 + 32 * ii4 - i2 + 32; i4 += 1) {
                        B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                      }
                      if (32 * ii2 + 4 == _PB_N && i2 + 2 == _PB_N && i3 + 3 >= _PB_N) {
                        B[_PB_N - 2][i3][32 * ii4 + 31] = ((((SCALAR_VAL(0.125) * ((A[_PB_N - 1][i3][32 * ii4 + 31] - (SCALAR_VAL(2.0) * A[_PB_N - 2][i3][32 * ii4 + 31])) + A[_PB_N - 3][i3][32 * ii4 + 31])) + (SCALAR_VAL(0.125) * ((A[_PB_N - 2][i3 + 1][32 * ii4 + 31] - (SCALAR_VAL(2.0) * A[_PB_N - 2][i3][32 * ii4 + 31])) + A[_PB_N - 2][i3 - 1][32 * ii4 + 31]))) + (SCALAR_VAL(0.125) * ((A[_PB_N - 2][i3][32 * ii4 + 32] - (SCALAR_VAL(2.0) * A[_PB_N - 2][i3][32 * ii4 + 31])) + A[_PB_N - 2][i3][32 * ii4 + 30]))) + A[_PB_N - 2][i3][32 * ii4 + 31]);
                      } else if (i3 == 32 * ii2) {
                        for (int i4 = 32 * ii2 + 32 * ii4 - i2 + 33; i4 <= 32 * ii4 + 32; i4 += 1) {
                          B[i2][32 * ii2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii2][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2 - 1][32 * ii2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii2 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2][32 * ii2 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2][32 * ii2][i4 - 1]))) + A[i2][32 * ii2][i4]);
                        }
                      }
                    }
                  }
                }
                for (int i2 = 64 * ii0 + 32 * ii2 - 2 * i0 + 3; i2 <= 32 * ii2; i2 += 1) {
                  for (int i3 = 64 * ii0 + 64 * ii2 - 2 * i0 - i2 + 4; i3 < _PB_N - 1; i3 += 1) {
                    for (int i4 = max(max(1, 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 4), 64 * ii0 + 64 * ii2 + 32 * ii4 - 2 * i0 - i2 - i3 + 5); i4 <= 64 * ii0 + 64 * ii2 + 32 * ii4 - 2 * i0 - i2 - i3 + 36; i4 += 1) {
                      A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                    }
                    for (int i4 = max(1, 64 * ii0 + 64 * ii2 + 32 * ii4 - 2 * i0 - i2 - i3 + 37); i4 <= 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i2 + 35; i4 += 1) {
                      A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                    }
                  }
                }
                for (int i2 = 32 * ii2 + 1; i2 < _PB_N - 1; i2 += 1) {
                  for (int i3 = 64 * ii0 + 32 * ii2 - 2 * i0 + 3; i3 <= 32 * ii2; i3 += 1) {
                    for (int i4 = max(1, 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i3 + 4); i4 <= 64 * ii0 + 32 * ii2 + 32 * ii4 - 2 * i0 - i3 + 35; i4 += 1) {
                      A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                    }
                  }
                  for (int i3 = 32 * ii2 + 1; i3 < _PB_N - 1; i3 += 1) {
                    for (int i4 = max(1, 64 * ii0 + 32 * ii4 - 2 * i0 + 3); i4 <= 64 * ii0 + 32 * ii4 - 2 * i0 + 34; i4 += 1) {
                      A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                    }
                  }
                }
              }
              if (ii4 == 0) {
                for (int i0 = 32 * ii0 + 17; i0 <= min(_PB_TSTEPS, 32 * ii0 + 32); i0 += 1) {
                  for (int i2 = max(1, 64 * ii0 + 32 * ii2 - 2 * i0 + 4); i2 < _PB_N - 1; i2 += 1) {
                    if (32 * ii0 + 16 * ii2 + 1 >= i0 && i2 >= 32 * ii2 + 1 && 32 * ii2 + 30 >= i2) {
                      for (int i4 = 1; i4 <= 32; i4 += 1) {
                        B[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 4][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii0 + 32 * ii2 - 2 * i0 + 4][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 4][i4])) + A[i2 - 1][64 * ii0 + 32 * ii2 - 2 * i0 + 4][i4])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 5][i4] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 4][i4])) + A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 3][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 4][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 4][i4])) + A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 4][i4 - 1]))) + A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 4][i4]);
                      }
                    }
                    if (32 * ii2 + 30 >= i2) {
                      for (int i3 = max(max(1, 64 * ii0 + 32 * ii2 - 2 * i0 + 5), 64 * ii0 + 64 * ii2 - 2 * i0 - i2 + 5); i3 <= min(32 * ii2 - 1, 64 * ii0 + 64 * ii2 - 2 * i0 - i2 + 36); i3 += 1) {
                        if (2 * i0 + i3 == 64 * ii0 + 32 * ii2 + 5) {
                          B[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 5][1] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][64 * ii0 + 32 * ii2 - 2 * i0 + 5][1] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 5][1])) + A[i2 - 1][64 * ii0 + 32 * ii2 - 2 * i0 + 5][1])) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 6][1] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 5][1])) + A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 4][1]))) + (SCALAR_VAL(0.125) * ((A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 5][2] - (SCALAR_VAL(2.0) * A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 5][1])) + A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 5][0]))) + A[i2][64 * ii0 + 32 * ii2 - 2 * i0 + 5][1]);
                        }
                        for (int i4 = max(1, 64 * ii0 + 32 * ii2 - 2 * i0 - i3 + 7); i4 <= 64 * ii0 + 64 * ii2 - 2 * i0 - i2 - i3 + 37; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = 64 * ii0 + 64 * ii2 - 2 * i0 - i2 - i3 + 38; i4 <= 64 * ii0 + 32 * ii2 - 2 * i0 - i3 + 36; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                      for (int i3 = max(1, 64 * ii0 + 64 * ii2 - 2 * i0 - i2 + 37); i3 <= min(32 * ii2 - 1, 64 * ii0 + 32 * ii2 - 2 * i0 + 35); i3 += 1) {
                        for (int i4 = 1; i4 <= 64 * ii0 + 32 * ii2 - 2 * i0 - i3 + 36; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                      if (2 * i0 + i2 >= 64 * ii0 + 32 * ii2 + 5 && 32 * ii2 >= i2) {
                        for (int i4 = 1; i4 <= 64 * ii0 + 32 * ii2 - 2 * i0 - i2 + 37; i4 += 1) {
                          B[i2][32 * ii2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii2][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2 - 1][32 * ii2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii2 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2][32 * ii2 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2][32 * ii2][i4 - 1]))) + A[i2][32 * ii2][i4]);
                        }
                      }
                    } else {
                      for (int i3 = max(1, 64 * ii0 + 32 * ii2 - 2 * i0 + 4); i3 <= 64 * ii0 + 32 * ii2 - 2 * i0 + 5; i3 += 1) {
                        for (int i4 = 1; i4 <= 64 * ii0 + 32 * ii2 - 2 * i0 - i3 + 36; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                      for (int i3 = max(1, 64 * ii0 + 32 * ii2 - 2 * i0 + 6); i3 <= min(32 * ii2 - 1, 64 * ii0 + 32 * ii2 - 2 * i0 + 35); i3 += 1) {
                        for (int i4 = 1; i4 <= 64 * ii0 + 32 * ii2 - 2 * i0 - i3 + 36; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    }
                    if (i0 == 32 * ii0 + 17 && i2 >= 32 * ii2 + 1) {
                      for (int i4 = 1; i4 <= 2; i4 += 1) {
                        B[i2][32 * ii2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii2][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2 - 1][32 * ii2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii2 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2][32 * ii2 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2][32 * ii2][i4 - 1]))) + A[i2][32 * ii2][i4]);
                      }
                      if (i2 >= 32 * ii2 + 2) {
                        for (int i3 = 32 * ii2 + 1; i3 < _PB_N - 1; i3 += 1) {
                          B[i2][i3][1] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][1] - (SCALAR_VAL(2.0) * A[i2][i3][1])) + A[i2 - 1][i3][1])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][1] - (SCALAR_VAL(2.0) * A[i2][i3][1])) + A[i2][i3 - 1][1]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][2] - (SCALAR_VAL(2.0) * A[i2][i3][1])) + A[i2][i3][0]))) + A[i2][i3][1]);
                        }
                      }
                    }
                    for (int i3 = 32 * ii2 + 1; i3 < _PB_N - 1; i3 += 1) {
                      for (int i4 = 1; i4 <= 64 * ii0 + 32 * ii2 - 2 * i0 - i2 + 36; i4 += 1) {
                        B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                      }
                    }
                  }
                  for (int i2 = max(1, 64 * ii0 + 32 * ii2 - 2 * i0 + 3); i2 <= 32 * ii2; i2 += 1) {
                    for (int i3 = max(1, 64 * ii0 + 64 * ii2 - 2 * i0 - i2 + 4); i3 <= min(32 * ii2 + 1, 64 * ii0 + 64 * ii2 - 2 * i0 - i2 + 35); i3 += 1) {
                      for (int i4 = 1; i4 <= 64 * ii0 + 64 * ii2 - 2 * i0 - i2 - i3 + 36; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                    for (int i3 = 32 * ii2 + 2; i3 < _PB_N - 1; i3 += 1) {
                      for (int i4 = 1; i4 <= 64 * ii0 + 32 * ii2 - 2 * i0 - i2 + 35; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                  }
                  for (int i2 = 32 * ii2 + 1; i2 < _PB_N - 1; i2 += 1) {
                    for (int i3 = max(1, 64 * ii0 + 32 * ii2 - 2 * i0 + 3); i3 <= 64 * ii0 + 32 * ii2 - 2 * i0 + 34; i3 += 1) {
                      for (int i4 = 1; i4 <= 64 * ii0 + 32 * ii2 - 2 * i0 - i3 + 35; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                  }
                }
              }
            }
            for (int i0 = 32 * ii0 + 1; i0 <= min(_PB_TSTEPS, 32 * ii0 + 32); i0 += 1) {
              if (i0 >= 32 * ii0 + 2) {
                for (int i2 = max(1, 64 * ii0 + 32 * ii2 - 2 * i0 + 4); i2 < _PB_N - 1; i2 += 1) {
                  if (i0 >= 32 * ii0 + 3) {
                    if (32 * ii2 + 3 == _PB_N && i2 + 2 == _PB_N) {
                      for (int i3 = max(-64 * ii0 + 2 * i0 - 1, _PB_N + 64 * ii0 - 2 * i0 + 1); i3 <= min(_PB_N - 4, _PB_N + 64 * ii0 - 2 * i0 + 3); i3 += 1) {
                        for (int i4 = 2 * _PB_N + 64 * ii0 - 2 * i0 - i3 - 1; i4 < _PB_N - 1; i4 += 1) {
                          B[_PB_N - 2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[_PB_N - 1][i3][i4] - (SCALAR_VAL(2.0) * A[_PB_N - 2][i3][i4])) + A[_PB_N - 3][i3][i4])) + (SCALAR_VAL(0.125) * ((A[_PB_N - 2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[_PB_N - 2][i3][i4])) + A[_PB_N - 2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[_PB_N - 2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[_PB_N - 2][i3][i4])) + A[_PB_N - 2][i3][i4 - 1]))) + A[_PB_N - 2][i3][i4]);
                        }
                      }
                      for (int i3 = _PB_N + 64 * ii0 - 2 * i0 + 1; i3 <= min(-64 * ii0 + 2 * i0 - 2, _PB_N + 64 * ii0 - 2 * i0 + 3); i3 += 1) {
                        for (int i4 = 2 * _PB_N + 64 * ii0 - 2 * i0 - i3 - 1; i4 < _PB_N - 1; i4 += 1) {
                          B[_PB_N - 2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[_PB_N - 1][i3][i4] - (SCALAR_VAL(2.0) * A[_PB_N - 2][i3][i4])) + A[_PB_N - 3][i3][i4])) + (SCALAR_VAL(0.125) * ((A[_PB_N - 2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[_PB_N - 2][i3][i4])) + A[_PB_N - 2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[_PB_N - 2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[_PB_N - 2][i3][i4])) + A[_PB_N - 2][i3][i4 - 1]))) + A[_PB_N - 2][i3][i4]);
                        }
                      }
                    } else {
                      for (int i3 = max(max(1, 64 * ii0 + 32 * ii2 - 2 * i0 + 4), 64 * ii0 + 64 * ii2 - 2 * i0 - i2 + 5); i3 <= min(min(min(min(_PB_N - 2, -64 * ii0 + 2 * i0 + i2 - 3), 2 * _PB_N + 64 * ii0 + 64 * ii2 - 2 * i0 - 3 * i2 + 3), _PB_N - 64 * ii0 - 53 * ii2 + 2 * i0 + 2 * i2 - (ii2 + i2 + 1) / 3 - 7), -64 * ii0 + 11 * ii2 + 2 * i0 - (ii2 + i2 + 1) / 3 - 2); i3 += 1) {
                        for (int i4 = max(max(1, 64 * ii0 + 64 * ii2 - 2 * i0 - i2 + 5), 64 * ii0 + 96 * ii2 - 2 * i0 - i2 - i3 + 6); i4 <= 64 * ii0 + 53 * ii2 - 2 * i0 - 2 * i2 + i3 + (ii2 + i2 + 1) / 3 + 4; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = max(1, 64 * ii0 + 64 * ii2 - 2 * i0 - i3 + 5); i4 <= 64 * ii0 + 60 * ii2 - 2 * i0 - i3 + (4 * ii2 + i2 - 2) / 9 + 5; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = max(max(max(1, 64 * ii0 + 96 * ii2 - 2 * i0 - i2 - i3 + 6), 64 * ii0 + 53 * ii2 - 2 * i0 - 2 * i2 + i3 + (ii2 + i2 + 1) / 3 + 5), 64 * ii0 + 60 * ii2 - 2 * i0 - i3 + (4 * ii2 + i2 - 2) / 9 + 6); i4 <= min(64 * ii0 + 96 * ii2 - 2 * i0 - i2 - i3 + 8, 11 * ii2 + i2 - (ii2 + i2 + 1) / 3); i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = max(max(1, 64 * ii0 + 96 * ii2 - 2 * i0 - i2 - i3 + 9), 64 * ii0 + 60 * ii2 - 2 * i0 - i3 + (4 * ii2 + i2 - 2) / 9 + 6); i4 <= min(64 * ii0 + 64 * ii2 - 2 * i0 - i3 + 7, 11 * ii2 + i2 - (ii2 + i2 + 1) / 3); i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        if (ii2 == 1) {
                          for (int i4 = max(max(max(1, 64 * ii0 - 2 * i0 - i3 + 72), 64 * ii0 - 2 * i0 - i2 - i3 + 105), 64 * ii0 - 2 * i0 - 2 * i2 + i3 + (i2 - 1) / 3 + 59); i4 <= min(-i2 + 64, -32 * ii0 + i0 + i2 + (i2 + i3) / 2 - 36); i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (i2 <= 31 && i2 + i3 >= 65) {
                            for (int i4 = -i2 + 65; i4 < -32 * ii0 + i0 + i2 + (i2 + i3) / 2 - 35; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                          if (i2 + i3 >= 65) {
                            for (int i4 = max(-i2 + 65, 64 * ii0 - 2 * i0 - i3 + 72); i4 <= min(-2 * i3 + 96, -32 * ii0 + i0 + i2 + (i2 + i3) / 2 - 36); i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            if (i3 <= 31) {
                              for (int i4 = -2 * i3 + 97; i4 < -32 * ii0 + i0 + i2 + (i2 + i3) / 2 - 35; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                          }
                          for (int i4 = max(-i2 + 65, -2 * i3 + 97); i4 <= 31; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (i2 + i3 <= 64) {
                            for (int i4 = max(-i2 + 65, 64 * ii0 - 2 * i0 - i3 + 72); i4 < -32 * ii0 + i0 + i2 + (i2 + i3) / 2 - 35; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                          if (i2 >= 32 && i2 + i3 >= 65 && i3 >= 32) {
                            for (int i4 = max(max(32, -i2 + 65), -2 * i3 + 97); i4 < -32 * ii0 + i0 + i2 + (i2 + i3) / 2 - 35; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                        } else {
                          for (int i4 = max(max(64 * ii0 + 64 * ii2 - 2 * i0 - i3 + 8, 64 * ii0 + 96 * ii2 - 2 * i0 - i2 - i3 + 9), 64 * ii0 + 53 * ii2 - 2 * i0 - 2 * i2 + i3 + (ii2 + i2 + 1) / 3 + 5); i4 < -32 * ii0 - 32 * ii2 + i0 + i2 + (i2 + i3) / 2 - 3; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                        for (int i4 = max(max(max(64 * ii0 + 96 * ii2 - 2 * i0 - i2 - i3 + 6, 11 * ii2 + i2 - (ii2 + i2 + 1) / 3 + 1), 64 * ii0 + 53 * ii2 - 2 * i0 - 2 * i2 + i3 + (ii2 + i2 + 1) / 3 + 5), -32 * ii0 - 32 * ii2 + i0 + i2 + (i2 + i3) / 2 - 3); i4 < _PB_N - 1; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    }
                    if (32 * ii2 + 1 >= i2) {
                      for (int i3 = max(1, 2 * _PB_N + 64 * ii0 + 64 * ii2 - 2 * i0 - 3 * i2 + 4); i3 <= min(min(_PB_N - 2, 64 * ii2 - i2), -64 * ii0 + 11 * ii2 + 2 * i0 - (ii2 + i2 + 1) / 3 - 2); i3 += 1) {
                        for (int i4 = max(max(1, 64 * ii0 + 64 * ii2 - 2 * i0 - i2 + 5), 64 * ii0 + 96 * ii2 - 2 * i0 - i2 - i3 + 6); i4 < _PB_N - 1; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                      if (3 * _PB_N >= 32 * ii2 + 2 * i2 + 8) {
                        for (int i3 = max(64 * ii0 + 64 * ii2 - 2 * i0 - i2 + 5, -64 * ii0 + 11 * ii2 + 2 * i0 - (ii2 + i2 + 1) / 3 - 1); i3 <= min(min(min(_PB_N - 2, -64 * ii0 + 2 * i0 + i2 - 3), 2 * _PB_N + 64 * ii0 + 64 * ii2 - 2 * i0 - 3 * i2 + 3), _PB_N - 64 * ii0 - 53 * ii2 + 2 * i0 + 2 * i2 - (ii2 + i2 + 1) / 3 - 7); i3 += 1) {
                          for (int i4 = max(max(1, 64 * ii0 + 64 * ii2 - 2 * i0 - i2 + 5), 64 * ii0 + 96 * ii2 - 2 * i0 - i2 - i3 + 6); i4 <= min(64 * ii2 - i2, 64 * ii0 + 53 * ii2 - 2 * i0 - 2 * i2 + i3 + (ii2 + i2 + 1) / 3 + 4); i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          for (int i4 = max(64 * ii0 + 96 * ii2 - 2 * i0 - i2 - i3 + 6, 64 * ii0 + 53 * ii2 - 2 * i0 - 2 * i2 + i3 + (ii2 + i2 + 1) / 3 + 5); i4 <= min(min(64 * ii2 - i2, 64 * ii0 + 96 * ii2 - 2 * i0 - i2 - i3 + 8), 11 * ii2 + i2 - (ii2 + i2 + 1) / 3); i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          for (int i4 = max(64 * ii0 + 96 * ii2 - 2 * i0 - i2 - i3 + 9, 64 * ii0 + 53 * ii2 - 2 * i0 - 2 * i2 + i3 + (ii2 + i2 + 1) / 3 + 5); i4 <= min(64 * ii2 - i2, -32 * ii0 - 32 * ii2 + i0 + i2 + (i2 + i3) / 2 - 4); i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (i2 + 1 == 32 * ii2 && 64 * ii0 + i3 + 4 == 32 * ii2 + 2 * i0) {
                            B[32 * ii2 - 1][-64 * ii0 + 32 * ii2 + 2 * i0 - 4][32 * ii2 + 2] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2][-64 * ii0 + 32 * ii2 + 2 * i0 - 4][32 * ii2 + 2] - (SCALAR_VAL(2.0) * A[32 * ii2 - 1][-64 * ii0 + 32 * ii2 + 2 * i0 - 4][32 * ii2 + 2])) + A[32 * ii2 - 2][-64 * ii0 + 32 * ii2 + 2 * i0 - 4][32 * ii2 + 2])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 - 1][-64 * ii0 + 32 * ii2 + 2 * i0 - 3][32 * ii2 + 2] - (SCALAR_VAL(2.0) * A[32 * ii2 - 1][-64 * ii0 + 32 * ii2 + 2 * i0 - 4][32 * ii2 + 2])) + A[32 * ii2 - 1][-64 * ii0 + 32 * ii2 + 2 * i0 - 5][32 * ii2 + 2]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 - 1][-64 * ii0 + 32 * ii2 + 2 * i0 - 4][32 * ii2 + 3] - (SCALAR_VAL(2.0) * A[32 * ii2 - 1][-64 * ii0 + 32 * ii2 + 2 * i0 - 4][32 * ii2 + 2])) + A[32 * ii2 - 1][-64 * ii0 + 32 * ii2 + 2 * i0 - 4][32 * ii2 + 1]))) + A[32 * ii2 - 1][-64 * ii0 + 32 * ii2 + 2 * i0 - 4][32 * ii2 + 2]);
                          }
                          if (32 * ii2 >= i2 + 1 && i2 + i3 >= 64 * ii2 + 1) {
                            for (int i4 = max(64 * ii2 - i2 + 1, 64 * ii0 + 53 * ii2 - 2 * i0 - 2 * i2 + i3 + (ii2 + i2 + 1) / 3 + 5); i4 < -32 * ii0 - 32 * ii2 + i0 + i2 + (i2 + i3) / 2 - 3; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          } else {
                            if (i0 >= 32 * ii0 + 4 && i2 == 32 * ii2 + 1 && i3 == 32 * ii2) {
                              B[32 * ii2 + 1][32 * ii2][32 * ii2] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][32 * ii2][32 * ii2] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][32 * ii2][32 * ii2])) + A[32 * ii2][32 * ii2][32 * ii2])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][32 * ii2 + 1][32 * ii2] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][32 * ii2][32 * ii2])) + A[32 * ii2 + 1][32 * ii2 - 1][32 * ii2]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][32 * ii2][32 * ii2 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][32 * ii2][32 * ii2])) + A[32 * ii2 + 1][32 * ii2][32 * ii2 - 1]))) + A[32 * ii2 + 1][32 * ii2][32 * ii2]);
                            }
                            if (i2 == 32 * ii2 + 1) {
                              for (int i4 = max(32 * ii2, 64 * ii0 + 64 * ii2 - 2 * i0 - i3 + 5); i4 <= min(min(32 * ii2 + 1, 96 * ii2 - 2 * i3), 64 * ii0 + 64 * ii2 - 2 * i0 - i3 + 7); i4 += 1) {
                                B[32 * ii2 + 1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][i3][i4])) + A[32 * ii2 + 1][i3][i4 - 1]))) + A[32 * ii2 + 1][i3][i4]);
                              }
                            }
                            if (64 * ii2 >= i2 + i3) {
                              for (int i4 = max(64 * ii2 - i2 + 1, 64 * ii0 + 96 * ii2 - 2 * i0 - i2 - i3 + 9); i4 < -32 * ii0 - 32 * ii2 + i0 + i2 + (i2 + i3) / 2 - 3; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            } else {
                              if (i0 == 32 * ii0 + 3 && i2 == 32 * ii2 && i3 == 32 * ii2 + 3) {
                                B[32 * ii2][32 * ii2 + 3][32 * ii2 + 1] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][32 * ii2 + 3][32 * ii2 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2][32 * ii2 + 3][32 * ii2 + 1])) + A[32 * ii2 - 1][32 * ii2 + 3][32 * ii2 + 1])) + (SCALAR_VAL(0.125) * ((A[32 * ii2][32 * ii2 + 4][32 * ii2 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2][32 * ii2 + 3][32 * ii2 + 1])) + A[32 * ii2][32 * ii2 + 2][32 * ii2 + 1]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2][32 * ii2 + 3][32 * ii2 + 2] - (SCALAR_VAL(2.0) * A[32 * ii2][32 * ii2 + 3][32 * ii2 + 1])) + A[32 * ii2][32 * ii2 + 3][32 * ii2]))) + A[32 * ii2][32 * ii2 + 3][32 * ii2 + 1]);
                              } else if (i0 == 32 * ii0 + 3 && i2 == 32 * ii2 + 1 && i3 == 32 * ii2) {
                                B[32 * ii2 + 1][32 * ii2][32 * ii2 + 1] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][32 * ii2][32 * ii2 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][32 * ii2][32 * ii2 + 1])) + A[32 * ii2][32 * ii2][32 * ii2 + 1])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][32 * ii2 + 1][32 * ii2 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][32 * ii2][32 * ii2 + 1])) + A[32 * ii2 + 1][32 * ii2 - 1][32 * ii2 + 1]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 1][32 * ii2][32 * ii2 + 2] - (SCALAR_VAL(2.0) * A[32 * ii2 + 1][32 * ii2][32 * ii2 + 1])) + A[32 * ii2 + 1][32 * ii2][32 * ii2]))) + A[32 * ii2 + 1][32 * ii2][32 * ii2 + 1]);
                              }
                              for (int i4 = max(64 * ii2 - i2 + 1, 96 * ii2 - 2 * i3 + 1); i4 < -32 * ii0 - 32 * ii2 + i0 + i2 + (i2 + i3) / 2 - 3; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                          }
                          for (int i4 = max(max(max(64 * ii0 + 96 * ii2 - 2 * i0 - i2 - i3 + 6, 11 * ii2 + i2 - (ii2 + i2 + 1) / 3 + 1), 64 * ii0 + 53 * ii2 - 2 * i0 - 2 * i2 + i3 + (ii2 + i2 + 1) / 3 + 5), -32 * ii0 - 32 * ii2 + i0 + i2 + (i2 + i3) / 2 - 3); i4 < _PB_N - 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                      for (int i3 = max(2 * _PB_N + 64 * ii0 + 64 * ii2 - 2 * i0 - 3 * i2 + 4, -64 * ii0 + 11 * ii2 + 2 * i0 - (ii2 + i2 + 1) / 3 - 1); i3 <= min(_PB_N - 2, 64 * ii2 - i2); i3 += 1) {
                        for (int i4 = max(64 * ii0 + 64 * ii2 - 2 * i0 - i2 + 5, 64 * ii0 + 96 * ii2 - 2 * i0 - i2 - i3 + 6); i4 < _PB_N - 1; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                      if (ii2 == 1 && i2 <= 31) {
                        for (int i3 = max(-i2 + 65, 2 * _PB_N + 64 * ii0 - 2 * i0 - 3 * i2 + 68); i3 <= min(_PB_N - 2, -64 * ii0 + 2 * i0 - (i2 + 2) / 3 + 9); i3 += 1) {
                          for (int i4 = max(1, 64 * ii0 - 2 * i0 - i2 + 69); i4 < _PB_N - 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                      if (32 * ii2 >= i2 + 1) {
                        for (int i3 = max(max(64 * ii2 - i2 + 1, 2 * _PB_N + 64 * ii0 + 64 * ii2 - 2 * i0 - 3 * i2 + 4), -64 * ii0 + 11 * ii2 + 2 * i0 - (ii2 + i2 + 1) / 3 - 1); i3 < _PB_N - 1; i3 += 1) {
                          for (int i4 = max(1, 64 * ii0 + 64 * ii2 - 2 * i0 - i2 + 5); i4 < _PB_N - 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                    }
                    for (int i3 = max(64 * ii0 + 64 * ii2 - 2 * i0 - i2 + 5, _PB_N - 64 * ii0 - 53 * ii2 + 2 * i0 + 2 * i2 - (ii2 + i2 + 1) / 3 - 6); i3 < min(_PB_N - 1, -64 * ii0 + 2 * i0 + i2 - 2); i3 += 1) {
                      for (int i4 = max(64 * ii0 + 64 * ii2 - 2 * i0 - i2 + 5, 64 * ii0 + 96 * ii2 - 2 * i0 - i2 - i3 + 6); i4 < _PB_N - 1; i4 += 1) {
                        B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                      }
                    }
                    if (i2 >= 32 * ii2 + 2) {
                      for (int i3 = max(64 * ii0 + 32 * ii2 - 2 * i0 + 4, -64 * ii0 + 11 * ii2 + 2 * i0 - (ii2 + i2 + 1) / 3 - 1); i3 <= min(2 * _PB_N + 64 * ii0 + 64 * ii2 - 2 * i0 - 3 * i2 + 3, 36 * ii2 - (4 * ii2 + i2 + 6) / 9 + 1); i3 += 1) {
                        for (int i4 = 64 * ii0 + 64 * ii2 - 2 * i0 - i3 + 5; i4 <= 64 * ii0 + 60 * ii2 - 2 * i0 - i3 + (4 * ii2 + i2 - 2) / 9 + 5; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        if (i2 == 32 * ii2 + 2 && 2 * i0 + i3 >= 64 * ii0 + 32 * ii2 + 5) {
                          B[32 * ii2 + 2][i3][64 * ii0 + 64 * ii2 - 2 * i0 - i3 + 6] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 3][i3][64 * ii0 + 64 * ii2 - 2 * i0 - i3 + 6] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][i3][64 * ii0 + 64 * ii2 - 2 * i0 - i3 + 6])) + A[32 * ii2 + 1][i3][64 * ii0 + 64 * ii2 - 2 * i0 - i3 + 6])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3 + 1][64 * ii0 + 64 * ii2 - 2 * i0 - i3 + 6] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][i3][64 * ii0 + 64 * ii2 - 2 * i0 - i3 + 6])) + A[32 * ii2 + 2][i3 - 1][64 * ii0 + 64 * ii2 - 2 * i0 - i3 + 6]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3][64 * ii0 + 64 * ii2 - 2 * i0 - i3 + 7] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][i3][64 * ii0 + 64 * ii2 - 2 * i0 - i3 + 6])) + A[32 * ii2 + 2][i3][64 * ii0 + 64 * ii2 - 2 * i0 - i3 + 5]))) + A[32 * ii2 + 2][i3][64 * ii0 + 64 * ii2 - 2 * i0 - i3 + 6]);
                        }
                        for (int i4 = max(64 * ii0 + 96 * ii2 - 2 * i0 - i2 - i3 + 9, 64 * ii0 + 60 * ii2 - 2 * i0 - i3 + (4 * ii2 + i2 - 2) / 9 + 6); i4 <= min(min(96 * ii2 - 2 * i3, 64 * ii0 + 64 * ii2 - 2 * i0 - i3 + 7), 11 * ii2 + i2 - (ii2 + i2 + 1) / 3); i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = 64 * ii0 + 64 * ii2 - 2 * i0 - i3 + 8; i4 <= min(96 * ii2 - 2 * i3, -32 * ii0 - 32 * ii2 + i0 + i2 + (i2 + i3) / 2 - 4); i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        if (32 * ii2 >= i3 + 1) {
                          for (int i4 = 96 * ii2 - 2 * i3 + 1; i4 < -32 * ii0 - 32 * ii2 + i0 + i2 + (i2 + i3) / 2 - 3; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        } else {
                          if (i0 >= 32 * ii0 + 4 && i2 == 32 * ii2 + 2 && i3 == 32 * ii2 + 1) {
                            B[32 * ii2 + 2][32 * ii2 + 1][32 * ii2 - 1] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 3][32 * ii2 + 1][32 * ii2 - 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][32 * ii2 + 1][32 * ii2 - 1])) + A[32 * ii2 + 1][32 * ii2 + 1][32 * ii2 - 1])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][32 * ii2 + 2][32 * ii2 - 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][32 * ii2 + 1][32 * ii2 - 1])) + A[32 * ii2 + 2][32 * ii2][32 * ii2 - 1]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][32 * ii2 + 1][32 * ii2] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][32 * ii2 + 1][32 * ii2 - 1])) + A[32 * ii2 + 2][32 * ii2 + 1][32 * ii2 - 2]))) + A[32 * ii2 + 2][32 * ii2 + 1][32 * ii2 - 1]);
                          }
                          for (int i4 = max(32 * ii2, 96 * ii2 - 2 * i3 + 1); i4 < -32 * ii0 - 32 * ii2 + i0 + i2 + (i2 + i3) / 2 - 3; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                        for (int i4 = -32 * ii0 - 32 * ii2 + i0 + i2 + (i2 + i3) / 2 - 3; i4 < _PB_N - 1; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                      for (int i3 = max(max(1, 64 * ii0 + 32 * ii2 - 2 * i0 + 4), 2 * _PB_N + 64 * ii0 + 64 * ii2 - 2 * i0 - 3 * i2 + 4); i3 <= min(min(32 * ii2 - 1, 64 * ii0 + 60 * ii2 - 2 * i0 + (4 * ii2 + i2 - 2) / 9 + 4), 36 * ii2 - (4 * ii2 + i2 + 6) / 9 + 1); i3 += 1) {
                        for (int i4 = max(1, 64 * ii0 + 64 * ii2 - 2 * i0 - i3 + 5); i4 <= 64 * ii0 + 60 * ii2 - 2 * i0 - i3 + (4 * ii2 + i2 - 2) / 9 + 5; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        if (i2 == 32 * ii2 + 2) {
                          B[32 * ii2 + 2][i3][64 * ii0 + 64 * ii2 - 2 * i0 - i3 + 6] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 3][i3][64 * ii0 + 64 * ii2 - 2 * i0 - i3 + 6] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][i3][64 * ii0 + 64 * ii2 - 2 * i0 - i3 + 6])) + A[32 * ii2 + 1][i3][64 * ii0 + 64 * ii2 - 2 * i0 - i3 + 6])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3 + 1][64 * ii0 + 64 * ii2 - 2 * i0 - i3 + 6] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][i3][64 * ii0 + 64 * ii2 - 2 * i0 - i3 + 6])) + A[32 * ii2 + 2][i3 - 1][64 * ii0 + 64 * ii2 - 2 * i0 - i3 + 6]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3][64 * ii0 + 64 * ii2 - 2 * i0 - i3 + 7] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][i3][64 * ii0 + 64 * ii2 - 2 * i0 - i3 + 6])) + A[32 * ii2 + 2][i3][64 * ii0 + 64 * ii2 - 2 * i0 - i3 + 5]))) + A[32 * ii2 + 2][i3][64 * ii0 + 64 * ii2 - 2 * i0 - i3 + 6]);
                        }
                        for (int i4 = max(64 * ii0 + 96 * ii2 - 2 * i0 - i2 - i3 + 9, 64 * ii0 + 60 * ii2 - 2 * i0 - i3 + (4 * ii2 + i2 - 2) / 9 + 6); i4 < _PB_N - 1; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                      for (int i3 = max(max(32 * ii2, 2 * _PB_N + 64 * ii0 + 64 * ii2 - 2 * i0 - 3 * i2 + 4), -64 * ii0 + 11 * ii2 + 2 * i0 - (ii2 + i2 + 1) / 3 - 1); i3 <= 36 * ii2 - (4 * ii2 + i2 + 6) / 9 + 1; i3 += 1) {
                        for (int i4 = 64 * ii0 + 64 * ii2 - 2 * i0 - i3 + 5; i4 <= 64 * ii0 + 60 * ii2 - 2 * i0 - i3 + (4 * ii2 + i2 - 2) / 9 + 5; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        if (i2 == 32 * ii2 + 2) {
                          B[32 * ii2 + 2][i3][64 * ii0 + 64 * ii2 - 2 * i0 - i3 + 6] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 3][i3][64 * ii0 + 64 * ii2 - 2 * i0 - i3 + 6] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][i3][64 * ii0 + 64 * ii2 - 2 * i0 - i3 + 6])) + A[32 * ii2 + 1][i3][64 * ii0 + 64 * ii2 - 2 * i0 - i3 + 6])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3 + 1][64 * ii0 + 64 * ii2 - 2 * i0 - i3 + 6] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][i3][64 * ii0 + 64 * ii2 - 2 * i0 - i3 + 6])) + A[32 * ii2 + 2][i3 - 1][64 * ii0 + 64 * ii2 - 2 * i0 - i3 + 6]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3][64 * ii0 + 64 * ii2 - 2 * i0 - i3 + 7] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][i3][64 * ii0 + 64 * ii2 - 2 * i0 - i3 + 6])) + A[32 * ii2 + 2][i3][64 * ii0 + 64 * ii2 - 2 * i0 - i3 + 5]))) + A[32 * ii2 + 2][i3][64 * ii0 + 64 * ii2 - 2 * i0 - i3 + 6]);
                        }
                        for (int i4 = max(64 * ii0 + 96 * ii2 - 2 * i0 - i2 - i3 + 9, 64 * ii0 + 60 * ii2 - 2 * i0 - i3 + (4 * ii2 + i2 - 2) / 9 + 6); i4 < _PB_N - 1; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                      if (ii2 == 1) {
                        for (int i3 = max(2 * _PB_N + 64 * ii0 - 2 * i0 - 3 * i2 + 68, 64 * ii0 - 2 * i0 + (i2 + 2) / 9 + 65); i3 <= -i2 + 64; i3 += 1) {
                          if (i2 == 34 && 2 * i0 + i3 == 64 * ii0 + 69) {
                            B[34][64 * ii0 - 2 * i0 + 69][1] = ((((SCALAR_VAL(0.125) * ((A[35][64 * ii0 - 2 * i0 + 69][1] - (SCALAR_VAL(2.0) * A[34][64 * ii0 - 2 * i0 + 69][1])) + A[33][64 * ii0 - 2 * i0 + 69][1])) + (SCALAR_VAL(0.125) * ((A[34][64 * ii0 - 2 * i0 + 70][1] - (SCALAR_VAL(2.0) * A[34][64 * ii0 - 2 * i0 + 69][1])) + A[34][64 * ii0 - 2 * i0 + 68][1]))) + (SCALAR_VAL(0.125) * ((A[34][64 * ii0 - 2 * i0 + 69][2] - (SCALAR_VAL(2.0) * A[34][64 * ii0 - 2 * i0 + 69][1])) + A[34][64 * ii0 - 2 * i0 + 69][0]))) + A[34][64 * ii0 - 2 * i0 + 69][1]);
                          }
                          for (int i4 = max(1, 64 * ii0 - 2 * i0 - i2 - i3 + 105); i4 < _PB_N - 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                    }
                    if (ii2 == 1) {
                      for (int i3 = max(max(-i2 + 65, 2 * _PB_N + 64 * ii0 - 2 * i0 - 3 * i2 + 68), 64 * ii0 - 2 * i0 + (i2 + 2) / 9 + 65); i3 <= min(31, -64 * ii0 + 2 * i0 - (i2 + 2) / 3 + 9); i3 += 1) {
                        if (i0 == 32 * ii0 + 19 && i2 == 34 && i3 == 31) {
                          B[34][31][1] = ((((SCALAR_VAL(0.125) * ((A[35][31][1] - (SCALAR_VAL(2.0) * A[34][31][1])) + A[33][31][1])) + (SCALAR_VAL(0.125) * ((A[34][32][1] - (SCALAR_VAL(2.0) * A[34][31][1])) + A[34][30][1]))) + (SCALAR_VAL(0.125) * ((A[34][31][2] - (SCALAR_VAL(2.0) * A[34][31][1])) + A[34][31][0]))) + A[34][31][1]);
                        }
                        for (int i4 = max(1, 64 * ii0 - 2 * i0 - i2 - i3 + 105); i4 < _PB_N - 1; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                      if (i2 >= 32) {
                        for (int i3 = max(max(32, -i2 + 65), 2 * _PB_N + 64 * ii0 - 2 * i0 - 3 * i2 + 68); i3 <= min(_PB_N - 2, -64 * ii0 + 2 * i0 - (i2 + 2) / 3 + 9); i3 += 1) {
                          if (i0 == 32 * ii0 + 18 && i2 >= 34 && i3 == 32) {
                            B[i2][32][1] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32][1] - (SCALAR_VAL(2.0) * A[i2][32][1])) + A[i2 - 1][32][1])) + (SCALAR_VAL(0.125) * ((A[i2][33][1] - (SCALAR_VAL(2.0) * A[i2][32][1])) + A[i2][31][1]))) + (SCALAR_VAL(0.125) * ((A[i2][32][2] - (SCALAR_VAL(2.0) * A[i2][32][1])) + A[i2][32][0]))) + A[i2][32][1]);
                          } else if (i2 == 32) {
                            for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 + i3 + 4; i4 += 1) {
                              B[32][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[33][i3][i4] - (SCALAR_VAL(2.0) * A[32][i3][i4])) + A[31][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32][i3][i4])) + A[32][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32][i3][i4])) + A[32][i3][i4 - 1]))) + A[32][i3][i4]);
                            }
                          }
                          for (int i4 = max(max(1, 64 * ii0 - 2 * i0 - 2 * i2 + i3 + (i2 - 1) / 3 + 59), 64 * ii0 - 2 * i0 - i3 + (i2 + 2) / 9 + 66); i4 < _PB_N - 1; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                    }
                    if (i2 >= 32 * ii2 + 2) {
                      for (int i3 = max(-64 * ii0 + 11 * ii2 + 2 * i0 - (ii2 + i2 + 1) / 3 - 1, 36 * ii2 - (4 * ii2 + i2 + 6) / 9 + 2); i3 < min(_PB_N - 1, -64 * ii0 + 2 * i0 + i2 - 2); i3 += 1) {
                        for (int i4 = max(max(1, 64 * ii0 + 32 * ii2 - 2 * i0 + 4), 64 * ii0 + 64 * ii2 - 2 * i0 - i3 + 5); i4 <= min(_PB_N - 2, 64 * ii0 - 11 * ii2 - 2 * i0 + i3 + (ii2 + i2 + 1) / 3 + 2); i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = 64 * ii0 - 11 * ii2 - 2 * i0 + i3 + (ii2 + i2 + 1) / 3 + 3; i4 < _PB_N - 1; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    }
                    for (int i3 = -64 * ii0 + 2 * i0 + i2 - 2; i3 < min(min(_PB_N - 1, _PB_N - 64 * ii0 - 53 * ii2 + 2 * i0 + 2 * i2 - (ii2 + i2 + 1) / 3 - 6), _PB_N - 64 * ii0 + 11 * ii2 + 2 * i0 - (ii2 + i2 + 1) / 3 - 4); i3 += 1) {
                      for (int i4 = max(64 * ii0 + 32 * ii2 - 2 * i0 + 4, 64 * ii0 + 64 * ii2 - 2 * i0 - i2 + 5); i4 < _PB_N - 1; i4 += 1) {
                        B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                      }
                    }
                    for (int i3 = max(-64 * ii0 + 2 * i0 + i2 - 2, _PB_N - 64 * ii0 - 53 * ii2 + 2 * i0 + 2 * i2 - (ii2 + i2 + 1) / 3 - 6); i3 < _PB_N - 1; i3 += 1) {
                      for (int i4 = 64 * ii0 + 64 * ii2 - 2 * i0 - i2 + 5; i4 < _PB_N - 1; i4 += 1) {
                        B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                      }
                    }
                    if (32 * ii2 + 3 == _PB_N && i0 == 32 * ii0 + 3 && i2 + 2 == _PB_N) {
                      for (int i4 = _PB_N - 4; i4 < _PB_N - 1; i4 += 1) {
                        B[_PB_N - 2][_PB_N - 3][i4] = ((((SCALAR_VAL(0.125) * ((A[_PB_N - 1][_PB_N - 3][i4] - (SCALAR_VAL(2.0) * A[_PB_N - 2][_PB_N - 3][i4])) + A[_PB_N - 3][_PB_N - 3][i4])) + (SCALAR_VAL(0.125) * ((A[_PB_N - 2][_PB_N - 2][i4] - (SCALAR_VAL(2.0) * A[_PB_N - 2][_PB_N - 3][i4])) + A[_PB_N - 2][_PB_N - 4][i4]))) + (SCALAR_VAL(0.125) * ((A[_PB_N - 2][_PB_N - 3][i4 + 1] - (SCALAR_VAL(2.0) * A[_PB_N - 2][_PB_N - 3][i4])) + A[_PB_N - 2][_PB_N - 3][i4 - 1]))) + A[_PB_N - 2][_PB_N - 3][i4]);
                      }
                    }
                    if (i2 >= 32 * ii2 && 32 * ii2 + 1 >= i2) {
                      for (int i3 = max(max(-64 * ii0 + 2 * i0 - 1, 64 * ii2 - i2 + 1), 2 * _PB_N + 64 * ii0 + 64 * ii2 - 2 * i0 - 3 * i2 + 4); i3 < _PB_N - 1; i3 += 1) {
                        for (int i4 = max(max(1, 64 * ii0 + 64 * ii2 - 2 * i0 - i2 + 5), 64 * ii0 + 96 * ii2 - 2 * i0 - i2 - i3 + 6); i4 < _PB_N - 1; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    }
                    for (int i3 = max(-64 * ii0 + 2 * i0 + i2 - 2, _PB_N - 64 * ii0 + 11 * ii2 + 2 * i0 - (ii2 + i2 + 1) / 3 - 4); i3 < _PB_N - 1; i3 += 1) {
                      for (int i4 = 64 * ii0 + 32 * ii2 - 2 * i0 + 4; i4 < _PB_N - 1; i4 += 1) {
                        B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                      }
                    }
                  } else {
                    for (int i3 = max(32 * ii2, 64 * ii2 - i2 + 1); i3 < _PB_N - 1; i3 += 1) {
                      for (int i4 = max(max(32 * ii2, 64 * ii2 - i2 + 1), 96 * ii2 - 2 * i3 + 1); i4 < _PB_N - 1; i4 += 1) {
                        B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                      }
                    }
                  }
                }
              }
              for (int i2 = max(1, 64 * ii0 + 32 * ii2 - 2 * i0 + 3); i2 <= 32 * ii2; i2 += 1) {
                for (int i3 = max(1, 64 * ii0 + 64 * ii2 - 2 * i0 - i2 + 4); i3 < _PB_N - 1; i3 += 1) {
                  for (int i4 = max(max(1, 64 * ii0 + 64 * ii2 - 2 * i0 - i2 + 4), 64 * ii0 + 96 * ii2 - 2 * i0 - i2 - i3 + 5); i4 < _PB_N - 1; i4 += 1) {
                    A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                  }
                }
              }
              for (int i2 = 32 * ii2 + 1; i2 < _PB_N - 1; i2 += 1) {
                for (int i3 = max(1, 64 * ii0 + 32 * ii2 - 2 * i0 + 3); i3 <= 32 * ii2; i3 += 1) {
                  for (int i4 = max(1, 64 * ii0 + 64 * ii2 - 2 * i0 - i3 + 4); i4 < _PB_N - 1; i4 += 1) {
                    A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                  }
                }
                for (int i3 = 32 * ii2 + 1; i3 < _PB_N - 1; i3 += 1) {
                  for (int i4 = max(1, 64 * ii0 + 32 * ii2 - 2 * i0 + 3); i4 < _PB_N - 1; i4 += 1) {
                    A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                  }
                }
              }
            }
          }
        }
      }
      if (_PB_N <= 36) {
        for (int ii3 = 0; ii3 <= 1; ii3 += 1) {
          for (int ii4 = 0; ii4 <= 1; ii4 += 1) {
            if (_PB_N == 35) {
              if (ii3 == 0 && ii4 == 0) {
                for (int i0 = 32 * ii0 + 1; i0 <= 32 * ii0 + 3; i0 += 1) {
                  if (i0 >= 32 * ii0 + 2) {
                    if (i0 == 32 * ii0 + 3) {
                      for (int i2 = 30; i2 <= 33; i2 += 1) {
                        if (i2 >= 32) {
                          for (int i3 = 1; i3 <= -i2 + 62; i3 += 1) {
                            if (2 * i2 + i3 >= 67 || 1 || 1) {
                              for (int i4 = 1; i4 <= -i2 + 62; i4 += 1) {
                                if (2 * i2 + i3 >= 67 || 1 || 1) {
                                  B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                                }
                              }
                            }
                          }
                        } else {
                          for (int i3 = 1; i3 <= 2; i3 += 1) {
                            if (i3 == 1) {
                              B[i2][1][1] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][1][1] - (SCALAR_VAL(2.0) * A[i2][1][1])) + A[i2 - 1][1][1])) + (SCALAR_VAL(0.125) * ((A[i2][2][1] - (SCALAR_VAL(2.0) * A[i2][1][1])) + A[i2][0][1]))) + (SCALAR_VAL(0.125) * ((A[i2][1][2] - (SCALAR_VAL(2.0) * A[i2][1][1])) + A[i2][1][0]))) + A[i2][1][1]);
                            }
                            for (int i4 = -i3 + 3; i4 <= -i2 + 62; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                          for (int i3 = 3; i3 <= -i2 + 62; i3 += 1) {
                            for (int i4 = 1; i4 <= -2 * i2 + 92; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            if (i2 == 31) {
                              B[31][i3][31] = ((((SCALAR_VAL(0.125) * ((A[32][i3][31] - (SCALAR_VAL(2.0) * A[31][i3][31])) + A[30][i3][31])) + (SCALAR_VAL(0.125) * ((A[31][i3 + 1][31] - (SCALAR_VAL(2.0) * A[31][i3][31])) + A[31][i3 - 1][31]))) + (SCALAR_VAL(0.125) * ((A[31][i3][32] - (SCALAR_VAL(2.0) * A[31][i3][31])) + A[31][i3][30]))) + A[31][i3][31]);
                            }
                          }
                        }
                      }
                    } else {
                      for (int i2 = 32; i2 <= 33; i2 += 1) {
                        for (int i3 = 1; i3 <= -i2 + 64; i3 += 1) {
                          for (int i4 = 1; i4 <= -i2 + 64; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                    }
                  }
                  for (int i2 = 64 * ii0 - 2 * i0 + 35; i2 <= 33; i2 += 1) {
                    for (int i3 = 1; i3 <= 64 * ii0 - 2 * i0 - i2 + 67; i3 += 1) {
                      for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 - i2 + 67; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                  }
                }
              }
              for (int i0 = max(32 * ii0 + 1, 32 * ii0 - 3 * ii3 - 3 * ii4 + 4); i0 <= min(min(_PB_TSTEPS, 32 * ii0 + 32), 32 * ii0 + 16 * ii3 + 15 * ii4 + 17); i0 += 1) {
                if (3 * ii3 + ii4 + i0 >= 32 * ii0 + 4 && i0 >= 32 * ii0 + ii3 + ii4 + 1) {
                  if (ii3 == 0) {
                    if (ii4 == 0) {
                      for (int i3 = 1; i3 <= 32; i3 += 1) {
                        if (i3 == 1) {
                          B[64 * ii0 - 2 * i0 + 36][1][1] = ((((SCALAR_VAL(0.125) * ((A[64 * ii0 - 2 * i0 + 37][1][1] - (SCALAR_VAL(2.0) * A[64 * ii0 - 2 * i0 + 36][1][1])) + A[64 * ii0 - 2 * i0 + 35][1][1])) + (SCALAR_VAL(0.125) * ((A[64 * ii0 - 2 * i0 + 36][2][1] - (SCALAR_VAL(2.0) * A[64 * ii0 - 2 * i0 + 36][1][1])) + A[64 * ii0 - 2 * i0 + 36][0][1]))) + (SCALAR_VAL(0.125) * ((A[64 * ii0 - 2 * i0 + 36][1][2] - (SCALAR_VAL(2.0) * A[64 * ii0 - 2 * i0 + 36][1][1])) + A[64 * ii0 - 2 * i0 + 36][1][0]))) + A[64 * ii0 - 2 * i0 + 36][1][1]);
                        }
                        for (int i4 = max(1, -i3 + 3); i4 <= 32; i4 += 1) {
                          B[64 * ii0 - 2 * i0 + 36][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii0 - 2 * i0 + 37][i3][i4] - (SCALAR_VAL(2.0) * A[64 * ii0 - 2 * i0 + 36][i3][i4])) + A[64 * ii0 - 2 * i0 + 35][i3][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii0 - 2 * i0 + 36][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[64 * ii0 - 2 * i0 + 36][i3][i4])) + A[64 * ii0 - 2 * i0 + 36][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii0 - 2 * i0 + 36][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii0 - 2 * i0 + 36][i3][i4])) + A[64 * ii0 - 2 * i0 + 36][i3][i4 - 1]))) + A[64 * ii0 - 2 * i0 + 36][i3][i4]);
                        }
                      }
                      for (int i3 = 1; i3 <= 31; i3 += 1) {
                        if (i3 == 1) {
                          B[64 * ii0 - 2 * i0 + 37][1][1] = ((((SCALAR_VAL(0.125) * ((A[64 * ii0 - 2 * i0 + 38][1][1] - (SCALAR_VAL(2.0) * A[64 * ii0 - 2 * i0 + 37][1][1])) + A[64 * ii0 - 2 * i0 + 36][1][1])) + (SCALAR_VAL(0.125) * ((A[64 * ii0 - 2 * i0 + 37][2][1] - (SCALAR_VAL(2.0) * A[64 * ii0 - 2 * i0 + 37][1][1])) + A[64 * ii0 - 2 * i0 + 37][0][1]))) + (SCALAR_VAL(0.125) * ((A[64 * ii0 - 2 * i0 + 37][1][2] - (SCALAR_VAL(2.0) * A[64 * ii0 - 2 * i0 + 37][1][1])) + A[64 * ii0 - 2 * i0 + 37][1][0]))) + A[64 * ii0 - 2 * i0 + 37][1][1]);
                        }
                        for (int i4 = max(1, -i3 + 3); i4 <= 31; i4 += 1) {
                          B[64 * ii0 - 2 * i0 + 37][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii0 - 2 * i0 + 38][i3][i4] - (SCALAR_VAL(2.0) * A[64 * ii0 - 2 * i0 + 37][i3][i4])) + A[64 * ii0 - 2 * i0 + 36][i3][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii0 - 2 * i0 + 37][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[64 * ii0 - 2 * i0 + 37][i3][i4])) + A[64 * ii0 - 2 * i0 + 37][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii0 - 2 * i0 + 37][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii0 - 2 * i0 + 37][i3][i4])) + A[64 * ii0 - 2 * i0 + 37][i3][i4 - 1]))) + A[64 * ii0 - 2 * i0 + 37][i3][i4]);
                        }
                      }
                    }
                    for (int i2 = max(1, 64 * ii0 - 2 * ii4 - 2 * i0 + 38); i2 <= -30 * ii4 + 31; i2 += 1) {
                      if (ii4 == 0) {
                        for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 - i2 + 68; i4 += 1) {
                          B[i2][1][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][1][i4] - (SCALAR_VAL(2.0) * A[i2][1][i4])) + A[i2 - 1][1][i4])) + (SCALAR_VAL(0.125) * ((A[i2][2][i4] - (SCALAR_VAL(2.0) * A[i2][1][i4])) + A[i2][0][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][1][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][1][i4])) + A[i2][1][i4 - 1]))) + A[i2][1][i4]);
                        }
                        if (2 * i0 + i2 == 64 * ii0 + 38) {
                          for (int i4 = 1; i4 <= 30; i4 += 1) {
                            B[64 * ii0 - 2 * i0 + 38][2][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii0 - 2 * i0 + 39][2][i4] - (SCALAR_VAL(2.0) * A[64 * ii0 - 2 * i0 + 38][2][i4])) + A[64 * ii0 - 2 * i0 + 37][2][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii0 - 2 * i0 + 38][3][i4] - (SCALAR_VAL(2.0) * A[64 * ii0 - 2 * i0 + 38][2][i4])) + A[64 * ii0 - 2 * i0 + 38][1][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii0 - 2 * i0 + 38][2][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii0 - 2 * i0 + 38][2][i4])) + A[64 * ii0 - 2 * i0 + 38][2][i4 - 1]))) + A[64 * ii0 - 2 * i0 + 38][2][i4]);
                          }
                        }
                        for (int i3 = max(2, 64 * ii0 - 2 * i0 - i2 + 41); i3 <= 64 * ii0 - 2 * i0 - i2 + 68; i3 += 1) {
                          for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 - i2 + 68; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      } else {
                        for (int i3 = 1; i3 <= 64 * ii0 - 2 * i0 + 67; i3 += 1) {
                          for (int i4 = 64 * ii0 - 2 * i0 + 68; i4 <= 33; i4 += 1) {
                            B[1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[2][i3][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[0][i3][i4])) + (SCALAR_VAL(0.125) * ((A[1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3][i4 - 1]))) + A[1][i3][i4]);
                          }
                        }
                      }
                    }
                    if (ii4 == 0 && i0 == 32 * ii0 + 17) {
                      for (int i3 = 1; i3 <= 2; i3 += 1) {
                        for (int i4 = 1; i4 <= 2; i4 += 1) {
                          B[32][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[33][i3][i4] - (SCALAR_VAL(2.0) * A[32][i3][i4])) + A[31][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32][i3][i4])) + A[32][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32][i3][i4])) + A[32][i3][i4 - 1]))) + A[32][i3][i4]);
                        }
                      }
                      B[33][1][1] = ((((SCALAR_VAL(0.125) * ((A[34][1][1] - (SCALAR_VAL(2.0) * A[33][1][1])) + A[32][1][1])) + (SCALAR_VAL(0.125) * ((A[33][2][1] - (SCALAR_VAL(2.0) * A[33][1][1])) + A[33][0][1]))) + (SCALAR_VAL(0.125) * ((A[33][1][2] - (SCALAR_VAL(2.0) * A[33][1][1])) + A[33][1][0]))) + A[33][1][1]);
                    } else if (ii4 == 1) {
                      for (int i2 = max(2, 64 * ii0 - 2 * i0 + 36); i2 <= min(31, 64 * ii0 - 2 * i0 + 65); i2 += 1) {
                        for (int i3 = 1; i3 <= 64 * ii0 - 2 * i0 - i2 + 68; i3 += 1) {
                          if (i3 <= 2) {
                            for (int i4 = 64 * ii0 - 2 * i0 - i2 + 69; i4 <= 33; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          } else {
                            if (2 * i0 + i2 >= 64 * ii0 + 38) {
                              for (int i4 = 64 * ii0 - 2 * i0 - i2 + 69; i4 <= 33; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                            if (i3 <= 31) {
                              for (int i4 = max(-128 * ii0 + 4 * i0 + 2 * i2 - 42, 128 * ii0 - 4 * i0 - 2 * i2 + 105); i4 <= 33; i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            } else {
                              B[64 * ii0 - 2 * i0 + 36][32][33] = ((((SCALAR_VAL(0.125) * ((A[64 * ii0 - 2 * i0 + 37][32][33] - (SCALAR_VAL(2.0) * A[64 * ii0 - 2 * i0 + 36][32][33])) + A[64 * ii0 - 2 * i0 + 35][32][33])) + (SCALAR_VAL(0.125) * ((A[64 * ii0 - 2 * i0 + 36][33][33] - (SCALAR_VAL(2.0) * A[64 * ii0 - 2 * i0 + 36][32][33])) + A[64 * ii0 - 2 * i0 + 36][31][33]))) + (SCALAR_VAL(0.125) * ((A[64 * ii0 - 2 * i0 + 36][32][34] - (SCALAR_VAL(2.0) * A[64 * ii0 - 2 * i0 + 36][32][33])) + A[64 * ii0 - 2 * i0 + 36][32][32]))) + A[64 * ii0 - 2 * i0 + 36][32][33]);
                            }
                          }
                        }
                      }
                      for (int i2 = 64 * ii0 - 2 * i0 + 66; i2 <= min(31, 64 * ii0 - 2 * i0 + 67); i2 += 1) {
                        for (int i3 = 1; i3 <= 64 * ii0 - 2 * i0 - i2 + 68; i3 += 1) {
                          for (int i4 = 64 * ii0 - 2 * i0 - i2 + 69; i4 <= 33; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                      if (32 * ii0 + 17 >= i0) {
                        for (int i2 = 32; i2 <= 33; i2 += 1) {
                          for (int i3 = 1; i3 <= 64 * ii0 - 2 * i0 - i2 + 68; i3 += 1) {
                            if (i2 == 33 && i3 == 1) {
                              B[33][1][64 * ii0 - 2 * i0 + 36] = ((((SCALAR_VAL(0.125) * ((A[34][1][64 * ii0 - 2 * i0 + 36] - (SCALAR_VAL(2.0) * A[33][1][64 * ii0 - 2 * i0 + 36])) + A[32][1][64 * ii0 - 2 * i0 + 36])) + (SCALAR_VAL(0.125) * ((A[33][2][64 * ii0 - 2 * i0 + 36] - (SCALAR_VAL(2.0) * A[33][1][64 * ii0 - 2 * i0 + 36])) + A[33][0][64 * ii0 - 2 * i0 + 36]))) + (SCALAR_VAL(0.125) * ((A[33][1][64 * ii0 - 2 * i0 + 37] - (SCALAR_VAL(2.0) * A[33][1][64 * ii0 - 2 * i0 + 36])) + A[33][1][64 * ii0 - 2 * i0 + 35]))) + A[33][1][64 * ii0 - 2 * i0 + 36]);
                            }
                            for (int i4 = max(64 * ii0 - 2 * i0 - i2 + 69, 64 * ii0 - 2 * i0 + i2 - i3 + 5); i4 <= 33; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                        }
                      }
                    }
                    if (ii4 == 0) {
                      for (int i2 = 32; i2 <= min(33, 64 * ii0 - 2 * i0 + 65); i2 += 1) {
                        for (int i3 = 1; i3 <= 64 * ii0 - 2 * i0 - i2 + 68; i3 += 1) {
                          for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 - i2 + 68; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                    }
                  } else if (ii4 == 1) {
                    for (int i2 = max(1, 64 * ii0 - 2 * i0 + 36); i2 <= 48 * ii0 - 2 * i0 + i0 / 2 + 35; i2 += 1) {
                      for (int i3 = 64 * ii0 - 2 * i0 - i2 + 69; i3 <= 33; i3 += 1) {
                        for (int i4 = 64 * ii0 - 2 * i0 - i2 - i3 + 102; i4 <= 33; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    }
                    for (int i2 = max(1, 48 * ii0 - 2 * i0 + i0 / 2 + 36); i2 <= min(32, 21 * ii0 - i0 + (ii0 + i0 + 1) / 3 + 45); i2 += 1) {
                      for (int i3 = max(1, 64 * ii0 - 2 * i0 - i2 + 69); i3 <= min(33, -64 * ii0 + 2 * i0 + 2 * i2 - (i2 + 2) / 3 - 25); i3 += 1) {
                        for (int i4 = max(1, 64 * ii0 - 2 * i0 - i2 - i3 + 102); i4 <= 64 * ii0 - 2 * i0 - 2 * i2 + i3 + (i2 - 1) / 3 + 58; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = max(max(1, 64 * ii0 - 2 * i0 - i2 - i3 + 102), 64 * ii0 - 2 * i0 - 2 * i2 + i3 + (i2 - 1) / 3 + 59); i4 <= min(64 * ii0 - 2 * i0 - i2 - i3 + 104, i2 - (i2 + 2) / 3 + 11); i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = max(max(1, 64 * ii0 - 2 * i0 - i2 - i3 + 105), 64 * ii0 - 2 * i0 - 2 * i2 + i3 + (i2 - 1) / 3 + 59); i4 <= min(32, -32 * ii0 + i0 + i2 + (i2 + i3) / 2 - 36); i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        if (i2 <= 31 && 2 * i0 + 3 * i2 >= 64 * ii0 + 105 && i3 == 33) {
                          B[i2][33][33] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][33][33] - (SCALAR_VAL(2.0) * A[i2][33][33])) + A[i2 - 1][33][33])) + (SCALAR_VAL(0.125) * ((A[i2][34][33] - (SCALAR_VAL(2.0) * A[i2][33][33])) + A[i2][32][33]))) + (SCALAR_VAL(0.125) * ((A[i2][33][34] - (SCALAR_VAL(2.0) * A[i2][33][33])) + A[i2][33][32]))) + A[i2][33][33]);
                        } else if (2 * i0 + 3 * i2 + i3 >= 64 * ii0 + 138 && i3 <= 32) {
                          B[i2][i3][33] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][33] - (SCALAR_VAL(2.0) * A[i2][i3][33])) + A[i2 - 1][i3][33])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][33] - (SCALAR_VAL(2.0) * A[i2][i3][33])) + A[i2][i3 - 1][33]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][34] - (SCALAR_VAL(2.0) * A[i2][i3][33])) + A[i2][i3][32]))) + A[i2][i3][33]);
                        } else if (i0 >= 32 * ii0 + 5 && i2 == 32 && i3 == 33) {
                          B[32][33][33] = ((((SCALAR_VAL(0.125) * ((A[33][33][33] - (SCALAR_VAL(2.0) * A[32][33][33])) + A[31][33][33])) + (SCALAR_VAL(0.125) * ((A[32][34][33] - (SCALAR_VAL(2.0) * A[32][33][33])) + A[32][32][33]))) + (SCALAR_VAL(0.125) * ((A[32][33][34] - (SCALAR_VAL(2.0) * A[32][33][33])) + A[32][33][32]))) + A[32][33][33]);
                        }
                        for (int i4 = max(max(max(64 * ii0 - 2 * i0 - i2 - i3 + 102, 64 * ii0 - 2 * i0 - 2 * i2 + i3 + (i2 - 1) / 3 + 59), i2 - (i2 + 2) / 3 + 12), -32 * ii0 + i0 + i2 + (i2 + i3) / 2 - 35); i4 <= 33; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                      for (int i3 = -64 * ii0 + 2 * i0 + 2 * i2 - (i2 + 2) / 3 - 24; i3 <= 33; i3 += 1) {
                        for (int i4 = 64 * ii0 - 2 * i0 - i2 - i3 + 102; i4 <= 33; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    }
                    for (int i2 = 21 * ii0 - i0 + (ii0 + i0 + 1) / 3 + 46; i2 <= 32; i2 += 1) {
                      for (int i3 = 1; i3 <= 33; i3 += 1) {
                        for (int i4 = max(1, 64 * ii0 - 2 * i0 - i2 - i3 + 102); i4 <= 33; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    }
                    for (int i3 = max(1, 64 * ii0 - 2 * i0 + 36); i3 <= 33; i3 += 1) {
                      for (int i4 = max(1, 64 * ii0 - 2 * i0 - i3 + 69); i4 <= 33; i4 += 1) {
                        B[33][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[34][i3][i4] - (SCALAR_VAL(2.0) * A[33][i3][i4])) + A[32][i3][i4])) + (SCALAR_VAL(0.125) * ((A[33][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[33][i3][i4])) + A[33][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[33][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[33][i3][i4])) + A[33][i3][i4 - 1]))) + A[33][i3][i4]);
                      }
                    }
                  } else if (i0 >= 32 * ii0 + 3) {
                    for (int i2 = max(1, 64 * ii0 - 2 * i0 + 36); i2 <= min(33, 64 * ii0 - 2 * i0 + 67); i2 += 1) {
                      for (int i3 = 64 * ii0 - 2 * i0 - i2 + 69; i3 <= 32; i3 += 1) {
                        if (64 * ii0 + 37 >= 2 * i0 + i3) {
                          B[i2][i3][1] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][1] - (SCALAR_VAL(2.0) * A[i2][i3][1])) + A[i2 - 1][i3][1])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][1] - (SCALAR_VAL(2.0) * A[i2][i3][1])) + A[i2][i3 - 1][1]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][2] - (SCALAR_VAL(2.0) * A[i2][i3][1])) + A[i2][i3][0]))) + A[i2][i3][1]);
                          if (i2 == 33 && 2 * i0 + i3 == 64 * ii0 + 36) {
                            for (int i4 = 2; i4 <= 32; i4 += 1) {
                              B[33][64 * ii0 - 2 * i0 + 36][i4] = ((((SCALAR_VAL(0.125) * ((A[34][64 * ii0 - 2 * i0 + 36][i4] - (SCALAR_VAL(2.0) * A[33][64 * ii0 - 2 * i0 + 36][i4])) + A[32][64 * ii0 - 2 * i0 + 36][i4])) + (SCALAR_VAL(0.125) * ((A[33][64 * ii0 - 2 * i0 + 37][i4] - (SCALAR_VAL(2.0) * A[33][64 * ii0 - 2 * i0 + 36][i4])) + A[33][64 * ii0 - 2 * i0 + 35][i4]))) + (SCALAR_VAL(0.125) * ((A[33][64 * ii0 - 2 * i0 + 36][i4 + 1] - (SCALAR_VAL(2.0) * A[33][64 * ii0 - 2 * i0 + 36][i4])) + A[33][64 * ii0 - 2 * i0 + 36][i4 - 1]))) + A[33][64 * ii0 - 2 * i0 + 36][i4]);
                            }
                          }
                        }
                        if (2 * i0 + i3 >= 64 * ii0 + 37) {
                          for (int i4 = max(1, 64 * ii0 - 2 * i0 - i3 + 39); i4 <= 64 * ii0 - 2 * i0 - i2 - i3 + 101; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                      for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 - i2 + 68; i4 += 1) {
                        B[i2][33][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][33][i4] - (SCALAR_VAL(2.0) * A[i2][33][i4])) + A[i2 - 1][33][i4])) + (SCALAR_VAL(0.125) * ((A[i2][34][i4] - (SCALAR_VAL(2.0) * A[i2][33][i4])) + A[i2][32][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][33][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][33][i4])) + A[i2][33][i4 - 1]))) + A[i2][33][i4]);
                      }
                    }
                    if (i0 >= 32 * ii0 + 19) {
                      for (int i2 = 64 * ii0 - 2 * i0 + 68; i2 <= 33; i2 += 1) {
                        for (int i3 = 1; i3 <= 64 * ii0 - 2 * i0 - i2 + 100; i3 += 1) {
                          for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 - i2 - i3 + 101; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                    } else if (i0 == 32 * ii0 + 18) {
                      for (int i2 = 32; i2 <= 33; i2 += 1) {
                        for (int i3 = 1; i3 <= -i2 + 64; i3 += 1) {
                          if (i3 == 1) {
                            B[i2][1][1] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][1][1] - (SCALAR_VAL(2.0) * A[i2][1][1])) + A[i2 - 1][1][1])) + (SCALAR_VAL(0.125) * ((A[i2][2][1] - (SCALAR_VAL(2.0) * A[i2][1][1])) + A[i2][0][1]))) + (SCALAR_VAL(0.125) * ((A[i2][1][2] - (SCALAR_VAL(2.0) * A[i2][1][1])) + A[i2][1][0]))) + A[i2][1][1]);
                          }
                          for (int i4 = max(1, -i3 + 3); i4 <= -i2 - i3 + 65; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                    }
                  } else {
                    for (int i2 = 32; i2 <= 33; i2 += 1) {
                      for (int i3 = -i2 + 65; i3 <= 33; i3 += 1) {
                        for (int i4 = 1; i4 <= -i2 + 64; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        if (i2 == 33 && i3 == 32) {
                          B[33][32][32] = ((((SCALAR_VAL(0.125) * ((A[34][32][32] - (SCALAR_VAL(2.0) * A[33][32][32])) + A[32][32][32])) + (SCALAR_VAL(0.125) * ((A[33][33][32] - (SCALAR_VAL(2.0) * A[33][32][32])) + A[33][31][32]))) + (SCALAR_VAL(0.125) * ((A[33][32][33] - (SCALAR_VAL(2.0) * A[33][32][32])) + A[33][32][31]))) + A[33][32][32]);
                        }
                      }
                    }
                  }
                } else if (ii3 == 1 && ii4 == 1 && i0 == 32 * ii0 + 2) {
                  for (int i2 = 32; i2 <= 33; i2 += 1) {
                    for (int i3 = -i2 + 65; i3 <= 33; i3 += 1) {
                      for (int i4 = max(-i2 + 65, -2 * i3 + 97); i4 <= 33; i4 += 1) {
                        B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                      }
                    }
                  }
                } else if (ii3 == 0 && ii4 == 1 && i0 == 32 * ii0 + 2) {
                  for (int i2 = 32; i2 <= 33; i2 += 1) {
                    for (int i3 = 1; i3 <= -i2 + 64; i3 += 1) {
                      for (int i4 = -i2 + 65; i4 <= 33; i4 += 1) {
                        B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                      }
                    }
                  }
                }
                if (ii3 == 1) {
                  for (int i2 = max(1, 64 * ii0 - 2 * i0 + 35); i2 <= 32; i2 += 1) {
                    for (int i3 = max(1, 64 * ii0 - 2 * i0 - i2 + 68); i3 <= min(33, 64 * ii0 + 30 * ii4 - 2 * i0 - i2 + 99); i3 += 1) {
                      for (int i4 = max(1, 64 * ii0 + 32 * ii4 - 2 * i0 - i2 - i3 + 69); i4 <= min(33, 64 * ii0 + 62 * ii4 - 2 * i0 - i2 - i3 + 100); i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                  }
                  if (ii4 == 0 && i0 >= 32 * ii0 + 17) {
                    for (int i3 = 1; i3 <= 64 * ii0 - 2 * i0 + 66; i3 += 1) {
                      for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 - i3 + 67; i4 += 1) {
                        A[33][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[34][i3][i4] - (SCALAR_VAL(2.0) * B[33][i3][i4])) + B[32][i3][i4])) + (SCALAR_VAL(0.125) * ((B[33][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[33][i3][i4])) + B[33][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[33][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[33][i3][i4])) + B[33][i3][i4 - 1]))) + B[33][i3][i4]);
                      }
                    }
                  }
                } else if (ii4 == 1) {
                  for (int i2 = max(1, 64 * ii0 - 2 * i0 + 35); i2 <= min(32, 64 * ii0 - 2 * i0 + 66); i2 += 1) {
                    for (int i3 = 1; i3 <= 64 * ii0 - 2 * i0 - i2 + 67; i3 += 1) {
                      for (int i4 = 64 * ii0 - 2 * i0 - i2 + 68; i4 <= 33; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                  }
                }
                if (ii4 == ii3) {
                  for (int i2 = max(32 * ii3 + 1, 64 * ii0 - 2 * i0 + 35); i2 <= ii3 + 32; i2 += 1) {
                    for (int i3 = max(1, 64 * ii0 + 26 * ii3 - 2 * i0 + 9); i3 <= min(32, 64 * ii0 + 63 * ii3 - 2 * i0 - i2 + 67); i3 += 1) {
                      if (ii3 == 0) {
                        for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 - i2 + 67; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      } else {
                        for (int i4 = max(1, 64 * ii0 - 2 * i0 - i3 + 68); i4 <= 33; i4 += 1) {
                          A[33][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[34][i3][i4] - (SCALAR_VAL(2.0) * B[33][i3][i4])) + B[32][i3][i4])) + (SCALAR_VAL(0.125) * ((B[33][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[33][i3][i4])) + B[33][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[33][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[33][i3][i4])) + B[33][i3][i4 - 1]))) + B[33][i3][i4]);
                        }
                      }
                    }
                    if (ii3 == 1 && i2 == 33) {
                      for (int i4 = max(1, 64 * ii0 - 2 * i0 + 35); i4 <= 33; i4 += 1) {
                        A[33][33][i4] = ((((SCALAR_VAL(0.125) * ((B[34][33][i4] - (SCALAR_VAL(2.0) * B[33][33][i4])) + B[32][33][i4])) + (SCALAR_VAL(0.125) * ((B[33][34][i4] - (SCALAR_VAL(2.0) * B[33][33][i4])) + B[33][32][i4]))) + (SCALAR_VAL(0.125) * ((B[33][33][i4 + 1] - (SCALAR_VAL(2.0) * B[33][33][i4])) + B[33][33][i4 - 1]))) + B[33][33][i4]);
                      }
                    }
                  }
                }
                if (ii3 + ii4 <= 1 && 32 * ii0 + 16 >= i0) {
                  if (ii4 == 0) {
                    for (int i3 = max(1, 64 * ii0 + 30 * ii3 - 2 * i0 + 5); i3 <= min(32, 64 * ii0 + 30 * ii3 - 2 * i0 + 34); i3 += 1) {
                      for (int i4 = 1; i4 <= min(64 * ii0 + 30 * ii3 - 2 * i0 + 34, 64 * ii0 - 2 * i0 - i3 + 67); i4 += 1) {
                        A[33][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[34][i3][i4] - (SCALAR_VAL(2.0) * B[33][i3][i4])) + B[32][i3][i4])) + (SCALAR_VAL(0.125) * ((B[33][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[33][i3][i4])) + B[33][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[33][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[33][i3][i4])) + B[33][i3][i4 - 1]))) + B[33][i3][i4]);
                      }
                    }
                    if (ii3 == 1) {
                      for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 + 34; i4 += 1) {
                        A[33][33][i4] = ((((SCALAR_VAL(0.125) * ((B[34][33][i4] - (SCALAR_VAL(2.0) * B[33][33][i4])) + B[32][33][i4])) + (SCALAR_VAL(0.125) * ((B[33][34][i4] - (SCALAR_VAL(2.0) * B[33][33][i4])) + B[33][32][i4]))) + (SCALAR_VAL(0.125) * ((B[33][33][i4 + 1] - (SCALAR_VAL(2.0) * B[33][33][i4])) + B[33][33][i4 - 1]))) + B[33][33][i4]);
                      }
                    }
                  } else {
                    for (int i3 = 1; i3 <= 64 * ii0 - 2 * i0 + 34; i3 += 1) {
                      for (int i4 = 64 * ii0 - 2 * i0 + 35; i4 <= 33; i4 += 1) {
                        A[33][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[34][i3][i4] - (SCALAR_VAL(2.0) * B[33][i3][i4])) + B[32][i3][i4])) + (SCALAR_VAL(0.125) * ((B[33][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[33][i3][i4])) + B[33][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[33][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[33][i3][i4])) + B[33][i3][i4 - 1]))) + B[33][i3][i4]);
                      }
                    }
                  }
                }
              }
              if (ii3 == 0 && ii4 == 0) {
                for (int i0 = 32 * ii0 + 18; i0 <= min(_PB_TSTEPS, 32 * ii0 + 32); i0 += 1) {
                  for (int i2 = 1; i2 <= 64 * ii0 - 2 * i0 + 67; i2 += 1) {
                    for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 - i2 + 68; i4 += 1) {
                      B[i2][1][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][1][i4] - (SCALAR_VAL(2.0) * A[i2][1][i4])) + A[i2 - 1][1][i4])) + (SCALAR_VAL(0.125) * ((A[i2][2][i4] - (SCALAR_VAL(2.0) * A[i2][1][i4])) + A[i2][0][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][1][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][1][i4])) + A[i2][1][i4 - 1]))) + A[i2][1][i4]);
                    }
                    if (i0 == 32 * ii0 + 18 && i2 <= 2) {
                      for (int i4 = 1; i4 <= -i2 + 32; i4 += 1) {
                        B[i2][2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][2][i4] - (SCALAR_VAL(2.0) * A[i2][2][i4])) + A[i2 - 1][2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][3][i4] - (SCALAR_VAL(2.0) * A[i2][2][i4])) + A[i2][1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][2][i4])) + A[i2][2][i4 - 1]))) + A[i2][2][i4]);
                      }
                    }
                    if (2 * i0 + i2 >= 64 * ii0 + 38) {
                      for (int i3 = max(2, 64 * ii0 - 2 * i0 - i2 + 41); i3 <= 64 * ii0 - 2 * i0 - i2 + 68; i3 += 1) {
                        for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 - i2 + 68; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    } else {
                      for (int i3 = 3; i3 <= 31; i3 += 1) {
                        for (int i4 = 1; i4 <= 31; i4 += 1) {
                          B[1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[2][i3][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[0][i3][i4])) + (SCALAR_VAL(0.125) * ((A[1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3][i4 - 1]))) + A[1][i3][i4]);
                        }
                      }
                    }
                  }
                  for (int i2 = 1; i2 <= 64 * ii0 - 2 * i0 + 66; i2 += 1) {
                    for (int i3 = 1; i3 <= 64 * ii0 - 2 * i0 - i2 + 67; i3 += 1) {
                      for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 - i2 + 67; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                  }
                }
              }
            } else if (ii3 >= ii4) {
              if (ii3 == 0 && ii4 == 0) {
                for (int i0 = 32 * ii0 + 1; i0 <= 32 * ii0 + 2; i0 += 1) {
                  if (i0 == 32 * ii0 + 2) {
                    for (int i2 = 32; i2 <= 34; i2 += 1) {
                      for (int i3 = 1; i3 <= -i2 + 64; i3 += 1) {
                        for (int i4 = 1; i4 <= -i2 + 64; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        if (i2 == 34) {
                          B[34][i3][31] = ((((SCALAR_VAL(0.125) * ((A[35][i3][31] - (SCALAR_VAL(2.0) * A[34][i3][31])) + A[33][i3][31])) + (SCALAR_VAL(0.125) * ((A[34][i3 + 1][31] - (SCALAR_VAL(2.0) * A[34][i3][31])) + A[34][i3 - 1][31]))) + (SCALAR_VAL(0.125) * ((A[34][i3][32] - (SCALAR_VAL(2.0) * A[34][i3][31])) + A[34][i3][30]))) + A[34][i3][31]);
                        }
                      }
                      if (i2 == 34) {
                        for (int i4 = 1; i4 <= 31; i4 += 1) {
                          B[34][31][i4] = ((((SCALAR_VAL(0.125) * ((A[35][31][i4] - (SCALAR_VAL(2.0) * A[34][31][i4])) + A[33][31][i4])) + (SCALAR_VAL(0.125) * ((A[34][32][i4] - (SCALAR_VAL(2.0) * A[34][31][i4])) + A[34][30][i4]))) + (SCALAR_VAL(0.125) * ((A[34][31][i4 + 1] - (SCALAR_VAL(2.0) * A[34][31][i4])) + A[34][31][i4 - 1]))) + A[34][31][i4]);
                        }
                      }
                    }
                  }
                  if (i0 == 32 * ii0 + 2) {
                    for (int i2 = 31; i2 <= 32; i2 += 1) {
                      for (int i3 = 1; i3 <= -i2 + 63; i3 += 1) {
                        for (int i4 = 1; i4 <= -i2 + 63; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                  }
                  for (int i2 = 33; i2 <= 34; i2 += 1) {
                    for (int i3 = 1; i3 <= 64 * ii0 - 2 * i0 + 34; i3 += 1) {
                      for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 + 34; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                  }
                }
              } else {
                for (int i0 = 32 * ii0 + 1; i0 <= 32 * ii0 - ii4 + 3; i0 += 1) {
                  if (i0 >= 32 * ii0 + 2) {
                    if (ii4 == 0 && i0 == 32 * ii0 + 3) {
                      for (int i3 = 33; i3 <= 34; i3 += 1) {
                        for (int i4 = 1; i4 <= 32; i4 += 1) {
                          B[30][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[31][i3][i4] - (SCALAR_VAL(2.0) * A[30][i3][i4])) + A[29][i3][i4])) + (SCALAR_VAL(0.125) * ((A[30][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[30][i3][i4])) + A[30][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[30][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[30][i3][i4])) + A[30][i3][i4 - 1]))) + A[30][i3][i4]);
                        }
                      }
                    }
                    for (int i2 = max(ii4 + 31, 32 * ii0 - ii4 - i0 + 34); i2 <= 34; i2 += 1) {
                      if (i0 == 32 * ii0 + 2) {
                        for (int i3 = max(32, -i2 + 65); i3 <= 34; i3 += 1) {
                          if (ii4 == 0) {
                            for (int i4 = 1; i4 <= -i2 + 64; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                          if (2 * ii4 + 32 >= i3) {
                            for (int i4 = max(max(ii4 + 31, -i2 + 65), 2 * ii4 - 2 * i3 + 95); i4 <= 2 * ii4 + 32; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          } else if (i2 == 34) {
                            B[34][i3][31] = ((((SCALAR_VAL(0.125) * ((A[35][i3][31] - (SCALAR_VAL(2.0) * A[34][i3][31])) + A[33][i3][31])) + (SCALAR_VAL(0.125) * ((A[34][i3 + 1][31] - (SCALAR_VAL(2.0) * A[34][i3][31])) + A[34][i3 - 1][31]))) + (SCALAR_VAL(0.125) * ((A[34][i3][32] - (SCALAR_VAL(2.0) * A[34][i3][31])) + A[34][i3][30]))) + A[34][i3][31]);
                          }
                        }
                      } else {
                        for (int i3 = max(30, -i2 + 63); i3 <= 32; i3 += 1) {
                          if (i3 <= 31) {
                            B[i2][i3][1] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][1] - (SCALAR_VAL(2.0) * A[i2][i3][1])) + A[i2 - 1][i3][1])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][1] - (SCALAR_VAL(2.0) * A[i2][i3][1])) + A[i2][i3 - 1][1]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][2] - (SCALAR_VAL(2.0) * A[i2][i3][1])) + A[i2][i3][0]))) + A[i2][i3][1]);
                          }
                          if (i3 >= 31) {
                            for (int i4 = -i3 + 33; i4 <= 29; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            if (i2 <= 32 && i3 == 32) {
                              for (int i4 = 30; i4 <= -i2 + 63; i4 += 1) {
                                B[i2][32][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32][i4] - (SCALAR_VAL(2.0) * A[i2][32][i4])) + A[i2 - 1][32][i4])) + (SCALAR_VAL(0.125) * ((A[i2][33][i4] - (SCALAR_VAL(2.0) * A[i2][32][i4])) + A[i2][31][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32][i4])) + A[i2][32][i4 - 1]))) + A[i2][32][i4]);
                              }
                            } else if (i3 == 31) {
                              for (int i4 = 30; i4 <= -i2 + 64; i4 += 1) {
                                B[i2][31][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][31][i4] - (SCALAR_VAL(2.0) * A[i2][31][i4])) + A[i2 - 1][31][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32][i4] - (SCALAR_VAL(2.0) * A[i2][31][i4])) + A[i2][30][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][31][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][31][i4])) + A[i2][31][i4 - 1]))) + A[i2][31][i4]);
                              }
                            }
                            if (i2 == 34 && i3 == 31) {
                              B[34][31][31] = ((((SCALAR_VAL(0.125) * ((A[35][31][31] - (SCALAR_VAL(2.0) * A[34][31][31])) + A[33][31][31])) + (SCALAR_VAL(0.125) * ((A[34][32][31] - (SCALAR_VAL(2.0) * A[34][31][31])) + A[34][30][31]))) + (SCALAR_VAL(0.125) * ((A[34][31][32] - (SCALAR_VAL(2.0) * A[34][31][31])) + A[34][31][30]))) + A[34][31][31]);
                            } else if (i2 >= 33 && i3 == 32) {
                              B[i2][32][30] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32][30] - (SCALAR_VAL(2.0) * A[i2][32][30])) + A[i2 - 1][32][30])) + (SCALAR_VAL(0.125) * ((A[i2][33][30] - (SCALAR_VAL(2.0) * A[i2][32][30])) + A[i2][31][30]))) + (SCALAR_VAL(0.125) * ((A[i2][32][31] - (SCALAR_VAL(2.0) * A[i2][32][30])) + A[i2][32][29]))) + A[i2][32][30]);
                            }
                          } else {
                            for (int i4 = 2; i4 <= 32; i4 += 1) {
                              B[i2][30][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][30][i4] - (SCALAR_VAL(2.0) * A[i2][30][i4])) + A[i2 - 1][30][i4])) + (SCALAR_VAL(0.125) * ((A[i2][31][i4] - (SCALAR_VAL(2.0) * A[i2][30][i4])) + A[i2][29][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][30][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][30][i4])) + A[i2][30][i4 - 1]))) + A[i2][30][i4]);
                            }
                          }
                        }
                        for (int i3 = 33; i3 <= 34; i3 += 1) {
                          if (i2 >= 33) {
                            for (int i4 = 1; i4 <= 29; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          } else {
                            for (int i4 = 1; i4 <= -i2 + 62; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                        }
                      }
                    }
                  }
                  for (int i2 = 64 * ii0 - 2 * i0 + 35; i2 <= 32; i2 += 1) {
                    for (int i3 = 64 * ii0 - 2 * i0 - i2 + 68; i3 <= 34; i3 += 1) {
                      for (int i4 = max(max(1, 34 * ii4 - i2 + 30), 34 * ii4 - i2 - i3 + 63); i4 <= min(34, 64 * ii0 + 4 * ii4 - 2 * i0 - i2 - i3 + 100); i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                      if (ii4 == 0 && i3 == 34) {
                        A[i2][34][64 * ii0 - 2 * i0 - i2 + 67] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][34][64 * ii0 - 2 * i0 - i2 + 67] - (SCALAR_VAL(2.0) * B[i2][34][64 * ii0 - 2 * i0 - i2 + 67])) + B[i2 - 1][34][64 * ii0 - 2 * i0 - i2 + 67])) + (SCALAR_VAL(0.125) * ((B[i2][35][64 * ii0 - 2 * i0 - i2 + 67] - (SCALAR_VAL(2.0) * B[i2][34][64 * ii0 - 2 * i0 - i2 + 67])) + B[i2][33][64 * ii0 - 2 * i0 - i2 + 67]))) + (SCALAR_VAL(0.125) * ((B[i2][34][64 * ii0 - 2 * i0 - i2 + 68] - (SCALAR_VAL(2.0) * B[i2][34][64 * ii0 - 2 * i0 - i2 + 67])) + B[i2][34][64 * ii0 - 2 * i0 - i2 + 66]))) + B[i2][34][64 * ii0 - 2 * i0 - i2 + 67]);
                      }
                    }
                  }
                  for (int i2 = 33; i2 <= 34; i2 += 1) {
                    if (ii4 == 1 && i0 == 32 * ii0 + 2) {
                      for (int i3 = 31; i3 <= 32; i3 += 1) {
                        for (int i4 = -i3 + 64; i4 <= 34; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    } else if (ii4 == 0) {
                      for (int i3 = 64 * ii0 - 2 * i0 + 35; i3 <= 32; i3 += 1) {
                        for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 - i3 + 67; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                    for (int i3 = 33; i3 <= 34; i3 += 1) {
                      if (i0 >= 32 * ii0 + ii4 + 1) {
                        for (int i4 = 30 * ii4 + 1; i4 <= 64 * ii0 + 2 * ii4 - 2 * i0 + 34; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                      if (ii4 == 1) {
                        for (int i4 = 33; i4 <= 34; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                  }
                }
              }
              for (int i0 = 32 * ii0 + ii3 - ii4 + 3; i0 <= min(_PB_TSTEPS, 32 * ii0 + 15 * ii4 + 17); i0 += 1) {
                if (ii3 == 1 && ii4 == 1 && i0 >= 32 * ii0 + 18 && 32 * ii0 + 22 >= i0) {
                  for (int i3 = 64 * ii0 - 2 * i0 + 68; i3 <= 34; i3 += 1) {
                    for (int i4 = max(64 * ii0 - 2 * i0 + 68, 64 * ii0 - 2 * i0 - i3 + 101); i4 <= 34; i4 += 1) {
                      B[1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[2][i3][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[0][i3][i4])) + (SCALAR_VAL(0.125) * ((A[1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3][i4 - 1]))) + A[1][i3][i4]);
                    }
                  }
                }
                if (ii3 == 1 && ii4 == 1) {
                  for (int i2 = max(2, 64 * ii0 - 2 * i0 + 36); i2 <= 48 * ii0 - 2 * i0 + (i0 + 1) / 2 + 34; i2 += 1) {
                    for (int i3 = 64 * ii0 - 2 * i0 - i2 + 69; i3 <= 34; i3 += 1) {
                      for (int i4 = max(64 * ii0 - 2 * i0 - i2 + 69, 64 * ii0 - 2 * i0 - i2 - i3 + 102); i4 <= 34; i4 += 1) {
                        B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                      }
                    }
                  }
                  for (int i2 = max(1, 48 * ii0 - 2 * i0 + (i0 + 1) / 2 + 35); i2 <= 32 * ii0 - i0 + 35; i2 += 1) {
                    for (int i3 = 64 * ii0 - 2 * i0 - i2 + 69; i3 <= 34; i3 += 1) {
                      for (int i4 = max(64 * ii0 - 2 * i0 - i2 + 69, 64 * ii0 - 2 * i0 - i2 - i3 + 102); i4 <= 34; i4 += 1) {
                        B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                      }
                    }
                  }
                } else {
                  if (ii3 == 0) {
                    for (int i2 = 64 * ii0 - 2 * i0 + 36; i2 <= 33; i2 += 1) {
                      if (2 * i0 + i2 >= 64 * ii0 + 39 && i2 <= 31) {
                        for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 - i2 + 68; i4 += 1) {
                          B[i2][1][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][1][i4] - (SCALAR_VAL(2.0) * A[i2][1][i4])) + A[i2 - 1][1][i4])) + (SCALAR_VAL(0.125) * ((A[i2][2][i4] - (SCALAR_VAL(2.0) * A[i2][1][i4])) + A[i2][0][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][1][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][1][i4])) + A[i2][1][i4 - 1]))) + A[i2][1][i4]);
                        }
                      } else if (2 * i0 + i2 >= 64 * ii0 + 39 && i2 >= 32) {
                        for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 - i2 + 68; i4 += 1) {
                          B[i2][1][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][1][i4] - (SCALAR_VAL(2.0) * A[i2][1][i4])) + A[i2 - 1][1][i4])) + (SCALAR_VAL(0.125) * ((A[i2][2][i4] - (SCALAR_VAL(2.0) * A[i2][1][i4])) + A[i2][0][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][1][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][1][i4])) + A[i2][1][i4 - 1]))) + A[i2][1][i4]);
                        }
                      } else {
                        for (int i3 = 1; i3 <= 2; i3 += 1) {
                          if (i3 == 2) {
                            B[i2][2][1] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][2][1] - (SCALAR_VAL(2.0) * A[i2][2][1])) + A[i2 - 1][2][1])) + (SCALAR_VAL(0.125) * ((A[i2][3][1] - (SCALAR_VAL(2.0) * A[i2][2][1])) + A[i2][1][1]))) + (SCALAR_VAL(0.125) * ((A[i2][2][2] - (SCALAR_VAL(2.0) * A[i2][2][1])) + A[i2][2][0]))) + A[i2][2][1]);
                          } else {
                            B[i2][1][1] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][1][1] - (SCALAR_VAL(2.0) * A[i2][1][1])) + A[i2 - 1][1][1])) + (SCALAR_VAL(0.125) * ((A[i2][2][1] - (SCALAR_VAL(2.0) * A[i2][1][1])) + A[i2][0][1]))) + (SCALAR_VAL(0.125) * ((A[i2][1][2] - (SCALAR_VAL(2.0) * A[i2][1][1])) + A[i2][1][0]))) + A[i2][1][1]);
                          }
                          for (int i4 = 2; i4 <= 64 * ii0 - 2 * i0 - i2 + 68; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                      if (2 * i0 + i2 >= 64 * ii0 + 38) {
                        for (int i3 = max(2, 64 * ii0 - 2 * i0 - i2 + 41); i3 <= 64 * ii0 - 2 * i0 - i2 + 68; i3 += 1) {
                          for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 - i2 + 68; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      } else {
                        for (int i3 = 3; i3 <= 64 * ii0 - 2 * i0 - i2 + 68; i3 += 1) {
                          for (int i4 = 1; i4 <= 128 * ii0 - 4 * i0 - 2 * i2 + 104; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (2 * i0 + i2 == 64 * ii0 + 37) {
                            B[64 * ii0 - 2 * i0 + 37][i3][31] = ((((SCALAR_VAL(0.125) * ((A[64 * ii0 - 2 * i0 + 38][i3][31] - (SCALAR_VAL(2.0) * A[64 * ii0 - 2 * i0 + 37][i3][31])) + A[64 * ii0 - 2 * i0 + 36][i3][31])) + (SCALAR_VAL(0.125) * ((A[64 * ii0 - 2 * i0 + 37][i3 + 1][31] - (SCALAR_VAL(2.0) * A[64 * ii0 - 2 * i0 + 37][i3][31])) + A[64 * ii0 - 2 * i0 + 37][i3 - 1][31]))) + (SCALAR_VAL(0.125) * ((A[64 * ii0 - 2 * i0 + 37][i3][32] - (SCALAR_VAL(2.0) * A[64 * ii0 - 2 * i0 + 37][i3][31])) + A[64 * ii0 - 2 * i0 + 37][i3][30]))) + A[64 * ii0 - 2 * i0 + 37][i3][31]);
                          }
                        }
                      }
                    }
                  }
                  if (ii3 == 0 && 32 * ii0 + 9 >= i0) {
                    for (int i3 = 1; i3 <= 64 * ii0 - 2 * i0 + 35; i3 += 1) {
                      for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 + 35; i4 += 1) {
                        B[34][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[35][i3][i4] - (SCALAR_VAL(2.0) * A[34][i3][i4])) + A[33][i3][i4])) + (SCALAR_VAL(0.125) * ((A[34][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[34][i3][i4])) + A[34][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[34][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[34][i3][i4])) + A[34][i3][i4 - 1]))) + A[34][i3][i4]);
                      }
                    }
                  }
                }
                if (ii3 == 1) {
                  for (int i2 = max(32 * ii0 + 17 * ii4 - i0 + 19, 64 * ii0 - 2 * i0 + 36); i2 <= 33; i2 += 1) {
                    if (ii4 == 1) {
                      if (32 * ii0 + 18 >= i0 && i2 == 32) {
                        for (int i4 = 33; i4 <= 34; i4 += 1) {
                          B[32][64 * ii0 - 2 * i0 + 37][i4] = ((((SCALAR_VAL(0.125) * ((A[33][64 * ii0 - 2 * i0 + 37][i4] - (SCALAR_VAL(2.0) * A[32][64 * ii0 - 2 * i0 + 37][i4])) + A[31][64 * ii0 - 2 * i0 + 37][i4])) + (SCALAR_VAL(0.125) * ((A[32][64 * ii0 - 2 * i0 + 38][i4] - (SCALAR_VAL(2.0) * A[32][64 * ii0 - 2 * i0 + 37][i4])) + A[32][64 * ii0 - 2 * i0 + 36][i4]))) + (SCALAR_VAL(0.125) * ((A[32][64 * ii0 - 2 * i0 + 37][i4 + 1] - (SCALAR_VAL(2.0) * A[32][64 * ii0 - 2 * i0 + 37][i4])) + A[32][64 * ii0 - 2 * i0 + 37][i4 - 1]))) + A[32][64 * ii0 - 2 * i0 + 37][i4]);
                        }
                      } else if (i2 <= 31) {
                        for (int i3 = max(1, 64 * ii0 - 2 * i0 - i2 + 69); i3 <= 64 * ii0 - 2 * i0 - 2 * i2 + (i2 - 1) / 3 + 91; i3 += 1) {
                          for (int i4 = 64 * ii0 - 2 * i0 - i2 - i3 + 102; i4 <= 34; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                      for (int i3 = max(1, 64 * ii0 - 2 * i0 - 2 * i2 + (i2 - 1) / 3 + 92); i3 <= min(min(-i2 + 64, 64 * ii0 - 2 * i0 - i2 + 103), i2 / 3 + 23); i3 += 1) {
                        for (int i4 = max(1, 64 * ii0 - 2 * i0 - i2 - i3 + 102); i4 <= min(64 * ii0 - 2 * i0 - i2 - i3 + 104, i2 - (i2 + 2) / 3 + 11); i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = 64 * ii0 - 2 * i0 - i2 - i3 + 105; i4 <= min(34, -32 * ii0 + i0 + i2 + (i2 + i3) / 2 - 36); i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = max(i2 - (i2 + 2) / 3 + 12, -32 * ii0 + i0 + i2 + (i2 + i3) / 2 - 35); i4 <= 34; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                      if (i2 >= 32) {
                        for (int i3 = 64 * ii0 - 2 * i0 - i2 + 104; i3 <= -i2 + 64; i3 += 1) {
                          for (int i4 = 1; i4 <= 34; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                        for (int i3 = -i2 + 65; i3 <= 34; i3 += 1) {
                          for (int i4 = max(max(1, 64 * ii0 - 2 * i0 - i2 + 69), 64 * ii0 - 2 * i0 - i2 - i3 + 102); i4 <= min(-i2 + 64, 64 * ii0 - 2 * i0 - i2 - i3 + 104); i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (i0 == 32 * ii0 + 3 && i2 == 33 && i3 == 32) {
                            B[33][32][32] = ((((SCALAR_VAL(0.125) * ((A[34][32][32] - (SCALAR_VAL(2.0) * A[33][32][32])) + A[32][32][32])) + (SCALAR_VAL(0.125) * ((A[33][33][32] - (SCALAR_VAL(2.0) * A[33][32][32])) + A[33][31][32]))) + (SCALAR_VAL(0.125) * ((A[33][32][33] - (SCALAR_VAL(2.0) * A[33][32][32])) + A[33][32][31]))) + A[33][32][32]);
                          }
                          if (i0 == 32 * ii0 + 3 && i2 == 33 && i3 <= 33) {
                            B[33][i3][-i3 + 65] = ((((SCALAR_VAL(0.125) * ((A[34][i3][-i3 + 65] - (SCALAR_VAL(2.0) * A[33][i3][-i3 + 65])) + A[32][i3][-i3 + 65])) + (SCALAR_VAL(0.125) * ((A[33][i3 + 1][-i3 + 65] - (SCALAR_VAL(2.0) * A[33][i3][-i3 + 65])) + A[33][i3 - 1][-i3 + 65]))) + (SCALAR_VAL(0.125) * ((A[33][i3][-i3 + 66] - (SCALAR_VAL(2.0) * A[33][i3][-i3 + 65])) + A[33][i3][-i3 + 64]))) + A[33][i3][-i3 + 65]);
                          }
                          for (int i4 = max(1, 64 * ii0 - 2 * i0 - i2 - i3 + 105); i4 <= min(34, -32 * ii0 + i0 + i2 + (i2 + i3) / 2 - 36); i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          for (int i4 = max(64 * ii0 - 2 * i0 - i2 - i3 + 105, -32 * ii0 + i0 + i2 + (i2 + i3) / 2 - 35); i4 <= 34; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      } else {
                        for (int i3 = 64 * ii0 - 2 * i0 - i2 + 104; i3 <= i2 / 3 + 23; i3 += 1) {
                          for (int i4 = 1; i4 <= 34; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                        for (int i3 = i2 / 3 + 24; i3 <= 34; i3 += 1) {
                          for (int i4 = max(max(1, 64 * ii0 - 2 * i0 - i2 + 69), 64 * ii0 - 2 * i0 - i2 - i3 + 102); i4 <= 34; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      }
                    } else {
                      for (int i3 = 64 * ii0 - 2 * i0 - i2 + 69; i3 <= 34; i3 += 1) {
                        if (i2 + 31 * i3 >= 1025) {
                          for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 - i2 + 68; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (i2 == 33 && i3 == 32) {
                            B[33][32][64 * ii0 - 2 * i0 + 36] = ((((SCALAR_VAL(0.125) * ((A[34][32][64 * ii0 - 2 * i0 + 36] - (SCALAR_VAL(2.0) * A[33][32][64 * ii0 - 2 * i0 + 36])) + A[32][32][64 * ii0 - 2 * i0 + 36])) + (SCALAR_VAL(0.125) * ((A[33][33][64 * ii0 - 2 * i0 + 36] - (SCALAR_VAL(2.0) * A[33][32][64 * ii0 - 2 * i0 + 36])) + A[33][31][64 * ii0 - 2 * i0 + 36]))) + (SCALAR_VAL(0.125) * ((A[33][32][64 * ii0 - 2 * i0 + 37] - (SCALAR_VAL(2.0) * A[33][32][64 * ii0 - 2 * i0 + 36])) + A[33][32][64 * ii0 - 2 * i0 + 35]))) + A[33][32][64 * ii0 - 2 * i0 + 36]);
                          }
                        } else {
                          if (64 * ii0 + 37 >= 2 * i0 + i3) {
                            B[i2][i3][1] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][1] - (SCALAR_VAL(2.0) * A[i2][i3][1])) + A[i2 - 1][i3][1])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][1] - (SCALAR_VAL(2.0) * A[i2][i3][1])) + A[i2][i3 - 1][1]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][2] - (SCALAR_VAL(2.0) * A[i2][i3][1])) + A[i2][i3][0]))) + A[i2][i3][1]);
                            if (i2 == 33 && 2 * i0 + i3 == 64 * ii0 + 36) {
                              for (int i4 = 2; i4 <= 32; i4 += 1) {
                                B[33][64 * ii0 - 2 * i0 + 36][i4] = ((((SCALAR_VAL(0.125) * ((A[34][64 * ii0 - 2 * i0 + 36][i4] - (SCALAR_VAL(2.0) * A[33][64 * ii0 - 2 * i0 + 36][i4])) + A[32][64 * ii0 - 2 * i0 + 36][i4])) + (SCALAR_VAL(0.125) * ((A[33][64 * ii0 - 2 * i0 + 37][i4] - (SCALAR_VAL(2.0) * A[33][64 * ii0 - 2 * i0 + 36][i4])) + A[33][64 * ii0 - 2 * i0 + 35][i4]))) + (SCALAR_VAL(0.125) * ((A[33][64 * ii0 - 2 * i0 + 36][i4 + 1] - (SCALAR_VAL(2.0) * A[33][64 * ii0 - 2 * i0 + 36][i4])) + A[33][64 * ii0 - 2 * i0 + 36][i4 - 1]))) + A[33][64 * ii0 - 2 * i0 + 36][i4]);
                              }
                            }
                          }
                          if (2 * i0 + i3 >= 64 * ii0 + 37) {
                            for (int i4 = max(1, 64 * ii0 - 2 * i0 - i3 + 39); i4 <= 64 * ii0 - 2 * i0 - i2 - i3 + 101; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          }
                        }
                      }
                    }
                  }
                }
                if (6 * ii3 + ii4 + i0 >= 32 * ii0 + 10) {
                  if (ii3 == 0 && ii4 == 0) {
                    for (int i3 = 1; i3 <= 64 * ii0 - 2 * i0 + 34; i3 += 1) {
                      for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 + 35; i4 += 1) {
                        B[34][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[35][i3][i4] - (SCALAR_VAL(2.0) * A[34][i3][i4])) + A[33][i3][i4])) + (SCALAR_VAL(0.125) * ((A[34][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[34][i3][i4])) + A[34][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[34][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[34][i3][i4])) + A[34][i3][i4 - 1]))) + A[34][i3][i4]);
                      }
                    }
                  }
                  for (int i3 = max(1, 64 * ii0 + ii3 - 2 * i0 + 35); i3 <= min(64 * ii0 + 2 * ii3 + 31 * ii4 - 2 * i0 + 35, 64 * ii0 + ii4 - 2 * i0 + 36); i3 += 1) {
                    if (ii3 == 1 && ii4 == 0 && 2 * i0 + i3 == 64 * ii0 + 36) {
                      B[34][64 * ii0 - 2 * i0 + 36][1] = ((((SCALAR_VAL(0.125) * ((A[35][64 * ii0 - 2 * i0 + 36][1] - (SCALAR_VAL(2.0) * A[34][64 * ii0 - 2 * i0 + 36][1])) + A[33][64 * ii0 - 2 * i0 + 36][1])) + (SCALAR_VAL(0.125) * ((A[34][64 * ii0 - 2 * i0 + 37][1] - (SCALAR_VAL(2.0) * A[34][64 * ii0 - 2 * i0 + 36][1])) + A[34][64 * ii0 - 2 * i0 + 35][1]))) + (SCALAR_VAL(0.125) * ((A[34][64 * ii0 - 2 * i0 + 36][2] - (SCALAR_VAL(2.0) * A[34][64 * ii0 - 2 * i0 + 36][1])) + A[34][64 * ii0 - 2 * i0 + 36][0]))) + A[34][64 * ii0 - 2 * i0 + 36][1]);
                    }
                    for (int i4 = 64 * ii0 + 2 * ii3 + 31 * ii4 - 2 * i0 - i3 + 36; i4 <= min(min(ii4 + 32, 64 * ii0 + 34 * ii3 - 2 * i0 + 34), 64 * ii0 - 2 * i0 - i3 + 69); i4 += 1) {
                      if (31 * ii3 + 2 >= i4 || 1) {
                        B[34][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[35][i3][i4] - (SCALAR_VAL(2.0) * A[34][i3][i4])) + A[33][i3][i4])) + (SCALAR_VAL(0.125) * ((A[34][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[34][i3][i4])) + A[34][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[34][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[34][i3][i4])) + A[34][i3][i4 - 1]))) + A[34][i3][i4]);
                      }
                    }
                    if (ii3 == 1 && ii4 == 1) {
                      for (int i4 = 64 * ii0 - 2 * i0 - i3 + 70; i4 <= 34; i4 += 1) {
                        B[34][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[35][i3][i4] - (SCALAR_VAL(2.0) * A[34][i3][i4])) + A[33][i3][i4])) + (SCALAR_VAL(0.125) * ((A[34][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[34][i3][i4])) + A[34][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[34][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[34][i3][i4])) + A[34][i3][i4 - 1]))) + A[34][i3][i4]);
                      }
                    } else if (ii3 == 0) {
                      B[34][64 * ii0 - 2 * i0 + 35][64 * ii0 - 2 * i0 + 35] = ((((SCALAR_VAL(0.125) * ((A[35][64 * ii0 - 2 * i0 + 35][64 * ii0 - 2 * i0 + 35] - (SCALAR_VAL(2.0) * A[34][64 * ii0 - 2 * i0 + 35][64 * ii0 - 2 * i0 + 35])) + A[33][64 * ii0 - 2 * i0 + 35][64 * ii0 - 2 * i0 + 35])) + (SCALAR_VAL(0.125) * ((A[34][64 * ii0 - 2 * i0 + 36][64 * ii0 - 2 * i0 + 35] - (SCALAR_VAL(2.0) * A[34][64 * ii0 - 2 * i0 + 35][64 * ii0 - 2 * i0 + 35])) + A[34][64 * ii0 - 2 * i0 + 34][64 * ii0 - 2 * i0 + 35]))) + (SCALAR_VAL(0.125) * ((A[34][64 * ii0 - 2 * i0 + 35][64 * ii0 - 2 * i0 + 36] - (SCALAR_VAL(2.0) * A[34][64 * ii0 - 2 * i0 + 35][64 * ii0 - 2 * i0 + 35])) + A[34][64 * ii0 - 2 * i0 + 35][64 * ii0 - 2 * i0 + 34]))) + A[34][64 * ii0 - 2 * i0 + 35][64 * ii0 - 2 * i0 + 35]);
                    }
                  }
                  if (ii3 == 1) {
                    for (int i3 = max(1, 64 * ii0 + ii4 - 2 * i0 + 37); i3 <= min(-ii4 + 34, 64 * ii0 - 2 * i0 + 70); i3 += 1) {
                      if (ii4 == 0 && 2 * i0 + i3 == 64 * ii0 + 37) {
                        B[34][64 * ii0 - 2 * i0 + 37][1] = ((((SCALAR_VAL(0.125) * ((A[35][64 * ii0 - 2 * i0 + 37][1] - (SCALAR_VAL(2.0) * A[34][64 * ii0 - 2 * i0 + 37][1])) + A[33][64 * ii0 - 2 * i0 + 37][1])) + (SCALAR_VAL(0.125) * ((A[34][64 * ii0 - 2 * i0 + 38][1] - (SCALAR_VAL(2.0) * A[34][64 * ii0 - 2 * i0 + 37][1])) + A[34][64 * ii0 - 2 * i0 + 36][1]))) + (SCALAR_VAL(0.125) * ((A[34][64 * ii0 - 2 * i0 + 37][2] - (SCALAR_VAL(2.0) * A[34][64 * ii0 - 2 * i0 + 37][1])) + A[34][64 * ii0 - 2 * i0 + 37][0]))) + A[34][64 * ii0 - 2 * i0 + 37][1]);
                      }
                      if (ii4 == 0 && i3 <= 31) {
                        for (int i4 = max(1, 64 * ii0 - 2 * i0 - i3 + 39); i4 <= 64 * ii0 - 2 * i0 + 34; i4 += 1) {
                          B[34][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[35][i3][i4] - (SCALAR_VAL(2.0) * A[34][i3][i4])) + A[33][i3][i4])) + (SCALAR_VAL(0.125) * ((A[34][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[34][i3][i4])) + A[34][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[34][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[34][i3][i4])) + A[34][i3][i4 - 1]))) + A[34][i3][i4]);
                        }
                        if (i3 >= 4) {
                          B[34][i3][64 * ii0 - 2 * i0 + 35] = ((((SCALAR_VAL(0.125) * ((A[35][i3][64 * ii0 - 2 * i0 + 35] - (SCALAR_VAL(2.0) * A[34][i3][64 * ii0 - 2 * i0 + 35])) + A[33][i3][64 * ii0 - 2 * i0 + 35])) + (SCALAR_VAL(0.125) * ((A[34][i3 + 1][64 * ii0 - 2 * i0 + 35] - (SCALAR_VAL(2.0) * A[34][i3][64 * ii0 - 2 * i0 + 35])) + A[34][i3 - 1][64 * ii0 - 2 * i0 + 35]))) + (SCALAR_VAL(0.125) * ((A[34][i3][64 * ii0 - 2 * i0 + 36] - (SCALAR_VAL(2.0) * A[34][i3][64 * ii0 - 2 * i0 + 35])) + A[34][i3][64 * ii0 - 2 * i0 + 34]))) + A[34][i3][64 * ii0 - 2 * i0 + 35]);
                        }
                      } else if (ii4 == 0 && i3 >= 32) {
                        for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 + 34; i4 += 1) {
                          B[34][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[35][i3][i4] - (SCALAR_VAL(2.0) * A[34][i3][i4])) + A[33][i3][i4])) + (SCALAR_VAL(0.125) * ((A[34][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[34][i3][i4])) + A[34][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[34][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[34][i3][i4])) + A[34][i3][i4 - 1]))) + A[34][i3][i4]);
                        }
                        if (i3 == 32) {
                          B[34][32][64 * ii0 - 2 * i0 + 35] = ((((SCALAR_VAL(0.125) * ((A[35][32][64 * ii0 - 2 * i0 + 35] - (SCALAR_VAL(2.0) * A[34][32][64 * ii0 - 2 * i0 + 35])) + A[33][32][64 * ii0 - 2 * i0 + 35])) + (SCALAR_VAL(0.125) * ((A[34][33][64 * ii0 - 2 * i0 + 35] - (SCALAR_VAL(2.0) * A[34][32][64 * ii0 - 2 * i0 + 35])) + A[34][31][64 * ii0 - 2 * i0 + 35]))) + (SCALAR_VAL(0.125) * ((A[34][32][64 * ii0 - 2 * i0 + 36] - (SCALAR_VAL(2.0) * A[34][32][64 * ii0 - 2 * i0 + 35])) + A[34][32][64 * ii0 - 2 * i0 + 34]))) + A[34][32][64 * ii0 - 2 * i0 + 35]);
                        }
                      }
                      if (ii4 == 0) {
                        for (int i4 = 64 * ii0 - 2 * i0 + 36; i4 <= 64 * ii0 - 2 * i0 - i3 + 68; i4 += 1) {
                          B[34][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[35][i3][i4] - (SCALAR_VAL(2.0) * A[34][i3][i4])) + A[33][i3][i4])) + (SCALAR_VAL(0.125) * ((A[34][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[34][i3][i4])) + A[34][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[34][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[34][i3][i4])) + A[34][i3][i4 - 1]))) + A[34][i3][i4]);
                        }
                      }
                      if (32 * ii4 + i3 >= 33 && 64 * ii0 + 68 >= 2 * i0 + i3) {
                        int i4 = i3 >= 33 ? 64 * ii0 + ii4 - 2 * i0 + 35 : 64 * ii0 - 2 * i0 - i3 + 69;
                        B[34][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[35][i3][i4] - (SCALAR_VAL(2.0) * A[34][i3][i4])) + A[33][i3][i4])) + (SCALAR_VAL(0.125) * ((A[34][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[34][i3][i4])) + A[34][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[34][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[34][i3][i4])) + A[34][i3][i4 - 1]))) + A[34][i3][i4]);
                      }
                      if (ii4 == 1) {
                        for (int i4 = max(1, 64 * ii0 - 2 * i0 - i3 + 70); i4 <= 34; i4 += 1) {
                          B[34][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[35][i3][i4] - (SCALAR_VAL(2.0) * A[34][i3][i4])) + A[33][i3][i4])) + (SCALAR_VAL(0.125) * ((A[34][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[34][i3][i4])) + A[34][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[34][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[34][i3][i4])) + A[34][i3][i4 - 1]))) + A[34][i3][i4]);
                        }
                      }
                    }
                    for (int i3 = 64 * ii0 - 2 * i0 + 71; i3 <= 33; i3 += 1) {
                      for (int i4 = 1; i4 <= 34; i4 += 1) {
                        B[34][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[35][i3][i4] - (SCALAR_VAL(2.0) * A[34][i3][i4])) + A[33][i3][i4])) + (SCALAR_VAL(0.125) * ((A[34][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[34][i3][i4])) + A[34][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[34][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[34][i3][i4])) + A[34][i3][i4 - 1]))) + A[34][i3][i4]);
                      }
                    }
                    if (ii4 == 1) {
                      for (int i4 = max(1, 64 * ii0 - 2 * i0 + 36); i4 <= 34; i4 += 1) {
                        B[34][34][i4] = ((((SCALAR_VAL(0.125) * ((A[35][34][i4] - (SCALAR_VAL(2.0) * A[34][34][i4])) + A[33][34][i4])) + (SCALAR_VAL(0.125) * ((A[34][35][i4] - (SCALAR_VAL(2.0) * A[34][34][i4])) + A[34][33][i4]))) + (SCALAR_VAL(0.125) * ((A[34][34][i4 + 1] - (SCALAR_VAL(2.0) * A[34][34][i4])) + A[34][34][i4 - 1]))) + A[34][34][i4]);
                      }
                    }
                  }
                }
                if (ii4 == 0) {
                  for (int i2 = 64 * ii0 - 2 * i0 + 35; i2 <= 32; i2 += 1) {
                    for (int i3 = max(1, 64 * ii0 + 32 * ii3 - 2 * i0 - i2 + 36); i3 <= min(34, 64 * ii0 + 33 * ii3 - 2 * i0 - i2 + 67); i3 += 1) {
                      for (int i4 = 1; i4 <= min(64 * ii0 + 31 * ii3 - 2 * i0 - i2 + 67, 64 * ii0 - 2 * i0 - i2 - i3 + 100); i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                      if (ii3 == 1 && i3 == 34) {
                        A[i2][34][64 * ii0 - 2 * i0 - i2 + 67] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][34][64 * ii0 - 2 * i0 - i2 + 67] - (SCALAR_VAL(2.0) * B[i2][34][64 * ii0 - 2 * i0 - i2 + 67])) + B[i2 - 1][34][64 * ii0 - 2 * i0 - i2 + 67])) + (SCALAR_VAL(0.125) * ((B[i2][35][64 * ii0 - 2 * i0 - i2 + 67] - (SCALAR_VAL(2.0) * B[i2][34][64 * ii0 - 2 * i0 - i2 + 67])) + B[i2][33][64 * ii0 - 2 * i0 - i2 + 67]))) + (SCALAR_VAL(0.125) * ((B[i2][34][64 * ii0 - 2 * i0 - i2 + 68] - (SCALAR_VAL(2.0) * B[i2][34][64 * ii0 - 2 * i0 - i2 + 67])) + B[i2][34][64 * ii0 - 2 * i0 - i2 + 66]))) + B[i2][34][64 * ii0 - 2 * i0 - i2 + 67]);
                      }
                    }
                  }
                  if (ii3 == 1 && i0 == 32 * ii0 + 17) {
                    for (int i2 = 33; i2 <= 34; i2 += 1) {
                      for (int i3 = 1; i3 <= 32; i3 += 1) {
                        for (int i4 = 1; i4 <= -i3 + 33; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                  } else if (32 * ii0 + 16 >= i0) {
                    for (int i2 = 33; i2 <= 34; i2 += 1) {
                      for (int i3 = max(1, 64 * ii0 + 28 * ii3 - 2 * i0 + 7); i3 <= min(32, 64 * ii0 + 30 * ii3 - 2 * i0 + 34); i3 += 1) {
                        for (int i4 = 1; i4 <= min(64 * ii0 + 30 * ii3 - 2 * i0 + 34, 64 * ii0 - 2 * i0 - i3 + 67); i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                      if (ii3 == 1) {
                        for (int i3 = 33; i3 <= 34; i3 += 1) {
                          for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 + 34; i4 += 1) {
                            A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                          }
                        }
                      }
                    }
                  }
                } else {
                  for (int i2 = max(1, 64 * ii0 - 2 * i0 + 35); i2 <= 34; i2 += 1) {
                    for (int i3 = max(max(1, 64 * ii0 - 2 * i0 + 35), 64 * ii0 - 2 * i0 - i2 + 68); i3 <= 34; i3 += 1) {
                      for (int i4 = max(max(max(max(1, 64 * ii0 - 2 * i0 + 35), 64 * ii0 - 2 * i0 - i2 + 68), 64 * ii0 - 2 * i0 - i3 + 68), 64 * ii0 - 2 * i0 - i2 - i3 + 101); i4 <= 34; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                  }
                }
              }
              if (ii3 == 0 && ii4 == 0) {
                for (int i0 = 32 * ii0 + 18; i0 <= min(_PB_TSTEPS, 32 * ii0 + 32); i0 += 1) {
                  for (int i2 = 1; i2 <= 64 * ii0 - 2 * i0 + 67; i2 += 1) {
                    for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 - i2 + 68; i4 += 1) {
                      B[i2][1][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][1][i4] - (SCALAR_VAL(2.0) * A[i2][1][i4])) + A[i2 - 1][1][i4])) + (SCALAR_VAL(0.125) * ((A[i2][2][i4] - (SCALAR_VAL(2.0) * A[i2][1][i4])) + A[i2][0][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][1][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][1][i4])) + A[i2][1][i4 - 1]))) + A[i2][1][i4]);
                    }
                    if (i0 == 32 * ii0 + 18 && i2 <= 2) {
                      for (int i4 = 1; i4 <= -i2 + 32; i4 += 1) {
                        B[i2][2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][2][i4] - (SCALAR_VAL(2.0) * A[i2][2][i4])) + A[i2 - 1][2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][3][i4] - (SCALAR_VAL(2.0) * A[i2][2][i4])) + A[i2][1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][2][i4])) + A[i2][2][i4 - 1]))) + A[i2][2][i4]);
                      }
                    }
                    if (2 * i0 + i2 >= 64 * ii0 + 38) {
                      for (int i3 = max(2, 64 * ii0 - 2 * i0 - i2 + 41); i3 <= 64 * ii0 - 2 * i0 - i2 + 68; i3 += 1) {
                        for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 - i2 + 68; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    } else {
                      for (int i3 = 3; i3 <= 31; i3 += 1) {
                        for (int i4 = 1; i4 <= 31; i4 += 1) {
                          B[1][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[2][i3][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[0][i3][i4])) + (SCALAR_VAL(0.125) * ((A[1][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[1][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[1][i3][i4])) + A[1][i3][i4 - 1]))) + A[1][i3][i4]);
                        }
                      }
                    }
                  }
                  for (int i2 = 1; i2 <= 64 * ii0 - 2 * i0 + 66; i2 += 1) {
                    for (int i3 = 1; i3 <= 64 * ii0 - 2 * i0 - i2 + 67; i3 += 1) {
                      for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 - i2 + 67; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                  }
                }
              } else if (ii4 == 0) {
                for (int i0 = 32 * ii0 + 18; i0 <= min(_PB_TSTEPS, 32 * ii0 + 32); i0 += 1) {
                  for (int i2 = 1; i2 <= 34; i2 += 1) {
                    for (int i3 = max(1, 64 * ii0 - 2 * i0 - i2 + 69); i3 <= min(32, 64 * ii0 - 2 * i0 - i2 + 100); i3 += 1) {
                      if (i0 == 32 * ii0 + 18 && i3 == 1) {
                        B[i2][1][1] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][1][1] - (SCALAR_VAL(2.0) * A[i2][1][1])) + A[i2 - 1][1][1])) + (SCALAR_VAL(0.125) * ((A[i2][2][1] - (SCALAR_VAL(2.0) * A[i2][1][1])) + A[i2][0][1]))) + (SCALAR_VAL(0.125) * ((A[i2][1][2] - (SCALAR_VAL(2.0) * A[i2][1][1])) + A[i2][1][0]))) + A[i2][1][1]);
                      }
                      for (int i4 = max(1, 64 * ii0 - 2 * i0 - i3 + 39); i4 <= 64 * ii0 - 2 * i0 - i2 - i3 + 101; i4 += 1) {
                        B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                      }
                      if (i2 == 34) {
                        B[34][i3][64 * ii0 - 2 * i0 - i3 + 68] = ((((SCALAR_VAL(0.125) * ((A[35][i3][64 * ii0 - 2 * i0 - i3 + 68] - (SCALAR_VAL(2.0) * A[34][i3][64 * ii0 - 2 * i0 - i3 + 68])) + A[33][i3][64 * ii0 - 2 * i0 - i3 + 68])) + (SCALAR_VAL(0.125) * ((A[34][i3 + 1][64 * ii0 - 2 * i0 - i3 + 68] - (SCALAR_VAL(2.0) * A[34][i3][64 * ii0 - 2 * i0 - i3 + 68])) + A[34][i3 - 1][64 * ii0 - 2 * i0 - i3 + 68]))) + (SCALAR_VAL(0.125) * ((A[34][i3][64 * ii0 - 2 * i0 - i3 + 69] - (SCALAR_VAL(2.0) * A[34][i3][64 * ii0 - 2 * i0 - i3 + 68])) + A[34][i3][64 * ii0 - 2 * i0 - i3 + 67]))) + A[34][i3][64 * ii0 - 2 * i0 - i3 + 68]);
                      }
                    }
                    if (i2 == 34) {
                      B[34][64 * ii0 - 2 * i0 + 67][1] = ((((SCALAR_VAL(0.125) * ((A[35][64 * ii0 - 2 * i0 + 67][1] - (SCALAR_VAL(2.0) * A[34][64 * ii0 - 2 * i0 + 67][1])) + A[33][64 * ii0 - 2 * i0 + 67][1])) + (SCALAR_VAL(0.125) * ((A[34][64 * ii0 - 2 * i0 + 68][1] - (SCALAR_VAL(2.0) * A[34][64 * ii0 - 2 * i0 + 67][1])) + A[34][64 * ii0 - 2 * i0 + 66][1]))) + (SCALAR_VAL(0.125) * ((A[34][64 * ii0 - 2 * i0 + 67][2] - (SCALAR_VAL(2.0) * A[34][64 * ii0 - 2 * i0 + 67][1])) + A[34][64 * ii0 - 2 * i0 + 67][0]))) + A[34][64 * ii0 - 2 * i0 + 67][1]);
                    }
                    for (int i3 = 33; i3 <= 34; i3 += 1) {
                      for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 - i2 + 68; i4 += 1) {
                        B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                      }
                    }
                  }
                  for (int i2 = 1; i2 <= 34; i2 += 1) {
                    if (i2 >= 33) {
                      for (int i3 = 1; i3 <= 64 * ii0 - 2 * i0 + 66; i3 += 1) {
                        for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 - i3 + 67; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    } else {
                      for (int i3 = max(1, 64 * ii0 - 2 * i0 - i2 + 68); i3 <= min(33, 64 * ii0 - 2 * i0 - i2 + 99); i3 += 1) {
                        for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 - i2 - i3 + 100; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                      for (int i4 = 1; i4 <= 64 * ii0 - 2 * i0 - i2 + 67; i4 += 1) {
                        A[i2][34][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][34][i4] - (SCALAR_VAL(2.0) * B[i2][34][i4])) + B[i2 - 1][34][i4])) + (SCALAR_VAL(0.125) * ((B[i2][35][i4] - (SCALAR_VAL(2.0) * B[i2][34][i4])) + B[i2][33][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][34][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][34][i4])) + B[i2][34][i4 - 1]))) + B[i2][34][i4]);
                      }
                    }
                  }
                }
              }
            } else {
              for (int i0 = 32 * ii0 + 1; i0 <= min(_PB_TSTEPS, 32 * ii0 + 32); i0 += 1) {
                if (i0 >= 32 * ii0 + 2) {
                  if (i0 >= 32 * ii0 + 3) {
                    for (int i2 = max(1, 64 * ii0 - 2 * i0 + 36); i2 <= min(34, 64 * ii0 - 2 * i0 + 67); i2 += 1) {
                      for (int i3 = 1; i3 < i2 - 31; i3 += 1) {
                        for (int i4 = 64 * ii0 - 2 * i0 + 36; i4 <= 64 * ii0 - 2 * i0 + i2 - i3 + 4; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                        for (int i4 = 64 * ii0 - 2 * i0 + i2 - i3 + 5; i4 <= 34; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                      if (i2 <= 33) {
                        for (int i3 = max(1, i2 - 31); i3 <= 64 * ii0 - 2 * i0 - i2 + 68; i3 += 1) {
                          if (i3 <= 2) {
                            for (int i4 = 64 * ii0 - 2 * i0 - i2 + 69; i4 <= 34; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          } else if (i3 <= 31) {
                            if (2 * i0 + i2 >= 64 * ii0 + 38) {
                              for (int i4 = 64 * ii0 - 2 * i0 - i2 + 69; i4 <= min(33, -i2 + 64); i4 += 1) {
                                B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                              }
                            }
                            for (int i4 = -i2 + 65; i4 <= 33; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            for (int i4 = max(-128 * ii0 + 4 * i0 + 2 * i2 - 42, 128 * ii0 - 4 * i0 - 2 * i2 + 105); i4 <= 33; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            B[i2][i3][34] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][34] - (SCALAR_VAL(2.0) * A[i2][i3][34])) + A[i2 - 1][i3][34])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][34] - (SCALAR_VAL(2.0) * A[i2][i3][34])) + A[i2][i3 - 1][34]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][35] - (SCALAR_VAL(2.0) * A[i2][i3][34])) + A[i2][i3][33]))) + A[i2][i3][34]);
                          } else {
                            for (int i4 = 33; i4 <= 34; i4 += 1) {
                              B[64 * ii0 - 2 * i0 + 36][32][i4] = ((((SCALAR_VAL(0.125) * ((A[64 * ii0 - 2 * i0 + 37][32][i4] - (SCALAR_VAL(2.0) * A[64 * ii0 - 2 * i0 + 36][32][i4])) + A[64 * ii0 - 2 * i0 + 35][32][i4])) + (SCALAR_VAL(0.125) * ((A[64 * ii0 - 2 * i0 + 36][33][i4] - (SCALAR_VAL(2.0) * A[64 * ii0 - 2 * i0 + 36][32][i4])) + A[64 * ii0 - 2 * i0 + 36][31][i4]))) + (SCALAR_VAL(0.125) * ((A[64 * ii0 - 2 * i0 + 36][32][i4 + 1] - (SCALAR_VAL(2.0) * A[64 * ii0 - 2 * i0 + 36][32][i4])) + A[64 * ii0 - 2 * i0 + 36][32][i4 - 1]))) + A[64 * ii0 - 2 * i0 + 36][32][i4]);
                            }
                          }
                        }
                      } else {
                        for (int i3 = 3; i3 <= 64 * ii0 - 2 * i0 + 35; i3 += 1) {
                          for (int i4 = 64 * ii0 - 2 * i0 + 36; i4 <= 34; i4 += 1) {
                            B[34][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[35][i3][i4] - (SCALAR_VAL(2.0) * A[34][i3][i4])) + A[33][i3][i4])) + (SCALAR_VAL(0.125) * ((A[34][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[34][i3][i4])) + A[34][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[34][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[34][i3][i4])) + A[34][i3][i4 - 1]))) + A[34][i3][i4]);
                          }
                        }
                      }
                    }
                    if (i0 == 32 * ii0 + 17) {
                      for (int i4 = 2; i4 <= 34; i4 += 1) {
                        B[34][1][i4] = ((((SCALAR_VAL(0.125) * ((A[35][1][i4] - (SCALAR_VAL(2.0) * A[34][1][i4])) + A[33][1][i4])) + (SCALAR_VAL(0.125) * ((A[34][2][i4] - (SCALAR_VAL(2.0) * A[34][1][i4])) + A[34][0][i4]))) + (SCALAR_VAL(0.125) * ((A[34][1][i4 + 1] - (SCALAR_VAL(2.0) * A[34][1][i4])) + A[34][1][i4 - 1]))) + A[34][1][i4]);
                      }
                    }
                  } else {
                    for (int i2 = 32; i2 <= 34; i2 += 1) {
                      for (int i3 = 1; i3 <= -i2 + 64; i3 += 1) {
                        for (int i4 = max(32, -i2 + 65); i4 <= 34; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                      if (i2 == 34) {
                        for (int i4 = 32; i4 <= 34; i4 += 1) {
                          B[34][31][i4] = ((((SCALAR_VAL(0.125) * ((A[35][31][i4] - (SCALAR_VAL(2.0) * A[34][31][i4])) + A[33][31][i4])) + (SCALAR_VAL(0.125) * ((A[34][32][i4] - (SCALAR_VAL(2.0) * A[34][31][i4])) + A[34][30][i4]))) + (SCALAR_VAL(0.125) * ((A[34][31][i4 + 1] - (SCALAR_VAL(2.0) * A[34][31][i4])) + A[34][31][i4 - 1]))) + A[34][31][i4]);
                        }
                      }
                    }
                  }
                }
                for (int i2 = max(1, 64 * ii0 - 2 * i0 + 35); i2 <= min(32, 64 * ii0 - 2 * i0 + 66); i2 += 1) {
                  for (int i3 = 1; i3 <= 64 * ii0 - 2 * i0 - i2 + 67; i3 += 1) {
                    for (int i4 = 64 * ii0 - 2 * i0 - i2 + 68; i4 <= 34; i4 += 1) {
                      A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                    }
                  }
                }
                for (int i2 = 33; i2 <= 34; i2 += 1) {
                  for (int i3 = 1; i3 <= 64 * ii0 - 2 * i0 + 34; i3 += 1) {
                    for (int i4 = 64 * ii0 - 2 * i0 + 35; i4 <= 34; i4 += 1) {
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
      for (int ii2 = 0; ii2 <= (_PB_N - 3) / 32; ii2 += 1) {
        for (int ii3 = 0; ii3 <= (_PB_N - 3) / 32; ii3 += 1) {
          if (_PB_N >= 32 * ii3 + 35) {
            for (int ii4 = 0; ii4 <= (_PB_N - 3) / 32; ii4 += 1) {
              for (int i0 = _PB_TSTEPS - 1; i0 <= _PB_TSTEPS; i0 += 1) {
                if (_PB_N >= 32 * ii2 + 35 && i0 == _PB_TSTEPS) {
                  for (int i2 = max(1, 32 * ii2); i2 <= 32 * ii2 + 31; i2 += 1) {
                    if (_PB_N >= 32 * ii4 + 35) {
                      if (ii3 >= 1 && i2 >= 32 * ii2 + 1 && 32 * ii2 + 2 >= i2) {
                        for (int i4 = 32 * ii4 + 1; i4 <= 32 * ii4 + 32; i4 += 1) {
                          B[i2][32 * ii3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii3][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][i4])) + A[i2 - 1][32 * ii3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][i4])) + A[i2][32 * ii3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][i4])) + A[i2][32 * ii3][i4 - 1]))) + A[i2][32 * ii3][i4]);
                        }
                      }
                      if (32 * ii2 + 1 >= i2) {
                        for (int i3 = 32 * ii3 + 1; i3 <= 32 * ii2 + 32 * ii3 - i2 + 32; i3 += 1) {
                          for (int i4 = max(1, 32 * ii2 + 32 * ii4 - i2 + 1); i4 <= 32 * ii2 + 32 * ii4 - i2 + 32; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        }
                      } else {
                        if (ii3 >= 1 && i2 >= 32 * ii2 + 3) {
                          for (int i4 = 32 * ii4 + 1; i4 <= 32 * ii4 + 32; i4 += 1) {
                            B[i2][32 * ii3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii3][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][i4])) + A[i2 - 1][32 * ii3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][i4])) + A[i2][32 * ii3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][i4])) + A[i2][32 * ii3][i4 - 1]))) + A[i2][32 * ii3][i4]);
                          }
                        }
                        for (int i3 = 32 * ii3 + 1; i3 <= 32 * ii3 + 31; i3 += 1) {
                          if (i2 >= 32 * ii2 + 3 && i2 + i3 >= 32 * ii2 + 32 * ii3 + 33) {
                            for (int i4 = max(1, 32 * ii4); i4 <= 32 * ii4 + 31; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          } else if (i2 >= 32 * ii2 + 3 && 32 * ii2 + 32 * ii3 + 32 >= i2 + i3) {
                            for (int i4 = max(1, 32 * ii4); i4 <= min(32 * ii4 + 31, 32 * ii3 + 32 * ii4 - i3 + 33); i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                            for (int i4 = 32 * ii3 + 32 * ii4 - i3 + 34; i4 <= 32 * ii4 + 31; i4 += 1) {
                              B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                            }
                          } else {
                            for (int i4 = max(1, 32 * ii4); i4 <= 32 * ii4 + 31; i4 += 1) {
                              B[32 * ii2 + 2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 3][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][i3][i4])) + A[32 * ii2 + 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][i3][i4])) + A[32 * ii2 + 2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][i3][i4])) + A[32 * ii2 + 2][i3][i4 - 1]))) + A[32 * ii2 + 2][i3][i4]);
                            }
                          }
                        }
                      }
                    } else {
                      for (int i3 = max(max(1, 32 * ii3), 32 * ii2 + 32 * ii3 - i2 + 1); i3 <= 32 * ii2 + 32 * ii3 - i2 + 32; i3 += 1) {
                        for (int i4 = max(max(32 * ii4, 64 * ii2 + 32 * ii4 - 2 * i2 + 1), 16 * ii3 + 32 * ii4 - i3 + i3 / 2 + 1); i4 < _PB_N - 1; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                      for (int i3 = 32 * ii2 + 32 * ii3 - i2 + 33; i3 <= 32 * ii3 + 31; i3 += 1) {
                        for (int i4 = 32 * ii4; i4 < _PB_N - 1; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    }
                  }
                } else if (32 * ii2 + 34 >= _PB_N && _PB_N >= 32 * ii4 + 35 && i0 == _PB_TSTEPS) {
                  for (int i2 = 32 * ii2; i2 < _PB_N - 1; i2 += 1) {
                    for (int i3 = max(max(1, 32 * ii3), 32 * ii2 + 32 * ii3 - i2 + 1); i3 <= 32 * ii2 + 32 * ii3 - i2 + 2; i3 += 1) {
                      for (int i4 = max(max(1, 32 * ii2 + 32 * ii4 - i2 + 1), 64 * ii3 + 32 * ii4 - 2 * i3 + 1); i4 <= 32 * ii2 + 32 * ii4 - i2 + 32; i4 += 1) {
                        B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                      }
                      if (i3 == 32 * ii3) {
                        for (int i4 = 32 * ii2 + 32 * ii4 - i2 + 33; i4 <= 32 * ii4 + 32; i4 += 1) {
                          B[i2][32 * ii3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii3][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][i4])) + A[i2 - 1][32 * ii3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][i4])) + A[i2][32 * ii3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][i4])) + A[i2][32 * ii3][i4 - 1]))) + A[i2][32 * ii3][i4]);
                        }
                      }
                    }
                    if (32 * ii2 + 1 >= i2) {
                      for (int i3 = 32 * ii2 + 32 * ii3 - i2 + 3; i3 <= 32 * ii2 + 32 * ii3 - i2 + 32; i3 += 1) {
                        for (int i4 = max(1, 32 * ii2 + 32 * ii4 - i2 + 1); i4 <= 32 * ii2 + 32 * ii4 - i2 + 32; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    } else {
                      for (int i3 = max(max(1, 32 * ii3), 32 * ii2 + 32 * ii3 - i2 + 3); i3 <= 32 * ii3 + 31; i3 += 1) {
                        if (i2 >= 32 * ii2 + 3) {
                          for (int i4 = max(max(1, 32 * ii4), 64 * ii3 + 32 * ii4 - 2 * i3 + 1); i4 <= 32 * ii4 + 31; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (i3 == 32 * ii3) {
                            B[i2][32 * ii3][32 * ii4 + 32] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii3][32 * ii4 + 32] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][32 * ii4 + 32])) + A[i2 - 1][32 * ii3][32 * ii4 + 32])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 1][32 * ii4 + 32] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][32 * ii4 + 32])) + A[i2][32 * ii3 - 1][32 * ii4 + 32]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3][32 * ii4 + 33] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][32 * ii4 + 32])) + A[i2][32 * ii3][32 * ii4 + 31]))) + A[i2][32 * ii3][32 * ii4 + 32]);
                          }
                        } else {
                          for (int i4 = max(1, 32 * ii4); i4 <= 32 * ii4 + 31; i4 += 1) {
                            B[32 * ii2 + 2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 3][i3][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][i3][i4])) + A[32 * ii2 + 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][i3][i4])) + A[32 * ii2 + 2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][i3][i4])) + A[32 * ii2 + 2][i3][i4 - 1]))) + A[32 * ii2 + 2][i3][i4]);
                          }
                        }
                      }
                    }
                  }
                }
                if (_PB_N >= 32 * ii4 + 35) {
                  if (i0 == _PB_TSTEPS) {
                    for (int i2 = max(1, 32 * ii2 - 1); i2 <= 32 * ii2; i2 += 1) {
                      for (int i3 = max(1, 32 * ii2 + 32 * ii3 - i2); i3 <= 32 * ii2 + 32 * ii3 - i2 + 31; i3 += 1) {
                        for (int i4 = max(max(1, 32 * ii2 + 32 * ii4 - i2), 32 * ii2 + 32 * ii3 + 32 * ii4 - i2 - i3 + 1); i4 <= 32 * ii2 + 32 * ii3 + 32 * ii4 - i2 - i3 + 32; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                        for (int i4 = 32 * ii2 + 32 * ii3 + 32 * ii4 - i2 - i3 + 33; i4 <= 32 * ii2 + 32 * ii4 - i2 + 31; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                  }
                  if (_PB_N >= 32 * ii2 + 35) {
                    for (int i2 = 32 * ii2 + 1; i2 <= 2 * _PB_TSTEPS + 32 * ii2 - 2 * i0 + 30; i2 += 1) {
                      if (i0 == _PB_TSTEPS) {
                        for (int i3 = max(1, 32 * ii3 - 1); i3 <= 32 * ii3; i3 += 1) {
                          for (int i4 = max(1, 32 * ii3 + 32 * ii4 - i3); i4 <= 32 * ii3 + 32 * ii4 - i3 + 31; i4 += 1) {
                            A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                          }
                        }
                      }
                      for (int i3 = 32 * ii3 + 1; i3 <= 2 * _PB_TSTEPS + 32 * ii3 - 2 * i0 + 30; i3 += 1) {
                        for (int i4 = max(1, 2 * _PB_TSTEPS + 32 * ii4 - 2 * i0 - 1); i4 <= 2 * _PB_TSTEPS + 32 * ii4 - 2 * i0 + 30; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                  } else {
                    for (int i2 = 32 * ii2 + 1; i2 < _PB_N - 1; i2 += 1) {
                      if (i0 == _PB_TSTEPS) {
                        for (int i3 = max(1, 32 * ii3 - 1); i3 <= 32 * ii3; i3 += 1) {
                          for (int i4 = max(1, 32 * ii3 + 32 * ii4 - i3); i4 <= 32 * ii3 + 32 * ii4 - i3 + 31; i4 += 1) {
                            A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                          }
                        }
                      }
                      for (int i3 = 32 * ii3 + 1; i3 <= 2 * _PB_TSTEPS + 32 * ii3 - 2 * i0 + 30; i3 += 1) {
                        if (i0 == _PB_TSTEPS) {
                          for (int i4 = max(1, 32 * ii4 - 1); i4 <= 32 * ii4; i4 += 1) {
                            A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                          }
                        }
                        for (int i4 = 32 * ii4 + 1; i4 <= 2 * _PB_TSTEPS + 32 * ii4 - 2 * i0 + 30; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                  }
                } else if (_PB_N >= 32 * ii2 + 35) {
                  if (i0 == _PB_TSTEPS) {
                    for (int i2 = max(1, 32 * ii2 - 1); i2 <= 32 * ii2; i2 += 1) {
                      for (int i3 = max(1, 32 * ii2 + 32 * ii3 - i2); i3 <= 32 * ii2 + 32 * ii3 - i2 + 31; i3 += 1) {
                        for (int i4 = max(32 * ii2 + 32 * ii4 - i2, 32 * ii2 + 32 * ii3 + 32 * ii4 - i2 - i3 + 1); i4 < _PB_N - 1; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                  }
                  for (int i2 = 32 * ii2 + 1; i2 <= 2 * _PB_TSTEPS + 32 * ii2 - 2 * i0 + 30; i2 += 1) {
                    if (i0 == _PB_TSTEPS) {
                      for (int i3 = max(1, 32 * ii3 - 1); i3 <= 32 * ii3; i3 += 1) {
                        for (int i4 = 32 * ii3 + 32 * ii4 - i3; i4 < _PB_N - 1; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                    for (int i3 = 32 * ii3 + 1; i3 <= 2 * _PB_TSTEPS + 32 * ii3 - 2 * i0 + 30; i3 += 1) {
                      for (int i4 = 2 * _PB_TSTEPS + 32 * ii4 - 2 * i0 - 1; i4 < _PB_N - 1; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                  }
                } else {
                  if (i0 == _PB_TSTEPS) {
                    for (int i2 = 32 * ii2; i2 < _PB_N - 1; i2 += 1) {
                      for (int i3 = max(max(1, 32 * ii3), 32 * ii2 + 32 * ii3 - i2 + 1); i3 <= 32 * ii2 + 32 * ii3 - i2 + 32; i3 += 1) {
                        for (int i4 = max(max(32 * ii2, 64 * ii2 - i2 + 1), 32 * ii2 + 16 * ii3 - i3 + i3 / 2 + 1); i4 < _PB_N - 1; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                      for (int i3 = 32 * ii2 + 32 * ii3 - i2 + 33; i3 <= 32 * ii3 + 31; i3 += 1) {
                        for (int i4 = 32 * ii2; i4 < _PB_N - 1; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    }
                  }
                  if (i0 == _PB_TSTEPS) {
                    for (int i2 = 32 * ii2 - 1; i2 <= 32 * ii2; i2 += 1) {
                      for (int i3 = max(1, 32 * ii2 + 32 * ii3 - i2); i3 <= 32 * ii2 + 32 * ii3 - i2 + 31; i3 += 1) {
                        for (int i4 = max(64 * ii2 - i2, 64 * ii2 + 32 * ii3 - i2 - i3 + 1); i4 < _PB_N - 1; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                  }
                  for (int i2 = 32 * ii2 + 1; i2 < _PB_N - 1; i2 += 1) {
                    if (i0 == _PB_TSTEPS) {
                      for (int i3 = max(1, 32 * ii3 - 1); i3 <= 32 * ii3; i3 += 1) {
                        for (int i4 = 32 * ii2 + 32 * ii3 - i3; i4 < _PB_N - 1; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                    for (int i3 = 32 * ii3 + 1; i3 <= 2 * _PB_TSTEPS + 32 * ii3 - 2 * i0 + 30; i3 += 1) {
                      if (i0 == _PB_TSTEPS) {
                        for (int i4 = 32 * ii2 - 1; i4 <= 32 * ii2; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                      for (int i4 = 32 * ii2 + 1; i4 < _PB_N - 1; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                  }
                }
              }
            }
          } else if (_PB_N >= 32 * ii2 + 35) {
            for (int ii4 = 0; ii4 <= ii3; ii4 += 1) {
              for (int i0 = _PB_TSTEPS - 1; i0 <= _PB_TSTEPS; i0 += 1) {
                if (i0 == _PB_TSTEPS) {
                  if (_PB_N >= 32 * ii4 + 35) {
                    for (int i2 = max(1, 32 * ii2); i2 <= 32 * ii2 + 31; i2 += 1) {
                      if (i2 >= 32 * ii2 + 3) {
                        for (int i3 = 32 * ii3; i3 < _PB_N - 1; i3 += 1) {
                          for (int i4 = max(max(1, 32 * ii4), 16 * ii3 + 32 * ii4 - i3 + i3 / 2 + 1); i4 <= 32 * ii4 + 31; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (i3 == 32 * ii3) {
                            B[i2][32 * ii3][32 * ii4 + 32] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii3][32 * ii4 + 32] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][32 * ii4 + 32])) + A[i2 - 1][32 * ii3][32 * ii4 + 32])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 1][32 * ii4 + 32] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][32 * ii4 + 32])) + A[i2][32 * ii3 - 1][32 * ii4 + 32]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3][32 * ii4 + 33] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][32 * ii4 + 32])) + A[i2][32 * ii3][32 * ii4 + 31]))) + A[i2][32 * ii3][32 * ii4 + 32]);
                          }
                        }
                      } else {
                        for (int i3 = max(32 * ii3, 32 * ii2 + 32 * ii3 - i2 + 1); i3 < _PB_N - 1; i3 += 1) {
                          for (int i4 = max(max(max(1, 32 * ii4), 32 * ii2 + 32 * ii4 - i2 + 1), 16 * ii3 + 32 * ii4 - i3 + i3 / 2 + 1); i4 <= 32 * ii2 + 32 * ii4 - i2 + 32; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (i3 == 32 * ii3) {
                            for (int i4 = 32 * ii2 + 32 * ii4 - i2 + 33; i4 <= 32 * ii4 + 32; i4 += 1) {
                              B[i2][32 * ii3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii3][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][i4])) + A[i2 - 1][32 * ii3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][i4])) + A[i2][32 * ii3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii3][i4])) + A[i2][32 * ii3][i4 - 1]))) + A[i2][32 * ii3][i4]);
                            }
                          } else if (i2 == 32 * ii2 + 2) {
                            B[32 * ii2 + 2][i3][32 * ii4 + 31] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 3][i3][32 * ii4 + 31] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][i3][32 * ii4 + 31])) + A[32 * ii2 + 1][i3][32 * ii4 + 31])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3 + 1][32 * ii4 + 31] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][i3][32 * ii4 + 31])) + A[32 * ii2 + 2][i3 - 1][32 * ii4 + 31]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3][32 * ii4 + 32] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][i3][32 * ii4 + 31])) + A[32 * ii2 + 2][i3][32 * ii4 + 30]))) + A[32 * ii2 + 2][i3][32 * ii4 + 31]);
                          }
                        }
                      }
                    }
                  } else {
                    for (int i2 = max(1, 32 * ii2); i2 <= 32 * ii2 + 31; i2 += 1) {
                      for (int i3 = max(32 * ii3, 32 * ii2 + 32 * ii3 - i2 + 1); i3 < _PB_N - 1; i3 += 1) {
                        for (int i4 = max(max(32 * ii3, 32 * ii2 + 32 * ii3 - i2 + 1), 96 * ii3 - 2 * i3 + 1); i4 < _PB_N - 1; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    }
                  }
                }
                if (_PB_N >= 32 * ii4 + 35 && i0 == _PB_TSTEPS) {
                  for (int i2 = max(1, 32 * ii2 - 1); i2 <= 32 * ii2; i2 += 1) {
                    for (int i3 = 32 * ii2 + 32 * ii3 - i2; i3 <= 32 * ii3 + 1; i3 += 1) {
                      for (int i4 = max(1, 32 * ii2 + 32 * ii3 + 32 * ii4 - i2 - i3 + 1); i4 <= 32 * ii2 + 32 * ii3 + 32 * ii4 - i2 - i3 + 32; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                    for (int i3 = 32 * ii3 + 2; i3 < _PB_N - 1; i3 += 1) {
                      for (int i4 = max(1, 32 * ii2 + 32 * ii4 - i2); i4 <= 32 * ii2 + 32 * ii4 - i2 + 31; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                  }
                } else if (ii4 == ii3 && i0 == _PB_TSTEPS) {
                  for (int i2 = max(1, 32 * ii2 - 1); i2 <= 32 * ii2; i2 += 1) {
                    for (int i3 = 32 * ii2 + 32 * ii3 - i2; i3 < _PB_N - 1; i3 += 1) {
                      for (int i4 = max(32 * ii2 + 32 * ii3 - i2, 32 * ii2 + 64 * ii3 - i2 - i3 + 1); i4 < _PB_N - 1; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                  }
                }
                for (int i2 = 32 * ii2 + 1; i2 <= 2 * _PB_TSTEPS + 32 * ii2 - 2 * i0 + 30; i2 += 1) {
                  if (_PB_N >= 32 * ii4 + 35 && i0 == _PB_TSTEPS) {
                    for (int i3 = 32 * ii3 - 1; i3 <= 32 * ii3; i3 += 1) {
                      for (int i4 = max(1, 32 * ii3 + 32 * ii4 - i3); i4 <= 32 * ii3 + 32 * ii4 - i3 + 31; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                  } else if (ii4 == ii3 && i0 == _PB_TSTEPS) {
                    for (int i3 = 32 * ii3 - 1; i3 <= 32 * ii3; i3 += 1) {
                      for (int i4 = 64 * ii3 - i3; i4 < _PB_N - 1; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                  }
                  for (int i3 = 32 * ii3 + 1; i3 < _PB_N - 1; i3 += 1) {
                    if (_PB_N >= 32 * ii4 + 35) {
                      for (int i4 = max(1, 2 * _PB_TSTEPS + 32 * ii4 - 2 * i0 - 1); i4 <= 2 * _PB_TSTEPS + 32 * ii4 - 2 * i0 + 30; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    } else {
                      if (i0 == _PB_TSTEPS) {
                        for (int i4 = 32 * ii3 - 1; i4 <= 32 * ii3; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                      for (int i4 = 32 * ii3 + 1; i4 < _PB_N - 1; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                  }
                }
              }
            }
          } else {
            for (int ii4 = 0; ii4 <= ii2; ii4 += 1) {
              if (_PB_N >= 32 * ii4 + 35) {
                for (int i0 = _PB_TSTEPS - 1; i0 <= _PB_TSTEPS; i0 += 1) {
                  if (i0 == _PB_TSTEPS) {
                    for (int i2 = 32 * ii2; i2 < _PB_N - 1; i2 += 1) {
                      for (int i3 = max(32 * ii2, 64 * ii2 - i2 + 1); i3 < _PB_N - 1; i3 += 1) {
                        if (i2 >= 32 * ii2 + 3) {
                          for (int i4 = max(max(1, 32 * ii4), 64 * ii2 + 32 * ii4 - 2 * i3 + 1); i4 <= 64 * ii2 + 32 * ii4 - 2 * i3 + 32; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          for (int i4 = max(max(1, 32 * ii4), 64 * ii2 + 32 * ii4 - 2 * i3 + 33); i4 <= 32 * ii4 + 31; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                        } else {
                          for (int i4 = max(max(max(1, 32 * ii4), 32 * ii2 + 32 * ii4 - i2 + 1), 64 * ii2 + 32 * ii4 - 2 * i3 + 1); i4 <= 32 * ii2 + 32 * ii4 - i2 + 32; i4 += 1) {
                            B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                          }
                          if (i2 == 32 * ii2 + 2 && i3 >= 32 * ii2 + 1) {
                            B[32 * ii2 + 2][i3][32 * ii4 + 31] = ((((SCALAR_VAL(0.125) * ((A[32 * ii2 + 3][i3][32 * ii4 + 31] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][i3][32 * ii4 + 31])) + A[32 * ii2 + 1][i3][32 * ii4 + 31])) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3 + 1][32 * ii4 + 31] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][i3][32 * ii4 + 31])) + A[32 * ii2 + 2][i3 - 1][32 * ii4 + 31]))) + (SCALAR_VAL(0.125) * ((A[32 * ii2 + 2][i3][32 * ii4 + 32] - (SCALAR_VAL(2.0) * A[32 * ii2 + 2][i3][32 * ii4 + 31])) + A[32 * ii2 + 2][i3][32 * ii4 + 30]))) + A[32 * ii2 + 2][i3][32 * ii4 + 31]);
                          } else if (i3 == 32 * ii2) {
                            for (int i4 = 32 * ii2 + 32 * ii4 - i2 + 33; i4 <= 32 * ii4 + 32; i4 += 1) {
                              B[i2][32 * ii2][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][32 * ii2][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2 - 1][32 * ii2][i4])) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii2 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2][32 * ii2 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][32 * ii2][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][32 * ii2][i4])) + A[i2][32 * ii2][i4 - 1]))) + A[i2][32 * ii2][i4]);
                            }
                          }
                        }
                      }
                    }
                  }
                  if (i0 == _PB_TSTEPS) {
                    for (int i2 = 32 * ii2 - 1; i2 <= 32 * ii2; i2 += 1) {
                      for (int i3 = 64 * ii2 - i2; i3 < _PB_N - 1; i3 += 1) {
                        for (int i4 = max(max(1, 32 * ii2 + 32 * ii4 - i2), 64 * ii2 + 32 * ii4 - i2 - i3 + 1); i4 <= 64 * ii2 + 32 * ii4 - i2 - i3 + 32; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                        for (int i4 = 64 * ii2 + 32 * ii4 - i2 - i3 + 33; i4 <= 32 * ii2 + 32 * ii4 - i2 + 31; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                  }
                  for (int i2 = 32 * ii2 + 1; i2 < _PB_N - 1; i2 += 1) {
                    if (i0 == _PB_TSTEPS) {
                      for (int i3 = 32 * ii2 - 1; i3 <= 32 * ii2; i3 += 1) {
                        for (int i4 = max(1, 32 * ii2 + 32 * ii4 - i3); i4 <= 32 * ii2 + 32 * ii4 - i3 + 31; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                    for (int i3 = 32 * ii2 + 1; i3 < _PB_N - 1; i3 += 1) {
                      for (int i4 = max(1, 2 * _PB_TSTEPS + 32 * ii4 - 2 * i0 - 1); i4 <= 2 * _PB_TSTEPS + 32 * ii4 - 2 * i0 + 30; i4 += 1) {
                        A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                      }
                    }
                  }
                }
              } else {
                for (int i0 = _PB_TSTEPS - 1; i0 <= _PB_TSTEPS; i0 += 1) {
                  if (i0 == _PB_TSTEPS) {
                    for (int i2 = 32 * ii2; i2 < _PB_N - 1; i2 += 1) {
                      for (int i3 = max(32 * ii2, 64 * ii2 - i2 + 1); i3 < _PB_N - 1; i3 += 1) {
                        for (int i4 = max(max(32 * ii2, 64 * ii2 - i2 + 1), 96 * ii2 - 2 * i3 + 1); i4 < _PB_N - 1; i4 += 1) {
                          B[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((A[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((A[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((A[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * A[i2][i3][i4])) + A[i2][i3][i4 - 1]))) + A[i2][i3][i4]);
                        }
                      }
                    }
                  }
                  if (i0 == _PB_TSTEPS) {
                    for (int i2 = 32 * ii2 - 1; i2 <= 32 * ii2; i2 += 1) {
                      for (int i3 = 64 * ii2 - i2; i3 < _PB_N - 1; i3 += 1) {
                        for (int i4 = max(64 * ii2 - i2, 96 * ii2 - i2 - i3 + 1); i4 < _PB_N - 1; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                  }
                  for (int i2 = 32 * ii2 + 1; i2 < _PB_N - 1; i2 += 1) {
                    if (i0 == _PB_TSTEPS) {
                      for (int i3 = 32 * ii2 - 1; i3 <= 32 * ii2; i3 += 1) {
                        for (int i4 = 64 * ii2 - i3; i4 < _PB_N - 1; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                    }
                    for (int i3 = 32 * ii2 + 1; i3 < _PB_N - 1; i3 += 1) {
                      if (i0 == _PB_TSTEPS) {
                        for (int i4 = 32 * ii2 - 1; i4 <= 32 * ii2; i4 += 1) {
                          A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                        }
                      }
                      for (int i4 = 32 * ii2 + 1; i4 < _PB_N - 1; i4 += 1) {
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
  if (_PB_N >= 35 && 32 * ii0 + 1 == _PB_TSTEPS) {
    for (int ii2 = 0; ii2 <= (_PB_N - 3) / 32; ii2 += 1) {
      if (_PB_N >= 32 * ii2 + 35) {
        for (int ii3 = 0; ii3 <= (_PB_N - 3) / 32; ii3 += 1) {
          for (int ii4 = 0; ii4 < (_PB_N - 3) / 32; ii4 += 1) {
            for (int i2 = 32 * ii2 + 1; i2 <= 32 * ii2 + 32; i2 += 1) {
              for (int i3 = 32 * ii3 + 1; i3 <= min(_PB_N - 2, 32 * ii3 + 32); i3 += 1) {
                for (int i4 = 32 * ii4 + 1; i4 <= 32 * ii4 + 32; i4 += 1) {
                  A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                }
              }
            }
          }
          if (_PB_N >= 32 * ii3 + 35) {
            for (int i2 = 32 * ii2 + 1; i2 <= 32 * ii2 + 32; i2 += 1) {
              for (int i3 = 32 * ii3 + 1; i3 <= 32 * ii3 + 32; i3 += 1) {
                for (int i4 = -((_PB_N - 3) % 32) + _PB_N - 2; i4 < _PB_N - 1; i4 += 1) {
                  A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                }
              }
            }
          } else {
            for (int i2 = 32 * ii2 + 1; i2 <= 32 * ii2 + 32; i2 += 1) {
              for (int i3 = 32 * ii3 + 1; i3 < _PB_N - 1; i3 += 1) {
                for (int i4 = 32 * ii3 + 1; i4 < _PB_N - 1; i4 += 1) {
                  A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                }
              }
            }
          }
        }
      } else {
        for (int ii3 = 0; ii3 <= ii2; ii3 += 1) {
          if (_PB_N >= 32 * ii3 + 35 || 1) {
            for (int ii4 = 0; ii4 < ii2; ii4 += 1) {
              for (int i2 = 32 * ii2 + 1; i2 < _PB_N - 1; i2 += 1) {
                for (int i3 = 32 * ii3 + 1; i3 <= min(_PB_N - 2, 32 * ii3 + 32); i3 += 1) {
                  for (int i4 = 32 * ii4 + 1; i4 <= 32 * ii4 + 32; i4 += 1) {
                    A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                  }
                }
              }
            }
            if (ii3 == ii2) {
              for (int i2 = 32 * ii2 + 1; i2 < _PB_N - 1; i2 += 1) {
                for (int i3 = 32 * ii2 + 1; i3 < _PB_N - 1; i3 += 1) {
                  for (int i4 = 32 * ii2 + 1; i4 < _PB_N - 1; i4 += 1) {
                    A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                  }
                }
              }
            } else {
              for (int i2 = 32 * ii2 + 1; i2 < _PB_N - 1; i2 += 1) {
                for (int i3 = 32 * ii3 + 1; i3 <= 32 * ii3 + 32; i3 += 1) {
                  for (int i4 = 32 * ii2 + 1; i4 < _PB_N - 1; i4 += 1) {
                    A[i2][i3][i4] = ((((SCALAR_VAL(0.125) * ((B[i2 + 1][i3][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2 - 1][i3][i4])) + (SCALAR_VAL(0.125) * ((B[i2][i3 + 1][i4] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3 - 1][i4]))) + (SCALAR_VAL(0.125) * ((B[i2][i3][i4 + 1] - (SCALAR_VAL(2.0) * B[i2][i3][i4])) + B[i2][i3][i4 - 1]))) + B[i2][i3][i4]);
                  }
                }
              }
            }
          }
        }
      }
    }
  } else if (_PB_N <= 34) {
    for (int i0 = 32 * ii0 + 1; i0 <= min(_PB_TSTEPS, 32 * ii0 + 32); i0 += 1) {
      if (i0 >= 32 * ii0 + 2) {
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
