/**
 * This version is stamped on Apr. 14, 2015
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* jacobi-2d.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "jacobi-2d.h"


/* Array initialization. */
static
void init_array (int n,
		 DATA_TYPE POLYBENCH_2D(A,N,N,n,n),
		 DATA_TYPE POLYBENCH_2D(B,N,N,n,n))
{
  int i, j;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      {
	A[i][j] = ((DATA_TYPE) i*(j+2) + 2) / n;
	B[i][j] = ((DATA_TYPE) i*(j+3) + 3) / n;
      }
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
      if ((i * n + j) % 20 == 0) fprintf(POLYBENCH_DUMP_TARGET, "\n");
      fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, A[i][j]);
    }
  POLYBENCH_DUMP_END("A");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_jacobi_2d(int tsteps,
			    int n,
			    DATA_TYPE POLYBENCH_2D(A,N,N,n,n),
			    DATA_TYPE POLYBENCH_2D(B,N,N,n,n))
{
  int t, i, j;

/* TC Optimizing Compiler 0.2.26 */
/* ./tc ../examples/polybench/jacobi-2d.scop.c --correction-tiling --lex-scheduling --serial-codegen -b 32 --debug */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#pragma scop
for (int ii0 = 0; ii0 <= floord(_PB_TSTEPS - 1, 32); ii0 += 1) {
  for (int ii2 = 0; ii2 <= floord(_PB_N - 3, 32); ii2 += 1) {
    for (int ii3 = 0; ii3 <= (_PB_N - 3) / 32; ii3 += 1) {
      for (int i2 = 32 * ii2 + 1; i2 <= min(_PB_N - 2, 32 * ii2 + 32); i2 += 1) {
        for (int i3 = 32 * ii3 + 1; i3 <= min(_PB_N - 2, 32 * ii3 + 32); i3 += 1) {
          B[i2][i3] = (SCALAR_VAL(0.2) * ((((A[i2][i3] + A[i2][i3 - 1]) + A[i2][i3 + 1]) + A[i2 + 1][i3]) + A[i2 - 1][i3]));
        }
      }
    }
  }
  for (int ii2 = 0; ii2 <= floord(_PB_N - 3, 32); ii2 += 1) {
    if (_PB_N >= 32 * ii2 + 35) {
      for (int i0 = 32 * ii0; i0 <= min(min(_PB_TSTEPS - 1, 32 * ii0 + 31), 32 * ii0 + 15 * ii2 + 16); i0 += 1) {
        if (i0 >= 32 * ii0 + 1) {
          for (int i2 = max(1, 64 * ii0 + 32 * ii2 - 2 * i0 + 2); i2 <= 64 * ii0 + 32 * ii2 - 2 * i0 + 33; i2 += 1) {
            if (i0 >= 32 * ii0 + 2 && 32 * ii2 >= i2) {
              for (int i3 = 1; i3 <= 64 * ii0 + 32 * ii2 - 2 * i0 - i2 + 34; i3 += 1) {
                B[i2][i3] = (SCALAR_VAL(0.2) * ((((A[i2][i3] + A[i2][i3 - 1]) + A[i2][i3 + 1]) + A[i2 + 1][i3]) + A[i2 - 1][i3]));
              }
            } else if (i2 >= 32 * ii2 + 1) {
              if (i0 == 32 * ii0 + 1 && i2 == 32 * ii2 + 1) {
                B[32 * ii2 + 1][1] = (SCALAR_VAL(0.2) * ((((A[32 * ii2 + 1][1] + A[32 * ii2 + 1][0]) + A[32 * ii2 + 1][2]) + A[32 * ii2 + 2][1]) + A[32 * ii2][1]));
              }
              for (int i3 = max(1, 32 * ii0 + 32 * ii2 - i0 - i2 + 4); i3 <= 64 * ii0 - 2 * i0 + 33; i3 += 1) {
                B[i2][i3] = (SCALAR_VAL(0.2) * ((((A[i2][i3] + A[i2][i3 - 1]) + A[i2][i3 + 1]) + A[i2 + 1][i3]) + A[i2 - 1][i3]));
              }
            } else {
              for (int i3 = 1; i3 <= 32; i3 += 1) {
                B[32 * ii2][i3] = (SCALAR_VAL(0.2) * ((((A[32 * ii2][i3] + A[32 * ii2][i3 - 1]) + A[32 * ii2][i3 + 1]) + A[32 * ii2 + 1][i3]) + A[32 * ii2 - 1][i3]));
              }
            }
          }
        }
        for (int i2 = max(1, 64 * ii0 + 32 * ii2 - 2 * i0 + 1); i2 <= 64 * ii0 + 32 * ii2 - 2 * i0 + 32; i2 += 1) {
          if (i2 >= 32 * ii2 + 1) {
            for (int i3 = 1; i3 <= 64 * ii0 - 2 * i0 + 32; i3 += 1) {
              A[i2][i3] = (SCALAR_VAL(0.2) * ((((B[i2][i3] + B[i2][i3 - 1]) + B[i2][i3 + 1]) + B[i2 + 1][i3]) + B[i2 - 1][i3]));
            }
          } else {
            for (int i3 = 1; i3 <= 64 * ii0 + 32 * ii2 - 2 * i0 - i2 + 33; i3 += 1) {
              A[i2][i3] = (SCALAR_VAL(0.2) * ((((B[i2][i3] + B[i2][i3 - 1]) + B[i2][i3 + 1]) + B[i2 + 1][i3]) + B[i2 - 1][i3]));
            }
          }
        }
      }
      if (32 * ii0 + 1 == _PB_TSTEPS) {
        for (int ii3 = 1; ii3 <= (_PB_N - 3) / 32; ii3 += 1) {
          for (int i2 = 32 * ii2 + 1; i2 <= 32 * ii2 + 32; i2 += 1) {
            for (int i3 = 32 * ii3 + 1; i3 <= min(_PB_N - 2, 32 * ii3 + 32); i3 += 1) {
              A[i2][i3] = (SCALAR_VAL(0.2) * ((((B[i2][i3] + B[i2][i3 - 1]) + B[i2][i3 + 1]) + B[i2 + 1][i3]) + B[i2 - 1][i3]));
            }
          }
        }
      } else {
        for (int ii3 = 1; ii3 <= (_PB_N - 3) / 32; ii3 += 1) {
          for (int i0 = 32 * ii0; i0 <= min(min(_PB_TSTEPS - 1, 32 * ii0 + 31), 32 * ii0 + 16 * ii2 + 16); i0 += 1) {
            if (i0 >= 32 * ii0 + 1) {
              if (_PB_N >= 32 * ii3 + 35) {
                for (int i2 = max(1, 64 * ii0 + 32 * ii2 - 2 * i0 + 2); i2 <= 64 * ii0 + 32 * ii2 - 2 * i0 + 33; i2 += 1) {
                  if (i0 >= 32 * ii0 + 2 && i2 >= 32 * ii2 + 1) {
                    for (int i3 = 64 * ii0 + 32 * ii3 - 2 * i0 + 2; i3 <= 64 * ii0 + 32 * ii3 - 2 * i0 + 33; i3 += 1) {
                      B[i2][i3] = (SCALAR_VAL(0.2) * ((((A[i2][i3] + A[i2][i3 - 1]) + A[i2][i3 + 1]) + A[i2 + 1][i3]) + A[i2 - 1][i3]));
                    }
                  } else if (i0 >= 32 * ii0 + 2 && 32 * ii2 >= i2) {
                    for (int i3 = 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 3; i3 <= 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 34; i3 += 1) {
                      B[i2][i3] = (SCALAR_VAL(0.2) * ((((A[i2][i3] + A[i2][i3 - 1]) + A[i2][i3 + 1]) + A[i2 + 1][i3]) + A[i2 - 1][i3]));
                    }
                  } else if (i2 >= 32 * ii2 + 1) {
                    B[i2][32 * ii3] = (SCALAR_VAL(0.2) * ((((A[i2][32 * ii3] + A[i2][32 * ii3 - 1]) + A[i2][32 * ii3 + 1]) + A[i2 + 1][32 * ii3]) + A[i2 - 1][32 * ii3]));
                    if (i2 == 32 * ii2 + 1) {
                      B[32 * ii2 + 1][32 * ii3 + 1] = (SCALAR_VAL(0.2) * ((((A[32 * ii2 + 1][32 * ii3 + 1] + A[32 * ii2 + 1][32 * ii3]) + A[32 * ii2 + 1][32 * ii3 + 2]) + A[32 * ii2 + 2][32 * ii3 + 1]) + A[32 * ii2][32 * ii3 + 1]));
                    }
                    for (int i3 = max(32 * ii3 + 1, 32 * ii2 + 32 * ii3 - i2 + 3); i3 <= 32 * ii3 + 31; i3 += 1) {
                      B[i2][i3] = (SCALAR_VAL(0.2) * ((((A[i2][i3] + A[i2][i3 - 1]) + A[i2][i3 + 1]) + A[i2 + 1][i3]) + A[i2 - 1][i3]));
                    }
                  } else {
                    for (int i3 = 32 * ii3 + 1; i3 <= 32 * ii3 + 32; i3 += 1) {
                      B[32 * ii2][i3] = (SCALAR_VAL(0.2) * ((((A[32 * ii2][i3] + A[32 * ii2][i3 - 1]) + A[32 * ii2][i3 + 1]) + A[32 * ii2 + 1][i3]) + A[32 * ii2 - 1][i3]));
                    }
                  }
                }
              } else {
                for (int i2 = max(1, 64 * ii0 + 32 * ii2 - 2 * i0 + 2); i2 <= 64 * ii0 + 32 * ii2 - 2 * i0 + 33; i2 += 1) {
                  for (int i3 = max(64 * ii0 + 32 * ii3 - 2 * i0 + 2, 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 3); i3 < _PB_N - 1; i3 += 1) {
                    B[i2][i3] = (SCALAR_VAL(0.2) * ((((A[i2][i3] + A[i2][i3 - 1]) + A[i2][i3 + 1]) + A[i2 + 1][i3]) + A[i2 - 1][i3]));
                  }
                }
              }
            }
            if (32 * ii0 + 16 * ii2 + 15 >= i0) {
              for (int i2 = max(1, 64 * ii0 + 32 * ii2 - 2 * i0 + 1); i2 <= min(32 * ii2, 64 * ii0 + 32 * ii2 - 2 * i0 + 32); i2 += 1) {
                if (_PB_N >= 32 * ii3 + 35) {
                  for (int i3 = 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 2; i3 <= 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 33; i3 += 1) {
                    A[i2][i3] = (SCALAR_VAL(0.2) * ((((B[i2][i3] + B[i2][i3 - 1]) + B[i2][i3 + 1]) + B[i2 + 1][i3]) + B[i2 - 1][i3]));
                  }
                } else {
                  for (int i3 = 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 2; i3 < _PB_N - 1; i3 += 1) {
                    A[i2][i3] = (SCALAR_VAL(0.2) * ((((B[i2][i3] + B[i2][i3 - 1]) + B[i2][i3 + 1]) + B[i2 + 1][i3]) + B[i2 - 1][i3]));
                  }
                }
              }
              for (int i2 = 32 * ii2 + 1; i2 <= 64 * ii0 + 32 * ii2 - 2 * i0 + 32; i2 += 1) {
                for (int i3 = 64 * ii0 + 32 * ii3 - 2 * i0 + 1; i3 <= 32 * ii3; i3 += 1) {
                  A[i2][i3] = (SCALAR_VAL(0.2) * ((((B[i2][i3] + B[i2][i3 - 1]) + B[i2][i3 + 1]) + B[i2 + 1][i3]) + B[i2 - 1][i3]));
                }
                if (32 * ii3 + 34 >= _PB_N) {
                  for (int i3 = 32 * ii3 + 1; i3 < _PB_N - 1; i3 += 1) {
                    A[i2][i3] = (SCALAR_VAL(0.2) * ((((B[i2][i3] + B[i2][i3 - 1]) + B[i2][i3 + 1]) + B[i2 + 1][i3]) + B[i2 - 1][i3]));
                  }
                } else {
                  for (int i3 = 32 * ii3 + 1; i3 <= 64 * ii0 + 32 * ii3 - 2 * i0 + 32; i3 += 1) {
                    A[i2][i3] = (SCALAR_VAL(0.2) * ((((B[i2][i3] + B[i2][i3 - 1]) + B[i2][i3 + 1]) + B[i2 + 1][i3]) + B[i2 - 1][i3]));
                  }
                }
              }
            }
          }
        }
      }
    } else {
      for (int ii3 = 0; ii3 <= ii2; ii3 += 1) {
        if (_PB_N >= 32 * ii3 + 35) {
          for (int i0 = 32 * ii0; i0 <= min(min(_PB_TSTEPS - 1, 32 * ii0 + 31), 32 * ii0 + 16 * ii3 + 15); i0 += 1) {
            if (i0 >= 32 * ii0 + 1) {
              if (i0 >= 32 * ii0 + 2) {
                for (int i2 = 64 * ii0 + 32 * ii2 - 2 * i0 + 2; i2 <= 32 * ii2; i2 += 1) {
                  for (int i3 = max(1, 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 3); i3 <= 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 34; i3 += 1) {
                    B[i2][i3] = (SCALAR_VAL(0.2) * ((((A[i2][i3] + A[i2][i3 - 1]) + A[i2][i3 + 1]) + A[i2 + 1][i3]) + A[i2 - 1][i3]));
                  }
                }
                for (int i2 = 32 * ii2 + 1; i2 < _PB_N - 1; i2 += 1) {
                  for (int i3 = max(1, 64 * ii0 + 32 * ii3 - 2 * i0 + 2); i3 <= 64 * ii0 + 32 * ii3 - 2 * i0 + 33; i3 += 1) {
                    B[i2][i3] = (SCALAR_VAL(0.2) * ((((A[i2][i3] + A[i2][i3 - 1]) + A[i2][i3 + 1]) + A[i2 + 1][i3]) + A[i2 - 1][i3]));
                  }
                }
              } else {
                for (int i2 = 32 * ii2; i2 < _PB_N - 1; i2 += 1) {
                  if (ii3 >= 1 && i2 >= 32 * ii2 + 3) {
                    B[i2][32 * ii3] = (SCALAR_VAL(0.2) * ((((A[i2][32 * ii3] + A[i2][32 * ii3 - 1]) + A[i2][32 * ii3 + 1]) + A[i2 + 1][32 * ii3]) + A[i2 - 1][32 * ii3]));
                  }
                  for (int i3 = max(max(1, 32 * ii3), 32 * ii2 + 32 * ii3 - i2 + 1); i3 <= 32 * ii2 + 32 * ii3 - i2 + 2; i3 += 1) {
                    B[i2][i3] = (SCALAR_VAL(0.2) * ((((A[i2][i3] + A[i2][i3 - 1]) + A[i2][i3 + 1]) + A[i2 + 1][i3]) + A[i2 - 1][i3]));
                  }
                  if (i2 == 32 * ii2) {
                    for (int i3 = 32 * ii3 + 3; i3 <= 32 * ii3 + 32; i3 += 1) {
                      B[32 * ii2][i3] = (SCALAR_VAL(0.2) * ((((A[32 * ii2][i3] + A[32 * ii2][i3 - 1]) + A[32 * ii2][i3 + 1]) + A[32 * ii2 + 1][i3]) + A[32 * ii2 - 1][i3]));
                    }
                  } else {
                    for (int i3 = max(32 * ii3 + 1, 32 * ii2 + 32 * ii3 - i2 + 3); i3 <= 32 * ii3 + 31; i3 += 1) {
                      B[i2][i3] = (SCALAR_VAL(0.2) * ((((A[i2][i3] + A[i2][i3 - 1]) + A[i2][i3 + 1]) + A[i2 + 1][i3]) + A[i2 - 1][i3]));
                    }
                  }
                }
              }
            }
            for (int i2 = 64 * ii0 + 32 * ii2 - 2 * i0 + 1; i2 <= 32 * ii2; i2 += 1) {
              for (int i3 = max(1, 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 2); i3 <= 64 * ii0 + 32 * ii2 + 32 * ii3 - 2 * i0 - i2 + 33; i3 += 1) {
                A[i2][i3] = (SCALAR_VAL(0.2) * ((((B[i2][i3] + B[i2][i3 - 1]) + B[i2][i3 + 1]) + B[i2 + 1][i3]) + B[i2 - 1][i3]));
              }
            }
            for (int i2 = 32 * ii2 + 1; i2 < _PB_N - 1; i2 += 1) {
              for (int i3 = max(1, 64 * ii0 + 32 * ii3 - 2 * i0 + 1); i3 <= 64 * ii0 + 32 * ii3 - 2 * i0 + 32; i3 += 1) {
                A[i2][i3] = (SCALAR_VAL(0.2) * ((((B[i2][i3] + B[i2][i3 - 1]) + B[i2][i3 + 1]) + B[i2 + 1][i3]) + B[i2 - 1][i3]));
              }
            }
          }
          if (ii3 == 0) {
            for (int i0 = 32 * ii0 + 16; i0 <= min(_PB_TSTEPS - 1, 32 * ii0 + 31); i0 += 1) {
              for (int i2 = max(1, 64 * ii0 + 32 * ii2 - 2 * i0 + 2); i2 <= min(32 * ii2, 64 * ii0 + 32 * ii2 - 2 * i0 + 33); i2 += 1) {
                for (int i3 = 1; i3 <= 64 * ii0 + 32 * ii2 - 2 * i0 - i2 + 34; i3 += 1) {
                  B[i2][i3] = (SCALAR_VAL(0.2) * ((((A[i2][i3] + A[i2][i3 - 1]) + A[i2][i3 + 1]) + A[i2 + 1][i3]) + A[i2 - 1][i3]));
                }
              }
              if (i0 == 32 * ii0 + 16) {
                for (int i2 = 32 * ii2 + 1; i2 < _PB_N - 1; i2 += 1) {
                  B[i2][1] = (SCALAR_VAL(0.2) * ((((A[i2][1] + A[i2][0]) + A[i2][2]) + A[i2 + 1][1]) + A[i2 - 1][1]));
                }
              }
              for (int i2 = max(1, 64 * ii0 + 32 * ii2 - 2 * i0 + 1); i2 <= 64 * ii0 + 32 * ii2 - 2 * i0 + 32; i2 += 1) {
                for (int i3 = 1; i3 <= 64 * ii0 + 32 * ii2 - 2 * i0 - i2 + 33; i3 += 1) {
                  A[i2][i3] = (SCALAR_VAL(0.2) * ((((B[i2][i3] + B[i2][i3 - 1]) + B[i2][i3 + 1]) + B[i2 + 1][i3]) + B[i2 - 1][i3]));
                }
              }
            }
          }
        } else {
          for (int i0 = 32 * ii0; i0 <= min(_PB_TSTEPS - 1, 32 * ii0 + 31); i0 += 1) {
            if (i0 >= 32 * ii0 + 1) {
              for (int i2 = max(1, 64 * ii0 + 32 * ii2 - 2 * i0 + 2); i2 < _PB_N - 1; i2 += 1) {
                for (int i3 = max(max(1, 64 * ii0 + 32 * ii2 - 2 * i0 + 2), 64 * ii0 + 64 * ii2 - 2 * i0 - i2 + 3); i3 < _PB_N - 1; i3 += 1) {
                  B[i2][i3] = (SCALAR_VAL(0.2) * ((((A[i2][i3] + A[i2][i3 - 1]) + A[i2][i3 + 1]) + A[i2 + 1][i3]) + A[i2 - 1][i3]));
                }
              }
            }
            for (int i2 = max(1, 64 * ii0 + 32 * ii2 - 2 * i0 + 1); i2 < _PB_N - 1; i2 += 1) {
              for (int i3 = max(max(1, 64 * ii0 + 32 * ii2 - 2 * i0 + 1), 64 * ii0 + 64 * ii2 - 2 * i0 - i2 + 2); i3 < _PB_N - 1; i3 += 1) {
                A[i2][i3] = (SCALAR_VAL(0.2) * ((((B[i2][i3] + B[i2][i3 - 1]) + B[i2][i3 + 1]) + B[i2 + 1][i3]) + B[i2 - 1][i3]));
              }
            }
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
  POLYBENCH_2D_ARRAY_DECL(A, DATA_TYPE, N, N, n, n);
  POLYBENCH_2D_ARRAY_DECL(B, DATA_TYPE, N, N, n, n);


  /* Initialize array(s). */
  init_array (n, POLYBENCH_ARRAY(A), POLYBENCH_ARRAY(B));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_jacobi_2d(tsteps, n, POLYBENCH_ARRAY(A), POLYBENCH_ARRAY(B));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(A)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(B);

  return 0;
}
