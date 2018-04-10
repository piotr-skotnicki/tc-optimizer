/**
 * This version is stamped on Apr. 14, 2015
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* floyd-warshall.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "floyd-warshall.h"


/* Array initialization. */
static
void init_array (int n,
		 DATA_TYPE POLYBENCH_2D(path,N,N,n,n))
{
  int i, j;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
      path[i][j] = i*j%7+1;
      if ((i+j)%13 == 0 || (i+j)%7==0 || (i+j)%11 == 0)
         path[i][j] = 999;
    }
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_2D(path,N,N,n,n))

{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("path");
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
      if ((i * n + j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
      fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, path[i][j]);
    }
  POLYBENCH_DUMP_END("path");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_floyd_warshall(int n,
			   DATA_TYPE POLYBENCH_2D(path,N,N,n,n))
{
  int i, j, k;

/* TC Optimizing Compiler 0.2.26 */
/* ./tc ../examples/polybench/floyd-warshall.scop.c --correction-tiling --lex-scheduling --serial-codegen -b 64 --debug --floyd-warshall-tc */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#pragma scop
{
  for (int ii0 = 0; ii0 < _PB_N / 64; ii0 += 1) {
    for (int ii1 = 0; ii1 < ii0; ii1 += 1) {
      for (int ii2 = 0; ii2 <= (_PB_N - 1) / 64; ii2 += 1) {
        for (int i1 = 64 * ii1; i1 <= 64 * ii1 + 63; i1 += 1) {
          for (int i2 = 64 * ii2; i2 <= min(_PB_N - 1, 64 * ii2 + 63); i2 += 1) {
            path[i1][i2] = ((path[i1][i2] < (path[i1][64 * ii0] + path[64 * ii0][i2])) ? path[i1][i2] : (path[i1][64 * ii0] + path[64 * ii0][i2]));
          }
        }
      }
    }
    for (int ii1 = ii0; ii1 < _PB_N / 64; ii1 += 1) {
      for (int ii2 = 0; ii2 <= (_PB_N - 1) / 64; ii2 += 1) {
        if (_PB_N >= 64 * ii2 + 65) {
          for (int i1 = 64 * ii1; i1 <= 64 * ii1 + 63; i1 += 1) {
            for (int i2 = 64 * ii2; i2 <= 64 * ii2 + 63; i2 += 1) {
              path[i1][i2] = ((path[i1][i2] < (path[i1][64 * ii0] + path[64 * ii0][i2])) ? path[i1][i2] : (path[i1][64 * ii0] + path[64 * ii0][i2]));
            }
          }
        } else {
          for (int i1 = 64 * ii1; i1 <= 64 * ii1 + 63; i1 += 1) {
            if (64 * ii0 + 64 == _PB_N && 64 * ii1 + 64 == _PB_N && 64 * ii2 + 64 == _PB_N && i1 + 64 == _PB_N) {
              path[_PB_N - 64][_PB_N - 64] = ((path[_PB_N - 64][_PB_N - 64] < (path[_PB_N - 64][_PB_N - 64] + path[_PB_N - 64][_PB_N - 64])) ? path[_PB_N - 64][_PB_N - 64] : (path[_PB_N - 64][_PB_N - 64] + path[_PB_N - 64][_PB_N - 64]));
            }
            if (64 * ii1 + 64 == _PB_N && 64 * ii2 + 64 == _PB_N) {
              for (int i2 = _PB_N - 64; i2 < i1; i2 += 1) {
                path[i1][i2] = ((path[i1][i2] < (path[i1][64 * ii0] + path[64 * ii0][i2])) ? path[i1][i2] : (path[i1][64 * ii0] + path[64 * ii0][i2]));
              }
            }
            for (int i2 = max(max(64 * ii0 + 1, 64 * ii2), i1); i2 < _PB_N; i2 += 1) {
              path[i1][i2] = ((path[i1][i2] < (path[i1][64 * ii0] + path[64 * ii0][i2])) ? path[i1][i2] : (path[i1][64 * ii0] + path[64 * ii0][i2]));
            }
          }
        }
        if (_PB_N >= 64 * ii0 + 128 && ii1 == ii0 && 64 * ii2 + 64 == _PB_N) {
          for (int i1 = 0; i1 < 64 * ii0; i1 += 1) {
            for (int i2 = _PB_N - 64; i2 < _PB_N; i2 += 1) {
              path[i1][i2] = ((path[i1][i2] < (path[i1][64 * ii0 + 1] + path[64 * ii0 + 1][i2])) ? path[i1][i2] : (path[i1][64 * ii0 + 1] + path[64 * ii0 + 1][i2]));
            }
          }
        } else if (ii1 == ii0 && 64 * ii2 + 63 >= _PB_N) {
          for (int i1 = 0; i1 < 64 * ii0; i1 += 1) {
            for (int i2 = 64 * ii2; i2 < _PB_N; i2 += 1) {
              path[i1][i2] = ((path[i1][i2] < (path[i1][64 * ii0 + 1] + path[64 * ii0 + 1][i2])) ? path[i1][i2] : (path[i1][64 * ii0 + 1] + path[64 * ii0 + 1][i2]));
            }
          }
        } else if (ii1 == ii0 && _PB_N >= 64 * ii2 + 65) {
          for (int i1 = 0; i1 < 64 * ii0; i1 += 1) {
            for (int i2 = 64 * ii2; i2 <= 64 * ii2 + 63; i2 += 1) {
              path[i1][i2] = ((path[i1][i2] < (path[i1][64 * ii0 + 1] + path[64 * ii0 + 1][i2])) ? path[i1][i2] : (path[i1][64 * ii0 + 1] + path[64 * ii0 + 1][i2]));
            }
          }
        } else if (_PB_N >= 64 * ii0 + 128 && 64 * ii1 + 64 == _PB_N && _PB_N >= 64 * ii2 + 128) {
          for (int i1 = 64 * ii0; i1 < _PB_N - 64; i1 += 1) {
            for (int i2 = 64 * ii2; i2 <= 64 * ii2 + 63; i2 += 1) {
              path[i1][i2] = ((path[i1][i2] < (path[i1][64 * ii0 + 1] + path[64 * ii0 + 1][i2])) ? path[i1][i2] : (path[i1][64 * ii0 + 1] + path[64 * ii0 + 1][i2]));
            }
          }
          if (ii2 == ii0) {
            for (int i1 = _PB_N - 64; i1 < _PB_N; i1 += 1) {
              for (int i2 = 0; i2 < 64 * ii0; i2 += 1) {
                path[i1][i2] = ((path[i1][i2] < (path[i1][64 * ii0 + 1] + path[64 * ii0 + 1][i2])) ? path[i1][i2] : (path[i1][64 * ii0 + 1] + path[64 * ii0 + 1][i2]));
              }
            }
          }
          if (ii0 == 0 && ii2 == 0) {
            path[0][0] = ((path[0][0] < (path[0][2] + path[2][0])) ? path[0][0] : (path[0][2] + path[2][0]));
          } else if (ii0 >= 1) {
            for (int i1 = 0; i1 < 64 * ii0; i1 += 1) {
              for (int i2 = 64 * ii2; i2 <= 64 * ii2 + 63; i2 += 1) {
                path[i1][i2] = ((path[i1][i2] < (path[i1][64 * ii0 + 2] + path[64 * ii0 + 2][i2])) ? path[i1][i2] : (path[i1][64 * ii0 + 2] + path[64 * ii0 + 2][i2]));
              }
            }
            if (ii2 == ii0) {
              for (int i1 = 64 * ii0; i1 < _PB_N - 64; i1 += 1) {
                if (i1 >= 64 * ii0 + 1) {
                  for (int i2 = 0; i2 < 64 * ii0; i2 += 1) {
                    path[i1][i2] = ((path[i1][i2] < (path[i1][64 * ii0 + 2] + path[64 * ii0 + 2][i2])) ? path[i1][i2] : (path[i1][64 * ii0 + 2] + path[64 * ii0 + 2][i2]));
                  }
                } else {
                  for (int i2 = 0; i2 <= 64 * ii0; i2 += 1) {
                    path[64 * ii0][i2] = ((path[64 * ii0][i2] < (path[64 * ii0][64 * ii0 + 2] + path[64 * ii0 + 2][i2])) ? path[64 * ii0][i2] : (path[64 * ii0][64 * ii0 + 2] + path[64 * ii0 + 2][i2]));
                  }
                }
              }
            }
          }
          if (ii2 == ii0) {
            for (int i1 = 0; i1 < 64 * ii0; i1 += 1) {
              for (int i2 = 0; i2 < 64 * ii0; i2 += 1) {
                path[i1][i2] = ((path[i1][i2] < (path[i1][64 * ii0 + 3] + path[64 * ii0 + 3][i2])) ? path[i1][i2] : (path[i1][64 * ii0 + 3] + path[64 * ii0 + 3][i2]));
              }
            }
          }
        } else if (64 * ii1 + 64 == _PB_N && 64 * ii2 + 64 == _PB_N) {
          for (int i0 = 64 * ii0 + 1; i0 <= min(_PB_N - 2, 64 * ii0 + 63); i0 += 1) {
            if (i0 >= 64 * ii0 + 2 && 64 * ii0 + 3 >= i0) {
              for (int i1 = 0; i1 < 64 * ii0; i1 += 1) {
                if (i0 == 64 * ii0 + 3) {
                  if (64 * ii0 + 64 == _PB_N) {
                    for (int i2 = 0; i2 < _PB_N - 64; i2 += 1) {
                      path[i1][i2] = ((path[i1][i2] < (path[i1][_PB_N - 61] + path[_PB_N - 61][i2])) ? path[i1][i2] : (path[i1][_PB_N - 61] + path[_PB_N - 61][i2]));
                    }
                  }
                  for (int i2 = 64 * ii0; i2 < _PB_N; i2 += 1) {
                    path[i1][i2] = ((path[i1][i2] < (path[i1][64 * ii0 + 3] + path[64 * ii0 + 3][i2])) ? path[i1][i2] : (path[i1][64 * ii0 + 3] + path[64 * ii0 + 3][i2]));
                  }
                } else if (_PB_N >= 64 * ii0 + 128) {
                  for (int i2 = _PB_N - 64; i2 < _PB_N; i2 += 1) {
                    path[i1][i2] = ((path[i1][i2] < (path[i1][64 * ii0 + 2] + path[64 * ii0 + 2][i2])) ? path[i1][i2] : (path[i1][64 * ii0 + 2] + path[64 * ii0 + 2][i2]));
                  }
                } else {
                  for (int i2 = 0; i2 < _PB_N; i2 += 1) {
                    path[i1][i2] = ((path[i1][i2] < (path[i1][_PB_N - 62] + path[_PB_N - 62][i2])) ? path[i1][i2] : (path[i1][_PB_N - 62] + path[_PB_N - 62][i2]));
                  }
                }
              }
              if (_PB_N >= 64 * ii0 + 128 && i0 == 64 * ii0 + 2) {
                for (int i2 = 64 * ii0 + 1; i2 < _PB_N; i2 += 1) {
                  path[64 * ii0][i2] = ((path[64 * ii0][i2] < (path[64 * ii0][64 * ii0 + 2] + path[64 * ii0 + 2][i2])) ? path[64 * ii0][i2] : (path[64 * ii0][64 * ii0 + 2] + path[64 * ii0 + 2][i2]));
                }
              }
              if (_PB_N >= 64 * ii0 + 128 && i0 == 64 * ii0 + 3) {
                for (int i2 = 0; i2 < _PB_N; i2 += 1) {
                  path[64 * ii0][i2] = ((path[64 * ii0][i2] < (path[64 * ii0][64 * ii0 + 3] + path[64 * ii0 + 3][i2])) ? path[64 * ii0][i2] : (path[64 * ii0][64 * ii0 + 3] + path[64 * ii0 + 3][i2]));
                }
              }
            } else {
              if (i0 == 64 * ii0 + 1) {
                for (int i1 = 64 * ii0; i1 < _PB_N - 64; i1 += 1) {
                  for (int i2 = _PB_N - 64; i2 < _PB_N; i2 += 1) {
                    path[i1][i2] = ((path[i1][i2] < (path[i1][64 * ii0 + 1] + path[64 * ii0 + 1][i2])) ? path[i1][i2] : (path[i1][64 * ii0 + 1] + path[64 * ii0 + 1][i2]));
                  }
                }
              }
              if (_PB_N >= 64 * ii0 + 128 && i0 >= 64 * ii0 + 4) {
                for (int i1 = 0; i1 <= 64 * ii0; i1 += 1) {
                  for (int i2 = 0; i2 < _PB_N; i2 += 1) {
                    path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
                  }
                }
              } else if (_PB_N == 64 && ii0 == 0 && i0 == 1) {
                for (int i2 = 0; i2 <= 63; i2 += 1) {
                  path[0][i2] = ((path[0][i2] < (path[0][1] + path[1][i2])) ? path[0][i2] : (path[0][1] + path[1][i2]));
                }
              } else if (_PB_N >= 64 * ii0 + 128 && i0 == 64 * ii0 + 1) {
                for (int i2 = 64 * ii0; i2 < _PB_N; i2 += 1) {
                  path[_PB_N - 64][i2] = ((path[_PB_N - 64][i2] < (path[_PB_N - 64][64 * ii0 + 1] + path[64 * ii0 + 1][i2])) ? path[_PB_N - 64][i2] : (path[_PB_N - 64][64 * ii0 + 1] + path[64 * ii0 + 1][i2]));
                }
              }
            }
            if (i0 >= 64 * ii0 + 2) {
              for (int i1 = 64 * ii0 + 1; i1 < _PB_N - 64; i1 += 1) {
                if (i0 >= 64 * ii0 + 3) {
                  for (int i2 = 0; i2 < _PB_N; i2 += 1) {
                    path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
                  }
                } else {
                  for (int i2 = 64 * ii0; i2 < _PB_N; i2 += 1) {
                    path[i1][i2] = ((path[i1][i2] < (path[i1][64 * ii0 + 2] + path[64 * ii0 + 2][i2])) ? path[i1][i2] : (path[i1][64 * ii0 + 2] + path[64 * ii0 + 2][i2]));
                  }
                }
              }
              if (_PB_N >= 64 * ii0 + 128 && i0 == 64 * ii0 + 2) {
                for (int i2 = 0; i2 < _PB_N; i2 += 1) {
                  path[_PB_N - 64][i2] = ((path[_PB_N - 64][i2] < (path[_PB_N - 64][64 * ii0 + 2] + path[64 * ii0 + 2][i2])) ? path[_PB_N - 64][i2] : (path[_PB_N - 64][64 * ii0 + 2] + path[64 * ii0 + 2][i2]));
                }
              }
              if (64 * ii0 + 64 == _PB_N && i0 + 60 >= _PB_N) {
                for (int i1 = 0; i1 < _PB_N - 64; i1 += 1) {
                  for (int i2 = 0; i2 < _PB_N; i2 += 1) {
                    path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
                  }
                }
              }
            } else if (64 * ii0 + 64 == _PB_N) {
              for (int i1 = 0; i1 < _PB_N - 64; i1 += 1) {
                for (int i2 = _PB_N - 64; i2 < _PB_N; i2 += 1) {
                  path[i1][i2] = ((path[i1][i2] < (path[i1][_PB_N - 63] + path[_PB_N - 63][i2])) ? path[i1][i2] : (path[i1][_PB_N - 63] + path[_PB_N - 63][i2]));
                }
              }
            }
            if (64 * ii0 + 64 == _PB_N && _PB_N >= i0 + 62 && i0 >= 2) {
              for (int i2 = 0; i2 < i0 - 1; i2 += 1) {
                path[_PB_N - 64][i2] = ((path[_PB_N - 64][i2] < (path[_PB_N - 64][i0] + path[i0][i2])) ? path[_PB_N - 64][i2] : (path[_PB_N - 64][i0] + path[i0][i2]));
              }
              for (int i2 = i0 - 1; i2 < 62 * _PB_N - 61 * i0 - 3843; i2 += 1) {
                path[_PB_N - 64][i2] = ((path[_PB_N - 64][i2] < (path[_PB_N - 64][i0] + path[i0][i2])) ? path[_PB_N - 64][i2] : (path[_PB_N - 64][i0] + path[i0][i2]));
              }
              if (i0 + 62 == _PB_N) {
                for (int i2 = _PB_N - 61; i2 < _PB_N; i2 += 1) {
                  path[_PB_N - 64][i2] = ((path[_PB_N - 64][i2] < (path[_PB_N - 64][_PB_N - 62] + path[_PB_N - 62][i2])) ? path[_PB_N - 64][i2] : (path[_PB_N - 64][_PB_N - 62] + path[_PB_N - 62][i2]));
                }
              }
            }
            if (64 * ii0 + 64 == _PB_N && _PB_N >= i0 + 62) {
              for (int i2 = 0; i2 < _PB_N - 63; i2 += 1) {
                path[_PB_N - 63][i2] = ((path[_PB_N - 63][i2] < (path[_PB_N - 63][i0] + path[i0][i2])) ? path[_PB_N - 63][i2] : (path[_PB_N - 63][i0] + path[i0][i2]));
              }
              for (int i2 = _PB_N - 63; i2 < _PB_N; i2 += 1) {
                path[_PB_N - 63][i2] = ((path[_PB_N - 63][i2] < (path[_PB_N - 63][i0] + path[i0][i2])) ? path[_PB_N - 63][i2] : (path[_PB_N - 63][i0] + path[i0][i2]));
              }
              if (i0 + 62 == _PB_N) {
                for (int i2 = 0; i2 < _PB_N; i2 += 1) {
                  path[_PB_N - 62][i2] = ((path[_PB_N - 62][i2] < (path[_PB_N - 62][_PB_N - 62] + path[_PB_N - 62][i2])) ? path[_PB_N - 62][i2] : (path[_PB_N - 62][_PB_N - 62] + path[_PB_N - 62][i2]));
                }
              }
            } else if (i0 >= 64 * ii0 + 3) {
              for (int i2 = 0; i2 < _PB_N - 64; i2 += 1) {
                path[_PB_N - 64][i2] = ((path[_PB_N - 64][i2] < (path[_PB_N - 64][i0] + path[i0][i2])) ? path[_PB_N - 64][i2] : (path[_PB_N - 64][i0] + path[i0][i2]));
              }
              if (64 * ii0 + 64 == _PB_N) {
                for (int i2 = _PB_N - 64; i2 <= i0; i2 += 1) {
                  path[_PB_N - 64][i2] = ((path[_PB_N - 64][i2] < (path[_PB_N - 64][i0] + path[i0][i2])) ? path[_PB_N - 64][i2] : (path[_PB_N - 64][i0] + path[i0][i2]));
                }
              }
              for (int i2 = max(_PB_N - 64, i0 + 1); i2 < _PB_N; i2 += 1) {
                path[_PB_N - 64][i2] = ((path[_PB_N - 64][i2] < (path[_PB_N - 64][i0] + path[i0][i2])) ? path[_PB_N - 64][i2] : (path[_PB_N - 64][i0] + path[i0][i2]));
              }
              if (64 * ii0 + 64 == _PB_N) {
                for (int i1 = _PB_N - 63; i1 <= i0; i1 += 1) {
                  for (int i2 = 0; i2 < _PB_N - 64; i2 += 1) {
                    path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
                  }
                  for (int i2 = _PB_N - 64; i2 < _PB_N; i2 += 1) {
                    path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
                  }
                }
              }
            }
            for (int i1 = max(_PB_N - 63, i0 + 1); i1 < _PB_N; i1 += 1) {
              if (64 * ii0 + 64 == _PB_N && i0 + 63 == _PB_N) {
                for (int i2 = 0; i2 < _PB_N - 63; i2 += 1) {
                  path[i1][i2] = ((path[i1][i2] < (path[i1][_PB_N - 63] + path[_PB_N - 63][i2])) ? path[i1][i2] : (path[i1][_PB_N - 63] + path[_PB_N - 63][i2]));
                }
              }
              if (i0 == 64 * ii0 + 1) {
                for (int i2 = 64 * ii0; i2 < _PB_N - 64; i2 += 1) {
                  path[i1][i2] = ((path[i1][i2] < (path[i1][64 * ii0 + 1] + path[64 * ii0 + 1][i2])) ? path[i1][i2] : (path[i1][64 * ii0 + 1] + path[64 * ii0 + 1][i2]));
                }
              } else {
                for (int i2 = 0; i2 < _PB_N - 64; i2 += 1) {
                  path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
                }
                if (64 * ii0 + 64 == _PB_N) {
                  for (int i2 = _PB_N - 64; i2 < i0; i2 += 1) {
                    path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
                  }
                }
              }
              for (int i2 = max(_PB_N - 64, i0); i2 < i1; i2 += 1) {
                path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
              }
              for (int i2 = i1; i2 < _PB_N; i2 += 1) {
                path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
              }
            }
          }
          if (64 * ii0 + 64 == _PB_N) {
            for (int i1 = 0; i1 < min(3845 * _PB_N - 249860, _PB_N); i1 += 1) {
              for (int i2 = 0; i2 < _PB_N; i2 += 1) {
                path[i1][i2] = ((path[i1][i2] < (path[i1][_PB_N - 1] + path[_PB_N - 1][i2])) ? path[i1][i2] : (path[i1][_PB_N - 1] + path[_PB_N - 1][i2]));
              }
            }
            if (_PB_N == 64) {
              for (int i1 = 0; i1 <= 63; i1 += 1) {
                for (int i2 = 0; i2 <= 63; i2 += 1) {
                  path[i1][i2] = ((path[i1][i2] < (path[i1][63] + path[63][i2])) ? path[i1][i2] : (path[i1][63] + path[63][i2]));
                }
              }
            }
          }
        }
      }
    }
    if (_PB_N % 64 >= 1) {
      for (int ii2 = 0; ii2 < (_PB_N - 1) / 64; ii2 += 1) {
        for (int i1 = -((_PB_N - 1) % 64) + _PB_N - 1; i1 < _PB_N; i1 += 1) {
          for (int i2 = 64 * ii2; i2 <= 64 * ii2 + 63; i2 += 1) {
            path[i1][i2] = ((path[i1][i2] < (path[i1][64 * ii0] + path[64 * ii0][i2])) ? path[i1][i2] : (path[i1][64 * ii0] + path[64 * ii0][i2]));
          }
        }
        if (ii0 >= ii2 + 1) {
          for (int i1 = 64 * ii0; i1 < -((_PB_N - 1) % 64) + _PB_N - 1; i1 += 1) {
            for (int i2 = 64 * ii2; i2 <= 64 * ii2 + 63; i2 += 1) {
              path[i1][i2] = ((path[i1][i2] < (path[i1][64 * ii0 + 1] + path[64 * ii0 + 1][i2])) ? path[i1][i2] : (path[i1][64 * ii0 + 1] + path[64 * ii0 + 1][i2]));
            }
          }
        } else {
          for (int i1 = 64 * ii0; i1 < -((_PB_N - 1) % 64) + _PB_N - 1; i1 += 1) {
            for (int i2 = 64 * ii2; i2 <= min(64 * ii2 + 63, i1 - 1); i2 += 1) {
              if (i2 >= 64 * ii0 + 1 || 1) {
                path[i1][i2] = ((path[i1][i2] < (path[i1][64 * ii0 + 1] + path[64 * ii0 + 1][i2])) ? path[i1][i2] : (path[i1][64 * ii0 + 1] + path[64 * ii0 + 1][i2]));
              }
            }
            for (int i2 = max(max(64 * ii0 + 64, 64 * ii2), i1); i2 <= 64 * ii2 + 63; i2 += 1) {
              path[i1][i2] = ((path[i1][i2] < (path[i1][64 * ii0 + 1] + path[64 * ii0 + 1][i2])) ? path[i1][i2] : (path[i1][64 * ii0 + 1] + path[64 * ii0 + 1][i2]));
            }
            if (ii2 == ii0) {
              for (int i2 = i1; i2 <= 64 * ii0 + 63; i2 += 1) {
                path[i1][i2] = ((path[i1][i2] < (path[i1][64 * ii0 + 1] + path[64 * ii0 + 1][i2])) ? path[i1][i2] : (path[i1][64 * ii0 + 1] + path[64 * ii0 + 1][i2]));
              }
            }
          }
          if (ii2 == ii0) {
            for (int i1 = -((_PB_N - 1) % 64) + _PB_N - 1; i1 < _PB_N; i1 += 1) {
              for (int i2 = 0; i2 < 64 * ii0; i2 += 1) {
                path[i1][i2] = ((path[i1][i2] < (path[i1][64 * ii0 + 1] + path[64 * ii0 + 1][i2])) ? path[i1][i2] : (path[i1][64 * ii0 + 1] + path[64 * ii0 + 1][i2]));
              }
            }
          }
        }
        if (ii0 == 0 && ii2 == 0) {
          path[0][0] = ((path[0][0] < (path[0][2] + path[2][0])) ? path[0][0] : (path[0][2] + path[2][0]));
        } else if (ii0 >= 1) {
          for (int i1 = 0; i1 < 64 * ii0; i1 += 1) {
            for (int i2 = 64 * ii2; i2 <= 64 * ii2 + 63; i2 += 1) {
              path[i1][i2] = ((path[i1][i2] < (path[i1][64 * ii0 + 2] + path[64 * ii0 + 2][i2])) ? path[i1][i2] : (path[i1][64 * ii0 + 2] + path[64 * ii0 + 2][i2]));
            }
          }
          if (ii2 == ii0) {
            for (int i2 = 0; i2 <= 64 * ii0; i2 += 1) {
              path[64 * ii0][i2] = ((path[64 * ii0][i2] < (path[64 * ii0][64 * ii0 + 2] + path[64 * ii0 + 2][i2])) ? path[64 * ii0][i2] : (path[64 * ii0][64 * ii0 + 2] + path[64 * ii0 + 2][i2]));
            }
            for (int i1 = 64 * ii0 + 1; i1 < -((_PB_N - 1) % 64) + _PB_N - 1; i1 += 1) {
              for (int i2 = 0; i2 < 64 * ii0; i2 += 1) {
                path[i1][i2] = ((path[i1][i2] < (path[i1][64 * ii0 + 2] + path[64 * ii0 + 2][i2])) ? path[i1][i2] : (path[i1][64 * ii0 + 2] + path[64 * ii0 + 2][i2]));
              }
            }
          }
        }
        if (ii2 == ii0) {
          for (int i1 = 0; i1 < 64 * ii0; i1 += 1) {
            for (int i2 = 0; i2 < 64 * ii0; i2 += 1) {
              path[i1][i2] = ((path[i1][i2] < (path[i1][64 * ii0 + 3] + path[64 * ii0 + 3][i2])) ? path[i1][i2] : (path[i1][64 * ii0 + 3] + path[64 * ii0 + 3][i2]));
            }
          }
        }
      }
      for (int i0 = 64 * ii0; i0 <= 64 * ii0 + 63; i0 += 1) {
        if (i0 == 64 * ii0 + 1) {
          for (int i1 = 64 * ii0; i1 < -((_PB_N - 1) % 64) + _PB_N - 1; i1 += 1) {
            for (int i2 = -((_PB_N - 1) % 64) + _PB_N - 1; i2 < _PB_N; i2 += 1) {
              path[i1][i2] = ((path[i1][i2] < (path[i1][64 * ii0 + 1] + path[64 * ii0 + 1][i2])) ? path[i1][i2] : (path[i1][64 * ii0 + 1] + path[64 * ii0 + 1][i2]));
            }
          }
        } else if (i0 >= 64 * ii0 + 2) {
          for (int i1 = 0; i1 < -((_PB_N - 1) % 64) + _PB_N - 1; i1 += 1) {
            if (i0 == 64 * ii0 + 3 && i1 >= 64 * ii0 && 64 * ii0 + 2 >= i1) {
              for (int i2 = 0; i2 < 64 * ii0; i2 += 1) {
                path[i1][i2] = ((path[i1][i2] < (path[i1][64 * ii0 + 3] + path[64 * ii0 + 3][i2])) ? path[i1][i2] : (path[i1][64 * ii0 + 3] + path[64 * ii0 + 3][i2]));
              }
              if (i1 == 64 * ii0) {
                path[64 * ii0][64 * ii0] = ((path[64 * ii0][64 * ii0] < (path[64 * ii0][64 * ii0 + 3] + path[64 * ii0 + 3][64 * ii0])) ? path[64 * ii0][64 * ii0] : (path[64 * ii0][64 * ii0 + 3] + path[64 * ii0 + 3][64 * ii0]));
              }
            } else if (i0 >= 64 * ii0 + 4 && i0 >= i1 + 1) {
              for (int i2 = 0; i2 < i1; i2 += 1) {
                path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
              }
              if (64 * ii0 >= i1 + 1) {
                for (int i2 = i1; i2 <= i0; i2 += 1) {
                  path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
                }
                for (int i2 = i0 + 1; i2 < -((_PB_N - 1) % 64) + _PB_N - 1; i2 += 1) {
                  path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
                }
              }
              if (i1 == 64 * ii0) {
                path[64 * ii0][64 * ii0] = ((path[64 * ii0][64 * ii0] < (path[64 * ii0][i0] + path[i0][64 * ii0])) ? path[64 * ii0][64 * ii0] : (path[64 * ii0][i0] + path[i0][64 * ii0]));
              }
            } else if (i0 == 64 * ii0 + 3 && 64 * ii0 >= i1 + 1) {
              for (int i2 = 64 * ii0; i2 < -((_PB_N - 1) % 64) + _PB_N - 1; i2 += 1) {
                path[i1][i2] = ((path[i1][i2] < (path[i1][64 * ii0 + 3] + path[64 * ii0 + 3][i2])) ? path[i1][i2] : (path[i1][64 * ii0 + 3] + path[64 * ii0 + 3][i2]));
              }
            } else if (i0 >= 64 * ii0 + 3 && i1 >= i0) {
              for (int i2 = 0; i2 < 64 * ii0; i2 += 1) {
                path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
              }
            }
            if (64 * ii0 + 3 >= i0 && i0 >= i1 + 1) {
              for (int i2 = 64 * ii0; i2 < i1; i2 += 1) {
                path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
              }
            } else if (i1 >= i0) {
              for (int i2 = 64 * ii0; i2 < i1; i2 += 1) {
                path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
              }
            }
            if (i1 >= 64 * ii0) {
              for (int i2 = max(64 * ii0 + 1, i1); i2 < -((_PB_N - 1) % 64) + _PB_N - 1; i2 += 1) {
                path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
              }
            }
            for (int i2 = -((_PB_N - 1) % 64) + _PB_N - 1; i2 < _PB_N; i2 += 1) {
              path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
            }
          }
        }
        for (int i1 = -((_PB_N - 1) % 64) + _PB_N - 1; i1 < _PB_N; i1 += 1) {
          if (i0 >= 64 * ii0 + 2) {
            for (int i2 = 0; i2 < 64 * ii0; i2 += 1) {
              path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
            }
          }
          if (i0 >= 64 * ii0 + 1) {
            for (int i2 = 64 * ii0; i2 < -((_PB_N - 1) % 64) + _PB_N - 1; i2 += 1) {
              path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
            }
          }
          for (int i2 = -((_PB_N - 1) % 64) + _PB_N - 1; i2 < _PB_N; i2 += 1) {
            path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
          }
        }
      }
    }
  }
  if (_PB_N % 64 >= 1) {
    for (int ii1 = 0; ii1 < (_PB_N - 1) / 64; ii1 += 1) {
      for (int ii2 = 0; ii2 < (_PB_N - 1) / 64; ii2 += 1) {
        for (int i1 = 64 * ii1; i1 <= 64 * ii1 + 63; i1 += 1) {
          for (int i2 = 64 * ii2; i2 <= 64 * ii2 + 63; i2 += 1) {
            path[i1][i2] = ((path[i1][i2] < (path[i1][-((_PB_N + 64) % 64) + _PB_N] + path[-((_PB_N + 64) % 64) + _PB_N][i2])) ? path[i1][i2] : (path[i1][-((_PB_N + 64) % 64) + _PB_N] + path[-((_PB_N + 64) % 64) + _PB_N][i2]));
          }
        }
      }
      for (int i1 = 64 * ii1; i1 <= 64 * ii1 + 63; i1 += 1) {
        for (int i2 = -((_PB_N + 64) % 64) + _PB_N; i2 < _PB_N; i2 += 1) {
          path[i1][i2] = ((path[i1][i2] < (path[i1][-((_PB_N + 64) % 64) + _PB_N] + path[-((_PB_N + 64) % 64) + _PB_N][i2])) ? path[i1][i2] : (path[i1][-((_PB_N + 64) % 64) + _PB_N] + path[-((_PB_N + 64) % 64) + _PB_N][i2]));
        }
      }
    }
    for (int ii2 = 0; ii2 < (_PB_N - 1) / 64; ii2 += 1) {
      for (int i0 = _PB_N - 63; i0 < _PB_N; i0 += 1) {
        if (i0 % 64 == 0) {
          for (int i1 = i0; i1 < _PB_N; i1 += 1) {
            for (int i2 = 64 * ii2; i2 <= 64 * ii2 + 63; i2 += 1) {
              path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
            }
          }
        } else if ((i0 - 1) % 64 == 0) {
          for (int i1 = 0; i1 < i0 - 1; i1 += 1) {
            for (int i2 = 64 * ii2; i2 <= 64 * ii2 + 63; i2 += 1) {
              path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
            }
          }
        }
      }
    }
    for (int i0 = -((_PB_N + 63) % 64) + _PB_N - 1; i0 < _PB_N; i0 += 1) {
      if ((_PB_N % 64) + i0 >= _PB_N + 1) {
        for (int i1 = 0; i1 < -((_PB_N + 63) % 64) + _PB_N - 1; i1 += 1) {
          if ((_PB_N % 64) + i0 >= _PB_N + 2) {
            for (int i2 = 0; i2 < -((_PB_N - 1) % 64) + _PB_N - 1; i2 += 1) {
              path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
            }
          }
          for (int i2 = -((_PB_N - 1) % 64) + _PB_N - 1; i2 < _PB_N; i2 += 1) {
            path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
          }
        }
      }
      for (int i1 = -((_PB_N + 63) % 64) + _PB_N - 1; i1 < _PB_N; i1 += 1) {
        if ((_PB_N % 64) + i0 >= _PB_N + 1) {
          for (int i2 = 0; i2 < -((_PB_N + 63) % 64) + _PB_N - 1; i2 += 1) {
            path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
          }
        }
        for (int i2 = -((_PB_N + 63) % 64) + _PB_N - 1; i2 < _PB_N; i2 += 1) {
          path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
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

  /* Variable declaration/allocation. */
  POLYBENCH_2D_ARRAY_DECL(path, DATA_TYPE, N, N, n, n);


  /* Initialize array(s). */
  init_array (n, POLYBENCH_ARRAY(path));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_floyd_warshall (n, POLYBENCH_ARRAY(path));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(path)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(path);

  return 0;
}
