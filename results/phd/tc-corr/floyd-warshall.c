/**
 * This version is stamped on May 10, 2016
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

/* TC Optimizing Compiler 0.5.1 */
/* ./tc ../examples/polybench/floyd-warshall.scop.c --correction-tiling --free-scheduling --omp-for-codegen --floyd-warshall-tc --debug -b 16 */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#pragma scop
for (register int k = 0; k <= floord(_PB_N - 1, 16); k += 1) {
  if (k >= 2) {
    for (register int i1 = 0; i1 <= 15; i1 += 1) {
      for (register int i2 = 0; i2 <= 15; i2 += 1) {
        path[i1][i2] = ((path[i1][i2] < (path[i1][16 * k] + path[16 * k][i2])) ? path[i1][i2] : (path[i1][16 * k] + path[16 * k][i2]));
      }
    }
  } else if (k == 1) {
    #pragma omp parallel for
    for (register int ii0 = 0; ii0 <= (_PB_N - 2) / 16; ii0 += 1) {
      if (ii0 == 1) {
        for (register int i1 = 0; i1 <= 15; i1 += 1) {
          for (register int i2 = 0; i2 <= 15; i2 += 1) {
            path[i1][i2] = ((path[i1][i2] < (path[i1][16] + path[16][i2])) ? path[i1][i2] : (path[i1][16] + path[16][i2]));
          }
        }
      }
      if (_PB_N >= 16 * ii0 + 17) {
        if (ii0 >= 1) {
          for (register int i0 = 16 * ii0; i0 <= 16 * ii0 + 1; i0 += 1) {
            if (i0 == 16 * ii0 + 1) {
              for (register int i1 = 0; i1 < 16 * ii0; i1 += 1) {
                for (register int i2 = 0; i2 <= 15; i2 += 1) {
                  path[i1][i2] = ((path[i1][i2] < (path[i1][16 * ii0 + 1] + path[16 * ii0 + 1][i2])) ? path[i1][i2] : (path[i1][16 * ii0 + 1] + path[16 * ii0 + 1][i2]));
                }
              }
            } else {
              for (register int i1 = 16 * ii0; i1 <= 16 * ii0 + 15; i1 += 1) {
                for (register int i2 = 0; i2 <= 15; i2 += 1) {
                  path[i1][i2] = ((path[i1][i2] < (path[i1][16 * ii0] + path[16 * ii0][i2])) ? path[i1][i2] : (path[i1][16 * ii0] + path[16 * ii0][i2]));
                }
              }
            }
          }
        } else {
          for (register int ii2 = 1; ii2 <= (_PB_N - 1) / 16; ii2 += 1) {
            for (register int i1 = 0; i1 <= 15; i1 += 1) {
              for (register int i2 = 16 * ii2; i2 <= min(_PB_N - 1, 16 * ii2 + 15); i2 += 1) {
                path[i1][i2] = ((path[i1][i2] < (path[i1][0] + path[0][i2])) ? path[i1][i2] : (path[i1][0] + path[0][i2]));
              }
            }
          }
        }
        if (ii0 == 0) {
          for (register int ii1 = 1; ii1 < (_PB_N - 1) / 16; ii1 += 1) {
            for (register int i1 = 16 * ii1; i1 <= 16 * ii1 + 15; i1 += 1) {
              for (register int i2 = 0; i2 <= 15; i2 += 1) {
                path[i1][i2] = ((path[i1][i2] < (path[i1][0] + path[0][i2])) ? path[i1][i2] : (path[i1][0] + path[0][i2]));
              }
            }
          }
        }
      } else if (16 * ii0 + 3 >= _PB_N) {
        for (register int i0 = 16 * ii0; i0 <= 16 * ii0 + 1; i0 += 1) {
          if (i0 == 16 * ii0 + 1) {
            for (register int i1 = 0; i1 < 16 * ii0; i1 += 1) {
              for (register int i2 = 0; i2 <= 15; i2 += 1) {
                path[i1][i2] = ((path[i1][i2] < (path[i1][16 * ii0 + 1] + path[16 * ii0 + 1][i2])) ? path[i1][i2] : (path[i1][16 * ii0 + 1] + path[16 * ii0 + 1][i2]));
              }
            }
          } else {
            for (register int i1 = 16 * ii0; i1 < _PB_N; i1 += 1) {
              for (register int i2 = 0; i2 <= 15; i2 += 1) {
                path[i1][i2] = ((path[i1][i2] < (path[i1][16 * ii0] + path[16 * ii0][i2])) ? path[i1][i2] : (path[i1][16 * ii0] + path[16 * ii0][i2]));
              }
            }
          }
        }
        if (16 * ii0 + 3 == _PB_N) {
          for (register int i0 = _PB_N - 3; i0 < _PB_N; i0 += 1) {
            if (i0 + 2 >= _PB_N) {
              for (register int i1 = 0; i1 < _PB_N - 3; i1 += 1) {
                if (i0 + 1 == _PB_N) {
                  for (register int i2 = 0; i2 < _PB_N - 3; i2 += 1) {
                    path[i1][i2] = ((path[i1][i2] < (path[i1][_PB_N - 1] + path[_PB_N - 1][i2])) ? path[i1][i2] : (path[i1][_PB_N - 1] + path[_PB_N - 1][i2]));
                  }
                }
                for (register int i2 = _PB_N - 3; i2 < _PB_N; i2 += 1) {
                  path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
                }
              }
            }
            for (register int i1 = _PB_N - 3; i1 < _PB_N; i1 += 1) {
              if (i0 + 2 >= _PB_N) {
                for (register int i2 = 0; i2 < _PB_N - 3; i2 += 1) {
                  path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
                }
              }
              for (register int i2 = _PB_N - 3; i2 < _PB_N; i2 += 1) {
                path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
              }
            }
          }
        }
      }
      if (_PB_N >= 16 * ii0 + 4) {
        if (ii0 >= 1) {
          for (register int i0 = 16 * ii0; i0 <= 16 * ii0 + 1; i0 += 1) {
            if (i0 == 16 * ii0 + 1) {
              for (register int i1 = 16 * ii0; i1 < -((_PB_N - 1) % 16) + _PB_N - 1; i1 += 1) {
                for (register int i2 = 0; i2 <= 15; i2 += 1) {
                  path[i1][i2] = ((path[i1][i2] < (path[i1][16 * ii0 + 1] + path[16 * ii0 + 1][i2])) ? path[i1][i2] : (path[i1][16 * ii0 + 1] + path[16 * ii0 + 1][i2]));
                }
              }
              if (16 * ii0 + 16 >= _PB_N) {
                for (register int i1 = 0; i1 < 16 * ii0; i1 += 1) {
                  for (register int i2 = 0; i2 <= 15; i2 += 1) {
                    path[i1][i2] = ((path[i1][i2] < (path[i1][16 * ii0 + 1] + path[16 * ii0 + 1][i2])) ? path[i1][i2] : (path[i1][16 * ii0 + 1] + path[16 * ii0 + 1][i2]));
                  }
                }
              }
            } else if (_PB_N >= 16 * ii0 + 17) {
              for (register int i1 = -((_PB_N - 1) % 16) + _PB_N - 1; i1 < _PB_N; i1 += 1) {
                for (register int i2 = 0; i2 <= 15; i2 += 1) {
                  path[i1][i2] = ((path[i1][i2] < (path[i1][16 * ii0] + path[16 * ii0][i2])) ? path[i1][i2] : (path[i1][16 * ii0] + path[16 * ii0][i2]));
                }
              }
            } else {
              for (register int i1 = 16 * ii0; i1 < _PB_N; i1 += 1) {
                for (register int i2 = 0; i2 <= 15; i2 += 1) {
                  path[i1][i2] = ((path[i1][i2] < (path[i1][16 * ii0] + path[16 * ii0][i2])) ? path[i1][i2] : (path[i1][16 * ii0] + path[16 * ii0][i2]));
                }
              }
            }
          }
          if (_PB_N >= 16 * ii0 + 17) {
            for (register int i1 = 0; i1 < 16 * ii0; i1 += 1) {
              for (register int i2 = 0; i2 <= 15; i2 += 1) {
                path[i1][i2] = ((path[i1][i2] < (path[i1][16 * ii0 + 2] + path[16 * ii0 + 2][i2])) ? path[i1][i2] : (path[i1][16 * ii0 + 2] + path[16 * ii0 + 2][i2]));
              }
            }
          }
        }
        if (_PB_N >= 16 * ii0 + 17) {
          for (register int i0 = 16 * ii0; i0 <= 16 * ii0 + 2; i0 += 1) {
            if (i0 == 16 * ii0 + 2) {
              for (register int i1 = 0; i1 <= 16 * ii0; i1 += 1) {
                if (16 * ii0 >= i1 + 1) {
                  for (register int i2 = 16 * ii0; i2 <= 16 * ii0 + 15; i2 += 1) {
                    path[i1][i2] = ((path[i1][i2] < (path[i1][16 * ii0 + 2] + path[16 * ii0 + 2][i2])) ? path[i1][i2] : (path[i1][16 * ii0 + 2] + path[16 * ii0 + 2][i2]));
                  }
                } else {
                  for (register int i2 = 0; i2 <= 16 * ii0; i2 += 1) {
                    path[16 * ii0][i2] = ((path[16 * ii0][i2] < (path[16 * ii0][16 * ii0 + 2] + path[16 * ii0 + 2][i2])) ? path[16 * ii0][i2] : (path[16 * ii0][16 * ii0 + 2] + path[16 * ii0 + 2][i2]));
                  }
                }
              }
              for (register int i1 = 16 * ii0 + 1; i1 < -((_PB_N - 1) % 16) + _PB_N - 1; i1 += 1) {
                for (register int i2 = 0; i2 < 16 * ii0; i2 += 1) {
                  path[i1][i2] = ((path[i1][i2] < (path[i1][16 * ii0 + 2] + path[16 * ii0 + 2][i2])) ? path[i1][i2] : (path[i1][16 * ii0 + 2] + path[16 * ii0 + 2][i2]));
                }
              }
            } else if (i0 == 16 * ii0 + 1) {
              for (register int i1 = 16 * ii0; i1 < min(_PB_N, -((_PB_N - 1) % 16) + _PB_N + 16 * ii0 - 1); i1 += 1) {
                if (((_PB_N - 1) % 16) + i1 + 1 >= _PB_N) {
                  for (register int i2 = 0; i2 < 16 * ii0; i2 += 1) {
                    path[i1][i2] = ((path[i1][i2] < (path[i1][16 * ii0 + 1] + path[16 * ii0 + 1][i2])) ? path[i1][i2] : (path[i1][16 * ii0 + 1] + path[16 * ii0 + 1][i2]));
                  }
                } else if (ii0 >= 1) {
                  for (register int i2 = 16 * ii0; i2 <= 16 * ii0 + 15; i2 += 1) {
                    path[i1][i2] = ((path[i1][i2] < (path[i1][16 * ii0 + 1] + path[16 * ii0 + 1][i2])) ? path[i1][i2] : (path[i1][16 * ii0 + 1] + path[16 * ii0 + 1][i2]));
                  }
                } else {
                  for (register int i2 = 0; i2 <= 15; i2 += 1) {
                    path[i1][i2] = ((path[i1][i2] < (path[i1][1] + path[1][i2])) ? path[i1][i2] : (path[i1][1] + path[1][i2]));
                  }
                }
              }
            } else {
              for (register int i1 = -((_PB_N - 1) % 16) + _PB_N - 1; i1 < _PB_N; i1 += 1) {
                if (ii0 >= 1) {
                  for (register int i2 = 16 * ii0; i2 <= 16 * ii0 + 15; i2 += 1) {
                    path[i1][i2] = ((path[i1][i2] < (path[i1][16 * ii0] + path[16 * ii0][i2])) ? path[i1][i2] : (path[i1][16 * ii0] + path[16 * ii0][i2]));
                  }
                } else {
                  for (register int i2 = 0; i2 <= 15; i2 += 1) {
                    path[i1][i2] = ((path[i1][i2] < (path[i1][0] + path[0][i2])) ? path[i1][i2] : (path[i1][0] + path[0][i2]));
                  }
                }
              }
            }
          }
          for (register int i1 = 0; i1 < 16 * ii0; i1 += 1) {
            for (register int i2 = 0; i2 < 16 * ii0; i2 += 1) {
              path[i1][i2] = ((path[i1][i2] < (path[i1][16 * ii0 + 3] + path[16 * ii0 + 3][i2])) ? path[i1][i2] : (path[i1][16 * ii0 + 3] + path[16 * ii0 + 3][i2]));
            }
          }
          if (ii0 == 0) {
            for (register int ii2 = 1; ii2 < (_PB_N - 1) / 16; ii2 += 1) {
              for (register int i0 = 0; i0 <= 1; i0 += 1) {
                if (i0 == 1) {
                  for (register int i1 = 0; i1 < -((_PB_N - 1) % 16) + _PB_N - 1; i1 += 1) {
                    for (register int i2 = 16 * ii2; i2 <= 16 * ii2 + 15; i2 += 1) {
                      path[i1][i2] = ((path[i1][i2] < (path[i1][1] + path[1][i2])) ? path[i1][i2] : (path[i1][1] + path[1][i2]));
                    }
                  }
                } else {
                  for (register int i1 = -((_PB_N - 1) % 16) + _PB_N - 1; i1 < _PB_N; i1 += 1) {
                    for (register int i2 = 16 * ii2; i2 <= 16 * ii2 + 15; i2 += 1) {
                      path[i1][i2] = ((path[i1][i2] < (path[i1][0] + path[0][i2])) ? path[i1][i2] : (path[i1][0] + path[0][i2]));
                    }
                  }
                }
              }
            }
          }
        }
        for (register int i0 = 16 * ii0; i0 <= min(_PB_N - 1, 16 * ii0 + 15); i0 += 1) {
          if (_PB_N >= 16 * ii0 + 17 && i0 >= 16 * ii0 + 2) {
            for (register int i1 = 0; i1 < 16 * ii0; i1 += 1) {
              if (i0 == 16 * ii0 + 3) {
                for (register int i2 = 16 * ii0; i2 < -((_PB_N - 1) % 16) + _PB_N - 1; i2 += 1) {
                  path[i1][i2] = ((path[i1][i2] < (path[i1][16 * ii0 + 3] + path[16 * ii0 + 3][i2])) ? path[i1][i2] : (path[i1][16 * ii0 + 3] + path[16 * ii0 + 3][i2]));
                }
              } else if (i0 >= 16 * ii0 + 4) {
                for (register int i2 = 0; i2 < -((_PB_N - 1) % 16) + _PB_N - 1; i2 += 1) {
                  path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
                }
              }
              for (register int i2 = -((_PB_N - 1) % 16) + _PB_N - 1; i2 < _PB_N; i2 += 1) {
                path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
              }
            }
          }
          if (i0 >= 16 * ii0 + 1) {
            for (register int i1 = 16 * ii0; i1 < -((_PB_N - 1) % 16) + _PB_N - 1; i1 += 1) {
              if (i0 >= 16 * ii0 + 2) {
                if (i0 >= 16 * ii0 + 3) {
                  for (register int i2 = 0; i2 < 16 * ii0; i2 += 1) {
                    path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
                  }
                  if (i1 == 16 * ii0) {
                    path[16 * ii0][16 * ii0] = ((path[16 * ii0][16 * ii0] < (path[16 * ii0][i0] + path[i0][16 * ii0])) ? path[16 * ii0][16 * ii0] : (path[16 * ii0][i0] + path[i0][16 * ii0]));
                  }
                }
                if (i1 >= 16 * ii0 + 1) {
                  path[i1][16 * ii0] = ((path[i1][16 * ii0] < (path[i1][i0] + path[i0][16 * ii0])) ? path[i1][16 * ii0] : (path[i1][i0] + path[i0][16 * ii0]));
                }
                for (register int i2 = 16 * ii0 + 1; i2 < _PB_N; i2 += 1) {
                  path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
                }
              } else {
                for (register int i2 = -((_PB_N - 1) % 16) + _PB_N - 1; i2 < _PB_N; i2 += 1) {
                  path[i1][i2] = ((path[i1][i2] < (path[i1][16 * ii0 + 1] + path[16 * ii0 + 1][i2])) ? path[i1][i2] : (path[i1][16 * ii0 + 1] + path[16 * ii0 + 1][i2]));
                }
              }
            }
            if (16 * ii0 + 16 >= _PB_N) {
              for (register int i1 = 0; i1 < 16 * ii0; i1 += 1) {
                if (i0 >= 16 * ii0 + 2) {
                  for (register int i2 = 0; i2 < 16 * ii0; i2 += 1) {
                    path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
                  }
                }
                for (register int i2 = 16 * ii0; i2 < _PB_N; i2 += 1) {
                  path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
                }
              }
            }
          }
          for (register int i1 = -((_PB_N - 1) % 16) + _PB_N - 1; i1 < _PB_N; i1 += 1) {
            if (_PB_N >= 16 * ii0 + 17 && i0 >= 16 * ii0 + 2) {
              for (register int i2 = 0; i2 < 16 * ii0; i2 += 1) {
                path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
              }
            } else if (16 * ii0 + 16 >= _PB_N && i0 >= 16 * ii0 + 1) {
              for (register int i2 = 0; i2 < 16 * ii0; i2 += 1) {
                path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
              }
            }
            if (i0 >= 16 * ii0 + 1) {
              for (register int i2 = 16 * ii0; i2 < -((_PB_N - 1) % 16) + _PB_N - 1; i2 += 1) {
                path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
              }
            }
            for (register int i2 = -((_PB_N - 1) % 16) + _PB_N - 1; i2 < _PB_N; i2 += 1) {
              path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
            }
          }
        }
      }
    }
    if (_PB_N == 17) {
      for (register int i1 = 0; i1 <= 15; i1 += 1) {
        for (register int i2 = 0; i2 <= 15; i2 += 1) {
          path[i1][i2] = ((path[i1][i2] < (path[i1][16] + path[16][i2])) ? path[i1][i2] : (path[i1][16] + path[16][i2]));
        }
      }
    }
  } else if (_PB_N >= 17) {
    for (register int i1 = 0; i1 <= 15; i1 += 1) {
      for (register int i2 = 0; i2 <= 15; i2 += 1) {
        path[i1][i2] = ((path[i1][i2] < (path[i1][0] + path[0][i2])) ? path[i1][i2] : (path[i1][0] + path[0][i2]));
      }
    }
  } else {
    for (register int i0 = 0; i0 < _PB_N; i0 += 1) {
      for (register int i1 = 0; i1 < _PB_N; i1 += 1) {
        for (register int i2 = 0; i2 < _PB_N; i2 += 1) {
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
