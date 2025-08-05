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

/* TC Optimizing Compiler 0.4.2 */
/* ./tc ../examples/polybench/floyd-warshall.scop.c --correction-tiling --free-scheduling --omp-for-codegen --floyd-warshall-tc --debug --inline -b 16 */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#pragma scop
for (register int k = 0; k <= floord(_PB_N - 1, 16); k += 1) {
  if (k == 1) {
    #pragma omp parallel for
    for (register int ii1 = 0; ii1 <= (_PB_N - 1) / 16; ii1 += 1) {
      if (ii1 >= 1) {
        for (register int i1 = 16 * ii1; i1 <= min(_PB_N - 1, 16 * ii1 + 15); i1 += 1) {
          for (register int i2 = 0; i2 <= 15; i2 += 1) {
            path[i1][i2] = ((path[i1][i2] < (path[i1][0] + path[0][i2])) ? path[i1][i2] : (path[i1][0] + path[0][i2]));
          }
        }
        if (16 * ii1 + 16 >= _PB_N) {
          for (register int i0 = 1; i0 <= 2; i0 += 1) {
            if (i0 == 2) {
              path[0][0] = ((path[0][0] < (path[0][2] + path[2][0])) ? path[0][0] : (path[0][2] + path[2][0]));
            } else {
              for (register int i1 = 0; i1 < 16 * ii1; i1 += 1) {
                for (register int i2 = 0; i2 <= 15; i2 += 1) {
                  path[i1][i2] = ((path[i1][i2] < (path[i1][1] + path[1][i2])) ? path[i1][i2] : (path[i1][1] + path[1][i2]));
                }
              }
            }
          }
        }
        if (16 * ii1 + 16 >= _PB_N) {
          for (register int ii2 = 1; ii2 < ii1; ii2 += 1) {
            for (register int i0 = 0; i0 <= 1; i0 += 1) {
              if (i0 == 1) {
                for (register int i1 = 0; i1 < 16 * ii1; i1 += 1) {
                  for (register int i2 = 16 * ii2; i2 <= 16 * ii2 + 15; i2 += 1) {
                    path[i1][i2] = ((path[i1][i2] < (path[i1][1] + path[1][i2])) ? path[i1][i2] : (path[i1][1] + path[1][i2]));
                  }
                }
              } else {
                for (register int i1 = 16 * ii1; i1 < _PB_N; i1 += 1) {
                  for (register int i2 = 16 * ii2; i2 <= 16 * ii2 + 15; i2 += 1) {
                    path[i1][i2] = ((path[i1][i2] < (path[i1][0] + path[0][i2])) ? path[i1][i2] : (path[i1][0] + path[0][i2]));
                  }
                }
              }
            }
          }
          for (register int i0 = 0; i0 <= 15; i0 += 1) {
            if (i0 >= 1) {
              for (register int i1 = 0; i1 < 16 * ii1; i1 += 1) {
                if (i0 >= 2) {
                  for (register int i2 = max(0, -i0 - i1 + 3); i2 < _PB_N; i2 += 1) {
                    path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
                  }
                } else {
                  for (register int i2 = 16 * ii1; i2 < _PB_N; i2 += 1) {
                    path[i1][i2] = ((path[i1][i2] < (path[i1][1] + path[1][i2])) ? path[i1][i2] : (path[i1][1] + path[1][i2]));
                  }
                }
              }
            }
            for (register int i1 = 16 * ii1; i1 < _PB_N; i1 += 1) {
              if (i0 >= 1) {
                for (register int i2 = 0; i2 < 16 * ii1; i2 += 1) {
                  path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
                }
              }
              for (register int i2 = 16 * ii1; i2 < _PB_N; i2 += 1) {
                path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
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
    }
  }
  if (k >= 1) {
    for (register int i1 = 0; i1 <= 15; i1 += 1) {
      for (register int i2 = 0; i2 <= 15; i2 += 1) {
        path[i1][i2] = ((path[i1][i2] < (path[i1][16 * k] + path[16 * k][i2])) ? path[i1][i2] : (path[i1][16 * k] + path[16 * k][i2]));
      }
    }
    if (_PB_N <= 20 && k == 1) {
      #pragma omp parallel for
      for (register int ii2 = 0; ii2 <= min(1, _PB_N - 18); ii2 += 1) {
        if (ii2 == 1) {
          for (register int i0 = 16; i0 < _PB_N; i0 += 1) {
            if (i0 >= 17) {
              for (register int i1 = 0; i1 <= 15; i1 += 1) {
                if (i0 >= 18) {
                  for (register int i2 = 0; i2 <= 15; i2 += 1) {
                    path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
                  }
                }
                for (register int i2 = 16; i2 < _PB_N; i2 += 1) {
                  path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
                }
              }
            }
            for (register int i1 = 16; i1 < _PB_N; i1 += 1) {
              if (i0 >= 17) {
                for (register int i2 = 0; i2 <= 15; i2 += 1) {
                  path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
                }
              }
              for (register int i2 = 16; i2 < _PB_N; i2 += 1) {
                path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
              }
            }
          }
        } else {
          for (register int i0 = 16; i0 <= 17; i0 += 1) {
            if (i0 == 17) {
              for (register int i1 = 0; i1 <= 15; i1 += 1) {
                for (register int i2 = 0; i2 <= 15; i2 += 1) {
                  path[i1][i2] = ((path[i1][i2] < (path[i1][17] + path[17][i2])) ? path[i1][i2] : (path[i1][17] + path[17][i2]));
                }
              }
            } else {
              for (register int i1 = 16; i1 < _PB_N; i1 += 1) {
                for (register int i2 = 0; i2 <= 15; i2 += 1) {
                  path[i1][i2] = ((path[i1][i2] < (path[i1][16] + path[16][i2])) ? path[i1][i2] : (path[i1][16] + path[16][i2]));
                }
              }
            }
          }
        }
      }
    } else {
      if (_PB_N >= 16 * k + 17) {
        for (register int i0 = 16 * k; i0 <= 16 * k + 1; i0 += 1) {
          if (i0 == 16 * k + 1) {
            for (register int i1 = 0; i1 < 16 * k; i1 += 1) {
              for (register int i2 = 0; i2 <= 15; i2 += 1) {
                path[i1][i2] = ((path[i1][i2] < (path[i1][16 * k + 1] + path[16 * k + 1][i2])) ? path[i1][i2] : (path[i1][16 * k + 1] + path[16 * k + 1][i2]));
              }
            }
          } else {
            for (register int i1 = 16 * k; i1 <= 16 * k + 15; i1 += 1) {
              for (register int i2 = 0; i2 <= 15; i2 += 1) {
                path[i1][i2] = ((path[i1][i2] < (path[i1][16 * k] + path[16 * k][i2])) ? path[i1][i2] : (path[i1][16 * k] + path[16 * k][i2]));
              }
            }
          }
        }
        if (k >= 2) {
          for (register int i0 = 16 * k; i0 <= 16 * k + 2; i0 += 1) {
            if (i0 == 16 * k + 1) {
              for (register int i1 = 16 * k; i1 < -((_PB_N - 1) % 16) + _PB_N - 1; i1 += 1) {
                for (register int i2 = 0; i2 <= 15; i2 += 1) {
                  path[i1][i2] = ((path[i1][i2] < (path[i1][16 * k + 1] + path[16 * k + 1][i2])) ? path[i1][i2] : (path[i1][16 * k + 1] + path[16 * k + 1][i2]));
                }
              }
            } else if (i0 == 16 * k + 2) {
              for (register int i1 = 0; i1 < 16 * k; i1 += 1) {
                for (register int i2 = 0; i2 <= 15; i2 += 1) {
                  path[i1][i2] = ((path[i1][i2] < (path[i1][16 * k + 2] + path[16 * k + 2][i2])) ? path[i1][i2] : (path[i1][16 * k + 2] + path[16 * k + 2][i2]));
                }
              }
            } else {
              for (register int i1 = -((_PB_N - 1) % 16) + _PB_N - 1; i1 < _PB_N; i1 += 1) {
                for (register int i2 = 0; i2 <= 15; i2 += 1) {
                  path[i1][i2] = ((path[i1][i2] < (path[i1][16 * k] + path[16 * k][i2])) ? path[i1][i2] : (path[i1][16 * k] + path[16 * k][i2]));
                }
              }
            }
          }
        }
      }
      if (k == 1) {
        #pragma omp parallel for
        for (register int ii2 = 0; ii2 <= min(1, (_PB_N - 1) / 16 - 1); ii2 += 1) {
          if (ii2 == 1) {
            for (register int i0 = 16; i0 <= 19; i0 += 1) {
              if (i0 >= 17) {
                if (i0 >= 18) {
                  for (register int i1 = 0; i1 <= 15; i1 += 1) {
                    if (i0 == 19) {
                      for (register int i2 = 0; i2 <= 15; i2 += 1) {
                        path[i1][i2] = ((path[i1][i2] < (path[i1][19] + path[19][i2])) ? path[i1][i2] : (path[i1][19] + path[19][i2]));
                      }
                    } else {
                      for (register int i2 = 16; i2 <= 31; i2 += 1) {
                        path[i1][i2] = ((path[i1][i2] < (path[i1][18] + path[18][i2])) ? path[i1][i2] : (path[i1][18] + path[18][i2]));
                      }
                    }
                  }
                }
                if (i0 <= 18) {
                  for (register int i1 = 16; i1 <= min(_PB_N - 1, -((_PB_N - 1) % 16) + _PB_N - 16 * i0 + 286); i1 += 1) {
                    if (i0 == 18 && i1 >= 17) {
                      for (register int i2 = 0; i2 <= 15; i2 += 1) {
                        path[i1][i2] = ((path[i1][i2] < (path[i1][18] + path[18][i2])) ? path[i1][i2] : (path[i1][18] + path[18][i2]));
                      }
                    } else if (i0 == 17 && _PB_N >= ((_PB_N - 1) % 16) + i1 + 2) {
                      for (register int i2 = 16; i2 <= 31; i2 += 1) {
                        path[i1][i2] = ((path[i1][i2] < (path[i1][17] + path[17][i2])) ? path[i1][i2] : (path[i1][17] + path[17][i2]));
                      }
                    } else if (i0 == 17 && ((_PB_N - 1) % 16) + i1 + 1 >= _PB_N) {
                      for (register int i2 = 0; i2 <= 15; i2 += 1) {
                        path[i1][i2] = ((path[i1][i2] < (path[i1][17] + path[17][i2])) ? path[i1][i2] : (path[i1][17] + path[17][i2]));
                      }
                    } else {
                      for (register int i2 = 0; i2 <= 16; i2 += 1) {
                        path[16][i2] = ((path[16][i2] < (path[16][18] + path[18][i2])) ? path[16][i2] : (path[16][18] + path[18][i2]));
                      }
                    }
                  }
                }
              } else {
                for (register int i1 = -((_PB_N - 1) % 16) + _PB_N - 1; i1 < _PB_N; i1 += 1) {
                  for (register int i2 = 16; i2 <= 31; i2 += 1) {
                    path[i1][i2] = ((path[i1][i2] < (path[i1][16] + path[16][i2])) ? path[i1][i2] : (path[i1][16] + path[16][i2]));
                  }
                }
              }
            }
          } else {
            for (register int i0 = 16; i0 <= 17; i0 += 1) {
              if (i0 == 17) {
                for (register int i1 = 16; i1 < -((_PB_N - 1) % 16) + _PB_N - 1; i1 += 1) {
                  for (register int i2 = 0; i2 <= 15; i2 += 1) {
                    path[i1][i2] = ((path[i1][i2] < (path[i1][17] + path[17][i2])) ? path[i1][i2] : (path[i1][17] + path[17][i2]));
                  }
                }
                if (_PB_N <= 32) {
                  for (register int i1 = 0; i1 <= 15; i1 += 1) {
                    for (register int i2 = 0; i2 <= 15; i2 += 1) {
                      path[i1][i2] = ((path[i1][i2] < (path[i1][17] + path[17][i2])) ? path[i1][i2] : (path[i1][17] + path[17][i2]));
                    }
                  }
                }
              } else {
                for (register int i1 = -((_PB_N - 1) % 16) + _PB_N - 1; i1 < _PB_N; i1 += 1) {
                  for (register int i2 = 0; i2 <= 15; i2 += 1) {
                    path[i1][i2] = ((path[i1][i2] < (path[i1][16] + path[16][i2])) ? path[i1][i2] : (path[i1][16] + path[16][i2]));
                  }
                }
              }
            }
            if (_PB_N >= 33) {
              for (register int i1 = 0; i1 <= 15; i1 += 1) {
                for (register int i2 = 0; i2 <= 15; i2 += 1) {
                  path[i1][i2] = ((path[i1][i2] < (path[i1][18] + path[18][i2])) ? path[i1][i2] : (path[i1][18] + path[18][i2]));
                }
              }
            }
          }
        }
        for (register int i0 = 16; i0 <= min(31, _PB_N - 1); i0 += 1) {
          if (_PB_N <= 32 && i0 == 17) {
            for (register int i1 = 0; i1 <= 15; i1 += 1) {
              for (register int i2 = 16; i2 < _PB_N; i2 += 1) {
                path[i1][i2] = ((path[i1][i2] < (path[i1][17] + path[17][i2])) ? path[i1][i2] : (path[i1][17] + path[17][i2]));
              }
            }
          }
          if (i0 >= 17) {
            for (register int i1 = max(0, -16 * i0 + 288); i1 < -((_PB_N - 1) % 16) + _PB_N - 1; i1 += 1) {
              if (_PB_N >= 33 && i1 <= 15) {
                if (i0 == 19) {
                  for (register int i2 = 16; i2 < -((_PB_N - 1) % 16) + _PB_N - 1; i2 += 1) {
                    path[i1][i2] = ((path[i1][i2] < (path[i1][19] + path[19][i2])) ? path[i1][i2] : (path[i1][19] + path[19][i2]));
                  }
                } else if (i0 >= 20) {
                  for (register int i2 = 0; i2 < -((_PB_N - 1) % 16) + _PB_N - 1; i2 += 1) {
                    path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
                  }
                }
                for (register int i2 = -((_PB_N - 1) % 16) + _PB_N - 1; i2 < _PB_N; i2 += 1) {
                  path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
                }
              } else if (_PB_N <= 32) {
                for (register int i2 = 0; i2 < _PB_N; i2 += 1) {
                  path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
                }
              } else if (i0 >= 18) {
                if (i0 == 19) {
                  for (register int i2 = 0; i2 <= 15; i2 += 1) {
                    path[i1][i2] = ((path[i1][i2] < (path[i1][19] + path[19][i2])) ? path[i1][i2] : (path[i1][19] + path[19][i2]));
                  }
                  if (i1 == 16) {
                    path[16][16] = ((path[16][16] < (path[16][19] + path[19][16])) ? path[16][16] : (path[16][19] + path[19][16]));
                  }
                } else if (i0 >= 20) {
                  for (register int i2 = 0; i2 <= 15; i2 += 1) {
                    path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
                  }
                  if (i1 == 16) {
                    path[16][16] = ((path[16][16] < (path[16][i0] + path[i0][16])) ? path[16][16] : (path[16][i0] + path[i0][16]));
                  }
                }
                for (register int i2 = max(16, -i1 + 33); i2 < _PB_N; i2 += 1) {
                  path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
                }
              } else {
                for (register int i2 = -((_PB_N - 1) % 16) + _PB_N - 1; i2 < _PB_N; i2 += 1) {
                  path[i1][i2] = ((path[i1][i2] < (path[i1][17] + path[17][i2])) ? path[i1][i2] : (path[i1][17] + path[17][i2]));
                }
              }
            }
          }
          for (register int i1 = -((_PB_N - 1) % 16) + _PB_N - 1; i1 < _PB_N; i1 += 1) {
            if (_PB_N <= 32 && i0 == 17) {
              for (register int i2 = 0; i2 <= 15; i2 += 1) {
                path[i1][i2] = ((path[i1][i2] < (path[i1][17] + path[17][i2])) ? path[i1][i2] : (path[i1][17] + path[17][i2]));
              }
            }
            if (i0 == 17) {
              for (register int i2 = 16; i2 < -((_PB_N - 1) % 16) + _PB_N - 1; i2 += 1) {
                path[i1][i2] = ((path[i1][i2] < (path[i1][17] + path[17][i2])) ? path[i1][i2] : (path[i1][17] + path[17][i2]));
              }
            } else if (i0 >= 18) {
              for (register int i2 = 0; i2 < -((_PB_N - 1) % 16) + _PB_N - 1; i2 += 1) {
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
  } else if (_PB_N >= 21) {
    for (register int i1 = 0; i1 <= 15; i1 += 1) {
      for (register int i2 = 0; i2 <= 15; i2 += 1) {
        path[i1][i2] = ((path[i1][i2] < (path[i1][0] + path[0][i2])) ? path[i1][i2] : (path[i1][0] + path[0][i2]));
      }
    }
  } else if (_PB_N <= 16) {
    for (register int i0 = 0; i0 < _PB_N; i0 += 1) {
      for (register int i1 = 0; i1 < _PB_N; i1 += 1) {
        for (register int i2 = 0; i2 < _PB_N; i2 += 1) {
          path[i1][i2] = ((path[i1][i2] < (path[i1][i0] + path[i0][i2])) ? path[i1][i2] : (path[i1][i0] + path[i0][i2]));
        }
      }
    }
  } else {
    for (register int i1 = 0; i1 <= 15; i1 += 1) {
      for (register int i2 = 0; i2 <= 15; i2 += 1) {
        path[i1][i2] = ((path[i1][i2] < (path[i1][0] + path[0][i2])) ? path[i1][i2] : (path[i1][0] + path[0][i2]));
      }
    }
  }
  if (k == 1) {
    #pragma omp parallel for
    for (register int ii0 = 2; ii0 <= (_PB_N - 2) / 16; ii0 += 1) {
      if (_PB_N >= 16 * ii0 + 5) {
        if (_PB_N >= 16 * ii0 + 17) {
          for (register int i0 = 16 * ii0; i0 <= 16 * ii0 + 3; i0 += 1) {
            if (i0 >= 16 * ii0 + 2) {
              for (register int i1 = 0; i1 < 16 * ii0; i1 += 1) {
                if (i0 == 16 * ii0 + 3) {
                  for (register int i2 = 0; i2 < 16 * ii0; i2 += 1) {
                    path[i1][i2] = ((path[i1][i2] < (path[i1][16 * ii0 + 3] + path[16 * ii0 + 3][i2])) ? path[i1][i2] : (path[i1][16 * ii0 + 3] + path[16 * ii0 + 3][i2]));
                  }
                } else {
                  for (register int i2 = 16 * ii0; i2 <= 16 * ii0 + 15; i2 += 1) {
                    path[i1][i2] = ((path[i1][i2] < (path[i1][16 * ii0 + 2] + path[16 * ii0 + 2][i2])) ? path[i1][i2] : (path[i1][16 * ii0 + 2] + path[16 * ii0 + 2][i2]));
                  }
                }
              }
              if (i0 == 16 * ii0 + 2) {
                for (register int i1 = 16 * ii0; i1 < -((_PB_N - 1) % 16) + _PB_N - 1; i1 += 1) {
                  if (i1 >= 16 * ii0 + 1) {
                    for (register int i2 = 0; i2 < 16 * ii0; i2 += 1) {
                      path[i1][i2] = ((path[i1][i2] < (path[i1][16 * ii0 + 2] + path[16 * ii0 + 2][i2])) ? path[i1][i2] : (path[i1][16 * ii0 + 2] + path[16 * ii0 + 2][i2]));
                    }
                  } else {
                    for (register int i2 = 0; i2 <= 16 * ii0; i2 += 1) {
                      path[16 * ii0][i2] = ((path[16 * ii0][i2] < (path[16 * ii0][16 * ii0 + 2] + path[16 * ii0 + 2][i2])) ? path[16 * ii0][i2] : (path[16 * ii0][16 * ii0 + 2] + path[16 * ii0 + 2][i2]));
                    }
                  }
                }
              }
            } else if (i0 == 16 * ii0 + 1) {
              for (register int i1 = 16 * ii0; i1 < _PB_N; i1 += 1) {
                if (_PB_N >= ((_PB_N - 1) % 16) + i1 + 2) {
                  for (register int i2 = 16 * ii0; i2 <= 16 * ii0 + 15; i2 += 1) {
                    path[i1][i2] = ((path[i1][i2] < (path[i1][16 * ii0 + 1] + path[16 * ii0 + 1][i2])) ? path[i1][i2] : (path[i1][16 * ii0 + 1] + path[16 * ii0 + 1][i2]));
                  }
                } else {
                  for (register int i2 = 0; i2 < 16 * ii0; i2 += 1) {
                    path[i1][i2] = ((path[i1][i2] < (path[i1][16 * ii0 + 1] + path[16 * ii0 + 1][i2])) ? path[i1][i2] : (path[i1][16 * ii0 + 1] + path[16 * ii0 + 1][i2]));
                  }
                }
              }
            } else {
              for (register int i1 = -((_PB_N - 1) % 16) + _PB_N - 1; i1 < _PB_N; i1 += 1) {
                for (register int i2 = 16 * ii0; i2 <= 16 * ii0 + 15; i2 += 1) {
                  path[i1][i2] = ((path[i1][i2] < (path[i1][16 * ii0] + path[16 * ii0][i2])) ? path[i1][i2] : (path[i1][16 * ii0] + path[16 * ii0][i2]));
                }
              }
            }
          }
        } else {
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
      } else {
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
        if (_PB_N >= 16 * ii0 + 3) {
          for (register int i0 = 16 * ii0; i0 < _PB_N; i0 += 1) {
            if (i0 >= 16 * ii0 + 1) {
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
            for (register int i1 = 16 * ii0; i1 < _PB_N; i1 += 1) {
              if (i0 >= 16 * ii0 + 1) {
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
