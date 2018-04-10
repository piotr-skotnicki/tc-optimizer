/**
 * This version is stamped on Apr. 14, 2015
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* 2mm.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "2mm.h"


/* Array initialization. */
static
void init_array(int ni, int nj, int nk, int nl,
		DATA_TYPE *alpha,
		DATA_TYPE *beta,
		DATA_TYPE POLYBENCH_2D(A,NI,NK,ni,nk),
		DATA_TYPE POLYBENCH_2D(B,NK,NJ,nk,nj),
		DATA_TYPE POLYBENCH_2D(C,NJ,NL,nj,nl),
		DATA_TYPE POLYBENCH_2D(D,NI,NL,ni,nl))
{
  int i, j;

  *alpha = 1.5;
  *beta = 1.2;
  for (i = 0; i < ni; i++)
    for (j = 0; j < nk; j++)
      A[i][j] = (DATA_TYPE) (i*j % ni) / ni;
  for (i = 0; i < nk; i++)
    for (j = 0; j < nj; j++)
      B[i][j] = (DATA_TYPE) (i*(j+1) % nj) / nj;
  for (i = 0; i < nj; i++)
    for (j = 0; j < nl; j++)
      C[i][j] = (DATA_TYPE) (i*(j+3) % nl) / nl;
  for (i = 0; i < ni; i++)
    for (j = 0; j < nl; j++)
      D[i][j] = (DATA_TYPE) (i*(j+2) % nk) / nk;
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int ni, int nl,
		 DATA_TYPE POLYBENCH_2D(D,NI,NL,ni,nl))
{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("D");
  for (i = 0; i < ni; i++)
    for (j = 0; j < nl; j++) {
	if ((i * ni + j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
	fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, D[i][j]);
    }
  POLYBENCH_DUMP_END("D");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_2mm(int ni, int nj, int nk, int nl,
		DATA_TYPE alpha,
		DATA_TYPE beta,
		DATA_TYPE POLYBENCH_2D(tmp,NI,NJ,ni,nj),
		DATA_TYPE POLYBENCH_2D(A,NI,NK,ni,nk),
		DATA_TYPE POLYBENCH_2D(B,NK,NJ,nk,nj),
		DATA_TYPE POLYBENCH_2D(C,NJ,NL,nj,nl),
		DATA_TYPE POLYBENCH_2D(D,NI,NL,ni,nl))
{
  int i, j, k;

/* TC Optimizing Compiler 0.2.26 */
/* ./tc ../examples/polybench/2mm.scop.c --correction-tiling --lex-scheduling --serial-codegen -b 32 --debug */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#pragma scop
{
  if (_PB_NL >= 1 && _PB_NJ >= _PB_NL && _PB_NJ <= 31) {
    for (int ii1 = 0; ii1 <= floord(_PB_NI - 1, 32); ii1 += 1) {
      for (int i1 = 32 * ii1; i1 <= min(_PB_NI - 1, 32 * ii1 + 31); i1 += 1) {
        for (int i2 = 0; i2 < _PB_NJ; i2 += 1) {
          tmp[i1][i2] = SCALAR_VAL(0.0);
        }
      }
      if (_PB_NJ >= _PB_NK + 1 && _PB_NI >= 32 * ii1 + 32) {
        for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
          for (int i2 = 0; i2 < _PB_NJ; i2 += 1) {
            for (int i4 = 0; i4 < _PB_NK; i4 += 1) {
              tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
            }
          }
        }
      } else if (32 * ii1 + 31 >= _PB_NI && _PB_NK % 32 >= 1) {
        for (int ii4 = 0; ii4 < _PB_NK / 32; ii4 += 1) {
          for (int i1 = 32 * ii1; i1 < _PB_NI; i1 += 1) {
            for (int i2 = 0; i2 < _PB_NJ; i2 += 1) {
              for (int i4 = 32 * ii4; i4 <= 32 * ii4 + 31; i4 += 1) {
                tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
              }
            }
          }
        }
        for (int i1 = 32 * ii1; i1 < _PB_NI; i1 += 1) {
          for (int i2 = 0; i2 < _PB_NJ; i2 += 1) {
            for (int i4 = -((_PB_NK + 31) % 32) + _PB_NK - 1; i4 < _PB_NK; i4 += 1) {
              tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
            }
          }
        }
      } else if (_PB_NK >= _PB_NJ && _PB_NI >= 32 * ii1 + 32 && _PB_NK % 32 >= 1) {
        for (int ii4 = 0; ii4 < _PB_NK / 32; ii4 += 1) {
          for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
            for (int i2 = 0; i2 < _PB_NJ; i2 += 1) {
              for (int i4 = 32 * ii4; i4 <= 32 * ii4 + 31; i4 += 1) {
                tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
              }
            }
          }
        }
        for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
          for (int i2 = 0; i2 < _PB_NJ; i2 += 1) {
            for (int i4 = -((_PB_NK - 1) % 32) + _PB_NK - 1; i4 < _PB_NK; i4 += 1) {
              tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
            }
          }
        }
      }
      if (_PB_NI >= 32 * ii1 + 32 && _PB_NK % 32 == 0) {
        for (int ii4 = 0; ii4 < _PB_NK / 32; ii4 += 1) {
          for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
            for (int i2 = 0; i2 < _PB_NJ; i2 += 1) {
              for (int i4 = 32 * ii4; i4 <= 32 * ii4 + 31; i4 += 1) {
                tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
              }
            }
          }
        }
      } else if (32 * ii1 + 31 >= _PB_NI && _PB_NK % 32 == 0) {
        for (int ii4 = 0; ii4 < _PB_NK / 32; ii4 += 1) {
          for (int i1 = 32 * ii1; i1 < _PB_NI; i1 += 1) {
            for (int i2 = 0; i2 < _PB_NJ; i2 += 1) {
              for (int i4 = 32 * ii4; i4 <= 32 * ii4 + 31; i4 += 1) {
                tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
              }
            }
          }
        }
      }
    }
  } else if (_PB_NL >= 1 && _PB_NJ >= 32 && _PB_NJ >= _PB_NL && _PB_NJ + _PB_NK >= _PB_NL + 1) {
    for (int ii1 = 0; ii1 <= floord(_PB_NI - 1, 32); ii1 += 1) {
      for (int ii2 = 0; ii2 < _PB_NJ / 32; ii2 += 1) {
        for (int i1 = 32 * ii1; i1 <= min(_PB_NI - 1, 32 * ii1 + 31); i1 += 1) {
          for (int i2 = 32 * ii2; i2 <= 32 * ii2 + 31; i2 += 1) {
            tmp[i1][i2] = SCALAR_VAL(0.0);
          }
        }
        if (_PB_NK >= 1) {
          if (_PB_NJ >= _PB_NK + 1 && 32 * ii1 + 31 >= _PB_NI) {
            for (int ii4 = 0; ii4 <= min(_PB_NJ / 32 - 1, floord(_PB_NK - 1, 32)); ii4 += 1) {
              for (int i1 = 32 * ii1; i1 < _PB_NI; i1 += 1) {
                for (int i2 = 32 * ii2; i2 <= 32 * ii2 + 31; i2 += 1) {
                  for (int i4 = 32 * ii4; i4 <= min(_PB_NK - 1, 32 * ii4 + 31); i4 += 1) {
                    tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
                  }
                }
              }
            }
          } else if (_PB_NK >= _PB_NJ) {
            for (int ii4 = 0; ii4 < _PB_NK / 32; ii4 += 1) {
              for (int i1 = 32 * ii1; i1 <= min(_PB_NI - 1, 32 * ii1 + 31); i1 += 1) {
                for (int i2 = 32 * ii2; i2 <= 32 * ii2 + 31; i2 += 1) {
                  for (int i4 = 32 * ii4; i4 <= 32 * ii4 + 31; i4 += 1) {
                    tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
                  }
                }
              }
            }
          }
          if (32 * ii1 + 31 >= _PB_NI && _PB_NK % 32 >= 1 && _PB_NK + 31 >= (_PB_NK % 32) + _PB_NJ) {
            for (int i1 = 32 * ii1; i1 < _PB_NI; i1 += 1) {
              for (int i2 = 32 * ii2; i2 <= 32 * ii2 + 31; i2 += 1) {
                for (int i4 = -((_PB_NK - 1) % 32) + _PB_NK - 1; i4 < _PB_NK; i4 += 1) {
                  tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
                }
              }
            }
          } else if (_PB_NK >= _PB_NJ && _PB_NI >= 32 * ii1 + 32 && _PB_NK % 32 >= 1) {
            for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
              for (int i2 = 32 * ii2; i2 <= 32 * ii2 + 31; i2 += 1) {
                for (int i4 = -((_PB_NK - 1) % 32) + _PB_NK - 1; i4 < _PB_NK; i4 += 1) {
                  tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
                }
              }
            }
          }
          if (_PB_NJ >= _PB_NK + 1 && _PB_NI >= 32 * ii1 + 32) {
            for (int ii4 = 0; ii4 <= floord(_PB_NK - 1, 32); ii4 += 1) {
              for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
                for (int i2 = 32 * ii2; i2 <= 32 * ii2 + 31; i2 += 1) {
                  for (int i4 = 32 * ii4; i4 <= min(_PB_NK - 1, 32 * ii4 + 31); i4 += 1) {
                    tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
                  }
                }
              }
            }
          }
        }
      }
      if (_PB_NJ % 32 >= 1) {
        for (int i1 = 32 * ii1; i1 <= min(_PB_NI - 1, 32 * ii1 + 31); i1 += 1) {
          for (int i2 = -((_PB_NJ - 1) % 32) + _PB_NJ - 1; i2 < _PB_NJ; i2 += 1) {
            tmp[i1][i2] = SCALAR_VAL(0.0);
          }
        }
        if (_PB_NK >= 1) {
          if (_PB_NK >= _PB_NJ) {
            for (int ii4 = 0; ii4 < _PB_NK / 32; ii4 += 1) {
              for (int i1 = 32 * ii1; i1 <= min(_PB_NI - 1, 32 * ii1 + 31); i1 += 1) {
                for (int i2 = -((_PB_NJ - 1) % 32) + _PB_NJ - 1; i2 < _PB_NJ; i2 += 1) {
                  for (int i4 = 32 * ii4; i4 <= 32 * ii4 + 31; i4 += 1) {
                    tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
                  }
                }
              }
            }
          } else if (32 * ii1 + 31 >= _PB_NI) {
            for (int ii4 = 0; ii4 <= min((_PB_NJ - 1) / 32 - 1, floord(_PB_NK - 1, 32)); ii4 += 1) {
              for (int i1 = 32 * ii1; i1 < _PB_NI; i1 += 1) {
                for (int i2 = -((_PB_NJ - 1) % 32) + _PB_NJ - 1; i2 < _PB_NJ; i2 += 1) {
                  for (int i4 = 32 * ii4; i4 <= min(_PB_NK - 1, 32 * ii4 + 31); i4 += 1) {
                    tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
                  }
                }
              }
            }
          }
          if (32 * ii1 + 31 >= _PB_NI && _PB_NK % 32 >= 1 && _PB_NK + 31 >= (_PB_NK % 32) + _PB_NJ) {
            for (int i1 = 32 * ii1; i1 < _PB_NI; i1 += 1) {
              for (int i2 = -((_PB_NJ - 1) % 32) + _PB_NJ - 1; i2 < _PB_NJ; i2 += 1) {
                for (int i4 = -((_PB_NK - 1) % 32) + _PB_NK - 1; i4 < _PB_NK; i4 += 1) {
                  tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
                }
              }
            }
          } else if (_PB_NK >= _PB_NJ && _PB_NI >= 32 * ii1 + 32 && _PB_NK % 32 >= 1) {
            for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
              for (int i2 = -((_PB_NJ - 1) % 32) + _PB_NJ - 1; i2 < _PB_NJ; i2 += 1) {
                for (int i4 = -((_PB_NK - 1) % 32) + _PB_NK - 1; i4 < _PB_NK; i4 += 1) {
                  tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
                }
              }
            }
          }
          if (_PB_NJ >= _PB_NK + 1 && _PB_NI >= 32 * ii1 + 32) {
            for (int ii4 = 0; ii4 <= floord(_PB_NK - 1, 32); ii4 += 1) {
              for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
                for (int i2 = -((_PB_NJ - 1) % 32) + _PB_NJ - 1; i2 < _PB_NJ; i2 += 1) {
                  for (int i4 = 32 * ii4; i4 <= min(_PB_NK - 1, 32 * ii4 + 31); i4 += 1) {
                    tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
                  }
                }
              }
            }
          }
        }
      }
    }
  } else if (_PB_NK >= 1 && (_PB_NJ % 32) + _PB_NL >= _PB_NJ + 32) {
    for (int ii1 = 0; ii1 <= floord(_PB_NI - 1, 32); ii1 += 1) {
      for (int ii2 = 0; ii2 <= floord(_PB_NJ - 1, 32); ii2 += 1) {
        for (int i1 = 32 * ii1; i1 <= min(_PB_NI - 1, 32 * ii1 + 31); i1 += 1) {
          for (int i2 = 32 * ii2; i2 <= min(_PB_NJ - 1, 32 * ii2 + 31); i2 += 1) {
            tmp[i1][i2] = SCALAR_VAL(0.0);
          }
        }
        if (_PB_NJ >= _PB_NK + 31 && _PB_NI >= 32 * ii1 + 32) {
          for (int ii4 = 0; ii4 <= floord(_PB_NK - 1, 32); ii4 += 1) {
            for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
              for (int i2 = 32 * ii2; i2 <= min(_PB_NJ - 1, 32 * ii2 + 31); i2 += 1) {
                for (int i4 = 32 * ii4; i4 <= min(_PB_NK - 1, 32 * ii4 + 31); i4 += 1) {
                  tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
                }
              }
            }
          }
        } else if (_PB_NK >= _PB_NJ && _PB_NI >= 32 * ii1 + 32) {
          for (int ii4 = 0; ii4 <= (_PB_NK - 1) / 32; ii4 += 1) {
            for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
              for (int i2 = 32 * ii2; i2 <= min(_PB_NJ - 1, 32 * ii2 + 31); i2 += 1) {
                for (int i4 = 32 * ii4; i4 <= min(_PB_NK - 1, 32 * ii4 + 31); i4 += 1) {
                  tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
                }
              }
            }
          }
        } else if (_PB_NJ >= _PB_NK + 1 && _PB_NK + 30 >= _PB_NJ && _PB_NI >= 32 * ii1 + 32) {
          for (int ii4 = 0; ii4 < _PB_NJ / 32; ii4 += 1) {
            for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
              for (int i2 = 32 * ii2; i2 <= min(_PB_NJ - 1, 32 * ii2 + 31); i2 += 1) {
                for (int i4 = 32 * ii4; i4 <= min(_PB_NK - 1, 32 * ii4 + 31); i4 += 1) {
                  tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
                }
              }
            }
          }
          for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
            for (int i2 = 32 * ii2; i2 <= min(_PB_NJ - 1, 32 * ii2 + 31); i2 += 1) {
              for (int i4 = -(_PB_NJ % 32) + _PB_NJ; i4 < _PB_NK; i4 += 1) {
                tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
              }
            }
          }
        } else if ((_PB_NK % 32) + _PB_NJ >= _PB_NK + 32) {
          for (int ii4 = 0; ii4 <= floord(_PB_NK - 1, 32); ii4 += 1) {
            for (int i1 = 32 * ii1; i1 < _PB_NI; i1 += 1) {
              for (int i2 = 32 * ii2; i2 <= min(_PB_NJ - 1, 32 * ii2 + 31); i2 += 1) {
                for (int i4 = 32 * ii4; i4 <= min(_PB_NK - 1, 32 * ii4 + 31); i4 += 1) {
                  tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
                }
              }
            }
          }
        } else if (_PB_NK % 32 >= 1 && _PB_NK + 31 >= (_PB_NK % 32) + _PB_NJ) {
          if (_PB_NJ >= _PB_NK + 1) {
            for (int ii4 = 0; ii4 < _PB_NJ / 32; ii4 += 1) {
              for (int i1 = 32 * ii1; i1 < _PB_NI; i1 += 1) {
                for (int i2 = 32 * ii2; i2 <= min(_PB_NJ - 1, 32 * ii2 + 31); i2 += 1) {
                  for (int i4 = 32 * ii4; i4 <= 32 * ii4 + 31; i4 += 1) {
                    tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
                  }
                }
              }
            }
          } else {
            for (int ii4 = 0; ii4 < _PB_NK / 32; ii4 += 1) {
              for (int i1 = 32 * ii1; i1 < _PB_NI; i1 += 1) {
                for (int i2 = 32 * ii2; i2 <= min(_PB_NJ - 1, 32 * ii2 + 31); i2 += 1) {
                  for (int i4 = 32 * ii4; i4 <= 32 * ii4 + 31; i4 += 1) {
                    tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
                  }
                }
              }
            }
          }
          for (int i1 = 32 * ii1; i1 < _PB_NI; i1 += 1) {
            for (int i2 = 32 * ii2; i2 <= min(_PB_NJ - 1, 32 * ii2 + 31); i2 += 1) {
              for (int i4 = -((_PB_NK + 31) % 32) + _PB_NK - 1; i4 < _PB_NK; i4 += 1) {
                tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
              }
            }
          }
        } else if (_PB_NJ >= _PB_NK + 1 && _PB_NK % 32 == 0) {
          for (int ii4 = 0; ii4 < _PB_NK / 32; ii4 += 1) {
            for (int i1 = 32 * ii1; i1 < _PB_NI; i1 += 1) {
              for (int i2 = 32 * ii2; i2 <= min(_PB_NJ - 1, 32 * ii2 + 31); i2 += 1) {
                for (int i4 = 32 * ii4; i4 <= 32 * ii4 + 31; i4 += 1) {
                  tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
                }
              }
            }
          }
        } else {
          for (int ii4 = 0; ii4 < _PB_NK / 32; ii4 += 1) {
            for (int i1 = 32 * ii1; i1 < _PB_NI; i1 += 1) {
              for (int i2 = 32 * ii2; i2 <= min(_PB_NJ - 1, 32 * ii2 + 31); i2 += 1) {
                for (int i4 = 32 * ii4; i4 <= 32 * ii4 + 31; i4 += 1) {
                  tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
                }
              }
            }
          }
        }
      }
    }
  } else if (_PB_NL >= _PB_NJ + 1 && _PB_NJ % 32 >= 1 && _PB_NJ + 31 >= (_PB_NJ % 32) + _PB_NL) {
    for (int ii1 = 0; ii1 <= floord(_PB_NI - 1, 32); ii1 += 1) {
      if (_PB_NK + 31 >= _PB_NJ && 32 * ii1 + 31 >= _PB_NI && _PB_NK % 32 == 0) {
        for (int ii2 = 0; ii2 < _PB_NL / 32; ii2 += 1) {
          for (int i1 = 32 * ii1; i1 < _PB_NI; i1 += 1) {
            for (int i2 = 32 * ii2; i2 <= 32 * ii2 + 31; i2 += 1) {
              tmp[i1][i2] = SCALAR_VAL(0.0);
            }
          }
          if (_PB_NJ >= _PB_NK + 1) {
            for (int ii4 = 0; ii4 < _PB_NK / 32; ii4 += 1) {
              for (int i1 = 32 * ii1; i1 < _PB_NI; i1 += 1) {
                for (int i2 = 32 * ii2; i2 <= 32 * ii2 + 31; i2 += 1) {
                  for (int i4 = 32 * ii4; i4 <= 32 * ii4 + 31; i4 += 1) {
                    tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
                  }
                }
              }
            }
          } else {
            for (int ii4 = 0; ii4 < _PB_NK / 32; ii4 += 1) {
              for (int i1 = 32 * ii1; i1 < _PB_NI; i1 += 1) {
                for (int i2 = 32 * ii2; i2 <= 32 * ii2 + 31; i2 += 1) {
                  for (int i4 = 32 * ii4; i4 <= 32 * ii4 + 31; i4 += 1) {
                    tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
                  }
                }
              }
            }
          }
        }
      } else if (32 * ii1 + 31 >= _PB_NI && (_PB_NK % 32) + _PB_NJ >= _PB_NK + 32) {
        for (int ii2 = 0; ii2 < _PB_NL / 32; ii2 += 1) {
          if (_PB_NK >= 1) {
            for (int i1 = 32 * ii1; i1 < _PB_NI; i1 += 1) {
              for (int i2 = 32 * ii2; i2 <= 32 * ii2 + 31; i2 += 1) {
                tmp[i1][i2] = SCALAR_VAL(0.0);
              }
            }
          } else {
            for (int i1 = 32 * ii1; i1 < _PB_NI; i1 += 1) {
              for (int i2 = 32 * ii2; i2 <= 32 * ii2 + 31; i2 += 1) {
                tmp[i1][i2] = SCALAR_VAL(0.0);
              }
            }
          }
          for (int ii4 = 0; ii4 <= floord(_PB_NK - 1, 32); ii4 += 1) {
            for (int i1 = 32 * ii1; i1 < _PB_NI; i1 += 1) {
              for (int i2 = 32 * ii2; i2 <= 32 * ii2 + 31; i2 += 1) {
                for (int i4 = 32 * ii4; i4 <= min(_PB_NK - 1, 32 * ii4 + 31); i4 += 1) {
                  tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
                }
              }
            }
          }
        }
      } else if (_PB_NI >= 32 * ii1 + 32) {
        for (int ii2 = 0; ii2 < _PB_NL / 32; ii2 += 1) {
          for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
            for (int i2 = 32 * ii2; i2 <= 32 * ii2 + 31; i2 += 1) {
              tmp[i1][i2] = SCALAR_VAL(0.0);
            }
          }
          if (_PB_NL >= _PB_NK + 31) {
            for (int ii4 = 0; ii4 <= floord(_PB_NK - 1, 32); ii4 += 1) {
              for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
                for (int i2 = 32 * ii2; i2 <= 32 * ii2 + 31; i2 += 1) {
                  for (int i4 = 32 * ii4; i4 <= min(_PB_NK - 1, 32 * ii4 + 31); i4 += 1) {
                    tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
                  }
                }
              }
            }
          } else if (_PB_NK >= _PB_NJ) {
            for (int ii4 = 0; ii4 <= (_PB_NK - 1) / 32; ii4 += 1) {
              for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
                for (int i2 = 32 * ii2; i2 <= 32 * ii2 + 31; i2 += 1) {
                  for (int i4 = 32 * ii4; i4 <= min(_PB_NK - 1, 32 * ii4 + 31); i4 += 1) {
                    tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
                  }
                }
              }
            }
          } else {
            for (int ii4 = 0; ii4 < _PB_NJ / 32; ii4 += 1) {
              for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
                for (int i2 = 32 * ii2; i2 <= 32 * ii2 + 31; i2 += 1) {
                  for (int i4 = 32 * ii4; i4 <= min(_PB_NK - 1, 32 * ii4 + 31); i4 += 1) {
                    tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
                  }
                }
              }
            }
            for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
              for (int i2 = 32 * ii2; i2 <= 32 * ii2 + 31; i2 += 1) {
                for (int i4 = ((31 * _PB_NL + 31) % 32) + _PB_NL - 31; i4 < _PB_NK; i4 += 1) {
                  tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
                }
              }
            }
          }
        }
      } else {
        for (int ii2 = 0; ii2 < _PB_NL / 32; ii2 += 1) {
          for (int i1 = 32 * ii1; i1 < _PB_NI; i1 += 1) {
            for (int i2 = 32 * ii2; i2 <= 32 * ii2 + 31; i2 += 1) {
              tmp[i1][i2] = SCALAR_VAL(0.0);
            }
          }
          if (_PB_NJ >= _PB_NK + 1) {
            for (int ii4 = 0; ii4 < _PB_NJ / 32; ii4 += 1) {
              for (int i1 = 32 * ii1; i1 < _PB_NI; i1 += 1) {
                for (int i2 = 32 * ii2; i2 <= 32 * ii2 + 31; i2 += 1) {
                  for (int i4 = 32 * ii4; i4 <= 32 * ii4 + 31; i4 += 1) {
                    tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
                  }
                }
              }
            }
          } else {
            for (int ii4 = 0; ii4 < _PB_NK / 32; ii4 += 1) {
              for (int i1 = 32 * ii1; i1 < _PB_NI; i1 += 1) {
                for (int i2 = 32 * ii2; i2 <= 32 * ii2 + 31; i2 += 1) {
                  for (int i4 = 32 * ii4; i4 <= 32 * ii4 + 31; i4 += 1) {
                    tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
                  }
                }
              }
            }
          }
          for (int i1 = 32 * ii1; i1 < _PB_NI; i1 += 1) {
            for (int i2 = 32 * ii2; i2 <= 32 * ii2 + 31; i2 += 1) {
              for (int i4 = -((_PB_NK + 31) % 32) + _PB_NK - 1; i4 < _PB_NK; i4 += 1) {
                tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
              }
            }
          }
        }
      }
      for (int i1 = 32 * ii1; i1 <= min(_PB_NI - 1, 32 * ii1 + 31); i1 += 1) {
        if (_PB_NJ >= _PB_NK + 1) {
          for (int i2 = ((31 * _PB_NL + 31) % 32) + _PB_NL - 31; i2 < _PB_NJ; i2 += 1) {
            tmp[i1][i2] = SCALAR_VAL(0.0);
          }
        } else {
          for (int i2 = ((31 * _PB_NL + 31) % 32) + _PB_NL - 31; i2 < _PB_NJ; i2 += 1) {
            tmp[i1][i2] = SCALAR_VAL(0.0);
          }
        }
      }
      if (_PB_NI >= 32 * ii1 + 32 && (_PB_NK % 32) + _PB_NJ >= _PB_NK + 32) {
        for (int ii4 = 0; ii4 <= floord(_PB_NK - 1, 32); ii4 += 1) {
          for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
            for (int i2 = ((31 * _PB_NL + 31) % 32) + _PB_NL - 31; i2 < _PB_NJ; i2 += 1) {
              for (int i4 = 32 * ii4; i4 <= min(_PB_NK - 1, 32 * ii4 + 31); i4 += 1) {
                tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
              }
            }
          }
        }
      } else if (_PB_NI >= 32 * ii1 + 32 && _PB_NK % 32 >= 1 && _PB_NK + 31 >= (_PB_NK % 32) + _PB_NJ) {
        if (_PB_NK >= _PB_NJ) {
          for (int ii4 = 0; ii4 < _PB_NK / 32; ii4 += 1) {
            for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
              for (int i2 = ((31 * _PB_NL + 31) % 32) + _PB_NL - 31; i2 < _PB_NJ; i2 += 1) {
                for (int i4 = 32 * ii4; i4 <= 32 * ii4 + 31; i4 += 1) {
                  tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
                }
              }
            }
          }
        } else {
          for (int ii4 = 0; ii4 < _PB_NL - (31 * _PB_NL + 31) / 32; ii4 += 1) {
            for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
              for (int i2 = ((31 * _PB_NL + 31) % 32) + _PB_NL - 31; i2 < _PB_NJ; i2 += 1) {
                for (int i4 = 32 * ii4; i4 <= 32 * ii4 + 31; i4 += 1) {
                  tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
                }
              }
            }
          }
        }
        for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
          for (int i2 = ((31 * _PB_NL + 31) % 32) + _PB_NL - 31; i2 < _PB_NJ; i2 += 1) {
            for (int i4 = -((_PB_NK + 31) % 32) + _PB_NK - 1; i4 < _PB_NK; i4 += 1) {
              tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
            }
          }
        }
      } else if (32 * ii1 + 31 >= _PB_NI && (_PB_NK % 32) + _PB_NJ >= _PB_NK + 32) {
        for (int ii4 = 0; ii4 <= floord(_PB_NK - 1, 32); ii4 += 1) {
          for (int i1 = 32 * ii1; i1 < _PB_NI; i1 += 1) {
            for (int i2 = ((31 * _PB_NL + 31) % 32) + _PB_NL - 31; i2 < _PB_NJ; i2 += 1) {
              for (int i4 = 32 * ii4; i4 <= min(_PB_NK - 1, 32 * ii4 + 31); i4 += 1) {
                tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
              }
            }
          }
        }
      } else if (32 * ii1 + 31 >= _PB_NI && _PB_NK % 32 >= 1) {
        if (_PB_NK >= _PB_NJ) {
          for (int ii4 = 0; ii4 < _PB_NK / 32; ii4 += 1) {
            for (int i1 = 32 * ii1; i1 < _PB_NI; i1 += 1) {
              for (int i2 = ((31 * _PB_NL + 31) % 32) + _PB_NL - 31; i2 < _PB_NJ; i2 += 1) {
                for (int i4 = 32 * ii4; i4 <= 32 * ii4 + 31; i4 += 1) {
                  tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
                }
              }
            }
          }
        } else {
          for (int ii4 = 0; ii4 < _PB_NL - (31 * _PB_NL + 31) / 32; ii4 += 1) {
            for (int i1 = 32 * ii1; i1 < _PB_NI; i1 += 1) {
              for (int i2 = ((31 * _PB_NL + 31) % 32) + _PB_NL - 31; i2 < _PB_NJ; i2 += 1) {
                for (int i4 = 32 * ii4; i4 <= 32 * ii4 + 31; i4 += 1) {
                  tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
                }
              }
            }
          }
        }
        for (int i1 = 32 * ii1; i1 < _PB_NI; i1 += 1) {
          for (int i2 = ((31 * _PB_NL + 31) % 32) + _PB_NL - 31; i2 < _PB_NJ; i2 += 1) {
            for (int i4 = -((_PB_NK + 31) % 32) + _PB_NK - 1; i4 < _PB_NK; i4 += 1) {
              tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
            }
          }
        }
      } else if (_PB_NJ >= _PB_NK && _PB_NI >= 32 * ii1 + 32) {
        for (int ii4 = 0; ii4 < _PB_NK / 32; ii4 += 1) {
          for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
            for (int i2 = _PB_NK; i2 < _PB_NJ; i2 += 1) {
              for (int i4 = 32 * ii4; i4 <= 32 * ii4 + 31; i4 += 1) {
                tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
              }
            }
          }
        }
      } else if (_PB_NJ >= _PB_NK && 32 * ii1 + 31 >= _PB_NI) {
        for (int ii4 = 0; ii4 < _PB_NK / 32; ii4 += 1) {
          for (int i1 = 32 * ii1; i1 < _PB_NI; i1 += 1) {
            for (int i2 = _PB_NK; i2 < _PB_NJ; i2 += 1) {
              for (int i4 = 32 * ii4; i4 <= 32 * ii4 + 31; i4 += 1) {
                tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
              }
            }
          }
        }
      } else if (_PB_NI >= 32 * ii1 + 32) {
        for (int ii4 = 0; ii4 < _PB_NK / 32; ii4 += 1) {
          for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
            for (int i2 = ((31 * _PB_NL + 31) % 32) + _PB_NL - 31; i2 < _PB_NJ; i2 += 1) {
              for (int i4 = 32 * ii4; i4 <= 32 * ii4 + 31; i4 += 1) {
                tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
              }
            }
          }
        }
      } else {
        for (int ii4 = 0; ii4 < _PB_NK / 32; ii4 += 1) {
          for (int i1 = 32 * ii1; i1 < _PB_NI; i1 += 1) {
            for (int i2 = ((31 * _PB_NL + 31) % 32) + _PB_NL - 31; i2 < _PB_NJ; i2 += 1) {
              for (int i4 = 32 * ii4; i4 <= 32 * ii4 + 31; i4 += 1) {
                tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
              }
            }
          }
        }
      }
    }
  } else if (_PB_NL == 0) {
    for (int ii1 = 0; ii1 <= floord(_PB_NI - 1, 32); ii1 += 1) {
      for (int ii2 = 0; ii2 <= floord(_PB_NJ - 1, 32); ii2 += 1) {
        for (int i1 = 32 * ii1; i1 <= min(_PB_NI - 1, 32 * ii1 + 31); i1 += 1) {
          for (int i2 = 32 * ii2; i2 <= min(_PB_NJ - 1, 32 * ii2 + 31); i2 += 1) {
            tmp[i1][i2] = SCALAR_VAL(0.0);
          }
        }
        for (int ii4 = 0; ii4 <= floord(_PB_NK - 1, 32); ii4 += 1) {
          for (int i1 = 32 * ii1; i1 <= min(_PB_NI - 1, 32 * ii1 + 31); i1 += 1) {
            for (int i2 = 32 * ii2; i2 <= min(_PB_NJ - 1, 32 * ii2 + 31); i2 += 1) {
              for (int i4 = 32 * ii4; i4 <= min(_PB_NK - 1, 32 * ii4 + 31); i4 += 1) {
                tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
              }
            }
          }
        }
      }
    }
  }
  if (_PB_NI <= 31 && _PB_NJ >= _PB_NL && _PB_NJ <= 31 && _PB_NK >= 1) {
    for (int i1 = 0; i1 < _PB_NI; i1 += 1) {
      for (int i2 = 0; i2 < _PB_NL; i2 += 1) {
        D[i1][i2] *= beta;
      }
    }
    for (int i1 = 0; i1 < _PB_NI; i1 += 1) {
      for (int i2 = 0; i2 < _PB_NL; i2 += 1) {
        for (int i4 = 0; i4 < _PB_NJ; i4 += 1) {
          D[i1][i2] += (tmp[i1][i4] * C[i4][i2]);
        }
      }
    }
  } else if (_PB_NI <= 31 && _PB_NL <= 31 && _PB_NJ >= 1 && _PB_NL >= _PB_NJ + 1 && _PB_NK >= 1) {
    for (int i1 = 0; i1 < _PB_NI; i1 += 1) {
      for (int i2 = 0; i2 < _PB_NL; i2 += 1) {
        D[i1][i2] *= beta;
      }
    }
    for (int i1 = 0; i1 < _PB_NI; i1 += 1) {
      for (int i2 = 0; i2 < _PB_NL; i2 += 1) {
        for (int i4 = 0; i4 < _PB_NJ; i4 += 1) {
          D[i1][i2] += (tmp[i1][i4] * C[i4][i2]);
        }
      }
    }
  } else if (_PB_NJ >= 32 && _PB_NJ >= _PB_NL && _PB_NK >= 1 && _PB_NI % 32 >= 1) {
    for (int ii1 = 0; ii1 <= (_PB_NI - 1) / 32; ii1 += 1) {
      if (32 * ii1 + 31 >= _PB_NI) {
        for (int ii2 = 0; ii2 <= floord(_PB_NL - 1, 32); ii2 += 1) {
          for (int i1 = 32 * ii1; i1 < _PB_NI; i1 += 1) {
            for (int i2 = 32 * ii2; i2 <= min(_PB_NL - 1, 32 * ii2 + 31); i2 += 1) {
              D[i1][i2] *= beta;
            }
          }
          if (_PB_NK >= _PB_NJ) {
            for (int ii4 = 0; ii4 <= min((_PB_NJ - 1) / 32, _PB_NK / 32 - 1); ii4 += 1) {
              for (int i1 = 32 * ii1; i1 < _PB_NI; i1 += 1) {
                for (int i2 = 32 * ii2; i2 <= min(_PB_NL - 1, 32 * ii2 + 31); i2 += 1) {
                  for (int i4 = 32 * ii4; i4 <= min(_PB_NJ - 1, 32 * ii4 + 31); i4 += 1) {
                    D[i1][i2] += (tmp[i1][i4] * C[i4][i2]);
                  }
                }
              }
            }
          } else {
            for (int ii4 = 0; ii4 < _PB_NJ / 32; ii4 += 1) {
              for (int i1 = 32 * ii1; i1 < _PB_NI; i1 += 1) {
                for (int i2 = 32 * ii2; i2 <= min(_PB_NL - 1, 32 * ii2 + 31); i2 += 1) {
                  for (int i4 = 32 * ii4; i4 <= 32 * ii4 + 31; i4 += 1) {
                    D[i1][i2] += (tmp[i1][i4] * C[i4][i2]);
                  }
                }
              }
            }
          }
          if (_PB_NJ >= 32 * ii2 + 32 && _PB_NJ % 32 >= 1 && _PB_NJ + 31 >= (_PB_NJ % 32) + _PB_NK) {
            for (int i1 = 32 * ii1; i1 < _PB_NI; i1 += 1) {
              for (int i2 = 32 * ii2; i2 <= min(_PB_NL - 1, 32 * ii2 + 31); i2 += 1) {
                for (int i4 = -((_PB_NJ - 1) % 32) + _PB_NJ - 1; i4 < _PB_NJ; i4 += 1) {
                  D[i1][i2] += (tmp[i1][i4] * C[i4][i2]);
                }
              }
            }
          } else if (32 * ii2 + 31 >= _PB_NJ && 32 * ii2 + 31 >= _PB_NK) {
            for (int i1 = 32 * ii1; i1 < _PB_NI; i1 += 1) {
              for (int i2 = 32 * ii2; i2 < _PB_NL; i2 += 1) {
                for (int i4 = 32 * ii2; i4 < _PB_NJ; i4 += 1) {
                  D[i1][i2] += (tmp[i1][i4] * C[i4][i2]);
                }
              }
            }
          }
        }
      } else {
        for (int ii2 = 0; ii2 <= floord(_PB_NL - 1, 32); ii2 += 1) {
          for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
            for (int i2 = 32 * ii2; i2 <= min(_PB_NL - 1, 32 * ii2 + 31); i2 += 1) {
              D[i1][i2] *= beta;
            }
          }
          if (_PB_NK >= _PB_NJ) {
            for (int ii4 = 0; ii4 <= min((_PB_NJ - 1) / 32, _PB_NK / 32 - 1); ii4 += 1) {
              for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
                for (int i2 = 32 * ii2; i2 <= min(_PB_NL - 1, 32 * ii2 + 31); i2 += 1) {
                  for (int i4 = 32 * ii4; i4 <= min(_PB_NJ - 1, 32 * ii4 + 31); i4 += 1) {
                    D[i1][i2] += (tmp[i1][i4] * C[i4][i2]);
                  }
                }
              }
            }
          } else {
            for (int ii4 = 0; ii4 < _PB_NJ / 32; ii4 += 1) {
              for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
                for (int i2 = 32 * ii2; i2 <= min(_PB_NL - 1, 32 * ii2 + 31); i2 += 1) {
                  for (int i4 = 32 * ii4; i4 <= 32 * ii4 + 31; i4 += 1) {
                    D[i1][i2] += (tmp[i1][i4] * C[i4][i2]);
                  }
                }
              }
            }
          }
          if (_PB_NJ >= 32 * ii2 + 32 && _PB_NJ % 32 >= 1 && _PB_NJ + 31 >= (_PB_NJ % 32) + _PB_NK) {
            for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
              for (int i2 = 32 * ii2; i2 <= min(_PB_NL - 1, 32 * ii2 + 31); i2 += 1) {
                for (int i4 = -((_PB_NJ - 1) % 32) + _PB_NJ - 1; i4 < _PB_NJ; i4 += 1) {
                  D[i1][i2] += (tmp[i1][i4] * C[i4][i2]);
                }
              }
            }
          } else if (32 * ii2 + 31 >= _PB_NJ && 32 * ii2 + 31 >= _PB_NK) {
            for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
              for (int i2 = 32 * ii2; i2 < _PB_NL; i2 += 1) {
                for (int i4 = 32 * ii2; i4 < _PB_NJ; i4 += 1) {
                  D[i1][i2] += (tmp[i1][i4] * C[i4][i2]);
                }
              }
            }
          }
        }
      }
    }
  } else {
    if (_PB_NL >= _PB_NJ + 1 && _PB_NJ + 31 >= _PB_NL && _PB_NK >= 1 && _PB_NJ % 32 == 0) {
      for (int ii1 = 0; ii1 <= floord(_PB_NI - 1, 32); ii1 += 1) {
        for (int ii2 = 0; ii2 < _PB_NJ / 32; ii2 += 1) {
          for (int i1 = 32 * ii1; i1 <= min(_PB_NI - 1, 32 * ii1 + 31); i1 += 1) {
            for (int i2 = 32 * ii2; i2 <= 32 * ii2 + 31; i2 += 1) {
              tmp[i1][i2] = SCALAR_VAL(0.0);
            }
          }
          if (_PB_NI >= 32 * ii1 + 32) {
            for (int ii4 = 0; ii4 <= floord(_PB_NK - 1, 32); ii4 += 1) {
              for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
                for (int i2 = 32 * ii2; i2 <= 32 * ii2 + 31; i2 += 1) {
                  for (int i4 = 32 * ii4; i4 <= min(_PB_NK - 1, 32 * ii4 + 31); i4 += 1) {
                    tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
                  }
                }
              }
            }
          } else if (_PB_NJ >= _PB_NK + 1) {
            for (int ii4 = 0; ii4 <= floord(_PB_NK - 1, 32); ii4 += 1) {
              for (int i1 = 32 * ii1; i1 < _PB_NI; i1 += 1) {
                for (int i2 = 32 * ii2; i2 <= 32 * ii2 + 31; i2 += 1) {
                  for (int i4 = 32 * ii4; i4 <= min(_PB_NK - 1, 32 * ii4 + 31); i4 += 1) {
                    tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
                  }
                }
              }
            }
          } else if (_PB_NK % 32 >= 1) {
            for (int ii4 = 0; ii4 < _PB_NK / 32; ii4 += 1) {
              for (int i1 = 32 * ii1; i1 < _PB_NI; i1 += 1) {
                for (int i2 = 32 * ii2; i2 <= 32 * ii2 + 31; i2 += 1) {
                  for (int i4 = 32 * ii4; i4 <= 32 * ii4 + 31; i4 += 1) {
                    tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
                  }
                }
              }
            }
            for (int i1 = 32 * ii1; i1 < _PB_NI; i1 += 1) {
              for (int i2 = 32 * ii2; i2 <= 32 * ii2 + 31; i2 += 1) {
                for (int i4 = -((_PB_NK - 1) % 32) + _PB_NK - 1; i4 < _PB_NK; i4 += 1) {
                  tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
                }
              }
            }
          } else {
            for (int ii4 = 0; ii4 < _PB_NK / 32; ii4 += 1) {
              for (int i1 = 32 * ii1; i1 < _PB_NI; i1 += 1) {
                for (int i2 = 32 * ii2; i2 <= 32 * ii2 + 31; i2 += 1) {
                  for (int i4 = 32 * ii4; i4 <= 32 * ii4 + 31; i4 += 1) {
                    tmp[i1][i2] += ((alpha * A[i1][i4]) * B[i4][i2]);
                  }
                }
              }
            }
          }
        }
      }
    }
    if (_PB_NL >= 32 && _PB_NJ >= 1 && _PB_NL >= _PB_NJ + 1 && _PB_NK >= 1 && _PB_NI % 32 >= 1) {
      if (_PB_NK >= _PB_NJ && _PB_NK <= 31) {
        for (int ii1 = 0; ii1 < _PB_NI / 32; ii1 += 1) {
          for (int ii2 = 0; ii2 <= (_PB_NL - 1) / 32; ii2 += 1) {
            for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
              for (int i2 = 32 * ii2; i2 <= min(_PB_NL - 1, 32 * ii2 + 31); i2 += 1) {
                D[i1][i2] *= beta;
              }
            }
            for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
              for (int i2 = 32 * ii2; i2 <= min(_PB_NL - 1, 32 * ii2 + 31); i2 += 1) {
                for (int i4 = 0; i4 < _PB_NJ; i4 += 1) {
                  D[i1][i2] += (tmp[i1][i4] * C[i4][i2]);
                }
              }
            }
          }
        }
        for (int ii2 = 0; ii2 <= (_PB_NL - 1) / 32; ii2 += 1) {
          for (int i1 = -((_PB_NI + 31) % 32) + _PB_NI - 1; i1 < _PB_NI; i1 += 1) {
            for (int i2 = 32 * ii2; i2 <= min(_PB_NL - 1, 32 * ii2 + 31); i2 += 1) {
              D[i1][i2] *= beta;
            }
          }
          for (int i1 = -((_PB_NI + 31) % 32) + _PB_NI - 1; i1 < _PB_NI; i1 += 1) {
            for (int i2 = 32 * ii2; i2 <= min(_PB_NL - 1, 32 * ii2 + 31); i2 += 1) {
              for (int i4 = 0; i4 < _PB_NJ; i4 += 1) {
                D[i1][i2] += (tmp[i1][i4] * C[i4][i2]);
              }
            }
          }
        }
      } else if (_PB_NJ >= _PB_NK + 1) {
        for (int ii1 = 0; ii1 < _PB_NI / 32; ii1 += 1) {
          for (int ii2 = 0; ii2 <= (_PB_NL - 1) / 32; ii2 += 1) {
            for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
              for (int i2 = 32 * ii2; i2 <= min(_PB_NL - 1, 32 * ii2 + 31); i2 += 1) {
                D[i1][i2] *= beta;
              }
            }
            if (_PB_NJ >= 32) {
              for (int ii4 = 0; ii4 <= (_PB_NJ - 1) / 32; ii4 += 1) {
                for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
                  for (int i2 = 32 * ii2; i2 <= min(_PB_NL - 1, 32 * ii2 + 31); i2 += 1) {
                    for (int i4 = 32 * ii4; i4 <= min(_PB_NJ - 1, 32 * ii4 + 31); i4 += 1) {
                      D[i1][i2] += (tmp[i1][i4] * C[i4][i2]);
                    }
                  }
                }
              }
            } else {
              for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
                for (int i2 = 32 * ii2; i2 <= min(_PB_NL - 1, 32 * ii2 + 31); i2 += 1) {
                  for (int i4 = 0; i4 < _PB_NJ; i4 += 1) {
                    D[i1][i2] += (tmp[i1][i4] * C[i4][i2]);
                  }
                }
              }
            }
          }
        }
        for (int ii2 = 0; ii2 <= (_PB_NL - 1) / 32; ii2 += 1) {
          for (int i1 = -((_PB_NI + 31) % 32) + _PB_NI - 1; i1 < _PB_NI; i1 += 1) {
            for (int i2 = 32 * ii2; i2 <= min(_PB_NL - 1, 32 * ii2 + 31); i2 += 1) {
              D[i1][i2] *= beta;
            }
          }
          if (_PB_NJ >= 32) {
            for (int ii4 = 0; ii4 <= (_PB_NJ - 1) / 32; ii4 += 1) {
              for (int i1 = -((_PB_NI + 31) % 32) + _PB_NI - 1; i1 < _PB_NI; i1 += 1) {
                for (int i2 = 32 * ii2; i2 <= min(_PB_NL - 1, 32 * ii2 + 31); i2 += 1) {
                  for (int i4 = 32 * ii4; i4 <= min(_PB_NJ - 1, 32 * ii4 + 31); i4 += 1) {
                    D[i1][i2] += (tmp[i1][i4] * C[i4][i2]);
                  }
                }
              }
            }
          } else {
            for (int i1 = -((_PB_NI + 31) % 32) + _PB_NI - 1; i1 < _PB_NI; i1 += 1) {
              for (int i2 = 32 * ii2; i2 <= min(_PB_NL - 1, 32 * ii2 + 31); i2 += 1) {
                for (int i4 = 0; i4 < _PB_NJ; i4 += 1) {
                  D[i1][i2] += (tmp[i1][i4] * C[i4][i2]);
                }
              }
            }
          }
        }
      } else {
        for (int ii1 = 0; ii1 <= (_PB_NI - 1) / 32; ii1 += 1) {
          if (32 * ii1 + 31 >= _PB_NI) {
            for (int ii2 = 0; ii2 <= (_PB_NL - 1) / 32; ii2 += 1) {
              for (int i1 = 32 * ii1; i1 < _PB_NI; i1 += 1) {
                for (int i2 = 32 * ii2; i2 <= min(_PB_NL - 1, 32 * ii2 + 31); i2 += 1) {
                  D[i1][i2] *= beta;
                }
              }
              for (int ii4 = 0; ii4 <= floord(_PB_NJ - 1, 32); ii4 += 1) {
                for (int i1 = 32 * ii1; i1 < _PB_NI; i1 += 1) {
                  for (int i2 = 32 * ii2; i2 <= min(_PB_NL - 1, 32 * ii2 + 31); i2 += 1) {
                    for (int i4 = 32 * ii4; i4 <= min(_PB_NJ - 1, 32 * ii4 + 31); i4 += 1) {
                      D[i1][i2] += (tmp[i1][i4] * C[i4][i2]);
                    }
                  }
                }
              }
            }
          } else {
            for (int ii2 = 0; ii2 <= (_PB_NL - 1) / 32; ii2 += 1) {
              for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
                for (int i2 = 32 * ii2; i2 <= min(_PB_NL - 1, 32 * ii2 + 31); i2 += 1) {
                  D[i1][i2] *= beta;
                }
              }
              for (int ii4 = 0; ii4 <= floord(_PB_NJ - 1, 32); ii4 += 1) {
                for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
                  for (int i2 = 32 * ii2; i2 <= min(_PB_NL - 1, 32 * ii2 + 31); i2 += 1) {
                    for (int i4 = 32 * ii4; i4 <= min(_PB_NJ - 1, 32 * ii4 + 31); i4 += 1) {
                      D[i1][i2] += (tmp[i1][i4] * C[i4][i2]);
                    }
                  }
                }
              }
            }
          }
        }
      }
    } else if (_PB_NI >= 32 && _PB_NJ >= _PB_NL && _PB_NJ <= 31 && _PB_NK >= 1 && _PB_NI % 32 >= 1) {
      for (int ii1 = 0; ii1 < _PB_NI / 32; ii1 += 1) {
        for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
          for (int i2 = 0; i2 < _PB_NL; i2 += 1) {
            D[i1][i2] *= beta;
          }
        }
        for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
          for (int i2 = 0; i2 < _PB_NL; i2 += 1) {
            for (int i4 = 0; i4 < _PB_NJ; i4 += 1) {
              D[i1][i2] += (tmp[i1][i4] * C[i4][i2]);
            }
          }
        }
      }
      for (int i1 = -((_PB_NI - 1) % 32) + _PB_NI - 1; i1 < _PB_NI; i1 += 1) {
        for (int i2 = 0; i2 < _PB_NL; i2 += 1) {
          D[i1][i2] *= beta;
        }
      }
      for (int i1 = -((_PB_NI - 1) % 32) + _PB_NI - 1; i1 < _PB_NI; i1 += 1) {
        for (int i2 = 0; i2 < _PB_NL; i2 += 1) {
          for (int i4 = 0; i4 < _PB_NJ; i4 += 1) {
            D[i1][i2] += (tmp[i1][i4] * C[i4][i2]);
          }
        }
      }
    } else if (_PB_NI >= 32 && _PB_NL <= 31 && _PB_NJ >= 1 && _PB_NL >= _PB_NJ + 1 && _PB_NK >= 1 && _PB_NI % 32 >= 1) {
      for (int ii1 = 0; ii1 < _PB_NI / 32; ii1 += 1) {
        for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
          for (int i2 = 0; i2 < _PB_NL; i2 += 1) {
            D[i1][i2] *= beta;
          }
        }
        for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
          for (int i2 = 0; i2 < _PB_NL; i2 += 1) {
            for (int i4 = 0; i4 < _PB_NJ; i4 += 1) {
              D[i1][i2] += (tmp[i1][i4] * C[i4][i2]);
            }
          }
        }
      }
      for (int i1 = -((_PB_NI - 1) % 32) + _PB_NI - 1; i1 < _PB_NI; i1 += 1) {
        for (int i2 = 0; i2 < _PB_NL; i2 += 1) {
          D[i1][i2] *= beta;
        }
      }
      for (int i1 = -((_PB_NI - 1) % 32) + _PB_NI - 1; i1 < _PB_NI; i1 += 1) {
        for (int i2 = 0; i2 < _PB_NL; i2 += 1) {
          for (int i4 = 0; i4 < _PB_NJ; i4 += 1) {
            D[i1][i2] += (tmp[i1][i4] * C[i4][i2]);
          }
        }
      }
    } else if (_PB_NK == 0 && (_PB_NJ % 32) + _PB_NL >= _PB_NJ + 32) {
      for (int ii1 = 0; ii1 <= floord(_PB_NI - 1, 32); ii1 += 1) {
        for (int ii2 = 0; ii2 <= floord(_PB_NJ - 1, 32); ii2 += 1) {
          for (int i1 = 32 * ii1; i1 <= min(_PB_NI - 1, 32 * ii1 + 31); i1 += 1) {
            for (int i2 = 32 * ii2; i2 <= min(_PB_NJ - 1, 32 * ii2 + 31); i2 += 1) {
              tmp[i1][i2] = SCALAR_VAL(0.0);
            }
          }
        }
      }
    } else if (_PB_NL >= _PB_NJ && _PB_NJ + 31 >= _PB_NL && _PB_NK == 0 && _PB_NJ % 32 == 0) {
      for (int ii1 = 0; ii1 <= floord(_PB_NI - 1, 32); ii1 += 1) {
        for (int ii2 = 0; ii2 < _PB_NJ / 32; ii2 += 1) {
          for (int i1 = 32 * ii1; i1 <= min(_PB_NI - 1, 32 * ii1 + 31); i1 += 1) {
            for (int i2 = 32 * ii2; i2 <= 32 * ii2 + 31; i2 += 1) {
              tmp[i1][i2] = SCALAR_VAL(0.0);
            }
          }
        }
      }
    }
  }
  if (_PB_NJ == 0) {
    for (int ii1 = 0; ii1 <= floord(_PB_NI - 1, 32); ii1 += 1) {
      for (int ii2 = 0; ii2 <= floord(_PB_NL - 1, 32); ii2 += 1) {
        for (int i1 = 32 * ii1; i1 <= min(_PB_NI - 1, 32 * ii1 + 31); i1 += 1) {
          for (int i2 = 32 * ii2; i2 <= min(_PB_NL - 1, 32 * ii2 + 31); i2 += 1) {
            D[i1][i2] *= beta;
          }
        }
      }
    }
  } else if (_PB_NJ >= 32 && _PB_NJ >= _PB_NL + 1 && _PB_NK == 0) {
    for (int ii1 = 0; ii1 <= floord(_PB_NI - 1, 32); ii1 += 1) {
      for (int ii2 = 0; ii2 <= floord(_PB_NL - 1, 32); ii2 += 1) {
        for (int i1 = 32 * ii1; i1 <= min(_PB_NI - 1, 32 * ii1 + 31); i1 += 1) {
          for (int i2 = 32 * ii2; i2 <= min(_PB_NL - 1, 32 * ii2 + 31); i2 += 1) {
            D[i1][i2] *= beta;
          }
        }
        for (int ii4 = 0; ii4 <= (_PB_NJ - 1) / 32; ii4 += 1) {
          for (int i1 = 32 * ii1; i1 <= min(_PB_NI - 1, 32 * ii1 + 31); i1 += 1) {
            for (int i2 = 32 * ii2; i2 <= min(_PB_NL - 1, 32 * ii2 + 31); i2 += 1) {
              for (int i4 = 32 * ii4; i4 <= min(_PB_NJ - 1, 32 * ii4 + 31); i4 += 1) {
                D[i1][i2] += (tmp[i1][i4] * C[i4][i2]);
              }
            }
          }
        }
      }
    }
  } else if (_PB_NJ >= _PB_NL + 1 && _PB_NJ <= 31 && _PB_NK == 0) {
    for (int ii1 = 0; ii1 <= floord(_PB_NI - 1, 32); ii1 += 1) {
      for (int i1 = 32 * ii1; i1 <= min(_PB_NI - 1, 32 * ii1 + 31); i1 += 1) {
        for (int i2 = 0; i2 < _PB_NL; i2 += 1) {
          D[i1][i2] *= beta;
        }
      }
      for (int i1 = 32 * ii1; i1 <= min(_PB_NI - 1, 32 * ii1 + 31); i1 += 1) {
        for (int i2 = 0; i2 < _PB_NL; i2 += 1) {
          for (int i4 = 0; i4 < _PB_NJ; i4 += 1) {
            D[i1][i2] += (tmp[i1][i4] * C[i4][i2]);
          }
        }
      }
    }
  } else {
    if (_PB_NL >= 32 && _PB_NJ == _PB_NL && _PB_NK == 0 && _PB_NL % 32 >= 1) {
      for (int ii1 = 0; ii1 <= floord(_PB_NI - 1, 32); ii1 += 1) {
        for (int ii2 = 0; ii2 < _PB_NL / 32; ii2 += 1) {
          for (int i1 = 32 * ii1; i1 <= min(_PB_NI - 1, 32 * ii1 + 31); i1 += 1) {
            for (int i2 = 32 * ii2; i2 <= 32 * ii2 + 31; i2 += 1) {
              tmp[i1][i2] = SCALAR_VAL(0.0);
            }
          }
        }
        for (int i1 = 32 * ii1; i1 <= min(_PB_NI - 1, 32 * ii1 + 31); i1 += 1) {
          for (int i2 = -((_PB_NL - 1) % 32) + _PB_NL - 1; i2 < _PB_NL; i2 += 1) {
            tmp[i1][i2] = SCALAR_VAL(0.0);
          }
        }
      }
    }
    if (_PB_NL >= _PB_NJ && _PB_NK == 0) {
      for (int ii1 = 0; ii1 <= floord(_PB_NI - 1, 32); ii1 += 1) {
        for (int ii2 = 0; ii2 <= (_PB_NL - 1) / 32; ii2 += 1) {
          for (int i1 = 32 * ii1; i1 <= min(_PB_NI - 1, 32 * ii1 + 31); i1 += 1) {
            for (int i2 = 32 * ii2; i2 <= min(_PB_NL - 1, 32 * ii2 + 31); i2 += 1) {
              D[i1][i2] *= beta;
            }
          }
          for (int ii4 = 0; ii4 <= floord(_PB_NJ - 1, 32); ii4 += 1) {
            for (int i1 = 32 * ii1; i1 <= min(_PB_NI - 1, 32 * ii1 + 31); i1 += 1) {
              for (int i2 = 32 * ii2; i2 <= min(_PB_NL - 1, 32 * ii2 + 31); i2 += 1) {
                for (int i4 = 32 * ii4; i4 <= min(_PB_NJ - 1, 32 * ii4 + 31); i4 += 1) {
                  D[i1][i2] += (tmp[i1][i4] * C[i4][i2]);
                }
              }
            }
          }
        }
      }
    }
  }
  if (_PB_NJ >= _PB_NL && _PB_NK >= 1 && _PB_NI % 32 == 0) {
    for (int ii1 = 0; ii1 < _PB_NI / 32; ii1 += 1) {
      for (int ii2 = 0; ii2 <= floord(_PB_NL - 1, 32); ii2 += 1) {
        for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
          for (int i2 = 32 * ii2; i2 <= min(_PB_NL - 1, 32 * ii2 + 31); i2 += 1) {
            D[i1][i2] *= beta;
          }
        }
        if (_PB_NJ >= 32 * ii2 + 32) {
          if (_PB_NK >= _PB_NJ) {
            for (int ii4 = 0; ii4 <= min((_PB_NJ - 1) / 32, _PB_NK / 32 - 1); ii4 += 1) {
              for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
                for (int i2 = 32 * ii2; i2 <= min(_PB_NL - 1, 32 * ii2 + 31); i2 += 1) {
                  for (int i4 = 32 * ii4; i4 <= min(_PB_NJ - 1, 32 * ii4 + 31); i4 += 1) {
                    D[i1][i2] += (tmp[i1][i4] * C[i4][i2]);
                  }
                }
              }
            }
          } else {
            for (int ii4 = 0; ii4 < _PB_NJ / 32; ii4 += 1) {
              for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
                for (int i2 = 32 * ii2; i2 <= min(_PB_NL - 1, 32 * ii2 + 31); i2 += 1) {
                  for (int i4 = 32 * ii4; i4 <= 32 * ii4 + 31; i4 += 1) {
                    D[i1][i2] += (tmp[i1][i4] * C[i4][i2]);
                  }
                }
              }
            }
          }
          if (_PB_NJ % 32 >= 1 && _PB_NJ + 31 >= (_PB_NJ % 32) + _PB_NK) {
            for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
              for (int i2 = 32 * ii2; i2 <= min(_PB_NL - 1, 32 * ii2 + 31); i2 += 1) {
                for (int i4 = -((_PB_NJ - 1) % 32) + _PB_NJ - 1; i4 < _PB_NJ; i4 += 1) {
                  D[i1][i2] += (tmp[i1][i4] * C[i4][i2]);
                }
              }
            }
          }
        } else if (_PB_NK >= 32 * ii2 + 32) {
          for (int ii4 = 0; ii4 <= ii2; ii4 += 1) {
            for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
              for (int i2 = 32 * ii2; i2 < _PB_NL; i2 += 1) {
                for (int i4 = 32 * ii4; i4 <= min(_PB_NJ - 1, 32 * ii4 + 31); i4 += 1) {
                  D[i1][i2] += (tmp[i1][i4] * C[i4][i2]);
                }
              }
            }
          }
        } else {
          if (_PB_NK >= _PB_NJ) {
            for (int ii4 = 0; ii4 < ii2; ii4 += 1) {
              for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
                for (int i2 = 32 * ii2; i2 < _PB_NL; i2 += 1) {
                  for (int i4 = 32 * ii4; i4 <= 32 * ii4 + 31; i4 += 1) {
                    D[i1][i2] += (tmp[i1][i4] * C[i4][i2]);
                  }
                }
              }
            }
          } else {
            for (int ii4 = 0; ii4 < ii2; ii4 += 1) {
              for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
                for (int i2 = 32 * ii2; i2 < _PB_NL; i2 += 1) {
                  for (int i4 = 32 * ii4; i4 <= 32 * ii4 + 31; i4 += 1) {
                    D[i1][i2] += (tmp[i1][i4] * C[i4][i2]);
                  }
                }
              }
            }
          }
          for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
            for (int i2 = 32 * ii2; i2 < _PB_NL; i2 += 1) {
              if (_PB_NK >= _PB_NJ) {
                for (int i4 = 32 * ii2; i4 < _PB_NJ; i4 += 1) {
                  D[i1][i2] += (tmp[i1][i4] * C[i4][i2]);
                }
              } else {
                for (int i4 = 32 * ii2; i4 < _PB_NJ; i4 += 1) {
                  D[i1][i2] += (tmp[i1][i4] * C[i4][i2]);
                }
              }
            }
          }
        }
      }
    }
  } else if (_PB_NJ >= 1 && _PB_NL >= _PB_NJ + 1 && _PB_NK >= 1 && _PB_NI % 32 == 0) {
    for (int ii1 = 0; ii1 < _PB_NI / 32; ii1 += 1) {
      for (int ii2 = 0; ii2 <= (_PB_NL - 1) / 32; ii2 += 1) {
        for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
          for (int i2 = 32 * ii2; i2 <= min(_PB_NL - 1, 32 * ii2 + 31); i2 += 1) {
            D[i1][i2] *= beta;
          }
        }
        if (_PB_NK >= 32 && _PB_NK >= _PB_NJ) {
          for (int ii4 = 0; ii4 <= floord(_PB_NJ - 1, 32); ii4 += 1) {
            for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
              for (int i2 = 32 * ii2; i2 <= min(_PB_NL - 1, 32 * ii2 + 31); i2 += 1) {
                for (int i4 = 32 * ii4; i4 <= min(_PB_NJ - 1, 32 * ii4 + 31); i4 += 1) {
                  D[i1][i2] += (tmp[i1][i4] * C[i4][i2]);
                }
              }
            }
          }
        } else if (_PB_NK >= _PB_NJ && _PB_NK <= 31) {
          for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
            for (int i2 = 32 * ii2; i2 <= min(_PB_NL - 1, 32 * ii2 + 31); i2 += 1) {
              for (int i4 = 0; i4 < _PB_NJ; i4 += 1) {
                D[i1][i2] += (tmp[i1][i4] * C[i4][i2]);
              }
            }
          }
        } else if (_PB_NJ >= 32) {
          for (int ii4 = 0; ii4 <= (_PB_NJ - 1) / 32; ii4 += 1) {
            for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
              for (int i2 = 32 * ii2; i2 <= min(_PB_NL - 1, 32 * ii2 + 31); i2 += 1) {
                for (int i4 = 32 * ii4; i4 <= min(_PB_NJ - 1, 32 * ii4 + 31); i4 += 1) {
                  D[i1][i2] += (tmp[i1][i4] * C[i4][i2]);
                }
              }
            }
          }
        } else {
          for (int i1 = 32 * ii1; i1 <= 32 * ii1 + 31; i1 += 1) {
            for (int i2 = 32 * ii2; i2 <= min(_PB_NL - 1, 32 * ii2 + 31); i2 += 1) {
              for (int i4 = 0; i4 < _PB_NJ; i4 += 1) {
                D[i1][i2] += (tmp[i1][i4] * C[i4][i2]);
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
  int ni = NI;
  int nj = NJ;
  int nk = NK;
  int nl = NL;

  /* Variable declaration/allocation. */
  DATA_TYPE alpha;
  DATA_TYPE beta;
  POLYBENCH_2D_ARRAY_DECL(tmp,DATA_TYPE,NI,NJ,ni,nj);
  POLYBENCH_2D_ARRAY_DECL(A,DATA_TYPE,NI,NK,ni,nk);
  POLYBENCH_2D_ARRAY_DECL(B,DATA_TYPE,NK,NJ,nk,nj);
  POLYBENCH_2D_ARRAY_DECL(C,DATA_TYPE,NJ,NL,nj,nl);
  POLYBENCH_2D_ARRAY_DECL(D,DATA_TYPE,NI,NL,ni,nl);

  /* Initialize array(s). */
  init_array (ni, nj, nk, nl, &alpha, &beta,
	      POLYBENCH_ARRAY(A),
	      POLYBENCH_ARRAY(B),
	      POLYBENCH_ARRAY(C),
	      POLYBENCH_ARRAY(D));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_2mm (ni, nj, nk, nl,
	      alpha, beta,
	      POLYBENCH_ARRAY(tmp),
	      POLYBENCH_ARRAY(A),
	      POLYBENCH_ARRAY(B),
	      POLYBENCH_ARRAY(C),
	      POLYBENCH_ARRAY(D));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(ni, nl,  POLYBENCH_ARRAY(D)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(tmp);
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(B);
  POLYBENCH_FREE_ARRAY(C);
  POLYBENCH_FREE_ARRAY(D);

  return 0;
}
