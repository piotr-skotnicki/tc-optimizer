/**
 * This version is stamped on Apr. 14, 2015
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* ludcmp.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "ludcmp.h"


/* Array initialization. */
static
void init_array (int n,
		 DATA_TYPE POLYBENCH_2D(A,N,N,n,n),
		 DATA_TYPE POLYBENCH_1D(b,N,n),
		 DATA_TYPE POLYBENCH_1D(x,N,n),
		 DATA_TYPE POLYBENCH_1D(y,N,n))
{
  int i, j;
  DATA_TYPE fn = (DATA_TYPE)n;

  for (i = 0; i < n; i++)
    {
      x[i] = 0;
      y[i] = 0;
      b[i] = (i+1)/fn/2.0 + 4;
    }

  for (i = 0; i < n; i++)
    {
      for (j = 0; j <= i; j++)
	A[i][j] = (DATA_TYPE)(-j % n) / n + 1;
      for (j = i+1; j < n; j++) {
	A[i][j] = 0;
      }
      A[i][i] = 1;
    }

  /* Make the matrix positive semi-definite. */
  /* not necessary for LU, but using same code as cholesky */
  int r,s,t;
  POLYBENCH_2D_ARRAY_DECL(B, DATA_TYPE, N, N, n, n);
  for (r = 0; r < n; ++r)
    for (s = 0; s < n; ++s)
      (POLYBENCH_ARRAY(B))[r][s] = 0;
  for (t = 0; t < n; ++t)
    for (r = 0; r < n; ++r)
      for (s = 0; s < n; ++s)
	(POLYBENCH_ARRAY(B))[r][s] += A[r][t] * A[s][t];
    for (r = 0; r < n; ++r)
      for (s = 0; s < n; ++s)
	A[r][s] = (POLYBENCH_ARRAY(B))[r][s];
  POLYBENCH_FREE_ARRAY(B);

}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_1D(x,N,n))

{
  int i;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("x");
  for (i = 0; i < n; i++) {
    if (i % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
    fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, x[i]);
  }
  POLYBENCH_DUMP_END("x");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_ludcmp(int n,
		   DATA_TYPE POLYBENCH_2D(A,N,N,n,n),
		   DATA_TYPE POLYBENCH_1D(b,N,n),
		   DATA_TYPE POLYBENCH_1D(x,N,n),
		   DATA_TYPE POLYBENCH_1D(y,N,n))
{
  int i, j, k;

  DATA_TYPE w;

/* TC Optimizing Compiler 0.2.26 */
/* ./tc ../examples/polybench/ludcmp.scop.c --correction-tiling --lex-scheduling --serial-codegen -b 64 --debug */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#pragma scop
{
  if (_PB_N >= 1) {
    for (int ii1 = (_PB_N - 1) / 64; ii1 <= (_PB_N - 1) / 32; ii1 += 1) {
      for (int ii3 = 0; ii3 <= ii1 - (_PB_N + 63) / 64; ii3 += 1) {
        if (_PB_N >= 32 * ii1 + 33 && ii3 >= 1) {
          w = A[-_PB_N + 64 * ii1 + 1][64 * ii3];
        } else if (32 * ii1 + 32 >= _PB_N && 64 * ii1 >= _PB_N + 1) {
          w = A[-_PB_N + 64 * ii1 + 1][64 * ii3];
        } else if (_PB_N >= 32 * ii1 + 33 && 64 * ii1 >= _PB_N + 1) {
          w = A[-_PB_N + 64 * ii1 + 1][0];
        } else {
          w = A[1][0];
        }
        for (int ii4 = max(1, -ii3 + 2); ii4 <= 2; ii4 += 1) {
          if (ii4 == 2) {
            if (_PB_N >= 32 * ii1 + 33) {
              for (int i3 = max(64, 64 * ii3); i3 <= min(-_PB_N + 64 * ii1, 64 * ii3 + 63); i3 += 1) {
                if (i3 >= 64 * ii3 + 1) {
                  w = A[-_PB_N + 64 * ii1 + 1][i3];
                  for (int i5 = 0; i5 < i3; i5 += 1) {
                    w -= (A[-_PB_N + 64 * ii1 + 1][i5] * A[i5][i3]);
                  }
                }
                A[-_PB_N + 64 * ii1 + 1][i3] = (w / A[i3][i3]);
              }
              if (64 * ii1 >= _PB_N + 1 && ii3 == 0) {
                for (int i3 = 0; i3 <= min(63, -_PB_N + 64 * ii1); i3 += 1) {
                  if (i3 >= 1) {
                    w = A[-_PB_N + 64 * ii1 + 1][i3];
                  }
                  for (int i5 = 0; i5 < i3; i5 += 1) {
                    w -= (A[-_PB_N + 64 * ii1 + 1][i5] * A[i5][i3]);
                  }
                  A[-_PB_N + 64 * ii1 + 1][i3] = (w / A[i3][i3]);
                }
              }
            } else if (64 * ii1 >= _PB_N + 1) {
              for (int i3 = 64 * ii3; i3 <= min(-_PB_N + 64 * ii1, 64 * ii3 + 63); i3 += 1) {
                if (i3 >= 64 * ii3 + 1) {
                  w = A[-_PB_N + 64 * ii1 + 1][i3];
                  for (int i5 = 0; i5 < i3; i5 += 1) {
                    w -= (A[-_PB_N + 64 * ii1 + 1][i5] * A[i5][i3]);
                  }
                }
                A[-_PB_N + 64 * ii1 + 1][i3] = (w / A[i3][i3]);
              }
            }
            if (64 * ii1 == _PB_N && ii3 == 0) {
              A[1][0] = (w / A[0][0]);
            }
          } else {
            for (int ii5 = 0; ii5 < ii3; ii5 += 1) {
              for (int i5 = 64 * ii5; i5 <= 64 * ii5 + 63; i5 += 1) {
                w -= (A[-_PB_N + 64 * ii1 + 1][i5] * A[i5][64 * ii3]);
              }
            }
          }
        }
      }
      if (_PB_N >= 32 * ii1 + 33 && 64 * ii1 >= _PB_N + 1) {
        for (int ii3 = ii1 - (_PB_N + 62) / 64; ii3 <= (_PB_N - 1) / 64; ii3 += 1) {
          for (int ii4 = 0; ii4 <= 2; ii4 += 1) {
            if (ii4 == 2) {
              if (_PB_N >= 64 * ii3 + 64) {
                for (int i3 = max(-_PB_N + 64 * ii1 + 1, 64 * ii3); i3 <= 64 * ii3 + 63; i3 += 1) {
                  if (_PB_N + i3 >= 64 * ii1 + 2 && i3 >= 64 * ii3 + 1) {
                    w = A[-_PB_N + 64 * ii1 + 1][i3];
                    for (int i5 = 0; i5 <= -_PB_N + 64 * ii1; i5 += 1) {
                      w -= (A[-_PB_N + 64 * ii1 + 1][i5] * A[i5][i3]);
                    }
                  }
                  A[-_PB_N + 64 * ii1 + 1][i3] = w;
                }
                if (64 * ii3 + 64 == _PB_N) {
                  for (int i1 = -_PB_N + 64 * ii1 + 2; i1 <= -_PB_N + 64 * ii1 + 64; i1 += 1) {
                    for (int i3 = 0; i3 < i1; i3 += 1) {
                      w = A[i1][i3];
                      for (int i5 = 0; i5 < i3; i5 += 1) {
                        w -= (A[i1][i5] * A[i5][i3]);
                      }
                      A[i1][i3] = (w / A[i3][i3]);
                    }
                    if (_PB_N >= i1 + 64) {
                      for (int i3 = i1; i3 < _PB_N - 64; i3 += 1) {
                        w = A[i1][i3];
                        for (int i5 = 0; i5 < i1; i5 += 1) {
                          w -= (A[i1][i5] * A[i5][i3]);
                        }
                        A[i1][i3] = w;
                      }
                      for (int i3 = _PB_N - 64; i3 < _PB_N; i3 += 1) {
                        w = A[i1][i3];
                        if (32 * ii1 + 64 == _PB_N && i1 + 64 == _PB_N && i3 + 64 == _PB_N) {
                          w -= (A[_PB_N - 64][0] * A[0][_PB_N - 64]);
                        }
                        for (int i5 = max(0, i1 - i3 + 1); i5 < i1; i5 += 1) {
                          w -= (A[i1][i5] * A[i5][i3]);
                        }
                        A[i1][i3] = w;
                      }
                    }
                  }
                }
              } else if (2 * ii3 >= ii1 + 1) {
                for (int i1 = -_PB_N + 64 * ii1 + 1; i1 <= -_PB_N + 64 * ii1 + 64; i1 += 1) {
                  if (_PB_N + i1 >= 64 * ii1 + 2) {
                    for (int i3 = 0; i3 < i1; i3 += 1) {
                      w = A[i1][i3];
                      for (int i5 = 0; i5 < i3; i5 += 1) {
                        w -= (A[i1][i5] * A[i5][i3]);
                      }
                      A[i1][i3] = (w / A[i3][i3]);
                    }
                  }
                  if (_PB_N + i1 >= 64 * ii1 + 2) {
                    for (int i3 = i1; i3 < 64 * ii3; i3 += 1) {
                      w = A[i1][i3];
                      if (i3 == i1) {
                        w -= (A[i1][0] * A[0][i1]);
                      }
                      for (int i5 = max(0, i1 - i3 + 1); i5 < i1; i5 += 1) {
                        w -= (A[i1][i5] * A[i5][i3]);
                      }
                      A[i1][i3] = w;
                    }
                  }
                  for (int i3 = 64 * ii3; i3 < _PB_N; i3 += 1) {
                    if (_PB_N + i1 + i3 >= 64 * ii1 + 64 * ii3 + 2) {
                      w = A[i1][i3];
                      for (int i5 = 0; i5 < i1; i5 += 1) {
                        w -= (A[i1][i5] * A[i5][i3]);
                      }
                    }
                    A[i1][i3] = w;
                  }
                }
              } else {
                for (int i1 = -_PB_N + 64 * ii1 + 1; i1 <= 32 * ii1; i1 += 1) {
                  if (_PB_N + i1 >= 64 * ii1 + 2) {
                    for (int i3 = 0; i3 < i1; i3 += 1) {
                      w = A[i1][i3];
                      for (int i5 = 0; i5 < i3; i5 += 1) {
                        w -= (A[i1][i5] * A[i5][i3]);
                      }
                      A[i1][i3] = (w / A[i3][i3]);
                    }
                  }
                  if (_PB_N + i1 >= 64 * ii1 + 2) {
                    for (int i3 = i1; i3 < 32 * ii1; i3 += 1) {
                      w = A[i1][i3];
                      if (i3 == i1) {
                        w -= (A[i1][0] * A[0][i1]);
                      }
                      for (int i5 = max(0, i1 - i3 + 1); i5 < i1; i5 += 1) {
                        w -= (A[i1][i5] * A[i5][i3]);
                      }
                      A[i1][i3] = w;
                    }
                  }
                  for (int i3 = 32 * ii1; i3 < _PB_N; i3 += 1) {
                    if (_PB_N + i1 + i3 >= 96 * ii1 + 2) {
                      w = A[i1][i3];
                      for (int i5 = 0; i5 < i1; i5 += 1) {
                        w -= (A[i1][i5] * A[i5][i3]);
                      }
                    }
                    A[i1][i3] = w;
                  }
                }
                for (int i1 = 32 * ii1 + 1; i1 <= -_PB_N + 64 * ii1 + 64; i1 += 1) {
                  for (int i3 = 0; i3 < i1; i3 += 1) {
                    w = A[i1][i3];
                    if (i3 >= 32 * ii1) {
                      for (int i5 = 0; i5 < i3; i5 += 1) {
                        w -= (A[i1][i5] * A[i5][i3]);
                      }
                    } else {
                      for (int i5 = 0; i5 < i3; i5 += 1) {
                        w -= (A[i1][i5] * A[i5][i3]);
                      }
                    }
                    A[i1][i3] = (w / A[i3][i3]);
                  }
                  for (int i3 = i1; i3 < _PB_N; i3 += 1) {
                    w = A[i1][i3];
                    for (int i5 = 0; i5 < i1; i5 += 1) {
                      w -= (A[i1][i5] * A[i5][i3]);
                    }
                    A[i1][i3] = w;
                  }
                }
              }
            } else if (ii4 == 0) {
              if (_PB_N + 64 * ii3 >= 64 * ii1 + 2) {
                w = A[-_PB_N + 64 * ii1 + 1][64 * ii3];
              } else {
                w = A[-_PB_N + 64 * ii1 + 1][-_PB_N + 64 * ii1 + 1];
              }
            } else if (_PB_N >= 64 * ii3 + 64) {
              for (int ii5 = 0; ii5 <= ii1 - (_PB_N + 63) / 64; ii5 += 1) {
                if (64 * ii1 + 1 >= _PB_N + 64 * ii3) {
                  if (ii5 == 0) {
                    w -= (A[-_PB_N + 64 * ii1 + 1][0] * A[0][-_PB_N + 64 * ii1 + 1]);
                  }
                  for (int i5 = max(1, 64 * ii5); i5 <= min(-_PB_N + 64 * ii1, 64 * ii5 + 63); i5 += 1) {
                    w -= (A[-_PB_N + 64 * ii1 + 1][i5] * A[i5][-_PB_N + 64 * ii1 + 1]);
                  }
                } else {
                  for (int i5 = 64 * ii5; i5 <= min(-_PB_N + 64 * ii1, 64 * ii5 + 63); i5 += 1) {
                    w -= (A[-_PB_N + 64 * ii1 + 1][i5] * A[i5][64 * ii3]);
                  }
                }
              }
            } else {
              for (int ii5 = 0; ii5 < ii1 - ii3 - 1; ii5 += 1) {
                for (int i5 = 64 * ii5; i5 <= 64 * ii5 + 63; i5 += 1) {
                  w -= (A[-_PB_N + 64 * ii1 + 1][i5] * A[i5][64 * ii3]);
                }
              }
              for (int i5 = 64 * ii1 - 64 * ii3 - 64; i5 <= -_PB_N + 64 * ii1; i5 += 1) {
                w -= (A[-_PB_N + 64 * ii1 + 1][i5] * A[i5][64 * ii3]);
              }
            }
          }
        }
      } else if (ii1 >= 1 && 32 * ii1 + 32 >= _PB_N) {
        if ((_PB_N + 62) % 64 <= 62) {
          w = A[-_PB_N + 64 * ii1 + 1][-_PB_N + 64 * ii1 + 1];
          if (32 * ii1 + 31 >= _PB_N && (ii1 - 1) % 2 == 0) {
            for (int ii4 = 1; ii4 <= 2; ii4 += 1) {
              if (ii4 == 2) {
                for (int i1 = -_PB_N + 64 * ii1 + 1; i1 < _PB_N; i1 += 1) {
                  if (_PB_N + i1 >= 64 * ii1 + 2) {
                    for (int i3 = 0; i3 < i1; i3 += 1) {
                      w = A[i1][i3];
                      for (int i5 = 0; i5 < i3; i5 += 1) {
                        w -= (A[i1][i5] * A[i5][i3]);
                      }
                      A[i1][i3] = (w / A[i3][i3]);
                    }
                  }
                  for (int i3 = i1; i3 < _PB_N; i3 += 1) {
                    if (_PB_N + i3 >= 64 * ii1 + 2) {
                      w = A[i1][i3];
                      for (int i5 = 0; i5 < i1; i5 += 1) {
                        w -= (A[i1][i5] * A[i5][i3]);
                      }
                    }
                    A[i1][i3] = w;
                  }
                }
              } else {
                for (int ii5 = 0; ii5 <= (ii1 - 1) / 2; ii5 += 1) {
                  for (int i5 = 64 * ii5; i5 <= min(-_PB_N + 64 * ii1, 64 * ii5 + 63); i5 += 1) {
                    w -= (A[-_PB_N + 64 * ii1 + 1][i5] * A[i5][-_PB_N + 64 * ii1 + 1]);
                  }
                }
              }
            }
          } else if (2 * _PB_N >= ((_PB_N + 62) % 64) + 64 * ii1 + 2) {
            for (int ii4 = 1; ii4 <= 2; ii4 += 1) {
              if (64 * ii1 >= _PB_N + 1 && ii4 == 1) {
                for (int ii5 = 0; ii5 <= ii1 - (_PB_N + 63) / 64; ii5 += 1) {
                  if (ii5 == 0) {
                    w -= (A[-_PB_N + 64 * ii1 + 1][0] * A[0][-_PB_N + 64 * ii1 + 1]);
                  }
                  for (int i5 = max(1, 64 * ii5); i5 <= min(-_PB_N + 64 * ii1, 64 * ii5 + 63); i5 += 1) {
                    w -= (A[-_PB_N + 64 * ii1 + 1][i5] * A[i5][-_PB_N + 64 * ii1 + 1]);
                  }
                }
              } else if (ii4 == 2 && ii1 % 2 == 0) {
                for (int i3 = -_PB_N + 64 * ii1 + 1; i3 < 32 * ii1; i3 += 1) {
                  if (_PB_N + i3 >= 64 * ii1 + 2) {
                    w = A[-_PB_N + 64 * ii1 + 1][i3];
                    for (int i5 = 0; i5 <= -_PB_N + 64 * ii1; i5 += 1) {
                      w -= (A[-_PB_N + 64 * ii1 + 1][i5] * A[i5][i3]);
                    }
                  }
                  A[-_PB_N + 64 * ii1 + 1][i3] = w;
                }
              } else if (_PB_N == 64 && ii1 == 1 && ii4 == 1) {
                w -= (A[1][0] * A[0][1]);
              } else {
                for (int i1 = _PB_N - 63; i1 < _PB_N; i1 += 1) {
                  if (i1 + 62 >= _PB_N) {
                    for (int i3 = 0; i3 < i1; i3 += 1) {
                      w = A[i1][i3];
                      for (int i5 = 0; i5 < i3; i5 += 1) {
                        w -= (A[i1][i5] * A[i5][i3]);
                      }
                      A[i1][i3] = (w / A[i3][i3]);
                    }
                  }
                  for (int i3 = i1; i3 < _PB_N; i3 += 1) {
                    if (i3 + 62 >= _PB_N) {
                      w = A[i1][i3];
                      for (int i5 = 0; i5 < i1; i5 += 1) {
                        w -= (A[i1][i5] * A[i5][i3]);
                      }
                    }
                    A[i1][i3] = w;
                  }
                }
              }
            }
          }
        }
        if (ii1 % 2 == 0) {
          w = A[-_PB_N + 64 * ii1 + 1][32 * ii1];
          for (int ii4 = 1; ii4 <= 2; ii4 += 1) {
            if (ii4 == 1) {
              for (int ii5 = 0; ii5 <= (ii1 / 2) - 1; ii5 += 1) {
                for (int i5 = 64 * ii5; i5 <= min(-_PB_N + 64 * ii1, 64 * ii5 + 63); i5 += 1) {
                  w -= (A[-_PB_N + 64 * ii1 + 1][i5] * A[i5][32 * ii1]);
                }
              }
            } else {
              for (int i1 = -_PB_N + 64 * ii1 + 1; i1 < _PB_N; i1 += 1) {
                if (_PB_N + i1 >= 64 * ii1 + 2) {
                  for (int i3 = 0; i3 < i1; i3 += 1) {
                    w = A[i1][i3];
                    for (int i5 = 0; i5 < i3; i5 += 1) {
                      w -= (A[i1][i5] * A[i5][i3]);
                    }
                    A[i1][i3] = (w / A[i3][i3]);
                  }
                }
                if (_PB_N + i1 >= 64 * ii1 + 2) {
                  for (int i3 = i1; i3 < 32 * ii1; i3 += 1) {
                    w = A[i1][i3];
                    for (int i5 = 0; i5 < i1; i5 += 1) {
                      w -= (A[i1][i5] * A[i5][i3]);
                    }
                    A[i1][i3] = w;
                  }
                }
                for (int i3 = max(32 * ii1, i1); i3 < _PB_N; i3 += 1) {
                  if (_PB_N + i1 + i3 >= 96 * ii1 + 2) {
                    w = A[i1][i3];
                    for (int i5 = 0; i5 < i1; i5 += 1) {
                      w -= (A[i1][i5] * A[i5][i3]);
                    }
                  }
                  A[i1][i3] = w;
                }
              }
            }
          }
        }
      } else if (_PB_N >= 64 * ii1 + 1) {
        for (int ii3 = 0; ii3 < _PB_N / 64; ii3 += 1) {
          w = A[0][64 * ii3];
          for (int i3 = 64 * ii3; i3 <= 64 * ii3 + 63; i3 += 1) {
            if (i3 >= 64 * ii3 + 1) {
              w = A[0][i3];
            }
            A[0][i3] = w;
          }
        }
        if (_PB_N >= 128 && 64 * ii1 + 63 >= _PB_N) {
          w = A[0][64 * ii1];
          for (int i1 = 0; i1 <= -_PB_N + 64 * ii1 + 64; i1 += 1) {
            for (int i3 = 0; i3 < i1; i3 += 1) {
              w = A[i1][i3];
              for (int i5 = 0; i5 < i3; i5 += 1) {
                w -= (A[i1][i5] * A[i5][i3]);
              }
              A[i1][i3] = (w / A[i3][i3]);
            }
            if (i1 >= 1) {
              for (int i3 = i1; i3 < _PB_N; i3 += 1) {
                if (i3 >= 64 * ii1) {
                  w = A[i1][i3];
                  for (int i5 = 0; i5 < i1; i5 += 1) {
                    w -= (A[i1][i5] * A[i5][i3]);
                  }
                  A[i1][i3] = w;
                } else {
                  w = A[i1][i3];
                  if (i3 >= 2) {
                    if (i3 == i1) {
                      w -= (A[i1][0] * A[0][i1]);
                    }
                    for (int i5 = max(0, i1 - i3 + 1); i5 < i1; i5 += 1) {
                      w -= (A[i1][i5] * A[i5][i3]);
                    }
                  } else {
                    w -= (A[1][0] * A[0][1]);
                  }
                  A[i1][i3] = w;
                }
              }
            } else {
              for (int i3 = 64 * ii1; i3 < _PB_N; i3 += 1) {
                if (i3 >= 64 * ii1 + 1) {
                  w = A[0][i3];
                }
                A[0][i3] = w;
              }
            }
          }
        } else if (_PB_N <= 127 && ii1 == 1) {
          w = A[0][64];
          if (_PB_N == 65) {
            A[0][64] = w;
          }
          for (int i1 = max(0, -_PB_N + 66); i1 <= 1; i1 += 1) {
            if (i1 == 1) {
              w = A[1][0];
              A[1][0] = (w / A[0][0]);
            }
            if (i1 == 1) {
              for (int i3 = 1; i3 < _PB_N; i3 += 1) {
                w = A[1][i3];
                w -= (A[1][0] * A[0][i3]);
                A[1][i3] = w;
              }
            } else {
              for (int i3 = 64; i3 < _PB_N; i3 += 1) {
                if (i3 >= 65) {
                  w = A[0][i3];
                }
                A[0][i3] = w;
              }
            }
          }
          for (int i1 = 2; i1 <= -_PB_N + 128; i1 += 1) {
            for (int i3 = 0; i3 <= 1; i3 += 1) {
              w = A[i1][i3];
              if (i3 == 1) {
                w -= (A[i1][0] * A[0][1]);
              }
              A[i1][i3] = (w / A[i3][i3]);
            }
            for (int i3 = 2; i3 < i1; i3 += 1) {
              w = A[i1][i3];
              for (int i5 = 0; i5 < i3; i5 += 1) {
                w -= (A[i1][i5] * A[i5][i3]);
              }
              A[i1][i3] = (w / A[i3][i3]);
            }
            for (int i3 = i1; i3 < _PB_N; i3 += 1) {
              w = A[i1][i3];
              if (i3 == i1) {
                w -= (A[i1][0] * A[0][i1]);
              }
              for (int i5 = max(0, i1 - i3 + 1); i5 < i1; i5 += 1) {
                w -= (A[i1][i5] * A[i5][i3]);
              }
              A[i1][i3] = w;
            }
          }
        }
        if (_PB_N <= 63 && ii1 == 0) {
          w = A[0][0];
          for (int i1 = 0; i1 <= min(_PB_N - 1, -_PB_N + 64); i1 += 1) {
            for (int i3 = 0; i3 < i1; i3 += 1) {
              w = A[i1][i3];
              for (int i5 = 0; i5 < i3; i5 += 1) {
                w -= (A[i1][i5] * A[i5][i3]);
              }
              A[i1][i3] = (w / A[i3][i3]);
            }
            for (int i3 = i1; i3 < _PB_N; i3 += 1) {
              if (_PB_N >= 33 && i1 >= 2) {
                w = A[i1][i3];
              } else if (_PB_N <= 32 && i3 >= 1) {
                w = A[i1][i3];
              } else if (_PB_N >= 33 && i3 >= i1 + 1) {
                w = A[i1][i3];
              } else if (_PB_N >= 33 && i1 == 1) {
                w = A[1][1];
              }
              for (int i5 = 0; i5 < i1; i5 += 1) {
                w -= (A[i1][i5] * A[i5][i3]);
              }
              A[i1][i3] = w;
            }
          }
        }
      } else {
        for (int ii3 = 0; ii3 < _PB_N / 64; ii3 += 1) {
          for (int ii4 = 0; ii4 <= 2; ii4 += 1) {
            if (ii4 <= 1) {
              if (ii4 == 1) {
                if (ii3 >= 1) {
                  w -= (A[1][0] * A[0][64 * ii3]);
                } else {
                  w -= (A[1][0] * A[0][1]);
                }
              } else if (ii3 >= 1) {
                w = A[1][64 * ii3];
              } else {
                w = A[1][1];
              }
            } else if (_PB_N >= 64 * ii3 + 128) {
              for (int i3 = max(1, 64 * ii3); i3 <= 64 * ii3 + 63; i3 += 1) {
                if (i3 >= 2 && i3 >= 64 * ii3 + 1) {
                  w = A[1][i3];
                  w -= (A[1][0] * A[0][i3]);
                }
                A[1][i3] = w;
              }
            } else {
              for (int i1 = 1; i1 <= 64; i1 += 1) {
                if (i1 >= 2) {
                  for (int i3 = 0; i3 < i1; i3 += 1) {
                    w = A[i1][i3];
                    for (int i5 = 0; i5 < i3; i5 += 1) {
                      w -= (A[i1][i5] * A[i5][i3]);
                    }
                    A[i1][i3] = (w / A[i3][i3]);
                  }
                }
                if (i1 >= 2) {
                  for (int i3 = i1; i3 < _PB_N - 64; i3 += 1) {
                    w = A[i1][i3];
                    for (int i5 = 0; i5 < i1; i5 += 1) {
                      w -= (A[i1][i5] * A[i5][i3]);
                    }
                    A[i1][i3] = w;
                  }
                  for (int i3 = _PB_N - 64; i3 < _PB_N; i3 += 1) {
                    w = A[i1][i3];
                    for (int i5 = 0; i5 < i1; i5 += 1) {
                      if ((i3 >= i1 + 1 && i5 >= 1) || (i1 + i3 + 61 >= _PB_N && i3 + i5 >= i1 + 1) || i5 == 0) {
                        w -= (A[i1][i5] * A[i5][i3]);
                      }
                    }
                    A[i1][i3] = w;
                  }
                } else {
                  for (int i3 = _PB_N - 64; i3 < _PB_N; i3 += 1) {
                    if (i3 + 63 >= _PB_N) {
                      w = A[1][i3];
                      w -= (A[1][0] * A[0][i3]);
                    }
                    A[1][i3] = w;
                  }
                }
              }
            }
          }
        }
      }
    }
    for (int ii1 = (_PB_N - 1) / 64; ii1 <= (_PB_N - 1) / 32; ii1 += 1) {
      if (64 * ii1 >= _PB_N) {
        w = b[-_PB_N + 64 * ii1 + 1];
      } else {
        w = b[0];
      }
      if (64 * ii1 + 64 == _PB_N) {
        y[0] = w;
      } else if (32 * ii1 + 32 >= _PB_N) {
        for (int ii2 = max(1, -ii1 + (_PB_N - 2 * ii1 - 1) / 62 + 2); ii2 <= min(2, _PB_N - 32 * ii1); ii2 += 1) {
          if (ii2 == 2) {
            for (int i1 = max(0, -_PB_N + 64 * ii1 + 1); i1 < _PB_N; i1 += 1) {
              if (_PB_N + i1 >= 64 * ii1 + 2 && i1 >= 1) {
                w = b[i1];
              }
              if (_PB_N + i1 >= 64 * ii1 + 2) {
                for (int i3 = 0; i3 < i1; i3 += 1) {
                  w -= (A[i1][i3] * y[i3]);
                }
              }
              y[i1] = w;
            }
          } else {
            for (int ii3 = 0; ii3 <= ii1 - (_PB_N + 63) / 64; ii3 += 1) {
              for (int i3 = 64 * ii3; i3 <= min(-_PB_N + 64 * ii1, 64 * ii3 + 63); i3 += 1) {
                w -= (A[-_PB_N + 64 * ii1 + 1][i3] * y[i3]);
              }
            }
          }
        }
        if (32 * ii1 + 1 == _PB_N) {
          y[_PB_N - 1] = w;
        }
      } else if (_PB_N >= 64 * ii1) {
        int ii2 = _PB_N >= 64 * ii1 + 1 ? 2 : 1;
        if (ii2 == 2) {
          for (int i1 = 0; i1 <= 1; i1 += 1) {
            if (i1 == 1) {
              w = b[1];
              w -= (A[1][0] * y[0]);
              y[1] = w;
            } else {
              y[0] = w;
            }
          }
          for (int i1 = 2; i1 <= -_PB_N + 64 * ii1 + 64; i1 += 1) {
            w = b[i1];
            w -= (A[i1][0] * y[0]);
            for (int i3 = 1; i3 < i1; i3 += 1) {
              w -= (A[i1][i3] * y[i3]);
            }
            y[i1] = w;
          }
        } else {
          w -= (A[1][0] * y[0]);
        }
      }
      if (_PB_N >= 32 * ii1 + 33) {
        for (int ii2 = max(1, _PB_N - 64 * ii1 + 2); ii2 <= 2; ii2 += 1) {
          if (ii2 == 2) {
            for (int i1 = -_PB_N + 64 * ii1 + 1; i1 <= -_PB_N + 64 * ii1 + 64; i1 += 1) {
              if (_PB_N + i1 >= 64 * ii1 + 2) {
                w = b[i1];
                w -= (A[i1][0] * y[0]);
                for (int i3 = 1; i3 < i1; i3 += 1) {
                  w -= (A[i1][i3] * y[i3]);
                }
              }
              y[i1] = w;
            }
          } else {
            for (int ii3 = 0; ii3 <= ii1 - (_PB_N + 63) / 64; ii3 += 1) {
              if (ii3 == 0) {
                w -= (A[-_PB_N + 64 * ii1 + 1][0] * y[0]);
              }
              for (int i3 = max(1, 64 * ii3); i3 <= min(-_PB_N + 64 * ii1, 64 * ii3 + 63); i3 += 1) {
                w -= (A[-_PB_N + 64 * ii1 + 1][i3] * y[i3]);
              }
            }
          }
        }
      }
    }
  }
  for (int ii1 = 0; ii1 <= floord(_PB_N - 1, 64); ii1 += 1) {
    w = y[_PB_N - 64 * ii1 - 1];
    for (int ii2 = max(1, -ii1 + 2); ii2 <= 2; ii2 += 1) {
      if (ii2 == 2) {
        for (int i1 = -_PB_N + 64 * ii1 + 1; i1 <= min(0, -_PB_N + 64 * ii1 + 64); i1 += 1) {
          if (_PB_N + i1 >= 64 * ii1 + 2) {
            w = y[-i1];
            for (int i3 = -i1 + 1; i3 < _PB_N; i3 += 1) {
              w -= (A[-i1][i3] * x[i3]);
            }
          }
          x[-i1] = (w / A[-i1][-i1]);
        }
      } else {
        if ((_PB_N - 63) % 64 == 0) {
          w -= (A[_PB_N - 64 * ii1 - 1][_PB_N - 64 * ii1] * x[_PB_N - 64 * ii1]);
        }
        for (int ii3 = -ii1 + (_PB_N + 1) / 64; ii3 < _PB_N / 64; ii3 += 1) {
          if (_PB_N >= 64 * ii1 + 64 * ii3) {
            w -= (A[_PB_N - 64 * ii1 - 1][_PB_N - 64 * ii1] * x[_PB_N - 64 * ii1]);
          }
          for (int i3 = max(_PB_N - 64 * ii1 + 1, 64 * ii3); i3 <= 64 * ii3 + 63; i3 += 1) {
            w -= (A[_PB_N - 64 * ii1 - 1][i3] * x[i3]);
          }
        }
        if (_PB_N % 64 >= 1) {
          for (int i3 = -((_PB_N - 1) % 64) + _PB_N - 1; i3 < _PB_N; i3 += 1) {
            w -= (A[_PB_N - 64 * ii1 - 1][i3] * x[i3]);
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
  POLYBENCH_2D_ARRAY_DECL(A, DATA_TYPE, N, N, n, n);
  POLYBENCH_1D_ARRAY_DECL(b, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(x, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(y, DATA_TYPE, N, n);


  /* Initialize array(s). */
  init_array (n,
	      POLYBENCH_ARRAY(A),
	      POLYBENCH_ARRAY(b),
	      POLYBENCH_ARRAY(x),
	      POLYBENCH_ARRAY(y));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_ludcmp (n,
		 POLYBENCH_ARRAY(A),
		 POLYBENCH_ARRAY(b),
		 POLYBENCH_ARRAY(x),
		 POLYBENCH_ARRAY(y));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(x)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(b);
  POLYBENCH_FREE_ARRAY(x);
  POLYBENCH_FREE_ARRAY(y);

  return 0;
}
