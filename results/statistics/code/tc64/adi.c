/**
 * This version is stamped on Apr. 14, 2015
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* adi.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "adi.h"


/* Array initialization. */
static
void init_array (int n,
		 DATA_TYPE POLYBENCH_2D(u,N,N,n,n))
{
  int i, j;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      {
	u[i][j] =  (DATA_TYPE)(i + n-j) / n;
      }
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_2D(u,N,N,n,n))

{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("u");
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
      if ((i * n + j) % 20 == 0) fprintf(POLYBENCH_DUMP_TARGET, "\n");
      fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, u[i][j]);
    }
  POLYBENCH_DUMP_END("u");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
/* Based on a Fortran code fragment from Figure 5 of
 * "Automatic Data and Computation Decomposition on Distributed Memory Parallel Computers"
 * by Peizong Lee and Zvi Meir Kedem, TOPLAS, 2002
 */
static
void kernel_adi(int tsteps, int n,
		DATA_TYPE POLYBENCH_2D(u,N,N,n,n),
		DATA_TYPE POLYBENCH_2D(v,N,N,n,n),
		DATA_TYPE POLYBENCH_2D(p,N,N,n,n),
		DATA_TYPE POLYBENCH_2D(q,N,N,n,n))
{
  int t, i, j;
  DATA_TYPE DX, DY, DT;
  DATA_TYPE B1, B2;
  DATA_TYPE mul1, mul2;
  DATA_TYPE a, b, c, d, e, f;


  DX = SCALAR_VAL(1.0)/(DATA_TYPE)_PB_N;
  DY = SCALAR_VAL(1.0)/(DATA_TYPE)_PB_N;
  DT = SCALAR_VAL(1.0)/(DATA_TYPE)_PB_TSTEPS;
  B1 = SCALAR_VAL(2.0);
  B2 = SCALAR_VAL(1.0);
  mul1 = B1 * DT / (DX * DX);
  mul2 = B2 * DT / (DY * DY);

  a = -mul1 /  SCALAR_VAL(2.0);
  b = SCALAR_VAL(1.0)+mul1;
  c = a;
  d = -mul2 / SCALAR_VAL(2.0);
  e = SCALAR_VAL(1.0)+mul2;
  f = d;

/* TC Optimizing Compiler 0.2.26 */
/* ./tc ../examples/polybench/adi.scop.c --correction-tiling --lex-scheduling --serial-codegen -b 64 --debug --floyd-warshall-tc */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#pragma scop
for (int ii0 = 0; ii0 <= floord(_PB_TSTEPS - 1, 64); ii0 += 1) {
  for (int ii1 = 0; ii1 <= 1; ii1 += 1) {
    for (int ii2 = 0; ii2 <= floord(_PB_N - 3, 64); ii2 += 1) {
      for (int ii3 = 0; ii3 <= 4; ii3 += 1) {
        if (ii1 == 0 && ii3 <= 2) {
          if (ii3 == 2) {
            for (int i2 = 64 * ii2 + 1; i2 <= min(_PB_N - 2, 64 * ii2 + 64); i2 += 1) {
              q[i2][0] = v[0][i2];
            }
          } else if (ii3 == 1) {
            for (int i2 = 64 * ii2 + 1; i2 <= min(_PB_N - 2, 64 * ii2 + 64); i2 += 1) {
              p[i2][0] = SCALAR_VAL(0.0);
            }
          } else {
            for (int i2 = 64 * ii2 + 1; i2 <= min(_PB_N - 2, 64 * ii2 + 64); i2 += 1) {
              v[0][i2] = SCALAR_VAL(1.0);
            }
          }
        } else if (ii1 == 1 && ii3 == 4) {
          for (int i2 = 64 * ii2 + 1; i2 <= min(_PB_N - 2, 64 * ii2 + 64); i2 += 1) {
            u[i2][_PB_N - 1] = SCALAR_VAL(1.0);
          }
        } else if (ii1 == 0 && ii3 == 4) {
          for (int i2 = 64 * ii2 + 1; i2 <= min(_PB_N - 2, 64 * ii2 + 64); i2 += 1) {
            v[_PB_N - 1][i2] = SCALAR_VAL(1.0);
          }
        } else if (ii1 == 1 && ii3 == 3) {
          for (int ii4 = (_PB_N - 1) / 64; ii4 <= (_PB_N - 2) / 32; ii4 += 1) {
            for (int i2 = 64 * ii2 + 1; i2 <= min(_PB_N - 2, 64 * ii2 + 64); i2 += 1) {
              for (int i4 = max(1, -_PB_N + 64 * ii4 + 2); i4 <= min(_PB_N - 2, -_PB_N + 64 * ii4 + 65); i4 += 1) {
                p[i2][i4] = ((-f) / ((d * p[i2][i4 - 1]) + e));
              }
            }
            for (int i2 = 64 * ii2 + 1; i2 <= min(_PB_N - 2, 64 * ii2 + 64); i2 += 1) {
              for (int i4 = max(1, -_PB_N + 64 * ii4 + 2); i4 <= min(_PB_N - 2, -_PB_N + 64 * ii4 + 65); i4 += 1) {
                q[i2][i4] = ((((((-a) * v[i2 - 1][i4]) + ((SCALAR_VAL(1.0) + (SCALAR_VAL(2.0) * a)) * v[i2][i4])) - (c * v[i2 + 1][i4])) - (d * q[i2][i4 - 1])) / ((d * p[i2][i4 - 1]) + e));
              }
            }
            if (_PB_TSTEPS >= 64 * ii0 + 2 && 64 * ii2 + 66 >= _PB_N) {
              for (int i2 = max(1, -_PB_N + 64 * ii4 + 2); i2 <= min(_PB_N - 2, -_PB_N + 64 * ii4 + 65); i2 += 1) {
                if (ii2 == 0) {
                  v[0][i2] = SCALAR_VAL(1.0);
                  if (_PB_N >= 64 * ii4 + 1) {
                    p[i2][0] = SCALAR_VAL(0.0);
                    if (_PB_N >= 32 * ii4 + 34) {
                      q[i2][0] = v[0][i2];
                    }
                  }
                  if (32 * ii4 + 33 >= _PB_N) {
                    q[i2][0] = v[0][i2];
                  }
                }
                v[_PB_N - 1][i2] = SCALAR_VAL(1.0);
              }
              if (_PB_N >= 64 * ii4 + 1) {
                for (int i2 = max(64 * ii2 + 1, -_PB_N + 64 * ii4 + 66); i2 < _PB_N - 1; i2 += 1) {
                  p[i2][0] = SCALAR_VAL(0.0);
                  if (ii2 >= 1) {
                    q[i2][0] = v[0][i2];
                  }
                }
              }
            } else if (_PB_TSTEPS >= 64 * ii0 + 2 && _PB_N >= 64 * ii2 + 67 && _PB_N >= 64 * ii4 + 1) {
              for (int i2 = 64 * ii2 + 1; i2 <= 64 * ii2 + 64; i2 += 1) {
                if (ii2 == 0 && 64 * ii4 + 65 >= _PB_N + i2) {
                  v[0][i2] = SCALAR_VAL(1.0);
                }
                p[i2][0] = SCALAR_VAL(0.0);
                if (ii2 >= 1) {
                  q[i2][0] = v[0][i2];
                } else if (64 * ii4 + 65 >= _PB_N + i2) {
                  q[i2][0] = v[0][i2];
                }
              }
            } else if (_PB_N >= 67 && _PB_TSTEPS >= 64 * ii0 + 2 && ii2 == 0 && 64 * ii4 >= _PB_N) {
              for (int i2 = -_PB_N + 64 * ii4 + 2; i2 <= min(_PB_N - 2, -_PB_N + 64 * ii4 + 65); i2 += 1) {
                v[0][i2] = SCALAR_VAL(1.0);
                if (i2 <= 64) {
                  q[i2][0] = v[0][i2];
                }
              }
            }
          }
        } else if (ii1 == 0 && ii3 == 3) {
          for (int ii4 = (_PB_N - 1) / 64; ii4 <= (_PB_N - 2) / 32; ii4 += 1) {
            for (int i2 = 64 * ii2 + 1; i2 <= min(_PB_N - 2, 64 * ii2 + 64); i2 += 1) {
              for (int i4 = max(1, -_PB_N + 64 * ii4 + 2); i4 <= min(_PB_N - 2, -_PB_N + 64 * ii4 + 65); i4 += 1) {
                p[i2][i4] = ((-c) / ((a * p[i2][i4 - 1]) + b));
              }
            }
            for (int i2 = 64 * ii2 + 1; i2 <= min(_PB_N - 2, 64 * ii2 + 64); i2 += 1) {
              for (int i4 = max(1, -_PB_N + 64 * ii4 + 2); i4 <= min(_PB_N - 2, -_PB_N + 64 * ii4 + 65); i4 += 1) {
                q[i2][i4] = ((((((-d) * u[i4][i2 - 1]) + ((SCALAR_VAL(1.0) + (SCALAR_VAL(2.0) * d)) * u[i4][i2])) - (f * u[i4][i2 + 1])) - (a * q[i2][i4 - 1])) / ((a * p[i2][i4 - 1]) + b));
              }
            }
          }
        } else if (ii3 == 2) {
          for (int i2 = 64 * ii2 + 1; i2 <= min(_PB_N - 2, 64 * ii2 + 64); i2 += 1) {
            q[i2][0] = u[i2][0];
          }
        } else if (ii3 == 1) {
          for (int i2 = 64 * ii2 + 1; i2 <= min(_PB_N - 2, 64 * ii2 + 64); i2 += 1) {
            p[i2][0] = SCALAR_VAL(0.0);
          }
        } else {
          for (int i2 = 64 * ii2 + 1; i2 <= min(_PB_N - 2, 64 * ii2 + 64); i2 += 1) {
            u[i2][0] = SCALAR_VAL(1.0);
          }
        }
      }
      if (64 * ii0 + 1 == _PB_TSTEPS && ii1 == 1 && _PB_N >= 64 * ii2 + 67) {
        for (int ii4 = 0; ii4 <= (_PB_N - 3) / 64; ii4 += 1) {
          for (int i2 = 64 * ii2 + 1; i2 <= 64 * ii2 + 64; i2 += 1) {
            for (int i4 = -_PB_N + 64 * ii4 + 2; i4 <= min(-1, -_PB_N + 64 * ii4 + 65); i4 += 1) {
              u[i2][-i4] = ((p[i2][-i4] * u[i2][-i4 + 1]) + q[i2][-i4]);
            }
          }
        }
      } else if (64 * ii0 + 1 == _PB_TSTEPS && ii1 == 1 && 64 * ii2 + 66 >= _PB_N) {
        for (int ii4 = 0; ii4 <= ii2; ii4 += 1) {
          for (int i2 = 64 * ii2 + 1; i2 < _PB_N - 1; i2 += 1) {
            for (int i4 = -_PB_N + 64 * ii4 + 2; i4 <= min(-1, -_PB_N + 64 * ii4 + 65); i4 += 1) {
              u[i2][-i4] = ((p[i2][-i4] * u[i2][-i4 + 1]) + q[i2][-i4]);
            }
          }
        }
      } else if (ii1 == 0) {
        for (int ii4 = 0; ii4 <= (_PB_N - 3) / 64; ii4 += 1) {
          for (int i2 = 64 * ii2 + 1; i2 <= min(_PB_N - 2, 64 * ii2 + 64); i2 += 1) {
            for (int i4 = -_PB_N + 64 * ii4 + 2; i4 <= min(-1, -_PB_N + 64 * ii4 + 65); i4 += 1) {
              v[-i4][i2] = ((p[i2][-i4] * v[-i4 + 1][i2]) + q[i2][-i4]);
            }
          }
        }
      } else {
        if (_PB_N >= 67 && _PB_N >= 64 * ii2 + 65 && 64 * ii2 + 66 >= _PB_N) {
          for (int i2 = 64 * ii2 + 1; i2 < _PB_N - 1; i2 += 1) {
            for (int i4 = -_PB_N + 2; i4 <= -_PB_N + 65; i4 += 1) {
              u[i2][-i4] = ((p[i2][-i4] * u[i2][-i4 + 1]) + q[i2][-i4]);
            }
          }
        } else if (_PB_N >= 64 * ii2 + 67) {
          for (int ii4 = 0; ii4 < min((_PB_N - 3) / 64, -ii2 + (_PB_N - 1) / 64); ii4 += 1) {
            for (int i2 = 64 * ii2 + 1; i2 <= 64 * ii2 + 64; i2 += 1) {
              for (int i4 = -_PB_N + 64 * ii4 + 2; i4 <= -_PB_N + 64 * ii4 + 65; i4 += 1) {
                u[i2][-i4] = ((p[i2][-i4] * u[i2][-i4 + 1]) + q[i2][-i4]);
              }
            }
          }
        }
        for (int ii4 = -ii2 + (_PB_N - 1) / 64; ii4 < (_PB_N - 3) / 64; ii4 += 1) {
          for (int i2 = 64 * ii2 + 1; i2 <= min(_PB_N - 2, 64 * ii2 + 64); i2 += 1) {
            for (int i4 = -_PB_N + 64 * ii4 + 2; i4 <= -_PB_N + 64 * ii4 + 65; i4 += 1) {
              u[i2][-i4] = ((p[i2][-i4] * u[i2][-i4 + 1]) + q[i2][-i4]);
            }
          }
          for (int i2 = _PB_N - 64 * ii4 - 64; i2 <= min(64 * ii2, _PB_N - 64 * ii4 - 1); i2 += 1) {
            for (int i4 = 64 * ii2 + 1; i4 <= min(_PB_N - 2, 64 * ii2 + 64); i4 += 1) {
              q[i2][i4] = ((((((-d) * u[i4][i2 - 1]) + ((SCALAR_VAL(1.0) + (SCALAR_VAL(2.0) * d)) * u[i4][i2])) - (f * u[i4][i2 + 1])) - (a * q[i2][i4 - 1])) / ((a * p[i2][i4 - 1]) + b));
            }
            if (64 * ii2 + 66 >= _PB_N) {
              for (int i4 = -_PB_N + 2; i4 < 0; i4 += 1) {
                v[-i4][i2] = ((p[i2][-i4] * v[-i4 + 1][i2]) + q[i2][-i4]);
              }
            }
          }
          if (64 * ii2 + 66 >= _PB_N) {
            for (int i2 = _PB_N - 64 * ii4 - 64; i2 <= min(64 * ii2, _PB_N - 64 * ii4 - 1); i2 += 1) {
              for (int i4 = 1; i4 < _PB_N - 1; i4 += 1) {
                p[i2][i4] = ((-f) / ((d * p[i2][i4 - 1]) + e));
              }
            }
          }
        }
        if (_PB_N >= 64 * ii2 + 67) {
          for (int i2 = 64 * ii2 + 1; i2 <= 64 * ii2 + 64; i2 += 1) {
            for (int i4 = -((_PB_N - 3) % 64) - 1; i4 < 0; i4 += 1) {
              u[i2][-i4] = ((p[i2][-i4] * u[i2][-i4 + 1]) + q[i2][-i4]);
            }
          }
          for (int i2 = 1; i2 <= min(64 * ii2, ((_PB_N - 3) % 64) + 2); i2 += 1) {
            for (int i4 = 64 * ii2 + 1; i4 <= 64 * ii2 + 64; i4 += 1) {
              q[i2][i4] = ((((((-d) * u[i4][i2 - 1]) + ((SCALAR_VAL(1.0) + (SCALAR_VAL(2.0) * d)) * u[i4][i2])) - (f * u[i4][i2 + 1])) - (a * q[i2][i4 - 1])) / ((a * p[i2][i4 - 1]) + b));
            }
          }
          for (int i2 = 64 * ii2 + 1; i2 <= 64 * ii2 + 64; i2 += 1) {
            for (int i4 = 1; i4 < _PB_N - 1; i4 += 1) {
              p[i2][i4] = ((-c) / ((a * p[i2][i4 - 1]) + b));
              if (64 * ii2 + 64 >= i4) {
                q[i2][i4] = ((((((-d) * u[i4][i2 - 1]) + ((SCALAR_VAL(1.0) + (SCALAR_VAL(2.0) * d)) * u[i4][i2])) - (f * u[i4][i2 + 1])) - (a * q[i2][i4 - 1])) / ((a * p[i2][i4 - 1]) + b));
              }
            }
          }
          for (int i2 = 64 * ii2 + 1; i2 <= 64 * ii2 + 64; i2 += 1) {
            u[i2][0] = SCALAR_VAL(1.0);
            p[i2][0] = SCALAR_VAL(0.0);
            q[i2][0] = u[i2][0];
          }
        } else {
          for (int i0 = 64 * ii0 + 1; i0 <= min(_PB_TSTEPS, 64 * ii0 + 64); i0 += 1) {
            if (i0 >= 64 * ii0 + 2) {
              if (i0 == 64 * ii0 + 2) {
                for (int i2 = 1; i2 <= min(_PB_N - 64 * ii2 - 1, 64 * ii2); i2 += 1) {
                  for (int i4 = 64 * ii2 + 1; i4 < _PB_N - 1; i4 += 1) {
                    q[i2][i4] = ((((((-d) * u[i4][i2 - 1]) + ((SCALAR_VAL(1.0) + (SCALAR_VAL(2.0) * d)) * u[i4][i2])) - (f * u[i4][i2 + 1])) - (a * q[i2][i4 - 1])) / ((a * p[i2][i4 - 1]) + b));
                  }
                  for (int i4 = -_PB_N + 2; i4 < 0; i4 += 1) {
                    v[-i4][i2] = ((p[i2][-i4] * v[-i4 + 1][i2]) + q[i2][-i4]);
                  }
                }
              } else {
                for (int i2 = 1; i2 <= 64 * ii2; i2 += 1) {
                  v[0][i2] = SCALAR_VAL(1.0);
                  p[i2][0] = SCALAR_VAL(0.0);
                  q[i2][0] = v[0][i2];
                  for (int i4 = 1; i4 < _PB_N - 1; i4 += 1) {
                    p[i2][i4] = ((-c) / ((a * p[i2][i4 - 1]) + b));
                    q[i2][i4] = ((((((-d) * u[i4][i2 - 1]) + ((SCALAR_VAL(1.0) + (SCALAR_VAL(2.0) * d)) * u[i4][i2])) - (f * u[i4][i2 + 1])) - (a * q[i2][i4 - 1])) / ((a * p[i2][i4 - 1]) + b));
                  }
                  v[_PB_N - 1][i2] = SCALAR_VAL(1.0);
                  for (int i4 = -_PB_N + 2; i4 < 0; i4 += 1) {
                    v[-i4][i2] = ((p[i2][-i4] * v[-i4 + 1][i2]) + q[i2][-i4]);
                  }
                }
              }
              for (int i2 = 64 * ii2 + 1; i2 < _PB_N - 1; i2 += 1) {
                if (i0 >= 64 * ii0 + 3) {
                  v[0][i2] = SCALAR_VAL(1.0);
                  p[i2][0] = SCALAR_VAL(0.0);
                  q[i2][0] = v[0][i2];
                }
                for (int i4 = 1; i4 < _PB_N - 1; i4 += 1) {
                  if (i0 >= 64 * ii0 + 3) {
                    p[i2][i4] = ((-c) / ((a * p[i2][i4 - 1]) + b));
                  } else {
                    p[i2][i4] = ((-c) / ((a * p[i2][i4 - 1]) + b));
                  }
                  q[i2][i4] = ((((((-d) * u[i4][i2 - 1]) + ((SCALAR_VAL(1.0) + (SCALAR_VAL(2.0) * d)) * u[i4][i2])) - (f * u[i4][i2 + 1])) - (a * q[i2][i4 - 1])) / ((a * p[i2][i4 - 1]) + b));
                }
                if (i0 >= 64 * ii0 + 3) {
                  v[_PB_N - 1][i2] = SCALAR_VAL(1.0);
                }
                for (int i4 = -_PB_N + 2; i4 < 0; i4 += 1) {
                  v[-i4][i2] = ((p[i2][-i4] * v[-i4 + 1][i2]) + q[i2][-i4]);
                }
              }
            }
            if (i0 >= 64 * ii0 + 2) {
              for (int i2 = 1; i2 <= 64 * ii2; i2 += 1) {
                if (i0 >= 64 * ii0 + 3) {
                  u[i2][0] = SCALAR_VAL(1.0);
                  p[i2][0] = SCALAR_VAL(0.0);
                  q[i2][0] = u[i2][0];
                }
                for (int i4 = 1; i4 < _PB_N - 1; i4 += 1) {
                  if (i0 >= 64 * ii0 + 3 && 64 * ii2 + i2 >= _PB_N) {
                    p[i2][i4] = ((-f) / ((d * p[i2][i4 - 1]) + e));
                  } else if (_PB_N >= 64 * ii2 + i2 + 1) {
                    p[i2][i4] = ((-f) / ((d * p[i2][i4 - 1]) + e));
                  }
                  q[i2][i4] = ((((((-a) * v[i2 - 1][i4]) + ((SCALAR_VAL(1.0) + (SCALAR_VAL(2.0) * a)) * v[i2][i4])) - (c * v[i2 + 1][i4])) - (d * q[i2][i4 - 1])) / ((d * p[i2][i4 - 1]) + e));
                }
                u[i2][_PB_N - 1] = SCALAR_VAL(1.0);
                for (int i4 = -_PB_N + 2; i4 < 0; i4 += 1) {
                  u[i2][-i4] = ((p[i2][-i4] * u[i2][-i4 + 1]) + q[i2][-i4]);
                }
              }
            }
            for (int i2 = 64 * ii2 + 1; i2 < _PB_N - 1; i2 += 1) {
              if (i0 >= 64 * ii0 + 2) {
                u[i2][0] = SCALAR_VAL(1.0);
                p[i2][0] = SCALAR_VAL(0.0);
                q[i2][0] = u[i2][0];
                for (int i4 = 1; i4 < _PB_N - 1; i4 += 1) {
                  p[i2][i4] = ((-f) / ((d * p[i2][i4 - 1]) + e));
                  q[i2][i4] = ((((((-a) * v[i2 - 1][i4]) + ((SCALAR_VAL(1.0) + (SCALAR_VAL(2.0) * a)) * v[i2][i4])) - (c * v[i2 + 1][i4])) - (d * q[i2][i4 - 1])) / ((d * p[i2][i4 - 1]) + e));
                }
                u[i2][_PB_N - 1] = SCALAR_VAL(1.0);
              }
              if (i0 >= 64 * ii0 + 2) {
                for (int i4 = -_PB_N + 2; i4 <= -_PB_N + 64 * ii2 + 1; i4 += 1) {
                  u[i2][-i4] = ((p[i2][-i4] * u[i2][-i4 + 1]) + q[i2][-i4]);
                }
              }
              for (int i4 = -_PB_N + 64 * ii2 + 2; i4 < 0; i4 += 1) {
                u[i2][-i4] = ((p[i2][-i4] * u[i2][-i4 + 1]) + q[i2][-i4]);
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
  POLYBENCH_2D_ARRAY_DECL(u, DATA_TYPE, N, N, n, n);
  POLYBENCH_2D_ARRAY_DECL(v, DATA_TYPE, N, N, n, n);
  POLYBENCH_2D_ARRAY_DECL(p, DATA_TYPE, N, N, n, n);
  POLYBENCH_2D_ARRAY_DECL(q, DATA_TYPE, N, N, n, n);


  /* Initialize array(s). */
  init_array (n, POLYBENCH_ARRAY(u));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_adi (tsteps, n, POLYBENCH_ARRAY(u), POLYBENCH_ARRAY(v), POLYBENCH_ARRAY(p), POLYBENCH_ARRAY(q));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(u)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(u);
  POLYBENCH_FREE_ARRAY(v);
  POLYBENCH_FREE_ARRAY(p);
  POLYBENCH_FREE_ARRAY(q);

  return 0;
}
