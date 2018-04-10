/**
 * This version is stamped on Apr. 14, 2015
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* gramschmidt.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "gramschmidt.h"


/* Array initialization. */
static
void init_array(int m, int n,
		DATA_TYPE POLYBENCH_2D(A,M,N,m,n),
		DATA_TYPE POLYBENCH_2D(R,N,N,n,n),
		DATA_TYPE POLYBENCH_2D(Q,M,N,m,n))
{
  int i, j;

  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++) {
      A[i][j] = (((DATA_TYPE) ((i*j) % m) / m )*100) + 10;
      Q[i][j] = 0.0;
    }
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      R[i][j] = 0.0;
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int m, int n,
		 DATA_TYPE POLYBENCH_2D(A,M,N,m,n),
		 DATA_TYPE POLYBENCH_2D(R,N,N,n,n),
		 DATA_TYPE POLYBENCH_2D(Q,M,N,m,n))
{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("R");
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
	if ((i*n+j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
	fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, R[i][j]);
    }
  POLYBENCH_DUMP_END("R");

  POLYBENCH_DUMP_BEGIN("Q");
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++) {
	if ((i*n+j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
	fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, Q[i][j]);
    }
  POLYBENCH_DUMP_END("Q");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
/* QR Decomposition with Modified Gram Schmidt:
 http://www.inf.ethz.ch/personal/gander/ */
static
void kernel_gramschmidt(int m, int n,
			DATA_TYPE POLYBENCH_2D(A,M,N,m,n),
			DATA_TYPE POLYBENCH_2D(R,N,N,n,n),
			DATA_TYPE POLYBENCH_2D(Q,M,N,m,n))
{
  int i, j, k;

  DATA_TYPE nrm;

/* TC Optimizing Compiler 0.2.26 */
/* ./tc ../examples/polybench/gramschmidt.scop.c --correction-tiling --lex-scheduling --serial-codegen -b 32 --debug --isl-union-map-tc */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#pragma scop
for (int ii0 = 0; ii0 <= floord(_PB_N - 1, 32); ii0 += 1) {
  nrm = SCALAR_VAL(0.0);
  if (_PB_N == 1 && _PB_M <= 31 && ii0 == 0) {
    for (int i2 = 0; i2 < _PB_M; i2 += 1) {
      nrm += (A[i2][0] * A[i2][0]);
    }
  } else if (_PB_N >= 2 && _PB_M >= _PB_N + 1 && _PB_M <= 31 && ii0 == 0) {
    for (int i2 = 0; i2 < _PB_M; i2 += 1) {
      nrm += (A[i2][0] * A[i2][0]);
    }
  } else if (_PB_N >= 2 && _PB_N >= _PB_M) {
    for (int ii2 = 0; ii2 <= floord(_PB_M - 1, 32); ii2 += 1) {
      if (32 * ii0 + 31 >= _PB_N) {
        for (int i2 = 32 * ii2; i2 <= min(_PB_M - 1, 32 * ii2 + 31); i2 += 1) {
          nrm += (A[i2][32 * ii0] * A[i2][32 * ii0]);
        }
      } else if (_PB_N >= 32 * ii2 + 32) {
        for (int i2 = 32 * ii2; i2 <= min(_PB_M - 1, 32 * ii2 + 31); i2 += 1) {
          nrm += (A[i2][32 * ii0] * A[i2][32 * ii0]);
        }
      } else {
        for (int i2 = 32 * ii2; i2 < _PB_M; i2 += 1) {
          nrm += (A[i2][32 * ii0] * A[i2][32 * ii0]);
        }
      }
    }
  } else {
    for (int ii2 = 0; ii2 <= (_PB_M - 1) / 32; ii2 += 1) {
      for (int i2 = 32 * ii2; i2 <= min(_PB_M - 1, 32 * ii2 + 31); i2 += 1) {
        nrm += (A[i2][32 * ii0] * A[i2][32 * ii0]);
      }
    }
  }
  if (_PB_M >= 1) {
    for (int i0 = 32 * ii0; i0 <= min(_PB_N - 1, 32 * ii0 + 1); i0 += 1) {
      if (i0 == 32 * ii0 + 1) {
        nrm = SCALAR_VAL(0.0);
      } else {
        R[32 * ii0][32 * ii0] = SQRT_FUN(nrm);
      }
    }
  } else {
    for (int i0 = 32 * ii0; i0 <= min(_PB_N - 1, 32 * ii0 + 31); i0 += 1) {
      if (i0 >= 32 * ii0 + 1) {
        nrm = SCALAR_VAL(0.0);
      }
      R[i0][i0] = SQRT_FUN(nrm);
    }
  }
  if (_PB_N >= _PB_M && _PB_M + 31 >= _PB_N && _PB_N >= 32 * ii0 + 32 && _PB_M % 32 == 0) {
    for (int ii2 = 0; ii2 < _PB_M / 32; ii2 += 1) {
      for (int i2 = 32 * ii2; i2 <= 32 * ii2 + 31; i2 += 1) {
        Q[i2][32 * ii0] = (A[i2][32 * ii0] / R[32 * ii0][32 * ii0]);
      }
    }
  } else if (_PB_N >= 32 * ii0 + 32 && (_PB_M % 32) + _PB_N >= _PB_M + 32) {
    for (int ii2 = 0; ii2 <= floord(_PB_M - 1, 32); ii2 += 1) {
      for (int i2 = 32 * ii2; i2 <= min(_PB_M - 1, 32 * ii2 + 31); i2 += 1) {
        Q[i2][32 * ii0] = (A[i2][32 * ii0] / R[32 * ii0][32 * ii0]);
      }
    }
    if (_PB_M == 0) {
      for (int ii2 = ii0; ii2 <= (_PB_N - 1) / 32; ii2 += 1) {
        if (_PB_N >= 32 * ii2 + 32) {
          for (int i0 = 32 * ii0; i0 <= min(32 * ii0 + 31, 32 * ii2 + 30); i0 += 1) {
            for (int i2 = max(32 * ii2, i0 + 1); i2 <= 32 * ii2 + 31; i2 += 1) {
              R[i0][i2] = SCALAR_VAL(0.0);
            }
          }
        } else {
          for (int i0 = 32 * ii0; i0 <= 32 * ii0 + 31; i0 += 1) {
            for (int i2 = 32 * ii2; i2 < _PB_N; i2 += 1) {
              R[i0][i2] = SCALAR_VAL(0.0);
            }
          }
        }
      }
    }
  } else {
    if (_PB_N >= 2 && _PB_N >= _PB_M && 32 * ii0 + 31 >= _PB_N) {
      for (int ii2 = 0; ii2 <= floord(_PB_M - 1, 32); ii2 += 1) {
        for (int i2 = 32 * ii2; i2 <= min(_PB_M - 1, 32 * ii2 + 31); i2 += 1) {
          Q[i2][32 * ii0] = (A[i2][32 * ii0] / R[32 * ii0][32 * ii0]);
        }
      }
      if (_PB_M >= 1 && 32 * ii0 + 2 == _PB_N) {
        R[_PB_N - 2][_PB_N - 1] = SCALAR_VAL(0.0);
        for (int ii4 = 0; ii4 <= floord(_PB_M - 1, 32); ii4 += 1) {
          for (int i4 = 32 * ii4; i4 <= min(_PB_M - 1, 32 * ii4 + 31); i4 += 1) {
            R[_PB_N - 2][_PB_N - 1] += (Q[i4][_PB_N - 2] * A[i4][_PB_N - 1]);
          }
        }
        for (int ii4 = 0; ii4 <= floord(_PB_M - 1, 32); ii4 += 1) {
          for (int i4 = 32 * ii4; i4 <= min(_PB_M - 1, 32 * ii4 + 31); i4 += 1) {
            A[i4][_PB_N - 1] = (A[i4][_PB_N - 1] - (Q[i4][_PB_N - 2] * R[_PB_N - 2][_PB_N - 1]));
          }
          for (int i2 = 32 * ii4; i2 <= min(_PB_M - 1, 32 * ii4 + 31); i2 += 1) {
            nrm += (A[i2][_PB_N - 1] * A[i2][_PB_N - 1]);
          }
        }
        R[_PB_N - 1][_PB_N - 1] = SQRT_FUN(nrm);
        for (int i2 = 0; i2 < _PB_M; i2 += 1) {
          Q[i2][_PB_N - 1] = (A[i2][_PB_N - 1] / R[_PB_N - 1][_PB_N - 1]);
        }
      }
    } else if (_PB_N >= _PB_M && _PB_N >= 32 * ii0 + 32 && _PB_M % 32 >= 1 && _PB_M + 31 >= (_PB_M % 32) + _PB_N) {
      for (int ii2 = 0; ii2 < _PB_N / 32; ii2 += 1) {
        for (int i2 = 32 * ii2; i2 <= 32 * ii2 + 31; i2 += 1) {
          Q[i2][32 * ii0] = (A[i2][32 * ii0] / R[32 * ii0][32 * ii0]);
        }
      }
      for (int i2 = ((31 * _PB_N + 31) % 32) + _PB_N - 31; i2 < _PB_M; i2 += 1) {
        Q[i2][32 * ii0] = (A[i2][32 * ii0] / R[32 * ii0][32 * ii0]);
      }
    } else if (_PB_N == 1 && ii0 == 0 && _PB_M % 32 == 0) {
      for (int ii2 = 0; ii2 < _PB_M / 32; ii2 += 1) {
        for (int i2 = 32 * ii2; i2 <= 32 * ii2 + 31; i2 += 1) {
          Q[i2][0] = (A[i2][0] / R[0][0]);
        }
      }
    } else if (_PB_N >= 2 && _PB_M >= _PB_N + 1 && 32 * ii0 + 31 >= _PB_N && _PB_M % 32 == 0) {
      for (int ii2 = 0; ii2 < _PB_M / 32; ii2 += 1) {
        for (int i2 = 32 * ii2; i2 <= 32 * ii2 + 31; i2 += 1) {
          Q[i2][32 * ii0] = (A[i2][32 * ii0] / R[32 * ii0][32 * ii0]);
        }
      }
    } else if (_PB_M >= _PB_N + 1 && _PB_N >= 32 * ii0 + 32) {
      for (int ii2 = 0; ii2 <= (_PB_M - 1) / 32; ii2 += 1) {
        for (int i2 = 32 * ii2; i2 <= min(_PB_M - 1, 32 * ii2 + 31); i2 += 1) {
          Q[i2][32 * ii0] = (A[i2][32 * ii0] / R[32 * ii0][32 * ii0]);
        }
      }
    } else if (_PB_N >= 2 && _PB_M >= _PB_N + 1) {
      for (int ii2 = 0; ii2 < _PB_M / 32; ii2 += 1) {
        for (int i2 = 32 * ii2; i2 <= 32 * ii2 + 31; i2 += 1) {
          Q[i2][32 * ii0] = (A[i2][32 * ii0] / R[32 * ii0][32 * ii0]);
        }
      }
      for (int i2 = -((_PB_M - 1) % 32) + _PB_M - 1; i2 < _PB_M; i2 += 1) {
        Q[i2][32 * ii0] = (A[i2][32 * ii0] / R[32 * ii0][32 * ii0]);
      }
    } else {
      for (int ii2 = 0; ii2 < _PB_M / 32; ii2 += 1) {
        for (int i2 = 32 * ii2; i2 <= 32 * ii2 + 31; i2 += 1) {
          Q[i2][0] = (A[i2][0] / R[0][0]);
        }
      }
      for (int i2 = -((_PB_M + 31) % 32) + _PB_M - 1; i2 < _PB_M; i2 += 1) {
        Q[i2][0] = (A[i2][0] / R[0][0]);
      }
    }
    if (_PB_M == 0 && 32 * ii0 + 31 >= _PB_N) {
      for (int i0 = 32 * ii0; i0 < _PB_N - 1; i0 += 1) {
        for (int i2 = i0 + 1; i2 < _PB_N; i2 += 1) {
          R[i0][i2] = SCALAR_VAL(0.0);
        }
      }
    }
  }
  if (_PB_M >= 1 && _PB_N >= _PB_M && _PB_N >= 32 * ii0 + 3) {
    for (int ii2 = ii0; ii2 <= (_PB_N - 1) / 32; ii2 += 1) {
      if (_PB_N >= 32 * ii2 + 32) {
        for (int i0 = 32 * ii0; i0 <= min(32 * ii0 + 31, 32 * ii2 + 30); i0 += 1) {
          for (int i2 = max(32 * ii2, i0 + 1); i2 <= 32 * ii2 + 31; i2 += 1) {
            R[i0][i2] = SCALAR_VAL(0.0);
          }
        }
      } else {
        for (int i0 = 32 * ii0; i0 <= min(_PB_N - 2, 32 * ii0 + 31); i0 += 1) {
          for (int i2 = max(32 * ii2, i0 + 1); i2 < _PB_N; i2 += 1) {
            R[i0][i2] = SCALAR_VAL(0.0);
          }
        }
      }
      if (_PB_N >= 32 * ii2 + 32) {
        for (int ii4 = 0; ii4 <= floord(_PB_M - 1, 32); ii4 += 1) {
          for (int i2 = max(32 * ii0 + 1, 32 * ii2); i2 <= 32 * ii2 + 31; i2 += 1) {
            for (int i4 = 32 * ii4; i4 <= min(_PB_M - 1, 32 * ii4 + 31); i4 += 1) {
              R[32 * ii0][i2] += (Q[i4][32 * ii0] * A[i4][i2]);
            }
          }
        }
        for (int ii4 = 0; ii4 <= floord(_PB_M - 1, 32); ii4 += 1) {
          for (int i2 = max(32 * ii0 + 1, 32 * ii2); i2 <= 32 * ii2 + 31; i2 += 1) {
            for (int i4 = 32 * ii4; i4 <= min(_PB_M - 1, 32 * ii4 + 31); i4 += 1) {
              A[i4][i2] = (A[i4][i2] - (Q[i4][32 * ii0] * R[32 * ii0][i2]));
            }
          }
          if (32 * ii4 + 32 >= _PB_M) {
            for (int i0 = 32 * ii0 + 1; i0 <= min(32 * ii0 + 31, 32 * ii2 + 30); i0 += 1) {
              if (ii2 == ii0 && i0 >= 32 * ii0 + 2) {
                nrm = SCALAR_VAL(0.0);
              }
              if (ii2 == ii0) {
                if (i0 >= 32 * ii0 + 2) {
                  for (int i2 = 0; i2 < _PB_M; i2 += 1) {
                    nrm += (A[i2][i0] * A[i2][i0]);
                  }
                } else {
                  for (int i2 = 32 * ii4; i2 < _PB_M; i2 += 1) {
                    nrm += (A[i2][32 * ii0 + 1] * A[i2][32 * ii0 + 1]);
                  }
                }
                R[i0][i0] = SQRT_FUN(nrm);
                for (int i2 = 0; i2 < _PB_M; i2 += 1) {
                  Q[i2][i0] = (A[i2][i0] / R[i0][i0]);
                }
              }
              for (int i2 = max(32 * ii2, i0 + 1); i2 <= 32 * ii2 + 31; i2 += 1) {
                if (i0 >= 32 * ii0 + 2) {
                  for (int i4 = 0; i4 < _PB_M; i4 += 1) {
                    R[i0][i2] += (Q[i4][i0] * A[i4][i2]);
                  }
                } else {
                  if (ii2 == ii0) {
                    for (int i4 = 0; i4 < 32 * ii4; i4 += 1) {
                      R[32 * ii0 + 1][i2] += (Q[i4][32 * ii0 + 1] * A[i4][i2]);
                    }
                  }
                  for (int i4 = 32 * ii4; i4 < _PB_M; i4 += 1) {
                    R[32 * ii0 + 1][i2] += (Q[i4][32 * ii0 + 1] * A[i4][i2]);
                  }
                }
                for (int i4 = 0; i4 < _PB_M; i4 += 1) {
                  A[i4][i2] = (A[i4][i2] - (Q[i4][i0] * R[i0][i2]));
                }
              }
            }
            if (ii2 == ii0) {
              nrm = SCALAR_VAL(0.0);
              for (int i2 = 0; i2 < _PB_M; i2 += 1) {
                nrm += (A[i2][32 * ii0 + 31] * A[i2][32 * ii0 + 31]);
              }
              R[32 * ii0 + 31][32 * ii0 + 31] = SQRT_FUN(nrm);
              for (int i2 = 0; i2 < _PB_M; i2 += 1) {
                Q[i2][32 * ii0 + 31] = (A[i2][32 * ii0 + 31] / R[32 * ii0 + 31][32 * ii0 + 31]);
              }
            }
          } else if (ii2 >= ii0 + 1) {
            for (int i2 = 32 * ii2; i2 <= 32 * ii2 + 31; i2 += 1) {
              for (int i4 = 32 * ii4; i4 <= 32 * ii4 + 31; i4 += 1) {
                R[32 * ii0 + 1][i2] += (Q[i4][32 * ii0 + 1] * A[i4][i2]);
              }
            }
          } else {
            for (int i2 = 32 * ii4; i2 <= 32 * ii4 + 31; i2 += 1) {
              nrm += (A[i2][32 * ii0 + 1] * A[i2][32 * ii0 + 1]);
            }
          }
        }
      } else {
        for (int ii4 = 0; ii4 <= floord(_PB_M - 1, 32); ii4 += 1) {
          for (int i2 = max(32 * ii0 + 1, 32 * ii2); i2 < _PB_N; i2 += 1) {
            for (int i4 = 32 * ii4; i4 <= min(_PB_M - 1, 32 * ii4 + 31); i4 += 1) {
              R[32 * ii0][i2] += (Q[i4][32 * ii0] * A[i4][i2]);
            }
          }
        }
        if (_PB_M <= 32 && ii2 >= ii0 + 1) {
          for (int i0 = 32 * ii0; i0 <= 32 * ii0 + 31; i0 += 1) {
            for (int i2 = 32 * ii2; i2 < _PB_N; i2 += 1) {
              if (i0 >= 32 * ii0 + 1) {
                for (int i4 = 0; i4 < _PB_M; i4 += 1) {
                  R[i0][i2] += (Q[i4][i0] * A[i4][i2]);
                }
              }
              for (int i4 = 0; i4 < _PB_M; i4 += 1) {
                A[i4][i2] = (A[i4][i2] - (Q[i4][i0] * R[i0][i2]));
              }
            }
          }
        } else if (_PB_M >= 33 && ii2 >= ii0 + 1) {
          for (int ii4 = 0; ii4 <= (_PB_M - 1) / 32; ii4 += 1) {
            if (_PB_M >= 32 * ii4 + 33) {
              for (int i0 = 32 * ii0; i0 <= 32 * ii0 + 1; i0 += 1) {
                if (i0 == 32 * ii0 + 1) {
                  for (int i2 = 32 * ii2; i2 < _PB_N; i2 += 1) {
                    for (int i4 = 32 * ii4; i4 <= 32 * ii4 + 31; i4 += 1) {
                      R[32 * ii0 + 1][i2] += (Q[i4][32 * ii0 + 1] * A[i4][i2]);
                    }
                  }
                } else {
                  for (int i2 = 32 * ii2; i2 < _PB_N; i2 += 1) {
                    for (int i4 = 32 * ii4; i4 <= 32 * ii4 + 31; i4 += 1) {
                      A[i4][i2] = (A[i4][i2] - (Q[i4][32 * ii0] * R[32 * ii0][i2]));
                    }
                  }
                }
              }
            } else {
              for (int i0 = 32 * ii0; i0 <= 32 * ii0 + 31; i0 += 1) {
                if (32 * ii2 >= i0 + 1) {
                  for (int i2 = 32 * ii2; i2 < _PB_N; i2 += 1) {
                    if (i0 >= 32 * ii0 + 1) {
                      if (i0 >= 32 * ii0 + 2) {
                        for (int i4 = 0; i4 < 32 * ii4; i4 += 1) {
                          R[i0][i2] += (Q[i4][i0] * A[i4][i2]);
                        }
                      }
                      for (int i4 = 32 * ii4; i4 < _PB_M; i4 += 1) {
                        R[i0][i2] += (Q[i4][i0] * A[i4][i2]);
                      }
                    }
                    if (i0 >= 32 * ii0 + 1) {
                      for (int i4 = 0; i4 < 32 * ii4; i4 += 1) {
                        A[i4][i2] = (A[i4][i2] - (Q[i4][i0] * R[i0][i2]));
                      }
                    }
                    for (int i4 = 32 * ii4; i4 < _PB_M; i4 += 1) {
                      A[i4][i2] = (A[i4][i2] - (Q[i4][i0] * R[i0][i2]));
                    }
                  }
                }
              }
            }
          }
        } else {
          for (int ii4 = 0; ii4 <= floord(_PB_M - 1, 32); ii4 += 1) {
            if (_PB_M >= 32 * ii4 + 33) {
              for (int i0 = 32 * ii0; i0 <= 32 * ii0 + 1; i0 += 1) {
                if (i0 == 32 * ii0 + 1) {
                  for (int i2 = 32 * ii4; i2 <= 32 * ii4 + 31; i2 += 1) {
                    nrm += (A[i2][32 * ii0 + 1] * A[i2][32 * ii0 + 1]);
                  }
                } else {
                  for (int i2 = 32 * ii0 + 1; i2 < _PB_N; i2 += 1) {
                    for (int i4 = 32 * ii4; i4 <= 32 * ii4 + 31; i4 += 1) {
                      A[i4][i2] = (A[i4][i2] - (Q[i4][32 * ii0] * R[32 * ii0][i2]));
                    }
                  }
                }
              }
            } else {
              for (int i0 = 32 * ii0; i0 < _PB_N; i0 += 1) {
                if (i0 >= 32 * ii0 + 2) {
                  nrm = SCALAR_VAL(0.0);
                }
                if (i0 >= 32 * ii0 + 1) {
                  if (i0 >= 32 * ii0 + 2) {
                    for (int i2 = 0; i2 < 32 * ii4; i2 += 1) {
                      nrm += (A[i2][i0] * A[i2][i0]);
                    }
                  }
                  for (int i2 = 32 * ii4; i2 < _PB_M; i2 += 1) {
                    nrm += (A[i2][i0] * A[i2][i0]);
                  }
                  R[i0][i0] = SQRT_FUN(nrm);
                  for (int i2 = 0; i2 < _PB_M; i2 += 1) {
                    Q[i2][i0] = (A[i2][i0] / R[i0][i0]);
                  }
                }
                for (int i2 = i0 + 1; i2 < _PB_N; i2 += 1) {
                  if (i0 >= 32 * ii0 + 1) {
                    for (int i4 = 0; i4 < _PB_M; i4 += 1) {
                      R[i0][i2] += (Q[i4][i0] * A[i4][i2]);
                    }
                  }
                  if (i0 >= 32 * ii0 + 1) {
                    for (int i4 = 0; i4 < 32 * ii4; i4 += 1) {
                      A[i4][i2] = (A[i4][i2] - (Q[i4][i0] * R[i0][i2]));
                    }
                  }
                  for (int i4 = 32 * ii4; i4 < _PB_M; i4 += 1) {
                    A[i4][i2] = (A[i4][i2] - (Q[i4][i0] * R[i0][i2]));
                  }
                }
              }
            }
          }
        }
      }
    }
  } else if (_PB_M >= _PB_N + 1 && _PB_N >= 32 * ii0 + 2) {
    for (int ii2 = ii0; ii2 <= (_PB_N - 1) / 32; ii2 += 1) {
      for (int i0 = 32 * ii0; i0 <= min(min(_PB_N - 2, 32 * ii0 + 31), 32 * ii2 + 30); i0 += 1) {
        for (int i2 = max(32 * ii2, i0 + 1); i2 <= min(_PB_N - 1, 32 * ii2 + 31); i2 += 1) {
          R[i0][i2] = SCALAR_VAL(0.0);
        }
      }
      if (_PB_M >= 32 * ii2 + 32) {
        for (int ii4 = 0; ii4 <= (_PB_M - 1) / 32; ii4 += 1) {
          for (int i2 = max(32 * ii0 + 1, 32 * ii2); i2 <= min(_PB_N - 1, 32 * ii2 + 31); i2 += 1) {
            for (int i4 = 32 * ii4; i4 <= min(_PB_M - 1, 32 * ii4 + 31); i4 += 1) {
              R[32 * ii0][i2] += (Q[i4][32 * ii0] * A[i4][i2]);
            }
          }
        }
        for (int ii4 = 0; ii4 <= (_PB_M - 1) / 32; ii4 += 1) {
          for (int i2 = max(32 * ii0 + 1, 32 * ii2); i2 <= min(_PB_N - 1, 32 * ii2 + 31); i2 += 1) {
            for (int i4 = 32 * ii4; i4 <= min(_PB_M - 1, 32 * ii4 + 31); i4 += 1) {
              A[i4][i2] = (A[i4][i2] - (Q[i4][32 * ii0] * R[32 * ii0][i2]));
            }
          }
          if (_PB_M >= 32 * ii4 + 33) {
            if (ii2 >= ii0 + 1) {
              for (int i2 = 32 * ii2; i2 <= min(_PB_N - 1, 32 * ii2 + 31); i2 += 1) {
                for (int i4 = 32 * ii4; i4 <= 32 * ii4 + 31; i4 += 1) {
                  R[32 * ii0 + 1][i2] += (Q[i4][32 * ii0 + 1] * A[i4][i2]);
                }
              }
            } else {
              for (int i2 = 32 * ii4; i2 <= 32 * ii4 + 31; i2 += 1) {
                nrm += (A[i2][32 * ii0 + 1] * A[i2][32 * ii0 + 1]);
              }
            }
          } else {
            if (32 * ii0 + 2 == _PB_N && 32 * ii2 + 2 == _PB_N && 32 * ii4 + 31 >= _PB_M) {
              for (int i2 = 32 * ii4; i2 < _PB_M; i2 += 1) {
                nrm += (A[i2][_PB_N - 1] * A[i2][_PB_N - 1]);
              }
              R[_PB_N - 1][_PB_N - 1] = SQRT_FUN(nrm);
              for (int i2 = 0; i2 < _PB_M; i2 += 1) {
                if (32 * ii4 >= i2 + 1) {
                  Q[i2][_PB_N - 1] = (A[i2][_PB_N - 1] / R[_PB_N - 1][_PB_N - 1]);
                } else {
                  Q[i2][_PB_N - 1] = (A[i2][_PB_N - 1] / R[_PB_N - 1][_PB_N - 1]);
                }
              }
            } else if (_PB_N >= 32 * ii0 + 3 && 32 * ii0 + 31 >= _PB_N && ii2 == ii0 && 32 * ii4 + 31 >= _PB_M) {
              for (int i2 = 32 * ii4; i2 < _PB_M; i2 += 1) {
                nrm += (A[i2][32 * ii0 + 1] * A[i2][32 * ii0 + 1]);
              }
              R[32 * ii0 + 1][32 * ii0 + 1] = SQRT_FUN(nrm);
              for (int i2 = 0; i2 < _PB_M; i2 += 1) {
                if (32 * ii4 >= i2 + 1) {
                  Q[i2][32 * ii0 + 1] = (A[i2][32 * ii0 + 1] / R[32 * ii0 + 1][32 * ii0 + 1]);
                } else {
                  Q[i2][32 * ii0 + 1] = (A[i2][32 * ii0 + 1] / R[32 * ii0 + 1][32 * ii0 + 1]);
                }
              }
              for (int i2 = 32 * ii0 + 2; i2 < _PB_N; i2 += 1) {
                for (int i4 = 0; i4 < _PB_M; i4 += 1) {
                  R[32 * ii0 + 1][i2] += (Q[i4][32 * ii0 + 1] * A[i4][i2]);
                }
                for (int i4 = 0; i4 < _PB_M; i4 += 1) {
                  A[i4][i2] = (A[i4][i2] - (Q[i4][32 * ii0 + 1] * R[32 * ii0 + 1][i2]));
                }
              }
            }
            if (32 * ii0 + 31 >= _PB_N && ii2 == ii0 && 32 * ii4 + 31 >= _PB_M) {
              for (int i0 = 32 * ii0 + 2; i0 < _PB_N; i0 += 1) {
                nrm = SCALAR_VAL(0.0);
                for (int i2 = 0; i2 < _PB_M; i2 += 1) {
                  if (32 * ii4 >= i2 + 1) {
                    nrm += (A[i2][i0] * A[i2][i0]);
                  } else {
                    nrm += (A[i2][i0] * A[i2][i0]);
                  }
                }
                R[i0][i0] = SQRT_FUN(nrm);
                for (int i2 = 0; i2 < _PB_M; i2 += 1) {
                  if (32 * ii4 >= i2 + 1) {
                    Q[i2][i0] = (A[i2][i0] / R[i0][i0]);
                  } else {
                    Q[i2][i0] = (A[i2][i0] / R[i0][i0]);
                  }
                }
                for (int i2 = i0 + 1; i2 < _PB_N; i2 += 1) {
                  for (int i4 = 0; i4 < _PB_M; i4 += 1) {
                    R[i0][i2] += (Q[i4][i0] * A[i4][i2]);
                  }
                  for (int i4 = 0; i4 < _PB_M; i4 += 1) {
                    A[i4][i2] = (A[i4][i2] - (Q[i4][i0] * R[i0][i2]));
                  }
                }
              }
            } else if (_PB_N >= 32 * ii0 + 32) {
              for (int i0 = 32 * ii0 + 1; i0 <= min(32 * ii0 + 31, 32 * ii2 + 30); i0 += 1) {
                if (ii2 == ii0 && i0 >= 32 * ii0 + 2) {
                  nrm = SCALAR_VAL(0.0);
                  for (int i2 = 0; i2 < _PB_M; i2 += 1) {
                    if (32 * ii4 + 31 >= _PB_M && i2 >= 32 * ii4) {
                      nrm += (A[i2][i0] * A[i2][i0]);
                    } else if ((i2 % 32) + _PB_M >= i2 + 32) {
                      nrm += (A[i2][i0] * A[i2][i0]);
                    }
                  }
                } else if (ii2 == ii0 && 32 * ii4 + 31 >= _PB_M && i0 == 32 * ii0 + 1) {
                  for (int i2 = 32 * ii4; i2 < _PB_M; i2 += 1) {
                    nrm += (A[i2][32 * ii0 + 1] * A[i2][32 * ii0 + 1]);
                  }
                } else if (ii2 == ii0 && 32 * ii4 + 32 == _PB_M && i0 == 32 * ii0 + 1) {
                  for (int i2 = _PB_M - 32; i2 < _PB_M; i2 += 1) {
                    nrm += (A[i2][32 * ii0 + 1] * A[i2][32 * ii0 + 1]);
                  }
                }
                if (ii2 == ii0) {
                  R[i0][i0] = SQRT_FUN(nrm);
                  for (int i2 = 0; i2 < _PB_M; i2 += 1) {
                    if (32 * ii4 + 31 >= _PB_M && i2 >= 32 * ii4) {
                      Q[i2][i0] = (A[i2][i0] / R[i0][i0]);
                    } else if ((i2 % 32) + _PB_M >= i2 + 32) {
                      Q[i2][i0] = (A[i2][i0] / R[i0][i0]);
                    }
                  }
                }
                for (int i2 = max(32 * ii2, i0 + 1); i2 <= min(_PB_N - 1, 32 * ii2 + 31); i2 += 1) {
                  if (ii2 == ii0 && i0 == 32 * ii0 + 1) {
                    for (int i4 = 0; i4 < 32 * ii4; i4 += 1) {
                      R[32 * ii0 + 1][i2] += (Q[i4][32 * ii0 + 1] * A[i4][i2]);
                    }
                  } else if (i0 >= 32 * ii0 + 2) {
                    for (int i4 = 0; i4 < 32 * ii4; i4 += 1) {
                      R[i0][i2] += (Q[i4][i0] * A[i4][i2]);
                    }
                  }
                  for (int i4 = 32 * ii4; i4 < _PB_M; i4 += 1) {
                    R[i0][i2] += (Q[i4][i0] * A[i4][i2]);
                  }
                  for (int i4 = 0; i4 < _PB_M; i4 += 1) {
                    A[i4][i2] = (A[i4][i2] - (Q[i4][i0] * R[i0][i2]));
                  }
                }
              }
              if (ii2 == ii0 && 32 * ii4 + 31 >= _PB_M) {
                nrm = SCALAR_VAL(0.0);
                for (int i2 = 0; i2 < _PB_M; i2 += 1) {
                  if (i2 >= 32 * ii4) {
                    nrm += (A[i2][32 * ii0 + 31] * A[i2][32 * ii0 + 31]);
                  } else {
                    nrm += (A[i2][32 * ii0 + 31] * A[i2][32 * ii0 + 31]);
                  }
                }
                R[32 * ii0 + 31][32 * ii0 + 31] = SQRT_FUN(nrm);
                for (int i2 = 0; i2 < _PB_M; i2 += 1) {
                  if (i2 >= 32 * ii4) {
                    Q[i2][32 * ii0 + 31] = (A[i2][32 * ii0 + 31] / R[32 * ii0 + 31][32 * ii0 + 31]);
                  } else {
                    Q[i2][32 * ii0 + 31] = (A[i2][32 * ii0 + 31] / R[32 * ii0 + 31][32 * ii0 + 31]);
                  }
                }
              }
              if (ii2 == ii0 && 32 * ii4 + 32 == _PB_M) {
                nrm = SCALAR_VAL(0.0);
                for (int i2 = 0; i2 < _PB_M; i2 += 1) {
                  nrm += (A[i2][32 * ii0 + 31] * A[i2][32 * ii0 + 31]);
                }
                R[32 * ii0 + 31][32 * ii0 + 31] = SQRT_FUN(nrm);
                for (int i2 = 0; i2 < _PB_M; i2 += 1) {
                  Q[i2][32 * ii0 + 31] = (A[i2][32 * ii0 + 31] / R[32 * ii0 + 31][32 * ii0 + 31]);
                }
              }
            } else {
              for (int i0 = 32 * ii0 + 1; i0 < _PB_N; i0 += 1) {
                if (i0 >= 32 * ii0 + 2) {
                  nrm = SCALAR_VAL(0.0);
                }
                if (i0 >= 32 * ii0 + 2) {
                  for (int i2 = 0; i2 < _PB_M; i2 += 1) {
                    nrm += (A[i2][i0] * A[i2][i0]);
                  }
                } else {
                  for (int i2 = _PB_M - 32; i2 < _PB_M; i2 += 1) {
                    nrm += (A[i2][32 * ii0 + 1] * A[i2][32 * ii0 + 1]);
                  }
                }
                R[i0][i0] = SQRT_FUN(nrm);
                for (int i2 = 0; i2 < _PB_M; i2 += 1) {
                  Q[i2][i0] = (A[i2][i0] / R[i0][i0]);
                }
                for (int i2 = i0 + 1; i2 < _PB_N; i2 += 1) {
                  for (int i4 = 0; i4 < _PB_M; i4 += 1) {
                    R[i0][i2] += (Q[i4][i0] * A[i4][i2]);
                  }
                  for (int i4 = 0; i4 < _PB_M; i4 += 1) {
                    A[i4][i2] = (A[i4][i2] - (Q[i4][i0] * R[i0][i2]));
                  }
                }
              }
            }
          }
        }
      } else if (ii2 >= ii0 + 1) {
        for (int ii3 = 1; ii3 <= 2; ii3 += 1) {
          for (int ii4 = 0; ii4 <= ii2 - ii3 + 1; ii4 += 1) {
            if (ii3 == 2) {
              for (int i2 = 32 * ii2; i2 < _PB_N; i2 += 1) {
                for (int i4 = 32 * ii4; i4 <= 32 * ii4 + 31; i4 += 1) {
                  A[i4][i2] = (A[i4][i2] - (Q[i4][32 * ii0] * R[32 * ii0][i2]));
                }
              }
              for (int i2 = 32 * ii2; i2 < _PB_N; i2 += 1) {
                for (int i4 = 32 * ii4; i4 <= 32 * ii4 + 31; i4 += 1) {
                  R[32 * ii0 + 1][i2] += (Q[i4][32 * ii0 + 1] * A[i4][i2]);
                }
              }
            } else {
              for (int i2 = 32 * ii2; i2 < _PB_N; i2 += 1) {
                for (int i4 = 32 * ii4; i4 <= min(_PB_M - 1, 32 * ii4 + 31); i4 += 1) {
                  R[32 * ii0][i2] += (Q[i4][32 * ii0] * A[i4][i2]);
                }
              }
            }
          }
          if (ii3 == 2) {
            for (int i0 = 32 * ii0; i0 <= 32 * ii0 + 31; i0 += 1) {
              if (32 * ii2 >= i0 + 1) {
                for (int i2 = 32 * ii2; i2 < _PB_N; i2 += 1) {
                  if (i0 >= 32 * ii0 + 1) {
                    if (i0 >= 32 * ii0 + 2) {
                      for (int i4 = 0; i4 < 32 * ii2; i4 += 1) {
                        R[i0][i2] += (Q[i4][i0] * A[i4][i2]);
                      }
                    }
                    for (int i4 = 32 * ii2; i4 < _PB_M; i4 += 1) {
                      R[i0][i2] += (Q[i4][i0] * A[i4][i2]);
                    }
                  }
                  if (i0 >= 32 * ii0 + 1) {
                    for (int i4 = 0; i4 < 32 * ii2; i4 += 1) {
                      A[i4][i2] = (A[i4][i2] - (Q[i4][i0] * R[i0][i2]));
                    }
                  }
                  for (int i4 = 32 * ii2; i4 < _PB_M; i4 += 1) {
                    A[i4][i2] = (A[i4][i2] - (Q[i4][i0] * R[i0][i2]));
                  }
                }
              }
            }
          }
        }
      } else {
        for (int ii3 = 1; ii3 <= 2; ii3 += 1) {
          for (int ii4 = 0; ii4 <= ii0 - ii3 + 1; ii4 += 1) {
            if (ii3 == 2) {
              for (int i2 = 32 * ii0 + 1; i2 < _PB_N; i2 += 1) {
                for (int i4 = 32 * ii4; i4 <= 32 * ii4 + 31; i4 += 1) {
                  A[i4][i2] = (A[i4][i2] - (Q[i4][32 * ii0] * R[32 * ii0][i2]));
                }
              }
              for (int i2 = 32 * ii4; i2 <= 32 * ii4 + 31; i2 += 1) {
                nrm += (A[i2][32 * ii0 + 1] * A[i2][32 * ii0 + 1]);
              }
            } else {
              for (int i2 = 32 * ii0 + 1; i2 < _PB_N; i2 += 1) {
                for (int i4 = 32 * ii4; i4 <= min(_PB_M - 1, 32 * ii4 + 31); i4 += 1) {
                  R[32 * ii0][i2] += (Q[i4][32 * ii0] * A[i4][i2]);
                }
              }
            }
          }
          if (ii3 == 2) {
            for (int i0 = 32 * ii0; i0 < _PB_N - 1; i0 += 1) {
              if (i0 >= 32 * ii0 + 2) {
                nrm = SCALAR_VAL(0.0);
              }
              if (i0 >= 32 * ii0 + 1) {
                if (i0 >= 32 * ii0 + 2) {
                  for (int i2 = 0; i2 < _PB_M; i2 += 1) {
                    nrm += (A[i2][i0] * A[i2][i0]);
                  }
                } else {
                  for (int i2 = 32 * ii0; i2 < _PB_M; i2 += 1) {
                    nrm += (A[i2][32 * ii0 + 1] * A[i2][32 * ii0 + 1]);
                  }
                }
                R[i0][i0] = SQRT_FUN(nrm);
                for (int i2 = 0; i2 < _PB_M; i2 += 1) {
                  Q[i2][i0] = (A[i2][i0] / R[i0][i0]);
                }
              }
              for (int i2 = i0 + 1; i2 < _PB_N; i2 += 1) {
                if (i0 >= 32 * ii0 + 1) {
                  for (int i4 = 0; i4 < _PB_M; i4 += 1) {
                    R[i0][i2] += (Q[i4][i0] * A[i4][i2]);
                  }
                }
                if (i0 >= 32 * ii0 + 1) {
                  for (int i4 = 0; i4 < 32 * ii0; i4 += 1) {
                    A[i4][i2] = (A[i4][i2] - (Q[i4][i0] * R[i0][i2]));
                  }
                }
                for (int i4 = 32 * ii0; i4 < _PB_M; i4 += 1) {
                  A[i4][i2] = (A[i4][i2] - (Q[i4][i0] * R[i0][i2]));
                }
              }
            }
            if (_PB_N >= 32 * ii0 + 3) {
              nrm = SCALAR_VAL(0.0);
            }
            if (_PB_N >= 32 * ii0 + 3) {
              for (int i2 = 0; i2 < _PB_M; i2 += 1) {
                nrm += (A[i2][_PB_N - 1] * A[i2][_PB_N - 1]);
              }
            } else {
              for (int i2 = _PB_N - 2; i2 < _PB_M; i2 += 1) {
                nrm += (A[i2][_PB_N - 1] * A[i2][_PB_N - 1]);
              }
            }
            R[_PB_N - 1][_PB_N - 1] = SQRT_FUN(nrm);
            for (int i2 = 0; i2 < _PB_M; i2 += 1) {
              Q[i2][_PB_N - 1] = (A[i2][_PB_N - 1] / R[_PB_N - 1][_PB_N - 1]);
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
  int m = M;
  int n = N;

  /* Variable declaration/allocation. */
  POLYBENCH_2D_ARRAY_DECL(A,DATA_TYPE,M,N,m,n);
  POLYBENCH_2D_ARRAY_DECL(R,DATA_TYPE,N,N,n,n);
  POLYBENCH_2D_ARRAY_DECL(Q,DATA_TYPE,M,N,m,n);

  /* Initialize array(s). */
  init_array (m, n,
	      POLYBENCH_ARRAY(A),
	      POLYBENCH_ARRAY(R),
	      POLYBENCH_ARRAY(Q));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_gramschmidt (m, n,
		      POLYBENCH_ARRAY(A),
		      POLYBENCH_ARRAY(R),
		      POLYBENCH_ARRAY(Q));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(m, n, POLYBENCH_ARRAY(A), POLYBENCH_ARRAY(R), POLYBENCH_ARRAY(Q)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(R);
  POLYBENCH_FREE_ARRAY(Q);

  return 0;
}
