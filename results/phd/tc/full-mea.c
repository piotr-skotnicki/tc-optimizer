/**
 * This version is stamped on Oct 25, 2025
 */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#define POLYBENCH_PADDING_FACTOR 2

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "full-mea.h"

static
int paired(char a, char b)
{
  return ((a == 'A' && b == 'U') || (a == 'U' && b == 'A') ||
          (a == 'G' && b == 'C') || (a == 'C' && b == 'G') ||
          (a == 'G' && b == 'U') || (a == 'U' && b == 'G'));
}

static
DATA_TYPE max_score(DATA_TYPE a, DATA_TYPE b)
{
  return a > b ? a : b;
}

/* Array initialization. */
static
void init_array (int n,
        int l,
        DATA_TYPE ERT,
        char POLYBENCH_1D(RNA,N,n),
        DATA_TYPE POLYBENCH_2D(Q,N,N,n,n),
        DATA_TYPE POLYBENCH_2D(Qbp,N,N,n,n),
        DATA_TYPE POLYBENCH_2D(Pbp,N,N,n,n),
        DATA_TYPE POLYBENCH_2D(Pu,N,N,n,n),
        DATA_TYPE POLYBENCH_2D(M,N,N,n,n))
{
  int i, j, k, p, q;

  srand(42);

  for (i = 0; i < _PB_N; i++) {
    RNA[i] = "ACGU"[rand() % 4];
  }

  //memcpy(RNA, "GGUCCAC", 7);

  for (i = 0; i < _PB_N + POLYBENCH_PADDING_FACTOR; i++) {
    for (j = 0; j < _PB_N + POLYBENCH_PADDING_FACTOR; j++) {
      Q[i][j] = SCALAR_VAL(1.0);
      Qbp[i][j] = SCALAR_VAL(0.0);
      Pbp[i][j] = SCALAR_VAL(0.0);
      Pu[i][j] = SCALAR_VAL(0.0);
      M[i][j] = SCALAR_VAL(0.0);
    }
  }
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
        DATA_TYPE POLYBENCH_2D(M,N,N,n,n))
{
  int i, j;
  int t = 0;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("M");
  for (i = 1; i <= _PB_N; i++) {
    for (j = 0; j <= _PB_N; j++) {
      if (t % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
      fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, M[i][j]);
      t++;
    }
  }
  POLYBENCH_DUMP_END("M");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_full_mea(int n,
        int l,
        DATA_TYPE ERT,
        DATA_TYPE gamma,
        char POLYBENCH_1D(RNA,N,n),
        DATA_TYPE POLYBENCH_2D(Q,N,N,n,n),
        DATA_TYPE POLYBENCH_2D(Qbp,N,N,n,n),
        DATA_TYPE POLYBENCH_2D(Pbp,N,N,n,n),
        DATA_TYPE POLYBENCH_2D(Pu,N,N,n,n),
        DATA_TYPE POLYBENCH_2D(M,N,N,n,n))
{
  int i, j, k, p, q;

/* TC Optimizing Compiler 0.5.1 */
/* ./tc ../examples/rnatools/mcc.scop.c --correction-tiling --isl-wave-scheduling --omp-for-codegen --floyd-warshall-tc --debug --align -b 16 */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#pragma scop
if (l >= 0 && l <= 5) {
  for (register int k0 = -1; k0 < _PB_N / 16; k0 += 1) {
    #pragma omp parallel for
    for (register int k1 = max(k0 - (_PB_N + 16) / 16 + 1, -((_PB_N + 14) / 16)); k1 < 0; k1 += 1) {
      {
        for (register int i0 = max(max(-_PB_N + 1, -16 * k0 + 16 * k1 - 14), 16 * k1); i0 <= 16 * k1 + 15; i0 += 1) {
          for (register int i1 = max(16 * k0 - 16 * k1, -i0 + 1); i1 <= min(min(_PB_N, 16 * k0 - 16 * k1 + 15), l - i0 + 1); i1 += 1) {
            Q[-i0][i1] = Q[-i0][i1 - 1];
          }
          if (16 * k0 + i0 >= l + 16 * k1 + 2) {
            Q[-i0][16 * k0 - 16 * k1] = Q[-i0][16 * k0 - 16 * k1 - 1];
          }
        }
        if (k0 >= 1 && k1 == -1) {
          for (register int k4 = 1; k4 <= 2; k4 += 1) {
            for (register int i0 = -15; i0 < 0; i0 += 1) {
              if (k4 == 2) {
                Q[-i0][16 * k0 + 16] += (Q[-i0][-i0 - 1] * Qbp[-i0][16 * k0 + 16]);
              }
              if (i0 + 16 >= k4) {
                if (k4 == 2) {
                  for (register int i3 = -i0 + 1; i3 <= 15; i3 += 1) {
                    Qbp[-i0][16 * k0 + 16] = ((paired(RNA[-i0 - 1], RNA[16 * k0 + 15]) * Q[-i0 + 1][16 * k0 + 15]) * ERT);
                  }
                } else {
                  Qbp[-i0][16 * k0 + 16] = ((paired(RNA[-i0 - 1], RNA[16 * k0 + 15]) * Q[-i0 + 1][16 * k0 + 15]) * ERT);
                }
              }
            }
          }
        } else if (_PB_N >= l + 2 && k0 <= 0 && k1 == -1) {
          for (register int k4 = 1; k4 <= 2; k4 += 1) {
            for (register int i0 = max(max(-15, -_PB_N + l + 1), l - 16 * k0 - 30); i0 < min(0, l - 16 * k0 + 28 * k4 - 41); i0 += 1) {
              for (register int i1 = max(16 * k0 + 16, l - i0 + 1); i1 <= min(min(min(_PB_N, l + 17), 16 * k0 + 31), l + 28 * k4 - i0 - 26); i1 += 1) {
                if (k4 == 2 && i1 >= 16 * k0 + 17 && i0 + i1 >= l + 2) {
                  Q[-i0][i1] = Q[-i0][i1 - 1];
                }
                if (l + i0 + 33 >= k4 + i1) {
                  if (k0 == 0 && k4 == 2 && i1 == l + 17) {
                    Qbp[-i0][l + 17] = ((paired(RNA[-i0 - 1], RNA[l + 16]) * Q[-i0 + 1][l + 16]) * ERT);
                    Q[-i0][l + 17] += (Q[-i0][-i0 - 1] * Qbp[-i0][l + 17]);
                  }
                  for (register int i3 = max(-i0, -l + k4 - i0 + i1 - 18); i3 <= min(min(15, 14 * k4 - i0 - 14), -l + i1 - 1); i3 += 1) {
                    if (k4 == 2 && i1 >= 16 * k0 + 17 && 15 * i0 + 4 * i1 + 11 * i3 >= 4 * l + 12) {
                      Qbp[-i0][i1] = ((paired(RNA[-i0 - 1], RNA[i1 - 1]) * Q[-i0 + 1][i1 - 1]) * ERT);
                    } else if (k4 == 1 && i0 + i3 == 0) {
                      Qbp[-i0][i1] = ((paired(RNA[-i0 - 1], RNA[i1 - 1]) * Q[-i0 + 1][i1 - 1]) * ERT);
                    } else if (k0 == 0 && i1 == 16 && i0 + i3 >= 1) {
                      Qbp[-i0][16] = ((paired(RNA[-i0 - 1], RNA[15]) * Q[-i0 + 1][15]) * ERT);
                    }
                    if (k4 == 2 && l + 16 >= i1) {
                      Q[-i0][i1] += (Q[-i0][i3 - 1] * Qbp[i3][i1]);
                    }
                  }
                } else {
                  Q[15][l + 17] += (Q[15][14] * Qbp[15][l + 17]);
                }
              }
            }
            if (k0 == 0 && k4 == 1) {
              for (register int i0 = l - 13; i0 < 0; i0 += 1) {
                Qbp[-i0][16] = ((paired(RNA[-i0 - 1], RNA[15]) * Q[-i0 + 1][15]) * ERT);
              }
            }
          }
        }
      }
      for (register int k3 = max(1, -k1 - 1); k3 < k0 - k1 - 1; k3 += 1) {
        for (register int k4 = 1; k4 <= 2; k4 += 1) {
          {
            if (k4 == 1) {
              for (register int i0 = 15 * k1 - k3; i0 <= k1 - 15 * k3; i0 += 1) {
                for (register int i1 = 16 * k0 - 16 * k1; i1 <= min(_PB_N, 16 * k0 - k1 + 15 * k3 + 15); i1 += 1) {
                  Qbp[-i0][i1] = ((paired(RNA[-i0 - 1], RNA[i1 - 1]) * Q[-i0 + 1][i1 - 1]) * ERT);
                }
              }
            }
            for (register int i0 = max(16 * k1, -16 * k3 - k4 + 2); i0 <= 16 * k1 - 15 * k4 + 30; i0 += 1) {
              for (register int i3 = 16 * k3; i3 <= 16 * k3 + 15; i3 += 1) {
                if (k4 == 1 && i0 >= 16 * k1 + 2) {
                  Qbp[-i0][16 * k0 - 16 * k1] = ((paired(RNA[-i0 - 1], RNA[16 * k0 - 16 * k1 - 1]) * Q[-i0 + 1][16 * k0 - 16 * k1 - 1]) * ERT);
                } else if (k4 == 1 && i0 == 16 * k1 + 1) {
                  Qbp[-16 * k1 - 1][16 * k0 - 16 * k1] = ((paired(RNA[-16 * k1 - 2], RNA[16 * k0 - 16 * k1 - 1]) * Q[-16 * k1][16 * k0 - 16 * k1 - 1]) * ERT);
                } else if (k4 == 1) {
                  Qbp[-16 * k1][16 * k0 - 16 * k1] = ((paired(RNA[-16 * k1 - 1], RNA[16 * k0 - 16 * k1 - 1]) * Q[-16 * k1 + 1][16 * k0 - 16 * k1 - 1]) * ERT);
                } else {
                  if (k1 + k3 == 0 && 16 * k1 + i3 >= 1) {
                    Qbp[-16 * k1][16 * k0 - 16 * k1] = ((paired(RNA[-16 * k1 - 1], RNA[16 * k0 - 16 * k1 - 1]) * Q[-16 * k1 + 1][16 * k0 - 16 * k1 - 1]) * ERT);
                  }
                  Q[-16 * k1][16 * k0 - 16 * k1] += (Q[-16 * k1][i3 - 1] * Qbp[i3][16 * k0 - 16 * k1]);
                }
              }
            }
          }
          if (k1 + k3 == -1 && k4 == 2) {
            for (register int i0 = 16 * k1 + 1; i0 <= 16 * k1 + 15; i0 += 1) {
              Q[-i0][16 * k0 - 16 * k1] += (Q[-i0][-i0 - 1] * Qbp[-i0][16 * k0 - 16 * k1]);
              for (register int i3 = -i0 + 1; i3 < -16 * k1; i3 += 1) {
                Qbp[-i0][16 * k0 - 16 * k1] = ((paired(RNA[-i0 - 1], RNA[16 * k0 - 16 * k1 - 1]) * Q[-i0 + 1][16 * k0 - 16 * k1 - 1]) * ERT);
              }
            }
          }
        }
      }
      if (_PB_N + 16 * k1 + 14 >= l) {
        for (register int k3 = max(max(1, -k1 - 1), k0 - k1 - 1); k3 <= min(k0 - k1, (_PB_N - l - 1) / 16); k3 += 1) {
          for (register int k4 = 1; k4 <= 2; k4 += 1) {
            if (k0 + 1 >= k1 + k3 + k4 && _PB_N + 15 >= l + 16 * k3 + 16 * k4) {
              if (k1 + k3 == 0 && k4 == 1) {
                for (register int i1 = max(16 * k0 - 16 * k1, l - 16 * k1 + 1); i1 <= min(_PB_N, 16 * k0 - 16 * k1 + 15); i1 += 1) {
                  Qbp[-16 * k1][i1] = ((paired(RNA[-16 * k1 - 1], RNA[i1 - 1]) * Q[-16 * k1 + 1][i1 - 1]) * ERT);
                }
              } else if (k1 + k3 == -1 && k4 == 1) {
                for (register int i0 = max(max(-_PB_N + l + 1, l - 16 * k0 + 16 * k1 - 14), 16 * k1 + 1); i0 <= l + 16 * k1 + 2; i0 += 1) {
                  for (register int i1 = max(16 * k0 - 16 * k1, l - i0 + 1); i1 <= min(min(_PB_N, 16 * k0 - 16 * k1 + 15), l - i0 + 2); i1 += 1) {
                    Qbp[-i0][i1] = ((paired(RNA[-i0 - 1], RNA[i1 - 1]) * Q[-i0 + 1][i1 - 1]) * ERT);
                  }
                }
              } else if (k0 >= 1 && k1 + k3 + 1 == k0 && k4 == 2) {
                for (register int i1 = 16 * k0 - 16 * k1; i1 <= l + 16 * k0 - 16 * k1 + 1; i1 += 1) {
                  if (16 * k1 + i1 >= 16 * k0 + 1) {
                    Q[-16 * k1][i1] = Q[-16 * k1][i1 - 1];
                  }
                  {
                    if (16 * k1 + i1 >= 16 * k0 + 1) {
                      for (register int i3 = -16 * k1; i3 < 16 * k0 - 16 * k1 - 16; i3 += 1) {
                        if (16 * k1 + i3 >= 1) {
                          Qbp[-16 * k1][i1] = ((paired(RNA[-16 * k1 - 1], RNA[i1 - 1]) * Q[-16 * k1 + 1][i1 - 1]) * ERT);
                        }
                        Q[-16 * k1][i1] += (Q[-16 * k1][i3 - 1] * Qbp[i3][i1]);
                      }
                    }
                    for (register int i3 = 16 * k0 - 16 * k1 - 16; i3 < min(16 * k0 - 16 * k1, -l + i1); i3 += 1) {
                      if (16 * k1 + i1 >= 16 * k0 + 1 && 16 * k1 + i3 >= 1) {
                        Qbp[-16 * k1][i1] = ((paired(RNA[-16 * k1 - 1], RNA[i1 - 1]) * Q[-16 * k1 + 1][i1 - 1]) * ERT);
                      } else if (k0 == 1 && 16 * k1 + i1 == 16 && 16 * k1 + i3 >= 1) {
                        Qbp[-16 * k1][-16 * k1 + 16] = ((paired(RNA[-16 * k1 - 1], RNA[-16 * k1 + 15]) * Q[-16 * k1 + 1][-16 * k1 + 15]) * ERT);
                      }
                      Q[-16 * k1][i1] += (Q[-16 * k1][i3 - 1] * Qbp[i3][i1]);
                    }
                  }
                }
              }
              if (k0 >= 0) {
                for (register int i0 = max(max(16 * k1 + k4 - 1, 15 * k1 - k3 - k4 + 2), l + 9 * k1 - 7 * k3 - 7 * k4 + 3); i0 <= 16 * k1 + 15; i0 += 1) {
                  if (k0 >= 1 && k4 == 1 && 16 * k3 + i0 >= 2 && i0 >= 16 * k1 + 1) {
                    for (register int i1 = l <= 4 && k1 + k3 + 1 == k0 ? 16 * k0 - 16 * k1 : l + 10 * k0 - 10 * k1 + 6 * k3 + 1; i1 <= (l <= 4 && k1 + k3 + 1 == k0 ? 16 * k0 - 16 * k1 : l + 10 * k0 - 10 * k1 + 6 * k3 + 1); i1 += 1) {
                      for (register int i3 = -i0 + 1; i3 < 16 * k3; i3 += 1) {
                        Q[-i0][i1] += (Q[-i0][i3 - 1] * Qbp[i3][i1]);
                      }
                      if (k1 + k3 == k0 && 16 * k1 + i1 == l + 16 * k0 + 1) {
                        Qbp[-i0][l + 16 * k0 - 16 * k1 + 1] = ((paired(RNA[-i0 - 1], RNA[l + 16 * k0 - 16 * k1]) * Q[-i0 + 1][l + 16 * k0 - 16 * k1]) * ERT);
                      } else {
                        for (register int i3 = 16 * k0 - 16 * k1 - 16; i3 < -l + 16 * k0 - 16 * k1; i3 += 1) {
                          Qbp[-i0][16 * k0 - 16 * k1] = ((paired(RNA[-i0 - 1], RNA[16 * k0 - 16 * k1 - 1]) * Q[-i0 + 1][16 * k0 - 16 * k1 - 1]) * ERT);
                        }
                      }
                    }
                  } else if (k1 + k3 + 1 == k0 && k4 == 2) {
                    for (register int i1 = max(16 * k0 - 16 * k1, l - i0 + 1); i1 <= l + 16 * k0 - 16 * k1; i1 += 1) {
                      if (i0 + i1 >= l + 2 && 16 * k1 + i1 >= 16 * k0 + 1) {
                        Q[-i0][i1] = Q[-i0][i1 - 1];
                      }
                      if (16 * k1 + i1 >= 16 * k0 + 1) {
                        if (k0 == 0 && l + 2 >= i0 + i1) {
                          Q[-i0][i1] += (Q[-i0][-i0 - 1] * Qbp[-i0][i1]);
                        }
                        for (register int i3 = max(-i0, -i0 - (-4 * l + 4 * i0 + 4 * i1 + 10) / 11 + 2); i3 < -l + i1; i3 += 1) {
                          Qbp[-i0][i1] = ((paired(RNA[-i0 - 1], RNA[i1 - 1]) * Q[-i0 + 1][i1 - 1]) * ERT);
                          Q[-i0][i1] += (Q[-i0][i3 - 1] * Qbp[i3][i1]);
                        }
                      } else {
                        for (register int i3 = max(16 * k0 - 16 * k1 - 16, -i0); i3 < -l + 16 * k0 - 16 * k1; i3 += 1) {
                          if (k0 == 0 && i0 + i3 >= 1) {
                            Qbp[-i0][-16 * k1] = ((paired(RNA[-i0 - 1], RNA[-16 * k1 - 1]) * Q[-i0 + 1][-16 * k1 - 1]) * ERT);
                          }
                          Q[-i0][16 * k0 - 16 * k1] += (Q[-i0][i3 - 1] * Qbp[i3][16 * k0 - 16 * k1]);
                        }
                      }
                    }
                    Q[-i0][l + 16 * k0 - 16 * k1 + 1] = Q[-i0][l + 16 * k0 - 16 * k1];
                    if (16 * k0 + i0 >= 16 * k1 + 2) {
                      Qbp[-i0][l + 16 * k0 - 16 * k1 + 1] = ((paired(RNA[-i0 - 1], RNA[l + 16 * k0 - 16 * k1]) * Q[-i0 + 1][l + 16 * k0 - 16 * k1]) * ERT);
                    }
                    Q[-i0][l + 16 * k0 - 16 * k1 + 1] += (Q[-i0][-i0 - 1] * Qbp[-i0][l + 16 * k0 - 16 * k1 + 1]);
                    for (register int i3 = -i0 + 1; i3 < 16 * k0 - 16 * k1; i3 += 1) {
                      Qbp[-i0][l + 16 * k0 - 16 * k1 + 1] = ((paired(RNA[-i0 - 1], RNA[l + 16 * k0 - 16 * k1]) * Q[-i0 + 1][l + 16 * k0 - 16 * k1]) * ERT);
                    }
                  } else if (k0 == 0 && k1 + k3 == 0) {
                    for (register int i3 = -i0 + 1; i3 < -16 * k1; i3 += 1) {
                      Q[-i0][l - 16 * k1 + 1] += (Q[-i0][i3 - 1] * Qbp[i3][l - 16 * k1 + 1]);
                    }
                    Qbp[-i0][l - 16 * k1 + 1] = ((paired(RNA[-i0 - 1], RNA[l - 16 * k1]) * Q[-i0 + 1][l - 16 * k1]) * ERT);
                  } else if (k0 == 0 && k1 + k3 == -1) {
                    Qbp[-i0][-16 * k1] = ((paired(RNA[-i0 - 1], RNA[-16 * k1 - 1]) * Q[-i0 + 1][-16 * k1 - 1]) * ERT);
                  } else if (k1 + k3 == k0 && i0 == 16 * k1) {
                    Qbp[-16 * k1][l + 16 * k0 - 16 * k1 + 1] = ((paired(RNA[-16 * k1 - 1], RNA[l + 16 * k0 - 16 * k1]) * Q[-16 * k1 + 1][l + 16 * k0 - 16 * k1]) * ERT);
                  } else if (i0 == 16 * k1) {
                    for (register int i3 = 16 * k0 - 16 * k1 - 16; i3 < -l + 16 * k0 - 16 * k1; i3 += 1) {
                      Qbp[-16 * k1][16 * k0 - 16 * k1] = ((paired(RNA[-16 * k1 - 1], RNA[16 * k0 - 16 * k1 - 1]) * Q[-16 * k1 + 1][16 * k0 - 16 * k1 - 1]) * ERT);
                    }
                  } else {
                    for (register int i3 = -16 * k1; i3 <= -l - 16 * k1 + 15; i3 += 1) {
                      Qbp[-16 * k1 - 1][-16 * k1 + 16] = ((paired(RNA[-16 * k1 - 2], RNA[-16 * k1 + 15]) * Q[-16 * k1][-16 * k1 + 15]) * ERT);
                    }
                  }
                }
              } else {
                for (register int i0 = max(-_PB_N + l + 1, l + 16 * k1 + 3); i0 <= 16 * k1 + 15; i0 += 1) {
                  for (register int i1 = l - i0 + 1; i1 <= min(_PB_N, l - i0 + 2); i1 += 1) {
                    Qbp[-i0][i1] = ((paired(RNA[-i0 - 1], RNA[i1 - 1]) * Q[-i0 + 1][i1 - 1]) * ERT);
                  }
                }
              }
            } else if (k1 + k3 == k0) {
              for (register int i0 = max(max(-_PB_N + l + 1, l - 16 * k0 + 16 * k1 - 14), 16 * k1); i0 <= 16 * k1 + 15; i0 += 1) {
                for (register int i1 = max(l + 16 * k0 - 16 * k1 + 1, l - i0 + 1); i1 <= min(_PB_N, 16 * k0 - 16 * k1 + 15); i1 += 1) {
                  if (16 * k1 + i1 >= l + 16 * k0 + 2 && i0 + i1 >= l + 2) {
                    Q[-i0][i1] = Q[-i0][i1 - 1];
                  }
                  {
                    if (16 * k1 + i1 >= l + 16 * k0 + 2) {
                      for (register int i3 = -i0; i3 < 16 * k0 - 16 * k1; i3 += 1) {
                        if (2 * i0 + i3 >= 16 * k1 + 1) {
                          Qbp[-i0][i1] = ((paired(RNA[-i0 - 1], RNA[i1 - 1]) * Q[-i0 + 1][i1 - 1]) * ERT);
                        }
                        Q[-i0][i1] += (Q[-i0][i3 - 1] * Qbp[i3][i1]);
                      }
                    }
                    for (register int i3 = max(16 * k0 - 16 * k1, -i0); i3 < -l + i1; i3 += 1) {
                      if (16 * k1 + i1 >= l + 16 * k0 + 2 && 2 * i0 + i3 >= 16 * k1 + 1 && 15 * i0 + 4 * i1 + 11 * i3 >= 4 * l + 12) {
                        Qbp[-i0][i1] = ((paired(RNA[-i0 - 1], RNA[i1 - 1]) * Q[-i0 + 1][i1 - 1]) * ERT);
                      }
                      Q[-i0][i1] += (Q[-i0][i3 - 1] * Qbp[i3][i1]);
                    }
                  }
                }
              }
            } else {
              for (register int i0 = max(-_PB_N + l + 1, 16 * k1); i0 <= 16 * k1 + 15; i0 += 1) {
                for (register int i1 = max(16 * k0 - 16 * k1, l - i0 + 1); i1 <= _PB_N; i1 += 1) {
                  if (i0 + i1 >= l + 2 && 16 * k1 + i1 >= 16 * k0 + 1) {
                    Q[-i0][i1] = Q[-i0][i1 - 1];
                  }
                  {
                    if (16 * k1 + i1 >= 16 * k0 + 1) {
                      for (register int i3 = -i0; i3 < 16 * k0 - 16 * k1 - 16; i3 += 1) {
                        if (2 * i0 + i3 >= 16 * k1 + 1) {
                          Qbp[-i0][i1] = ((paired(RNA[-i0 - 1], RNA[i1 - 1]) * Q[-i0 + 1][i1 - 1]) * ERT);
                        }
                        Q[-i0][i1] += (Q[-i0][i3 - 1] * Qbp[i3][i1]);
                      }
                    }
                    for (register int i3 = max(16 * k0 - 16 * k1 - 16, -i0); i3 < -l + i1; i3 += 1) {
                      if (16 * k1 + i1 >= 16 * k0 + 1 && 2 * i0 + i3 >= 16 * k1 + 1 && 15 * i0 + 4 * i1 + 11 * i3 >= 4 * l + 12) {
                        Qbp[-i0][i1] = ((paired(RNA[-i0 - 1], RNA[i1 - 1]) * Q[-i0 + 1][i1 - 1]) * ERT);
                      } else if (k0 == 0 && 16 * k1 + i1 == 0 && i0 + i3 >= 1) {
                        Qbp[-i0][-16 * k1] = ((paired(RNA[-i0 - 1], RNA[-16 * k1 - 1]) * Q[-i0 + 1][-16 * k1 - 1]) * ERT);
                      } else if (k0 == 1 && i0 == 16 * k1 && 16 * k1 + i1 == 16 && 16 * k1 + i3 >= 1) {
                        Qbp[-16 * k1][-16 * k1 + 16] = ((paired(RNA[-16 * k1 - 1], RNA[-16 * k1 + 15]) * Q[-16 * k1 + 1][-16 * k1 + 15]) * ERT);
                      }
                      Q[-i0][i1] += (Q[-i0][i3 - 1] * Qbp[i3][i1]);
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
#pragma endscop

/* TC Optimizing Compiler 0.5.1 */
/* ./tc ../examples/rnatools/pbp-fission-1.scop.c --rectangular-tiling --isl-wave-scheduling --omp-for-codegen --isl-map-tc --align --debug -b 16 */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#pragma scop
for (register int k0 = 0; k0 <= floord(_PB_N, 8); k0 += 1) {
  #pragma omp parallel for
  for (register int k1 = max(0, k0 - (_PB_N + 16) / 16 + 1); k1 <= min(k0, _PB_N / 16); k1 += 1) {
    for (register int i0 = max(1, 16 * k1); i0 <= min(_PB_N, 16 * k1 + 15); i0 += 1) {
      for (register int i1 = max(1, 16 * k0 - 16 * k1); i1 <= min(_PB_N, 16 * k0 - 16 * k1 + 15); i1 += 1) {
        Pbp[i0][i1] = (((Q[1][i0 - 1] * Qbp[i0][i1]) * Q[i1 + 1][_PB_N]) / Q[1][_PB_N]);
      }
    }
  }
}
#pragma endscop

/* TC Optimizing Compiler 0.5.1 */
/* ./tc ../examples/rnatools/pbp-fission-2.scop.c --merge-tiling --isl-wave-scheduling --omp-for-codegen --iterative-tc --align --debug -b 16 */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#pragma scop
{
  for (register int k0 = -((_PB_N + 14) / 16) - 1; k0 < floord(_PB_N - 1, 8) - 1; k0 += 1) {
    #pragma omp parallel for
    for (register int k1 = max(0, floord(k0 + 1, 2) + 1); k1 <= min(k0 + (_PB_N - 2) / 16 + 2, _PB_N / 16); k1 += 1) {
      if (16 * k1 + 1 == _PB_N) {
        for (register int k2 = max((-_PB_N + 1) / 16, -((3 * _PB_N - 16 * k0 + 29) / 32) + 1); k2 <= ((-_PB_N + 9) / 8) + k0; k2 += 1) {
          for (register int i1 = max(2 * _PB_N - 16 * k0 + 32 * k2 - 16, 16 * k2); i1 <= 16 * k2 + 15; i1 += 1) {
            for (register int i3 = max(-2 * _PB_N + 16 * k0 - 32 * k2 + 2, -i1 + 1); i3 <= min(_PB_N, -2 * _PB_N + 16 * k0 - 32 * k2 + 17); i3 += 1) {
              Pbp[_PB_N][-i1] += (paired(RNA[_PB_N - 2], RNA[i3 - 1]) * (((((Pbp[_PB_N - 1][i3] * ERT) * Q[_PB_N][_PB_N - 1]) * Qbp[_PB_N][-i1]) * Q[-i1 + 1][i3 - 1]) / ((Qbp[_PB_N - 1][i3] == 0.0) ? 1.0 : Qbp[_PB_N - 1][i3])));
            }
          }
        }
      }
      for (register int k2 = max(max(k0 - 2 * k1 + 1, -((_PB_N + 14) / 16)), k0 - k1 - (_PB_N + 14) / 16 + 2); k2 <= min(-1, k0 - k1 + 1); k2 += 1) {
        for (register int i0 = max(max(2, 16 * k1), 16 * k0 - 16 * k1 - 16 * k2 + 17); i0 <= min(_PB_N, 16 * k1 + 15); i0 += 1) {
          for (register int i1 = max(-_PB_N + 1, 16 * k2); i1 <= 16 * k2 + 15; i1 += 1) {
            for (register int i2 = max(1, 16 * k0 - 16 * k1 - 16 * k2 + 16); i2 <= min(16 * k0 - 16 * k1 - 16 * k2 + 31, i0 - 1); i2 += 1) {
              for (register int i3 = -i1 + 1; i3 <= _PB_N; i3 += 1) {
                Pbp[i0][-i1] += (paired(RNA[i2 - 1], RNA[i3 - 1]) * (((((Pbp[i2][i3] * ERT) * Q[i2 + 1][i0 - 1]) * Qbp[i0][-i1]) * Q[-i1 + 1][i3 - 1]) / ((Qbp[i2][i3] == 0.0) ? 1.0 : Qbp[i2][i3])));
              }
            }
          }
        }
      }
    }
  }
  if ((_PB_N - 1) % 16 == 0) {
    for (register int k0 = (_PB_N - 9) / 8; k0 < (3 * _PB_N - 19) / 16; k0 += 1) {
      for (register int k2 = -((3 * _PB_N - 16 * k0 + 29) / 32) + 1; k2 < 0; k2 += 1) {
        for (register int i1 = 16 * k2; i1 <= 16 * k2 + 15; i1 += 1) {
          for (register int i3 = max(-2 * _PB_N + 16 * k0 - 32 * k2 + 2, -i1 + 1); i3 <= min(_PB_N, -2 * _PB_N + 16 * k0 - 32 * k2 + 17); i3 += 1) {
            Pbp[_PB_N][-i1] += (paired(RNA[_PB_N - 2], RNA[i3 - 1]) * (((((Pbp[_PB_N - 1][i3] * ERT) * Q[_PB_N][_PB_N - 1]) * Qbp[_PB_N][-i1]) * Q[-i1 + 1][i3 - 1]) / ((Qbp[_PB_N - 1][i3] == 0.0) ? 1.0 : Qbp[_PB_N - 1][i3])));
          }
        }
      }
    }
  }
}
#pragma endscop

/* TC Optimizing Compiler 0.5.1 */
/* ./tc ../examples/rnatools/pu.scop.c --merge-tiling --sfs-single-scheduling --omp-for-codegen --floyd-warshall-tc --debug --align -b 16 */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#pragma scop
#pragma omp parallel for
for (register int ir0 = 1; ir0 <= _PB_N; ir0 += 1) {
  for (register int ir1 = 1; ir1 <= _PB_N; ir1 += 1) {
    Pu[ir0][ir1] = ((Q[1][ir0 - 1] * Q[ir1 + 1][_PB_N]) / Q[1][_PB_N]);
    if (ir0 >= 2) {
      for (register int ii3 = 0; ii3 <= (ir0 - 1) / 16; ii3 += 1) {
        if (_PB_N >= 16 * ii3 + 2) {
          for (register int i3 = max(1, 16 * ii3); i3 <= min(ir0 - 1, 16 * ii3 + 15); i3 += 1) {
            for (register int i4 = ir1 + 1; i4 <= _PB_N; i4 += 1) {
              Pu[ir0][ir1] += (paired(RNA[i3 - 1], RNA[i4 - 1]) * ((((Pbp[i3][i4] * ERT) * Q[i3 + 1][ir0 - 1]) * Q[ir1 + 1][i4 - 1]) / ((Qbp[i3][i4] == 0.0) ? 1.0 : Qbp[i3][i4])));
            }
          }
        } else {
          for (register int ii4 = (ir1 + 1) / 16; ii4 <= (_PB_N - 1) / 16; ii4 += 1) {
            for (register int i4 = max(ir1 + 1, 16 * ii4); i4 <= min(_PB_N, 16 * ii4 + 15); i4 += 1) {
              Pu[_PB_N][ir1] += (paired(RNA[_PB_N - 2], RNA[i4 - 1]) * ((((Pbp[_PB_N - 1][i4] * ERT) * Q[_PB_N][_PB_N - 1]) * Q[ir1 + 1][i4 - 1]) / ((Qbp[_PB_N - 1][i4] == 0.0) ? 1.0 : Qbp[_PB_N - 1][i4])));
            }
          }
        }
      }
    }
  }
}
#pragma endscop

/* TC Optimizing Compiler 0.5.1 */
/* ./tc ../examples/rnatools/mea.scop.c --correction-tiling --isl-wave-scheduling --omp-for-codegen --iterative-tc --debug --align -b 16 */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#pragma scop
if (l >= 0 && l <= 5) {
  for (register int k0 = -1; k0 < _PB_N / 16; k0 += 1) {
    #pragma omp parallel for
    for (register int k1 = max(k0 - (_PB_N + 16) / 16 + 1, -((_PB_N + 15) / 16)); k1 < 0; k1 += 1) {
      if (_PB_N + 16 * k1 + 14 >= l) {
        if (k1 <= -2) {
          for (register int i0 = max(max(-_PB_N, -16 * k0 + 16 * k1 - 15), 16 * k1); i0 <= 16 * k1 + 15; i0 += 1) {
            if (16 * k0 + i0 >= l + 16 * k1 + 2) {
              M[-i0][16 * k0 - 16 * k1] = (M[-i0][16 * k0 - 16 * k1 - 1] + Pu[16 * k0 - 16 * k1][16 * k0 - 16 * k1]);
            }
            for (register int i1 = max(16 * k0 - 16 * k1, -i0); i1 <= min(min(_PB_N, 16 * k0 - 16 * k1 + 15), l - i0 + 1); i1 += 1) {
              M[-i0][i1] = (M[-i0][i1 - 1] + Pu[i1][i1]);
            }
          }
        }
        for (register int k3 = -k1 - 1; k3 <= min(k0 - k1, (_PB_N - l - 1) / 16); k3 += 1) {
          if (k1 == -1 && k3 == 0) {
            for (register int i0 = max(max(-16, -_PB_N), -16 * k0 - 31); i0 < 0; i0 += 1) {
              for (register int i1 = max(16 * k0 + 16, -i0); i1 <= min(min(_PB_N, 16 * k0 + 31), l - i0 + 1); i1 += 1) {
                M[-i0][i1] = (M[-i0][i1 - 1] + Pu[i1][i1]);
              }
              if (16 * k0 + i0 + 14 >= l) {
                M[-i0][16 * k0 + 16] = (M[-i0][16 * k0 + 15] + Pu[16 * k0 + 16][16 * k0 + 16]);
              }
            }
          }
          for (register int i0 = max(max(max(-_PB_N + l + 1, l - 16 * k0 + 16 * k1 - 14), 16 * k1), -16 * k3 - 15); i0 <= 16 * k1 + 15; i0 += 1) {
            for (register int i1 = max(max(16 * k0 - 16 * k1, l + 16 * k3 + 1), l - i0 + 1); i1 <= min(min(_PB_N, 16 * k0 - 16 * k1 + 15), l + 16 * k3 + 17); i1 += 1) {
              if (i1 >= l + 16 * k3 + 2 && i0 + i1 >= l + 2 && 16 * k1 + i1 >= 16 * k0 + 1) {
                M[-i0][i1] = (M[-i0][i1 - 1] + Pu[i1][i1]);
              }
              {
                if (i1 >= l + 16 * k3 + 2 && 16 * k1 + i1 >= 16 * k0 + 1) {
                  for (register int i3 = -i0; i3 < 16 * k3; i3 += 1) {
                    M[-i0][i1] = max_score(M[-i0][i1], paired(RNA[i3 - 1], RNA[i1 - 1]) * ((M[-i0][i3 - 1] + M[i3 + 1][i1 - 1]) + (gamma * Pbp[i3][i1])));
                  }
                }
                for (register int i3 = max(16 * k3, -i0); i3 <= min(16 * k3 + 15, -l + i1 - 1); i3 += 1) {
                  M[-i0][i1] = max_score(M[-i0][i1], paired(RNA[i3 - 1], RNA[i1 - 1]) * ((M[-i0][i3 - 1] + M[i3 + 1][i1 - 1]) + (gamma * Pbp[i3][i1])));
                }
              }
            }
            if (k0 >= k1 + k3 + 2) {
              for (register int i3 = max(16 * k3, -i0); i3 <= 16 * k3 + 15; i3 += 1) {
                M[-i0][16 * k0 - 16 * k1] = max_score(M[-i0][16 * k0 - 16 * k1], paired(RNA[i3 - 1], RNA[16 * k0 - 16 * k1 - 1]) * ((M[-i0][i3 - 1] + M[i3 + 1][16 * k0 - 16 * k1 - 1]) + (gamma * Pbp[i3][16 * k0 - 16 * k1])));
              }
            }
          }
        }
      } else {
        for (register int i0 = -_PB_N; i0 <= 16 * k1 + 15; i0 += 1) {
          for (register int i1 = -i0; i1 <= _PB_N; i1 += 1) {
            M[-i0][i1] = (M[-i0][i1 - 1] + Pu[i1][i1]);
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

  int l           = 1;
  DATA_TYPE Ebp   = SCALAR_VAL(-1.0);
  DATA_TYPE RT    = SCALAR_VAL(1.0);
  DATA_TYPE gamma = SCALAR_VAL(2.0);
  DATA_TYPE ERT   = EXP_FUN(-Ebp/RT);

  /* Variable declaration/allocation. */
  POLYBENCH_1D_ARRAY_DECL(RNA, char, N, n);
  POLYBENCH_2D_ARRAY_DECL(Q, DATA_TYPE, N, N, n, n);
  POLYBENCH_2D_ARRAY_DECL(Qbp, DATA_TYPE, N, N, n, n);
  POLYBENCH_2D_ARRAY_DECL(Pbp, DATA_TYPE, N, N, n, n);
  POLYBENCH_2D_ARRAY_DECL(Pu, DATA_TYPE, N, N, n, n);
  POLYBENCH_2D_ARRAY_DECL(M, DATA_TYPE, N, N, n, n);

  /* Initialize array(s). */
  init_array (n, l, ERT, POLYBENCH_ARRAY(RNA), POLYBENCH_ARRAY(Q), POLYBENCH_ARRAY(Qbp), POLYBENCH_ARRAY(Pbp), POLYBENCH_ARRAY(Pu), POLYBENCH_ARRAY(M));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_full_mea (n, l, ERT, gamma, POLYBENCH_ARRAY(RNA), POLYBENCH_ARRAY(Q), POLYBENCH_ARRAY(Qbp), POLYBENCH_ARRAY(Pbp), POLYBENCH_ARRAY(Pu), POLYBENCH_ARRAY(M));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(M)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(RNA);
  POLYBENCH_FREE_ARRAY(Q);
  POLYBENCH_FREE_ARRAY(Qbp);
  POLYBENCH_FREE_ARRAY(Pbp);
  POLYBENCH_FREE_ARRAY(Pu);
  POLYBENCH_FREE_ARRAY(M);

  return 0;
}

