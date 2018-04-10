/**
 * This version is stamped on Apr. 14, 2015
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* nussinov.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "nussinov.h"

/* RNA bases represented as chars, range is [0,3] */
typedef char base;

#define match(b1, b2) (((b1)+(b2)) == 3 ? 1 : 0)
#define max_score(s1, s2) ((s1 >= s2) ? s1 : s2)

/* Array initialization. */
static
void init_array (int n,
                 base POLYBENCH_1D(seq,N,n),
		 DATA_TYPE POLYBENCH_2D(table,N,N,n,n))
{
  int i, j;

  //base is AGCT/0..3
  for (i=0; i <n; i++) {
     seq[i] = (base)((i+1)%4);
  }

  for (i=0; i <n; i++)
     for (j=0; j <n; j++)
       table[i][j] = 0;
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_2D(table,N,N,n,n))

{
  int i, j;
  int t = 0;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("table");
  for (i = 0; i < n; i++) {
    for (j = i; j < n; j++) {
      if (t % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
      fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, table[i][j]);
      t++;
    }
  }
  POLYBENCH_DUMP_END("table");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
/*
  Original version by Dave Wonnacott at Haverford College <davew@cs.haverford.edu>, 
  with help from Allison Lake, Ting Zhou, and Tian Jin, 
  based on algorithm by Nussinov, described in Allison Lake's senior thesis.
*/
static
void kernel_nussinov(int n, base POLYBENCH_1D(seq,N,n),
			   DATA_TYPE POLYBENCH_2D(table,N,N,n,n))
{
  int i, j, k;

/* TC Optimizing Compiler 0.2.26 */
/* ./tc ../examples/polybench/nussinov.scop.c --correction-tiling --lex-scheduling --serial-codegen -b 64 --debug --iterative-tc */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#pragma scop
if (_PB_N >= 2) {
  for (int ii0 = 0; ii0 <= (_PB_N - 2) / 64; ii0 += 1) {
    if (_PB_N >= 3) {
      if (ii0 >= 1 && _PB_N >= 64 * ii0 + 65) {
        for (int ii2 = 0; ii2 <= min(2, -ii0 + (_PB_N - 2 * ii0 - 4) / 62 + 1); ii2 += 1) {
          if (_PB_N >= 64 * ii0 + ii2 + 64 && ii2 >= 1) {
            for (int i0 = -63; i0 <= -_PB_N + 64 * ii0 + 65; i0 += 1) {
              if (ii2 == 2) {
                table[-i0][-i0 + 1] = max_score(table[-i0][-i0 + 1], table[-i0 + 1][-i0]);
                if (i0 >= -62) {
                  table[-i0][-i0 + 2] = max_score(table[-i0][-i0 + 2], table[-i0][-i0 + 1]);
                }
              }
              if (i0 + 64 >= ii2) {
                table[-i0][ii2 - i0] = max_score(table[-i0][ii2 - i0], table[-i0 + 1][ii2 - i0]);
                if (ii2 == 2) {
                  table[-i0][-i0 + 2] = max_score(table[-i0][-i0 + 2], table[-i0 + 1][-i0 + 1] + match(seq[-i0], seq[-i0 + 2]));
                }
              }
            }
          } else {
            for (int i0 = -63; i0 <= -_PB_N + 64 * ii0 + 65; i0 += 1) {
              table[-i0][-i0 + 1] = max_score(table[-i0][-i0 + 1], table[-i0][-i0]);
            }
          }
        }
        for (int ii2 = -ii0 + (_PB_N - 2 * ii0 - 4) / 62 + 2; ii2 <= 3; ii2 += 1) {
          if (ii2 == 3) {
            for (int i0 = -62; i0 <= -_PB_N + 64 * ii0 + 65; i0 += 1) {
              for (int i1 = -i0 + 2; i1 <= 64; i1 += 1) {
                if (_PB_N + i0 >= 64 * ii0 + 5 && i0 + i1 >= 3) {
                  table[-i0][i1] = max_score(table[-i0][i1], table[-i0][i1 - 1]);
                } else if (64 * ii0 + 65 == _PB_N && i0 == -61 && i1 == 64) {
                  table[61][64] = max_score(table[61][64], table[61][63]);
                }
                if (i0 + i1 >= 3) {
                  table[-i0][i1] = max_score(table[-i0][i1], table[-i0 + 1][i1]);
                  table[-i0][i1] = max_score(table[-i0][i1], table[-i0 + 1][i1 - 1] + match(seq[-i0], seq[i1]));
                }
                for (int i3 = -i0 + 1; i3 < i1; i3 += 1) {
                  table[-i0][i1] = max_score(table[-i0][i1], table[-i0][i3] + table[i3 + 1][i1]);
                }
              }
            }
          } else {
            for (int i0 = -63; i0 <= 0; i0 += 1) {
              table[-i0][-i0 + 1] = max_score(table[-i0][-i0 + 1], table[-i0 + 1][-i0]);
              if (i0 >= -62) {
                table[-i0][-i0 + 2] = max_score(table[-i0][-i0 + 2], table[-i0][-i0 + 1]);
                table[-i0][-i0 + 2] = max_score(table[-i0][-i0 + 2], table[-i0 + 1][-i0 + 2]);
                table[-i0][-i0 + 2] = max_score(table[-i0][-i0 + 2], table[-i0 + 1][-i0 + 1] + match(seq[-i0], seq[-i0 + 2]));
              }
            }
          }
        }
      } else if (ii0 >= 1 && 64 * ii0 + 64 >= _PB_N) {
        for (int ii2 = 0; ii2 <= 3; ii2 += 1) {
          if (ii2 >= 0 && ii2 <= 2) {
            int ii3 = ii2 <= 1 ? 0 : 1;
            if (ii2 == 2 && ii3 == 1) {
              table[_PB_N - 64 * ii0 - 2][_PB_N - 64 * ii0 - 1] = max_score(table[_PB_N - 64 * ii0 - 2][_PB_N - 64 * ii0 - 1], table[_PB_N - 64 * ii0 - 1][_PB_N - 64 * ii0 - 2]);
              table[_PB_N - 64 * ii0 - 2][_PB_N - 64 * ii0] = max_score(table[_PB_N - 64 * ii0 - 2][_PB_N - 64 * ii0], table[_PB_N - 64 * ii0 - 2][_PB_N - 64 * ii0 - 1]);
              table[_PB_N - 64 * ii0 - 2][_PB_N - 64 * ii0] = max_score(table[_PB_N - 64 * ii0 - 2][_PB_N - 64 * ii0], table[_PB_N - 64 * ii0 - 1][_PB_N - 64 * ii0]);
              table[_PB_N - 64 * ii0 - 2][_PB_N - 64 * ii0] = max_score(table[_PB_N - 64 * ii0 - 2][_PB_N - 64 * ii0], table[_PB_N - 64 * ii0 - 1][_PB_N - 64 * ii0 - 1] + match(seq[_PB_N - 64 * ii0 - 2], seq[_PB_N - 64 * ii0]));
            }
            for (int i0 = -_PB_N + 64 * ii0 + ii2 + 1; i0 <= min(0, -_PB_N + 64 * ii0 + 62 * ii2 - 60); i0 += 1) {
              if (ii2 == 2) {
                table[-i0][-i0 + 1] = max_score(table[-i0][-i0 + 1], table[-i0 + 1][-i0]);
                table[-i0][-i0 + 2] = max_score(table[-i0][-i0 + 2], table[-i0][-i0 + 1]);
              }
              table[-i0][ii2 - i0] = max_score(table[-i0][ii2 - i0], table[-i0 + 1][ii2 - i0]);
              if (ii2 == 2) {
                table[-i0][-i0 + 2] = max_score(table[-i0][-i0 + 2], table[-i0 + 1][-i0 + 1] + match(seq[-i0], seq[-i0 + 2]));
              }
            }
            if (ii3 == 0) {
              for (int i0 = max(-_PB_N + 64 * ii0 + 2, -_PB_N + 64 * ii0 + 62 * ii2 - 59); i0 <= 0; i0 += 1) {
                if (ii2 == 1) {
                  table[-i0][-i0 + 1] = max_score(table[-i0][-i0 + 1], table[-i0 + 1][-i0 + 1]);
                } else {
                  table[-i0][-i0 + 1] = max_score(table[-i0][-i0 + 1], table[-i0][-i0]);
                }
              }
            }
          } else {
            for (int i0 = -_PB_N + 64 * ii0 + 2; i0 <= 0; i0 += 1) {
              for (int i1 = -i0 + 2; i1 <= 64; i1 += 1) {
                if (_PB_N + i0 >= 64 * ii0 + 3 && i0 + i1 >= 3 && _PB_N + 4 * i0 + 3 * i1 >= 64 * ii0 + 15) {
                  table[-i0][i1] = max_score(table[-i0][i1], table[-i0][i1 - 1]);
                } else if (_PB_N + i0 == 64 * ii0 + 2 && 64 * ii0 + i1 >= _PB_N + 2) {
                  table[_PB_N - 64 * ii0 - 2][i1] = max_score(table[_PB_N - 64 * ii0 - 2][i1], table[_PB_N - 64 * ii0 - 2][i1 - 1]);
                }
                if (64 * ii0 + 3 >= _PB_N + i0 && i0 + i1 >= 5) {
                  table[-i0][i1] = max_score(table[-i0][i1], table[-i0 + 1][i1]);
                } else {
                  if (_PB_N + i0 >= 64 * ii0 + 3 && 64 * ii0 + 5 >= _PB_N + i0 && i0 + i1 == 3) {
                    table[-i0][-i0 + 3] = max_score(table[-i0][-i0 + 3], table[-i0][-i0 + 2]);
                  } else if (_PB_N + i0 == 64 * ii0 + 2 && 64 * ii0 + i1 == _PB_N + 1) {
                    table[_PB_N - 64 * ii0 - 2][_PB_N - 64 * ii0 + 1] = max_score(table[_PB_N - 64 * ii0 - 2][_PB_N - 64 * ii0 + 1], table[_PB_N - 64 * ii0 - 2][_PB_N - 64 * ii0]);
                  }
                  if (i0 + i1 >= 3 && i0 + i1 <= 4) {
                    table[-i0][i1] = max_score(table[-i0][i1], table[-i0 + 1][i1]);
                  } else if (_PB_N + i0 >= 64 * ii0 + 4 && i0 + i1 >= 5) {
                    table[-i0][i1] = max_score(table[-i0][i1], table[-i0 + 1][i1]);
                  }
                }
                if (i0 + i1 >= 3 && 64 * ii0 >= i1) {
                  table[-i0][i1] = max_score(table[-i0][i1], table[-i0 + 1][i1 - 1] + match(seq[-i0], seq[i1]));
                }
                for (int i3 = -i0 + 1; i3 < i1; i3 += 1) {
                  table[-i0][i1] = max_score(table[-i0][i1], table[-i0][i3] + table[i3 + 1][i1]);
                }
              }
            }
          }
        }
      }
      for (int ii1 = max(1, -ii0 + (_PB_N - 1) / 64 - 1); ii1 < ii0; ii1 += 1) {
        for (int ii2 = 0; ii2 <= min(3, -_PB_N + 32 * ii0 + 32 * ii1 + (_PB_N + 1) / 2 + 66); ii2 += 1) {
          if (_PB_N >= 64 * ii0 + 64 * ii1 + 2 && ii2 >= 1 && ii2 <= 2) {
            for (int i0 = max(-_PB_N + 64 * ii0 + 2, -64 * ii1 - 63); i0 <= min(-_PB_N + 64 * ii0 + 65, -64 * ii1 + ii2 - 1); i0 += 1) {
              if (ii2 == 2 && 64 * ii1 + i0 <= 0) {
                table[-i0][-i0 + 1] = max_score(table[-i0][-i0 + 1], table[-i0 + 1][-i0]);
                if (64 * ii1 + i0 >= -62) {
                  table[-i0][-i0 + 2] = max_score(table[-i0][-i0 + 2], table[-i0][-i0 + 1]);
                }
              }
              if (64 * ii1 + i0 + 64 >= ii2) {
                table[-i0][ii2 - i0] = max_score(table[-i0][ii2 - i0], table[-i0 + 1][ii2 - i0]);
                if (ii2 == 2) {
                  table[-i0][-i0 + 2] = max_score(table[-i0][-i0 + 2], table[-i0 + 1][-i0 + 1] + match(seq[-i0], seq[-i0 + 2]));
                }
              }
            }
          } else if (64 * ii0 + 64 * ii1 + 1 >= _PB_N && ii2 >= 1 && ii2 <= 2) {
            if (ii2 == 2) {
              table[_PB_N - 64 * ii0 - 2][64 * ii1 + 1] = max_score(table[_PB_N - 64 * ii0 - 2][64 * ii1 + 1], table[_PB_N - 64 * ii0 - 1][64 * ii1] + match(seq[_PB_N - 64 * ii0 - 2], seq[64 * ii1 + 1]));
            } else if (64 * ii0 + 64 >= _PB_N) {
              table[_PB_N - 64 * ii0 - 2][64 * ii1 + 1] = max_score(table[_PB_N - 64 * ii0 - 2][64 * ii1 + 1], table[_PB_N - 64 * ii0 - 1][64 * ii1 + 1]);
            } else {
              table[_PB_N - 64 * ii0 - 2][64 * ii1 + 1] = max_score(table[_PB_N - 64 * ii0 - 2][64 * ii1 + 1], table[_PB_N - 64 * ii0 - 1][64 * ii1 + 1]);
            }
          } else if (ii2 == 3) {
            for (int ii3 = -ii0 + (_PB_N - 1) / 64; ii3 < ii1; ii3 += 1) {
              for (int i3 = max(_PB_N - 64 * ii0 - 1, 64 * ii3); i3 <= 64 * ii3 + 63; i3 += 1) {
                table[_PB_N - 64 * ii0 - 2][64 * ii1 + 1] = max_score(table[_PB_N - 64 * ii0 - 2][64 * ii1 + 1], table[_PB_N - 64 * ii0 - 2][i3] + table[i3 + 1][64 * ii1 + 1]);
              }
            }
            if (_PB_N >= 64 * ii0 + 65) {
              for (int i0 = max(-_PB_N + 64 * ii0 + 2, -64 * ii1 - 62); i0 <= -_PB_N + 64 * ii0 + 65; i0 += 1) {
                for (int i1 = max(64 * ii1 + 1, -i0 + 2); i1 <= 64 * ii1 + 64; i1 += 1) {
                  if (i0 + i1 >= 3 && i1 >= 64 * ii1 + 2 && (_PB_N - 63 * i0 - 4) % 64 >= 1) {
                    table[-i0][i1] = max_score(table[-i0][i1], table[-i0][i1 - 1]);
                  }
                  if (_PB_N + i0 >= 64 * ii0 + 5 && i0 + i1 >= 5) {
                    table[-i0][i1] = max_score(table[-i0][i1], table[-i0 + 1][i1]);
                  } else if (64 * ii0 + 3 >= _PB_N + i0 && _PB_N + i0 + i1 >= 64 * ii0 + 64 * ii1 + 4 && 64 * ii0 + i1 >= _PB_N + 2) {
                    table[-i0][i1] = max_score(table[-i0][i1], table[-i0 + 1][i1]);
                  } else {
                    if (_PB_N + i0 == 64 * ii0 + 4 && i1 >= 64 * ii1 + 2 && 64 * ii0 + i1 + 1 >= _PB_N) {
                      table[_PB_N - 64 * ii0 - 4][i1] = max_score(table[_PB_N - 64 * ii0 - 4][i1], table[_PB_N - 64 * ii0 - 4][i1 - 1]);
                    }
                    if (_PB_N + i0 + i1 >= 64 * ii0 + 64 * ii1 + 4 && i0 + i1 >= 3 && _PB_N + 1 >= 64 * ii0 + i1 && i0 + i1 <= 4) {
                      table[-i0][i1] = max_score(table[-i0][i1], table[-i0 + 1][i1]);
                    } else if (_PB_N + i0 == 64 * ii0 + 4 && 64 * ii0 + i1 >= _PB_N + 1) {
                      table[_PB_N - 64 * ii0 - 4][i1] = max_score(table[_PB_N - 64 * ii0 - 4][i1], table[_PB_N - 64 * ii0 - 3][i1]);
                    }
                  }
                  if (_PB_N + i0 >= 64 * ii0 + 4 && i0 + i1 >= 3) {
                    table[-i0][i1] = max_score(table[-i0][i1], table[-i0 + 1][i1 - 1] + match(seq[-i0], seq[i1]));
                  } else if (64 * ii0 + 3 >= _PB_N + i0 && _PB_N + i0 + i1 >= 64 * ii0 + 64 * ii1 + 4 && i0 + i1 >= 3) {
                    if (64 * ii0 + 64 * ii1 >= _PB_N && _PB_N + i0 == 64 * ii0 + 3) {
                      table[_PB_N - 64 * ii0 - 3][i1] = max_score(table[_PB_N - 64 * ii0 - 3][i1], table[_PB_N - 64 * ii0 - 2][i1 - 1] + match(seq[_PB_N - 64 * ii0 - 3], seq[i1]));
                    } else if (_PB_N >= 64 * ii0 + 64 * ii1 + 1 && _PB_N + i0 == 64 * ii0 + 3 && 64 * ii0 + i1 >= _PB_N + 3) {
                      table[_PB_N - 64 * ii0 - 3][i1] = max_score(table[_PB_N - 64 * ii0 - 3][i1], table[_PB_N - 64 * ii0 - 2][i1 - 1] + match(seq[_PB_N - 64 * ii0 - 3], seq[i1]));
                    } else if (_PB_N >= 64 * ii0 + 64 * ii1 + 1 && _PB_N + i0 == 64 * ii0 + 3 && _PB_N + 2 >= 64 * ii0 + i1) {
                      table[_PB_N - 64 * ii0 - 3][i1] = max_score(table[_PB_N - 64 * ii0 - 3][i1], table[_PB_N - 64 * ii0 - 2][i1 - 1] + match(seq[_PB_N - 64 * ii0 - 3], seq[i1]));
                    } else {
                      table[_PB_N - 64 * ii0 - 2][i1] = max_score(table[_PB_N - 64 * ii0 - 2][i1], table[_PB_N - 64 * ii0 - 1][i1 - 1] + match(seq[_PB_N - 64 * ii0 - 2], seq[i1]));
                    }
                  }
                  if (_PB_N + i0 + i1 >= 64 * ii0 + 64 * ii1 + 4) {
                    for (int i3 = -i0 + 1; i3 < 64 * ii1; i3 += 1) {
                      table[-i0][i1] = max_score(table[-i0][i1], table[-i0][i3] + table[i3 + 1][i1]);
                    }
                  }
                  for (int i3 = max(64 * ii1, -i0 + 1); i3 < i1; i3 += 1) {
                    table[-i0][i1] = max_score(table[-i0][i1], table[-i0][i3] + table[i3 + 1][i1]);
                  }
                }
              }
            } else {
              for (int i0 = -_PB_N + 64 * ii0 + 2; i0 <= 0; i0 += 1) {
                for (int i1 = 64 * ii1 + 1; i1 <= 64 * ii1 + 64; i1 += 1) {
                  if (i1 >= 64 * ii1 + 2) {
                    table[-i0][i1] = max_score(table[-i0][i1], table[-i0][i1 - 1]);
                  }
                  if (_PB_N + i0 + i1 >= 64 * ii0 + 64 * ii1 + 4) {
                    table[-i0][i1] = max_score(table[-i0][i1], table[-i0 + 1][i1]);
                    table[-i0][i1] = max_score(table[-i0][i1], table[-i0 + 1][i1 - 1] + match(seq[-i0], seq[i1]));
                  }
                  if (_PB_N + i0 + i1 >= 64 * ii0 + 64 * ii1 + 4) {
                    for (int i3 = -i0 + 1; i3 < 64 * ii1; i3 += 1) {
                      table[-i0][i1] = max_score(table[-i0][i1], table[-i0][i3] + table[i3 + 1][i1]);
                    }
                  }
                  for (int i3 = 64 * ii1; i3 < i1; i3 += 1) {
                    table[-i0][i1] = max_score(table[-i0][i1], table[-i0][i3] + table[i3 + 1][i1]);
                  }
                }
              }
            }
          } else if (_PB_N >= 64 * ii0 + 65) {
            for (int i0 = max(-_PB_N + 64 * ii0 + 2, -64 * ii1 - 63); i0 <= min(-_PB_N + 64 * ii0 + 65, -64 * ii1); i0 += 1) {
              table[-i0][-i0 + 1] = max_score(table[-i0][-i0 + 1], table[-i0][-i0]);
            }
            for (int i0 = max(-_PB_N + 64 * ii0 + 2, -64 * ii1 + 1); i0 <= -_PB_N + 64 * ii0 + 65; i0 += 1) {
              table[-i0][64 * ii1 + 1] = max_score(table[-i0][64 * ii1 + 1], table[-i0][64 * ii1]);
            }
          } else {
            for (int i0 = -_PB_N + 64 * ii0 + 2; i0 <= 0; i0 += 1) {
              if (_PB_N + i0 >= 64 * ii0 + 3) {
                table[-i0][64 * ii1 + 1] = max_score(table[-i0][64 * ii1 + 1], table[-i0][64 * ii1]);
              } else {
                table[_PB_N - 64 * ii0 - 2][64 * ii1 + 1] = max_score(table[_PB_N - 64 * ii0 - 2][64 * ii1 + 1], table[_PB_N - 64 * ii0 - 2][64 * ii1]);
              }
            }
          }
        }
      }
      if (_PB_N >= 64 * ii0 + 65) {
        for (int ii1 = max(ii0, -ii0 + (_PB_N - 1) / 64 - 1); ii1 <= (_PB_N - 2) / 64; ii1 += 1) {
          for (int i0 = max(-_PB_N + 64 * ii0 + 2, -64 * ii1 - 63); i0 <= min(-_PB_N + 64 * ii0 + 65, -64 * ii1); i0 += 1) {
            table[-i0][-i0 + 1] = max_score(table[-i0][-i0 + 1], table[-i0][-i0]);
          }
          for (int i0 = max(-_PB_N + 64 * ii0 + 2, -64 * ii1 + 1); i0 <= -_PB_N + 64 * ii0 + 65; i0 += 1) {
            table[-i0][64 * ii1 + 1] = max_score(table[-i0][64 * ii1 + 1], table[-i0][64 * ii1]);
          }
          if (128 * ii0 + 128 == _PB_N && 128 * ii1 + 128 == _PB_N) {
            table[(_PB_N / 2) - 1][_PB_N / 2] = max_score(table[(_PB_N / 2) - 1][_PB_N / 2], table[_PB_N / 2][_PB_N / 2]);
          } else if (128 * ii0 + 2 == _PB_N && 128 * ii1 + 2 == _PB_N) {
            table[(_PB_N / 2) - 1][_PB_N / 2] = max_score(table[(_PB_N / 2) - 1][_PB_N / 2], table[_PB_N / 2][_PB_N / 2]);
          } else if (ii0 == 0 && 64 * ii1 + 64 >= _PB_N) {
            for (int ii2 = 1; ii2 <= 2; ii2 += 1) {
              if (64 * ii1 + 63 >= _PB_N && ii2 == 2) {
                table[_PB_N - 2][_PB_N - 1] = max_score(table[_PB_N - 2][_PB_N - 1], table[_PB_N - 1][_PB_N - 2]);
              }
              for (int i0 = -_PB_N + 2; i0 <= -64 * ii1 - 61 * ii2 + 60; i0 += 1) {
                if (ii2 == 1) {
                  table[-i0][-i0 + 1] = max_score(table[-i0][-i0 + 1], table[-i0 + 1][-i0 + 1]);
                } else {
                  table[_PB_N - 2][_PB_N - 1] = max_score(table[_PB_N - 2][_PB_N - 1], table[_PB_N - 1][_PB_N - 2]);
                }
              }
              for (int i0 = max(-_PB_N + ii2 + 1, -64 * ii1 - 61 * ii2 + 61); i0 < -64 * ii1 + ii2; i0 += 1) {
                if (ii2 == 2 && 64 * ii1 + i0 <= 0) {
                  table[-i0][-i0 + 1] = max_score(table[-i0][-i0 + 1], table[-i0 + 1][-i0]);
                  if (_PB_N + i0 >= 7) {
                    table[-i0][-i0 + 2] = max_score(table[-i0][-i0 + 2], table[-i0][-i0 + 1]);
                  } else if (_PB_N + i0 <= 5) {
                    table[-i0][-i0 + 2] = max_score(table[-i0][-i0 + 2], table[-i0][-i0 + 1]);
                  } else {
                    table[_PB_N - 6][_PB_N - 4] = max_score(table[_PB_N - 6][_PB_N - 4], table[_PB_N - 6][_PB_N - 5]);
                  }
                }
                table[-i0][ii2 - i0] = max_score(table[-i0][ii2 - i0], table[-i0 + 1][ii2 - i0]);
                if (ii2 == 2) {
                  table[-i0][-i0 + 2] = max_score(table[-i0][-i0 + 2], table[-i0 + 1][-i0 + 1] + match(seq[-i0], seq[-i0 + 2]));
                }
              }
            }
          } else if (ii0 >= 1 && _PB_N >= 128 * ii0 + 66 && 64 * ii0 + 64 * ii1 + 2 == _PB_N) {
            table[_PB_N - 64 * ii0 - 2][_PB_N - 64 * ii0 - 1] = max_score(table[_PB_N - 64 * ii0 - 2][_PB_N - 64 * ii0 - 1], table[_PB_N - 64 * ii0 - 1][_PB_N - 64 * ii0 - 1]);
          } else if (_PB_N >= 128 * ii0 + 3 && 128 * ii0 + 4 >= _PB_N && ii1 == ii0) {
            table[_PB_N - 64 * ii0 - 2][_PB_N - 64 * ii0 - 1] = max_score(table[_PB_N - 64 * ii0 - 2][_PB_N - 64 * ii0 - 1], table[_PB_N - 64 * ii0 - 1][_PB_N - 64 * ii0 - 1]);
            for (int i0 = -_PB_N + 64 * ii0 + 3; i0 <= -64 * ii0; i0 += 1) {
              table[-i0][-i0 + 1] = max_score(table[-i0][-i0 + 1], table[-i0 + 1][-i0 + 1]);
            }
          } else if (_PB_N >= 128 * ii0 + 5 && _PB_N >= 64 * ii1 + 65 && _PB_N >= 64 * ii0 + 64 * ii1 + 3) {
            if (64 * ii0 + 64 * ii1 + 65 >= _PB_N) {
              table[_PB_N - 64 * ii0 - 2][_PB_N - 64 * ii0 - 1] = max_score(table[_PB_N - 64 * ii0 - 2][_PB_N - 64 * ii0 - 1], table[_PB_N - 64 * ii0 - 1][_PB_N - 64 * ii0 - 1]);
            }
            for (int i0 = max(-_PB_N + 64 * ii0 + 3, -64 * ii1 - 63); i0 <= min(min(-_PB_N + 64 * ii0 + 65, -64 * ii1), -_PB_N + (_PB_N + 1) / 2); i0 += 1) {
              table[-i0][-i0 + 1] = max_score(table[-i0][-i0 + 1], table[-i0 + 1][-i0 + 1]);
            }
            if (ii1 == ii0) {
              for (int i0 = -_PB_N + (_PB_N + 1) / 2 + 1; i0 <= min(-64 * ii0, -_PB_N + 64 * ii0 + 65); i0 += 1) {
                table[-i0][-i0 + 1] = max_score(table[-i0][-i0 + 1], table[-i0 + 1][-i0 + 1]);
              }
            }
          }
          if (_PB_N >= 64 * ii1 + 65 && _PB_N >= 64 * ii0 + 64 * ii1 + 2) {
            for (int i0 = max(-_PB_N + 64 * ii0 + 2, -64 * ii1 - 63); i0 <= min(-_PB_N + 64 * ii0 + 65, -64 * ii1 + 1); i0 += 1) {
              if (64 * ii1 + i0 <= 0) {
                table[-i0][-i0 + 1] = max_score(table[-i0][-i0 + 1], table[-i0 + 1][-i0]);
                if (64 * ii1 + i0 >= -62) {
                  table[-i0][-i0 + 2] = max_score(table[-i0][-i0 + 2], table[-i0][-i0 + 1]);
                }
              }
              if (64 * ii1 + i0 >= -62) {
                table[-i0][-i0 + 2] = max_score(table[-i0][-i0 + 2], table[-i0 + 1][-i0 + 2]);
                table[-i0][-i0 + 2] = max_score(table[-i0][-i0 + 2], table[-i0 + 1][-i0 + 1] + match(seq[-i0], seq[-i0 + 2]));
              }
            }
            if (64 * ii0 + 64 * ii1 + 127 == _PB_N) {
              table[_PB_N - 64 * ii0 - 65][_PB_N - 64 * ii0 - 63] = max_score(table[_PB_N - 64 * ii0 - 65][_PB_N - 64 * ii0 - 63], table[_PB_N - 64 * ii0 - 65][_PB_N - 64 * ii0 - 64] + table[_PB_N - 64 * ii0 - 63][_PB_N - 64 * ii0 - 63]);
            }
          } else if (64 * ii0 + 64 * ii1 + 1 >= _PB_N) {
            for (int ii2 = 1; ii2 <= 2; ii2 += 1) {
              if (ii2 == 2) {
                table[_PB_N - 64 * ii0 - 2][64 * ii1 + 1] = max_score(table[_PB_N - 64 * ii0 - 2][64 * ii1 + 1], table[_PB_N - 64 * ii0 - 1][64 * ii1] + match(seq[_PB_N - 64 * ii0 - 2], seq[64 * ii1 + 1]));
              } else {
                table[_PB_N - 64 * ii0 - 2][64 * ii1 + 1] = max_score(table[_PB_N - 64 * ii0 - 2][64 * ii1 + 1], table[_PB_N - 64 * ii0 - 1][64 * ii1 + 1]);
              }
            }
          }
          if (64 * ii0 + 64 * ii1 + 126 >= _PB_N) {
            for (int ii3 = -ii0 + (_PB_N - 1) / 64; ii3 < ii1; ii3 += 1) {
              for (int i3 = max(_PB_N - 64 * ii0 - 1, 64 * ii3); i3 <= 64 * ii3 + 63; i3 += 1) {
                table[_PB_N - 64 * ii0 - 2][64 * ii1 + 1] = max_score(table[_PB_N - 64 * ii0 - 2][64 * ii1 + 1], table[_PB_N - 64 * ii0 - 2][i3] + table[i3 + 1][64 * ii1 + 1]);
              }
            }
            if (ii0 >= 1 && 64 * ii1 + 64 >= _PB_N) {
              for (int i0 = -_PB_N + 64 * ii0 + 2; i0 <= -_PB_N + 64 * ii0 + 65; i0 += 1) {
                if (64 * ii1 + 1 >= _PB_N + i0) {
                  for (int i1 = 64 * ii1 + 1; i1 < _PB_N; i1 += 1) {
                    if (i1 >= 64 * ii1 + 2 && (_PB_N - 63 * i0 - 4) % 64 >= 1) {
                      table[-i0][i1] = max_score(table[-i0][i1], table[-i0][i1 - 1]);
                    }
                    if (_PB_N + i0 >= 64 * ii0 + 5) {
                      table[-i0][i1] = max_score(table[-i0][i1], table[-i0 + 1][i1]);
                    } else {
                      if (_PB_N + i0 == 64 * ii0 + 4 && i1 >= 64 * ii1 + 2) {
                        table[_PB_N - 64 * ii0 - 4][i1] = max_score(table[_PB_N - 64 * ii0 - 4][i1], table[_PB_N - 64 * ii0 - 4][i1 - 1]);
                      }
                      if (_PB_N + i0 >= 64 * ii0 + 3 && _PB_N + 2 * i0 + i1 >= 64 * ii0 + 9 && _PB_N + i0 + i1 >= 64 * ii0 + 64 * ii1 + 5) {
                        table[-i0][i1] = max_score(table[-i0][i1], table[-i0 + 1][i1]);
                      } else if (_PB_N + i0 == 64 * ii0 + 2 && i1 >= 64 * ii1 + 2) {
                        table[_PB_N - 64 * ii0 - 2][i1] = max_score(table[_PB_N - 64 * ii0 - 2][i1], table[_PB_N - 64 * ii0 - 1][i1]);
                      } else if (ii0 == 1 && _PB_N + i0 == 67 && _PB_N >= i1 + 62) {
                        table[_PB_N - 67][i1] = max_score(table[_PB_N - 67][i1], table[_PB_N - 66][i1]);
                      } else if (64 * ii0 + 64 * ii1 >= _PB_N + 2 && _PB_N + i0 == 64 * ii0 + 3) {
                        table[_PB_N - 64 * ii0 - 3][64 * ii1 + 1] = max_score(table[_PB_N - 64 * ii0 - 3][64 * ii1 + 1], table[_PB_N - 64 * ii0 - 2][64 * ii1 + 1]);
                      }
                    }
                    if (_PB_N + i0 >= 64 * ii0 + 4 || (i0 + i1 >= 5 && _PB_N + i0 + i1 >= 64 * ii0 + 64 * ii1 + 4) || (ii0 == 1 && 64 * ii1 + 64 == _PB_N && i0 + i1 == 4)) {
                      table[-i0][i1] = max_score(table[-i0][i1], table[-i0 + 1][i1 - 1] + match(seq[-i0], seq[i1]));
                    }
                    if (_PB_N + i0 + i1 >= 64 * ii0 + 64 * ii1 + 4) {
                      for (int i3 = -i0 + 1; i3 < 64 * ii1; i3 += 1) {
                        table[-i0][i1] = max_score(table[-i0][i1], table[-i0][i3] + table[i3 + 1][i1]);
                      }
                    }
                    for (int i3 = 64 * ii1; i3 < i1; i3 += 1) {
                      table[-i0][i1] = max_score(table[-i0][i1], table[-i0][i3] + table[i3 + 1][i1]);
                    }
                  }
                }
              }
            } else if (_PB_N >= 64 * ii1 + 65) {
              if (_PB_N >= 64 * ii0 + 64 * ii1 + 64) {
                table[64 * ii1 + 62][64 * ii1 + 64] = max_score(table[64 * ii1 + 62][64 * ii1 + 64], table[64 * ii1 + 62][64 * ii1 + 63] + table[64 * ii1 + 64][64 * ii1 + 64]);
              }
              for (int i0 = max(-_PB_N + 64 * ii0 + 2, -64 * ii1 - 61); i0 <= -_PB_N + 64 * ii0 + 65; i0 += 1) {
                if (64 * ii1 + i0 <= 0) {
                  table[-i0][-i0 + 2] = max_score(table[-i0][-i0 + 2], table[-i0][-i0 + 1] + table[-i0 + 2][-i0 + 2]);
                } else {
                  if (_PB_N + i0 >= 64 * ii0 + 5 && 64 * ii1 + i0 >= 4) {
                    table[-i0][64 * ii1 + 1] = max_score(table[-i0][64 * ii1 + 1], table[-i0 + 1][64 * ii1 + 1]);
                  } else if (64 * ii1 + i0 >= 3 && _PB_N + i0 >= 64 * ii0 + 3 && 64 * ii0 + 4 >= _PB_N + i0) {
                    table[-i0][64 * ii1 + 1] = max_score(table[-i0][64 * ii1 + 1], table[-i0 + 1][64 * ii1 + 1]);
                  } else if (64 * ii1 + i0 >= 2 && _PB_N + 1 >= 64 * ii0 + 128 * ii1 + i0 && 64 * ii1 + i0 <= 3) {
                    table[-i0][64 * ii1 + 1] = max_score(table[-i0][64 * ii1 + 1], table[-i0 + 1][64 * ii1 + 1]);
                  }
                  if (_PB_N + i0 >= 64 * ii0 + 4 && 64 * ii1 + i0 >= 3) {
                    table[-i0][64 * ii1 + 1] = max_score(table[-i0][64 * ii1 + 1], table[-i0 + 1][64 * ii1] + match(seq[-i0], seq[64 * ii1 + 1]));
                  } else if (64 * ii0 + 64 * ii1 + 1 >= _PB_N && _PB_N + i0 == 64 * ii0 + 3) {
                    table[_PB_N - 64 * ii0 - 3][64 * ii1 + 1] = max_score(table[_PB_N - 64 * ii0 - 3][64 * ii1 + 1], table[_PB_N - 64 * ii0 - 2][64 * ii1] + match(seq[_PB_N - 64 * ii0 - 3], seq[64 * ii1 + 1]));
                  } else if (_PB_N >= 64 * ii0 + 64 * ii1 + 2 && 64 * ii1 + i0 == 2) {
                    table[64 * ii1 - 2][64 * ii1 + 1] = max_score(table[64 * ii1 - 2][64 * ii1 + 1], table[64 * ii1 - 1][64 * ii1] + match(seq[64 * ii1 - 2], seq[64 * ii1 + 1]));
                  }
                  if (_PB_N + i0 >= 64 * ii0 + 3) {
                    for (int i3 = -i0 + 1; i3 < 64 * ii1; i3 += 1) {
                      table[-i0][64 * ii1 + 1] = max_score(table[-i0][64 * ii1 + 1], table[-i0][i3] + table[i3 + 1][64 * ii1 + 1]);
                    }
                  }
                  table[-i0][64 * ii1 + 1] = max_score(table[-i0][64 * ii1 + 1], table[-i0][64 * ii1] + table[64 * ii1 + 1][64 * ii1 + 1]);
                }
                if (_PB_N + i0 == 64 * ii0 + 4) {
                  for (int i1 = max(_PB_N - 64 * ii0 - 1, 64 * ii1 + 2); i1 <= 64 * ii1 + 64; i1 += 1) {
                    table[_PB_N - 64 * ii0 - 4][i1] = max_score(table[_PB_N - 64 * ii0 - 4][i1], table[_PB_N - 64 * ii0 - 4][i1 - 1]);
                    if (64 * ii0 + i1 >= _PB_N + 1) {
                      table[_PB_N - 64 * ii0 - 4][i1] = max_score(table[_PB_N - 64 * ii0 - 4][i1], table[_PB_N - 64 * ii0 - 3][i1]);
                    } else if (64 * ii0 + i1 == _PB_N) {
                      table[_PB_N - 64 * ii0 - 4][_PB_N - 64 * ii0] = max_score(table[_PB_N - 64 * ii0 - 4][_PB_N - 64 * ii0], table[_PB_N - 64 * ii0 - 3][_PB_N - 64 * ii0]);
                    } else {
                      table[_PB_N - 64 * ii0 - 4][_PB_N - 64 * ii0 - 1] = max_score(table[_PB_N - 64 * ii0 - 4][_PB_N - 64 * ii0 - 1], table[_PB_N - 64 * ii0 - 3][_PB_N - 64 * ii0 - 1]);
                    }
                    table[_PB_N - 64 * ii0 - 4][i1] = max_score(table[_PB_N - 64 * ii0 - 4][i1], table[_PB_N - 64 * ii0 - 3][i1 - 1] + match(seq[_PB_N - 64 * ii0 - 4], seq[i1]));
                    for (int i3 = _PB_N - 64 * ii0 - 3; i3 < i1; i3 += 1) {
                      table[_PB_N - 64 * ii0 - 4][i1] = max_score(table[_PB_N - 64 * ii0 - 4][i1], table[_PB_N - 64 * ii0 - 4][i3] + table[i3 + 1][i1]);
                    }
                  }
                } else if ((_PB_N - 63 * i0 - 4) % 64 >= 1) {
                  for (int i1 = max(64 * ii1 + 2, -i0 + 3); i1 <= 64 * ii1 + 64; i1 += 1) {
                    table[-i0][i1] = max_score(table[-i0][i1], table[-i0][i1 - 1]);
                    if (_PB_N + i0 >= 64 * ii0 + 5 && i0 + i1 >= 4) {
                      table[-i0][i1] = max_score(table[-i0][i1], table[-i0 + 1][i1]);
                    } else if (64 * ii0 + 3 >= _PB_N + i0 && 64 * ii0 + i1 >= _PB_N + 2) {
                      table[-i0][i1] = max_score(table[-i0][i1], table[-i0 + 1][i1]);
                    } else if (64 * ii0 + 3 >= _PB_N + i0 && 64 * ii0 + i1 == _PB_N + 1) {
                      table[-i0][_PB_N - 64 * ii0 + 1] = max_score(table[-i0][_PB_N - 64 * ii0 + 1], table[-i0 + 1][_PB_N - 64 * ii0 + 1]);
                    } else {
                      table[-i0][-i0 + 3] = max_score(table[-i0][-i0 + 3], table[-i0 + 1][-i0 + 3]);
                    }
                    if (_PB_N + i0 >= 64 * ii0 + 5) {
                      table[-i0][i1] = max_score(table[-i0][i1], table[-i0 + 1][i1 - 1] + match(seq[-i0], seq[i1]));
                    } else if (_PB_N + i0 == 64 * ii0 + 2) {
                      table[_PB_N - 64 * ii0 - 2][i1] = max_score(table[_PB_N - 64 * ii0 - 2][i1], table[_PB_N - 64 * ii0 - 1][i1 - 1] + match(seq[_PB_N - 64 * ii0 - 2], seq[i1]));
                    } else if (64 * ii0 + 64 * ii1 >= _PB_N) {
                      table[_PB_N - 64 * ii0 - 3][i1] = max_score(table[_PB_N - 64 * ii0 - 3][i1], table[_PB_N - 64 * ii0 - 2][i1 - 1] + match(seq[_PB_N - 64 * ii0 - 3], seq[i1]));
                    } else if (64 * ii0 + i1 >= _PB_N + 3) {
                      table[_PB_N - 64 * ii0 - 3][i1] = max_score(table[_PB_N - 64 * ii0 - 3][i1], table[_PB_N - 64 * ii0 - 2][i1 - 1] + match(seq[_PB_N - 64 * ii0 - 3], seq[i1]));
                    } else {
                      table[_PB_N - 64 * ii0 - 3][i1] = max_score(table[_PB_N - 64 * ii0 - 3][i1], table[_PB_N - 64 * ii0 - 2][i1 - 1] + match(seq[_PB_N - 64 * ii0 - 3], seq[i1]));
                    }
                    for (int i3 = -i0 + 1; i3 < i1; i3 += 1) {
                      table[-i0][i1] = max_score(table[-i0][i1], table[-i0][i3] + table[i3 + 1][i1]);
                    }
                  }
                }
              }
            }
            if (ii0 == 0 && 64 * ii1 + 64 >= _PB_N) {
              table[_PB_N - 3][_PB_N - 1] = max_score(table[_PB_N - 3][_PB_N - 1], table[_PB_N - 3][_PB_N - 2] + table[_PB_N - 1][_PB_N - 1]);
              if (64 * ii1 + 2 == _PB_N) {
                table[_PB_N - 4][_PB_N - 1] = max_score(table[_PB_N - 4][_PB_N - 1], table[_PB_N - 3][_PB_N - 1]);
                table[_PB_N - 4][_PB_N - 1] = max_score(table[_PB_N - 4][_PB_N - 1], table[_PB_N - 3][_PB_N - 2] + match(seq[_PB_N - 4], seq[_PB_N - 1]));
                for (int i3 = _PB_N - 3; i3 < _PB_N - 1; i3 += 1) {
                  table[_PB_N - 4][_PB_N - 1] = max_score(table[_PB_N - 4][_PB_N - 1], table[_PB_N - 4][i3] + table[i3 + 1][_PB_N - 1]);
                }
              } else {
                for (int i1 = _PB_N - 2; i1 < _PB_N; i1 += 1) {
                  if (i1 + 1 == _PB_N) {
                    table[_PB_N - 4][_PB_N - 1] = max_score(table[_PB_N - 4][_PB_N - 1], table[_PB_N - 4][_PB_N - 2]);
                    table[_PB_N - 4][_PB_N - 1] = max_score(table[_PB_N - 4][_PB_N - 1], table[_PB_N - 3][_PB_N - 1]);
                    table[_PB_N - 4][_PB_N - 1] = max_score(table[_PB_N - 4][_PB_N - 1], table[_PB_N - 3][_PB_N - 2] + match(seq[_PB_N - 4], seq[_PB_N - 1]));
                  }
                  for (int i3 = _PB_N - 3; i3 < i1; i3 += 1) {
                    table[_PB_N - 4][i1] = max_score(table[_PB_N - 4][i1], table[_PB_N - 4][i3] + table[i3 + 1][i1]);
                  }
                }
              }
              for (int i0 = -_PB_N + 5; i0 <= -_PB_N + 65; i0 += 1) {
                if (64 * ii1 + 1 >= _PB_N + i0) {
                  for (int i1 = max(64 * ii1 + 1, -i0 + 2); i1 < _PB_N; i1 += 1) {
                    if (i0 + i1 >= 3 && i1 >= 64 * ii1 + 2) {
                      table[-i0][i1] = max_score(table[-i0][i1], table[-i0][i1 - 1]);
                    }
                    if (i0 + i1 >= 3) {
                      table[-i0][i1] = max_score(table[-i0][i1], table[-i0 + 1][i1]);
                      table[-i0][i1] = max_score(table[-i0][i1], table[-i0 + 1][i1 - 1] + match(seq[-i0], seq[i1]));
                    }
                    for (int i3 = -i0 + 1; i3 < i1; i3 += 1) {
                      table[-i0][i1] = max_score(table[-i0][i1], table[-i0][i3] + table[i3 + 1][i1]);
                    }
                  }
                }
              }
            }
          }
        }
      } else {
        if (ii0 == 0) {
          for (int ii2 = 0; ii2 <= 2; ii2 += 1) {
            if (ii2 >= 1) {
              for (int i0 = -_PB_N + 2; i0 <= 0; i0 += 1) {
                if (ii2 == 2) {
                  table[-i0][-i0 + 1] = max_score(table[-i0][-i0 + 1], table[-i0 + 1][-i0]);
                  if (_PB_N + i0 >= 3) {
                    table[-i0][-i0 + 2] = max_score(table[-i0][-i0 + 2], table[-i0][-i0 + 1]);
                  }
                }
                if (_PB_N + i0 >= ii2 + 1) {
                  table[-i0][ii2 - i0] = max_score(table[-i0][ii2 - i0], table[-i0 + 1][ii2 - i0]);
                  if (ii2 == 2) {
                    table[-i0][-i0 + 2] = max_score(table[-i0][-i0 + 2], table[-i0 + 1][-i0 + 1] + match(seq[-i0], seq[-i0 + 2]));
                  }
                }
              }
            } else {
              for (int i0 = -_PB_N + 2; i0 <= 0; i0 += 1) {
                table[-i0][-i0 + 1] = max_score(table[-i0][-i0 + 1], table[-i0][-i0]);
              }
            }
          }
        }
        for (int ii2 = 0; ii2 <= min(2, -_PB_N + 128 * ii0 + 1); ii2 += 1) {
          if (ii2 == 2) {
            table[_PB_N - 64 * ii0 - 2][64 * ii0 + 1] = max_score(table[_PB_N - 64 * ii0 - 2][64 * ii0 + 1], table[_PB_N - 64 * ii0 - 1][64 * ii0] + match(seq[_PB_N - 64 * ii0 - 2], seq[64 * ii0 + 1]));
          } else if (ii2 == 1) {
            table[_PB_N - 64 * ii0 - 2][64 * ii0 + 1] = max_score(table[_PB_N - 64 * ii0 - 2][64 * ii0 + 1], table[_PB_N - 64 * ii0 - 1][64 * ii0 + 1]);
          } else {
            for (int i0 = -_PB_N + 64 * ii0 + 2; i0 <= 0; i0 += 1) {
              table[-i0][64 * ii0 + 1] = max_score(table[-i0][64 * ii0 + 1], table[-i0][64 * ii0]);
            }
          }
        }
        if (_PB_N == 128 && ii0 == 1) {
          table[62][65] = max_score(table[62][65], table[63][64] + match(seq[62], seq[65]));
        }
        for (int ii3 = 0; ii3 < ii0; ii3 += 1) {
          for (int i3 = max(_PB_N - 64 * ii0 - 1, 64 * ii3); i3 <= 64 * ii3 + 63; i3 += 1) {
            table[_PB_N - 64 * ii0 - 2][64 * ii0 + 1] = max_score(table[_PB_N - 64 * ii0 - 2][64 * ii0 + 1], table[_PB_N - 64 * ii0 - 2][i3] + table[i3 + 1][64 * ii0 + 1]);
          }
        }
        for (int i0 = max(-_PB_N + 3, -_PB_N + 64 * ii0 + 2); i0 <= 0; i0 += 1) {
          for (int i1 = max(64 * ii0 + 1, -i0 + 2); i1 < _PB_N; i1 += 1) {
            if (i0 + i1 >= 3 && i1 >= 64 * ii0 + 2) {
              table[-i0][i1] = max_score(table[-i0][i1], table[-i0][i1 - 1]);
            }
            if (_PB_N + i0 + i1 >= 128 * ii0 + 4 && i0 + i1 >= 3) {
              table[-i0][i1] = max_score(table[-i0][i1], table[-i0 + 1][i1]);
              table[-i0][i1] = max_score(table[-i0][i1], table[-i0 + 1][i1 - 1] + match(seq[-i0], seq[i1]));
            }
            if (_PB_N + i0 + i1 >= 128 * ii0 + 4) {
              for (int i3 = -i0 + 1; i3 < 64 * ii0; i3 += 1) {
                table[-i0][i1] = max_score(table[-i0][i1], table[-i0][i3] + table[i3 + 1][i1]);
              }
            }
            for (int i3 = max(64 * ii0, -i0 + 1); i3 < i1; i3 += 1) {
              table[-i0][i1] = max_score(table[-i0][i1], table[-i0][i3] + table[i3 + 1][i1]);
            }
          }
        }
      }
    } else {
      table[0][1] = max_score(table[0][1], table[0][0]);
      table[0][1] = max_score(table[0][1], table[1][1]);
      table[0][1] = max_score(table[0][1], table[1][0]);
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
  POLYBENCH_1D_ARRAY_DECL(seq, base, N, n);
  POLYBENCH_2D_ARRAY_DECL(table, DATA_TYPE, N, N, n, n);

  /* Initialize array(s). */
  init_array (n, POLYBENCH_ARRAY(seq), POLYBENCH_ARRAY(table));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_nussinov (n, POLYBENCH_ARRAY(seq), POLYBENCH_ARRAY(table));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(table)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(seq);
  POLYBENCH_FREE_ARRAY(table);

  return 0;
}
