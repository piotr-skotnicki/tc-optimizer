/**
 * This version is stamped on Apr. 14, 2015
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* doitgen.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "doitgen.h"


/* Array initialization. */
static
void init_array(int nr, int nq, int np,
		DATA_TYPE POLYBENCH_3D(A,NR,NQ,NP,nr,nq,np),
		DATA_TYPE POLYBENCH_2D(C4,NP,NP,np,np))
{
  int i, j, k;

  for (i = 0; i < nr; i++)
    for (j = 0; j < nq; j++)
      for (k = 0; k < np; k++)
	A[i][j][k] = (DATA_TYPE) ((i*j + k)%np) / np;
  for (i = 0; i < np; i++)
    for (j = 0; j < np; j++)
      C4[i][j] = (DATA_TYPE) (i*j % np) / np;
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int nr, int nq, int np,
		 DATA_TYPE POLYBENCH_3D(A,NR,NQ,NP,nr,nq,np))
{
  int i, j, k;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("A");
  for (i = 0; i < nr; i++)
    for (j = 0; j < nq; j++)
      for (k = 0; k < np; k++) {
	if ((i*nq*np+j*np+k) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
	fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, A[i][j][k]);
      }
  POLYBENCH_DUMP_END("A");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
void kernel_doitgen(int nr, int nq, int np,
		    DATA_TYPE POLYBENCH_3D(A,NR,NQ,NP,nr,nq,np),
		    DATA_TYPE POLYBENCH_2D(C4,NP,NP,np,np),
		    DATA_TYPE POLYBENCH_1D(sum,NP,np))
{
  int r, q, p, s;


/* TC Optimizing Compiler 0.2.26 */
/* ./tc ../examples/polybench/doitgen.scop.c --correction-tiling --lex-scheduling --serial-codegen -b 64 --debug --isl-union-map-tc */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#pragma scop
for (int ii0 = 0; ii0 <= floord(_PB_NR - 1, 64); ii0 += 1) {
  for (int ii1 = 0; ii1 <= floord(_PB_NQ - 1, 64); ii1 += 1) {
    for (int ii3 = 0; ii3 <= floord(_PB_NP - 1, 64); ii3 += 1) {
      for (int i3 = 64 * ii3; i3 <= min(_PB_NP - 1, 64 * ii3 + 63); i3 += 1) {
        sum[i3] = SCALAR_VAL(0.0);
      }
      for (int ii5 = 0; ii5 <= (_PB_NP - 1) / 64; ii5 += 1) {
        if (_PB_NP >= 64 * ii3 + 64) {
          for (int i3 = 64 * ii3; i3 <= 64 * ii3 + 63; i3 += 1) {
            for (int i5 = 64 * ii5; i5 <= min(_PB_NP - 1, 64 * ii5 + 63); i5 += 1) {
              sum[i3] += (A[64 * ii0][64 * ii1][i5] * C4[i5][i3]);
            }
          }
        } else {
          for (int i3 = 64 * ii3; i3 < _PB_NP; i3 += 1) {
            for (int i5 = 64 * ii5; i5 <= min(_PB_NP - 1, 64 * ii5 + 63); i5 += 1) {
              sum[i3] += (A[64 * ii0][64 * ii1][i5] * C4[i5][i3]);
            }
          }
        }
      }
    }
    for (int ii3 = 0; ii3 <= floord(_PB_NP - 1, 64); ii3 += 1) {
      if (_PB_NP >= 64 * ii3 + 65) {
        for (int i3 = 64 * ii3; i3 <= 64 * ii3 + 63; i3 += 1) {
          A[64 * ii0][64 * ii1][i3] = sum[i3];
        }
        if (_PB_NQ >= 64 * ii1 + 2) {
          for (int i3 = 64 * ii3; i3 <= 64 * ii3 + 63; i3 += 1) {
            sum[i3] = SCALAR_VAL(0.0);
            for (int i5 = 0; i5 < _PB_NP; i5 += 1) {
              sum[i3] += (A[64 * ii0][64 * ii1 + 1][i5] * C4[i5][i3]);
            }
          }
        }
        if (_PB_NR >= 64 * ii0 + 2 && 64 * ii1 + 1 == _PB_NQ) {
          for (int i3 = 64 * ii3; i3 <= 64 * ii3 + 63; i3 += 1) {
            sum[i3] = SCALAR_VAL(0.0);
            for (int i5 = 0; i5 < _PB_NP; i5 += 1) {
              sum[i3] += (A[64 * ii0 + 1][0][i5] * C4[i5][i3]);
            }
          }
        }
      } else if (_PB_NQ >= 64 * ii1 + 65) {
        for (int i1 = 64 * ii1; i1 <= 64 * ii1 + 63; i1 += 1) {
          if (i1 >= 64 * ii1 + 1) {
            if (i1 >= 64 * ii1 + 2) {
              for (int i3 = 0; i3 < _PB_NP; i3 += 1) {
                if (64 * ii3 >= i3 + 1) {
                  sum[i3] = SCALAR_VAL(0.0);
                } else {
                  sum[i3] = SCALAR_VAL(0.0);
                }
                if (64 * ii3 >= i3 + 1) {
                  for (int i5 = 0; i5 < _PB_NP; i5 += 1) {
                    sum[i3] += (A[64 * ii0][i1][i5] * C4[i5][i3]);
                  }
                } else {
                  for (int i5 = 0; i5 < _PB_NP; i5 += 1) {
                    sum[i3] += (A[64 * ii0][i1][i5] * C4[i5][i3]);
                  }
                }
              }
            } else {
              for (int i3 = 64 * ii3; i3 < _PB_NP; i3 += 1) {
                sum[i3] = SCALAR_VAL(0.0);
                if (64 * ii3 + 63 >= _PB_NP) {
                  for (int i5 = 0; i5 < _PB_NP; i5 += 1) {
                    sum[i3] += (A[64 * ii0][64 * ii1 + 1][i5] * C4[i5][i3]);
                  }
                } else {
                  for (int i5 = 0; i5 < _PB_NP; i5 += 1) {
                    sum[i3] += (A[64 * ii0][64 * ii1 + 1][i5] * C4[i5][i3]);
                  }
                }
              }
            }
          }
          if (i1 >= 64 * ii1 + 1) {
            for (int i3 = 0; i3 < 64 * ii3; i3 += 1) {
              A[64 * ii0][i1][i3] = sum[i3];
            }
          }
          for (int i3 = 64 * ii3; i3 < _PB_NP; i3 += 1) {
            A[64 * ii0][i1][i3] = sum[i3];
          }
        }
      } else {
        for (int i0 = 64 * ii0; i0 <= min(_PB_NR - 1, 64 * ii0 + 63); i0 += 1) {
          if (i0 >= 64 * ii0 + 1) {
            if (_PB_NQ >= 65 && 64 * ii1 + 1 == _PB_NQ && i0 == 64 * ii0 + 1) {
              if (64 * ii3 + 63 >= _PB_NP) {
                for (int i3 = 64 * ii3; i3 < _PB_NP; i3 += 1) {
                  sum[i3] = SCALAR_VAL(0.0);
                  for (int i5 = 0; i5 < _PB_NP; i5 += 1) {
                    sum[i3] += (A[64 * ii0 + 1][0][i5] * C4[i5][i3]);
                  }
                }
              } else {
                for (int i3 = _PB_NP - 64; i3 < _PB_NP; i3 += 1) {
                  sum[i3] = SCALAR_VAL(0.0);
                  for (int i5 = 0; i5 < _PB_NP; i5 += 1) {
                    sum[i3] += (A[64 * ii0 + 1][0][i5] * C4[i5][i3]);
                  }
                }
              }
              for (int i3 = 0; i3 < _PB_NP; i3 += 1) {
                A[64 * ii0 + 1][0][i3] = sum[i3];
              }
            }
            for (int i1 = max(0, -_PB_NQ + 64 * ii0 + 64 * ii1 - i0 + 3); i1 < _PB_NQ; i1 += 1) {
              for (int i3 = 0; i3 < _PB_NP; i3 += 1) {
                if (64 * ii1 >= i1 + 1) {
                  sum[i3] = SCALAR_VAL(0.0);
                } else if (64 * ii3 >= i3 + 1) {
                  sum[i3] = SCALAR_VAL(0.0);
                } else {
                  sum[i3] = SCALAR_VAL(0.0);
                }
                if (64 * ii3 + 63 >= _PB_NP && i3 >= 64 * ii3) {
                  for (int i5 = 0; i5 < _PB_NP; i5 += 1) {
                    sum[i3] += (A[i0][i1][i5] * C4[i5][i3]);
                  }
                } else if ((i3 % 64) + _PB_NP >= i3 + 64) {
                  for (int i5 = 0; i5 < _PB_NP; i5 += 1) {
                    sum[i3] += (A[i0][i1][i5] * C4[i5][i3]);
                  }
                }
              }
              for (int i3 = 0; i3 < _PB_NP; i3 += 1) {
                A[i0][i1][i3] = sum[i3];
              }
            }
            if (_PB_NQ == 1 && ii1 == 0 && i0 == 64 * ii0 + 1) {
              for (int i3 = 64 * ii3; i3 < _PB_NP; i3 += 1) {
                sum[i3] = SCALAR_VAL(0.0);
                for (int i5 = 0; i5 < _PB_NP; i5 += 1) {
                  sum[i3] += (A[64 * ii0 + 1][0][i5] * C4[i5][i3]);
                }
              }
              for (int i3 = 0; i3 < 64 * ii3; i3 += 1) {
                A[64 * ii0 + 1][0][i3] = sum[i3];
              }
              for (int i3 = 64 * ii3; i3 < _PB_NP; i3 += 1) {
                A[64 * ii0 + 1][0][i3] = sum[i3];
              }
            }
          } else {
            for (int i1 = 64 * ii1; i1 < _PB_NQ; i1 += 1) {
              if (i1 >= 64 * ii1 + 1) {
                if (i1 >= 64 * ii1 + 2) {
                  for (int i3 = 0; i3 < 64 * ii3; i3 += 1) {
                    sum[i3] = SCALAR_VAL(0.0);
                    for (int i5 = 0; i5 < _PB_NP; i5 += 1) {
                      sum[i3] += (A[64 * ii0][i1][i5] * C4[i5][i3]);
                    }
                  }
                }
                for (int i3 = 64 * ii3; i3 < _PB_NP; i3 += 1) {
                  sum[i3] = SCALAR_VAL(0.0);
                  for (int i5 = 0; i5 < _PB_NP; i5 += 1) {
                    sum[i3] += (A[64 * ii0][i1][i5] * C4[i5][i3]);
                  }
                }
              }
              if (i1 >= 64 * ii1 + 1) {
                for (int i3 = 0; i3 < 64 * ii3; i3 += 1) {
                  A[64 * ii0][i1][i3] = sum[i3];
                }
              }
              for (int i3 = 64 * ii3; i3 < _PB_NP; i3 += 1) {
                A[64 * ii0][i1][i3] = sum[i3];
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
  int nr = NR;
  int nq = NQ;
  int np = NP;

  /* Variable declaration/allocation. */
  POLYBENCH_3D_ARRAY_DECL(A,DATA_TYPE,NR,NQ,NP,nr,nq,np);
  POLYBENCH_1D_ARRAY_DECL(sum,DATA_TYPE,NP,np);
  POLYBENCH_2D_ARRAY_DECL(C4,DATA_TYPE,NP,NP,np,np);

  /* Initialize array(s). */
  init_array (nr, nq, np,
	      POLYBENCH_ARRAY(A),
	      POLYBENCH_ARRAY(C4));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_doitgen (nr, nq, np,
		  POLYBENCH_ARRAY(A),
		  POLYBENCH_ARRAY(C4),
		  POLYBENCH_ARRAY(sum));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(nr, nq, np,  POLYBENCH_ARRAY(A)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(sum);
  POLYBENCH_FREE_ARRAY(C4);

  return 0;
}
