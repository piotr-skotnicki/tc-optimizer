/**
 * This version is stamped on Oct 25, 2025
 */

#include <omp.h>
#include <math.h>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

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

{
  int t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
  int lb, ub, lbp, ubp, lb2, ub2;
  register int lbv, ubv;
  /* Start of CLooG code */
  if ((_PB_N >= 2) && (l >= 0) && (l <= 5)) {
    for (t2=-1;t2<=floord(_PB_N-16,16);t2++) {
      lbp=t2+1;
      ubp=min(floord(_PB_N,16),floord(16*t2+_PB_N+14,16));
  #pragma omp parallel for private(lbv,ubv,t5,t6,t7,t8,t9,t10)
      for (t4=lbp;t4<=ubp;t4++) {
        if (t2 == -1) {
          for (t5=max(-_PB_N+1,-16*t4-14);t5<=-16*t4+l-15;t5++) {
            for (t7=-t5+1;t7<=min(_PB_N,16*t4+15);t7++) {
              Q[-t5][t7] = Q[-t5][t7-1];;
            }
          }
        }
        for (t5=max(max(-_PB_N+1,16*t2-16*t4),-16*t4+l-14);t5<=min(-_PB_N+l,16*t2-16*t4+15);t5++) {
          for (t7=max(16*t4,-t5+1);t7<=_PB_N;t7++) {
            Q[-t5][t7] = Q[-t5][t7-1];;
          }
        }
        if (l >= 1) {
          for (t5=max(max(16*t2-16*t4,-_PB_N+l+1),-16*t4+l-14);t5<=min(-16*t4+l,16*t2-16*t4+15);t5++) {
            for (t7=max(16*t4,-t5+1);t7<=-t5+l;t7++) {
              Q[-t5][t7] = Q[-t5][t7-1];;
            }
            for (t7=-t5+l+1;t7<=min(_PB_N,16*t4+15);t7++) {
              Q[-t5][t7] = Q[-t5][t7-1];;
              for (t9=-t5;t9<=t7-l-1;t9++) {
                Qbp[-t5][t7] = paired(RNA[-t5-1], RNA[t7-1]) * Q[-t5+1][t7-1] * ERT;;
                Q[-t5][t7] += Q[-t5][t9-1] * Qbp[t9][t7];;
              }
            }
          }
        }
        if (l >= 1) {
          for (t5=max(16*t2-16*t4,-16*t4+l+1);t5<=16*t2-16*t4+15;t5++) {
            for (t7=16*t4;t7<=min(_PB_N,16*t4+15);t7++) {
              Q[-t5][t7] = Q[-t5][t7-1];;
              for (t9=-t5;t9<=t7-l-1;t9++) {
                Qbp[-t5][t7] = paired(RNA[-t5-1], RNA[t7-1]) * Q[-t5+1][t7-1] * ERT;;
                Q[-t5][t7] += Q[-t5][t9-1] * Qbp[t9][t7];;
              }
            }
          }
        }
        if (l == 0) {
          for (t5=max(max(-_PB_N+1,16*t2-16*t4),-16*t4-14);t5<=16*t2-16*t4+15;t5++) {
            for (t7=max(16*t4,-t5+1);t7<=min(_PB_N,16*t4+15);t7++) {
              Q[-t5][t7] = Q[-t5][t7-1];;
              for (t9=-t5;t9<=t7-1;t9++) {
                Qbp[-t5][t7] = paired(RNA[-t5-1], RNA[t7-1]) * Q[-t5+1][t7-1] * ERT;;
                Q[-t5][t7] += Q[-t5][t9-1] * Qbp[t9][t7];;
              }
            }
          }
        }
      }
    }
  }
}

{
  int t1, t2, t3, t4, t5, t6, t7;
  int lb, ub, lbp, ubp, lb2, ub2;
  register int lbv, ubv;
  if (_PB_N >= 1) {
    lbp=0;
    ubp=floord(_PB_N,16);
  #pragma omp parallel for private(lbv,ubv,t3,t4,t5,t6,t7)
    for (t2=lbp;t2<=ubp;t2++) {
      for (t3=0;t3<=floord(_PB_N,16);t3++) {
        for (t4=max(1,16*t2);t4<=(min(_PB_N,16*t2+15))-7;t4+=8) {
          lbv=max(1,16*t3);
          ubv=min(_PB_N,16*t3+15);
  #pragma ivdep
  #pragma vector always
          for (t5=lbv;t5<=ubv;t5++) {
            Pbp[t4][t5] = (Q[1][t4-1] * Qbp[t4][t5] * Q[t5+1][_PB_N]) / Q[1][_PB_N];;
            Pbp[(t4+1)][t5] = (Q[1][(t4+1)-1] * Qbp[(t4+1)][t5] * Q[t5+1][_PB_N]) / Q[1][_PB_N];;
            Pbp[(t4+2)][t5] = (Q[1][(t4+2)-1] * Qbp[(t4+2)][t5] * Q[t5+1][_PB_N]) / Q[1][_PB_N];;
            Pbp[(t4+3)][t5] = (Q[1][(t4+3)-1] * Qbp[(t4+3)][t5] * Q[t5+1][_PB_N]) / Q[1][_PB_N];;
            Pbp[(t4+4)][t5] = (Q[1][(t4+4)-1] * Qbp[(t4+4)][t5] * Q[t5+1][_PB_N]) / Q[1][_PB_N];;
            Pbp[(t4+5)][t5] = (Q[1][(t4+5)-1] * Qbp[(t4+5)][t5] * Q[t5+1][_PB_N]) / Q[1][_PB_N];;
            Pbp[(t4+6)][t5] = (Q[1][(t4+6)-1] * Qbp[(t4+6)][t5] * Q[t5+1][_PB_N]) / Q[1][_PB_N];;
            Pbp[(t4+7)][t5] = (Q[1][(t4+7)-1] * Qbp[(t4+7)][t5] * Q[t5+1][_PB_N]) / Q[1][_PB_N];;
          }
        }
        for (;t4<=min(_PB_N,16*t2+15);t4++) {
          lbv=max(1,16*t3);
          ubv=min(_PB_N,16*t3+15);
  #pragma ivdep
  #pragma vector always
          for (t5=lbv;t5<=ubv;t5++) {
            Pbp[t4][t5] = (Q[1][t4-1] * Qbp[t4][t5] * Q[t5+1][_PB_N]) / Q[1][_PB_N];;
          }
        }
      }
    }
    if (_PB_N >= 2) {
      for (t2=0;t2<=floord(_PB_N-1,8);t2++) {
        lbp=max(0,ceild(16*t2-_PB_N,16));
        ubp=floord(t2,2);
  #pragma omp parallel for private(lbv,ubv,t4,t5,t6,t7)
        for (t3=lbp;t3<=ubp;t3++) {
          for (t4=max(max(2,16*t2-16*t3),16*t3+1);t4<=min(_PB_N,16*t2-16*t3+15);t4++) {
            for (t5=max(1,16*t3);t5<=min(16*t3+15,t4-1);t5++) {
              for (t6=2;t6<=_PB_N;t6++) {
                lbv=1;
                ubv=t6-1;
  #pragma ivdep
  #pragma vector always
                for (t7=lbv;t7<=ubv;t7++) {
                  Pbp[t4][t7] += paired(RNA[t5-1], RNA[t6-1]) * ((Pbp[t5][t6] * ERT * Q[t5+1][t4-1] * Qbp[t4][t7] * Q[t7+1][t6-1]) / (Qbp[t5][t6] == 0.0 ? 1.0 : Qbp[t5][t6]));;
                }
              }
            }
          }
        }
      }
    }
  }
}

{
  int t1, t2, t3, t4, t5, t6, t7, t8;
  int lb, ub, lbp, ubp, lb2, ub2;
  register int lbv, ubv;
  if (_PB_N >= 1) {
    lbp=0;
    ubp=floord(_PB_N,16);
  #pragma omp parallel for private(lbv,ubv,t3,t4,t5,t6,t7,t8)
    for (t2=lbp;t2<=ubp;t2++) {
      for (t3=0;t3<=floord(_PB_N,16);t3++) {
        for (t4=max(1,16*t2);t4<=(min(_PB_N,16*t2+15))-7;t4+=8) {
          lbv=max(1,16*t3);
          ubv=min(_PB_N,16*t3+15);
  #pragma ivdep
  #pragma vector always
          for (t5=lbv;t5<=ubv;t5++) {
            Pu[t4][t5] = (Q[1][t4-1] * Q[t5+1][_PB_N]) / Q[1][_PB_N];;
            Pu[(t4+1)][t5] = (Q[1][(t4+1)-1] * Q[t5+1][_PB_N]) / Q[1][_PB_N];;
            Pu[(t4+2)][t5] = (Q[1][(t4+2)-1] * Q[t5+1][_PB_N]) / Q[1][_PB_N];;
            Pu[(t4+3)][t5] = (Q[1][(t4+3)-1] * Q[t5+1][_PB_N]) / Q[1][_PB_N];;
            Pu[(t4+4)][t5] = (Q[1][(t4+4)-1] * Q[t5+1][_PB_N]) / Q[1][_PB_N];;
            Pu[(t4+5)][t5] = (Q[1][(t4+5)-1] * Q[t5+1][_PB_N]) / Q[1][_PB_N];;
            Pu[(t4+6)][t5] = (Q[1][(t4+6)-1] * Q[t5+1][_PB_N]) / Q[1][_PB_N];;
            Pu[(t4+7)][t5] = (Q[1][(t4+7)-1] * Q[t5+1][_PB_N]) / Q[1][_PB_N];;
          }
        }
        for (;t4<=min(_PB_N,16*t2+15);t4++) {
          lbv=max(1,16*t3);
          ubv=min(_PB_N,16*t3+15);
  #pragma ivdep
  #pragma vector always
          for (t5=lbv;t5<=ubv;t5++) {
            Pu[t4][t5] = (Q[1][t4-1] * Q[t5+1][_PB_N]) / Q[1][_PB_N];;
          }
        }
      }
    }
    if (_PB_N >= 2) {
      lbp=0;
      ubp=floord(_PB_N,16);
  #pragma omp parallel for private(lbv,ubv,t3,t4,t5,t6,t7,t8)
      for (t2=lbp;t2<=ubp;t2++) {
        for (t3=0;t3<=floord(_PB_N-1,16);t3++) {
          for (t4=0;t4<=min(floord(_PB_N-1,16),t2);t4++) {
            for (t5=max(max(2,16*t2),16*t4+1);t5<=min(_PB_N,16*t2+15);t5++) {
              for (t6=max(1,16*t4);t6<=min(16*t4+15,t5-1);t6++) {
                for (t7=max(2,16*t3+1);t7<=_PB_N;t7++) {
                  lbv=max(1,16*t3);
                  ubv=min(16*t3+15,t7-1);
  #pragma ivdep
  #pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    Pu[t5][t8] += paired(RNA[t6-1], RNA[t7-1]) * ((Pbp[t6][t7] * ERT * Q[t6+1][t5-1] * Q[t8+1][t7-1]) / (Qbp[t6][t7] == 0.0 ? 1.0 : Qbp[t6][t7]));;
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

{
  int t1, t2, t3, t4, t5;
  int lb, ub, lbp, ubp, lb2, ub2;
  register int lbv, ubv;
  if ((_PB_N >= 1) && (l >= 0) && (l <= 5)) {
    for (t1=1;t1<=_PB_N;t1++) {
      lbp=0;
      ubp=floord(t1,16);
  #pragma omp parallel for private(lbv,ubv,t3,t4,t5)
      for (t2=lbp;t2<=ubp;t2++) {
        for (t3=0;t3<=floord(t1-1,16);t3++) {
          if ((t2 == 0) && (t3 == 0)) {
            for (t4=1;t4<=min(15,t1-l-1);t4++) {
              M[t4][t1] = M[t4][t1-1] + Pu[t1][t1];;
              for (t5=t4;t5<=min(15,t1-l-1);t5++) {
                M[t4][t1] = max_score(M[t4][t1], paired(RNA[t5-1], RNA[t1-1]) * (M[t4][t5-1] + M[t5+1][t1-1] + gamma * Pbp[t5][t1]));;
              }
            }
          }
          if (t3 == 0) {
            for (t4=max(16,16*t2);t4<=min(16*t2+15,t1-l-1);t4++) {
              M[t4][t1] = M[t4][t1-1] + Pu[t1][t1];;
            }
          }
          if (t3 == 0) {
            for (t4=max(max(1,16*t2),t1-l);t4<=min(t1,16*t2+15);t4++) {
              M[t4][t1] = M[t4][t1-1] + Pu[t1][t1];;
            }
          }
          if ((t1 >= 16*t3+l+1) && (t3 >= 1)) {
            for (t4=max(1,16*t2);t4<=min(min(16*t2+15,16*t3+15),t1-l-1);t4++) {
              for (t5=max(16*t3,t4);t5<=min(16*t3+15,t1-l-1);t5++) {
                M[t4][t1] = max_score(M[t4][t1], paired(RNA[t5-1], RNA[t1-1]) * (M[t4][t5-1] + M[t5+1][t1-1] + gamma * Pbp[t5][t1]));;
              }
            }
          }
        }
      }
    }
  }
}

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
