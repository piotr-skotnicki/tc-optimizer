/**
 * This version is stamped on Oct 25, 2025
 */

#include <omp.h>
#include <math.h>
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
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
#include "pu.h"

static
int paired(char a, char b)
{
  return ((a == 'A' && b == 'U') || (a == 'U' && b == 'A') ||
          (a == 'G' && b == 'C') || (a == 'C' && b == 'G') ||
          (a == 'G' && b == 'U') || (a == 'U' && b == 'G'));
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
        DATA_TYPE POLYBENCH_2D(Pu,N,N,n,n))
{
  int i, j, k, p, q;

  srand(42);

  for (i = 0; i < _PB_N; i++) {
    RNA[i] = "ACGU"[rand() % 4];
  }

  for (i = 0; i < _PB_N + POLYBENCH_PADDING_FACTOR; i++) {
    for (j = 0; j < _PB_N + POLYBENCH_PADDING_FACTOR; j++) {
      Q[i][j] = SCALAR_VAL(1.0);
      Qbp[i][j] = SCALAR_VAL(0.0);
      Pbp[i][j] = SCALAR_VAL(0.0);
      Pu[i][j] = SCALAR_VAL(0.0);
    }
  }

  for (i = _PB_N; i >= 1; i--) {
    for (j = i+1; j <= _PB_N; j++) {
      Q[i][j] = Q[i][j-1];
      for (k = i; k < j-l; k++) {
        Qbp[i][j] = paired(RNA[i-1], RNA[j-1]) * Q[i+1][j-1] * ERT;
        Q[i][j] += Q[i][k-1] * Qbp[k][j];
      }
    }
  }

  for (i = 1; i <= _PB_N; i++) {
    for (j = 1; j <= _PB_N; j++) {
      Pbp[i][j] = (Q[1][i-1] * Qbp[i][j] * Q[j+1][_PB_N]) / Q[1][_PB_N];
      for (p = 1; p < i; p++) {
        for (q = j+1; q <= _PB_N; q++) {
          Pbp[i][j] += paired(RNA[p-1], RNA[q-1]) * ((Pbp[p][q] * ERT * Q[p+1][i-1] * Qbp[i][j] * Q[j+1][q-1]) / (Qbp[p][q] == 0.0 ? 1.0 : Qbp[p][q]));
        }
      }
    }
  }
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
        DATA_TYPE POLYBENCH_2D(Pu,N,N,n,n))
{
  int i;
  int t = 0;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("Pu");
  for (i = 1; i <= _PB_N; i++) {
    if (t % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
    fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, Pu[i][i]);
    t++;
  }
  POLYBENCH_DUMP_END("Pu");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_pu(int n,
        DATA_TYPE ERT,
        char POLYBENCH_1D(RNA,N,n),
        DATA_TYPE POLYBENCH_2D(Q,N,N,n,n),
        DATA_TYPE POLYBENCH_2D(Qbp,N,N,n,n),
        DATA_TYPE POLYBENCH_2D(Pbp,N,N,n,n),
        DATA_TYPE POLYBENCH_2D(Pu,N,N,n,n))
{
  int i, j, p, q;

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


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int n = N;

  int l         = 1;
  DATA_TYPE Ebp = SCALAR_VAL(-1.0);
  DATA_TYPE RT  = SCALAR_VAL(1.0);
  DATA_TYPE ERT = EXP_FUN(-Ebp/RT);

  /* Variable declaration/allocation. */
  POLYBENCH_1D_ARRAY_DECL(RNA, char, N, n);
  POLYBENCH_2D_ARRAY_DECL(Q, DATA_TYPE, N, N, n, n);
  POLYBENCH_2D_ARRAY_DECL(Qbp, DATA_TYPE, N, N, n, n);
  POLYBENCH_2D_ARRAY_DECL(Pbp, DATA_TYPE, N, N, n, n);
  POLYBENCH_2D_ARRAY_DECL(Pu, DATA_TYPE, N, N, n, n);

  /* Initialize array(s). */
  init_array (n, l, ERT, POLYBENCH_ARRAY(RNA), POLYBENCH_ARRAY(Q), POLYBENCH_ARRAY(Qbp), POLYBENCH_ARRAY(Pbp), POLYBENCH_ARRAY(Pu));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_pu (n, ERT, POLYBENCH_ARRAY(RNA), POLYBENCH_ARRAY(Q), POLYBENCH_ARRAY(Qbp), POLYBENCH_ARRAY(Pbp), POLYBENCH_ARRAY(Pu));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(Pu)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(RNA);
  POLYBENCH_FREE_ARRAY(Q);
  POLYBENCH_FREE_ARRAY(Qbp);
  POLYBENCH_FREE_ARRAY(Pbp);
  POLYBENCH_FREE_ARRAY(Pu);

  return 0;
}
