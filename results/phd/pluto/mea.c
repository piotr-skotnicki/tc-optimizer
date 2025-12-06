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
#include "mea.h"

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

  for (i = 0; i < _PB_N + POLYBENCH_PADDING_FACTOR; i++) {
    for (j = 0; j < _PB_N + POLYBENCH_PADDING_FACTOR; j++) {
      Q[i][j] = SCALAR_VAL(1.0);
      Qbp[i][j] = SCALAR_VAL(0.0);
      Pbp[i][j] = SCALAR_VAL(0.0);
      Pu[i][j] = SCALAR_VAL(0.0);
      M[i][j] = SCALAR_VAL(0.0);
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

  for (i = 1; i <= _PB_N; i++) {
    for (j = 1; j <= _PB_N; j++) {
      Pu[i][j] = (Q[1][i-1] * Q[j+1][_PB_N]) / Q[1][_PB_N];
      for (p = 1; p < i; p++) {
        for (q = j+1; q <= _PB_N; q++) {
          Pu[i][j] += paired(RNA[p-1], RNA[q-1]) * ((Pbp[p][q] * ERT * Q[p+1][i-1] * Q[j+1][q-1]) / (Qbp[p][q] == 0.0 ? 1.0 : Qbp[p][q]));
        }
      }
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
void kernel_mea(int n,
        int l,
        DATA_TYPE gamma,
        char POLYBENCH_1D(RNA,N,n),
        DATA_TYPE POLYBENCH_2D(Pbp,N,N,n,n),
        DATA_TYPE POLYBENCH_2D(Pu,N,N,n,n),
        DATA_TYPE POLYBENCH_2D(M,N,N,n,n))
{
  int i, j, k;

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
  kernel_mea (n, l, gamma, POLYBENCH_ARRAY(RNA), POLYBENCH_ARRAY(Pbp), POLYBENCH_ARRAY(Pu), POLYBENCH_ARRAY(M));

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
