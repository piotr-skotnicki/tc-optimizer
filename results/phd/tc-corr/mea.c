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
