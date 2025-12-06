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
#include "pbp.h"

static
int paired(char a, char b)
{
  return ((a == 'A' && b == 'U') || (a == 'U' && b == 'A') ||
          (a == 'G' && b == 'C') || (a == 'C' && b == 'G') ||
          (a == 'G' && b == 'U') || (a == 'U' && b == 'G'));
}

const int l        =  1;
const double Ebp   = -1.0;
const double RT    =  1.0;

/* Array initialization. */
static
void init_array (int n,
        int l,
        DATA_TYPE ERT,
        char POLYBENCH_1D(RNA,N,n),
        DATA_TYPE POLYBENCH_2D(Q,N,N,n,n),
        DATA_TYPE POLYBENCH_2D(Qbp,N,N,n,n),
        DATA_TYPE POLYBENCH_2D(Pbp,N,N,n,n))
{
  int i, j, k;

  srand(42);

  for (i = 0; i < _PB_N; i++) {
    RNA[i] = "ACGU"[rand() % 4];
  }

  for (i = 0; i < _PB_N + POLYBENCH_PADDING_FACTOR; i++) {
    for (j = 0; j < _PB_N + POLYBENCH_PADDING_FACTOR; j++) {
      Q[i][j] = SCALAR_VAL(1.0);
      Qbp[i][j] = SCALAR_VAL(0.0);
      Pbp[i][j] = SCALAR_VAL(0.0);
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
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
        DATA_TYPE POLYBENCH_2D(Pbp,N,N,n,n))
{
  int i, j;
  int t = 0;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("Pbp");
  for (i = 1; i <= _PB_N; i++) {
    for (j = 0; j <= _PB_N; j++) {
      if (t % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
      fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, Pbp[i][j]);
      t++;
    }
  }
  POLYBENCH_DUMP_END("Pbp");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_pbp(int n,
        DATA_TYPE ERT,
        char POLYBENCH_1D(RNA,N,n),
        DATA_TYPE POLYBENCH_2D(Q,N,N,n,n),
        DATA_TYPE POLYBENCH_2D(Qbp,N,N,n,n),
        DATA_TYPE POLYBENCH_2D(Pbp,N,N,n,n))
{
  int i, j, p, q;

/* TC Optimizing Compiler 0.5.1 */
/* ./tc ../examples/rnatools/pbp-fission-1.scop.c --rectangular-tiling --isl-wave-scheduling --omp-for-codegen --iterative-tc --align --debug -b 16 */
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

  /* Initialize array(s). */
  init_array (n, l, ERT, POLYBENCH_ARRAY(RNA), POLYBENCH_ARRAY(Q), POLYBENCH_ARRAY(Qbp), POLYBENCH_ARRAY(Pbp));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_pbp (n, ERT, POLYBENCH_ARRAY(RNA), POLYBENCH_ARRAY(Q), POLYBENCH_ARRAY(Qbp), POLYBENCH_ARRAY(Pbp));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(Pbp)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(RNA);
  POLYBENCH_FREE_ARRAY(Q);
  POLYBENCH_FREE_ARRAY(Qbp);
  POLYBENCH_FREE_ARRAY(Pbp);

  return 0;
}
