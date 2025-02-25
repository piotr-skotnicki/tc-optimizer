/*
 * Discretized 2D heat equation stencil with non periodic boundary conditions
 * Adapted from Pochoir test bench
 *
 * Irshad Pananilath: irshad@csa.iisc.ernet.in
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

/*
 * N is the number of points
 * T is the number of timesteps
 */
#ifdef HAS_DECLS
#include "decls.h"
#else
#define N 4000L
#define T 1000L
#endif

#define NUM_FP_OPS 10

/* Define our arrays */
double A[2][N + 2][N + 2];
double total = 0;
double sum_err_sqr = 0;
int chtotal = 0;
int timeval_subtract(struct timeval *result, struct timeval *x,
                     struct timeval *y) {
  if (x->tv_usec < y->tv_usec) {
    int nsec = (y->tv_usec - x->tv_usec) / 1000000 + 1;

    y->tv_usec -= 1000000 * nsec;
    y->tv_sec += nsec;
  }

  if (x->tv_usec - y->tv_usec > 1000000) {
    int nsec = (x->tv_usec - y->tv_usec) / 1000000;

    y->tv_usec += 1000000 * nsec;
    y->tv_sec -= nsec;
  }

  result->tv_sec = x->tv_sec - y->tv_sec;
  result->tv_usec = x->tv_usec - y->tv_usec;

  return x->tv_sec < y->tv_sec;
}

int main(int argc, char *argv[]) {
  long int i, j;
  const int BASE = 1024;

  // for timekeeping
  int ts_return = -1;
  struct timeval start, end, result;
  double tdiff = 0.0;

  printf("Number of points = %ld\t|Number of timesteps = %ld\t", N * N, T);

  /* Initialization */
  srand(42); // seed with a constant value to verify results

  for (i = 1; i < N + 1; i++) {
    for (j = 1; j < N + 1; j++) {
      A[0][i][j] = 1.0 * (rand() % BASE);
    }
  }

#ifdef TIME
  gettimeofday(&start, 0);
#endif

/* TC Optimizing Compiler 0.4.1 */
/* ./tc ../examples/pluto/heat-2d-t.scop.c --semi-diamond-tiling --omp-for-codegen --iterative-tc --debug -b 16 --drop-bounds --use-macros */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define S1_I(t,i,j) A[(t + 1)%2][i][j] = (((0.125 * ((A[t%2][i + 1][j] - (2.0 * A[t%2][i][j])) + A[t%2][i - 1][j])) + (0.125 * ((A[t%2][i][j + 1] - (2.0 * A[t%2][i][j])) + A[t%2][i][j - 1]))) + A[t%2][i][j]);
#define S1(t,i,j) S1_I((t),(i),(j))
#pragma scop
if (N >= 1) {
  const int ii0_60485_lb = 0, ii0_83148_ub = floord(T - 1, 16);
  for (register int ii0 = ii0_60485_lb; ii0 <= ii0_83148_ub; ii0 += 1) {
    const int ii1_92980_lb = 0, ii1_8075_ub = floord(N - 1, 16);
    #pragma omp parallel for
    for (register int ii1 = ii1_92980_lb; ii1 < ii1_8075_ub; ii1 += 1) {
      const int ii2_91851_lb = 0, ii2_56022_ub = min(-ii0 + (T + N - 1) / 16, (N + 6) / 16);
      for (register int ii2 = ii2_91851_lb; ii2 <= ii2_56022_ub; ii2 += 1) {
        const int i0_10899_lb = max(16 * ii0, -N + 16 * ii0 + 16 * ii2), i0_33307_ub = min(min(T - 1, 16 * ii0 + 6), N + 16 * ii0 - 16 * ii1 - 17);
        for (register int i0 = i0_10899_lb; i0 <= i0_33307_ub; i0 += 1) {
          const int i1_32766_lb = -16 * ii0 + 16 * ii1 + i0 + 17, i1_28974_ub = min(N, 16 * ii0 + 16 * ii1 - i0 + 30);
          for (register int i1 = i1_32766_lb; i1 <= i1_28974_ub; i1 += 1) {
            const int i2_66006_lb = max(1, 16 * ii0 + 16 * ii2 - i0), i2_7979_ub = min(N, 16 * ii0 + 16 * ii2 - i0 + 15);
            for (register int i2 = i2_66006_lb; i2 <= i2_7979_ub; i2 += 1) {
              S1(i0, i1, i2);
            }
          }
        }
      }
    }
    const int k_69126_lb = 1, k_90962_ub = min(2, -ii0 + (T - ii0 + 13) / 15 + 1);
    for (register int k = k_69126_lb; k <= k_90962_ub; k += 1) {
      if (k == 2) {
        const int ii2_57393_lb = 0, ii2_20855_ub = min((N + 14) / 16, -ii0 + floord(T + N - 2, 16));
        for (register int ii2 = ii2_57393_lb; ii2 <= ii2_20855_ub; ii2 += 1) {
          const int i0_16648_lb = max(16 * ii0 + 1, -N + 16 * ii0 + 16 * ii2 + 1), i0_22412_ub = min(T - 1, 16 * ii0 + 15);
          for (register int i0 = i0_16648_lb; i0 <= i0_22412_ub; i0 += 1) {
            const int i1_26538_lb = 1, i1_85039_ub = min(min(N, -16 * ii0 + i0), N - 16 * ii0 - 16 * ii2 + i0);
            for (register int i1 = i1_26538_lb; i1 <= i1_85039_ub; i1 += 1) {
              const int i2_61735_lb = max(1, 16 * ii0 + 16 * ii2 - i0 + i1), i2_60888_ub = min(N, 16 * ii0 + 16 * ii2 - i0 + i1 + 15);
              for (register int i2 = i2_61735_lb; i2 <= i2_60888_ub; i2 += 1) {
                S1(i0, i1, i2);
              }
            }
          }
        }
      } else {
        const int ii2_66945_lb = max(0, -((N + 14) / 16) + 1), ii2_36108_ub = min(min(floord(N - 1, 8), (N + 14) / 16), -ii0 + floord(T + N - 1, 16));
        for (register int ii2 = ii2_66945_lb; ii2 <= ii2_36108_ub; ii2 += 1) {
          const int i0_2242_lb = max(16 * ii0, -N + 16 * ii0 + 16 * ii2), i0_50782_ub = min(min(T - 1, N + 16 * ii0 - 1), 16 * ii0 + 14);
          for (register int i0 = i0_2242_lb; i0 <= i0_50782_ub; i0 += 1) {
            if (ii2 == 0) {
              const int i2_22608_lb = 1, i2_64506_ub = min(N, 16 * ii0 - i0 + 15);
              for (register int i2 = i2_22608_lb; i2 <= i2_64506_ub; i2 += 1) {
                S1(i0, -16 * ii0 + i0 + 1, i2);
              }
            }
            const int i1_22325_lb = max(-16 * ii0 + i0 + 1, -16 * ii0 - 16 * ii2 + i0 + 2), i1_74159_ub = min(min(min(N, -16 * ii0 + i0 + 16), N - 16 * ii0 - 16 * ii2 + i0 + 16), 16 * ii0 - i0 + 30);
            for (register int i1 = i1_22325_lb; i1 <= i1_74159_ub; i1 += 1) {
              const int i2_58192_lb = max(max(1, 16 * ii0 + 16 * ii2 - i0), 16 * ii0 + 16 * ii2 - i0 + i1 - 16), i2_99162_ub = min(N, 16 * ii0 + 16 * ii2 - i0 + i1 - 1);
              for (register int i2 = i2_58192_lb; i2 <= i2_99162_ub; i2 += 1) {
                S1(i0, i1, i2);
              }
              const int i2_57307_lb = 16 * ii0 + 16 * ii2 - i0 + i1, i2_67524_ub = min(N, 16 * ii0 + 16 * ii2 - i0 + 15);
              for (register int i2 = i2_57307_lb; i2 <= i2_67524_ub; i2 += 1) {
                S1(i0, i1, i2);
              }
            }
          }
        }
        if (N == 1) {
          S1(16 * ii0, 1, 1);
        }
      }
      if (T + 7 >= 16 * ii0 + 8 * k) {
        const int ii1_7238_lb = 1, ii1_49158_ub = min(-ii0 + floord(T + N, 16), (N + 8 * k) / 16);
        #pragma omp parallel for
        for (register int ii1 = ii1_7238_lb; ii1 < ii1_49158_ub; ii1 += 1) {
          const int ii2_23547_lb = 0, ii2_34489_ub = min(min(min(k - 2 * ii1 + 3 * N / 16 - 2, -ii0 - ii1 + (T + 2 * N) / 16 - 1), (N - 7 * k + 5) / 16 + 1), -ii0 - k + (T + N + 8 * k + 7) / 16);
          for (register int ii2 = ii2_23547_lb; ii2 <= ii2_34489_ub; ii2 += 1) {
            if (k == 2) {
              const int i0_98817_lb = max(max(16 * ii0 + 8, -N + 16 * ii0 + 16 * ii1 + 15), -N + 16 * ii0 + 16 * ii2 + 8), i0_72665_ub = min(T - 1, 16 * ii0 + 15);
              for (register int i0 = i0_98817_lb; i0 <= i0_72665_ub; i0 += 1) {
                const int i1_63463_lb = max(-N + 16 * ii0 + 16 * ii1 + 16 * ii2 - i0 + 15, 16 * ii0 + 16 * ii1 - i0 + 15), i1_64824_ub = min(min(N, -16 * ii0 + 16 * ii1 + i0), N - 16 * ii0 + 16 * ii1 - 16 * ii2 + i0);
                for (register int i1 = i1_63463_lb; i1 <= i1_64824_ub; i1 += 1) {
                  const int i2_80644_lb = max(max(1, 16 * ii0 - 16 * ii1 + 16 * ii2 - i0 + i1), 16 * ii0 + 16 * ii1 + 16 * ii2 - i0 - i1 + 15), i2_48941_ub = min(N, 16 * ii0 + 16 * ii1 + 16 * ii2 - i0 - i1 + 30);
                  for (register int i2 = i2_80644_lb; i2 <= i2_48941_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  const int i2_72138_lb = 16 * ii0 + 16 * ii1 + 16 * ii2 - i0 - i1 + 31, i2_54389_ub = min(N, 16 * ii0 - 16 * ii1 + 16 * ii2 - i0 + i1 + 15);
                  for (register int i2 = i2_72138_lb; i2 <= i2_54389_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
              }
            } else {
              const int i0_86149_lb = max(max(max(16 * ii0, -N + 16 * ii0 + 16 * ii1 + 15), -N + 16 * ii0 + 16 * ii2), -2 * N + 16 * ii0 + 16 * ii1 + 16 * ii2 + 15), i0_5139_ub = min(min(T - 1, 16 * ii0 + 14), N + 16 * ii0 - 16 * ii1 - 1);
              for (register int i0 = i0_86149_lb; i0 <= i0_5139_ub; i0 += 1) {
                const int i1_93154_lb = max(max(-16 * ii0 + 16 * ii1 + i0 + 1, -N + 16 * ii0 + 16 * ii1 + 16 * ii2 - i0 + 15), 16 * ii0 + 16 * ii1 - i0 + 15), i1_12687_ub = min(N, 16 * ii1 + 15);
                for (register int i1 = i1_93154_lb; i1 <= i1_12687_ub; i1 += 1) {
                  const int i2_90178_lb = max(1, 16 * ii0 + 16 * ii1 + 16 * ii2 - i0 - i1 + 15), i2_71241_ub = min(N, 16 * ii0 + 16 * ii1 + 16 * ii2 - i0 - i1 + 30);
                  for (register int i2 = i2_90178_lb; i2 <= i2_71241_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
                const int i1_73575_lb = 16 * ii1 + 16, i1_57123_ub = min(min(min(min(N, -16 * ii0 + 16 * ii1 + i0 + 16), N - 16 * ii0 + 16 * ii1 - 16 * ii2 + i0 + 16), -N + 16 * ii0 + 16 * ii1 + 16 * ii2 - i0 + 30), 16 * ii0 + 16 * ii1 - i0 + 30);
                for (register int i1 = i1_73575_lb; i1 <= i1_57123_ub; i1 += 1) {
                  const int i2_7349_lb = 16 * ii0 - 16 * ii1 + 16 * ii2 - i0 + i1 - 16, i2_92169_ub = N;
                  for (register int i2 = i2_7349_lb; i2 <= i2_92169_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
                const int i1_24257_lb = max(16 * ii1 + 16, -N + 16 * ii0 + 16 * ii1 + 16 * ii2 - i0 + 31), i1_46309_ub = min(min(N, -16 * ii0 + 16 * ii1 + i0 + 16), 16 * ii0 + 16 * ii1 - i0 + 30);
                for (register int i1 = i1_24257_lb; i1 <= i1_46309_ub; i1 += 1) {
                  const int i2_56676_lb = max(1, 16 * ii0 - 16 * ii1 + 16 * ii2 - i0 + i1 - 16), i2_62934_ub = min(N, 16 * ii0 - 16 * ii1 + 16 * ii2 - i0 + i1 - 1);
                  for (register int i2 = i2_56676_lb; i2 <= i2_62934_ub; i2 += 1) {
                    S1(i0, i1, i2);
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


#ifdef TIME
  gettimeofday(&end, 0);

  ts_return = timeval_subtract(&result, &end, &start);
  tdiff = (double)(result.tv_sec + result.tv_usec * 1.0e-6);

  printf("|Time taken =  %7.5lfms\t", tdiff * 1.0e3);
  printf("|MFLOPS =  %f\n",
         ((((double)NUM_FP_OPS * N * N * T) / tdiff) / 1000000L));
#endif

  if (fopen(".test", "r")) {
    for (int i = 1; i < N + 1; i++) {
      for (int j = 1; j < N + 1; j++) {
        total += A[T % 2][i][j];
      }
    }
    fprintf(stderr, "|sum: %e\t", total);
    for (int i = 1; i < N + 1; i++) {
      for (int j = 1; j < N + 1; j++) {
        sum_err_sqr +=
            (A[T % 2][i][j] - (total / N)) * (A[T % 2][i][j] - (total / N));
      }
    }
    fprintf(stderr, "|rms(A) = %7.2f\t", sqrt(sum_err_sqr));
    for (int i = 1; i < N + 1; i++) {
      for (int j = 1; j < N + 1; j++) {
        chtotal += ((char *)A[T % 2][i])[j];
      }
    }
    fprintf(stderr, "|sum(rep(A)) = %d\n", chtotal);
  }

  return 0;
}

// icc -O3 -fp-model precise heat_1d_np.c -o op-heat-1d-np -lm
// /* @ begin PrimeTile (num_tiling_levels=1; first_depth=1; last_depth=-1;
// boundary_tiling_level=-1;) @*/
// /* @ begin PrimeRegTile (scalar_replacement=0; T1t3=8; T1t4=8; ) @*/
// /* @ end @*/
