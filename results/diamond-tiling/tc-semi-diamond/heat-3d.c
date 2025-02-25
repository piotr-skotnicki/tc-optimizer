/*
 * Discretized 3D heat equation stencil with non periodic boundary conditions
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
#define N 150L
#define T 100L
#endif

#define NUM_FP_OPS 15

/* Define our arrays */

double A[2][N + 2][N + 2][N + 2];
double total = 0;
double sum_err_sqr = 0;
int chtotal = 0;
/* Subtract the `struct timeval' values X and Y,
 * storing the result in RESULT.
 *
 * Return 1 if the difference is negative, otherwise 0.
 */
int timeval_subtract(struct timeval *result, struct timeval *x,
                     struct timeval *y) {
  /* Perform the carry for the later subtraction by updating y. */
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

  /* Compute the time remaining to wait.
   * tv_usec is certainly positive.
   */
  result->tv_sec = x->tv_sec - y->tv_sec;
  result->tv_usec = x->tv_usec - y->tv_usec;

  /* Return 1 if result is negative. */
  return x->tv_sec < y->tv_sec;
}

int main(int argc, char *argv[]) {
  long int t, i, j, k;
  const int BASE = 1024;
  long count = 0;
  // for timekeeping
  int ts_return = -1;
  struct timeval start, end, result;
  double tdiff = 0.0;

  printf("Number of points = %ld\t|Number of timesteps = %ld\t", N, T);

  /* Initialization */
  srand(42); // seed with a constant value to verify results

  for (i = 1; i < N + 1; i++) {
    for (j = 1; j < N + 1; j++) {
      for (k = 1; k < N + 1; k++) {
        A[0][i][j][k] = 1.0 * (rand() % BASE);
      }
    }
  }

#ifdef TIME
  gettimeofday(&start, 0);
#endif

/* TC Optimizing Compiler 0.4.1 */
/* ./tc ../examples/pluto/heat-3d-t.scop.c --semi-diamond-tiling --omp-for-codegen --omega-map-tc --debug -b 8 --drop-bounds --use-macros */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define S1_I(t,i,j,k) A[(t + 1) % 2][i][j][k] = ((((0.125 * ((A[t % 2][i + 1][j][k] - (2.0 * A[t % 2][i][j][k])) + A[t % 2][i - 1][j][k])) + (0.125 * ((A[t % 2][i][j + 1][k] - (2.0 * A[t % 2][i][j][k])) + A[t % 2][i][j - 1][k]))) + (0.125 * ((A[t % 2][i][j][k - 1] - (2.0 * A[t % 2][i][j][k])) + A[t % 2][i][j][k + 1]))) + A[t % 2][i][j][k]);
#define S1(t,i,j,k) S1_I((t),(i),(j),(k))
#pragma scop
const int ii0_96143_lb = 0, ii0_29352_ub = floord(T - 2, 8);
for (register int ii0 = ii0_96143_lb; ii0 <= ii0_29352_ub; ii0 += 1) {
  const int ii1_22408_lb = 0, ii1_4804_ub = floord(N - 1, 8);
  #pragma omp parallel for
  for (register int ii1 = ii1_22408_lb; ii1 < ii1_4804_ub; ii1 += 1) {
    const int ii2_28134_lb = 0, ii2_75485_ub = min(-ii0 + (T + N - 2) / 8, (N + 2) / 8);
    for (register int ii2 = ii2_28134_lb; ii2 <= ii2_75485_ub; ii2 += 1) {
      const int ii3_700_lb = 0, ii3_53484_ub = min(-ii0 + (T + N - 2) / 8, (N + 2) / 8);
      for (register int ii3 = ii3_700_lb; ii3 <= ii3_53484_ub; ii3 += 1) {
        const int i0_82164_lb = max(max(8 * ii0, -N + 8 * ii0 + 8 * ii2), -N + 8 * ii0 + 8 * ii3), i0_94144_ub = min(min(T - 2, 8 * ii0 + 2), N + 8 * ii0 - 8 * ii1 - 9);
        for (register int i0 = i0_82164_lb; i0 <= i0_94144_ub; i0 += 1) {
          const int i1_81675_lb = -8 * ii0 + 8 * ii1 + i0 + 9, i1_90644_ub = min(N, 8 * ii0 + 8 * ii1 - i0 + 14);
          for (register int i1 = i1_81675_lb; i1 <= i1_90644_ub; i1 += 1) {
            const int i2_96983_lb = max(1, 8 * ii0 + 8 * ii2 - i0), i2_9999_ub = min(N, 8 * ii0 + 8 * ii2 - i0 + 7);
            for (register int i2 = i2_96983_lb; i2 <= i2_9999_ub; i2 += 1) {
              const int i3_80412_lb = max(1, 8 * ii0 + 8 * ii3 - i0), i3_94329_ub = min(N, 8 * ii0 + 8 * ii3 - i0 + 7);
              for (register int i3 = i3_80412_lb; i3 <= i3_94329_ub; i3 += 1) {
                S1(i0, i1, i2, i3);
              }
            }
          }
        }
      }
    }
  }
  const int k_23748_lb = 1, k_34954_ub = min(2, -ii0 + (T + 5) / 8 + 1);
  for (register int k = k_23748_lb; k <= k_34954_ub; k += 1) {
    const int ii2_99360_lb = 0, ii2_35645_ub = min(min(min(k + floord(N - k, 4) - 1, -ii0 + floord(T + N - k - 1, 8)), (N - k + 8) / 8), (2 * N + 5 * k + 8) / 16);
    for (register int ii2 = ii2_99360_lb; ii2 <= ii2_35645_ub; ii2 += 1) {
      const int ii3_87596_lb = 0, ii3_36056_ub = min(min(min(k + floord(N - k, 4) - 1, -ii0 + (T + N - k - 1) / 8), (N - k + 8) / 8), (2 * N + 5 * k + 8) / 16);
      for (register int ii3 = ii3_87596_lb; ii3 <= ii3_36056_ub; ii3 += 1) {
        const int i0_86612_lb = max(max(max(8 * ii0, 8 * ii0 + 8 * k - 15), -N + 8 * ii0 + k + 8 * ii2 - 1), -N + 8 * ii0 + k + 8 * ii3 - 1), i0_90156_ub = min(min(min(T - 2, 8 * ii0 + 7), N + 8 * ii0 + 7 * k - 8), 8 * ii0 + 4 * k + 2);
        for (register int i0 = i0_86612_lb; i0 <= i0_90156_ub; i0 += 1) {
          const int i1_29953_lb = max(1, -8 * ii0 - 7 * k + i0 + 8), i1_56473_ub = min(min(N, 8 * ii0 - i0 + 14), -8 * ii0 - 8 * k + i0 + 16);
          for (register int i1 = i1_29953_lb; i1 <= i1_56473_ub; i1 += 1) {
            const int i2_2756_lb = max(1, 8 * ii0 + k + 8 * ii2 - i0 - 1), i2_4057_ub = min(N, 8 * ii0 + k + 8 * ii2 - i0 + 6);
            for (register int i2 = i2_2756_lb; i2 <= i2_4057_ub; i2 += 1) {
              const int i3_81945_lb = max(1, 8 * ii0 + k + 8 * ii3 - i0 - 1), i3_49005_ub = min(N, 8 * ii0 + k + 8 * ii3 - i0 + 6);
              for (register int i3 = i3_81945_lb; i3 <= i3_49005_ub; i3 += 1) {
                S1(i0, i1, i2, i3);
              }
            }
          }
        }
      }
    }
    if (T + 2 >= 8 * ii0 + 4 * k) {
      const int ii1_73410_lb = 1, ii1_78088_ub = min(-ii0 + floord(T + N - 1, 8), (N + 4 * k) / 8);
      #pragma omp parallel for
      for (register int ii1 = ii1_73410_lb; ii1 < ii1_78088_ub; ii1 += 1) {
        const int ii2_94709_lb = 0, ii2_95819_ub = min(min((N - 3 * k + 1) / 8 + 1, k - ii1 + (2 * N - 3 * k + 1) / 8 - 1), -ii0 - k + (T + N + 4 * k + 2) / 8);
        for (register int ii2 = ii2_94709_lb; ii2 <= ii2_95819_ub; ii2 += 1) {
          const int ii3_82892_lb = 0, ii3_22844_ub = min(min((N - 3 * k + 1) / 8 + 1, k - ii1 + (2 * N - 3 * k + 1) / 8 - 1), -ii0 - k + (T + N + 4 * k + 2) / 8);
          for (register int ii3 = ii3_82892_lb; ii3 <= ii3_22844_ub; ii3 += 1) {
            if (k == 2) {
              const int i0_71304_lb = max(max(max(8 * ii0 + 4, -N + 8 * ii0 + 8 * ii1 + 7), -N + 8 * ii0 + 8 * ii2 + 4), -N + 8 * ii0 + 8 * ii3 + 4), i0_99944_ub = min(T - 2, 8 * ii0 + 7);
              for (register int i0 = i0_71304_lb; i0 <= i0_99944_ub; i0 += 1) {
                const int i1_92680_lb = 8 * ii0 + 8 * ii1 - i0 + 7, i1_53468_ub = min(N, -8 * ii0 + 8 * ii1 + i0);
                for (register int i1 = i1_92680_lb; i1 <= i1_53468_ub; i1 += 1) {
                  const int i2_94088_lb = max(1, 8 * ii0 + 8 * ii2 - i0 + 4), i2_90707_ub = min(N, 8 * ii0 + 8 * ii2 - i0 + 11);
                  for (register int i2 = i2_94088_lb; i2 <= i2_90707_ub; i2 += 1) {
                    const int i3_60464_lb = max(1, 8 * ii0 + 8 * ii3 - i0 + 4), i3_7423_ub = min(N, 8 * ii0 + 8 * ii3 - i0 + 11);
                    for (register int i3 = i3_60464_lb; i3 <= i3_7423_ub; i3 += 1) {
                      S1(i0, i1, i2, i3);
                    }
                  }
                }
              }
            } else {
              const int i0_17059_lb = max(max(max(8 * ii0, -N + 8 * ii0 + 8 * ii1 + 7), -N + 8 * ii0 + 8 * ii2), -N + 8 * ii0 + 8 * ii3), i0_40877_ub = min(min(T - 2, 8 * ii0 + 6), N + 8 * ii0 - 8 * ii1 - 1);
              for (register int i0 = i0_17059_lb; i0 <= i0_40877_ub; i0 += 1) {
                const int i1_1753_lb = max(-8 * ii0 + 8 * ii1 + i0 + 1, 8 * ii0 + 8 * ii1 - i0 + 7), i1_40807_ub = min(min(N, -8 * ii0 + 8 * ii1 + i0 + 8), 8 * ii0 + 8 * ii1 - i0 + 14);
                for (register int i1 = i1_1753_lb; i1 <= i1_40807_ub; i1 += 1) {
                  const int i2_92183_lb = max(1, 8 * ii0 + 8 * ii2 - i0), i2_1113_ub = min(N, 8 * ii0 + 8 * ii2 - i0 + 7);
                  for (register int i2 = i2_92183_lb; i2 <= i2_1113_ub; i2 += 1) {
                    const int i3_76452_lb = max(1, 8 * ii0 + 8 * ii3 - i0), i3_79779_ub = min(N, 8 * ii0 + 8 * ii3 - i0 + 7);
                    for (register int i3 = i3_76452_lb; i3 <= i3_79779_ub; i3 += 1) {
                      S1(i0, i1, i2, i3);
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

#ifdef TIME
  gettimeofday(&end, 0);

  ts_return = timeval_subtract(&result, &end, &start);
  tdiff = (double)(result.tv_sec + result.tv_usec * 1.0e-6);

  printf("|Time taken: %7.5lfms\t", tdiff * 1.0e3);
  printf("|MFLOPS: %f\n",
         ((((double)NUM_FP_OPS * N * N * N * (T - 1)) / tdiff) / 1000000L));
#endif

#ifdef VERIFY
  for (i = 1; i < N + 1; i++) {
    for (j = 1; j < N + 1; j++) {
      for (k = 1; k < N + 1; k++) {
        total += A[T % 2][i][j][k];
      }
    }
  }
  fprintf(stderr, "|sum: %e\t", total);
  for (i = 1; i < N + 1; i++) {
    for (j = 1; j < N + 1; j++) {
      for (k = 1; k < N + 1; k++) {
        sum_err_sqr += (A[T % 2][i][j][k] - (total / N)) *
                       (A[T % 2][i][j][k] - (total / N));
      }
    }
  }
  fprintf(stderr, "|rms(A) = %7.2f\t", sqrt(sum_err_sqr));
  for (i = 1; i < N + 1; i++) {
    for (j = 1; j < N + 1; j++) {
      for (k = 1; k < N + 1; k++) {
        chtotal += ((char *)A[T % 2][i][j])[k];
      }
    }
  }
  fprintf(stderr, "|sum(rep(A)) = %d\n", chtotal);
#endif

  return 0;
}

// icc -O3 -fp-model precise heat_1d_np.c -o op-heat-1d-np -lm
// /* @ begin PrimeTile (num_tiling_levels=1; first_depth=1; last_depth=-1;
// boundary_tiling_level=-1;) @*/
// /* @ begin PrimeRegTile (scalar_replacement=0; T1t5=4; T1t6=4; T1t7=4;
// T1t8=4; ) @*/
// /* @ end @*/
