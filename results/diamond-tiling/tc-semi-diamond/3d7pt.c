/*
 * Order-1, 3D 7 point stencil
 * Adapted from Pochoir test bench
 *
 * Irshad Pananilath: irshad@csa.iisc.ernet.in
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define NUM_FP_OPS 8

#define MAX(a,b) ((a) > (b) ? a : b)
#define MIN(a,b) ((a) < (b) ? a : b)


#define N 256L
#define T 200L

double A[2][N+2][N+2][N+2];

/* Subtract the `struct timeval' values X and Y,
 * storing the result in RESULT.
 *
 * Return 1 if the difference is negative, otherwise 0.
 */
int timeval_subtract(struct timeval *result, struct timeval *x, struct timeval *y)
{
  /* Perform the carry for the later subtraction by updating y. */
  if (x->tv_usec < y->tv_usec)
  {
    int nsec = (y->tv_usec - x->tv_usec) / 1000000 + 1;

    y->tv_usec -= 1000000 * nsec;
    y->tv_sec += nsec;
  }

  if (x->tv_usec - y->tv_usec > 1000000)
  {
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

int main(int argc, char *argv[])
{
  // for timekeeping
  int ts_return = -1;
  struct timeval start, end, result;
  double tdiff = 0.0;

  int t, i, j, k;

  const int BASE = 1024;
  const double alpha = 0.0876;
  const double beta = 0.0765;
  double total;

    printf("Number of points = %ld\t|Number of timesteps = %ld\t", N, T);

  // initialize variables
  //
  srand(42);
    for (i = 1; i < N+1; i++) {
        for (j = 1; j < N+1; j++) {
            for (k = 1; k < N+1; k++) {
                A[0][i][j][k] = 1.0 * (rand() % BASE);
            }
        }
    }


#ifdef TIME
  gettimeofday(&start, 0);
#endif

/* TC Optimizing Compiler 0.4.1 */
/* ./tc ../examples/pluto/3d7pt-t.scop.c --semi-diamond-tiling --omp-for-codegen --omega-map-tc --debug --use-macros -b 16 --drop-bounds */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define S1_I(t,i,j,k) A[(t + 1) % 2][i][j][k] = ((alpha * A[t % 2][i][j][k]) + (beta * (((((A[t % 2][i - 1][j][k] + A[t % 2][i][j - 1][k]) + A[t % 2][i][j][k - 1]) + A[t % 2][i + 1][j][k]) + A[t % 2][i][j + 1][k]) + A[t % 2][i][j][k + 1])));
#define S1(t,i,j,k) S1_I((t),(i),(j),(k))
#pragma scop
const int ii0_48423_lb = 0, ii0_13252_ub = floord(T - 2, 16);
for (register int ii0 = ii0_48423_lb; ii0 <= ii0_13252_ub; ii0 += 1) {
  const int ii1_41645_lb = 0, ii1_97347_ub = floord(N - 1, 16);
  #pragma omp parallel for
  for (register int ii1 = ii1_41645_lb; ii1 < ii1_97347_ub; ii1 += 1) {
    const int ii2_54647_lb = 0, ii2_5366_ub = min(-ii0 + (T + N - 2) / 16, (N + 6) / 16);
    for (register int ii2 = ii2_54647_lb; ii2 <= ii2_5366_ub; ii2 += 1) {
      const int ii3_84595_lb = 0, ii3_19292_ub = min(-ii0 + (T + N - 2) / 16, (N + 6) / 16);
      for (register int ii3 = ii3_84595_lb; ii3 <= ii3_19292_ub; ii3 += 1) {
        const int i0_38989_lb = max(max(16 * ii0, -N + 16 * ii0 + 16 * ii2), -N + 16 * ii0 + 16 * ii3), i0_99019_ub = min(min(T - 2, 16 * ii0 + 6), N + 16 * ii0 - 16 * ii1 - 17);
        for (register int i0 = i0_38989_lb; i0 <= i0_99019_ub; i0 += 1) {
          const int i1_46224_lb = -16 * ii0 + 16 * ii1 + i0 + 17, i1_9780_ub = min(N, 16 * ii0 + 16 * ii1 - i0 + 30);
          for (register int i1 = i1_46224_lb; i1 <= i1_9780_ub; i1 += 1) {
            const int i2_55958_lb = max(1, 16 * ii0 + 16 * ii2 - i0), i2_60116_ub = min(N, 16 * ii0 + 16 * ii2 - i0 + 15);
            for (register int i2 = i2_55958_lb; i2 <= i2_60116_ub; i2 += 1) {
              const int i3_5454_lb = max(1, 16 * ii0 + 16 * ii3 - i0), i3_7840_ub = min(N, 16 * ii0 + 16 * ii3 - i0 + 15);
              for (register int i3 = i3_5454_lb; i3 <= i3_7840_ub; i3 += 1) {
                S1(i0, i1, i2, i3);
              }
            }
          }
        }
      }
    }
  }
  const int k_43404_lb = 1, k_80118_ub = min(2, -ii0 + (T + 13) / 16 + 1);
  for (register int k = k_43404_lb; k <= k_80118_ub; k += 1) {
    const int ii2_68291_lb = 0, ii2_20065_ub = min(min(min(k + floord(N - k, 8) - 1, -ii0 + floord(T + N - k - 1, 16)), (N - k + 16) / 16), (2 * N + 13 * k + 16) / 32);
    for (register int ii2 = ii2_68291_lb; ii2 <= ii2_20065_ub; ii2 += 1) {
      const int ii3_54939_lb = 0, ii3_75343_ub = min(min(min(k + floord(N - k, 8) - 1, -ii0 + (T + N - k - 1) / 16), (N - k + 16) / 16), (2 * N + 13 * k + 16) / 32);
      for (register int ii3 = ii3_54939_lb; ii3 <= ii3_75343_ub; ii3 += 1) {
        const int i0_75587_lb = max(max(max(16 * ii0, 16 * ii0 + 16 * k - 31), -N + 16 * ii0 + k + 16 * ii2 - 1), -N + 16 * ii0 + k + 16 * ii3 - 1), i0_83239_ub = min(min(min(T - 2, 16 * ii0 + 15), N + 16 * ii0 + 15 * k - 16), 16 * ii0 + 8 * k + 6);
        for (register int i0 = i0_75587_lb; i0 <= i0_83239_ub; i0 += 1) {
          const int i1_59189_lb = max(1, -16 * ii0 - 15 * k + i0 + 16), i1_86792_ub = min(min(N, 16 * ii0 - i0 + 30), -16 * ii0 - 16 * k + i0 + 32);
          for (register int i1 = i1_59189_lb; i1 <= i1_86792_ub; i1 += 1) {
            const int i2_89196_lb = max(1, 16 * ii0 + k + 16 * ii2 - i0 - 1), i2_35070_ub = min(N, 16 * ii0 + k + 16 * ii2 - i0 + 14);
            for (register int i2 = i2_89196_lb; i2 <= i2_35070_ub; i2 += 1) {
              const int i3_88945_lb = max(1, 16 * ii0 + k + 16 * ii3 - i0 - 1), i3_20671_ub = min(N, 16 * ii0 + k + 16 * ii3 - i0 + 14);
              for (register int i3 = i3_88945_lb; i3 <= i3_20671_ub; i3 += 1) {
                S1(i0, i1, i2, i3);
              }
            }
          }
        }
      }
    }
    if (T + 6 >= 16 * ii0 + 8 * k) {
      const int ii1_61103_lb = 1, ii1_37369_ub = min(-ii0 + floord(T + N - 1, 16), (N + 8 * k) / 16);
      #pragma omp parallel for
      for (register int ii1 = ii1_61103_lb; ii1 < ii1_37369_ub; ii1 += 1) {
        const int ii2_50275_lb = 0, ii2_19100_ub = min(min((N - 7 * k + 5) / 16 + 1, k - ii1 + (2 * N - 7 * k + 5) / 16 - 1), -ii0 - k + (T + N + 8 * k + 6) / 16);
        for (register int ii2 = ii2_50275_lb; ii2 <= ii2_19100_ub; ii2 += 1) {
          const int ii3_34716_lb = 0, ii3_4923_ub = min(min((N - 7 * k + 5) / 16 + 1, k - ii1 + (2 * N - 7 * k + 5) / 16 - 1), -ii0 - k + (T + N + 8 * k + 6) / 16);
          for (register int ii3 = ii3_34716_lb; ii3 <= ii3_4923_ub; ii3 += 1) {
            if (k == 2) {
              const int i0_40818_lb = max(max(max(16 * ii0 + 8, -N + 16 * ii0 + 16 * ii1 + 15), -N + 16 * ii0 + 16 * ii2 + 8), -N + 16 * ii0 + 16 * ii3 + 8), i0_35663_ub = min(T - 2, 16 * ii0 + 15);
              for (register int i0 = i0_40818_lb; i0 <= i0_35663_ub; i0 += 1) {
                const int i1_24215_lb = 16 * ii0 + 16 * ii1 - i0 + 15, i1_96160_ub = min(N, -16 * ii0 + 16 * ii1 + i0);
                for (register int i1 = i1_24215_lb; i1 <= i1_96160_ub; i1 += 1) {
                  const int i2_34683_lb = max(1, 16 * ii0 + 16 * ii2 - i0 + 8), i2_86791_ub = min(N, 16 * ii0 + 16 * ii2 - i0 + 23);
                  for (register int i2 = i2_34683_lb; i2 <= i2_86791_ub; i2 += 1) {
                    const int i3_5940_lb = max(1, 16 * ii0 + 16 * ii3 - i0 + 8), i3_6993_ub = min(N, 16 * ii0 + 16 * ii3 - i0 + 23);
                    for (register int i3 = i3_5940_lb; i3 <= i3_6993_ub; i3 += 1) {
                      S1(i0, i1, i2, i3);
                    }
                  }
                }
              }
            } else {
              const int i0_63259_lb = max(max(max(16 * ii0, -N + 16 * ii0 + 16 * ii1 + 15), -N + 16 * ii0 + 16 * ii2), -N + 16 * ii0 + 16 * ii3), i0_11394_ub = min(min(T - 2, 16 * ii0 + 14), N + 16 * ii0 - 16 * ii1 - 1);
              for (register int i0 = i0_63259_lb; i0 <= i0_11394_ub; i0 += 1) {
                const int i1_14833_lb = max(-16 * ii0 + 16 * ii1 + i0 + 1, 16 * ii0 + 16 * ii1 - i0 + 15), i1_23015_ub = min(min(N, -16 * ii0 + 16 * ii1 + i0 + 16), 16 * ii0 + 16 * ii1 - i0 + 30);
                for (register int i1 = i1_14833_lb; i1 <= i1_23015_ub; i1 += 1) {
                  const int i2_91512_lb = max(1, 16 * ii0 + 16 * ii2 - i0), i2_99476_ub = min(N, 16 * ii0 + 16 * ii2 - i0 + 15);
                  for (register int i2 = i2_91512_lb; i2 <= i2_99476_ub; i2 += 1) {
                    const int i3_43081_lb = max(1, 16 * ii0 + 16 * ii3 - i0), i3_62804_ub = min(N, 16 * ii0 + 16 * ii3 - i0 + 15);
                    for (register int i3 = i3_43081_lb; i3 <= i3_62804_ub; i3 += 1) {
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
  tdiff = (double) (result.tv_sec + result.tv_usec * 1.0e-6);

  printf("Time taken: %7.5lfms\n", tdiff * 1.0e3);
  printf("|MFLOPS: %f\t", ((((double)NUM_FP_OPS * N *N * N * (T-1)) / tdiff) / 1000000L));
#endif

#ifdef VERIFY
  // display the initial setting
  total = 0.0;
  for (i = 0; i < N; ++i) {
    for (j = 0; j < N; ++j) {
      for (k = 0; k < N; ++k) {
        total += A[T%2][i][j][k];
      }
    }
  }
  fprintf(stderr, "Sum(final): %e\n", total);
#endif

  return 0;
}

// icc -O3 -DTIME -DDEBUG 3d7pt.c -o op-3d7pt
