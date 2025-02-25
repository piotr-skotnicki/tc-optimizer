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
/* ./tc ../examples/pluto/heat-2d-t.scop.c --diamond-tiling --omp-for-codegen --iterative-tc --debug -b 16 --drop-bounds --use-macros */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define S1_I(t,i,j) A[(t + 1)%2][i][j] = (((0.125 * ((A[t%2][i + 1][j] - (2.0 * A[t%2][i][j])) + A[t%2][i - 1][j])) + (0.125 * ((A[t%2][i][j + 1] - (2.0 * A[t%2][i][j])) + A[t%2][i][j - 1]))) + A[t%2][i][j]);
#define S1(t,i,j) S1_I((t),(i),(j))
#pragma scop
if (T >= 9 && T + N <= 31) {
  const int ii2_76461_lb = 0, ii2_86940_ub = floord(N - 1, 8);
  for (register int ii2 = ii2_76461_lb; ii2 < ii2_86940_ub; ii2 += 1) {
    const int i0_56425_lb = 0, i0_4658_ub = N - 16;
    for (register int i0 = i0_56425_lb; i0 < i0_4658_ub; i0 += 1) {
      const int i1_87318_lb = i0 + 17, i1_16508_ub = N;
      for (register int i1 = i1_87318_lb; i1 <= i1_16508_ub; i1 += 1) {
        const int i2_45718_lb = max(1, 16 * ii2 - i0), i2_46330_ub = min(N, 16 * ii2 - i0 + 15);
        for (register int i2 = i2_45718_lb; i2 <= i2_46330_ub; i2 += 1) {
          S1(i0, i1, i2);
        }
      }
    }
  }
  for (register int k = 1; k <= 2; k += 1) {
    if (k == 2) {
      const int ii2_16891_lb = 0, ii2_58939_ub = (T + N - 2) / 16;
      for (register int ii2 = ii2_16891_lb; ii2 <= ii2_58939_ub; ii2 += 1) {
        const int i0_6878_lb = max(1, -N + 16 * ii2 + 1), i0_31196_ub = min(T - 1, N + 16 * ii2 + 14);
        for (register int i0 = i0_6878_lb; i0 <= i0_31196_ub; i0 += 1) {
          const int i1_96943_lb = max(1, -16 * ii2 + i0 - 14), i1_4467_ub = min(min(N, i0), N - 16 * ii2 + i0);
          for (register int i1 = i1_96943_lb; i1 <= i1_4467_ub; i1 += 1) {
            const int i2_54556_lb = max(1, 16 * ii2 - i0 + i1), i2_37862_ub = min(N, 16 * ii2 - i0 + i1 + 15);
            for (register int i2 = i2_54556_lb; i2 <= i2_37862_ub; i2 += 1) {
              S1(i0, i1, i2);
            }
          }
        }
      }
    } else {
      const int i0_44695_lb = 0, i0_24155_ub = min(T, N);
      for (register int i0 = i0_44695_lb; i0 < i0_24155_ub; i0 += 1) {
        const int i1_90411_lb = i0 + 1, i1_6633_ub = min(N, i0 + 16);
        for (register int i1 = i1_90411_lb; i1 <= i1_6633_ub; i1 += 1) {
          const int i2_12377_lb = 1, i2_49097_ub = -i0 + i1;
          for (register int i2 = i2_12377_lb; i2 < i2_49097_ub; i2 += 1) {
            S1(i0, i1, i2);
          }
          const int i2_1800_lb = -i0 + i1, i2_62285_ub = min(N, -i0 + 15);
          for (register int i2 = i2_1800_lb; i2 <= i2_62285_ub; i2 += 1) {
            S1(i0, i1, i2);
          }
        }
      }
      const int i0_5153_lb = max(0, -N + 16), i0_81346_ub = min(T, N);
      for (register int i0 = i0_5153_lb; i0 < i0_81346_ub; i0 += 1) {
        const int i1_18073_lb = i0 + 1, i1_18637_ub = min(N, i0 + 16);
        for (register int i1 = i1_18073_lb; i1 <= i1_18637_ub; i1 += 1) {
          const int i2_49994_lb = max(-i0 + 16, -i0 + i1), i2_94534_ub = N;
          for (register int i2 = i2_49994_lb; i2 <= i2_94534_ub; i2 += 1) {
            S1(i0, i1, i2);
          }
        }
      }
    }
  }
} else {
  if (T + N >= 32 && N >= 1) {
    if (N <= 15) {
      for (register int k = 1; k <= 2; k += 1) {
        if (k == 2) {
          const int ii2_15544_lb = 0, ii2_9247_ub = min(-((-N + 19) / 16) + 2, (T + N - 2) / 16);
          for (register int ii2 = ii2_15544_lb; ii2 <= ii2_9247_ub; ii2 += 1) {
            const int i0_39280_lb = max(1, -N + 16 * ii2 + 1), i0_61263_ub = min(min(min(29, T - 1), N + 16 * ii2 + 14), 8 * ii2 + 22);
            for (register int i0 = i0_39280_lb; i0 <= i0_61263_ub; i0 += 1) {
              const int i1_55578_lb = max(1, -16 * ii2 + i0 - 14), i1_32268_ub = min(min(min(N, i0), N - 16 * ii2 + i0), -i0 + 30);
              for (register int i1 = i1_55578_lb; i1 <= i1_32268_ub; i1 += 1) {
                const int i2_65504_lb = max(1, 16 * ii2 - i0 + i1), i2_88821_ub = min(N, 16 * ii2 - i0 + i1 + 15);
                for (register int i2 = i2_65504_lb; i2 <= i2_88821_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
            }
          }
        } else {
          for (register int i0 = 0; i0 < N; i0 += 1) {
            const int i1_20017_lb = i0 + 1, i1_4503_ub = N;
            for (register int i1 = i1_20017_lb; i1 <= i1_4503_ub; i1 += 1) {
              const int i2_93202_lb = 1, i2_90925_ub = min(N, -i0 + 15);
              for (register int i2 = i2_93202_lb; i2 <= i2_90925_ub; i2 += 1) {
                S1(i0, i1, i2);
              }
            }
          }
          const int i0_42365_lb = -N + 16, i0_37897_ub = N;
          for (register int i0 = i0_42365_lb; i0 < i0_37897_ub; i0 += 1) {
            const int i1_15081_lb = i0 + 1, i1_49129_ub = N;
            for (register int i1 = i1_15081_lb; i1 <= i1_49129_ub; i1 += 1) {
              const int i2_60882_lb = -i0 + 16, i2_43810_ub = N;
              for (register int i2 = i2_60882_lb; i2 <= i2_43810_ub; i2 += 1) {
                S1(i0, i1, i2);
              }
            }
          }
        }
      }
    }
    const int ii0_98226_lb = max(0, -((N + 14) / 15) + 2), ii0_62682_ub = min(min(floord(T - 1, 8) - 1, floord(T - 1, 16)), (T + N) / 16 - 1);
    for (register int ii0 = ii0_98226_lb; ii0 <= ii0_62682_ub; ii0 += 1) {
      if (N >= 24 && ii0 == 0) {
        const int ii1_22447_lb = 0, ii1_19731_ub = (N - 1) / 16;
        #pragma omp parallel for
        for (register int ii1 = ii1_22447_lb; ii1 < ii1_19731_ub; ii1 += 1) {
          const int ii2_44028_lb = 0, ii2_40521_ub = (N + 6) / 16;
          for (register int ii2 = ii2_44028_lb; ii2 <= ii2_40521_ub; ii2 += 1) {
            const int i0_38369_lb = max(0, -N + 16 * ii2), i0_10375_ub = min(6, N - 16 * ii1 - 17);
            for (register int i0 = i0_38369_lb; i0 <= i0_10375_ub; i0 += 1) {
              const int i1_51407_lb = 16 * ii1 + i0 + 17, i1_60298_ub = min(N, 16 * ii1 - i0 + 30);
              for (register int i1 = i1_51407_lb; i1 <= i1_60298_ub; i1 += 1) {
                const int i2_33147_lb = max(1, 16 * ii2 - i0), i2_66952_ub = min(N, 16 * ii2 - i0 + 15);
                for (register int i2 = i2_33147_lb; i2 <= i2_66952_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
            }
          }
        }
        const int ii1_69546_lb = 0, ii1_72427_ub = (N + 8) / 16;
        #pragma omp parallel for
        for (register int ii1 = ii1_69546_lb; ii1 < ii1_72427_ub; ii1 += 1) {
          const int ii2_28215_lb = 0, ii2_41476_ub = min(min(min((N - 2) / 16 + 1, (T + N - 1) / 16), -2 * ii1 + 3 * N / 16 - 1), -ii1 + (T + 2 * N) / 16 - 1);
          for (register int ii2 = ii2_28215_lb; ii2 <= ii2_41476_ub; ii2 += 1) {
            const int i0_4696_lb = max(max(max(0, -N + 16 * ii1 + 15), -N + 16 * ii2), -2 * N + 16 * ii1 + 16 * ii2 + 15), i0_10071_ub = min(min(14, T - 1), N - 16 * ii1 - 1);
            for (register int i0 = i0_4696_lb; i0 <= i0_10071_ub; i0 += 1) {
              {
                const int i1_30297_lb = max(max(max(16 * ii1 + i0 + 1, -N + 16 * ii1 + 16 * ii2 - i0 + 15), 16 * ii1 - i0 + 15), -i0 + 31), i1_12255_ub = min(N, 16 * ii1 + 15);
                for (register int i1 = i1_30297_lb; i1 <= i1_12255_ub; i1 += 1) {
                  const int i2_82454_lb = max(1, 16 * ii1 + 16 * ii2 - i0 - i1 + 15), i2_66667_ub = min(N, 16 * ii1 + 16 * ii2 - i0 - i1 + 30);
                  for (register int i2 = i2_82454_lb; i2 <= i2_66667_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
                const int i1_16758_lb = max(16 * ii1 + 16, -i0 + 31), i1_75656_ub = min(min(min(min(N, 16 * ii1 + i0 + 16), N + 16 * ii1 - 16 * ii2 + i0 + 16), -N + 16 * ii1 + 16 * ii2 - i0 + 30), 16 * ii1 - i0 + 30);
                for (register int i1 = i1_16758_lb; i1 <= i1_75656_ub; i1 += 1) {
                  const int i2_57592_lb = -16 * ii1 + 16 * ii2 - i0 + i1 - 16, i2_75476_ub = N;
                  for (register int i2 = i2_57592_lb; i2 <= i2_75476_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
                const int i1_29905_lb = max(max(16 * ii1 + 16, -N + 16 * ii1 + 16 * ii2 - i0 + 31), -i0 + 31), i1_72673_ub = min(min(N, 16 * ii1 + i0 + 16), 16 * ii1 - i0 + 30);
                for (register int i1 = i1_29905_lb; i1 <= i1_72673_ub; i1 += 1) {
                  const int i2_40957_lb = max(1, -16 * ii1 + 16 * ii2 - i0 + i1 - 16), i2_90788_ub = min(N, -16 * ii1 + 16 * ii2 - i0 + i1 - 1);
                  for (register int i2 = i2_40957_lb; i2 <= i2_90788_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
              }
              if (ii1 == 0) {
                if (ii2 == 0) {
                  const int i2_32835_lb = 1, i2_39183_ub = -i0 + 15;
                  for (register int i2 = i2_32835_lb; i2 <= i2_39183_ub; i2 += 1) {
                    S1(i0, i0 + 1, i2);
                  }
                }
                const int i1_53470_lb = max(i0 + 1, -16 * ii2 + i0 + 2), i1_55283_ub = min(min(i0 + 16, N - 16 * ii2 + i0 + 16), -i0 + 30);
                for (register int i1 = i1_53470_lb; i1 <= i1_55283_ub; i1 += 1) {
                  const int i2_75267_lb = max(max(1, 16 * ii2 - i0), 16 * ii2 - i0 + i1 - 16), i2_13851_ub = min(N, 16 * ii2 - i0 + i1 - 1);
                  for (register int i2 = i2_75267_lb; i2 <= i2_13851_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  const int i2_95804_lb = 16 * ii2 - i0 + i1, i2_29988_ub = min(N, 16 * ii2 - i0 + 15);
                  for (register int i2 = i2_95804_lb; i2 <= i2_29988_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
              }
            }
          }
        }
      } else if (N <= 23 && ii0 == 0) {
        const int k_24226_lb = max(0, -N + 17), k_63563_ub = 1;
        for (register int k = k_24226_lb; k <= k_63563_ub; k += 1) {
          const int ii2_90286_lb = 0, ii2_57373_ub = min(min(k + (N - 1) / 8 - 1, (N - 2) / 16 + 1), (T + N - 1) / 16);
          for (register int ii2 = ii2_90286_lb; ii2 <= ii2_57373_ub; ii2 += 1) {
            if (k == 1) {
              const int i0_30515_lb = max(0, -N + 16 * ii2), i0_76184_ub = min(14, T - 1);
              for (register int i0 = i0_30515_lb; i0 <= i0_76184_ub; i0 += 1) {
                if (ii2 == 0) {
                  const int i2_46152_lb = 1, i2_75082_ub = -i0 + 15;
                  for (register int i2 = i2_46152_lb; i2 <= i2_75082_ub; i2 += 1) {
                    S1(i0, i0 + 1, i2);
                  }
                }
                const int i1_34012_lb = max(i0 + 1, -16 * ii2 + i0 + 2), i1_67200_ub = min(min(min(N, i0 + 16), N - 16 * ii2 + i0 + 16), -i0 + 30);
                for (register int i1 = i1_34012_lb; i1 <= i1_67200_ub; i1 += 1) {
                  const int i2_85154_lb = max(max(1, 16 * ii2 - i0), 16 * ii2 - i0 + i1 - 16), i2_64310_ub = min(N, 16 * ii2 - i0 + i1 - 1);
                  for (register int i2 = i2_85154_lb; i2 <= i2_64310_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  const int i2_79456_lb = 16 * ii2 - i0 + i1, i2_83960_ub = min(N, 16 * ii2 - i0 + 15);
                  for (register int i2 = i2_79456_lb; i2 <= i2_83960_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
              }
            } else {
              const int i0_47329_lb = 0, i0_96214_ub = N - 16;
              for (register int i0 = i0_47329_lb; i0 < i0_96214_ub; i0 += 1) {
                const int i1_59617_lb = i0 + 17, i1_4921_ub = N;
                for (register int i1 = i1_59617_lb; i1 <= i1_4921_ub; i1 += 1) {
                  const int i2_71690_lb = max(1, 16 * ii2 - i0), i2_5874_ub = min(N, 16 * ii2 - i0 + 15);
                  for (register int i2 = i2_71690_lb; i2 <= i2_5874_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
              }
            }
          }
        }
      }
      const int k_93947_lb = max(max(1, -ii0 + 2), -((N + 8) / 8) + 3), k_12647_ub = min(2, -2 * ii0 + (2 * T - 2 * ii0 - 2) / 15 + 1);
      for (register int k = k_93947_lb; k <= k_12647_ub; k += 1) {
        if (ii0 == 0 && k == 2) {
          const int ii2_96662_lb = 0, ii2_26782_ub = min((N - 4) / 16 + 2, (T + N - 2) / 16);
          for (register int ii2 = ii2_96662_lb; ii2 <= ii2_26782_ub; ii2 += 1) {
            const int i0_68183_lb = max(1, -N + 16 * ii2 + 1), i0_66485_ub = min(min(29, T - 1), 8 * ii2 + 22);
            for (register int i0 = i0_68183_lb; i0 <= i0_66485_ub; i0 += 1) {
              const int i1_82065_lb = max(1, -16 * ii2 + i0 - 14), i1_59802_ub = min(min(i0, N - 16 * ii2 + i0), -i0 + 30);
              for (register int i1 = i1_82065_lb; i1 <= i1_59802_ub; i1 += 1) {
                const int i2_80336_lb = max(1, 16 * ii2 - i0 + i1), i2_94221_ub = min(N, 16 * ii2 - i0 + i1 + 15);
                for (register int i2 = i2_80336_lb; i2 <= i2_94221_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
            }
          }
        }
        const int ii1_89790_lb = max(0, -ii0 + 1), ii1_4562_ub = min(-ii0 + (T + N) / 16, (N + 8 * k) / 16);
        #pragma omp parallel for
        for (register int ii1 = ii1_89790_lb; ii1 < ii1_4562_ub; ii1 += 1) {
          if (ii1 >= 1) {
            const int ii2_57785_lb = 0, ii2_96428_ub = min(min(min((N - 2) / 16 + 1, k - 2 * ii1 + 3 * N / 16 - 2), -ii0 - ii1 + (T + 2 * N) / 16 - 1), -ii0 - k + (T + N + 8 * k + 7) / 16);
            for (register int ii2 = ii2_57785_lb; ii2 <= ii2_96428_ub; ii2 += 1) {
              if (k == 2) {
                const int i0_78287_lb = max(max(max(16 * ii0 + 8, -N + 16 * ii0 + 16 * ii1 + 15), -N + 16 * ii0 + 16 * ii2 + 8), -2 * N + 16 * ii0 + 16 * ii1 + 16 * ii2 + 15), i0_4652_ub = min(min(T - 1, 16 * ii0 + 22), N + 16 * ii0 - 16 * ii1 + 15);
                for (register int i0 = i0_78287_lb; i0 <= i0_4652_ub; i0 += 1) {
                  const int i1_88965_lb = max(max(-16 * ii0 + 16 * ii1 + i0 - 15, -N + 16 * ii0 + 16 * ii1 + 16 * ii2 - i0 + 15), 16 * ii0 + 16 * ii1 - i0 + 15), i1_24439_ub = min(N, 16 * ii1 + 7);
                  for (register int i1 = i1_88965_lb; i1 <= i1_24439_ub; i1 += 1) {
                    const int i2_96087_lb = max(1, 16 * ii0 + 16 * ii1 + 16 * ii2 - i0 - i1 + 15), i2_22977_ub = min(N, 16 * ii0 + 16 * ii1 + 16 * ii2 - i0 - i1 + 30);
                    for (register int i2 = i2_96087_lb; i2 <= i2_22977_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                  const int i1_91640_lb = 16 * ii1 + 8, i1_97593_ub = min(min(min(min(N, -16 * ii0 + 16 * ii1 + i0), N - 16 * ii0 + 16 * ii1 - 16 * ii2 + i0), -N + 16 * ii0 + 16 * ii1 + 16 * ii2 - i0 + 30), 16 * ii0 + 16 * ii1 - i0 + 30);
                  for (register int i1 = i1_91640_lb; i1 <= i1_97593_ub; i1 += 1) {
                    const int i2_3639_lb = 16 * ii0 - 16 * ii1 + 16 * ii2 - i0 + i1, i2_87448_ub = N;
                    for (register int i2 = i2_3639_lb; i2 <= i2_87448_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                  const int i1_81553_lb = max(16 * ii1 + 8, -N + 16 * ii0 + 16 * ii1 + 16 * ii2 - i0 + 31), i1_50968_ub = min(min(N, -16 * ii0 + 16 * ii1 + i0), 16 * ii0 + 16 * ii1 - i0 + 30);
                  for (register int i1 = i1_81553_lb; i1 <= i1_50968_ub; i1 += 1) {
                    const int i2_83662_lb = max(1, 16 * ii0 - 16 * ii1 + 16 * ii2 - i0 + i1), i2_57522_ub = min(N, 16 * ii0 - 16 * ii1 + 16 * ii2 - i0 + i1 + 15);
                    for (register int i2 = i2_83662_lb; i2 <= i2_57522_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                }
              } else {
                const int i0_72242_lb = max(max(max(16 * ii0, -N + 16 * ii0 + 16 * ii1 + 15), -N + 16 * ii0 + 16 * ii2), -2 * N + 16 * ii0 + 16 * ii1 + 16 * ii2 + 15), i0_71705_ub = min(min(T - 1, 16 * ii0 + 14), N + 16 * ii0 - 16 * ii1 - 1);
                for (register int i0 = i0_72242_lb; i0 <= i0_71705_ub; i0 += 1) {
                  const int i1_63397_lb = max(max(-16 * ii0 + 16 * ii1 + i0 + 1, -N + 16 * ii0 + 16 * ii1 + 16 * ii2 - i0 + 15), 16 * ii0 + 16 * ii1 - i0 + 15), i1_66189_ub = min(N, 16 * ii1 + 15);
                  for (register int i1 = i1_63397_lb; i1 <= i1_66189_ub; i1 += 1) {
                    const int i2_704_lb = max(1, 16 * ii0 + 16 * ii1 + 16 * ii2 - i0 - i1 + 15), i2_76411_ub = min(N, 16 * ii0 + 16 * ii1 + 16 * ii2 - i0 - i1 + 30);
                    for (register int i2 = i2_704_lb; i2 <= i2_76411_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                  const int i1_9323_lb = 16 * ii1 + 16, i1_68887_ub = min(min(min(min(N, -16 * ii0 + 16 * ii1 + i0 + 16), N - 16 * ii0 + 16 * ii1 - 16 * ii2 + i0 + 16), -N + 16 * ii0 + 16 * ii1 + 16 * ii2 - i0 + 30), 16 * ii0 + 16 * ii1 - i0 + 30);
                  for (register int i1 = i1_9323_lb; i1 <= i1_68887_ub; i1 += 1) {
                    const int i2_42896_lb = 16 * ii0 - 16 * ii1 + 16 * ii2 - i0 + i1 - 16, i2_91389_ub = N;
                    for (register int i2 = i2_42896_lb; i2 <= i2_91389_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                  const int i1_45041_lb = max(16 * ii1 + 16, -N + 16 * ii0 + 16 * ii1 + 16 * ii2 - i0 + 31), i1_23232_ub = min(min(N, -16 * ii0 + 16 * ii1 + i0 + 16), 16 * ii0 + 16 * ii1 - i0 + 30);
                  for (register int i1 = i1_45041_lb; i1 <= i1_23232_ub; i1 += 1) {
                    const int i2_1962_lb = max(1, 16 * ii0 - 16 * ii1 + 16 * ii2 - i0 + i1 - 16), i2_34831_ub = min(N, 16 * ii0 - 16 * ii1 + 16 * ii2 - i0 + i1 - 1);
                    for (register int i2 = i2_1962_lb; i2 <= i2_34831_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                }
              }
            }
          } else if (k == 1) {
            const int ii2_44146_lb = 0, ii2_59747_ub = min((N - 1) / 16, -2 * ii0 + (2 * T + N - 1) / 16 - 1);
            for (register int ii2 = ii2_44146_lb; ii2 <= ii2_59747_ub; ii2 += 1) {
              const int i0_31260_lb = max(16 * ii0, -N + 16 * ii0 + 15), i0_22433_ub = min(min(T - 1, N + 16 * ii0 - 1), 16 * ii0 + 14);
              for (register int i0 = i0_31260_lb; i0 <= i0_22433_ub; i0 += 1) {
                const int i1_80752_lb = max(-16 * ii0 + i0 + 1, 16 * ii0 - i0 + 15), i1_36577_ub = min(min(15, N), N - 16 * ii0 - 16 * ii2 + i0);
                for (register int i1 = i1_80752_lb; i1 <= i1_36577_ub; i1 += 1) {
                  const int i2_63225_lb = max(1, 16 * ii0 + 16 * ii2 - i0 - i1 + 15), i2_76839_ub = min(N, 16 * ii0 + 16 * ii2 - i0 - i1 + 30);
                  for (register int i2 = i2_63225_lb; i2 <= i2_76839_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
                const int i1_59554_lb = max(N - 16 * ii0 - 16 * ii2 + i0 + 1, 16 * ii0 - i0 + 15), i1_71217_ub = 15;
                for (register int i1 = i1_59554_lb; i1 <= i1_71217_ub; i1 += 1) {
                  const int i2_90784_lb = 16 * ii0 + 16 * ii2 - i0 - i1 + 15, i2_63194_ub = N;
                  for (register int i2 = i2_90784_lb; i2 <= i2_63194_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
                const int i1_58665_lb = 16, i1_88689_ub = min(min(N, -16 * ii0 + i0 + 16), 16 * ii0 - i0 + 30);
                for (register int i1 = i1_58665_lb; i1 <= i1_88689_ub; i1 += 1) {
                  const int i2_30514_lb = max(1, 16 * ii0 + 16 * ii2 - i0 + i1 - 16), i2_42327_ub = min(N, 16 * ii0 + 16 * ii2 - i0 + i1 - 1);
                  for (register int i2 = i2_30514_lb; i2 <= i2_42327_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
              }
            }
            if (((N - 1) % 16) + 2 * T >= 32 * ii0 + 32 && (N - 1) % 16 >= 1) {
              const int i0_46212_lb = max(-((N - 2) % 16) + 16 * ii0 + 14, -((N - 2) % 16) - N + 16 * ii0 + 29), i0_2756_ub = min(min(T - 1, N + 16 * ii0 - 1), 16 * ii0 + 14);
              for (register int i0 = i0_46212_lb; i0 <= i0_2756_ub; i0 += 1) {
                const int i1_30384_lb = max(-16 * ii0 + i0 + 1, -((N - 2) % 16) + 16 * ii0 - i0 + 29), i1_25961_ub = min(min(N, 16 * ii0 - i0 + 30), ((N - 2) % 16) - 16 * ii0 + i0 + 2);
                for (register int i1 = i1_30384_lb; i1 <= i1_25961_ub; i1 += 1) {
                  const int i2_85297_lb = max(-((N - 2) % 16) + N + 16 * ii0 - i0 + i1 - 2, -((N - 2) % 16) + N + 16 * ii0 - i0 - i1 + 29), i2_31089_ub = N;
                  for (register int i2 = i2_85297_lb; i2 <= i2_31089_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
              }
            }
            const int i0_18724_lb = max(max(16 * ii0, -((2 * T + N + 15) % 16) + 2 * T - 16 * ii0 - 1), -((2 * T + N + 15) % 16) + 2 * T - N - 16 * ii0 + 14), i0_94621_ub = T;
            for (register int i0 = i0_18724_lb; i0 < i0_94621_ub; i0 += 1) {
              const int i1_16328_lb = max(16 * ii0 - i0 + 15, -((2 * T + N + 15) % 16) + 2 * T - 16 * ii0 - i0 + 14), i1_77973_ub = min(min(N, -16 * ii0 + i0 + 16), ((2 * T + N + 15) % 16) - 2 * T + 16 * ii0 + i0 + 17);
              for (register int i1 = i1_16328_lb; i1 <= i1_77973_ub; i1 += 1) {
                const int i2_2362_lb = max(-((2 * T + N + 15) % 16) + 2 * T + N - 16 * ii0 - i0 + i1 - 17, -((2 * T + N + 15) % 16) + 2 * T + N - 16 * ii0 - i0 - i1 + 14), i2_61370_ub = N;
                for (register int i2 = i2_2362_lb; i2 <= i2_61370_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
            }
          } else if (T >= 16 * ii0 + 10) {
            const int i0_1205_lb = max(16 * ii0 + 8, -N + 16 * ii0 + 15), i0_4324_ub = min(min(T - 1, N + 16 * ii0 + 14), 16 * ii0 + 22);
            for (register int i0 = i0_1205_lb; i0 <= i0_4324_ub; i0 += 1) {
              const int i1_96201_lb = 1, i1_61704_ub = -16 * ii0 + i0 - 14;
              for (register int i1 = i1_96201_lb; i1 < i1_61704_ub; i1 += 1) {
                const int i2_80424_lb = 1, i2_43813_ub = min(N, 16 * ii0 - i0 - i1 + 30);
                for (register int i2 = i2_80424_lb; i2 <= i2_43813_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
              const int i1_489_lb = max(-16 * ii0 + i0 - 14, 16 * ii0 - i0 + 15), i1_77528_ub = min(min(N, -16 * ii0 + i0), 16 * ii0 - i0 + 30);
              for (register int i1 = i1_489_lb; i1 <= i1_77528_ub; i1 += 1) {
                const int i2_80390_lb = 1, i2_63714_ub = min(N, 16 * ii0 - i0 + i1 + 15);
                for (register int i2 = i2_80390_lb; i2 <= i2_63714_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
                const int i2_70719_lb = 16 * ii0 - i0 + i1 + 16, i2_56297_ub = min(N, 16 * ii0 - i0 - i1 + 30);
                for (register int i2 = i2_70719_lb; i2 <= i2_56297_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
            }
            const int i0_34931_lb = N + 16 * ii0 + 15, i0_61503_ub = min(T - 1, 16 * ii0 + 22);
            for (register int i0 = i0_34931_lb; i0 <= i0_61503_ub; i0 += 1) {
              for (register int i1 = 1; i1 <= N; i1 += 1) {
                const int i2_50192_lb = 1, i2_66357_ub = min(N, 16 * ii0 - i0 - i1 + 30);
                for (register int i2 = i2_50192_lb; i2 <= i2_66357_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
            }
            const int i0_52276_lb = 16 * ii0 + 23, i0_12756_ub = min(T - 1, 16 * ii0 + 28);
            for (register int i0 = i0_52276_lb; i0 <= i0_12756_ub; i0 += 1) {
              const int i1_85466_lb = 1, i1_82660_ub = min(N, 16 * ii0 - i0 + 29);
              for (register int i1 = i1_85466_lb; i1 <= i1_82660_ub; i1 += 1) {
                const int i2_55069_lb = 1, i2_87115_ub = min(N, 16 * ii0 - i0 - i1 + 30);
                for (register int i2 = i2_55069_lb; i2 <= i2_87115_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
            }
            const int ii2_13749_lb = 1, ii2_73794_ub = min(min((N + 15) / 16, -ii0 + (T + 2 * N) / 16 - 1), -ii0 + (T + N + 7) / 16 - 1);
            for (register int ii2 = ii2_13749_lb; ii2 <= ii2_73794_ub; ii2 += 1) {
              const int i0_81736_lb = max(max(16 * ii0 + 8, -N + 16 * ii0 + 16 * ii2 + 8), -2 * N + 16 * ii0 + 16 * ii2 + 15), i0_30078_ub = min(T - 1, 16 * ii0 + 29);
              for (register int i0 = i0_81736_lb; i0 <= i0_30078_ub; i0 += 1) {
                const int i1_51767_lb = max(max(1, -N + 16 * ii0 + 16 * ii2 - i0 + 15), 16 * ii0 - i0 + 15), i1_84098_ub = min(min(min(N, -16 * ii0 + i0), N - 16 * ii0 - 16 * ii2 + i0), 16 * ii0 - i0 + 30);
                for (register int i1 = i1_51767_lb; i1 <= i1_84098_ub; i1 += 1) {
                  const int i2_7800_lb = max(16 * ii0 + 16 * ii2 - i0 + i1, 16 * ii0 + 16 * ii2 - i0 - i1 + 15), i2_69324_ub = min(N, 16 * ii2 - 1);
                  for (register int i2 = i2_7800_lb; i2 <= i2_69324_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  const int i2_4775_lb = 16 * ii2, i2_4001_ub = min(N, 16 * ii0 + 16 * ii2 - i0 + i1 + 15);
                  for (register int i2 = i2_4775_lb; i2 <= i2_4001_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  const int i2_47380_lb = max(16 * ii2, 16 * ii0 + 16 * ii2 - i0 + i1 + 16), i2_1551_ub = min(N, 16 * ii0 + 16 * ii2 - i0 - i1 + 30);
                  for (register int i2 = i2_47380_lb; i2 <= i2_1551_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
              }
            }
          } else {
            const int ii2_64167_lb = 0, ii2_47870_ub = N / 16;
            for (register int ii2 = ii2_64167_lb; ii2 <= ii2_47870_ub; ii2 += 1) {
              const int i1_79079_lb = 7, i1_60909_ub = min(8, N);
              for (register int i1 = i1_79079_lb; i1 <= i1_60909_ub; i1 += 1) {
                const int i2_27936_lb = max(1, 16 * ii2), i2_49798_ub = min(N, 16 * ii2 + 15);
                for (register int i2 = i2_27936_lb; i2 <= i2_49798_ub; i2 += 1) {
                  S1(T - 1, i1, i2);
                }
              }
            }
          }
        }
      }
    }
  }
  if (T <= 8 && N >= 1) {
    const int ii1_33558_lb = 0, ii1_62868_ub = floord(N - 1, 16);
    #pragma omp parallel for
    for (register int ii1 = ii1_33558_lb; ii1 < ii1_62868_ub; ii1 += 1) {
      const int ii2_11301_lb = 0, ii2_69401_ub = min((T + N - 1) / 16, (N + 6) / 16);
      for (register int ii2 = ii2_11301_lb; ii2 <= ii2_69401_ub; ii2 += 1) {
        const int i0_72816_lb = max(0, -N + 16 * ii2), i0_77845_ub = min(min(6, T - 1), N - 16 * ii1 - 17);
        for (register int i0 = i0_72816_lb; i0 <= i0_77845_ub; i0 += 1) {
          const int i1_52111_lb = 16 * ii1 + i0 + 17, i1_25092_ub = min(N, 16 * ii1 - i0 + 30);
          for (register int i1 = i1_52111_lb; i1 <= i1_25092_ub; i1 += 1) {
            const int i2_90602_lb = max(1, 16 * ii2 - i0), i2_53929_ub = min(N, 16 * ii2 - i0 + 15);
            for (register int i2 = i2_90602_lb; i2 <= i2_53929_ub; i2 += 1) {
              S1(i0, i1, i2);
            }
          }
        }
      }
    }
    const int k_24105_lb = 1, k_62023_ub = min(2, T);
    for (register int k = k_24105_lb; k <= k_62023_ub; k += 1) {
      const int ii2_41044_lb = 0, ii2_37854_ub = floord(T + N - k, 16);
      for (register int ii2 = ii2_41044_lb; ii2 <= ii2_37854_ub; ii2 += 1) {
        const int i0_52169_lb = max(k - 1, -N + k + 16 * ii2 - 1), i0_22781_ub = min(T, N + 7 * k - 7);
        for (register int i0 = i0_52169_lb; i0 < i0_22781_ub; i0 += 1) {
          const int i1_84284_lb = max(1, -7 * k + i0 + 8), i1_20288_ub = min(min(N, -16 * k + i0 + 32), N - 16 * k - 16 * ii2 + i0 + 32);
          for (register int i1 = i1_84284_lb; i1 <= i1_20288_ub; i1 += 1) {
            const int i2_23231_lb = max(max(1, 16 * ii2 - i0), 16 * k + 16 * ii2 - i0 + i1 - 32), i2_92084_ub = min(N, 16 * k + 16 * ii2 - i0 + i1 - 17);
            for (register int i2 = i2_23231_lb; i2 <= i2_92084_ub; i2 += 1) {
              S1(i0, i1, i2);
            }
            if (k == 1) {
              const int i2_89613_lb = 16 * ii2 - i0 + i1, i2_44358_ub = min(N, 16 * ii2 - i0 + 15);
              for (register int i2 = i2_89613_lb; i2 <= i2_44358_ub; i2 += 1) {
                S1(i0, i1, i2);
              }
            }
          }
        }
      }
      if (k == 1) {
        const int ii1_12438_lb = 1, ii1_53345_ub = floord(T + N, 16);
        #pragma omp parallel for
        for (register int ii1 = ii1_12438_lb; ii1 < ii1_53345_ub; ii1 += 1) {
          const int ii2_45909_lb = 0, ii2_92957_ub = min((T + N - 1) / 16, -ii1 + (T + 2 * N) / 16 - 1);
          for (register int ii2 = ii2_45909_lb; ii2 <= ii2_92957_ub; ii2 += 1) {
            const int i0_1215_lb = max(max(max(0, -N + 16 * ii1 + 15), -N + 16 * ii2), -2 * N + 16 * ii1 + 16 * ii2 + 15), i0_41340_ub = T;
            for (register int i0 = i0_1215_lb; i0 < i0_41340_ub; i0 += 1) {
              const int i1_53866_lb = max(-N + 16 * ii1 + 16 * ii2 - i0 + 15, 16 * ii1 - i0 + 15), i1_29152_ub = min(min(N, 16 * ii1 + i0 + 16), N + 16 * ii1 - 16 * ii2 + i0 + 16);
              for (register int i1 = i1_53866_lb; i1 <= i1_29152_ub; i1 += 1) {
                const int i2_91138_lb = max(max(1, -16 * ii1 + 16 * ii2 - i0 + i1 - 16), 16 * ii1 + 16 * ii2 - i0 - i1 + 15), i2_3777_ub = min(N, 16 * ii1 + 16 * ii2 - i0 - i1 + 30);
                for (register int i2 = i2_91138_lb; i2 <= i2_3777_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
                const int i2_8372_lb = 16 * ii1 + 16 * ii2 - i0 - i1 + 31, i2_18791_ub = min(N, -16 * ii1 + 16 * ii2 - i0 + i1 - 1);
                for (register int i2 = i2_8372_lb; i2 <= i2_18791_ub; i2 += 1) {
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
