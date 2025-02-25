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
/* ./tc ../examples/pluto/3d7pt-t.scop.c --diamond-tiling --omp-for-codegen --omega-map-tc --debug --use-macros -b 16 --drop-bounds */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define S1_I(t,i,j,k) A[(t + 1)%2][i][j][k] = ((alpha * A[t%2][i][j][k]) + (beta * (((((A[t%2][i - 1][j][k] + A[t%2][i][j - 1][k]) + A[t%2][i][j][k - 1]) + A[t%2][i + 1][j][k]) + A[t%2][i][j + 1][k]) + A[t%2][i][j][k + 1])));
#define S1(t,i,j,k) S1_I((t),(i),(j),(k))
#pragma scop
if (T + N >= 33) {
  if (N <= 15) {
    for (register int k = 1; k <= 2; k += 1) {
      if (k == 2) {
        const int ii2_24661_lb = 0, ii2_60279_ub = min(-((-N + 19) / 16) + 2, (T + N - 3) / 16);
        for (register int ii2 = ii2_24661_lb; ii2 <= ii2_60279_ub; ii2 += 1) {
          const int ii3_34280_lb = max(0, ii2 - (N + 14) / 16), ii3_96311_ub = min(min(-((-N + 19) / 16) + 2, (T + N - 3) / 16), ii2 - (-N + 17) / 16 + 1);
          for (register int ii3 = ii3_34280_lb; ii3 <= ii3_96311_ub; ii3 += 1) {
            const int i0_50038_lb = max(max(1, -N + 16 * ii2 + 1), -N + 16 * ii3 + 1), i0_97259_ub = min(min(min(29, T - 2), 16 * ii2 + 15), 16 * ii3 + 15);
            for (register int i0 = i0_50038_lb; i0 <= i0_97259_ub; i0 += 1) {
              const int i1_72667_lb = 1, i1_71528_ub = min(min(N, i0), -i0 + 30);
              for (register int i1 = i1_72667_lb; i1 <= i1_71528_ub; i1 += 1) {
                const int i2_8094_lb = max(1, 16 * ii2 - i0 + 1), i2_18896_ub = min(N, 16 * ii2 - i0 + 16);
                for (register int i2 = i2_8094_lb; i2 <= i2_18896_ub; i2 += 1) {
                  const int i3_75938_lb = max(1, 16 * ii3 - i0 + 1), i3_61824_ub = min(N, 16 * ii3 - i0 + 16);
                  for (register int i3 = i3_75938_lb; i3 <= i3_61824_ub; i3 += 1) {
                    S1(i0, i1, i2, i3);
                  }
                }
              }
            }
          }
        }
      } else {
        const int ii2_79319_lb = 0, ii2_82539_ub = floord(N - 1, 8);
        for (register int ii2 = ii2_79319_lb; ii2 <= ii2_82539_ub; ii2 += 1) {
          const int ii3_63304_lb = 0, ii3_1392_ub = (N - 1) / 8;
          for (register int ii3 = ii3_63304_lb; ii3 <= ii3_1392_ub; ii3 += 1) {
            const int i0_28558_lb = max(max(0, -N + 16 * ii2), -N + 16 * ii3), i0_4007_ub = N;
            for (register int i0 = i0_28558_lb; i0 < i0_4007_ub; i0 += 1) {
              const int i1_40697_lb = i0 + 1, i1_36997_ub = N;
              for (register int i1 = i1_40697_lb; i1 <= i1_36997_ub; i1 += 1) {
                const int i2_67067_lb = max(1, 16 * ii2 - i0), i2_81092_ub = min(N, 16 * ii2 - i0 + 15);
                for (register int i2 = i2_67067_lb; i2 <= i2_81092_ub; i2 += 1) {
                  const int i3_16860_lb = max(1, 16 * ii3 - i0), i3_84951_ub = min(N, 16 * ii3 - i0 + 15);
                  for (register int i3 = i3_16860_lb; i3 <= i3_84951_ub; i3 += 1) {
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
  const int ii0_37120_lb = max(0, -((N + 14) / 15) + 2), ii0_3933_ub = min((T + 6) / 16, (T + N - 1) / 16);
  for (register int ii0 = ii0_37120_lb; ii0 < ii0_3933_ub; ii0 += 1) {
    if (ii0 == 0) {
      const int ii1_12982_lb = 0, ii1_94154_ub = (N - 1) / 16;
      #pragma omp parallel for
      for (register int ii1 = ii1_12982_lb; ii1 < ii1_94154_ub; ii1 += 1) {
        const int ii2_76147_lb = 0, ii2_56962_ub = (N + 6) / 16;
        for (register int ii2 = ii2_76147_lb; ii2 <= ii2_56962_ub; ii2 += 1) {
          const int ii3_48709_lb = 0, ii3_808_ub = (N + 6) / 16;
          for (register int ii3 = ii3_48709_lb; ii3 <= ii3_808_ub; ii3 += 1) {
            const int i0_17242_lb = max(max(0, -N + 16 * ii2), -N + 16 * ii3), i0_82989_ub = min(6, N - 16 * ii1 - 17);
            for (register int i0 = i0_17242_lb; i0 <= i0_82989_ub; i0 += 1) {
              const int i1_97120_lb = 16 * ii1 + i0 + 17, i1_67280_ub = min(N, 16 * ii1 - i0 + 30);
              for (register int i1 = i1_97120_lb; i1 <= i1_67280_ub; i1 += 1) {
                const int i2_96600_lb = max(1, 16 * ii2 - i0), i2_86139_ub = min(N, 16 * ii2 - i0 + 15);
                for (register int i2 = i2_96600_lb; i2 <= i2_86139_ub; i2 += 1) {
                  const int i3_55160_lb = max(1, 16 * ii3 - i0), i3_4695_ub = min(N, 16 * ii3 - i0 + 15);
                  for (register int i3 = i3_55160_lb; i3 <= i3_4695_ub; i3 += 1) {
                    S1(i0, i1, i2, i3);
                  }
                }
              }
            }
          }
        }
      }
      if (N <= 23) {
        const int ii2_5035_lb = 0, ii2_47450_ub = min((N - 2) / 16 + 1, (T + N - 2) / 16);
        for (register int ii2 = ii2_5035_lb; ii2 <= ii2_47450_ub; ii2 += 1) {
          const int ii3_66519_lb = 0, ii3_707_ub = min((N - 2) / 16 + 1, (T + N - 2) / 16);
          for (register int ii3 = ii3_66519_lb; ii3 <= ii3_707_ub; ii3 += 1) {
            const int i0_46341_lb = max(max(0, -N + 16 * ii2), -N + 16 * ii3), i0_46175_ub = min(14, T - 2);
            for (register int i0 = i0_46341_lb; i0 <= i0_46175_ub; i0 += 1) {
              const int i1_2099_lb = i0 + 1, i1_74899_ub = min(min(N, i0 + 16), -i0 + 30);
              for (register int i1 = i1_2099_lb; i1 <= i1_74899_ub; i1 += 1) {
                const int i2_50183_lb = max(1, 16 * ii2 - i0), i2_42796_ub = min(N, 16 * ii2 - i0 + 15);
                for (register int i2 = i2_50183_lb; i2 <= i2_42796_ub; i2 += 1) {
                  const int i3_28249_lb = max(1, 16 * ii3 - i0), i3_33602_ub = min(N, 16 * ii3 - i0 + 15);
                  for (register int i3 = i3_28249_lb; i3 <= i3_33602_ub; i3 += 1) {
                    S1(i0, i1, i2, i3);
                  }
                }
              }
            }
          }
        }
      }
    }
    if (N >= 8 && N + 16 * ii0 >= 24) {
      if (ii0 == 0) {
        const int ii2_40241_lb = 0, ii2_61461_ub = min((N - 2) / 16 + 1, (T + N - 2) / 16);
        for (register int ii2 = ii2_40241_lb; ii2 <= ii2_61461_ub; ii2 += 1) {
          const int ii3_18553_lb = 0, ii3_93713_ub = min((N - 2) / 16 + 1, (T + N - 2) / 16);
          for (register int ii3 = ii3_18553_lb; ii3 <= ii3_93713_ub; ii3 += 1) {
            const int i0_81746_lb = max(max(0, -N + 16 * ii2), -N + 16 * ii3), i0_31535_ub = min(14, T - 2);
            for (register int i0 = i0_81746_lb; i0 <= i0_31535_ub; i0 += 1) {
              const int i1_4219_lb = i0 + 1, i1_57893_ub = min(i0 + 16, -i0 + 30);
              for (register int i1 = i1_4219_lb; i1 <= i1_57893_ub; i1 += 1) {
                const int i2_4850_lb = max(1, 16 * ii2 - i0), i2_52928_ub = min(N, 16 * ii2 - i0 + 15);
                for (register int i2 = i2_4850_lb; i2 <= i2_52928_ub; i2 += 1) {
                  const int i3_75054_lb = max(1, 16 * ii3 - i0), i3_22092_ub = min(N, 16 * ii3 - i0 + 15);
                  for (register int i3 = i3_75054_lb; i3 <= i3_22092_ub; i3 += 1) {
                    S1(i0, i1, i2, i3);
                  }
                }
              }
            }
          }
        }
      }
      const int ii1_52270_lb = max(0, -ii0 + 1), ii1_88526_ub = (N + 8) / 16;
      #pragma omp parallel for
      for (register int ii1 = ii1_52270_lb; ii1 < ii1_88526_ub; ii1 += 1) {
        const int ii2_5724_lb = 0, ii2_48870_ub = min(min(-ii1 + (N - 1) / 8, (N - 2) / 16 + 1), -ii0 + (T + N - 2) / 16);
        for (register int ii2 = ii2_5724_lb; ii2 <= ii2_48870_ub; ii2 += 1) {
          const int ii3_74665_lb = 0, ii3_60884_ub = min(min(-ii1 + (N - 1) / 8, (N - 2) / 16 + 1), -ii0 + (T + N - 2) / 16);
          for (register int ii3 = ii3_74665_lb; ii3 <= ii3_60884_ub; ii3 += 1) {
            const int i0_69917_lb = max(max(max(16 * ii0, -N + 16 * ii0 + 16 * ii1 + 15), -N + 16 * ii0 + 16 * ii2), -N + 16 * ii0 + 16 * ii3), i0_96053_ub = min(min(T - 2, 16 * ii0 + 14), N + 16 * ii0 - 16 * ii1 - 1);
            for (register int i0 = i0_69917_lb; i0 <= i0_96053_ub; i0 += 1) {
              const int i1_24686_lb = max(-16 * ii0 + 16 * ii1 + i0 + 1, 16 * ii0 + 16 * ii1 - i0 + 15), i1_36437_ub = min(min(N, -16 * ii0 + 16 * ii1 + i0 + 16), 16 * ii0 + 16 * ii1 - i0 + 30);
              for (register int i1 = i1_24686_lb; i1 <= i1_36437_ub; i1 += 1) {
                const int i2_96760_lb = max(1, 16 * ii0 + 16 * ii2 - i0), i2_71028_ub = min(N, 16 * ii0 + 16 * ii2 - i0 + 15);
                for (register int i2 = i2_96760_lb; i2 <= i2_71028_ub; i2 += 1) {
                  const int i3_82612_lb = max(1, 16 * ii0 + 16 * ii3 - i0), i3_15211_ub = min(N, 16 * ii0 + 16 * ii3 - i0 + 15);
                  for (register int i3 = i3_82612_lb; i3 <= i3_15211_ub; i3 += 1) {
                    S1(i0, i1, i2, i3);
                  }
                }
              }
            }
          }
        }
      }
      if (N >= 16 && 16 * ii0 + 32 >= T + N) {
        for (register int ii2 = 0; ii2 < 2; ii2 += 1) {
          for (register int ii3 = 0; ii3 < 2; ii3 += 1) {
            const int i0_82749_lb = 16 * ii0 + 8, i0_14600_ub = T - 1;
            for (register int i0 = i0_82749_lb; i0 < i0_14600_ub; i0 += 1) {
              const int i1_68341_lb = 16 * ii0 - i0 + 15, i1_1302_ub = -16 * ii0 + i0;
              for (register int i1 = i1_68341_lb; i1 <= i1_1302_ub; i1 += 1) {
                const int i2_24665_lb = max(1, 16 * ii0 + 16 * ii2 - i0 + 8), i2_50088_ub = min(N, 16 * ii0 + 16 * ii2 - i0 + 23);
                for (register int i2 = i2_24665_lb; i2 <= i2_50088_ub; i2 += 1) {
                  const int i3_49190_lb = max(1, 16 * ii0 + 16 * ii3 - i0 + 8), i3_45237_ub = min(N, 16 * ii0 + 16 * ii3 - i0 + 23);
                  for (register int i3 = i3_49190_lb; i3 <= i3_45237_ub; i3 += 1) {
                    S1(i0, i1, i2, i3);
                  }
                }
              }
            }
          }
        }
      }
    }
    if (N <= 15) {
      const int ii2_7981_lb = 0, ii2_70392_ub = min((N + 5) / 16 + 1, -ii0 + (T + N + 6) / 16 - 1);
      for (register int ii2 = ii2_7981_lb; ii2 <= ii2_70392_ub; ii2 += 1) {
        const int ii3_98165_lb = max(0, ii2 - (N + 14) / 16), ii3_83035_ub = min(min(ii2 - (-N + 17) / 16 + 1, (N + 5) / 16 + 1), -ii0 + (T + N + 6) / 16 - 1);
        for (register int ii3 = ii3_98165_lb; ii3 <= ii3_83035_ub; ii3 += 1) {
          const int i0_92484_lb = max(max(max(16 * ii0 + 8, -N + 16 * ii0 + 15), -N + 16 * ii0 + 16 * ii2 + 8), -N + 16 * ii0 + 16 * ii3 + 8), i0_66787_ub = min(min(min(T - 2, 16 * ii0 + 29), 16 * ii0 + 16 * ii2 + 22), 16 * ii0 + 16 * ii3 + 22);
          for (register int i0 = i0_92484_lb; i0 <= i0_66787_ub; i0 += 1) {
            const int i1_71561_lb = max(1, 16 * ii0 - i0 + 15), i1_14560_ub = min(min(N, -16 * ii0 + i0), 16 * ii0 - i0 + 30);
            for (register int i1 = i1_71561_lb; i1 <= i1_14560_ub; i1 += 1) {
              const int i2_32010_lb = max(1, 16 * ii0 + 16 * ii2 - i0 + 8), i2_62579_ub = min(N, 16 * ii0 + 16 * ii2 - i0 + 23);
              for (register int i2 = i2_32010_lb; i2 <= i2_62579_ub; i2 += 1) {
                const int i3_91796_lb = max(1, 16 * ii0 + 16 * ii3 - i0 + 8), i3_1927_ub = min(N, 16 * ii0 + 16 * ii3 - i0 + 23);
                for (register int i3 = i3_91796_lb; i3 <= i3_1927_ub; i3 += 1) {
                  S1(i0, i1, i2, i3);
                }
              }
            }
          }
        }
      }
    } else if (T + N >= 16 * ii0 + 33) {
      const int ii1_58632_lb = 0, ii1_16482_ub = min(-ii0 + (T + N - 1) / 16 - 1, N / 16);
      #pragma omp parallel for
      for (register int ii1 = ii1_58632_lb; ii1 <= ii1_16482_ub; ii1 += 1) {
        if (ii1 >= 1) {
          const int ii2_38364_lb = 0, ii2_55392_ub = min(min(-ii1 + (N + 3) / 8, (N - 2) / 16 + 1), -ii0 + (T + N + 6) / 16 - 1);
          for (register int ii2 = ii2_38364_lb; ii2 <= ii2_55392_ub; ii2 += 1) {
            const int ii3_87510_lb = 0, ii3_37329_ub = min(min(-ii1 + (N + 3) / 8, (N - 2) / 16 + 1), -ii0 + (T + N + 6) / 16 - 1);
            for (register int ii3 = ii3_87510_lb; ii3 <= ii3_37329_ub; ii3 += 1) {
              const int i0_70603_lb = max(max(max(16 * ii0 + 8, -N + 16 * ii0 + 16 * ii1 + 15), -N + 16 * ii0 + 16 * ii2 + 8), -N + 16 * ii0 + 16 * ii3 + 8), i0_49790_ub = min(min(T - 2, 16 * ii0 + 22), N + 16 * ii0 - 16 * ii1 + 15);
              for (register int i0 = i0_70603_lb; i0 <= i0_49790_ub; i0 += 1) {
                const int i1_2828_lb = max(-16 * ii0 + 16 * ii1 + i0 - 15, 16 * ii0 + 16 * ii1 - i0 + 15), i1_44962_ub = min(min(N, -16 * ii0 + 16 * ii1 + i0), 16 * ii0 + 16 * ii1 - i0 + 30);
                for (register int i1 = i1_2828_lb; i1 <= i1_44962_ub; i1 += 1) {
                  const int i2_56670_lb = max(1, 16 * ii0 + 16 * ii2 - i0 + 8), i2_85578_ub = min(N, 16 * ii0 + 16 * ii2 - i0 + 23);
                  for (register int i2 = i2_56670_lb; i2 <= i2_85578_ub; i2 += 1) {
                    const int i3_75915_lb = max(1, 16 * ii0 + 16 * ii3 - i0 + 8), i3_25012_ub = min(N, 16 * ii0 + 16 * ii3 - i0 + 23);
                    for (register int i3 = i3_75915_lb; i3 <= i3_25012_ub; i3 += 1) {
                      S1(i0, i1, i2, i3);
                    }
                  }
                }
              }
            }
          }
        } else if (ii0 >= 1) {
          const int ii2_3232_lb = 0, ii2_16932_ub = min((N + 5) / 16 + 1, -ii0 + (T + N + 6) / 16 - 1);
          for (register int ii2 = ii2_3232_lb; ii2 <= ii2_16932_ub; ii2 += 1) {
            const int ii3_75100_lb = max(0, ii2 - (N + 14) / 16), ii3_68774_ub = min(min(ii2 + (N - 2) / 16 + 1, (N + 5) / 16 + 1), -ii0 + (T + N + 6) / 16 - 1);
            for (register int ii3 = ii3_75100_lb; ii3 <= ii3_68774_ub; ii3 += 1) {
              const int i0_62169_lb = max(max(16 * ii0 + 8, -N + 16 * ii0 + 16 * ii2 + 8), -N + 16 * ii0 + 16 * ii3 + 8), i0_99433_ub = min(min(min(T - 2, 16 * ii0 + 29), 16 * ii0 + 16 * ii2 + 22), 16 * ii0 + 16 * ii3 + 22);
              for (register int i0 = i0_62169_lb; i0 <= i0_99433_ub; i0 += 1) {
                const int i1_39166_lb = max(1, 16 * ii0 - i0 + 15), i1_76687_ub = min(-16 * ii0 + i0, 16 * ii0 - i0 + 30);
                for (register int i1 = i1_39166_lb; i1 <= i1_76687_ub; i1 += 1) {
                  const int i2_82469_lb = max(1, 16 * ii0 + 16 * ii2 - i0 + 8), i2_48002_ub = min(N, 16 * ii0 + 16 * ii2 - i0 + 23);
                  for (register int i2 = i2_82469_lb; i2 <= i2_48002_ub; i2 += 1) {
                    const int i3_43474_lb = max(1, 16 * ii0 + 16 * ii3 - i0 + 8), i3_70382_ub = min(N, 16 * ii0 + 16 * ii3 - i0 + 23);
                    for (register int i3 = i3_43474_lb; i3 <= i3_70382_ub; i3 += 1) {
                      S1(i0, i1, i2, i3);
                    }
                  }
                }
              }
            }
          }
        } else {
          const int ii2_62562_lb = 0, ii2_75484_ub = min((N - 4) / 16 + 2, (T + N - 3) / 16);
          for (register int ii2 = ii2_62562_lb; ii2 <= ii2_75484_ub; ii2 += 1) {
            const int ii3_32961_lb = max(0, ii2 - (N + 14) / 16), ii3_54358_ub = min(min((N - 4) / 16 + 2, (T + N - 3) / 16), ii2 + (N - 2) / 16 + 1);
            for (register int ii3 = ii3_32961_lb; ii3 <= ii3_54358_ub; ii3 += 1) {
              const int i0_77412_lb = max(max(1, -N + 16 * ii2 + 1), -N + 16 * ii3 + 1), i0_91593_ub = min(min(min(29, T - 2), 16 * ii2 + 15), 16 * ii3 + 15);
              for (register int i0 = i0_77412_lb; i0 <= i0_91593_ub; i0 += 1) {
                const int i1_70841_lb = 1, i1_32128_ub = min(i0, -i0 + 30);
                for (register int i1 = i1_70841_lb; i1 <= i1_32128_ub; i1 += 1) {
                  const int i2_63337_lb = max(1, 16 * ii2 - i0 + 1), i2_58351_ub = min(N, 16 * ii2 - i0 + 16);
                  for (register int i2 = i2_63337_lb; i2 <= i2_58351_ub; i2 += 1) {
                    const int i3_85809_lb = max(1, 16 * ii3 - i0 + 1), i3_50292_ub = min(N, 16 * ii3 - i0 + 16);
                    for (register int i3 = i3_85809_lb; i3 <= i3_50292_ub; i3 += 1) {
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
  if (((T + 6) % 16) + N >= 23 && (T + 6) % 16 >= 8) {
    if (T <= 9) {
      const int ii1_24493_lb = 0, ii1_88638_ub = (N - 1) / 16;
      #pragma omp parallel for
      for (register int ii1 = ii1_24493_lb; ii1 < ii1_88638_ub; ii1 += 1) {
        const int ii2_11607_lb = 0, ii2_81164_ub = min((T + N - 2) / 16, (N + 6) / 16);
        for (register int ii2 = ii2_11607_lb; ii2 <= ii2_81164_ub; ii2 += 1) {
          const int ii3_90568_lb = 0, ii3_87522_ub = min((T + N - 2) / 16, (N + 6) / 16);
          for (register int ii3 = ii3_90568_lb; ii3 <= ii3_87522_ub; ii3 += 1) {
            const int i0_6176_lb = max(max(0, -N + 16 * ii2), -N + 16 * ii3), i0_93800_ub = min(min(6, T - 2), N - 16 * ii1 - 17);
            for (register int i0 = i0_6176_lb; i0 <= i0_93800_ub; i0 += 1) {
              const int i1_4454_lb = 16 * ii1 + i0 + 17, i1_81276_ub = min(N, 16 * ii1 - i0 + 30);
              for (register int i1 = i1_4454_lb; i1 <= i1_81276_ub; i1 += 1) {
                const int i2_62575_lb = max(1, 16 * ii2 - i0), i2_82976_ub = min(N, 16 * ii2 - i0 + 15);
                for (register int i2 = i2_62575_lb; i2 <= i2_82976_ub; i2 += 1) {
                  const int i3_80709_lb = max(1, 16 * ii3 - i0), i3_18093_ub = min(N, 16 * ii3 - i0 + 15);
                  for (register int i3 = i3_80709_lb; i3 <= i3_18093_ub; i3 += 1) {
                    S1(i0, i1, i2, i3);
                  }
                }
              }
            }
          }
        }
      }
    }
    const int ii1_59663_lb = 0, ii1_79530_ub = -((T - 2) / 16) + (T + N - 1) / 16;
    #pragma omp parallel for
    for (register int ii1 = ii1_59663_lb; ii1 < ii1_79530_ub; ii1 += 1) {
      const int ii2_66096_lb = 0, ii2_19489_ub = -((T - 2) / 16) + (T + N - 2) / 16;
      for (register int ii2 = ii2_66096_lb; ii2 <= ii2_19489_ub; ii2 += 1) {
        const int ii3_49913_lb = 0, ii3_28658_ub = -((T - 2) / 16) + (T + N - 2) / 16;
        for (register int ii3 = ii3_49913_lb; ii3 <= ii3_28658_ub; ii3 += 1) {
          const int i0_94974_lb = max(max(max(max(-((T + 14) % 16) + T - N + 16 * ii3 - 2, -((T + 14) % 16) + T - N + 16 * ii2 - 2), -((T + 14) % 16) + T - 2), -8 * ii1 + 8 * ((T - 2) / 16) + 8), -((T + 14) % 16) + T - N + 16 * ii1 + 13), i0_99226_ub = T - 1;
          for (register int i0 = i0_94974_lb; i0 < i0_99226_ub; i0 += 1) {
            const int i1_99369_lb = -((T + 14) % 16) + T + 16 * ii1 - i0 + 13, i1_88738_ub = min(N, ((T + 14) % 16) - T + 16 * ii1 + i0 + 18);
            for (register int i1 = i1_99369_lb; i1 <= i1_88738_ub; i1 += 1) {
              const int i2_7172_lb = max(1, -((T + 14) % 16) + T + 16 * ii2 - i0 - 2), i2_70210_ub = min(N, -((T + 14) % 16) + T + 16 * ii2 - i0 + 13);
              for (register int i2 = i2_7172_lb; i2 <= i2_70210_ub; i2 += 1) {
                const int i3_37218_lb = max(1, -((T + 14) % 16) + T + 16 * ii3 - i0 - 2), i3_70509_ub = min(N, -((T + 14) % 16) + T + 16 * ii3 - i0 + 13);
                for (register int i3 = i3_37218_lb; i3 <= i3_70509_ub; i3 += 1) {
                  S1(i0, i1, i2, i3);
                }
              }
            }
          }
          if (T <= 9 && ii1 == 0) {
            const int i0_44913_lb = max(max(0, -N + 16 * ii2), -N + 16 * ii3), i0_39380_ub = T - 1;
            for (register int i0 = i0_44913_lb; i0 < i0_39380_ub; i0 += 1) {
              const int i1_37154_lb = i0 + 1, i1_69407_ub = i0 + 16;
              for (register int i1 = i1_37154_lb; i1 <= i1_69407_ub; i1 += 1) {
                const int i2_28018_lb = max(1, 16 * ii2 - i0), i2_48761_ub = min(N, 16 * ii2 - i0 + 15);
                for (register int i2 = i2_28018_lb; i2 <= i2_48761_ub; i2 += 1) {
                  const int i3_50571_lb = max(1, 16 * ii3 - i0), i3_18586_ub = min(N, 16 * ii3 - i0 + 15);
                  for (register int i3 = i3_50571_lb; i3 <= i3_18586_ub; i3 += 1) {
                    S1(i0, i1, i2, i3);
                  }
                }
              }
            }
          }
        }
      }
    }
    if (T <= 9) {
      const int ii2_36283_lb = 0, ii2_56747_ub = (T + N - 3) / 16;
      for (register int ii2 = ii2_36283_lb; ii2 <= ii2_56747_ub; ii2 += 1) {
        const int ii3_28738_lb = 0, ii3_57089_ub = (T + N - 3) / 16;
        for (register int ii3 = ii3_28738_lb; ii3 <= ii3_57089_ub; ii3 += 1) {
          const int i0_54375_lb = max(max(1, -N + 16 * ii2 + 1), -N + 16 * ii3 + 1), i0_7665_ub = T - 1;
          for (register int i0 = i0_54375_lb; i0 < i0_7665_ub; i0 += 1) {
            for (register int i1 = 1; i1 <= i0; i1 += 1) {
              const int i2_25759_lb = max(1, 16 * ii2 - i0 + 1), i2_16080_ub = min(N, 16 * ii2 - i0 + 16);
              for (register int i2 = i2_25759_lb; i2 <= i2_16080_ub; i2 += 1) {
                const int i3_30967_lb = max(1, 16 * ii3 - i0 + 1), i3_8207_ub = min(N, 16 * ii3 - i0 + 16);
                for (register int i3 = i3_30967_lb; i3 <= i3_8207_ub; i3 += 1) {
                  S1(i0, i1, i2, i3);
                }
              }
            }
          }
        }
      }
    }
  }
} else {
  const int ii2_35570_lb = 0, ii2_97232_ub = floord(T + N - 2, 16);
  for (register int ii2 = ii2_35570_lb; ii2 <= ii2_97232_ub; ii2 += 1) {
    const int ii3_36865_lb = 0, ii3_46896_ub = (T + N - 2) / 16;
    for (register int ii3 = ii3_36865_lb; ii3 <= ii3_46896_ub; ii3 += 1) {
      const int i0_96458_lb = 0, i0_52586_ub = min(T - 1, N - 16);
      for (register int i0 = i0_96458_lb; i0 < i0_52586_ub; i0 += 1) {
        const int i1_51986_lb = i0 + 17, i1_19982_ub = N;
        for (register int i1 = i1_51986_lb; i1 <= i1_19982_ub; i1 += 1) {
          const int i2_22796_lb = max(1, 16 * ii2 - i0), i2_5556_ub = min(N, 16 * ii2 - i0 + 15);
          for (register int i2 = i2_22796_lb; i2 <= i2_5556_ub; i2 += 1) {
            const int i3_6844_lb = max(1, 16 * ii3 - i0), i3_67710_ub = min(N, 16 * ii3 - i0 + 15);
            for (register int i3 = i3_6844_lb; i3 <= i3_67710_ub; i3 += 1) {
              S1(i0, i1, i2, i3);
            }
          }
        }
      }
    }
  }
  const int k_44936_lb = 1, k_43998_ub = min(2, T - 1);
  for (register int k = k_44936_lb; k <= k_43998_ub; k += 1) {
    if (k == 2) {
      const int ii2_53469_lb = 0, ii2_89306_ub = floord(T + N - 3, 16);
      for (register int ii2 = ii2_53469_lb; ii2 <= ii2_89306_ub; ii2 += 1) {
        const int ii3_92759_lb = max(0, ii2 - (N + 14) / 16), ii3_4040_ub = min((T + N - 3) / 16, ii2 + (N + 14) / 16);
        for (register int ii3 = ii3_92759_lb; ii3 <= ii3_4040_ub; ii3 += 1) {
          const int i0_24244_lb = max(max(1, -N + 16 * ii2 + 1), -N + 16 * ii3 + 1), i0_45394_ub = min(min(T - 2, 16 * ii2 + 15), 16 * ii3 + 15);
          for (register int i0 = i0_24244_lb; i0 <= i0_45394_ub; i0 += 1) {
            const int i1_77139_lb = 1, i1_52983_ub = min(N, i0);
            for (register int i1 = i1_77139_lb; i1 <= i1_52983_ub; i1 += 1) {
              const int i2_18835_lb = max(1, 16 * ii2 - i0 + 1), i2_31514_ub = min(N, 16 * ii2 - i0 + 16);
              for (register int i2 = i2_18835_lb; i2 <= i2_31514_ub; i2 += 1) {
                const int i3_77000_lb = max(1, 16 * ii3 - i0 + 1), i3_75253_ub = min(N, 16 * ii3 - i0 + 16);
                for (register int i3 = i3_77000_lb; i3 <= i3_75253_ub; i3 += 1) {
                  S1(i0, i1, i2, i3);
                }
              }
            }
          }
        }
      }
    } else {
      const int ii2_82950_lb = 0, ii2_19111_ub = min(floord(N - 1, 8), floord(T + N - 2, 16));
      for (register int ii2 = ii2_82950_lb; ii2 <= ii2_19111_ub; ii2 += 1) {
        const int ii3_7685_lb = 0, ii3_13917_ub = min((N - 1) / 8, (T + N - 2) / 16);
        for (register int ii3 = ii3_7685_lb; ii3 <= ii3_13917_ub; ii3 += 1) {
          const int i0_27318_lb = max(max(0, -N + 16 * ii2), -N + 16 * ii3), i0_59607_ub = min(T - 1, N);
          for (register int i0 = i0_27318_lb; i0 < i0_59607_ub; i0 += 1) {
            const int i1_11149_lb = i0 + 1, i1_80536_ub = min(N, i0 + 16);
            for (register int i1 = i1_11149_lb; i1 <= i1_80536_ub; i1 += 1) {
              const int i2_6503_lb = max(1, 16 * ii2 - i0), i2_23960_ub = min(N, 16 * ii2 - i0 + 15);
              for (register int i2 = i2_6503_lb; i2 <= i2_23960_ub; i2 += 1) {
                const int i3_33122_lb = max(1, 16 * ii3 - i0), i3_74841_ub = min(N, 16 * ii3 - i0 + 15);
                for (register int i3 = i3_33122_lb; i3 <= i3_74841_ub; i3 += 1) {
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
