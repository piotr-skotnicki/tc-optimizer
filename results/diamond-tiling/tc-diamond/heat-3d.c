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
/* ./tc ../examples/pluto/heat-3d-t.scop.c --diamond-tiling --omp-for-codegen --omega-map-tc --debug -b 8 --drop-bounds --use-macros */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define S1_I(t,i,j,k) A[(t + 1)%2][i][j][k] = ((((0.125 * ((A[t%2][i + 1][j][k] - (2.0 * A[t%2][i][j][k])) + A[t%2][i - 1][j][k])) + (0.125 * ((A[t%2][i][j + 1][k] - (2.0 * A[t%2][i][j][k])) + A[t%2][i][j - 1][k]))) + (0.125 * ((A[t%2][i][j][k - 1] - (2.0 * A[t%2][i][j][k])) + A[t%2][i][j][k + 1]))) + A[t%2][i][j][k]);
#define S1(t,i,j,k) S1_I((t),(i),(j),(k))
#pragma scop
if (T + N >= 17) {
  if (N <= 7) {
    for (register int k = 1; k <= 2; k += 1) {
      if (k == 2) {
        const int ii2_73987_lb = 0, ii2_20304_ub = min((T + N - 3) / 8, (N + 4) / 8 + 1);
        for (register int ii2 = ii2_73987_lb; ii2 <= ii2_20304_ub; ii2 += 1) {
          const int ii3_44878_lb = max(0, ii2 - (N + 6) / 8), ii3_65030_ub = min(min((T + N - 3) / 8, ii2 - (-N + 9) / 8 + 1), (N + 4) / 8 + 1);
          for (register int ii3 = ii3_44878_lb; ii3 <= ii3_65030_ub; ii3 += 1) {
            const int i0_29528_lb = max(max(1, -N + 8 * ii2 + 1), -N + 8 * ii3 + 1), i0_92423_ub = min(min(min(13, T - 2), 8 * ii2 + 7), 8 * ii3 + 7);
            for (register int i0 = i0_29528_lb; i0 <= i0_92423_ub; i0 += 1) {
              const int i1_72085_lb = 1, i1_39154_ub = min(min(N, i0), -i0 + 14);
              for (register int i1 = i1_72085_lb; i1 <= i1_39154_ub; i1 += 1) {
                const int i2_50049_lb = max(1, 8 * ii2 - i0 + 1), i2_68710_ub = min(N, 8 * ii2 - i0 + 8);
                for (register int i2 = i2_50049_lb; i2 <= i2_68710_ub; i2 += 1) {
                  const int i3_40407_lb = max(1, 8 * ii3 - i0 + 1), i3_55466_ub = min(N, 8 * ii3 - i0 + 8);
                  for (register int i3 = i3_40407_lb; i3 <= i3_55466_ub; i3 += 1) {
                    S1(i0, i1, i2, i3);
                  }
                }
              }
            }
          }
        }
      } else {
        const int ii2_78227_lb = 0, ii2_22861_ub = floord(N - 1, 4);
        for (register int ii2 = ii2_78227_lb; ii2 <= ii2_22861_ub; ii2 += 1) {
          const int ii3_50594_lb = 0, ii3_97908_ub = (N - 1) / 4;
          for (register int ii3 = ii3_50594_lb; ii3 <= ii3_97908_ub; ii3 += 1) {
            const int i0_46213_lb = max(max(0, -N + 8 * ii2), -N + 8 * ii3), i0_84364_ub = N;
            for (register int i0 = i0_46213_lb; i0 < i0_84364_ub; i0 += 1) {
              const int i1_29754_lb = i0 + 1, i1_62490_ub = N;
              for (register int i1 = i1_29754_lb; i1 <= i1_62490_ub; i1 += 1) {
                const int i2_67806_lb = max(1, 8 * ii2 - i0), i2_80623_ub = min(N, 8 * ii2 - i0 + 7);
                for (register int i2 = i2_67806_lb; i2 <= i2_80623_ub; i2 += 1) {
                  const int i3_49161_lb = max(1, 8 * ii3 - i0), i3_57362_ub = min(N, 8 * ii3 - i0 + 7);
                  for (register int i3 = i3_49161_lb; i3 <= i3_57362_ub; i3 += 1) {
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
  const int ii0_70272_lb = max(0, -((N + 6) / 7) + 2), ii0_94410_ub = min((T + 2) / 8, (T + N - 1) / 8);
  for (register int ii0 = ii0_70272_lb; ii0 < ii0_94410_ub; ii0 += 1) {
    if (ii0 == 0) {
      const int ii1_70598_lb = 0, ii1_31224_ub = (N - 1) / 8;
      #pragma omp parallel for
      for (register int ii1 = ii1_70598_lb; ii1 < ii1_31224_ub; ii1 += 1) {
        const int ii2_55812_lb = 0, ii2_79035_ub = (N + 2) / 8;
        for (register int ii2 = ii2_55812_lb; ii2 <= ii2_79035_ub; ii2 += 1) {
          const int ii3_64618_lb = 0, ii3_46151_ub = (N + 2) / 8;
          for (register int ii3 = ii3_64618_lb; ii3 <= ii3_46151_ub; ii3 += 1) {
            const int i0_99339_lb = max(max(0, -N + 8 * ii2), -N + 8 * ii3), i0_25849_ub = min(2, N - 8 * ii1 - 9);
            for (register int i0 = i0_99339_lb; i0 <= i0_25849_ub; i0 += 1) {
              const int i1_11181_lb = 8 * ii1 + i0 + 9, i1_45219_ub = min(N, 8 * ii1 - i0 + 14);
              for (register int i1 = i1_11181_lb; i1 <= i1_45219_ub; i1 += 1) {
                const int i2_34624_lb = max(1, 8 * ii2 - i0), i2_99618_ub = min(N, 8 * ii2 - i0 + 7);
                for (register int i2 = i2_34624_lb; i2 <= i2_99618_ub; i2 += 1) {
                  const int i3_84373_lb = max(1, 8 * ii3 - i0), i3_84673_ub = min(N, 8 * ii3 - i0 + 7);
                  for (register int i3 = i3_84373_lb; i3 <= i3_84673_ub; i3 += 1) {
                    S1(i0, i1, i2, i3);
                  }
                }
              }
            }
          }
        }
      }
      if (N <= 11) {
        const int ii2_68328_lb = 0, ii2_41133_ub = min((N - 2) / 8 + 1, (T + N - 2) / 8);
        for (register int ii2 = ii2_68328_lb; ii2 <= ii2_41133_ub; ii2 += 1) {
          const int ii3_40140_lb = 0, ii3_62908_ub = min((N - 2) / 8 + 1, (T + N - 2) / 8);
          for (register int ii3 = ii3_40140_lb; ii3 <= ii3_62908_ub; ii3 += 1) {
            const int i0_80346_lb = max(max(0, -N + 8 * ii2), -N + 8 * ii3), i0_7086_ub = min(6, T - 2);
            for (register int i0 = i0_80346_lb; i0 <= i0_7086_ub; i0 += 1) {
              const int i1_60816_lb = i0 + 1, i1_26559_ub = min(min(N, i0 + 8), -i0 + 14);
              for (register int i1 = i1_60816_lb; i1 <= i1_26559_ub; i1 += 1) {
                const int i2_91450_lb = max(1, 8 * ii2 - i0), i2_6922_ub = min(N, 8 * ii2 - i0 + 7);
                for (register int i2 = i2_91450_lb; i2 <= i2_6922_ub; i2 += 1) {
                  const int i3_5402_lb = max(1, 8 * ii3 - i0), i3_75608_ub = min(N, 8 * ii3 - i0 + 7);
                  for (register int i3 = i3_5402_lb; i3 <= i3_75608_ub; i3 += 1) {
                    S1(i0, i1, i2, i3);
                  }
                }
              }
            }
          }
        }
      }
    }
    if (N >= 4 && N + 8 * ii0 >= 12) {
      if (ii0 == 0) {
        const int ii2_3897_lb = 0, ii2_54563_ub = min((N - 2) / 8 + 1, (T + N - 2) / 8);
        for (register int ii2 = ii2_3897_lb; ii2 <= ii2_54563_ub; ii2 += 1) {
          const int ii3_49323_lb = 0, ii3_74169_ub = min((N - 2) / 8 + 1, (T + N - 2) / 8);
          for (register int ii3 = ii3_49323_lb; ii3 <= ii3_74169_ub; ii3 += 1) {
            const int i0_65325_lb = max(max(0, -N + 8 * ii2), -N + 8 * ii3), i0_19921_ub = min(6, T - 2);
            for (register int i0 = i0_65325_lb; i0 <= i0_19921_ub; i0 += 1) {
              const int i1_5394_lb = i0 + 1, i1_37489_ub = min(i0 + 8, -i0 + 14);
              for (register int i1 = i1_5394_lb; i1 <= i1_37489_ub; i1 += 1) {
                const int i2_15308_lb = max(1, 8 * ii2 - i0), i2_86364_ub = min(N, 8 * ii2 - i0 + 7);
                for (register int i2 = i2_15308_lb; i2 <= i2_86364_ub; i2 += 1) {
                  const int i3_83640_lb = max(1, 8 * ii3 - i0), i3_30999_ub = min(N, 8 * ii3 - i0 + 7);
                  for (register int i3 = i3_83640_lb; i3 <= i3_30999_ub; i3 += 1) {
                    S1(i0, i1, i2, i3);
                  }
                }
              }
            }
          }
        }
      }
      const int ii1_12213_lb = max(0, -ii0 + 1), ii1_11173_ub = (N + 4) / 8;
      #pragma omp parallel for
      for (register int ii1 = ii1_12213_lb; ii1 < ii1_11173_ub; ii1 += 1) {
        const int ii2_76219_lb = 0, ii2_46837_ub = min(min(-ii1 + (N - 1) / 4, (N - 2) / 8 + 1), -ii0 + (T + N - 2) / 8);
        for (register int ii2 = ii2_76219_lb; ii2 <= ii2_46837_ub; ii2 += 1) {
          const int ii3_10791_lb = 0, ii3_60592_ub = min(min(-ii1 + (N - 1) / 4, (N - 2) / 8 + 1), -ii0 + (T + N - 2) / 8);
          for (register int ii3 = ii3_10791_lb; ii3 <= ii3_60592_ub; ii3 += 1) {
            const int i0_47863_lb = max(max(max(8 * ii0, -N + 8 * ii0 + 8 * ii1 + 7), -N + 8 * ii0 + 8 * ii2), -N + 8 * ii0 + 8 * ii3), i0_95472_ub = min(min(T - 2, 8 * ii0 + 6), N + 8 * ii0 - 8 * ii1 - 1);
            for (register int i0 = i0_47863_lb; i0 <= i0_95472_ub; i0 += 1) {
              const int i1_1725_lb = max(-8 * ii0 + 8 * ii1 + i0 + 1, 8 * ii0 + 8 * ii1 - i0 + 7), i1_4355_ub = min(min(N, -8 * ii0 + 8 * ii1 + i0 + 8), 8 * ii0 + 8 * ii1 - i0 + 14);
              for (register int i1 = i1_1725_lb; i1 <= i1_4355_ub; i1 += 1) {
                const int i2_74732_lb = max(1, 8 * ii0 + 8 * ii2 - i0), i2_82072_ub = min(N, 8 * ii0 + 8 * ii2 - i0 + 7);
                for (register int i2 = i2_74732_lb; i2 <= i2_82072_ub; i2 += 1) {
                  const int i3_27793_lb = max(1, 8 * ii0 + 8 * ii3 - i0), i3_51900_ub = min(N, 8 * ii0 + 8 * ii3 - i0 + 7);
                  for (register int i3 = i3_27793_lb; i3 <= i3_51900_ub; i3 += 1) {
                    S1(i0, i1, i2, i3);
                  }
                }
              }
            }
          }
        }
      }
      if (N >= 8 && 8 * ii0 + 16 >= T + N) {
        for (register int ii2 = 0; ii2 < 2; ii2 += 1) {
          for (register int ii3 = 0; ii3 < 2; ii3 += 1) {
            const int i0_11203_lb = 8 * ii0 + 4, i0_79071_ub = T - 1;
            for (register int i0 = i0_11203_lb; i0 < i0_79071_ub; i0 += 1) {
              const int i1_1300_lb = 8 * ii0 - i0 + 7, i1_76878_ub = -8 * ii0 + i0;
              for (register int i1 = i1_1300_lb; i1 <= i1_76878_ub; i1 += 1) {
                const int i2_53241_lb = max(1, 8 * ii0 + 8 * ii2 - i0 + 4), i2_82978_ub = min(N, 8 * ii0 + 8 * ii2 - i0 + 11);
                for (register int i2 = i2_53241_lb; i2 <= i2_82978_ub; i2 += 1) {
                  const int i3_96799_lb = max(1, 8 * ii0 + 8 * ii3 - i0 + 4), i3_74987_ub = min(N, 8 * ii0 + 8 * ii3 - i0 + 11);
                  for (register int i3 = i3_96799_lb; i3 <= i3_74987_ub; i3 += 1) {
                    S1(i0, i1, i2, i3);
                  }
                }
              }
            }
          }
        }
      }
    }
    if (N <= 7) {
      const int ii2_36819_lb = 0, ii2_28459_ub = min((N + 9) / 8, -ii0 + (T + N + 2) / 8 - 1);
      for (register int ii2 = ii2_36819_lb; ii2 <= ii2_28459_ub; ii2 += 1) {
        const int ii3_77703_lb = max(0, ii2 - (N + 6) / 8), ii3_20460_ub = min(min(ii2 - (-N + 9) / 8 + 1, (N + 9) / 8), -ii0 + (T + N + 2) / 8 - 1);
        for (register int ii3 = ii3_77703_lb; ii3 <= ii3_20460_ub; ii3 += 1) {
          const int i0_59459_lb = max(max(max(8 * ii0 + 4, -N + 8 * ii0 + 7), -N + 8 * ii0 + 8 * ii2 + 4), -N + 8 * ii0 + 8 * ii3 + 4), i0_6269_ub = min(min(min(T - 2, 8 * ii0 + 13), 8 * ii0 + 8 * ii2 + 10), 8 * ii0 + 8 * ii3 + 10);
          for (register int i0 = i0_59459_lb; i0 <= i0_6269_ub; i0 += 1) {
            const int i1_47985_lb = max(1, 8 * ii0 - i0 + 7), i1_35678_ub = min(min(N, -8 * ii0 + i0), 8 * ii0 - i0 + 14);
            for (register int i1 = i1_47985_lb; i1 <= i1_35678_ub; i1 += 1) {
              const int i2_69458_lb = max(1, 8 * ii0 + 8 * ii2 - i0 + 4), i2_75129_ub = min(N, 8 * ii0 + 8 * ii2 - i0 + 11);
              for (register int i2 = i2_69458_lb; i2 <= i2_75129_ub; i2 += 1) {
                const int i3_96270_lb = max(1, 8 * ii0 + 8 * ii3 - i0 + 4), i3_33673_ub = min(N, 8 * ii0 + 8 * ii3 - i0 + 11);
                for (register int i3 = i3_96270_lb; i3 <= i3_33673_ub; i3 += 1) {
                  S1(i0, i1, i2, i3);
                }
              }
            }
          }
        }
      }
    } else if (T + N >= 8 * ii0 + 17) {
      const int ii1_86953_lb = 0, ii1_14348_ub = min(-ii0 + (T + N - 1) / 8 - 1, N / 8);
      #pragma omp parallel for
      for (register int ii1 = ii1_86953_lb; ii1 <= ii1_14348_ub; ii1 += 1) {
        if (ii1 >= 1) {
          const int ii2_38028_lb = 0, ii2_61685_ub = min(min(-ii1 + (N + 1) / 4, (N - 2) / 8 + 1), -ii0 + (T + N + 2) / 8 - 1);
          for (register int ii2 = ii2_38028_lb; ii2 <= ii2_61685_ub; ii2 += 1) {
            const int ii3_12772_lb = 0, ii3_82173_ub = min(min(-ii1 + (N + 1) / 4, (N - 2) / 8 + 1), -ii0 + (T + N + 2) / 8 - 1);
            for (register int ii3 = ii3_12772_lb; ii3 <= ii3_82173_ub; ii3 += 1) {
              const int i0_29937_lb = max(max(max(8 * ii0 + 4, -N + 8 * ii0 + 8 * ii1 + 7), -N + 8 * ii0 + 8 * ii2 + 4), -N + 8 * ii0 + 8 * ii3 + 4), i0_37755_ub = min(min(T - 2, 8 * ii0 + 10), N + 8 * ii0 - 8 * ii1 + 7);
              for (register int i0 = i0_29937_lb; i0 <= i0_37755_ub; i0 += 1) {
                const int i1_17768_lb = max(-8 * ii0 + 8 * ii1 + i0 - 7, 8 * ii0 + 8 * ii1 - i0 + 7), i1_5111_ub = min(min(N, -8 * ii0 + 8 * ii1 + i0), 8 * ii0 + 8 * ii1 - i0 + 14);
                for (register int i1 = i1_17768_lb; i1 <= i1_5111_ub; i1 += 1) {
                  const int i2_84493_lb = max(1, 8 * ii0 + 8 * ii2 - i0 + 4), i2_45324_ub = min(N, 8 * ii0 + 8 * ii2 - i0 + 11);
                  for (register int i2 = i2_84493_lb; i2 <= i2_45324_ub; i2 += 1) {
                    const int i3_84182_lb = max(1, 8 * ii0 + 8 * ii3 - i0 + 4), i3_2145_ub = min(N, 8 * ii0 + 8 * ii3 - i0 + 11);
                    for (register int i3 = i3_84182_lb; i3 <= i3_2145_ub; i3 += 1) {
                      S1(i0, i1, i2, i3);
                    }
                  }
                }
              }
            }
          }
        } else if (ii0 >= 1) {
          const int ii2_22202_lb = 0, ii2_53775_ub = min((N + 1) / 8 + 1, -ii0 + (T + N + 2) / 8 - 1);
          for (register int ii2 = ii2_22202_lb; ii2 <= ii2_53775_ub; ii2 += 1) {
            const int ii3_85123_lb = max(0, ii2 - (N + 6) / 8), ii3_35354_ub = min(min(ii2 + (N - 2) / 8 + 1, (N + 1) / 8 + 1), -ii0 + (T + N + 2) / 8 - 1);
            for (register int ii3 = ii3_85123_lb; ii3 <= ii3_35354_ub; ii3 += 1) {
              const int i0_45114_lb = max(max(8 * ii0 + 4, -N + 8 * ii0 + 8 * ii2 + 4), -N + 8 * ii0 + 8 * ii3 + 4), i0_21943_ub = min(min(min(T - 2, 8 * ii0 + 13), 8 * ii0 + 8 * ii2 + 10), 8 * ii0 + 8 * ii3 + 10);
              for (register int i0 = i0_45114_lb; i0 <= i0_21943_ub; i0 += 1) {
                const int i1_63813_lb = max(1, 8 * ii0 - i0 + 7), i1_22818_ub = min(-8 * ii0 + i0, 8 * ii0 - i0 + 14);
                for (register int i1 = i1_63813_lb; i1 <= i1_22818_ub; i1 += 1) {
                  const int i2_58755_lb = max(1, 8 * ii0 + 8 * ii2 - i0 + 4), i2_23272_ub = min(N, 8 * ii0 + 8 * ii2 - i0 + 11);
                  for (register int i2 = i2_58755_lb; i2 <= i2_23272_ub; i2 += 1) {
                    const int i3_29087_lb = max(1, 8 * ii0 + 8 * ii3 - i0 + 4), i3_23092_ub = min(N, 8 * ii0 + 8 * ii3 - i0 + 11);
                    for (register int i3 = i3_29087_lb; i3 <= i3_23092_ub; i3 += 1) {
                      S1(i0, i1, i2, i3);
                    }
                  }
                }
              }
            }
          }
        } else {
          const int ii2_58950_lb = 0, ii2_14897_ub = min((T + N - 3) / 8, (N + 4) / 8 + 1);
          for (register int ii2 = ii2_58950_lb; ii2 <= ii2_14897_ub; ii2 += 1) {
            const int ii3_98221_lb = max(0, ii2 - (N + 6) / 8), ii3_71573_ub = min(min((T + N - 3) / 8, ii2 + (N - 2) / 8 + 1), (N + 4) / 8 + 1);
            for (register int ii3 = ii3_98221_lb; ii3 <= ii3_71573_ub; ii3 += 1) {
              const int i0_48571_lb = max(max(1, -N + 8 * ii2 + 1), -N + 8 * ii3 + 1), i0_1526_ub = min(min(min(13, T - 2), 8 * ii2 + 7), 8 * ii3 + 7);
              for (register int i0 = i0_48571_lb; i0 <= i0_1526_ub; i0 += 1) {
                const int i1_85921_lb = 1, i1_2951_ub = min(i0, -i0 + 14);
                for (register int i1 = i1_85921_lb; i1 <= i1_2951_ub; i1 += 1) {
                  const int i2_79563_lb = max(1, 8 * ii2 - i0 + 1), i2_98693_ub = min(N, 8 * ii2 - i0 + 8);
                  for (register int i2 = i2_79563_lb; i2 <= i2_98693_ub; i2 += 1) {
                    const int i3_85125_lb = max(1, 8 * ii3 - i0 + 1), i3_9500_ub = min(N, 8 * ii3 - i0 + 8);
                    for (register int i3 = i3_85125_lb; i3 <= i3_9500_ub; i3 += 1) {
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
  if (((T + 2) % 8) + N >= 11 && (T + 2) % 8 >= 4) {
    if (T <= 5) {
      const int ii1_52800_lb = 0, ii1_2893_ub = (N - 1) / 8;
      #pragma omp parallel for
      for (register int ii1 = ii1_52800_lb; ii1 < ii1_2893_ub; ii1 += 1) {
        const int ii2_14611_lb = 0, ii2_53645_ub = min((T + N - 2) / 8, (N + 2) / 8);
        for (register int ii2 = ii2_14611_lb; ii2 <= ii2_53645_ub; ii2 += 1) {
          const int ii3_64569_lb = 0, ii3_15146_ub = min((T + N - 2) / 8, (N + 2) / 8);
          for (register int ii3 = ii3_64569_lb; ii3 <= ii3_15146_ub; ii3 += 1) {
            const int i0_55791_lb = max(max(0, -N + 8 * ii2), -N + 8 * ii3), i0_86772_ub = min(min(2, T - 2), N - 8 * ii1 - 9);
            for (register int i0 = i0_55791_lb; i0 <= i0_86772_ub; i0 += 1) {
              const int i1_68921_lb = 8 * ii1 + i0 + 9, i1_57266_ub = min(N, 8 * ii1 - i0 + 14);
              for (register int i1 = i1_68921_lb; i1 <= i1_57266_ub; i1 += 1) {
                const int i2_38478_lb = max(1, 8 * ii2 - i0), i2_14036_ub = min(N, 8 * ii2 - i0 + 7);
                for (register int i2 = i2_38478_lb; i2 <= i2_14036_ub; i2 += 1) {
                  const int i3_95561_lb = max(1, 8 * ii3 - i0), i3_18643_ub = min(N, 8 * ii3 - i0 + 7);
                  for (register int i3 = i3_95561_lb; i3 <= i3_18643_ub; i3 += 1) {
                    S1(i0, i1, i2, i3);
                  }
                }
              }
            }
          }
        }
      }
    }
    const int ii1_53206_lb = 0, ii1_54316_ub = -((T - 2) / 8) + (T + N - 1) / 8;
    #pragma omp parallel for
    for (register int ii1 = ii1_53206_lb; ii1 < ii1_54316_ub; ii1 += 1) {
      const int ii2_41916_lb = 0, ii2_98645_ub = -((T - 2) / 8) + (T + N - 2) / 8;
      for (register int ii2 = ii2_41916_lb; ii2 <= ii2_98645_ub; ii2 += 1) {
        const int ii3_93761_lb = 0, ii3_17218_ub = -((T - 2) / 8) + (T + N - 2) / 8;
        for (register int ii3 = ii3_93761_lb; ii3 <= ii3_17218_ub; ii3 += 1) {
          const int i0_13542_lb = max(max(max(max(-((T + 6) % 8) + T - N + 8 * ii3 - 2, -((T + 6) % 8) + T - N + 8 * ii2 - 2), -((T + 6) % 8) + T - 2), -4 * ii1 + 4 * ((T - 2) / 8) + 4), -((T + 6) % 8) + T - N + 8 * ii1 + 5), i0_8334_ub = T - 1;
          for (register int i0 = i0_13542_lb; i0 < i0_8334_ub; i0 += 1) {
            const int i1_5143_lb = -((T + 6) % 8) + T + 8 * ii1 - i0 + 5, i1_62113_ub = min(N, ((T + 6) % 8) - T + 8 * ii1 + i0 + 10);
            for (register int i1 = i1_5143_lb; i1 <= i1_62113_ub; i1 += 1) {
              const int i2_9861_lb = max(1, -((T + 6) % 8) + T + 8 * ii2 - i0 - 2), i2_7416_ub = min(N, -((T + 6) % 8) + T + 8 * ii2 - i0 + 5);
              for (register int i2 = i2_9861_lb; i2 <= i2_7416_ub; i2 += 1) {
                const int i3_65065_lb = max(1, -((T + 6) % 8) + T + 8 * ii3 - i0 - 2), i3_89424_ub = min(N, -((T + 6) % 8) + T + 8 * ii3 - i0 + 5);
                for (register int i3 = i3_65065_lb; i3 <= i3_89424_ub; i3 += 1) {
                  S1(i0, i1, i2, i3);
                }
              }
            }
          }
          if (T <= 5 && ii1 == 0) {
            const int i0_22461_lb = max(max(0, -N + 8 * ii2), -N + 8 * ii3), i0_50190_ub = T - 1;
            for (register int i0 = i0_22461_lb; i0 < i0_50190_ub; i0 += 1) {
              const int i1_98925_lb = i0 + 1, i1_91614_ub = i0 + 8;
              for (register int i1 = i1_98925_lb; i1 <= i1_91614_ub; i1 += 1) {
                const int i2_69435_lb = max(1, 8 * ii2 - i0), i2_29888_ub = min(N, 8 * ii2 - i0 + 7);
                for (register int i2 = i2_69435_lb; i2 <= i2_29888_ub; i2 += 1) {
                  const int i3_45259_lb = max(1, 8 * ii3 - i0), i3_34005_ub = min(N, 8 * ii3 - i0 + 7);
                  for (register int i3 = i3_45259_lb; i3 <= i3_34005_ub; i3 += 1) {
                    S1(i0, i1, i2, i3);
                  }
                }
              }
            }
          }
        }
      }
    }
    if (T <= 5) {
      const int ii2_61386_lb = 0, ii2_1050_ub = (T + N - 3) / 8;
      for (register int ii2 = ii2_61386_lb; ii2 <= ii2_1050_ub; ii2 += 1) {
        const int ii3_37129_lb = 0, ii3_30308_ub = (T + N - 3) / 8;
        for (register int ii3 = ii3_37129_lb; ii3 <= ii3_30308_ub; ii3 += 1) {
          const int i0_58317_lb = max(max(1, -N + 8 * ii2 + 1), -N + 8 * ii3 + 1), i0_75607_ub = T - 1;
          for (register int i0 = i0_58317_lb; i0 < i0_75607_ub; i0 += 1) {
            for (register int i1 = 1; i1 <= i0; i1 += 1) {
              const int i2_94250_lb = max(1, 8 * ii2 - i0 + 1), i2_30254_ub = min(N, 8 * ii2 - i0 + 8);
              for (register int i2 = i2_94250_lb; i2 <= i2_30254_ub; i2 += 1) {
                const int i3_24547_lb = max(1, 8 * ii3 - i0 + 1), i3_52518_ub = min(N, 8 * ii3 - i0 + 8);
                for (register int i3 = i3_24547_lb; i3 <= i3_52518_ub; i3 += 1) {
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
  const int ii2_28899_lb = 0, ii2_34660_ub = floord(T + N - 2, 8);
  for (register int ii2 = ii2_28899_lb; ii2 <= ii2_34660_ub; ii2 += 1) {
    const int ii3_86089_lb = 0, ii3_42441_ub = (T + N - 2) / 8;
    for (register int ii3 = ii3_86089_lb; ii3 <= ii3_42441_ub; ii3 += 1) {
      const int i0_42994_lb = 0, i0_7584_ub = min(T - 1, N - 8);
      for (register int i0 = i0_42994_lb; i0 < i0_7584_ub; i0 += 1) {
        const int i1_4555_lb = i0 + 9, i1_52855_ub = N;
        for (register int i1 = i1_4555_lb; i1 <= i1_52855_ub; i1 += 1) {
          const int i2_15001_lb = max(1, 8 * ii2 - i0), i2_69620_ub = min(N, 8 * ii2 - i0 + 7);
          for (register int i2 = i2_15001_lb; i2 <= i2_69620_ub; i2 += 1) {
            const int i3_58632_lb = max(1, 8 * ii3 - i0), i3_53814_ub = min(N, 8 * ii3 - i0 + 7);
            for (register int i3 = i3_58632_lb; i3 <= i3_53814_ub; i3 += 1) {
              S1(i0, i1, i2, i3);
            }
          }
        }
      }
    }
  }
  const int k_36162_lb = 1, k_73909_ub = min(2, T - 1);
  for (register int k = k_36162_lb; k <= k_73909_ub; k += 1) {
    if (k == 2) {
      const int ii2_45428_lb = 0, ii2_5597_ub = floord(T + N - 3, 8);
      for (register int ii2 = ii2_45428_lb; ii2 <= ii2_5597_ub; ii2 += 1) {
        const int ii3_20149_lb = max(0, ii2 - (N + 6) / 8), ii3_90688_ub = min((T + N - 3) / 8, ii2 + (N + 6) / 8);
        for (register int ii3 = ii3_20149_lb; ii3 <= ii3_90688_ub; ii3 += 1) {
          const int i0_39602_lb = max(max(1, -N + 8 * ii2 + 1), -N + 8 * ii3 + 1), i0_81536_ub = min(min(T - 2, 8 * ii2 + 7), 8 * ii3 + 7);
          for (register int i0 = i0_39602_lb; i0 <= i0_81536_ub; i0 += 1) {
            const int i1_91738_lb = 1, i1_76731_ub = min(N, i0);
            for (register int i1 = i1_91738_lb; i1 <= i1_76731_ub; i1 += 1) {
              const int i2_28196_lb = max(1, 8 * ii2 - i0 + 1), i2_66407_ub = min(N, 8 * ii2 - i0 + 8);
              for (register int i2 = i2_28196_lb; i2 <= i2_66407_ub; i2 += 1) {
                const int i3_68690_lb = max(1, 8 * ii3 - i0 + 1), i3_5244_ub = min(N, 8 * ii3 - i0 + 8);
                for (register int i3 = i3_68690_lb; i3 <= i3_5244_ub; i3 += 1) {
                  S1(i0, i1, i2, i3);
                }
              }
            }
          }
        }
      }
    } else {
      const int ii2_36638_lb = 0, ii2_79293_ub = min(floord(N - 1, 4), floord(T + N - 2, 8));
      for (register int ii2 = ii2_36638_lb; ii2 <= ii2_79293_ub; ii2 += 1) {
        const int ii3_35498_lb = 0, ii3_77537_ub = min((N - 1) / 4, (T + N - 2) / 8);
        for (register int ii3 = ii3_35498_lb; ii3 <= ii3_77537_ub; ii3 += 1) {
          const int i0_48163_lb = max(max(0, -N + 8 * ii2), -N + 8 * ii3), i0_64397_ub = min(T - 1, N);
          for (register int i0 = i0_48163_lb; i0 < i0_64397_ub; i0 += 1) {
            const int i1_12197_lb = i0 + 1, i1_34252_ub = min(N, i0 + 8);
            for (register int i1 = i1_12197_lb; i1 <= i1_34252_ub; i1 += 1) {
              const int i2_23190_lb = max(1, 8 * ii2 - i0), i2_55191_ub = min(N, 8 * ii2 - i0 + 7);
              for (register int i2 = i2_23190_lb; i2 <= i2_55191_ub; i2 += 1) {
                const int i3_58189_lb = max(1, 8 * ii3 - i0), i3_27745_ub = min(N, 8 * ii3 - i0 + 7);
                for (register int i3 = i3_58189_lb; i3 <= i3_27745_ub; i3 += 1) {
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
 