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
#define N 1600L
#define T 500L
#endif

#define NUM_FP_OPS 10

/* Define our arrays */
double A[2][N][N];
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
  long int t, i, j, k;
  const int BASE = 1024;

  // for timekeeping
  int ts_return = -1;
  struct timeval start, end, result;
  double tdiff = 0.0;
  double total;

  printf("Number of points = %ld\t|Number of timesteps = %ld\t", N * N, T);

  /* Initialization */
  srand(42); // seed with a constant value to verify results

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      A[0][i][j] = 1.0 * (rand() % BASE);
    }
  }

#ifdef TIME
  gettimeofday(&start, 0);
#endif

  short _N = N - 1;
  int _T = T;

/* TC Optimizing Compiler 0.4.1 */
/* ./tc ../examples/pluto+/heat-2dp-nt-t.scop.c --diamond-tiling --omp-for-codegen --iterative-tc --debug -b 32 --use-macros */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define S1_I(t,i,j) A[(t + 1) % 2][i][j] = (((0.125 * ((A[t % 2][i == 1599 ? 0 : i + 1][j] - (2.0 * A[t % 2][i][j])) + A[t % 2][i == 0 ? 1599 : i - 1][j])) + (0.125 * ((A[t % 2][i][j == 1599 ? 0 : j + 1] - (2.0 * A[t % 2][i][j])) + A[t % 2][i][j == 0 ? 1599 : j - 1]))) + A[t % 2][i][j]);
#define S1(t,i,j) S1_I((t),(i),(j))
#pragma scop
for (register int ii0 = 0; ii0 <= 15; ii0 += 1) {
  if (ii0 >= 1) {
    for (register int k = 1; k <= 2; k += 1) {
      if (ii0 == 15) {
        const int ii1_37629_lb = 0, ii1_23576_ub = k + 47;
        #pragma omp parallel for
        for (register int ii1 = ii1_37629_lb; ii1 <= ii1_23576_ub; ii1 += 1) {
          if (ii1 >= k) {
            for (register int ii2 = 0; ii2 <= 49; ii2 += 1) {
              if (k == 2) {
                const int i0_31833_lb = max(16 * ii1 + 442, -16 * ii1 + 522), i0_84366_ub = min(499, 16 * ii1 + 8 * ii2 + 464);
                for (register int i0 = i0_31833_lb; i0 <= i0_84366_ub; i0 += 1) {
                  {
                    if (ii1 == 3 && ii2 == 0) {
                      const int i1_87371_lb = -i0 + 595, i1_83528_ub = 105;
                      for (register int i1 = i1_87371_lb; i1 <= i1_83528_ub; i1 += 1) {
                        const int i2_6841_lb = i0 + i1 - 595, i2_57957_ub = -i0 - i1 + 626;
                        for (register int i2 = i2_6841_lb; i2 <= i2_57957_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                      }
                    } else if (ii1 == 2) {
                      const int i1_61717_lb = -i0 + 563, i1_78557_ub = i0 - 415;
                      for (register int i1 = i1_61717_lb; i1 < i1_78557_ub; i1 += 1) {
                        if (ii2 <= 48) {
                          const int i2_75400_lb = max(i0 + i1 - 563, 32 * ii2 - i0 - i1 + 575), i2_84723_ub = 32 * ii2 - i0 + i1 + 415;
                          for (register int i2 = i2_75400_lb; i2 <= i2_84723_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                          if (ii2 == 0) {
                            const int i2_5298_lb = max(i0 + i1 - 575, i0 - i1 - 416), i2_30763_ub = i0 + i1 - 563;
                            for (register int i2 = i2_5298_lb; i2 < i2_30763_ub; i2 += 1) {
                              S1(i0, i1, i2);
                            }
                          }
                          const int i2_89516_lb = max(max(i0 + i1 - 563, 32 * ii2 - i0 + i1 + 416), 32 * ii2 - i0 - i1 + 563), i2_82723_ub = 32 * ii2 - i0 - i1 + 594;
                          for (register int i2 = i2_89516_lb; i2 <= i2_82723_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                          const int i2_73151_lb = 32 * ii2 - i0 - i1 + 595, i2_9110_ub = min(32 * ii2 - i0 + i1 + 447, 32 * ii2 - i0 - i1 + 606);
                          for (register int i2 = i2_73151_lb; i2 <= i2_9110_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                        } else {
                          const int i2_15951_lb = 0, i2_89079_ub = i0 + i1 - 575;
                          for (register int i2 = i2_15951_lb; i2 < i2_89079_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                          const int i2_95658_lb = -i0 - i1 + 2143, i2_32241_ub = -i0 + i1 + 1983;
                          for (register int i2 = i2_95658_lb; i2 <= i2_32241_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                          const int i2_93765_lb = max(0, i0 + i1 - 575), i2_11753_ub = min(i0 + i1 - 563, i0 - i1 - 416);
                          for (register int i2 = i2_93765_lb; i2 < i2_11753_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                          const int i2_37348_lb = max(-i0 + i1 + 1984, -i0 - i1 + 2131), i2_55512_ub = 1599;
                          for (register int i2 = i2_37348_lb; i2 <= i2_55512_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                        }
                      }
                      if (ii2 == 0) {
                        const int i1_87899_lb = max(i0 - 415, -i0 + 575), i1_95594_ub = 85;
                        for (register int i1 = i1_87899_lb; i1 <= i1_95594_ub; i1 += 1) {
                          const int i2_78592_lb = i0 + i1 - 575, i2_25528_ub = -i0 - i1 + 606;
                          for (register int i2 = i2_78592_lb; i2 <= i2_25528_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                        }
                        const int i1_19171_lb = 86, i1_2436_ub = min(89, i0 - 404);
                        for (register int i1 = i1_19171_lb; i1 <= i1_2436_ub; i1 += 1) {
                          const int i2_98735_lb = i0 - i1 - 404, i2_67356_ub = -i0 + i1 + 435;
                          for (register int i2 = i2_98735_lb; i2 <= i2_67356_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                        }
                      }
                    }
                    if (ii2 == 0) {
                      const int i1_86802_lb = 16 * ii1 + 58, i1_2459_ub = min(32 * ii1 + i0 - 467, i0 - 383);
                      for (register int i1 = i1_86802_lb; i1 < i1_2459_ub; i1 += 1) {
                        if (ii1 == 2) {
                          const int i2_67236_lb = i0 - i1 - 404, i2_93643_ub = i0 - i1 - 384;
                          for (register int i2 = i2_67236_lb; i2 < i2_93643_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                        }
                        const int i2_76768_lb = max(i0 + i1 - 607, i0 - i1 - 384), i2_28954_ub = -32 * ii1 + i0 + i1 - 499;
                        for (register int i2 = i2_76768_lb; i2 < i2_28954_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                        if (ii1 == 2) {
                          const int i2_88552_lb = i0 + i1 - 563, i2_68520_ub = -i0 + i1 + 435;
                          for (register int i2 = i2_88552_lb; i2 <= i2_68520_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                        } else {
                          const int i2_13677_lb = i0 + i1 - 595, i2_10202_ub = min(-i0 + i1 + 415, -i0 - i1 + 638);
                          for (register int i2 = i2_13677_lb; i2 <= i2_10202_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                        }
                      }
                    } else if (ii1 == 3 && ii2 <= 48) {
                      const int i1_99284_lb = -i0 + 595, i1_19545_ub = 105;
                      for (register int i1 = i1_99284_lb; i1 <= i1_19545_ub; i1 += 1) {
                        const int i2_92926_lb = 32 * ii2 - i0 - i1 + 595, i2_72435_ub = 32 * ii2 - i0 - i1 + 626;
                        for (register int i2 = i2_92926_lb; i2 <= i2_72435_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                      }
                      const int i1_28656_lb = 106, i1_8877_ub = i0 - 383;
                      for (register int i1 = i1_28656_lb; i1 < i1_8877_ub; i1 += 1) {
                        const int i2_77866_lb = 32 * ii2 - i0 - i1 + 607, i2_40666_ub = 32 * ii2 - i0 + i1 + 383;
                        for (register int i2 = i2_77866_lb; i2 <= i2_40666_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                        const int i2_57470_lb = 32 * ii2 - i0 + i1 + 384, i2_87983_ub = min(32 * ii2 - i0 + i1 + 415, 32 * ii2 - i0 - i1 + 638);
                        for (register int i2 = i2_57470_lb; i2 <= i2_87983_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                      }
                    }
                    if (ii1 == 3 && ii2 <= 48) {
                      const int i1_52420_lb = max(i0 - 383, -i0 + 607), i1_94818_ub = i0 - 373;
                      for (register int i1 = i1_52420_lb; i1 < i1_94818_ub; i1 += 1) {
                        if (ii2 == 0) {
                          const int i2_43495_lb = i0 + i1 - 607, i2_40319_ub = i0 - i1 - 372;
                          for (register int i2 = i2_43495_lb; i2 < i2_40319_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                        }
                        const int i2_90413_lb = max(max(i0 + i1 - 595, 32 * ii2 - i0 + i1 + 372), 32 * ii2 - i0 - i1 + 607), i2_38439_ub = 32 * ii2;
                        for (register int i2 = i2_90413_lb; i2 < i2_38439_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                        const int i2_82199_lb = max(32 * ii2, i0 - i1 - 372), i2_25936_ub = 32 * ii2 - i0 + i1 + 403;
                        for (register int i2 = i2_82199_lb; i2 <= i2_25936_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                        const int i2_40875_lb = 32 * ii2 - i0 + i1 + 404, i2_97286_ub = 32 * ii2 - i0 - i1 + 638;
                        for (register int i2 = i2_40875_lb; i2 <= i2_97286_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                      }
                      const int i1_93292_lb = i0 - 373, i1_27678_ub = i0 - 371;
                      for (register int i1 = i1_93292_lb; i1 < i1_27678_ub; i1 += 1) {
                        if (ii2 == 0 && i0 == 490 && i1 == 117) {
                          S1(490, 117, 0);
                        } else if (i0 >= 491 && 16 * ii2 + 483 >= i0 && i1 + 373 == i0) {
                          S1(i0, i0 - 373, 32 * ii2 - 1);
                        }
                        const int i2_99745_lb = max(32 * ii2, i0 - i1 - 372), i2_76880_ub = 32 * ii2 - i0 + i1 + 403;
                        for (register int i2 = i2_99745_lb; i2 <= i2_76880_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                        if (i0 == 490 && i1 == 117) {
                          S1(490, 117, 32 * ii2 + 31);
                        }
                      }
                    } else if (ii1 == 2 && ii2 >= 1 && ii2 <= 48) {
                      const int i1_21321_lb = max(i0 - 415, -i0 + 575), i1_92865_ub = 85;
                      for (register int i1 = i1_21321_lb; i1 <= i1_92865_ub; i1 += 1) {
                        const int i2_5834_lb = 32 * ii2 - i0 - i1 + 575, i2_26226_ub = 32 * ii2 - i0 - i1 + 606;
                        for (register int i2 = i2_5834_lb; i2 <= i2_26226_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                      }
                      const int i1_61386_lb = 86, i1_35863_ub = i0 - 403;
                      for (register int i1 = i1_61386_lb; i1 < i1_35863_ub; i1 += 1) {
                        const int i2_36428_lb = 32 * ii2 - i0 + i1 + 404, i2_60670_ub = 32 * ii2 - i0 + i1 + 435;
                        for (register int i2 = i2_36428_lb; i2 <= i2_60670_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                      }
                    } else if (ii1 == 2 && ii2 == 49) {
                      const int i1_55409_lb = max(i0 - 415, -i0 + 575), i1_29354_ub = i0 - 403;
                      for (register int i1 = i1_55409_lb; i1 < i1_29354_ub; i1 += 1) {
                        const int i2_49457_lb = 0, i2_417_ub = min(i0 + i1 - 575, i0 - i1 - 404);
                        for (register int i2 = i2_49457_lb; i2 < i2_417_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                        const int i2_54584_lb = max(-i0 + i1 + 1972, -i0 - i1 + 2143), i2_27324_ub = 1599;
                        for (register int i2 = i2_54584_lb; i2 <= i2_27324_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                      }
                    }
                  }
                  if (ii1 == 3 && ii2 == 49) {
                    const int i1_41083_lb = -i0 + 595, i1_28406_ub = i0 - 383;
                    for (register int i1 = i1_41083_lb; i1 < i1_28406_ub; i1 += 1) {
                      const int i2_31659_lb = 0, i2_9855_ub = min(i0 + i1 - 595, i0 - i1 - 384);
                      for (register int i2 = i2_31659_lb; i2 < i2_9855_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                      const int i2_23225_lb = i0 - i1 - 384, i2_75155_ub = i0 + i1 - 607;
                      for (register int i2 = i2_23225_lb; i2 < i2_75155_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                      const int i2_66526_lb = -i0 - i1 + 2175, i2_13638_ub = -i0 + i1 + 1951;
                      for (register int i2 = i2_66526_lb; i2 <= i2_13638_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                      const int i2_13594_lb = max(-i0 + i1 + 1952, -i0 - i1 + 2163), i2_48725_ub = 1599;
                      for (register int i2 = i2_13594_lb; i2 <= i2_48725_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                    const int i1_39574_lb = max(i0 - 383, -i0 + 607), i1_70822_ub = i0 - 371;
                    for (register int i1 = i1_39574_lb; i1 < i1_70822_ub; i1 += 1) {
                      const int i2_62364_lb = 0, i2_49218_ub = min(i0 + i1 - 607, i0 - i1 - 372);
                      for (register int i2 = i2_62364_lb; i2 < i2_49218_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                      const int i2_98500_lb = max(-i0 + i1 + 1940, -i0 - i1 + 2175), i2_78461_ub = 1599;
                      for (register int i2 = i2_98500_lb; i2 <= i2_78461_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                  }
                }
                if (ii1 == 2 && ii2 == 0) {
                  for (register int i0 = 497; i0 <= 499; i0 += 1) {
                    const int i1_71327_lb = -i0 + 563, i1_48285_ub = i0 - 403;
                    for (register int i1 = i1_71327_lb; i1 < i1_48285_ub; i1 += 1) {
                      const int i2_62399_lb = i0 - i1 - 404, i2_32713_ub = i0 + i1 - 575;
                      for (register int i2 = i2_62399_lb; i2 < i2_32713_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                      const int i2_84148_lb = max(i0 + i1 - 575, i0 - i1 - 416), i2_15180_ub = i0 + i1 - 563;
                      for (register int i2 = i2_84148_lb; i2 < i2_15180_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                      const int i2_9735_lb = i0 + i1 - 563, i2_55909_ub = -i0 - i1 + 594;
                      for (register int i2 = i2_9735_lb; i2 <= i2_55909_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                      const int i2_60886_lb = max(i0 + i1 - 563, -i0 - i1 + 595), i2_59192_ub = min(-i0 + i1 + 447, -i0 - i1 + 606);
                      for (register int i2 = i2_60886_lb; i2 <= i2_59192_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                      const int i2_56326_lb = max(i0 + i1 - 563, -i0 - i1 + 607), i2_31822_ub = -i0 + i1 + 435;
                      for (register int i2 = i2_56326_lb; i2 <= i2_31822_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                  }
                }
                if (ii1 >= 48 && ii2 >= 1) {
                  if (ii1 == 49 && ii2 == 49) {
                    for (register int i1 = 1577; i1 <= 1578; i1 += 1) {
                      for (register int i2 = 1568; i2 <= 1599; i2 += 1) {
                        S1(490, i1, i2);
                      }
                    }
                    for (register int i1 = 1589; i1 <= 1590; i1 += 1) {
                      for (register int i2 = 1568; i2 <= 1599; i2 += 1) {
                        S1(490, i1, i2);
                      }
                    }
                  } else if (ii1 == 48 && ii2 == 49) {
                    for (register int i1 = 1545; i1 <= 1546; i1 += 1) {
                      const int i2_71573_lb = max(i1 + 22, -i1 + 3113), i2_53018_ub = 1599;
                      for (register int i2 = i2_71573_lb; i2 <= i2_53018_ub; i2 += 1) {
                        S1(490, i1, i2);
                      }
                    }
                    for (register int i1 = 1557; i1 <= 1558; i1 += 1) {
                      for (register int i2 = 1568; i2 <= 1599; i2 += 1) {
                        S1(490, i1, i2);
                      }
                    }
                  }
                  const int i0_12399_lb = max(max(max(16 * ii1 - 294, -16 * ii1 + 1258), 16 * ii1 + ii2 - 342), -16 * ii1 + ii2 + 1210), i0_44686_ub = 499;
                  for (register int i0 = i0_12399_lb; i0 <= i0_44686_ub; i0 += 1) {
                    if (ii2 <= 48) {
                      if (ii1 == 48) {
                        const int i1_94181_lb = -i0 + 2035, i1_83726_ub = i0 + 1056;
                        for (register int i1 = i1_94181_lb; i1 <= i1_83726_ub; i1 += 1) {
                          const int i2_92971_lb = 32 * ii2 - i0 - i1 + 2047, i2_72932_ub = 32 * ii2 - i0 + i1 - 1056;
                          for (register int i2 = i2_92971_lb; i2 < i2_72932_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                          const int i2_32791_lb = max(32 * ii2 - i0 + i1 - 1056, 32 * ii2 - i0 - i1 + 2035), i2_77120_ub = 32 * ii2 - i0 - i1 + 2066;
                          for (register int i2 = i2_32791_lb; i2 <= i2_77120_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                          const int i2_88112_lb = 32 * ii2 - i0 - i1 + 2067, i2_42526_ub = min(32 * ii2 - i0 + i1 - 1025, 32 * ii2 - i0 - i1 + 2078);
                          for (register int i2 = i2_88112_lb; i2 <= i2_42526_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                        }
                      }
                      const int i1_49381_lb = max(max(32 * ii1 - i0 + 499, i0 + 1057), -i0 + 2047), i1_48999_ub = min(32 * ii1 + i0 - 468, i0 + 1088);
                      for (register int i1 = i1_49381_lb; i1 <= i1_48999_ub; i1 += 1) {
                        const int i2_1718_lb = 32 * ii2 - i0 - i1 + 2079, i2_5708_ub = 32 * ii2 - i0 + i1 - 1088;
                        for (register int i2 = i2_1718_lb; i2 < i2_5708_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                        const int i2_80821_lb = max(32 * ii2 - i0 + i1 - 1068, 32 * ii2 - i0 - i1 + 2047), i2_20939_ub = 32 * ii2;
                        for (register int i2 = i2_80821_lb; i2 < i2_20939_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                        const int i2_35822_lb = max(32 * ii2 - i0 + i1 - 1088, 32 * ii2 - i0 - i1 + 2067), i2_57402_ub = 32 * ii2;
                        for (register int i2 = i2_35822_lb; i2 < i2_57402_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                        const int i2_55467_lb = 32 * ii2, i2_75791_ub = 32 * ii1 + 32 * ii2 - i0 - i1 + 530;
                        for (register int i2 = i2_55467_lb; i2 <= i2_75791_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                        const int i2_40856_lb = 32 * ii1 + 32 * ii2 - i0 - i1 + 531, i2_97854_ub = min(32 * ii2 - i0 + i1 - 1057, 32 * ii2 - i0 - i1 + 2110);
                        for (register int i2 = i2_40856_lb; i2 <= i2_97854_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                        const int i2_82287_lb = max(32 * ii2 - i0 + i1 - 1056, 32 * ii2 - i0 - i1 + 2067), i2_54300_ub = 32 * ii2 - i0 - i1 + 2078;
                        for (register int i2 = i2_82287_lb; i2 <= i2_54300_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                        if (ii1 == 48) {
                          const int i2_53835_lb = 32 * ii2 - i0 - i1 + 2079, i2_70213_ub = 32 * ii2 - i0 + i1 - 1036;
                          for (register int i2 = i2_53835_lb; i2 < i2_70213_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                        }
                      }
                      if (ii1 == 49) {
                        const int i1_23670_lb = max(i0 + 1089, -i0 + 2079), i1_80639_ub = i0 + 1100;
                        for (register int i1 = i1_23670_lb; i1 <= i1_80639_ub; i1 += 1) {
                          const int i2_4150_lb = max(32 * ii2 - i0 + i1 - 1100, 32 * ii2 - i0 - i1 + 2079), i2_42258_ub = 32 * ii2 - i0 + i1 - 1068;
                          for (register int i2 = i2_4150_lb; i2 < i2_42258_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                          const int i2_38646_lb = 32 * ii2 - i0 + i1 - 1068, i2_16549_ub = 32 * ii2 - i0 - i1 + 2110;
                          for (register int i2 = i2_38646_lb; i2 <= i2_16549_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                        }
                      }
                    } else {
                      if (ii1 == 48) {
                        const int i1_3297_lb = -i0 + 2035, i1_32827_ub = min(i0 + 1056, -i0 + 2047);
                        for (register int i1 = i1_3297_lb; i1 <= i1_32827_ub; i1 += 1) {
                          const int i2_16627_lb = 0, i2_96268_ub = min(i0 + i1 - 2036, i0 - i1 + 1055);
                          for (register int i2 = i2_16627_lb; i2 <= i2_96268_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                          const int i2_5760_lb = max(-i0 + i1 + 512, -i0 - i1 + 3603), i2_49418_ub = 1599;
                          for (register int i2 = i2_5760_lb; i2 <= i2_49418_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                        }
                        const int i1_89740_lb = -i0 + 2048, i1_10224_ub = i0 + 1056;
                        for (register int i1 = i1_89740_lb; i1 <= i1_10224_ub; i1 += 1) {
                          const int i2_91944_lb = 0, i2_39122_ub = i0 - i1 + 1055;
                          for (register int i2 = i2_91944_lb; i2 <= i2_39122_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                          const int i2_75575_lb = i0 - i1 + 1056, i2_10015_ub = i0 + i1 - 2047;
                          for (register int i2 = i2_75575_lb; i2 < i2_10015_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                          const int i2_61182_lb = -i0 - i1 + 3615, i2_72749_ub = -i0 + i1 + 511;
                          for (register int i2 = i2_61182_lb; i2 <= i2_72749_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                          const int i2_30954_lb = -i0 + i1 + 512, i2_13356_ub = 1599;
                          for (register int i2 = i2_30954_lb; i2 <= i2_13356_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                        }
                      }
                      const int i1_30151_lb = max(max(32 * ii1 - i0 + 499, i0 + 1057), -i0 + 2047), i1_2773_ub = min(32 * ii1 + i0 - 468, i0 + 1088);
                      for (register int i1 = i1_30151_lb; i1 <= i1_2773_ub; i1 += 1) {
                        const int i2_89147_lb = 0, i2_71008_ub = i0 + i1 - 2079;
                        for (register int i2 = i2_89147_lb; i2 < i2_71008_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                        const int i2_627_lb = max(0, i0 + i1 - 2079), i2_87787_ub = min(min(min(i0 + i1 - 2048, -32 * ii1 + i0 + i1 - 500), 32 * ii1 + i0 - i1 - 469), i0 - i1 + 1087);
                        for (register int i2 = i2_627_lb; i2 <= i2_87787_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                        const int i2_41660_lb = -i0 - i1 + 3647, i2_54462_ub = -i0 + i1 + 479;
                        for (register int i2 = i2_41660_lb; i2 <= i2_54462_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                        const int i2_58000_lb = max(max(max(-i0 + i1 + 480, -32 * ii1 - i0 + i1 + 2036), 32 * ii1 - i0 - i1 + 2067), -i0 - i1 + 3615), i2_65331_ub = 1599;
                        for (register int i2 = i2_58000_lb; i2 <= i2_65331_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                      }
                      if (ii1 == 49) {
                        const int i1_51453_lb = max(i0 + 1089, -i0 + 2079), i1_78502_ub = i0 + 1100;
                        for (register int i1 = i1_51453_lb; i1 <= i1_78502_ub; i1 += 1) {
                          const int i2_23941_lb = 0, i2_90100_ub = min(i0 + i1 - 2080, i0 - i1 + 1099);
                          for (register int i2 = i2_23941_lb; i2 <= i2_90100_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                          const int i2_11404_lb = max(-i0 + i1 + 468, -i0 - i1 + 3647), i2_27238_ub = 1599;
                          for (register int i2 = i2_11404_lb; i2 <= i2_27238_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                        }
                      }
                    }
                  }
                } else if (ii1 >= 4 && ii1 <= 47 && ii2 >= 1 && ii2 <= 48) {
                  if (ii1 == 47 && ii2 >= 2) {
                    for (register int i1 = 1513; i1 <= 1514; i1 += 1) {
                      const int i2_39859_lb = 32 * ii2, i2_61391_ub = 32 * ii2 + 31;
                      for (register int i2 = i2_39859_lb; i2 <= i2_61391_ub; i2 += 1) {
                        S1(490, i1, i2);
                      }
                    }
                    for (register int i1 = 1525; i1 <= 1526; i1 += 1) {
                      const int i2_71616_lb = 32 * ii2, i2_2098_ub = 32 * ii2 + 31;
                      for (register int i2 = i2_71616_lb; i2 <= i2_2098_ub; i2 += 1) {
                        S1(490, i1, i2);
                      }
                    }
                  } else if (ii1 == 47 && ii2 == 1) {
                    for (register int i1 = 1513; i1 <= 1514; i1 += 1) {
                      for (register int i2 = 32; i2 <= 63; i2 += 1) {
                        S1(490, i1, i2);
                      }
                    }
                    for (register int i1 = 1525; i1 <= 1526; i1 += 1) {
                      for (register int i2 = 32; i2 <= 63; i2 += 1) {
                        S1(490, i1, i2);
                      }
                    }
                  }
                  const int i0_78544_lb = max(490, ii1 + 444), i0_65111_ub = 499;
                  for (register int i0 = i0_78544_lb; i0 <= i0_65111_ub; i0 += 1) {
                    const int i1_37452_lb = 32 * ii1 - i0 + 499, i1_79171_ub = min(32 * ii1 + i0 - 480, 32 * ii1 - i0 + 510);
                    for (register int i1 = i1_37452_lb; i1 <= i1_79171_ub; i1 += 1) {
                      const int i2_69250_lb = max(-32 * ii1 + 32 * ii2 - i0 + i1 + 480, 32 * ii1 + 32 * ii2 - i0 - i1 + 499), i2_95464_ub = 32 * ii1 + 32 * ii2 - i0 - i1 + 530;
                      for (register int i2 = i2_69250_lb; i2 <= i2_95464_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                      const int i2_49986_lb = 32 * ii1 + 32 * ii2 - i0 - i1 + 531, i2_43602_ub = -32 * ii1 + 32 * ii2 - i0 + i1 + 511;
                      for (register int i2 = i2_49986_lb; i2 <= i2_43602_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                    if (i0 == 490) {
                      const int i2_77147_lb = 32 * ii2, i2_17791_ub = 32 * ii2 + 31;
                      for (register int i2 = i2_77147_lb; i2 <= i2_17791_ub; i2 += 1) {
                        S1(490, 32 * ii1 + 21, i2);
                      }
                    }
                    const int i1_22104_lb = 32 * ii1 - i0 + 511, i1_1089_ub = min(32 * ii1 + i0 - 481, -i0 + 2014);
                    for (register int i1 = i1_22104_lb; i1 <= i1_1089_ub; i1 += 1) {
                      const int i2_7891_lb = -32 * ii1 + 32 * ii2 - i0 + i1 + 480, i2_33508_ub = 32 * ii1 + 32 * ii2 - i0 - i1 + 510;
                      for (register int i2 = i2_7891_lb; i2 <= i2_33508_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                      const int i2_44679_lb = 32 * ii1 + 32 * ii2 - i0 - i1 + 511, i2_63523_ub = 32 * ii2;
                      for (register int i2 = i2_44679_lb; i2 < i2_63523_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                      const int i2_77892_lb = 32 * ii2, i2_84538_ub = min(-32 * ii1 + 32 * ii2 - i0 + i1 + 511, 32 * ii1 + 32 * ii2 - i0 - i1 + 542);
                      for (register int i2 = i2_77892_lb; i2 <= i2_84538_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                    const int i1_24914_lb = max(32 * ii1 + i0 - 480, 32 * ii1 - i0 + 511), i1_71694_ub = min(min(32 * ii1 + 21, 32 * ii1 + i0 - 470), -i0 + 2014);
                    for (register int i1 = i1_24914_lb; i1 <= i1_71694_ub; i1 += 1) {
                      const int i2_14138_lb = 32 * ii1 + 32 * ii2 - i0 - i1 + 511, i2_96530_ub = 32 * ii1 + 32 * ii2 - i0 - i1 + 542;
                      for (register int i2 = i2_14138_lb; i2 <= i2_96530_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                    const int i1_73792_lb = 32 * ii1 + 22, i1_15563_ub = min(32 * ii1 + i0 - 468, -i0 + 2014);
                    for (register int i1 = i1_73792_lb; i1 <= i1_15563_ub; i1 += 1) {
                      const int i2_60074_lb = -32 * ii1 + 32 * ii2 - i0 + i1 + 468, i2_2258_ub = -32 * ii1 + 32 * ii2 - i0 + i1 + 499;
                      for (register int i2 = i2_60074_lb; i2 <= i2_2258_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                    if (ii1 == 47) {
                      const int i1_78171_lb = -i0 + 2015, i1_80014_ub = i0 + 1036;
                      for (register int i1 = i1_78171_lb; i1 <= i1_80014_ub; i1 += 1) {
                        const int i2_61677_lb = 32 * ii2 - i0 + i1 - 1024, i2_54134_ub = 32 * ii2 - i0 - i1 + 2014;
                        for (register int i2 = i2_61677_lb; i2 <= i2_54134_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                        const int i2_46458_lb = max(32 * ii2 - i0 + i1 - 1036, 32 * ii2 - i0 - i1 + 2015), i2_56574_ub = min(32 * ii2 - i0 + i1 - 993, 32 * ii2 - i0 - i1 + 2046);
                        for (register int i2 = i2_46458_lb; i2 <= i2_56574_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                        const int i2_35597_lb = 32 * ii2 - i0 - i1 + 2047, i2_262_ub = 32 * ii2 - i0 + i1 - 1004;
                        for (register int i2 = i2_35597_lb; i2 < i2_262_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                      }
                    }
                  }
                } else if (ii2 == 0) {
                  const int i0_35745_lb = max(max(485, 16 * ii1 - 262), -16 * ii1 + 538), i0_4847_ub = min(499, 16 * ii1 + 441);
                  for (register int i0 = i0_35745_lb; i0 <= i0_4847_ub; i0 += 1) {
                    const int i1_12079_lb = 32 * ii1 - i0 + 499, i1_2083_ub = 32 * ii1 + 9;
                    for (register int i1 = i1_12079_lb; i1 <= i1_2083_ub; i1 += 1) {
                      const int i2_64801_lb = -32 * ii1 + i0 + i1 - 499, i2_89226_ub = 32 * ii1 - i0 - i1 + 530;
                      for (register int i2 = i2_64801_lb; i2 <= i2_89226_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                    const int i1_19875_lb = max(32 * ii1 + 10, 32 * ii1 - i0 + 500), i1_86906_ub = min(min(min(32 * ii1 + 21, 32 * ii1 + i0 - 469), i0 + 1024), -i0 + 2014);
                    for (register int i1 = i1_19875_lb; i1 <= i1_86906_ub; i1 += 1) {
                      if (((-i0 - i1 + 2014) % 32) + 2 * i0 >= 990) {
                        const int i2_90315_lb = max(-32 * ii1 + i0 + i1 - 511, 32 * ii1 + i0 - i1 - 480), i2_44118_ub = min(-32 * ii1 + i0 + i1 - 499, 32 * ii1 + i0 - i1 - 468);
                        for (register int i2 = i2_90315_lb; i2 < i2_44118_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                        if (32 * ii1 + 510 >= i0 + i1) {
                          const int i2_36766_lb = -32 * ii1 + i0 + i1 - 499, i2_34995_ub = 32 * ii1 - i0 - i1 + 530;
                          for (register int i2 = i2_36766_lb; i2 <= i2_34995_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                        } else {
                          const int i2_7641_lb = 32 * ii1 + i0 - i1 - 468, i2_14658_ub = -32 * ii1 + i0 + i1 - 499;
                          for (register int i2 = i2_7641_lb; i2 < i2_14658_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                          const int i2_35885_lb = -32 * ii1 + i0 + i1 - 499, i2_32556_ub = 32 * ii1 - i0 - i1 + 530;
                          for (register int i2 = i2_35885_lb; i2 <= i2_32556_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                          const int i2_2704_lb = max(-32 * ii1 + i0 + i1 - 499, 32 * ii1 - i0 - i1 + 531), i2_66375_ub = min(32 * ii1 + i0 - i1 - 440, -32 * ii1 - i0 + i1 + 499);
                          for (register int i2 = i2_2704_lb; i2 <= i2_66375_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                        }
                        const int i2_45438_lb = max(-32 * ii1 - i0 + i1 + 500, 32 * ii1 - i0 - i1 + 531), i2_76497_ub = min(min(32 * ii1 + i0 - i1 - 440, -32 * ii1 - i0 + i1 + 511), 32 * ii1 - i0 - i1 + 542);
                        for (register int i2 = i2_45438_lb; i2 <= i2_76497_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                        if (i0 == 490 && i1 == 32 * ii1 + 21) {
                          for (register int i2 = 30; i2 <= 31; i2 += 1) {
                            S1(490, 32 * ii1 + 21, i2);
                          }
                        }
                      }
                    }
                    const int i1_95107_lb = 32 * ii1 + 22, i1_60110_ub = min(min(32 * ii1 + i0 - 468, 32 * ii1 - i0 + 514), i0 + 1024);
                    for (register int i1 = i1_95107_lb; i1 <= i1_60110_ub; i1 += 1) {
                      const int i2_1879_lb = 32 * ii1 + i0 - i1 - 468, i2_73136_ub = -32 * ii1 - i0 + i1 + 499;
                      for (register int i2 = i2_1879_lb; i2 <= i2_73136_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                    const int i1_14244_lb = max(32 * ii1 + 22, 32 * ii1 - i0 + 515), i1_64689_ub = min(32 * ii1 + i0 - 468, -i0 + 2014);
                    for (register int i1 = i1_14244_lb; i1 <= i1_64689_ub; i1 += 1) {
                      const int i2_29710_lb = 32 * ii1 + i0 - i1 - 468, i2_49842_ub = -32 * ii1 - i0 + i1 + 499;
                      for (register int i2 = i2_29710_lb; i2 <= i2_49842_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                    if (ii1 == 47) {
                      const int i1_81304_lb = -i0 + 2015, i1_65456_ub = i0 + 1036;
                      for (register int i1 = i1_81304_lb; i1 <= i1_65456_ub; i1 += 1) {
                        const int i2_71041_lb = i0 - i1 + 1036, i2_93383_ub = i0 + i1 - 2015;
                        for (register int i2 = i2_71041_lb; i2 < i2_93383_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                        const int i2_83891_lb = max(i0 + i1 - 2015, i0 - i1 + 1024), i2_35843_ub = i0 + i1 - 2003;
                        for (register int i2 = i2_83891_lb; i2 < i2_35843_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                        const int i2_82609_lb = i0 + i1 - 2003, i2_3766_ub = min(-i0 + i1 - 993, -i0 - i1 + 2046);
                        for (register int i2 = i2_82609_lb; i2 <= i2_3766_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                        const int i2_39101_lb = max(i0 + i1 - 2003, -i0 - i1 + 2047), i2_89277_ub = -i0 + i1 - 1004;
                        for (register int i2 = i2_39101_lb; i2 < i2_89277_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                      }
                    }
                  }
                  if (ii1 >= 48) {
                    if (ii1 == 48) {
                      for (register int i0 = 490; i0 <= 493; i0 += 1) {
                        const int i1_40624_lb = -i0 + 2035, i1_71878_ub = i0 + 1056;
                        for (register int i1 = i1_40624_lb; i1 <= i1_71878_ub; i1 += 1) {
                          const int i2_6878_lb = i0 - i1 + 1056, i2_76509_ub = i0 + i1 - 2035;
                          for (register int i2 = i2_6878_lb; i2 < i2_76509_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                          const int i2_20786_lb = i0 + i1 - 2035, i2_9582_ub = -i0 - i1 + 2066;
                          for (register int i2 = i2_20786_lb; i2 <= i2_9582_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                          const int i2_42885_lb = -i0 - i1 + 2067, i2_82577_ub = -i0 + i1 - 1024;
                          for (register int i2 = i2_42885_lb; i2 < i2_82577_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                        }
                        const int i1_2431_lb = -i0 + 2047, i1_24824_ub = i0 + 1068;
                        for (register int i1 = i1_2431_lb; i1 <= i1_24824_ub; i1 += 1) {
                          const int i2_88089_lb = i0 - i1 + 1068, i2_97538_ub = i0 + i1 - 2047;
                          for (register int i2 = i2_88089_lb; i2 < i2_97538_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                          const int i2_1286_lb = i0 + i1 - 2047, i2_6320_ub = -i0 - i1 + 2078;
                          for (register int i2 = i2_1286_lb; i2 <= i2_6320_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                          const int i2_70675_lb = -i0 - i1 + 2079, i2_15530_ub = -i0 + i1 - 1036;
                          for (register int i2 = i2_70675_lb; i2 < i2_15530_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                        }
                      }
                    }
                    const int i0_87362_lb = max(16 * ii1 - 294, -16 * ii1 + 1262), i0_385_ub = 499;
                    for (register int i0 = i0_87362_lb; i0 <= i0_385_ub; i0 += 1) {
                      if (ii1 == 48) {
                        const int i1_81724_lb = -i0 + 2035, i1_68666_ub = i0 + 1056;
                        for (register int i1 = i1_81724_lb; i1 <= i1_68666_ub; i1 += 1) {
                          const int i2_82193_lb = max(i0 + i1 - 2047, i0 - i1 + 1056), i2_52766_ub = i0 + i1 - 2035;
                          for (register int i2 = i2_82193_lb; i2 < i2_52766_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                          const int i2_62049_lb = i0 + i1 - 2035, i2_66085_ub = -i0 - i1 + 2066;
                          for (register int i2 = i2_62049_lb; i2 <= i2_66085_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                          const int i2_88609_lb = max(i0 + i1 - 2035, -i0 - i1 + 2067), i2_61010_ub = min(-i0 + i1 - 1025, -i0 - i1 + 2078);
                          for (register int i2 = i2_88609_lb; i2 <= i2_61010_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                        }
                        const int i1_86203_lb = max(i0 + 1057, -i0 + 2047), i1_27710_ub = 1561;
                        for (register int i1 = i1_86203_lb; i1 <= i1_27710_ub; i1 += 1) {
                          const int i2_66639_lb = i0 - i1 + 1068, i2_50440_ub = i0 + i1 - 2047;
                          for (register int i2 = i2_66639_lb; i2 < i2_50440_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                          const int i2_19929_lb = i0 + i1 - 2047, i2_7263_ub = i0 + i1 - 2035;
                          for (register int i2 = i2_19929_lb; i2 < i2_7263_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                          const int i2_38671_lb = i0 + i1 - 2035, i2_26807_ub = -i0 - i1 + 2078;
                          for (register int i2 = i2_38671_lb; i2 <= i2_26807_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                          const int i2_125_lb = max(i0 + i1 - 2035, -i0 - i1 + 2079), i2_59457_ub = -i0 + i1 - 1036;
                          for (register int i2 = i2_125_lb; i2 < i2_59457_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                        }
                      } else {
                        const int i1_36390_lb = -i0 + 2067, i1_43010_ub = 1577;
                        for (register int i1 = i1_36390_lb; i1 <= i1_43010_ub; i1 += 1) {
                          const int i2_42034_lb = i0 + i1 - 2067, i2_38821_ub = -i0 - i1 + 2098;
                          for (register int i2 = i2_42034_lb; i2 <= i2_38821_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                        }
                      }
                      const int i1_84186_lb = 16 * ii1 + 794, i1_46476_ub = min(32 * ii1 + i0 - 468, i0 + 1088);
                      for (register int i1 = i1_84186_lb; i1 <= i1_46476_ub; i1 += 1) {
                        if (ii1 == 48) {
                          const int i2_36360_lb = i0 - i1 + 1068, i2_85472_ub = i0 - i1 + 1087;
                          for (register int i2 = i2_36360_lb; i2 <= i2_85472_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                        }
                        const int i2_69148_lb = max(i0 + i1 - 2079, i0 - i1 + 1088), i2_23387_ub = -32 * ii1 + i0 + i1 - 499;
                        for (register int i2 = i2_69148_lb; i2 < i2_23387_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                        if (ii1 == 49) {
                          const int i2_1002_lb = i0 + i1 - 2067, i2_56510_ub = min(-i0 + i1 - 1057, -i0 - i1 + 2110);
                          for (register int i2 = i2_1002_lb; i2 <= i2_56510_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                        } else {
                          const int i2_40124_lb = i0 + i1 - 2035, i2_82727_ub = -i0 + i1 - 1036;
                          for (register int i2 = i2_40124_lb; i2 < i2_82727_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                        }
                      }
                      if (ii1 == 49) {
                        const int i1_25176_lb = max(i0 + 1089, -i0 + 2079), i1_22318_ub = i0 + 1100;
                        for (register int i1 = i1_25176_lb; i1 <= i1_22318_ub; i1 += 1) {
                          const int i2_51845_lb = i0 + i1 - 2079, i2_3577_ub = i0 - i1 + 1099;
                          for (register int i2 = i2_51845_lb; i2 <= i2_3577_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                          const int i2_4755_lb = i0 - i1 + 1100, i2_56806_ub = -i0 + i1 - 1068;
                          for (register int i2 = i2_4755_lb; i2 < i2_56806_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                          const int i2_80940_lb = -i0 + i1 - 1068, i2_7310_ub = -i0 - i1 + 2110;
                          for (register int i2 = i2_80940_lb; i2 <= i2_7310_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                        }
                      }
                    }
                  }
                } else if (ii2 == 49) {
                  const int i0_868_lb = max(max(484, 16 * ii1 - 262), -16 * ii1 + 538), i0_47579_ub = min(499, 16 * ii1 + 441);
                  for (register int i0 = i0_868_lb; i0 <= i0_47579_ub; i0 += 1) {
                    const int i1_57751_lb = 32 * ii1 - i0 + 499, i1_20797_ub = min(min(32 * ii1 + i0 - 468, i0 + 1024), -i0 + 2014);
                    for (register int i1 = i1_57751_lb; i1 <= i1_20797_ub; i1 += 1) {
                      if (((-i0 - i1 + 2014) % 32) + 2 * i0 >= 990) {
                        const int i2_71195_lb = 0, i2_96422_ub = min(-32 * ii1 + i0 + i1 - 499, 32 * ii1 + i0 - i1 - 468);
                        for (register int i2 = i2_71195_lb; i2 < i2_96422_ub; i2 += 1) {
                          if (2 * i0 >= ((i0 + i1 - i2) % 32) + 2 * i2 + 961) {
                            S1(i0, i1, i2);
                          }
                        }
                        const int i2_47605_lb = max(-32 * ii1 - i0 + i1 + 2036, 32 * ii1 - i0 - i1 + 2067), i2_71320_ub = 1599;
                        for (register int i2 = i2_47605_lb; i2 <= i2_71320_ub; i2 += 1) {
                          if (i2 >= 1568 || 2 * i0 + 2 * i2 >= ((i0 + i1 + i2 + 1) % 32) + 4095) {
                            S1(i0, i1, i2);
                          }
                        }
                      }
                    }
                    if (ii1 == 47) {
                      const int i1_72231_lb = -i0 + 2015, i1_83995_ub = i0 + 1036;
                      for (register int i1 = i1_72231_lb; i1 <= i1_83995_ub; i1 += 1) {
                        const int i2_14330_lb = 0, i2_30618_ub = min(i0 + i1 - 2016, i0 - i1 + 1035);
                        for (register int i2 = i2_14330_lb; i2 <= i2_30618_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                        const int i2_39168_lb = i0 + i1 - 2015, i2_14868_ub = i0 - i1 + 1023;
                        for (register int i2 = i2_39168_lb; i2 <= i2_14868_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                        const int i2_77094_lb = -i0 + i1 + 544, i2_91880_ub = -i0 - i1 + 3582;
                        for (register int i2 = i2_77094_lb; i2 <= i2_91880_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                        const int i2_340_lb = max(-i0 + i1 + 532, -i0 - i1 + 3583), i2_62594_ub = 1599;
                        for (register int i2 = i2_340_lb; i2 <= i2_62594_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                      }
                    }
                  }
                }
              } else if (ii1 >= 47) {
                if (ii1 == 47 && ii2 == 49) {
                  for (register int i1 = 1535; i1 <= 1536; i1 += 1) {
                    const int i2_19105_lb = max(i1 + 32, -i1 + 3103), i2_71744_ub = 1599;
                    for (register int i2 = i2_19105_lb; i2 <= i2_71744_ub; i2 += 1) {
                      S1(480, i1, i2);
                    }
                  }
                }
                const int i0_16773_lb = max(480, -ii1 + ii2 + 479), i0_60633_ub = min(498, 8 * ii2 + 496);
                for (register int i0 = i0_16773_lb; i0 <= i0_60633_ub; i0 += 1) {
                  if (ii1 == 48) {
                    if (ii2 == 0 && i0 == 496) {
                      for (register int i1 = 1565; i1 <= 1566; i1 += 1) {
                        const int i2_80563_lb = i1 - 1551, i2_15169_ub = -i1 + 1582;
                        for (register int i2 = i2_80563_lb; i2 <= i2_15169_ub; i2 += 1) {
                          S1(496, i1, i2);
                        }
                      }
                    } else if (ii2 == 49) {
                      const int i1_41776_lb = max(i0 + 1069, -i0 + 2047), i1_77855_ub = 1567;
                      for (register int i1 = i1_41776_lb; i1 <= i1_77855_ub; i1 += 1) {
                        const int i2_22479_lb = 0, i2_42644_ub = i0 + i1 - 2047;
                        for (register int i2 = i2_22479_lb; i2 < i2_42644_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                        const int i2_25434_lb = -i0 - i1 + 3615, i2_96582_ub = 1599;
                        for (register int i2 = i2_25434_lb; i2 <= i2_96582_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                      }
                    } else {
                      const int i1_63442_lb = max(i0 + 1069, -i0 + 2047), i1_96629_ub = 1567;
                      for (register int i1 = i1_63442_lb; i1 <= i1_96629_ub; i1 += 1) {
                        if (ii2 == 1) {
                          const int i2_9356_lb = max(i0 + i1 - 2047, -i0 - i1 + 2079), i2_11047_ub = i0 - i1 + 1087;
                          for (register int i2 = i2_9356_lb; i2 <= i2_11047_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                        } else if (ii2 == 0) {
                          const int i2_67949_lb = i0 + i1 - 2047, i2_81588_ub = i0 - i1 + 1087;
                          for (register int i2 = i2_67949_lb; i2 <= i2_81588_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                        }
                        const int i2_11394_lb = max(i0 - i1 + 1088, 32 * ii2 - i0 - i1 + 2047), i2_98631_ub = 32 * ii2 - i0 - i1 + 2078;
                        for (register int i2 = i2_11394_lb; i2 <= i2_98631_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                      }
                    }
                    const int i1_28558_lb = max(1568, -16 * ii2 + i0 + 1073), i1_66914_ub = min(i0 + 1088, -i0 + 2066);
                    for (register int i1 = i1_28558_lb; i1 <= i1_66914_ub; i1 += 1) {
                      if (ii2 <= 48) {
                        const int i2_13499_lb = max(32 * ii2 - i0 + i1 - 1088, i0 - i1 + 1088), i2_22004_ub = 32 * ii2 - i0 + i1 - 1056;
                        for (register int i2 = i2_13499_lb; i2 < i2_22004_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                      } else {
                        const int i2_58795_lb = 0, i2_30191_ub = i0 - i1 + 1087;
                        for (register int i2 = i2_58795_lb; i2 <= i2_30191_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                        const int i2_84598_lb = -i0 + i1 + 480, i2_90414_ub = 1599;
                        for (register int i2 = i2_84598_lb; i2 <= i2_90414_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                      }
                    }
                  } else if (ii2 <= 48) {
                    const int i1_47886_lb = max(i0 + 1037, -i0 + 2015), i1_20055_ub = min(1535, 16 * ii2 - i0 + 2030);
                    for (register int i1 = i1_47886_lb; i1 <= i1_20055_ub; i1 += 1) {
                      const int i2_62158_lb = max(i0 + i1 - 2015, 32 * ii2 - i0 - i1 + 2015), i2_64659_ub = 32 * ii2 - i0 - i1 + 2046;
                      for (register int i2 = i2_62158_lb; i2 <= i2_64659_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                    if (ii2 == 0 && i0 == 496) {
                      for (register int i1 = 1537; i1 <= 1538; i1 += 1) {
                        const int i2_49630_lb = -i1 + 1552, i2_77604_ub = i1 - 1520;
                        for (register int i2 = i2_49630_lb; i2 < i2_77604_ub; i2 += 1) {
                          S1(496, i1, i2);
                        }
                      }
                    } else {
                      const int i1_4093_lb = 1536, i1_91406_ub = min(i0 + 1056, -i0 + 2034);
                      for (register int i1 = i1_4093_lb; i1 <= i1_91406_ub; i1 += 1) {
                        if (ii2 == 0) {
                          const int i2_55459_lb = i0 - i1 + 1056, i2_42925_ub = i0 + i1 - 2015;
                          for (register int i2 = i2_55459_lb; i2 < i2_42925_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                        } else if (ii2 == 1) {
                          const int i2_34051_lb = max(-i0 + i1 - 1024, i0 - i1 + 1056), i2_80893_ub = i0 + i1 - 2015;
                          for (register int i2 = i2_34051_lb; i2 < i2_80893_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                        }
                        const int i2_39507_lb = max(i0 + i1 - 2015, 32 * ii2 - i0 + i1 - 1056), i2_13845_ub = 32 * ii2 - i0 + i1 - 1024;
                        for (register int i2 = i2_39507_lb; i2 < i2_13845_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                      }
                    }
                  } else {
                    const int i1_93875_lb = max(i0 + 1037, -i0 + 2015), i1_65216_ub = min(i0 + 1056, -i0 + 2034);
                    for (register int i1 = i1_93875_lb; i1 <= i1_65216_ub; i1 += 1) {
                      const int i2_24892_lb = 0, i2_78176_ub = min(i0 + i1 - 2016, i0 - i1 + 1055);
                      for (register int i2 = i2_24892_lb; i2 <= i2_78176_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                      const int i2_63156_lb = max(-i0 + i1 + 512, -i0 - i1 + 3583), i2_52638_ub = 1599;
                      for (register int i2 = i2_63156_lb; i2 <= i2_52638_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                  }
                }
              } else {
                const int i0_93160_lb = 480, i0_91714_ub = min(498, 8 * ii2 + 496);
                for (register int i0 = i0_93160_lb; i0 <= i0_91714_ub; i0 += 1) {
                  if (ii1 >= 3) {
                    const int i1_19552_lb = max(32 * ii1 + i0 - 467, 32 * ii1 - i0 + 511), i1_6659_ub = min(32 * ii1 + 31, 32 * ii1 + 16 * ii2 - i0 + 526);
                    for (register int i1 = i1_19552_lb; i1 <= i1_6659_ub; i1 += 1) {
                      if (ii2 <= 48) {
                        const int i2_13718_lb = max(-32 * ii1 + i0 + i1 - 511, 32 * ii1 + 32 * ii2 - i0 - i1 + 511), i2_94699_ub = 32 * ii1 + 32 * ii2 - i0 - i1 + 542;
                        for (register int i2 = i2_13718_lb; i2 <= i2_94699_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                      } else {
                        const int i2_53203_lb = 0, i2_14668_ub = -32 * ii1 + i0 + i1 - 511;
                        for (register int i2 = i2_53203_lb; i2 < i2_14668_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                        const int i2_1466_lb = 32 * ii1 - i0 - i1 + 2079, i2_1089_ub = 1599;
                        for (register int i2 = i2_1466_lb; i2 <= i2_1089_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                      }
                    }
                    if (ii2 == 0 && i0 == 496) {
                      const int i1_34724_lb = max(113, 32 * ii1 + 33), i1_79976_ub = 32 * ii1 + 34;
                      for (register int i1 = i1_34724_lb; i1 <= i1_79976_ub; i1 += 1) {
                        const int i2_82100_lb = 32 * ii1 - i1 + 48, i2_15413_ub = -32 * ii1 + i1 - 16;
                        for (register int i2 = i2_82100_lb; i2 < i2_15413_ub; i2 += 1) {
                          S1(496, i1, i2);
                        }
                      }
                    } else if (ii2 == 49) {
                      const int i1_68901_lb = 32 * ii1 + 32, i1_48082_ub = min(32 * ii1 + i0 - 448, 32 * ii1 - i0 + 530);
                      for (register int i1 = i1_68901_lb; i1 <= i1_48082_ub; i1 += 1) {
                        const int i2_93017_lb = 0, i2_72994_ub = 32 * ii1 + i0 - i1 - 448;
                        for (register int i2 = i2_93017_lb; i2 < i2_72994_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                        const int i2_55841_lb = -32 * ii1 - i0 + i1 + 2016, i2_64828_ub = 1599;
                        for (register int i2 = i2_55841_lb; i2 <= i2_64828_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                      }
                    } else {
                      const int i1_15919_lb = max(32 * ii1 + 32, i0 - 383), i1_6244_ub = min(32 * ii1 + i0 - 448, 32 * ii1 - i0 + 530);
                      for (register int i1 = i1_15919_lb; i1 <= i1_6244_ub; i1 += 1) {
                        const int i2_45721_lb = max(32 * ii1 + i0 - i1 - 448, -32 * ii1 + 32 * ii2 - i0 + i1 + 448), i2_71779_ub = -32 * ii1 + i0 + i1 - 511;
                        for (register int i2 = i2_45721_lb; i2 < i2_71779_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                        const int i2_20089_lb = max(-32 * ii1 + i0 + i1 - 511, -32 * ii1 + 32 * ii2 - i0 + i1 + 448), i2_55948_ub = -32 * ii1 + 32 * ii2 - i0 + i1 + 479;
                        for (register int i2 = i2_20089_lb; i2 <= i2_55948_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                      }
                    }
                  } else if (ii1 == 2 && ii2 <= 48) {
                    if (ii2 == 0 && i0 == 496) {
                      for (register int i1 = 93; i1 <= 94; i1 += 1) {
                        const int i2_50477_lb = i1 - 79, i2_16503_ub = -i1 + 110;
                        for (register int i2 = i2_50477_lb; i2 <= i2_16503_ub; i2 += 1) {
                          S1(496, i1, i2);
                        }
                      }
                    } else {
                      const int i1_30323_lb = max(i0 - 403, -i0 + 575), i1_43637_ub = 95;
                      for (register int i1 = i1_30323_lb; i1 <= i1_43637_ub; i1 += 1) {
                        if (ii2 == 0) {
                          const int i2_8217_lb = i0 + i1 - 575, i2_49875_ub = i0 - i1 - 384;
                          for (register int i2 = i2_8217_lb; i2 < i2_49875_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                        } else if (ii2 == 1) {
                          const int i2_66648_lb = max(i0 + i1 - 575, -i0 - i1 + 607), i2_38287_ub = i0 - i1 - 384;
                          for (register int i2 = i2_66648_lb; i2 < i2_38287_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                        }
                        const int i2_60927_lb = max(i0 - i1 - 384, 32 * ii2 - i0 - i1 + 575), i2_19851_ub = 32 * ii2 - i0 - i1 + 606;
                        for (register int i2 = i2_60927_lb; i2 <= i2_19851_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                      }
                    }
                    const int i1_52955_lb = max(96, -16 * ii2 + i0 - 399), i1_78745_ub = min(i0 - 384, -i0 + 594);
                    for (register int i1 = i1_52955_lb; i1 <= i1_78745_ub; i1 += 1) {
                      const int i2_20940_lb = max(i0 - i1 - 384, 32 * ii2 - i0 + i1 + 384), i2_87679_ub = 32 * ii2 - i0 + i1 + 415;
                      for (register int i2 = i2_20940_lb; i2 <= i2_87679_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                  } else if (ii1 == 1 && ii2 <= 48) {
                    if (ii2 == 0 && i0 == 496) {
                      for (register int i1 = 61; i1 <= 62; i1 += 1) {
                        const int i2_19444_lb = i1 - 47, i2_27622_ub = -i1 + 78;
                        for (register int i2 = i2_19444_lb; i2 <= i2_27622_ub; i2 += 1) {
                          S1(496, i1, i2);
                        }
                      }
                    } else {
                      const int i1_67475_lb = max(i0 - 435, -i0 + 543), i1_28813_ub = 63;
                      for (register int i1 = i1_67475_lb; i1 <= i1_28813_ub; i1 += 1) {
                        if (ii2 == 0) {
                          const int i2_16969_lb = i0 + i1 - 543, i2_39668_ub = i0 - i1 - 416;
                          for (register int i2 = i2_16969_lb; i2 < i2_39668_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                        } else if (ii2 == 1) {
                          const int i2_93641_lb = max(i0 + i1 - 543, -i0 - i1 + 575), i2_49240_ub = i0 - i1 - 416;
                          for (register int i2 = i2_93641_lb; i2 < i2_49240_ub; i2 += 1) {
                            S1(i0, i1, i2);
                          }
                        }
                        const int i2_45912_lb = max(i0 - i1 - 416, 32 * ii2 - i0 - i1 + 543), i2_55715_ub = 32 * ii2 - i0 - i1 + 574;
                        for (register int i2 = i2_45912_lb; i2 <= i2_55715_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                      }
                    }
                    const int i1_21019_lb = max(64, -16 * ii2 + i0 - 431), i1_66001_ub = min(i0 - 416, -i0 + 562);
                    for (register int i1 = i1_21019_lb; i1 <= i1_66001_ub; i1 += 1) {
                      const int i2_28015_lb = max(i0 - i1 - 416, 32 * ii2 - i0 + i1 + 416), i2_74366_ub = 32 * ii2 - i0 + i1 + 447;
                      for (register int i2 = i2_28015_lb; i2 <= i2_74366_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                  } else if (ii1 == 2) {
                    const int i1_43686_lb = max(i0 - 403, -i0 + 575), i1_78492_ub = min(i0 - 384, -i0 + 594);
                    for (register int i1 = i1_43686_lb; i1 <= i1_78492_ub; i1 += 1) {
                      const int i2_90869_lb = 0, i2_74009_ub = min(i0 + i1 - 575, i0 - i1 - 384);
                      for (register int i2 = i2_90869_lb; i2 < i2_74009_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                      const int i2_22129_lb = max(-i0 + i1 + 1952, -i0 - i1 + 2143), i2_15438_ub = 1599;
                      for (register int i2 = i2_22129_lb; i2 <= i2_15438_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                  } else {
                    const int i1_23885_lb = max(i0 - 435, -i0 + 543), i1_88778_ub = min(i0 - 416, -i0 + 562);
                    for (register int i1 = i1_23885_lb; i1 <= i1_88778_ub; i1 += 1) {
                      const int i2_53725_lb = 0, i2_1164_ub = min(i0 + i1 - 543, i0 - i1 - 416);
                      for (register int i2 = i2_53725_lb; i2 < i2_1164_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                      const int i2_8629_lb = max(-i0 + i1 + 1984, -i0 - i1 + 2111), i2_6681_ub = 1599;
                      for (register int i2 = i2_8629_lb; i2 <= i2_6681_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                  }
                }
              }
            }
          } else if (ii1 + 1 == k) {
            const int ii2_79909_lb = max(2 * k - 4, -2 * k + 2), ii2_45922_ub = 48;
            for (register int ii2 = ii2_79909_lb; ii2 <= ii2_45922_ub; ii2 += 1) {
              const int i0_10712_lb = 10 * k + 470, i0_54982_ub = min(min(499, 10 * k + 488), 20 * k + 8 * ii2 + 476);
              for (register int i0 = i0_10712_lb; i0 <= i0_54982_ub; i0 += 1) {
                if (k == 1 && ii2 == 0 && i0 == 496) {
                  for (register int i1 = 29; i1 <= 30; i1 += 1) {
                    const int i2_98957_lb = i1 - 15, i2_49142_ub = -i1 + 46;
                    for (register int i2 = i2_98957_lb; i2 <= i2_49142_ub; i2 += 1) {
                      S1(496, i1, i2);
                    }
                  }
                } else if (k == 1 && ii2 == 0 && i0 <= 495) {
                  const int i1_75322_lb = max(i0 - 467, -i0 + 511), i1_15926_ub = i0 - 463;
                  for (register int i1 = i1_75322_lb; i1 < i1_15926_ub; i1 += 1) {
                    const int i2_88811_lb = i0 + i1 - 511, i2_68964_ub = -i0 - i1 + 542;
                    for (register int i2 = i2_88811_lb; i2 <= i2_68964_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                }
                const int i1_65166_lb = max(max(i0 - 467, -19 * k - 16 * ii2 + i0 - 444), 20 * k - i0 + 491), i1_51075_ub = min(min(41, i0 - 448), 20 * k - i0 + 510);
                for (register int i1 = i1_65166_lb; i1 <= i1_51075_ub; i1 += 1) {
                  if (ii1 == 1) {
                    const int i2_41031_lb = max(i0 - i1 - 448, 32 * ii2 - i0 - i1 + 531), i2_2538_ub = 32 * ii2;
                    for (register int i2 = i2_41031_lb; i2 < i2_2538_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  } else if (ii2 == 0) {
                    const int i2_17077_lb = i0 + i1 - 511, i2_69046_ub = i0 - i1 - 448;
                    for (register int i2 = i2_17077_lb; i2 < i2_69046_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  } else if (ii2 == 1) {
                    const int i2_76904_lb = max(i0 + i1 - 511, -i0 - i1 + 543), i2_60763_ub = i0 - i1 - 448;
                    for (register int i2 = i2_76904_lb; i2 < i2_60763_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                  const int i2_63891_lb = max(max(max(max(19 * k + 32 * ii2 - 38, i0 + i1 - 531), -19 * k + i0 - i1 - 429), 32 * ii2 - i0 + i1 + 448), 32 * ii2 - i0 - i1 + 511), i2_84126_ub = min(19 * k + 32 * ii2 - i0 + i1 + 460, 32 * ii2 - i0 - i1 + 562);
                  for (register int i2 = i2_63891_lb; i2 <= i2_84126_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  if (ii1 == 0) {
                    const int i2_34773_lb = 32 * ii2 - i0 + i1 + 480, i2_86020_ub = 32 * ii2 - i0 - i1 + 542;
                    for (register int i2 = i2_34773_lb; i2 <= i2_86020_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                }
                if (k == 2) {
                  const int i1_99564_lb = 42, i1_75010_ub = min(i0 - 448, -i0 + 543);
                  for (register int i1 = i1_99564_lb; i1 <= i1_75010_ub; i1 += 1) {
                    const int i2_74798_lb = max(i0 - i1 - 448, 32 * ii2 - i0 + i1 + 448), i2_53290_ub = 32 * ii2 - i0 + i1 + 479;
                    for (register int i2 = i2_74798_lb; i2 <= i2_53290_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                  if (ii2 >= 1 && i0 <= 495) {
                    const int i2_92526_lb = 32 * ii2, i2_99780_ub = 32 * ii2 + 31;
                    for (register int i2 = i2_92526_lb; i2 <= i2_99780_ub; i2 += 1) {
                      S1(i0, -i0 + 543, i2);
                    }
                  }
                  const int i1_76323_lb = -i0 + 544, i1_88787_ub = min(i0 - 437, 32 * ii2 - i0 + 542);
                  for (register int i1 = i1_76323_lb; i1 <= i1_88787_ub; i1 += 1) {
                    const int i2_45702_lb = 32 * ii2 - i0 + i1 + 448, i2_87035_ub = 32 * ii2 - i0 - i1 + 542;
                    for (register int i2 = i2_45702_lb; i2 <= i2_87035_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                    const int i2_43769_lb = max(32 * ii2 - i0 + i1 + 436, 32 * ii2 - i0 - i1 + 543), i2_27369_ub = 32 * ii2;
                    for (register int i2 = i2_43769_lb; i2 < i2_27369_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                    const int i2_33544_lb = 32 * ii2, i2_42726_ub = min(32 * ii2 - i0 + i1 + 479, 32 * ii2 - i0 - i1 + 574);
                    for (register int i2 = i2_33544_lb; i2 <= i2_42726_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                    const int i2_92863_lb = 32 * ii2 - i0 - i1 + 575, i2_8867_ub = 32 * ii2 - i0 + i1 + 467;
                    for (register int i2 = i2_92863_lb; i2 <= i2_8867_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                  if (ii2 >= 1) {
                    const int i2_75004_lb = 32 * ii2, i2_81674_ub = 32 * ii2 + 31;
                    for (register int i2 = i2_75004_lb; i2 <= i2_81674_ub; i2 += 1) {
                      S1(i0, i0 - 436, i2);
                    }
                  } else if (ii2 == 0) {
                    const int i1_94183_lb = -i0 + 544, i1_40171_ub = i0 - 447;
                    for (register int i1 = i1_94183_lb; i1 < i1_40171_ub; i1 += 1) {
                      const int i2_32750_lb = max(i0 + i1 - 543, i0 - i1 - 448), i2_51566_ub = min(-i0 + i1 + 479, -i0 - i1 + 574);
                      for (register int i2 = i2_32750_lb; i2 <= i2_51566_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                    const int i1_42709_lb = max(i0 - 447, -i0 + 543), i1_66179_ub = i0 - 435;
                    for (register int i1 = i1_42709_lb; i1 < i1_66179_ub; i1 += 1) {
                      const int i2_36964_lb = i0 - i1 - 436, i2_35965_ub = i0 + i1 - 543;
                      for (register int i2 = i2_36964_lb; i2 < i2_35965_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                      const int i2_43294_lb = i0 + i1 - 543, i2_855_ub = i0 + i1 - 531;
                      for (register int i2 = i2_43294_lb; i2 < i2_855_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                      const int i2_20091_lb = i0 + i1 - 531, i2_94419_ub = -i0 - i1 + 574;
                      for (register int i2 = i2_20091_lb; i2 <= i2_94419_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                      const int i2_86876_lb = max(i0 + i1 - 531, -i0 - i1 + 575), i2_36008_ub = -i0 + i1 + 467;
                      for (register int i2 = i2_86876_lb; i2 <= i2_36008_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                  }
                }
              }
              if (k == 1) {
                const int i0_85781_lb = 480, i0_78026_ub = min(498, 8 * ii2 + 496);
                for (register int i0 = i0_85781_lb; i0 <= i0_78026_ub; i0 += 1) {
                  const int i1_5650_lb = max(0, -16 * ii2 + i0 - 495), i1_78307_ub = min(i0 - 480, -i0 + 498);
                  for (register int i1 = i1_5650_lb; i1 <= i1_78307_ub; i1 += 1) {
                    const int i2_94158_lb = max(i0 - i1 - 480, 32 * ii2 - i0 + i1 + 480), i2_81973_ub = 32 * ii2 - i0 + i1 + 511;
                    for (register int i2 = i2_94158_lb; i2 <= i2_81973_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                  const int i1_83446_lb = max(i0 + 1101, -i0 + 2079), i1_39860_ub = min(1599, 16 * ii2 - i0 + 2094);
                  for (register int i1 = i1_83446_lb; i1 <= i1_39860_ub; i1 += 1) {
                    const int i2_85360_lb = max(i0 + i1 - 2079, 32 * ii2 - i0 - i1 + 2079), i2_43568_ub = 32 * ii2 - i0 - i1 + 2110;
                    for (register int i2 = i2_85360_lb; i2 <= i2_43568_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                }
              }
            }
            if (k == 2) {
              for (register int i0 = 490; i0 <= 499; i0 += 1) {
                const int i1_86294_lb = -i0 + 531, i1_76445_ub = min(i0 - 448, -i0 + 543);
                for (register int i1 = i1_86294_lb; i1 <= i1_76445_ub; i1 += 1) {
                  const int i2_60476_lb = 0, i2_61299_ub = min(i0 + i1 - 531, i0 - i1 - 448);
                  for (register int i2 = i2_60476_lb; i2 < i2_61299_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  const int i2_58119_lb = max(-i0 + i1 + 2016, -i0 - i1 + 2099), i2_71011_ub = 1599;
                  for (register int i2 = i2_58119_lb; i2 <= i2_71011_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
                if (i0 <= 495) {
                  for (register int i2 = 1568; i2 <= 1599; i2 += 1) {
                    S1(i0, -i0 + 543, i2);
                  }
                }
                const int i1_38929_lb = -i0 + 544, i1_60531_ub = i0 - 436;
                for (register int i1 = i1_38929_lb; i1 < i1_60531_ub; i1 += 1) {
                  const int i2_89752_lb = 0, i2_75893_ub = i0 - i1 - 448;
                  for (register int i2 = i2_89752_lb; i2 < i2_75893_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  const int i2_12848_lb = -i0 + i1 + 2016, i2_49399_ub = -i0 - i1 + 2110;
                  for (register int i2 = i2_12848_lb; i2 <= i2_49399_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  const int i2_76749_lb = max(0, i0 - i1 - 448), i2_49292_ub = min(i0 + i1 - 543, i0 - i1 - 436);
                  for (register int i2 = i2_76749_lb; i2 < i2_49292_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  const int i2_43818_lb = max(-i0 + i1 + 2004, -i0 - i1 + 2111), i2_79977_ub = 1567;
                  for (register int i2 = i2_43818_lb; i2 <= i2_79977_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  if (i0 >= i1 + 448) {
                    for (register int i2 = 1568; i2 <= 1599; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  } else {
                    for (register int i2 = 1568; i2 <= 1599; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                }
                for (register int i2 = 1568; i2 <= 1599; i2 += 1) {
                  S1(i0, i0 - 436, i2);
                }
              }
            } else {
              for (register int i0 = 480; i0 <= 498; i0 += 1) {
                const int i1_8374_lb = max(i0 - 467, -i0 + 511), i1_90987_ub = min(i0 - 448, -i0 + 530);
                for (register int i1 = i1_8374_lb; i1 <= i1_90987_ub; i1 += 1) {
                  const int i2_67626_lb = 0, i2_8308_ub = min(i0 + i1 - 511, i0 - i1 - 448);
                  for (register int i2 = i2_67626_lb; i2 < i2_8308_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  const int i2_42596_lb = max(-i0 + i1 + 2016, -i0 - i1 + 2079), i2_53920_ub = 1599;
                  for (register int i2 = i2_42596_lb; i2 <= i2_53920_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
              }
            }
            if (k == 1) {
              for (register int i0 = 480; i0 <= 498; i0 += 1) {
                const int i1_31571_lb = 0, i1_59224_ub = min(i0 - 480, -i0 + 498);
                for (register int i1 = i1_31571_lb; i1 <= i1_59224_ub; i1 += 1) {
                  const int i2_90435_lb = 0, i2_49393_ub = i0 - i1 - 480;
                  for (register int i2 = i2_90435_lb; i2 < i2_49393_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  const int i2_66446_lb = -i0 + i1 + 2048, i2_29364_ub = 1599;
                  for (register int i2 = i2_66446_lb; i2 <= i2_29364_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
                const int i1_26276_lb = max(i0 + 1101, -i0 + 2079), i1_72550_ub = 1599;
                for (register int i1 = i1_26276_lb; i1 <= i1_72550_ub; i1 += 1) {
                  const int i2_5258_lb = 0, i2_39125_ub = i0 + i1 - 2079;
                  for (register int i2 = i2_5258_lb; i2 < i2_39125_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  const int i2_21949_lb = -i0 - i1 + 3647, i2_98359_ub = 1599;
                  for (register int i2 = i2_21949_lb; i2 <= i2_98359_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
              }
            }
          } else {
            for (register int ii2 = 0; ii2 <= 48; ii2 += 1) {
              for (register int i0 = 490; i0 <= 499; i0 += 1) {
                const int i1_11720_lb = -i0 + 499, i1_52691_ub = i0 - 479;
                for (register int i1 = i1_11720_lb; i1 < i1_52691_ub; i1 += 1) {
                  if (ii2 == 0) {
                    const int i2_97371_lb = i0 + i1 - 499, i2_52331_ub = i0 - i1 - 480;
                    for (register int i2 = i2_97371_lb; i2 < i2_52331_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                  const int i2_37557_lb = max(max(i0 - i1 - 480, 32 * ii2 - i0 + i1 + 480), 32 * ii2 - i0 - i1 + 499), i2_86646_ub = min(32 * ii2 - 1, 32 * ii2 - i0 - i1 + 510);
                  for (register int i2 = i2_37557_lb; i2 <= i2_86646_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  const int i2_76389_lb = max(i0 + i1 - 511, 32 * ii2 - i0 - i1 + 511), i2_45932_ub = 32 * ii2;
                  for (register int i2 = i2_76389_lb; i2 < i2_45932_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  const int i2_93985_lb = max(max(32 * ii2, i0 + i1 - 511), i0 - i1 - 480), i2_44015_ub = 32 * ii2 - i0 - i1 + 530;
                  for (register int i2 = i2_93985_lb; i2 <= i2_44015_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  const int i2_70592_lb = 32 * ii2 - i0 - i1 + 531, i2_36582_ub = min(32 * ii2 - i0 + i1 + 511, 32 * ii2 - i0 - i1 + 542);
                  for (register int i2 = i2_70592_lb; i2 <= i2_36582_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
                const int i1_14288_lb = max(i0 - 479, -i0 + 511), i1_55345_ub = i0 - 467;
                for (register int i1 = i1_14288_lb; i1 < i1_55345_ub; i1 += 1) {
                  if (ii2 == 0) {
                    const int i2_56006_lb = i0 + i1 - 511, i2_62211_ub = i0 - i1 - 468;
                    for (register int i2 = i2_56006_lb; i2 < i2_62211_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                  const int i2_14569_lb = max(max(i0 - i1 - 468, 32 * ii2 - i0 + i1 + 468), 32 * ii2 - i0 - i1 + 511), i2_46442_ub = 32 * ii2 - i0 + i1 + 499;
                  for (register int i2 = i2_14569_lb; i2 <= i2_46442_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  const int i2_27957_lb = 32 * ii2 - i0 + i1 + 500, i2_97367_ub = 32 * ii2 - i0 - i1 + 542;
                  for (register int i2 = i2_27957_lb; i2 <= i2_97367_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
              }
            }
            for (register int i0 = 490; i0 <= 499; i0 += 1) {
              const int i1_69918_lb = -i0 + 499, i1_97416_ub = i0 - 479;
              for (register int i1 = i1_69918_lb; i1 < i1_97416_ub; i1 += 1) {
                const int i2_9710_lb = 0, i2_8219_ub = min(i0 + i1 - 499, i0 - i1 - 480);
                for (register int i2 = i2_9710_lb; i2 < i2_8219_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
                const int i2_95775_lb = i0 - i1 - 480, i2_14479_ub = i0 + i1 - 511;
                for (register int i2 = i2_95775_lb; i2 < i2_14479_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
                const int i2_90339_lb = -i0 - i1 + 2079, i2_90463_ub = -i0 + i1 + 2047;
                for (register int i2 = i2_90339_lb; i2 <= i2_90463_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
                const int i2_20900_lb = max(-i0 + i1 + 2048, -i0 - i1 + 2067), i2_18411_ub = 1599;
                for (register int i2 = i2_20900_lb; i2 <= i2_18411_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
              const int i1_59507_lb = max(i0 - 479, -i0 + 511), i1_18271_ub = i0 - 467;
              for (register int i1 = i1_59507_lb; i1 < i1_18271_ub; i1 += 1) {
                const int i2_70743_lb = 0, i2_97064_ub = min(i0 + i1 - 511, i0 - i1 - 468);
                for (register int i2 = i2_70743_lb; i2 < i2_97064_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
                const int i2_21269_lb = max(-i0 + i1 + 2036, -i0 - i1 + 2079), i2_47132_ub = 1599;
                for (register int i2 = i2_21269_lb; i2 <= i2_47132_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
            }
          }
        }
      } else if (k == 2) {
        #pragma omp parallel for
        for (register int ii1 = 0; ii1 <= 1; ii1 += 1) {
          for (register int ii2 = 0; ii2 <= 49; ii2 += 1) {
            if (ii2 <= 48) {
              const int i0_68189_lb = 32 * ii0 + 16, i0_21788_ub = min(32 * ii0 + 15 * ii1 + 31, 32 * ii0 - 16 * ii1 + 16 * ii2 + 47);
              for (register int i0 = i0_68189_lb; i0 <= i0_21788_ub; i0 += 1) {
                if (2 * i0 >= 64 * ii0 + 33 * ii1 + 31) {
                  if (ii1 == 0) {
                    const int i1_85285_lb = 32 * ii0 - i0 + 31, i1_40547_ub = 15;
                    for (register int i1 = i1_85285_lb; i1 <= i1_40547_ub; i1 += 1) {
                      if (ii2 == 0) {
                        const int i2_83999_lb = -32 * ii0 + i0 + i1 - 31, i2_16207_ub = -32 * ii0 + i0 - i1;
                        for (register int i2 = i2_83999_lb; i2 < i2_16207_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                      }
                      const int i2_86989_lb = max(-32 * ii0 + i0 - i1, 32 * ii0 + 32 * ii2 - i0 - i1 + 31), i2_11956_ub = 32 * ii0 + 32 * ii2 - i0 + i1 + 31;
                      for (register int i2 = i2_86989_lb; i2 <= i2_11956_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                      const int i2_13574_lb = max(-32 * ii0 + i0 - i1, 32 * ii0 + 32 * ii2 - i0 + i1 + 32), i2_79148_ub = 32 * ii0 + 32 * ii2 - i0 - i1 + 62;
                      for (register int i2 = i2_13574_lb; i2 <= i2_79148_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                  }
                  const int i1_82542_lb = max(16, -32 * ii0 + 32 * ii1 + i0 - 31), i1_99844_ub = min(47, -32 * ii0 + 15 * ii1 + i0);
                  for (register int i1 = i1_82542_lb; i1 <= i1_99844_ub; i1 += 1) {
                    {
                      if (i0 >= 32 * ii0 + 32) {
                        const int i2_76564_lb = max(-32 * ii0 + i0 + i1 - 63, 32 * ii0 + 32 * ii2 - i0 - i1 + 63), i2_8604_ub = 32 * ii0 + 32 * ii2 - i0 + i1;
                        for (register int i2 = i2_76564_lb; i2 < i2_8604_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                      }
                      if (48 * ii1 + 498 >= i0 + i1) {
                        const int i2_8064_lb = max(-32 * ii0 + 32 * ii1 + i0 - i1, 32 * ii0 + 32 * ii2 - i0 + i1), i2_88692_ub = min(32 * ii0 + 32 * ii2 - i0 + i1 + 31, 32 * ii0 + 32 * ii2 - i0 - i1 + 94);
                        for (register int i2 = i2_8064_lb; i2 <= i2_88692_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                      }
                    }
                    if (ii0 == 14 && ii1 == 0 && i0 + i1 >= 499) {
                      const int i2_39436_lb = max(i0 - i1 - 448, 32 * ii2 - i0 + i1 + 448), i2_98403_ub = min(32 * ii2 - i0 + i1 + 479, 32 * ii2 - i0 - i1 + 524);
                      for (register int i2 = i2_39436_lb; i2 <= i2_98403_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                      const int i2_79155_lb = 32 * ii2 - i0 - i1 + 525, i2_60336_ub = 32 * ii2 - i0 + i1 + 479;
                      for (register int i2 = i2_79155_lb; i2 <= i2_60336_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                  }
                  if (ii1 == 1) {
                    const int i1_16815_lb = 48, i1_38662_ub = 32 * ii0 - i0 + 94;
                    for (register int i1 = i1_16815_lb; i1 <= i1_38662_ub; i1 += 1) {
                      const int i2_78608_lb = max(32 * ii0 + 32 * ii2 - i0 + i1 - 32, -32 * ii0 + i0 - i1 + 32), i2_3910_ub = 32 * ii0 + 32 * ii2 - i0 + i1;
                      for (register int i2 = i2_78608_lb; i2 < i2_3910_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                  }
                } else {
                  if (ii0 == 14) {
                    const int i1_52079_lb = -i0 + 511, i1_16229_ub = i0 - 435;
                    for (register int i1 = i1_52079_lb; i1 < i1_16229_ub; i1 += 1) {
                      const int i2_67394_lb = max(i0 + i1 - 511, 32 * ii2 - i0 - i1 + 511), i2_11427_ub = 32 * ii2 - i0 - i1 + 542;
                      for (register int i2 = i2_67394_lb; i2 <= i2_11427_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                  }
                  const int i1_47836_lb = max(i0 - 435, 32 * ii0 - i0 + 63), i1_74894_ub = 47;
                  for (register int i1 = i1_47836_lb; i1 <= i1_74894_ub; i1 += 1) {
                    if (ii2 == 0) {
                      const int i2_41368_lb = -32 * ii0 + i0 + i1 - 63, i2_16025_ub = -32 * ii0 + i0 - i1 + 31;
                      for (register int i2 = i2_41368_lb; i2 <= i2_16025_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                    const int i2_13034_lb = max(-32 * ii0 + i0 - i1 + 32, 32 * ii0 + 32 * ii2 - i0 - i1 + 63), i2_43005_ub = 32 * ii0 + 32 * ii2 - i0 + i1;
                    for (register int i2 = i2_13034_lb; i2 < i2_43005_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                    const int i2_56573_lb = max(32 * ii0 + 32 * ii2 - i0 + i1, -32 * ii0 + i0 - i1 + 32), i2_13386_ub = 32 * ii0 + 32 * ii2 - i0 - i1 + 94;
                    for (register int i2 = i2_56573_lb; i2 <= i2_13386_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                  const int i1_59212_lb = 48, i1_59914_ub = -32 * ii0 + i0 + 32;
                  for (register int i1 = i1_59212_lb; i1 <= i1_59914_ub; i1 += 1) {
                    const int i2_41694_lb = max(32 * ii0 + 32 * ii2 - i0 + i1 - 32, -32 * ii0 + i0 - i1 + 32), i2_72787_ub = 32 * ii0 + 32 * ii2 - i0 + i1;
                    for (register int i2 = i2_41694_lb; i2 < i2_72787_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                }
              }
              if (ii1 == 0 && ii2 == 0) {
                const int i0_55414_lb = 32 * ii0 + 32, i0_40588_ub = 32 * ii0 + 38;
                for (register int i0 = i0_55414_lb; i0 <= i0_40588_ub; i0 += 1) {
                  if (ii0 == 14) {
                    const int i1_72631_lb = i0 - 479, i1_48331_ub = -i0 + 492;
                    for (register int i1 = i1_72631_lb; i1 <= i1_48331_ub; i1 += 1) {
                      const int i2_49193_lb = i0 + i1 - 479, i2_97047_ub = -i0 - i1 + 510;
                      for (register int i2 = i2_49193_lb; i2 <= i2_97047_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                  }
                  const int i1_37023_lb = max(64 * ii0 - i0 - 403, -32 * ii0 + i0 - 31), i1_88629_ub = 32 * ii0 - i0 + 46;
                  for (register int i1 = i1_37023_lb; i1 <= i1_88629_ub; i1 += 1) {
                    const int i2_11803_lb = -32 * ii0 + i0 + i1 - 31, i2_32530_ub = 32 * ii0 - i0 - i1 + 62;
                    for (register int i2 = i2_11803_lb; i2 <= i2_32530_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                  const int i1_65317_lb = -32 * ii0 + i0 - 15, i1_28618_ub = 32 * ii0 - i0 + 62;
                  for (register int i1 = i1_65317_lb; i1 <= i1_28618_ub; i1 += 1) {
                    const int i2_87545_lb = -32 * ii0 + i0 - i1, i2_60277_ub = 32 * ii0 - i0 + i1 + 31;
                    for (register int i2 = i2_87545_lb; i2 <= i2_60277_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                }
              } else if (ii1 == 0 && ii2 >= 1) {
                const int i0_48880_lb = 32 * ii0 + 32, i0_39624_ub = 32 * ii0 + 46;
                for (register int i0 = i0_48880_lb; i0 <= i0_39624_ub; i0 += 1) {
                  if (ii0 == 14) {
                    const int i1_76507_lb = i0 - 479, i1_16274_ub = -i0 + 492;
                    for (register int i1 = i1_76507_lb; i1 <= i1_16274_ub; i1 += 1) {
                      const int i2_51051_lb = 32 * ii2 - i0 - i1 + 479, i2_24343_ub = 32 * ii2;
                      for (register int i2 = i2_51051_lb; i2 < i2_24343_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                      const int i2_7521_lb = 32 * ii2, i2_8771_ub = 32 * ii2 - i0 + i1 + 479;
                      for (register int i2 = i2_7521_lb; i2 <= i2_8771_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                      const int i2_40369_lb = 32 * ii2 - i0 + i1 + 480, i2_20555_ub = 32 * ii2 - i0 - i1 + 510;
                      for (register int i2 = i2_40369_lb; i2 <= i2_20555_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                  }
                  const int i1_51777_lb = max(64 * ii0 - i0 - 403, -32 * ii0 + i0 - 31), i1_13294_ub = 15;
                  for (register int i1 = i1_51777_lb; i1 <= i1_13294_ub; i1 += 1) {
                    const int i2_50293_lb = max(-32 * ii0 + i0 + i1 - 31, 32 * ii0 + 32 * ii2 - i0 - i1 + 31), i2_10989_ub = 32 * ii2;
                    for (register int i2 = i2_50293_lb; i2 < i2_10989_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                    const int i2_89560_lb = 32 * ii2, i2_91988_ub = 32 * ii0 + 32 * ii2 - i0 - i1 + 62;
                    for (register int i2 = i2_89560_lb; i2 <= i2_91988_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                  if (ii0 == 14) {
                    const int i1_83776_lb = 16, i1_44975_ub = -i0 + 510;
                    for (register int i1 = i1_83776_lb; i1 <= i1_44975_ub; i1 += 1) {
                      const int i2_48928_lb = max(i0 - i1 - 448, 32 * ii2 - i0 + i1 + 448), i2_72760_ub = 32 * ii2 - i0 + i1 + 479;
                      for (register int i2 = i2_48928_lb; i2 <= i2_72760_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                  } else {
                    const int i1_93306_lb = 16, i1_98121_ub = 32 * ii0 - i0 + 62;
                    for (register int i1 = i1_93306_lb; i1 <= i1_98121_ub; i1 += 1) {
                      const int i2_69807_lb = max(-32 * ii0 + i0 - i1, 32 * ii0 + 32 * ii2 - i0 + i1), i2_30329_ub = 32 * ii0 + 32 * ii2 - i0 + i1 + 31;
                      for (register int i2 = i2_69807_lb; i2 <= i2_30329_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                  }
                }
              } else if (ii2 == 0) {
                const int i0_3102_lb = 32 * ii0 + 32, i0_81610_ub = 32 * ii0 + 38;
                for (register int i0 = i0_3102_lb; i0 <= i0_81610_ub; i0 += 1) {
                  const int i1_79211_lb = -32 * ii0 + i0 + 1, i1_84772_ub = 32 * ii0 - i0 + 78;
                  for (register int i1 = i1_79211_lb; i1 <= i1_84772_ub; i1 += 1) {
                    const int i2_26580_lb = -32 * ii0 + i0 + i1 - 63, i2_66756_ub = 32 * ii0 - i0 - i1 + 94;
                    for (register int i2 = i2_26580_lb; i2 <= i2_66756_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                  const int i1_45049_lb = -32 * ii0 + i0 + 17, i1_91812_ub = 32 * ii0 - i0 + 94;
                  for (register int i1 = i1_45049_lb; i1 <= i1_91812_ub; i1 += 1) {
                    const int i2_6380_lb = -32 * ii0 + i0 - i1 + 32, i2_21556_ub = 32 * ii0 - i0 + i1;
                    for (register int i2 = i2_6380_lb; i2 < i2_21556_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                }
              }
            } else {
              const int i0_8087_lb = 32 * ii0 + 16, i0_73784_ub = 32 * ii0 + 46;
              for (register int i0 = i0_8087_lb; i0 <= i0_73784_ub; i0 += 1) {
                {
                  if (ii0 == 14 && 16 * ii1 + i0 >= 480 && 11 * ii1 + i0 <= 490) {
                    const int i1_45900_lb = max(i0 - 479, 30 * ii1 - i0 + 481), i1_15608_ub = min(i0 - 436, 30 * ii1 - i0 + 492);
                    for (register int i1 = i1_45900_lb; i1 <= i1_15608_ub; i1 += 1) {
                      if (ii1 == 0) {
                        const int i2_82555_lb = 0, i2_2621_ub = i0 + i1 - 479;
                        for (register int i2 = i2_82555_lb; i2 < i2_2621_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                        const int i2_52515_lb = -i0 - i1 + 2047, i2_34332_ub = 1567;
                        for (register int i2 = i2_52515_lb; i2 <= i2_34332_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                      } else {
                        const int i2_32267_lb = 0, i2_19161_ub = i0 + i1 - 511;
                        for (register int i2 = i2_32267_lb; i2 < i2_19161_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                      }
                      const int i2_45322_lb = max(-15 * ii1 + 1568, 30 * ii1 - i0 - i1 + 2049), i2_21827_ub = 1599;
                      for (register int i2 = i2_45322_lb; i2 <= i2_21827_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                  }
                  if (ii1 == 1 && 32 * ii0 + 31 >= i0) {
                    const int i1_27501_lb = max(i0 - 435, 32 * ii0 - i0 + 63), i1_45450_ub = -32 * ii0 + i0 + 32;
                    for (register int i1 = i1_27501_lb; i1 <= i1_45450_ub; i1 += 1) {
                      const int i2_83154_lb = 0, i2_76429_ub = min(-32 * ii0 + i0 + i1 - 64, -32 * ii0 + i0 - i1 + 31);
                      for (register int i2 = i2_83154_lb; i2 <= i2_76429_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                      const int i2_18210_lb = max(32 * ii0 - i0 + i1 + 1536, 32 * ii0 - i0 - i1 + 1631), i2_76460_ub = 1599;
                      for (register int i2 = i2_18210_lb; i2 <= i2_76460_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                  } else if (ii0 == 14 && ii1 == 0 && i0 >= 480) {
                    const int i1_90903_lb = max(i0 - 479, -i0 + 493), i1_4370_ub = -i0 + 510;
                    for (register int i1 = i1_90903_lb; i1 <= i1_4370_ub; i1 += 1) {
                      const int i2_23141_lb = 0, i2_94005_ub = min(i0 + i1 - 479, i0 - i1 - 448);
                      for (register int i2 = i2_23141_lb; i2 < i2_94005_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                      const int i2_2332_lb = max(-i0 + i1 + 2016, -i0 - i1 + 2047), i2_2353_ub = 1599;
                      for (register int i2 = i2_2332_lb; i2 <= i2_2353_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                  }
                }
                if (ii0 <= 13 && 2 * i0 >= 64 * ii0 + 33 * ii1 + 31 && 64 * ii0 + 31 * ii1 + 62 >= 2 * i0) {
                  const int i1_95129_lb = max(-32 * ii0 + 32 * ii1 + i0 - 31, 32 * ii0 - i0 + 31), i1_45265_ub = min(-32 * ii0 + 30 * ii1 + i0, 32 * ii0 - i0 + 94);
                  for (register int i1 = i1_95129_lb; i1 <= i1_45265_ub; i1 += 1) {
                    if (ii1 == 1) {
                      const int i2_85461_lb = 0, i2_40179_ub = min(-32 * ii0 + i0 + i1 - 64, -32 * ii0 + i0 - i1 + 31);
                      for (register int i2 = i2_85461_lb; i2 <= i2_40179_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                      const int i2_37077_lb = max(32 * ii0 - i0 + i1 + 1536, 32 * ii0 - i0 - i1 + 1631), i2_91842_ub = 1599;
                      for (register int i2 = i2_37077_lb; i2 <= i2_91842_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    } else {
                      const int i2_78087_lb = 0, i2_45164_ub = min(-32 * ii0 + i0 + i1 - 31, -32 * ii0 + i0 - i1);
                      for (register int i2 = i2_78087_lb; i2 < i2_45164_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                      const int i2_65626_lb = max(32 * ii0 - i0 + i1 + 1568, 32 * ii0 - i0 - i1 + 1599), i2_23987_ub = 1599;
                      for (register int i2 = i2_65626_lb; i2 <= i2_23987_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                  }
                } else if (ii0 == 14 && i0 >= 16 * ii1 + 464) {
                  if (ii1 == 0 && i0 >= 470 && i0 <= 479) {
                    const int i1_77124_lb = -i0 + 479, i1_48181_ub = -i0 + 480;
                    for (register int i1 = i1_77124_lb; i1 <= i1_48181_ub; i1 += 1) {
                      if (i0 + i1 == 480) {
                        S1(i0, -i0 + 480, 0);
                      }
                      const int i2_42960_lb = -i0 - i1 + 2047, i2_45992_ub = 1599;
                      for (register int i2 = i2_42960_lb; i2 <= i2_45992_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                  }
                  if (ii1 == 0 && i0 <= 479) {
                    const int i1_98866_lb = -i0 + 481, i1_91579_ub = i0 - 459;
                    for (register int i1 = i1_98866_lb; i1 < i1_91579_ub; i1 += 1) {
                      const int i2_65153_lb = 0, i2_44188_ub = min(i0 + i1 - 479, i0 - i1 - 448);
                      for (register int i2 = i2_65153_lb; i2 < i2_44188_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                      const int i2_29759_lb = max(-i0 + i1 + 2016, -i0 - i1 + 2047), i2_9006_ub = 1599;
                      for (register int i2 = i2_29759_lb; i2 <= i2_9006_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                  }
                  if (ii1 == 0) {
                    const int i1_89638_lb = max(i0 - 459, -i0 + 479), i1_12913_ub = min(i0 - 448, -i0 + 498);
                    for (register int i1 = i1_89638_lb; i1 <= i1_12913_ub; i1 += 1) {
                      const int i2_85435_lb = 0, i2_24201_ub = min(i0 + i1 - 479, i0 - i1 - 448);
                      for (register int i2 = i2_85435_lb; i2 < i2_24201_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                      const int i2_5726_lb = max(-i0 + i1 + 2016, -i0 - i1 + 2047), i2_76338_ub = 1599;
                      for (register int i2 = i2_5726_lb; i2 <= i2_76338_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                  } else {
                    const int i1_28571_lb = i0 - 447, i1_45219_ub = min(i0 - 437, -i0 + 542);
                    for (register int i1 = i1_28571_lb; i1 <= i1_45219_ub; i1 += 1) {
                      const int i2_86696_lb = 0, i2_30903_ub = min(i0 + i1 - 511, i0 - i1 - 416);
                      for (register int i2 = i2_86696_lb; i2 < i2_30903_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                      const int i2_47572_lb = max(-i0 + i1 + 1984, -i0 - i1 + 2079), i2_81825_ub = 1599;
                      for (register int i2 = i2_47572_lb; i2 <= i2_81825_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                    const int i1_76168_lb = i0 - 436, i1_49386_ub = min(i0 - 430, -i0 + 542);
                    for (register int i1 = i1_76168_lb; i1 <= i1_49386_ub; i1 += 1) {
                      const int i2_38356_lb = 0, i2_29598_ub = min(i0 + i1 - 511, i0 - i1 - 416);
                      for (register int i2 = i2_38356_lb; i2 < i2_29598_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                      const int i2_41228_lb = max(-i0 + i1 + 1984, -i0 - i1 + 2079), i2_16444_ub = 1599;
                      for (register int i2 = i2_41228_lb; i2 <= i2_16444_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                  }
                  if (11 * ii1 + 479 >= i0) {
                    const int i1_74762_lb = max(30 * ii1 + i0 - 459, -i0 + 499), i1_23206_ub = min(30 * ii1 + i0 - 448, -i0 + 542);
                    for (register int i1 = i1_74762_lb; i1 <= i1_23206_ub; i1 += 1) {
                      if (ii1 == 1) {
                        const int i2_56783_lb = 0, i2_68239_ub = i0 - 480;
                        for (register int i2 = i2_56783_lb; i2 < i2_68239_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                        const int i2_87739_lb = i0 - 480, i2_99744_ub = i0 - i1 - 416;
                        for (register int i2 = i2_87739_lb; i2 < i2_99744_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                        const int i2_30583_lb = -i0 + i1 + 1984, i2_86605_ub = -i0 + 2047;
                        for (register int i2 = i2_30583_lb; i2 <= i2_86605_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                      } else {
                        const int i2_7675_lb = 0, i2_12088_ub = i0 - i1 - 448;
                        for (register int i2 = i2_7675_lb; i2 < i2_12088_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                        const int i2_30793_lb = -i0 + i1 + 2016, i2_37434_ub = -i0 - i1 + 2066;
                        for (register int i2 = i2_30793_lb; i2 <= i2_37434_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                      }
                      const int i2_21094_lb = max(max(16 * ii1 - i0 + 2032, -30 * ii1 - i0 + i1 + 2016), -i0 - i1 + 2067), i2_36784_ub = 1599;
                      for (register int i2 = i2_21094_lb; i2 <= i2_36784_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                  }
                } else if (ii1 == 0 && i0 >= 32 * ii0 + 32) {
                  const int i1_50348_lb = -32 * ii0 + i0 - 31, i1_22881_ub = 32 * ii0 - i0 + 62;
                  for (register int i1 = i1_50348_lb; i1 <= i1_22881_ub; i1 += 1) {
                    const int i2_60985_lb = 0, i2_72426_ub = min(-32 * ii0 + i0 + i1 - 31, -32 * ii0 + i0 - i1);
                    for (register int i2 = i2_60985_lb; i2 < i2_72426_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                    const int i2_99220_lb = max(32 * ii0 - i0 + i1 + 1568, 32 * ii0 - i0 - i1 + 1599), i2_5908_ub = 1599;
                    for (register int i2 = i2_99220_lb; i2 <= i2_5908_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                }
              }
            }
          }
        }
        #pragma omp parallel for
        for (register int ii1 = 2; ii1 <= 49; ii1 += 1) {
          if (ii1 <= 47) {
            for (register int ii2 = 0; ii2 <= 48; ii2 += 1) {
              const int i0_84093_lb = 32 * ii0 + 16, i0_29332_ub = min(min(32 * ii0 + 46, 32 * ii0 + 15 * ii1 + 1), 32 * ii0 + 8 * ii2 + 38);
              for (register int i0 = i0_84093_lb; i0 <= i0_29332_ub; i0 += 1) {
                if (32 * ii0 + 31 >= i0) {
                  const int i1_30956_lb = 32 * ii0 + 32 * ii1 - i0 + 31, i1_22450_ub = 32 * ii1 + 15;
                  for (register int i1 = i1_30956_lb; i1 <= i1_22450_ub; i1 += 1) {
                    const int i2_58930_lb = max(-32 * ii0 - 32 * ii1 + i0 + i1 - 31, 32 * ii0 + 32 * ii1 + 32 * ii2 - i0 - i1 + 31), i2_88536_ub = 32 * ii0 + 32 * ii1 + 32 * ii2 - i0 - i1 + 62;
                    for (register int i2 = i2_58930_lb; i2 <= i2_88536_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                } else {
                  const int i1_38894_lb = -32 * ii0 + 32 * ii1 + i0 - 31, i1_50044_ub = min(32 * ii1 + 15, 32 * ii0 + 32 * ii1 + 16 * ii2 - i0 + 46);
                  for (register int i1 = i1_38894_lb; i1 <= i1_50044_ub; i1 += 1) {
                    const int i2_11742_lb = max(-32 * ii0 - 32 * ii1 + i0 + i1 - 31, 32 * ii0 + 32 * ii1 + 32 * ii2 - i0 - i1 + 31), i2_95677_ub = 32 * ii0 + 32 * ii1 + 32 * ii2 - i0 - i1 + 62;
                    for (register int i2 = i2_11742_lb; i2 <= i2_95677_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                }
                const int i1_34635_lb = max(32 * ii1 + 16, -32 * ii0 + 32 * ii1 - 16 * ii2 + i0 - 15), i1_15833_ub = min(-32 * ii0 + 32 * ii1 + i0, 32 * ii0 + 32 * ii1 - i0 + 62);
                for (register int i1 = i1_34635_lb; i1 <= i1_15833_ub; i1 += 1) {
                  if (ii2 == 0 && 32 * ii0 + 31 >= i0) {
                    const int i2_11773_lb = -32 * ii0 + 32 * ii1 + i0 - i1, i2_65218_ub = -32 * ii0 - 32 * ii1 + i0 + i1 - 31;
                    for (register int i2 = i2_11773_lb; i2 < i2_65218_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  } else if (ii2 == 1) {
                    const int i2_2439_lb = max(-32 * ii0 + 32 * ii1 + i0 - i1, 32 * ii0 - 32 * ii1 - i0 + i1 + 32), i2_19449_ub = -32 * ii0 - 32 * ii1 + i0 + i1 - 31;
                    for (register int i2 = i2_2439_lb; i2 < i2_19449_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                  const int i2_93658_lb = max(-32 * ii0 - 32 * ii1 + i0 + i1 - 31, 32 * ii0 - 32 * ii1 + 32 * ii2 - i0 + i1), i2_49584_ub = 32 * ii0 + 32 * ii1 + 32 * ii2 - i0 - i1 + 62;
                  for (register int i2 = i2_93658_lb; i2 <= i2_49584_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  if (32 * ii0 + 31 >= i0) {
                    const int i2_56883_lb = max(-32 * ii0 - 32 * ii1 + i0 + i1 - 31, 32 * ii0 + 32 * ii1 + 32 * ii2 - i0 - i1 + 63), i2_31104_ub = 32 * ii0 - 32 * ii1 + 32 * ii2 - i0 + i1 + 31;
                    for (register int i2 = i2_56883_lb; i2 <= i2_31104_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  } else {
                    const int i2_2720_lb = max(-32 * ii0 + 32 * ii1 + i0 - i1, 32 * ii0 + 32 * ii1 + 32 * ii2 - i0 - i1 + 63), i2_23583_ub = 32 * ii0 - 32 * ii1 + 32 * ii2 - i0 + i1 + 31;
                    for (register int i2 = i2_2720_lb; i2 <= i2_23583_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                }
              }
              if (ii1 == 2) {
                const int i0_53986_lb = 32 * ii0 + 32, i0_63705_ub = min(32 * ii0 + 46, 32 * ii0 + 7 * ii2 + (ii2 + 3) / 4 + 38);
                for (register int i0 = i0_53986_lb; i0 <= i0_63705_ub; i0 += 1) {
                  if (ii2 >= 2) {
                    const int i1_96009_lb = -32 * ii0 + i0 + 33, i1_53206_ub = 79;
                    for (register int i1 = i1_96009_lb; i1 <= i1_53206_ub; i1 += 1) {
                      const int i2_85965_lb = 32 * ii0 + 32 * ii2 - i0 - i1 + 95, i2_13655_ub = 32 * ii0 + 32 * ii2 - i0 - i1 + 126;
                      for (register int i2 = i2_85965_lb; i2 <= i2_13655_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                    const int i1_55474_lb = 80, i1_22777_ub = 32 * ii0 - i0 + 126;
                    for (register int i1 = i1_55474_lb; i1 <= i1_22777_ub; i1 += 1) {
                      const int i2_95225_lb = 32 * ii0 + 32 * ii2 - i0 + i1 - 64, i2_39567_ub = 32 * ii0 + 32 * ii2 - i0 + i1 - 32;
                      for (register int i2 = i2_95225_lb; i2 < i2_39567_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                  } else {
                    const int i1_68461_lb = -32 * ii0 + i0 + 33, i1_42533_ub = min(79, 32 * ii0 + 15 * ii2 - i0 + 110);
                    for (register int i1 = i1_68461_lb; i1 <= i1_42533_ub; i1 += 1) {
                      const int i2_78369_lb = max(-32 * ii0 + i0 + i1 - 95, 32 * ii0 + 28 * ii2 - i0 - i1 + 99), i2_27391_ub = min(-32 * ii0 + i0 - i1 + 63, 32 * ii0 + 29 * ii2 - i0 - i1 + 126);
                      for (register int i2 = i2_78369_lb; i2 <= i2_27391_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                      if (ii2 == 1) {
                        const int i2_31069_lb = -32 * ii0 + i0 - i1 + 64, i2_17263_ub = 32 * ii0 - i0 - i1 + 158;
                        for (register int i2 = i2_31069_lb; i2 <= i2_17263_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                      }
                    }
                    const int i1_93787_lb = max(80, -32 * ii0 - 16 * ii2 + i0 + 49), i1_59163_ub = 32 * ii0 - i0 + 126;
                    for (register int i1 = i1_93787_lb; i1 <= i1_59163_ub; i1 += 1) {
                      const int i2_12941_lb = max(32 * ii0 + 32 * ii2 - i0 + i1 - 64, -32 * ii0 + i0 - i1 + 64), i2_28423_ub = 32 * ii0 + 32 * ii2 - i0 + i1 - 32;
                      for (register int i2 = i2_12941_lb; i2 < i2_28423_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                  }
                }
              }
            }
          } else {
            for (register int ii2 = 0; ii2 <= 48; ii2 += 1) {
              const int i0_9993_lb = 32 * ii0 + 16, i0_93787_ub = min(32 * ii0 + 46, 32 * ii0 + 8 * ii2 + 38);
              for (register int i0 = i0_9993_lb; i0 <= i0_93787_ub; i0 += 1) {
                if (ii1 == 49 && i0 >= 32 * ii0 + 32) {
                  if (ii0 == 14 && ii2 == 0) {
                    const int i1_60515_lb = i0 + 1089, i1_20004_ub = min(i0 + 1099, -i0 + 2062);
                    for (register int i1 = i1_60515_lb; i1 <= i1_20004_ub; i1 += 1) {
                      const int i2_59724_lb = i0 + i1 - 2047, i2_33751_ub = -i0 - i1 + 2078;
                      for (register int i2 = i2_59724_lb; i2 <= i2_33751_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                  } else if (ii0 == 14 && ii2 >= 1) {
                    const int i1_51108_lb = i0 + 1089, i1_62444_ub = min(i0 + 1099, -i0 + 2078);
                    for (register int i1 = i1_51108_lb; i1 <= i1_62444_ub; i1 += 1) {
                      if (ii2 == 1) {
                        const int i2_57334_lb = max(i0 + i1 - 2047, -i0 - i1 + 2079), i2_5094_ub = min(i0 + i1 - 2030, i0 - i1 + 1119);
                        for (register int i2 = i2_57334_lb; i2 <= i2_5094_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                      }
                      const int i2_42502_lb = max(max(32 * ii2 - i0 + i1 - 1120, i0 - i1 + 1120), 32 * ii2 - i0 - i1 + 2047), i2_69696_ub = min(32 * ii2, i0 + i1 - 2029);
                      for (register int i2 = i2_42502_lb; i2 < i2_69696_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                      const int i2_74652_lb = max(max(i0 + i1 - 2029, 32 * ii2 - i0 + i1 - 1120), 32 * ii2 - i0 - i1 + 2047), i2_28467_ub = 32 * ii2 - i0 + 479;
                      for (register int i2 = i2_74652_lb; i2 <= i2_28467_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                      const int i2_99703_lb = max(32 * ii2 - i0 + 480, i0 + i1 - 2029), i2_46478_ub = min(32 * ii2 - 1, 32 * ii2 - i0 - i1 + 2060);
                      for (register int i2 = i2_99703_lb; i2 <= i2_46478_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                      if (ii2 == 1) {
                        const int i2_67596_lb = 32, i2_94928_ub = min(i0 + i1 - 2029, -i0 + i1 - 1056);
                        for (register int i2 = i2_67596_lb; i2 < i2_94928_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                      }
                      const int i2_86046_lb = 32 * ii2, i2_36057_ub = 32 * ii2 - i0 - i1 + 2060;
                      for (register int i2 = i2_86046_lb; i2 <= i2_36057_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                      const int i2_37461_lb = max(max(32 * ii2 - i0 + 480, i0 + i1 - 2029), 32 * ii2 - i0 - i1 + 2061), i2_64415_ub = 32 * ii2 - i0 + i1 - 1088;
                      for (register int i2 = i2_37461_lb; i2 < i2_64415_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                      const int i2_79800_lb = max(32 * ii2 - i0 + i1 - 1088, 32 * ii2 - i0 - i1 + 2061), i2_84882_ub = 32 * ii2 - i0 - i1 + 2078;
                      for (register int i2 = i2_79800_lb; i2 <= i2_84882_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                  }
                  const int i1_98031_lb = max(i0 + 1100, -32 * ii0 + i0 + 1537), i1_89940_ub = min(1583, 32 * ii0 + 16 * ii2 - i0 + 1614);
                  for (register int i1 = i1_98031_lb; i1 <= i1_89940_ub; i1 += 1) {
                    const int i2_44045_lb = max(-32 * ii0 + i0 + i1 - 1599, 32 * ii0 + 32 * ii2 - i0 - i1 + 1599), i2_27324_ub = 32 * ii2;
                    for (register int i2 = i2_44045_lb; i2 < i2_27324_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                    if (ii0 == 14 && ii2 >= 1 && i0 == 480 && i1 == 1580) {
                      S1(480, 1580, 32 * ii2);
                    }
                    const int i2_34715_lb = max(max(32 * ii2, 64 * ii0 + 32 * ii2 - i0 - i1 + 1165), -32 * ii0 + i0 - i1 + 1568), i2_19041_ub = 32 * ii0 + 32 * ii2 - i0 + i1 - 1536;
                    for (register int i2 = i2_34715_lb; i2 < i2_19041_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                    const int i2_68390_lb = max(-32 * ii0 + i0 + i1 - 1599, 32 * ii0 + 32 * ii2 - i0 + i1 - 1536), i2_44708_ub = 32 * ii0 + 32 * ii2 - i0 - i1 + 1630;
                    for (register int i2 = i2_68390_lb; i2 <= i2_44708_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                  const int i1_29181_lb = max(max(1584, i0 + 1100), -32 * ii0 - 16 * ii2 + i0 + 1553), i1_45258_ub = 32 * ii0 - i0 + 1630;
                  for (register int i1 = i1_29181_lb; i1 <= i1_45258_ub; i1 += 1) {
                    const int i2_81064_lb = max(32 * ii0 + 32 * ii2 - i0 + i1 - 1568, -32 * ii0 + i0 - i1 + 1568), i2_88905_ub = 32 * ii2;
                    for (register int i2 = i2_81064_lb; i2 < i2_88905_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                    const int i2_79009_lb = max(32 * ii2, -32 * ii0 + i0 - i1 + 1568), i2_32173_ub = 32 * ii0 + 32 * ii2 - i0 + i1 - 1536;
                    for (register int i2 = i2_79009_lb; i2 < i2_32173_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                } else if (ii1 == 49 && 32 * ii0 + 31 >= i0) {
                  if (ii0 == 14) {
                    const int i1_51349_lb = -i0 + 2047, i1_52695_ub = i0 + 1100;
                    for (register int i1 = i1_51349_lb; i1 <= i1_52695_ub; i1 += 1) {
                      const int i2_53619_lb = max(i0 + i1 - 2047, 32 * ii2 - i0 - i1 + 2047), i2_10203_ub = 32 * ii2 - i0 - i1 + 2078;
                      for (register int i2 = i2_53619_lb; i2 <= i2_10203_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                  }
                  const int i1_22391_lb = max(i0 + 1101, 32 * ii0 - i0 + 1599), i1_44624_ub = 1583;
                  for (register int i1 = i1_22391_lb; i1 <= i1_44624_ub; i1 += 1) {
                    const int i2_38671_lb = max(-32 * ii0 + i0 + i1 - 1599, 32 * ii0 + 32 * ii2 - i0 - i1 + 1599), i2_38446_ub = 32 * ii0 + 32 * ii2 - i0 - i1 + 1630;
                    for (register int i2 = i2_38671_lb; i2 <= i2_38446_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                  const int i1_91102_lb = 1584, i1_22619_ub = -32 * ii0 + i0 + 1568;
                  for (register int i1 = i1_91102_lb; i1 <= i1_22619_ub; i1 += 1) {
                    const int i2_49726_lb = max(32 * ii0 + 32 * ii2 - i0 + i1 - 1568, -32 * ii0 + i0 - i1 + 1568), i2_93500_ub = 32 * ii0 + 32 * ii2 - i0 + i1 - 1536;
                    for (register int i2 = i2_49726_lb; i2 < i2_93500_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                } else if (i0 >= 32 * ii0 + 16 * ii2 + 16) {
                  const int i1_75029_lb = max(-32 * ii0 + i0 + 1505, 32 * ii0 - i0 + 1567), i1_87187_ub = min(1551, 32 * ii0 + 15 * ii2 - i0 + 1582);
                  for (register int i1 = i1_75029_lb; i1 <= i1_87187_ub; i1 += 1) {
                    const int i2_74268_lb = max(-32 * ii0 + i0 + i1 - 1567, 32 * ii0 + 32 * ii2 - i0 - i1 + 1567), i2_54829_ub = min(-32 * ii0 + i0 - i1 + 1535, 32 * ii0 + 29 * ii2 - i0 - i1 + 1598);
                    for (register int i2 = i2_74268_lb; i2 <= i2_54829_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                    const int i2_88421_lb = -32 * ii0 + i0 - i1 + 1536, i2_72299_ub = 32 * ii0 + 32 * ii2 - i0 - i1 + 1598;
                    for (register int i2 = i2_88421_lb; i2 <= i2_72299_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                  const int i1_61121_lb = max(1552, -32 * ii0 - 16 * ii2 + i0 + 1521), i1_32466_ub = min(-32 * ii0 + i0 + 1536, 32 * ii0 - i0 + 1598);
                  for (register int i1 = i1_61121_lb; i1 <= i1_32466_ub; i1 += 1) {
                    const int i2_15975_lb = max(32 * ii0 + 32 * ii2 - i0 + i1 - 1536, -32 * ii0 + i0 - i1 + 1536), i2_95836_ub = 32 * ii0 + 32 * ii2 - i0 + i1 - 1504;
                    for (register int i2 = i2_15975_lb; i2 < i2_95836_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                } else {
                  const int i1_67860_lb = max(-32 * ii0 + i0 + 1505, 32 * ii0 - i0 + 1567), i1_84365_ub = min(-32 * ii0 + i0 + 1536, 32 * ii0 - i0 + 1598);
                  for (register int i1 = i1_67860_lb; i1 <= i1_84365_ub; i1 += 1) {
                    const int i2_56897_lb = max(32 * ii0 + 32 * ii2 - i0 + i1 - 1536, 32 * ii0 + 32 * ii2 - i0 - i1 + 1567), i2_97041_ub = 32 * ii0 + 32 * ii2 - i0 + i1 - 1504;
                    for (register int i2 = i2_56897_lb; i2 < i2_97041_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                    const int i2_45975_lb = 32 * ii0 + 32 * ii2 - i0 + i1 - 1504, i2_37961_ub = 32 * ii0 + 32 * ii2 - i0 - i1 + 1598;
                    for (register int i2 = i2_45975_lb; i2 <= i2_37961_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                }
              }
            }
          }
          const int i0_85946_lb = 32 * ii0 + 16, i0_24984_ub = 32 * ii0 + 46;
          for (register int i0 = i0_85946_lb; i0 <= i0_24984_ub; i0 += 1) {
            if (ii1 <= 47 && 32 * ii0 + 31 >= i0) {
              const int i1_86486_lb = 32 * ii0 + 32 * ii1 - i0 + 31, i1_53647_ub = -32 * ii0 + 32 * ii1 + i0;
              for (register int i1 = i1_86486_lb; i1 <= i1_53647_ub; i1 += 1) {
                const int i2_94032_lb = 0, i2_56458_ub = min(-32 * ii0 - 32 * ii1 + i0 + i1 - 31, -32 * ii0 + 32 * ii1 + i0 - i1);
                for (register int i2 = i2_94032_lb; i2 < i2_56458_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
                const int i2_63851_lb = max(32 * ii0 - 32 * ii1 - i0 + i1 + 1568, 32 * ii0 + 32 * ii1 - i0 - i1 + 1599), i2_32775_ub = 1599;
                for (register int i2 = i2_63851_lb; i2 <= i2_32775_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
            } else if (ii1 >= 3 && ii1 <= 47 && i0 >= 32 * ii0 + 32) {
              const int i1_1082_lb = -32 * ii0 + 32 * ii1 + i0 - 31, i1_18874_ub = 32 * ii0 + 32 * ii1 - i0 + 62;
              for (register int i1 = i1_1082_lb; i1 <= i1_18874_ub; i1 += 1) {
                const int i2_71222_lb = 0, i2_8536_ub = min(-32 * ii0 - 32 * ii1 + i0 + i1 - 31, -32 * ii0 + 32 * ii1 + i0 - i1);
                for (register int i2 = i2_71222_lb; i2 < i2_8536_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
                const int i2_41493_lb = max(32 * ii0 - 32 * ii1 - i0 + i1 + 1568, 32 * ii0 + 32 * ii1 - i0 - i1 + 1599), i2_20948_ub = 1599;
                for (register int i2 = i2_41493_lb; i2 <= i2_20948_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
            } else if (ii1 == 48) {
              const int i1_18389_lb = max(-32 * ii0 + i0 + 1505, 32 * ii0 - i0 + 1567), i1_16522_ub = min(-32 * ii0 + i0 + 1536, 32 * ii0 - i0 + 1598);
              for (register int i1 = i1_18389_lb; i1 <= i1_16522_ub; i1 += 1) {
                const int i2_24488_lb = 0, i2_92657_ub = min(-32 * ii0 + i0 + i1 - 1568, -32 * ii0 + i0 - i1 + 1535);
                for (register int i2 = i2_24488_lb; i2 <= i2_92657_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
                const int i2_87704_lb = max(32 * ii0 - i0 + i1 + 32, 32 * ii0 - i0 - i1 + 3135), i2_12909_ub = 1599;
                for (register int i2 = i2_87704_lb; i2 <= i2_12909_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
            } else if (ii1 == 49 && i0 >= 32 * ii0 + 32) {
              const int i1_81308_lb = -32 * ii0 + i0 + 1537, i1_48825_ub = 32 * ii0 - i0 + 1630;
              for (register int i1 = i1_81308_lb; i1 <= i1_48825_ub; i1 += 1) {
                const int i2_45376_lb = 0, i2_97283_ub = min(-32 * ii0 + i0 + i1 - 1600, -32 * ii0 + i0 - i1 + 1567);
                for (register int i2 = i2_45376_lb; i2 <= i2_97283_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
                const int i2_61014_lb = max(32 * ii0 - i0 + i1, 32 * ii0 - i0 - i1 + 3167), i2_13236_ub = 1599;
                for (register int i2 = i2_61014_lb; i2 <= i2_13236_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
            } else if (ii1 == 49 && 32 * ii0 + 31 >= i0) {
              if (ii0 == 14) {
                const int i1_98000_lb = -i0 + 2047, i1_17911_ub = i0 + 1100;
                for (register int i1 = i1_98000_lb; i1 <= i1_17911_ub; i1 += 1) {
                  const int i2_26629_lb = 0, i2_43976_ub = i0 + i1 - 2047;
                  for (register int i2 = i2_26629_lb; i2 < i2_43976_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  const int i2_55872_lb = -i0 - i1 + 3615, i2_12575_ub = 1599;
                  for (register int i2 = i2_55872_lb; i2 <= i2_12575_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
              }
              const int i1_85312_lb = max(i0 + 1101, 32 * ii0 - i0 + 1599), i1_58711_ub = -32 * ii0 + i0 + 1568;
              for (register int i1 = i1_85312_lb; i1 <= i1_58711_ub; i1 += 1) {
                const int i2_66222_lb = 0, i2_79344_ub = min(-32 * ii0 + i0 + i1 - 1600, -32 * ii0 + i0 - i1 + 1567);
                for (register int i2 = i2_66222_lb; i2 <= i2_79344_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
                const int i2_15169_lb = max(32 * ii0 - i0 + i1, 32 * ii0 - i0 - i1 + 3167), i2_30073_ub = 1599;
                for (register int i2 = i2_15169_lb; i2 <= i2_30073_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
            } else {
              const int i1_28472_lb = -32 * ii0 + i0 + 33, i1_32603_ub = 32 * ii0 - i0 + 126;
              for (register int i1 = i1_28472_lb; i1 <= i1_32603_ub; i1 += 1) {
                const int i2_48947_lb = 0, i2_99694_ub = min(-32 * ii0 + i0 + i1 - 96, -32 * ii0 + i0 - i1 + 63);
                for (register int i2 = i2_48947_lb; i2 <= i2_99694_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
                const int i2_57491_lb = max(32 * ii0 - i0 + i1 + 1504, 32 * ii0 - i0 - i1 + 1663), i2_6793_ub = 1599;
                for (register int i2 = i2_57491_lb; i2 <= i2_6793_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
            }
          }
        }
      } else {
        #pragma omp parallel for
        for (register int ii1 = 0; ii1 <= 47; ii1 += 1) {
          if (ii1 >= 3) {
            for (register int ii2 = 0; ii2 <= 49; ii2 += 1) {
              const int i0_84889_lb = 32 * ii0, i0_27371_ub = min(32 * ii0 + 30, 32 * ii0 + 8 * ii2 + 22);
              for (register int i0 = i0_84889_lb; i0 <= i0_27371_ub; i0 += 1) {
                const int i1_90744_lb = max(-32 * ii0 + 32 * ii1 + i0 + 1, 32 * ii0 + 32 * ii1 - i0 + 31), i1_82549_ub = min(min(-32 * ii0 + 32 * ii1 + i0 + 32, 32 * ii0 + 32 * ii1 + 16 * ii2 - i0 + 46), 32 * ii0 + 32 * ii1 - i0 + 62);
                for (register int i1 = i1_90744_lb; i1 <= i1_82549_ub; i1 += 1) {
                  if (ii2 <= 48) {
                    const int i2_92549_lb = max(32 * ii0 - 32 * ii1 + 32 * ii2 - i0 + i1 - 32, -32 * ii0 + 32 * ii1 + i0 - i1 + 32), i2_52472_ub = -32 * ii0 - 32 * ii1 + i0 + i1 - 31;
                    for (register int i2 = i2_92549_lb; i2 < i2_52472_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                    const int i2_96184_lb = max(max(32 * ii0 - 32 * ii1 + 32 * ii2 - i0 + i1 - 32, -32 * ii0 - 32 * ii1 + i0 + i1 - 31), 32 * ii0 + 32 * ii1 + 32 * ii2 - i0 - i1 + 31), i2_53563_ub = 32 * ii2;
                    for (register int i2 = i2_96184_lb; i2 < i2_53563_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                    const int i2_82060_lb = max(32 * ii2, -32 * ii0 - 32 * ii1 + i0 + i1 - 31), i2_94185_ub = 32 * ii0 + 32 * ii1 + 32 * ii2 - i0 - i1 + 62;
                    for (register int i2 = i2_82060_lb; i2 <= i2_94185_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                    const int i2_71474_lb = 32 * ii0 + 32 * ii1 + 32 * ii2 - i0 - i1 + 63, i2_8689_ub = 32 * ii0 - 32 * ii1 + 32 * ii2 - i0 + i1;
                    for (register int i2 = i2_71474_lb; i2 < i2_8689_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  } else {
                    const int i2_54513_lb = 0, i2_43698_ub = min(-32 * ii0 - 32 * ii1 + i0 + i1 - 32, -32 * ii0 + 32 * ii1 + i0 - i1 + 31);
                    for (register int i2 = i2_54513_lb; i2 <= i2_43698_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                    const int i2_21264_lb = max(32 * ii0 - 32 * ii1 - i0 + i1 + 1536, 32 * ii0 + 32 * ii1 - i0 - i1 + 1599), i2_56177_ub = 1599;
                    for (register int i2 = i2_21264_lb; i2 <= i2_56177_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                }
                if (ii2 == 0 && i0 >= 32 * ii0 + 16) {
                  const int i1_2409_lb = -32 * ii0 + 32 * ii1 + i0 + 17, i1_87486_ub = 32 * ii0 + 32 * ii1 - i0 + 62;
                  for (register int i1 = i1_2409_lb; i1 <= i1_87486_ub; i1 += 1) {
                    const int i2_51874_lb = -32 * ii0 + 32 * ii1 + i0 - i1 + 32, i2_33930_ub = 32 * ii0 - 32 * ii1 - i0 + i1;
                    for (register int i2 = i2_51874_lb; i2 < i2_33930_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                } else if (ii2 == 0 && 32 * ii0 + 15 >= i0) {
                  const int i1_33912_lb = 32 * ii0 + 32 * ii1 - i0 + 47, i1_80346_ub = -32 * ii0 + 32 * ii1 + i0 + 32;
                  for (register int i1 = i1_33912_lb; i1 <= i1_80346_ub; i1 += 1) {
                    const int i2_66533_lb = -32 * ii0 + 32 * ii1 + i0 - i1 + 32, i2_99211_ub = 32 * ii0 - 32 * ii1 - i0 + i1;
                    for (register int i2 = i2_66533_lb; i2 < i2_99211_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                }
              }
            }
          } else if (ii1 >= 1) {
            const int ii2_96392_lb = max(ii1 - 2, -ii1 + 1), ii2_24025_ub = 48;
            for (register int ii2 = ii2_96392_lb; ii2 <= ii2_24025_ub; ii2 += 1) {
              const int i0_6004_lb = 32 * ii0, i0_33386_ub = min(min(32 * ii0 + 30, 32 * ii0 - 15 * ii1 + 16 * ii2 + 45), 32 * ii0 + 15 * ii1 + 16 * ii2);
              for (register int i0 = i0_6004_lb; i0 <= i0_33386_ub; i0 += 1) {
                if (ii1 == 1) {
                  const int i1_16257_lb = max(-32 * ii0 + i0 + 33, 32 * ii0 - i0 + 63), i1_45672_ub = 63;
                  for (register int i1 = i1_16257_lb; i1 <= i1_45672_ub; i1 += 1) {
                    const int i2_94869_lb = max(-32 * ii0 + i0 + i1 - 63, 32 * ii0 + 32 * ii2 - i0 - i1 + 63), i2_1147_ub = 32 * ii0 + 32 * ii2 - i0 - i1 + 94;
                    for (register int i2 = i2_94869_lb; i2 <= i2_1147_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                }
                const int i1_89395_lb = max(max(64, -32 * ii0 + 31 * ii1 + i0 + 3), 32 * ii0 + 31 * ii1 - i0 + 33), i1_1965_ub = min(min(95, -32 * ii0 + 31 * ii1 + i0 + 33), 32 * ii0 + 31 * ii1 - i0 + 63);
                for (register int i1 = i1_89395_lb; i1 <= i1_1965_ub; i1 += 1) {
                  const int i2_48_lb = max(max(max(-32 * ii0 + i0 + i1 - 95, 32 * ii0 - 31 * ii1 + 32 * ii2 - i0 + i1 - 33), 32 * ii0 + 31 * ii1 + 32 * ii2 - i0 - i1 + 33), -32 * ii0 + i0 - i1 + 64), i2_81944_ub = -32 * ii0 + 31 * ii1 + i0 - i1 + 33;
                  for (register int i2 = i2_48_lb; i2 <= i2_81944_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  if (ii1 == 1) {
                    const int i2_54437_lb = max(32 * ii0 + 32 * ii2 - i0 + i1 - 64, -32 * ii0 + i0 - i1 + 65), i2_96233_ub = -32 * ii0 + i0 + i1 - 63;
                    for (register int i2 = i2_54437_lb; i2 < i2_96233_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                  const int i2_51859_lb = max(32 * ii0 + 32 * ii2 - i0 + i1 - 64, -32 * ii0 + i0 + i1 - 63), i2_36497_ub = 32 * ii0 + 32 * ii2 - i0 - i1 + 94;
                  for (register int i2 = i2_51859_lb; i2 <= i2_36497_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  const int i2_6770_lb = max(32 * ii0 + 32 * ii2 - i0 - i1 + 95, -32 * ii0 + i0 - i1 + 96), i2_39685_ub = 32 * ii0 + 32 * ii2 - i0 + i1 - 64;
                  for (register int i2 = i2_6770_lb; i2 < i2_39685_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  const int i2_45186_lb = max(max(max(32 * ii0 + 32 * ii2 - i0 + i1 - 64, -32 * ii0 - 31 * ii1 + i0 + i1 - 32), -32 * ii0 + 31 * ii1 + i0 - i1 + 34), 32 * ii0 + 32 * ii2 - i0 - i1 + 95), i2_61283_ub = min(32 * ii0 + 32 * ii2 - i0 + i1 - 33, 32 * ii0 + 32 * ii2 - i0 - i1 + 126);
                  for (register int i2 = i2_45186_lb; i2 <= i2_61283_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
                if (ii1 == 2) {
                  const int i1_83384_lb = 96, i1_82802_ub = min(-32 * ii0 + i0 + 96, 32 * ii0 - i0 + 126);
                  for (register int i1 = i1_83384_lb; i1 <= i1_82802_ub; i1 += 1) {
                    const int i2_33812_lb = max(32 * ii0 + 32 * ii2 - i0 + i1 - 96, -32 * ii0 + i0 - i1 + 96), i2_2145_ub = 32 * ii0 + 32 * ii2 - i0 + i1 - 64;
                    for (register int i2 = i2_33812_lb; i2 < i2_2145_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                }
              }
              if (ii1 == 1 && ii2 == 0) {
                const int i0_86640_lb = 32 * ii0 + 16, i0_85686_ub = 32 * ii0 + 22;
                for (register int i0 = i0_86640_lb; i0 <= i0_85686_ub; i0 += 1) {
                  const int i1_36076_lb = -32 * ii0 + i0 + 33, i1_20552_ub = 32 * ii0 - i0 + 78;
                  for (register int i1 = i1_36076_lb; i1 <= i1_20552_ub; i1 += 1) {
                    const int i2_66032_lb = -32 * ii0 + i0 + i1 - 63, i2_2609_ub = 32 * ii0 - i0 - i1 + 94;
                    for (register int i2 = i2_66032_lb; i2 <= i2_2609_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                  const int i1_19764_lb = -32 * ii0 + i0 + 49, i1_62424_ub = 32 * ii0 - i0 + 94;
                  for (register int i1 = i1_19764_lb; i1 <= i1_62424_ub; i1 += 1) {
                    const int i2_42986_lb = -32 * ii0 + i0 - i1 + 64, i2_25768_ub = 32 * ii0 - i0 + i1 - 32;
                    for (register int i2 = i2_42986_lb; i2 < i2_25768_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                }
              } else if (ii1 == 2 && ii2 == 0) {
                const int i0_95811_lb = 32 * ii0 + 16, i0_59244_ub = 32 * ii0 + 22;
                for (register int i0 = i0_95811_lb; i0 <= i0_59244_ub; i0 += 1) {
                  const int i1_87792_lb = -32 * ii0 + i0 + 65, i1_7032_ub = 32 * ii0 - i0 + 110;
                  for (register int i1 = i1_87792_lb; i1 <= i1_7032_ub; i1 += 1) {
                    const int i2_76743_lb = -32 * ii0 + i0 + i1 - 95, i2_77188_ub = 32 * ii0 - i0 - i1 + 126;
                    for (register int i2 = i2_76743_lb; i2 <= i2_77188_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                  const int i1_8997_lb = -32 * ii0 + i0 + 81, i1_76791_ub = 32 * ii0 - i0 + 126;
                  for (register int i1 = i1_8997_lb; i1 <= i1_76791_ub; i1 += 1) {
                    const int i2_75484_lb = -32 * ii0 + i0 - i1 + 96, i2_79786_ub = 32 * ii0 - i0 + i1 - 64;
                    for (register int i2 = i2_75484_lb; i2 < i2_79786_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                }
              }
            }
            const int i0_73024_lb = 32 * ii0, i0_27344_ub = 32 * ii0 + 30;
            for (register int i0 = i0_73024_lb; i0 <= i0_27344_ub; i0 += 1) {
              const int i1_16283_lb = max(-32 * ii0 + 32 * ii1 + i0 + 1, 32 * ii0 + 32 * ii1 - i0 + 31), i1_96146_ub = min(-32 * ii0 + 32 * ii1 + i0 + 32, 32 * ii0 + 32 * ii1 - i0 + 62);
              for (register int i1 = i1_16283_lb; i1 <= i1_96146_ub; i1 += 1) {
                const int i2_67029_lb = 0, i2_77821_ub = min(-32 * ii0 + i0 + i1 - 64, -32 * ii0 + i0 - i1 + 63);
                for (register int i2 = i2_67029_lb; i2 <= i2_77821_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
                const int i2_57429_lb = 0, i2_66765_ub = min(-32 * ii0 + i0 + i1 - 96, -32 * ii0 + i0 - i1 + 95);
                for (register int i2 = i2_57429_lb; i2 <= i2_66765_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
                const int i2_60623_lb = max(32 * ii0 - 32 * ii1 - i0 + i1 + 1536, 32 * ii0 + 32 * ii1 - i0 - i1 + 1599), i2_7594_ub = 1599;
                for (register int i2 = i2_60623_lb; i2 <= i2_7594_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
            }
          } else {
            for (register int ii2 = 0; ii2 <= 49; ii2 += 1) {
              const int i0_93280_lb = 32 * ii0, i0_4987_ub = min(32 * ii0 + 30, 32 * ii0 + 8 * ii2 + 22);
              for (register int i0 = i0_93280_lb; i0 <= i0_4987_ub; i0 += 1) {
                const int i1_84168_lb = max(-32 * ii0 + i0 + 1, 32 * ii0 - i0 + 31), i1_75665_ub = min(min(-32 * ii0 + i0 + 32, 32 * ii0 + 16 * ii2 - i0 + 46), 32 * ii0 - i0 + 62);
                for (register int i1 = i1_84168_lb; i1 <= i1_75665_ub; i1 += 1) {
                  if (ii2 <= 48) {
                    if (ii2 == 1) {
                      const int i2_23948_lb = max(32 * ii0 - i0 + i1, -32 * ii0 + i0 - i1 + 32), i2_3932_ub = -32 * ii0 + i0 + i1 - 31;
                      for (register int i2 = i2_23948_lb; i2 < i2_3932_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    } else if (ii2 == 0) {
                      const int i2_38089_lb = -32 * ii0 + i0 - i1 + 32, i2_66935_ub = -32 * ii0 + i0 + i1 - 31;
                      for (register int i2 = i2_38089_lb; i2 < i2_66935_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                    const int i2_46052_lb = max(max(32 * ii0 + 32 * ii2 - i0 + i1 - 32, -32 * ii0 + i0 + i1 - 31), 32 * ii0 + 32 * ii2 - i0 - i1 + 31), i2_50252_ub = 32 * ii0 + 32 * ii2 - i0 - i1 + 62;
                    for (register int i2 = i2_46052_lb; i2 <= i2_50252_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                    const int i2_42531_lb = 32 * ii0 + 32 * ii2 - i0 - i1 + 63, i2_50197_ub = 32 * ii0 + 32 * ii2 - i0 + i1;
                    for (register int i2 = i2_42531_lb; i2 < i2_50197_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  } else {
                    const int i2_73636_lb = 0, i2_35626_ub = min(-32 * ii0 + i0 + i1 - 32, -32 * ii0 + i0 - i1 + 31);
                    for (register int i2 = i2_73636_lb; i2 <= i2_35626_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                    const int i2_27385_lb = max(32 * ii0 - i0 + i1 + 1536, 32 * ii0 - i0 - i1 + 1599), i2_98985_ub = 1599;
                    for (register int i2 = i2_27385_lb; i2 <= i2_98985_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                }
                if (ii2 == 0) {
                  const int i1_12417_lb = max(-32 * ii0 + i0 + 17, 32 * ii0 - i0 + 47), i1_19221_ub = min(-32 * ii0 + i0 + 32, 32 * ii0 - i0 + 62);
                  for (register int i1 = i1_12417_lb; i1 <= i1_19221_ub; i1 += 1) {
                    const int i2_78771_lb = -32 * ii0 + i0 - i1 + 32, i2_1794_ub = min(-32 * ii0 + i0 + i1 - 31, 32 * ii0 - i0 + i1);
                    for (register int i2 = i2_78771_lb; i2 < i2_1794_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                    const int i2_46565_lb = -32 * ii0 + i0 + i1 - 31, i2_95054_ub = 32 * ii0 - i0 + i1;
                    for (register int i2 = i2_46565_lb; i2 < i2_95054_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                }
              }
              {
                const int i0_97940_lb = 32 * ii0, i0_29947_ub = min(min(32 * ii0 + 30, 16 * ii0 + 248), 32 * ii0 + 8 * ii2 + 22);
                for (register int i0 = i0_97940_lb; i0 <= i0_29947_ub; i0 += 1) {
                  if (ii0 == 14) {
                    const int i1_89227_lb = max(0, -16 * ii2 + i0 - 463), i1_71722_ub = i0 - 467;
                    for (register int i1 = i1_89227_lb; i1 < i1_71722_ub; i1 += 1) {
                      if (ii2 <= 48) {
                        const int i2_96712_lb = max(i0 - i1 - 448, 32 * ii2 - i0 + i1 + 448), i2_66202_ub = 32 * ii2 - i0 + i1 + 479;
                        for (register int i2 = i2_96712_lb; i2 <= i2_66202_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                      } else {
                        const int i2_79316_lb = 0, i2_81975_ub = i0 - i1 - 448;
                        for (register int i2 = i2_79316_lb; i2 < i2_81975_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                        const int i2_29818_lb = -i0 + i1 + 2016, i2_88948_ub = 1599;
                        for (register int i2 = i2_29818_lb; i2 <= i2_88948_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                      }
                    }
                  }
                  const int i1_86962_lb = max(max(0, i0 - 467), -32 * ii0 - 16 * ii2 + i0 - 15), i1_13986_ub = min(-32 * ii0 + i0, 32 * ii0 - i0 + 30);
                  for (register int i1 = i1_86962_lb; i1 <= i1_13986_ub; i1 += 1) {
                    if (ii2 <= 48) {
                      const int i2_64613_lb = max(-32 * ii0 + i0 - i1, 32 * ii0 + 32 * ii2 - i0 + i1), i2_27263_ub = 32 * ii0 + 32 * ii2 - i0 + i1 + 31;
                      for (register int i2 = i2_64613_lb; i2 <= i2_27263_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    } else {
                      const int i2_17918_lb = 0, i2_2703_ub = -32 * ii0 + i0 - i1;
                      for (register int i2 = i2_17918_lb; i2 < i2_2703_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                      const int i2_10550_lb = 32 * ii0 - i0 + i1 + 1568, i2_80322_ub = 1599;
                      for (register int i2 = i2_10550_lb; i2 <= i2_80322_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                  }
                  const int i1_69307_lb = max(-32 * ii0 + i0 + 1569, 32 * ii0 - i0 + 1599), i1_53081_ub = min(min(1599, 32 * ii0 + 16 * ii2 - i0 + 1614), -i0 + 2066);
                  for (register int i1 = i1_69307_lb; i1 <= i1_53081_ub; i1 += 1) {
                    if (ii2 <= 48) {
                      const int i2_30519_lb = max(-32 * ii0 + i0 + i1 - 1599, 32 * ii0 + 32 * ii2 - i0 - i1 + 1599), i2_42944_ub = 32 * ii0 + 32 * ii2 - i0 - i1 + 1630;
                      for (register int i2 = i2_30519_lb; i2 <= i2_42944_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    } else {
                      const int i2_5059_lb = 0, i2_74256_ub = -32 * ii0 + i0 + i1 - 1599;
                      for (register int i2 = i2_5059_lb; i2 < i2_74256_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                      const int i2_41929_lb = 32 * ii0 - i0 - i1 + 3167, i2_17476_ub = 1599;
                      for (register int i2 = i2_41929_lb; i2 <= i2_17476_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                  }
                  if (ii0 == 14) {
                    const int i1_93478_lb = -i0 + 2067, i1_20701_ub = min(1599, 16 * ii2 - i0 + 2062);
                    for (register int i1 = i1_93478_lb; i1 <= i1_20701_ub; i1 += 1) {
                      if (ii2 <= 48) {
                        const int i2_35622_lb = max(i0 + i1 - 2047, 32 * ii2 - i0 - i1 + 2047), i2_40043_ub = 32 * ii2 - i0 - i1 + 2078;
                        for (register int i2 = i2_35622_lb; i2 <= i2_40043_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                      } else {
                        const int i2_32107_lb = 0, i2_33563_ub = i0 + i1 - 2047;
                        for (register int i2 = i2_32107_lb; i2 < i2_33563_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                        const int i2_86342_lb = -i0 - i1 + 3615, i2_21335_ub = 1599;
                        for (register int i2 = i2_86342_lb; i2 <= i2_21335_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                      }
                    }
                  }
                }
                if (ii0 == 14 && ii2 >= 1) {
                  for (register int i0 = 473; i0 <= 478; i0 += 1) {
                    const int i1_3889_lb = 0, i1_17305_ub = -i0 + 478;
                    for (register int i1 = i1_3889_lb; i1 <= i1_17305_ub; i1 += 1) {
                      if (ii2 <= 48) {
                        const int i2_65030_lb = max(i0 - i1 - 448, 32 * ii2 - i0 + i1 + 448), i2_33707_ub = 32 * ii2 - i0 + i1 + 479;
                        for (register int i2 = i2_65030_lb; i2 <= i2_33707_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                      } else {
                        const int i2_6253_lb = 0, i2_68345_ub = i0 - i1 - 448;
                        for (register int i2 = i2_6253_lb; i2 < i2_68345_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                        const int i2_47693_lb = -i0 + i1 + 2016, i2_70867_ub = 1599;
                        for (register int i2 = i2_47693_lb; i2 <= i2_70867_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                      }
                    }
                    const int i1_95608_lb = i0 + 1121, i1_81963_ub = 1599;
                    for (register int i1 = i1_95608_lb; i1 <= i1_81963_ub; i1 += 1) {
                      if (ii2 <= 48) {
                        const int i2_73570_lb = max(i0 + i1 - 2047, 32 * ii2 - i0 - i1 + 2047), i2_22510_ub = 32 * ii2 - i0 - i1 + 2078;
                        for (register int i2 = i2_73570_lb; i2 <= i2_22510_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                      } else {
                        const int i2_78638_lb = 0, i2_59229_ub = i0 + i1 - 2047;
                        for (register int i2 = i2_78638_lb; i2 < i2_59229_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                        const int i2_75591_lb = -i0 - i1 + 3615, i2_9157_ub = 1599;
                        for (register int i2 = i2_75591_lb; i2 <= i2_9157_ub; i2 += 1) {
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
        for (register int ii2 = 0; ii2 <= 49; ii2 += 1) {
          const int i0_99766_lb = 32 * ii0, i0_60455_ub = min(32 * ii0 + 30, 32 * ii0 + 8 * ii2 + 22);
          for (register int i0 = i0_99766_lb; i0 <= i0_60455_ub; i0 += 1) {
            if (i0 >= 32 * ii0 + 16 * ii2) {
              const int i1_14478_lb = max(-32 * ii0 + i0 + 1537, 32 * ii0 - i0 + 1567), i1_93244_ub = min(1567, 32 * ii0 + 15 * ii2 - i0 + 1582);
              for (register int i1 = i1_14478_lb; i1 <= i1_93244_ub; i1 += 1) {
                const int i2_81156_lb = max(-32 * ii0 + i0 + i1 - 1567, 32 * ii0 + 32 * ii2 - i0 - i1 + 1567), i2_66453_ub = min(-32 * ii0 + i0 - i1 + 1567, 32 * ii0 + 29 * ii2 - i0 - i1 + 1598);
                for (register int i2 = i2_81156_lb; i2 <= i2_66453_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
                const int i2_49639_lb = -32 * ii0 + i0 - i1 + 1568, i2_29615_ub = 32 * ii0 + 32 * ii2 - i0 - i1 + 1598;
                for (register int i2 = i2_49639_lb; i2 <= i2_29615_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
              const int i1_16368_lb = max(1568, -32 * ii0 - 16 * ii2 + i0 + 1553), i1_35982_ub = min(-32 * ii0 + i0 + 1568, 32 * ii0 - i0 + 1598);
              for (register int i1 = i1_16368_lb; i1 <= i1_35982_ub; i1 += 1) {
                const int i2_67302_lb = max(32 * ii0 + 32 * ii2 - i0 + i1 - 1568, -32 * ii0 + i0 - i1 + 1568), i2_38005_ub = 32 * ii0 + 32 * ii2 - i0 + i1 - 1536;
                for (register int i2 = i2_67302_lb; i2 < i2_38005_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
            } else if (ii2 <= 48) {
              const int i1_35389_lb = max(-32 * ii0 + i0 + 1537, 32 * ii0 - i0 + 1567), i1_71192_ub = 1567;
              for (register int i1 = i1_35389_lb; i1 <= i1_71192_ub; i1 += 1) {
                const int i2_55310_lb = 32 * ii0 + 32 * ii2 - i0 - i1 + 1567, i2_16771_ub = 32 * ii0 + 32 * ii2 - i0 - i1 + 1598;
                for (register int i2 = i2_55310_lb; i2 <= i2_16771_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
              const int i1_4899_lb = 1568, i1_61563_ub = min(-32 * ii0 + i0 + 1568, 32 * ii0 - i0 + 1598);
              for (register int i1 = i1_4899_lb; i1 <= i1_61563_ub; i1 += 1) {
                const int i2_85116_lb = 32 * ii0 + 32 * ii2 - i0 + i1 - 1568, i2_68945_ub = 32 * ii0 + 32 * ii2 - i0 + i1 - 1536;
                for (register int i2 = i2_85116_lb; i2 < i2_68945_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
            } else {
              const int i1_48782_lb = max(-32 * ii0 + i0 + 1537, 32 * ii0 - i0 + 1567), i1_97076_ub = min(-32 * ii0 + i0 + 1568, 32 * ii0 - i0 + 1598);
              for (register int i1 = i1_48782_lb; i1 <= i1_97076_ub; i1 += 1) {
                const int i2_67260_lb = 0, i2_38704_ub = min(-32 * ii0 + i0 + i1 - 1568, -32 * ii0 + i0 - i1 + 1567);
                for (register int i2 = i2_67260_lb; i2 <= i2_38704_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
                const int i2_19586_lb = max(32 * ii0 - i0 + i1, 32 * ii0 - i0 - i1 + 3135), i2_45898_ub = 1599;
                for (register int i2 = i2_19586_lb; i2 <= i2_45898_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
            }
          }
        }
      }
    }
  } else {
    for (register int k = 0; k <= 1; k += 1) {
      const int ii1_71408_lb = 0, ii1_107_ub = k + 47;
      #pragma omp parallel for
      for (register int ii1 = ii1_71408_lb; ii1 <= ii1_107_ub; ii1 += 1) {
        for (register int ii2 = 0; ii2 <= 49; ii2 += 1) {
          if (k == 1) {
            if (ii1 >= 1) {
              const int i0_60562_lb = 0, i0_39362_ub = min(30, 8 * ii2 + 22);
              for (register int i0 = i0_60562_lb; i0 <= i0_39362_ub; i0 += 1) {
                if (ii1 <= 47 && ii2 <= 48) {
                  const int i1_64418_lb = max(32 * ii1 + i0 + 1, 32 * ii1 - i0 + 31), i1_58070_ub = min(min(32 * ii1 + i0 + 32, 32 * ii1 + 16 * ii2 - i0 + 46), 32 * ii1 - i0 + 62);
                  for (register int i1 = i1_64418_lb; i1 <= i1_58070_ub; i1 += 1) {
                    const int i2_5815_lb = max(-32 * ii1 + 32 * ii2 - i0 + i1 - 32, 32 * ii1 + i0 - i1 + 32), i2_14057_ub = -32 * ii1 + i0 + i1 - 31;
                    for (register int i2 = i2_5815_lb; i2 < i2_14057_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                    const int i2_4038_lb = max(max(-32 * ii1 + 32 * ii2 - i0 + i1 - 32, -32 * ii1 + i0 + i1 - 31), 32 * ii1 + 32 * ii2 - i0 - i1 + 31), i2_22183_ub = 32 * ii1 + 32 * ii2 - i0 - i1 + 62;
                    for (register int i2 = i2_4038_lb; i2 <= i2_22183_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                    const int i2_66391_lb = 32 * ii1 + 32 * ii2 - i0 - i1 + 63, i2_71340_ub = -32 * ii1 + 32 * ii2 - i0 + i1;
                    for (register int i2 = i2_66391_lb; i2 < i2_71340_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                  if (ii2 == 0 && i0 >= 16) {
                    const int i1_76540_lb = 32 * ii1 + i0 + 17, i1_1780_ub = 32 * ii1 - i0 + 62;
                    for (register int i1 = i1_76540_lb; i1 <= i1_1780_ub; i1 += 1) {
                      const int i2_58884_lb = 32 * ii1 + i0 - i1 + 32, i2_31850_ub = -32 * ii1 - i0 + i1;
                      for (register int i2 = i2_58884_lb; i2 < i2_31850_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                  } else if (ii2 == 0 && i0 <= 15) {
                    const int i1_18552_lb = 32 * ii1 - i0 + 47, i1_63784_ub = 32 * ii1 + i0 + 32;
                    for (register int i1 = i1_18552_lb; i1 <= i1_63784_ub; i1 += 1) {
                      const int i2_9765_lb = 32 * ii1 + i0 - i1 + 32, i2_20020_ub = -32 * ii1 - i0 + i1;
                      for (register int i2 = i2_9765_lb; i2 < i2_20020_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                  }
                } else if (ii2 == 49) {
                  const int i1_49081_lb = max(32 * ii1 + i0 + 1, 32 * ii1 - i0 + 31), i1_74900_ub = min(32 * ii1 + i0 + 32, 32 * ii1 - i0 + 62);
                  for (register int i1 = i1_49081_lb; i1 <= i1_74900_ub; i1 += 1) {
                    if (ii1 <= 47) {
                      const int i2_33449_lb = 0, i2_16341_ub = min(-32 * ii1 + i0 + i1 - 32, 32 * ii1 + i0 - i1 + 31);
                      for (register int i2 = i2_33449_lb; i2 <= i2_16341_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                    const int i2_13604_lb = 0, i2_69387_ub = min(i0 + i1 - 1568, i0 - i1 + 1567);
                    for (register int i2 = i2_13604_lb; i2 <= i2_69387_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                    const int i2_78592_lb = max(-32 * ii1 - i0 + i1 + 1536, 32 * ii1 - i0 - i1 + 1599), i2_11538_ub = 1599;
                    for (register int i2 = i2_78592_lb; i2 <= i2_11538_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                } else if (i0 >= 16 * ii2) {
                  const int i1_97269_lb = max(i0 + 1537, -i0 + 1567), i1_50000_ub = min(1567, 15 * ii2 - i0 + 1582);
                  for (register int i1 = i1_97269_lb; i1 <= i1_50000_ub; i1 += 1) {
                    const int i2_27998_lb = max(i0 + i1 - 1567, 32 * ii2 - i0 - i1 + 1567), i2_22152_ub = min(i0 - i1 + 1567, 29 * ii2 - i0 - i1 + 1598);
                    for (register int i2 = i2_27998_lb; i2 <= i2_22152_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                    const int i2_37526_lb = i0 - i1 + 1568, i2_4912_ub = 32 * ii2 - i0 - i1 + 1598;
                    for (register int i2 = i2_37526_lb; i2 <= i2_4912_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                  const int i1_61514_lb = max(1568, -16 * ii2 + i0 + 1553), i1_1944_ub = min(i0 + 1568, -i0 + 1598);
                  for (register int i1 = i1_61514_lb; i1 <= i1_1944_ub; i1 += 1) {
                    const int i2_79335_lb = max(32 * ii2 - i0 + i1 - 1568, i0 - i1 + 1568), i2_83681_ub = 32 * ii2 - i0 + i1 - 1536;
                    for (register int i2 = i2_79335_lb; i2 < i2_83681_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                } else {
                  const int i1_32353_lb = max(i0 + 1537, -i0 + 1567), i1_83373_ub = 1567;
                  for (register int i1 = i1_32353_lb; i1 <= i1_83373_ub; i1 += 1) {
                    const int i2_22216_lb = 32 * ii2 - i0 - i1 + 1567, i2_15097_ub = 32 * ii2 - i0 - i1 + 1598;
                    for (register int i2 = i2_22216_lb; i2 <= i2_15097_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                  const int i1_54713_lb = 1568, i1_98756_ub = min(i0 + 1568, -i0 + 1598);
                  for (register int i1 = i1_54713_lb; i1 <= i1_98756_ub; i1 += 1) {
                    const int i2_16877_lb = 32 * ii2 - i0 + i1 - 1568, i2_13598_ub = 32 * ii2 - i0 + i1 - 1536;
                    for (register int i2 = i2_16877_lb; i2 < i2_13598_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                }
              }
            } else if (ii2 <= 48) {
              const int i0_46958_lb = 0, i0_51781_ub = min(30, 8 * ii2 + 22);
              for (register int i0 = i0_46958_lb; i0 <= i0_51781_ub; i0 += 1) {
                const int i1_93734_lb = max(i0 + 1, -i0 + 31), i1_73076_ub = min(min(i0 + 32, 16 * ii2 - i0 + 46), -i0 + 62);
                for (register int i1 = i1_93734_lb; i1 <= i1_73076_ub; i1 += 1) {
                  if (ii2 == 0) {
                    const int i2_88154_lb = i0 - i1 + 32, i2_59167_ub = i0 + i1 - 31;
                    for (register int i2 = i2_88154_lb; i2 < i2_59167_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  } else if (ii2 == 1) {
                    const int i2_47976_lb = max(-i0 + i1, i0 - i1 + 32), i2_37955_ub = i0 + i1 - 31;
                    for (register int i2 = i2_47976_lb; i2 < i2_37955_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                  const int i2_91860_lb = max(max(32 * ii2 - i0 + i1 - 32, i0 + i1 - 31), 32 * ii2 - i0 - i1 + 31), i2_61580_ub = 32 * ii2 - i0 - i1 + 62;
                  for (register int i2 = i2_91860_lb; i2 <= i2_61580_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  const int i2_23694_lb = 32 * ii2 - i0 - i1 + 63, i2_70452_ub = 32 * ii2 - i0 + i1;
                  for (register int i2 = i2_23694_lb; i2 < i2_70452_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
                if (ii2 == 0) {
                  const int i1_73119_lb = max(i0 + 17, -i0 + 47), i1_20963_ub = min(i0 + 32, -i0 + 62);
                  for (register int i1 = i1_73119_lb; i1 <= i1_20963_ub; i1 += 1) {
                    const int i2_36804_lb = i0 - i1 + 32, i2_17469_ub = -i0 + i1;
                    for (register int i2 = i2_36804_lb; i2 < i2_17469_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                }
              }
            } else {
              for (register int i0 = 0; i0 <= 30; i0 += 1) {
                const int i1_22381_lb = max(i0 + 1, -i0 + 31), i1_20982_ub = min(i0 + 32, -i0 + 62);
                for (register int i1 = i1_22381_lb; i1 <= i1_20982_ub; i1 += 1) {
                  const int i2_92626_lb = 0, i2_1716_ub = min(i0 + i1 - 32, i0 - i1 + 31);
                  for (register int i2 = i2_92626_lb; i2 <= i2_1716_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  const int i2_21016_lb = max(-i0 + i1 + 1536, -i0 - i1 + 1599), i2_24980_ub = 1599;
                  for (register int i2 = i2_21016_lb; i2 <= i2_24980_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
              }
            }
          } else {
            for (register int i0 = 0; i0 <= 14; i0 += 1) {
              const int i1_56429_lb = 32 * ii1 + i0 + 33, i1_56155_ub = 32 * ii1 - i0 + 62;
              for (register int i1 = i1_56429_lb; i1 <= i1_56155_ub; i1 += 1) {
                if (ii2 <= 48) {
                  const int i2_41989_lb = max(0, 32 * ii2 - i0), i2_89658_ub = 32 * ii2;
                  for (register int i2 = i2_41989_lb; i2 < i2_89658_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  const int i2_69753_lb = max(32 * ii2, i0), i2_5299_ub = 32 * ii2 - i0 + 31;
                  for (register int i2 = i2_69753_lb; i2 <= i2_5299_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                } else {
                  for (register int i2 = 0; i2 < i0; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  const int i2_78375_lb = -i0 + 1568, i2_45946_ub = 1599;
                  for (register int i2 = i2_78375_lb; i2 <= i2_45946_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
              }
            }
            if (ii1 == 0) {
              for (register int i0 = 0; i0 <= 14; i0 += 1) {
                const int i1_253_lb = i0 + 1, i1_30866_ub = -i0 + 30;
                for (register int i1 = i1_253_lb; i1 <= i1_30866_ub; i1 += 1) {
                  if (ii2 <= 48) {
                    const int i2_4284_lb = max(0, 32 * ii2 - i0), i2_23947_ub = 32 * ii2;
                    for (register int i2 = i2_4284_lb; i2 < i2_23947_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                    const int i2_17671_lb = max(32 * ii2, i0), i2_93755_ub = 32 * ii2 - i0 + 31;
                    for (register int i2 = i2_17671_lb; i2 <= i2_93755_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  } else {
                    for (register int i2 = 0; i2 < i0; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                    const int i2_11224_lb = -i0 + 1568, i2_20731_ub = 1599;
                    for (register int i2 = i2_11224_lb; i2 <= i2_20731_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                }
              }
            }
          }
          if (ii1 == 0) {
            if (k == 1) {
              const int i0_28806_lb = 0, i0_33605_ub = min(30, 8 * ii2 + 22);
              for (register int i0 = i0_28806_lb; i0 <= i0_33605_ub; i0 += 1) {
                const int i1_58065_lb = max(0, -16 * ii2 + i0 - 15), i1_21432_ub = min(i0, -i0 + 30);
                for (register int i1 = i1_58065_lb; i1 <= i1_21432_ub; i1 += 1) {
                  if (ii2 <= 48) {
                    const int i2_35322_lb = max(i0 - i1, 32 * ii2 - i0 + i1), i2_79081_ub = 32 * ii2 - i0 + i1 + 31;
                    for (register int i2 = i2_35322_lb; i2 <= i2_79081_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  } else {
                    const int i2_62764_lb = 0, i2_36763_ub = i0 - i1;
                    for (register int i2 = i2_62764_lb; i2 < i2_36763_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                    const int i2_38666_lb = -i0 + i1 + 1568, i2_19193_ub = 1599;
                    for (register int i2 = i2_38666_lb; i2 <= i2_19193_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                }
                const int i1_92918_lb = max(i0 + 1569, -i0 + 1599), i1_97007_ub = min(1599, 16 * ii2 - i0 + 1614);
                for (register int i1 = i1_92918_lb; i1 <= i1_97007_ub; i1 += 1) {
                  if (ii2 <= 48) {
                    const int i2_8852_lb = max(i0 + i1 - 1599, 32 * ii2 - i0 - i1 + 1599), i2_62671_ub = 32 * ii2 - i0 - i1 + 1630;
                    for (register int i2 = i2_8852_lb; i2 <= i2_62671_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  } else {
                    const int i2_2306_lb = 0, i2_66644_ub = i0 + i1 - 1599;
                    for (register int i2 = i2_2306_lb; i2 < i2_66644_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                    const int i2_58862_lb = -i0 - i1 + 3167, i2_97034_ub = 1599;
                    for (register int i2 = i2_58862_lb; i2 <= i2_97034_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                }
              }
            } else {
              for (register int i0 = 0; i0 <= 14; i0 += 1) {
                const int i1_23385_lb = i0 + 1569, i1_45547_ub = -i0 + 1598;
                for (register int i1 = i1_23385_lb; i1 <= i1_45547_ub; i1 += 1) {
                  if (ii2 <= 48) {
                    const int i2_45087_lb = max(0, 32 * ii2 - i0), i2_44021_ub = 32 * ii2;
                    for (register int i2 = i2_45087_lb; i2 < i2_44021_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                    const int i2_85846_lb = max(32 * ii2, i0), i2_79110_ub = 32 * ii2 - i0 + 31;
                    for (register int i2 = i2_85846_lb; i2 <= i2_79110_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  } else {
                    for (register int i2 = 0; i2 < i0; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                    const int i2_33585_lb = -i0 + 1568, i2_49000_ub = 1599;
                    for (register int i2 = i2_33585_lb; i2 <= i2_49000_ub; i2 += 1) {
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
    #pragma omp parallel for
    for (register int ii1 = 0; ii1 <= 49; ii1 += 1) {
      for (register int ii2 = 0; ii2 <= 49; ii2 += 1) {
        if (ii1 >= 48) {
          const int i0_16528_lb = 16, i0_34280_ub = min(46, 8 * ii2 + 38);
          for (register int i0 = i0_16528_lb; i0 <= i0_34280_ub; i0 += 1) {
            if (ii1 == 49 && ii2 <= 48) {
              const int i1_37691_lb = max(i0 + 1537, -i0 + 1599), i1_95644_ub = min(1583, 16 * ii2 - i0 + 1614);
              for (register int i1 = i1_37691_lb; i1 <= i1_95644_ub; i1 += 1) {
                if (i0 >= 32) {
                  const int i2_71043_lb = max(i0 + i1 - 1599, 32 * ii2 - i0 - i1 + 1599), i2_92709_ub = 32 * ii2;
                  for (register int i2 = i2_71043_lb; i2 < i2_92709_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  const int i2_31190_lb = max(32 * ii2, i0 - i1 + 1568), i2_63962_ub = 32 * ii2 - i0 + i1 - 1536;
                  for (register int i2 = i2_31190_lb; i2 < i2_63962_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  const int i2_89716_lb = max(i0 + i1 - 1599, 32 * ii2 - i0 + i1 - 1536), i2_40042_ub = 32 * ii2 - i0 - i1 + 1630;
                  for (register int i2 = i2_89716_lb; i2 <= i2_40042_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                } else {
                  const int i2_42985_lb = max(i0 + i1 - 1599, 32 * ii2 - i0 - i1 + 1599), i2_8375_ub = 32 * ii2 - i0 - i1 + 1630;
                  for (register int i2 = i2_42985_lb; i2 <= i2_8375_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
              }
              const int i1_23038_lb = max(1584, -16 * ii2 + i0 + 1553), i1_1848_ub = min(i0 + 1568, -i0 + 1630);
              for (register int i1 = i1_23038_lb; i1 <= i1_1848_ub; i1 += 1) {
                if (i0 >= 32) {
                  const int i2_5409_lb = max(32 * ii2 - i0 + i1 - 1568, i0 - i1 + 1568), i2_68332_ub = 32 * ii2;
                  for (register int i2 = i2_5409_lb; i2 < i2_68332_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  const int i2_16068_lb = max(32 * ii2, i0 - i1 + 1568), i2_45146_ub = 32 * ii2 - i0 + i1 - 1536;
                  for (register int i2 = i2_16068_lb; i2 < i2_45146_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                } else {
                  const int i2_13879_lb = max(32 * ii2 - i0 + i1 - 1568, i0 - i1 + 1568), i2_77507_ub = 32 * ii2 - i0 + i1 - 1536;
                  for (register int i2 = i2_13879_lb; i2 < i2_77507_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
              }
            } else if (ii1 == 48 && i0 >= 16 * ii2 + 16) {
              const int i1_89168_lb = max(i0 + 1505, -i0 + 1567), i1_99725_ub = min(1551, 15 * ii2 - i0 + 1582);
              for (register int i1 = i1_89168_lb; i1 <= i1_99725_ub; i1 += 1) {
                const int i2_72969_lb = max(i0 + i1 - 1567, 32 * ii2 - i0 - i1 + 1567), i2_26944_ub = min(i0 - i1 + 1535, 29 * ii2 - i0 - i1 + 1598);
                for (register int i2 = i2_72969_lb; i2 <= i2_26944_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
                const int i2_63187_lb = i0 - i1 + 1536, i2_22907_ub = 32 * ii2 - i0 - i1 + 1598;
                for (register int i2 = i2_63187_lb; i2 <= i2_22907_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
              const int i1_92297_lb = max(1552, -16 * ii2 + i0 + 1521), i1_47379_ub = min(i0 + 1536, -i0 + 1598);
              for (register int i1 = i1_92297_lb; i1 <= i1_47379_ub; i1 += 1) {
                const int i2_1650_lb = max(32 * ii2 - i0 + i1 - 1536, i0 - i1 + 1536), i2_74903_ub = 32 * ii2 - i0 + i1 - 1504;
                for (register int i2 = i2_1650_lb; i2 < i2_74903_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
            } else if (ii1 == 48 && ii2 <= 48) {
              const int i1_5989_lb = max(i0 + 1505, -i0 + 1567), i1_34530_ub = 1551;
              for (register int i1 = i1_5989_lb; i1 <= i1_34530_ub; i1 += 1) {
                const int i2_9183_lb = 32 * ii2 - i0 - i1 + 1567, i2_60033_ub = 32 * ii2 - i0 - i1 + 1598;
                for (register int i2 = i2_9183_lb; i2 <= i2_60033_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
              const int i1_30175_lb = 1552, i1_96578_ub = min(i0 + 1536, -i0 + 1598);
              for (register int i1 = i1_30175_lb; i1 <= i1_96578_ub; i1 += 1) {
                const int i2_52742_lb = 32 * ii2 - i0 + i1 - 1536, i2_77717_ub = 32 * ii2 - i0 + i1 - 1504;
                for (register int i2 = i2_52742_lb; i2 < i2_77717_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
            } else if (ii1 == 49) {
              const int i1_60540_lb = max(i0 + 1537, -i0 + 1599), i1_58811_ub = min(i0 + 1568, -i0 + 1630);
              for (register int i1 = i1_60540_lb; i1 <= i1_58811_ub; i1 += 1) {
                const int i2_17759_lb = 0, i2_19878_ub = min(i0 + i1 - 1600, i0 - i1 + 1567);
                for (register int i2 = i2_17759_lb; i2 <= i2_19878_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
                const int i2_67186_lb = max(-i0 + i1, -i0 - i1 + 3167), i2_57149_ub = 1599;
                for (register int i2 = i2_67186_lb; i2 <= i2_57149_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
            } else {
              const int i1_38078_lb = max(i0 + 1505, -i0 + 1567), i1_88947_ub = min(i0 + 1536, -i0 + 1598);
              for (register int i1 = i1_38078_lb; i1 <= i1_88947_ub; i1 += 1) {
                const int i2_25481_lb = 0, i2_70498_ub = min(i0 + i1 - 1568, i0 - i1 + 1535);
                for (register int i2 = i2_25481_lb; i2 <= i2_70498_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
                const int i2_34093_lb = max(-i0 + i1 + 32, -i0 - i1 + 3135), i2_55712_ub = 1599;
                for (register int i2 = i2_34093_lb; i2 <= i2_55712_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
            }
          }
        } else if (ii1 >= 1 && ii2 <= 48) {
          for (register int i0 = 16; i0 <= 31; i0 += 1) {
            {
              const int i1_55437_lb = max(32 * ii1 - i0 + 31, i0 + 33), i1_37327_ub = 32 * ii1 + 15;
              for (register int i1 = i1_55437_lb; i1 <= i1_37327_ub; i1 += 1) {
                const int i2_66558_lb = max(-32 * ii1 + i0 + i1 - 31, 32 * ii1 + 32 * ii2 - i0 - i1 + 31), i2_18624_ub = 32 * ii1 + 32 * ii2 - i0 - i1 + 62;
                for (register int i2 = i2_66558_lb; i2 <= i2_18624_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
              const int i1_60234_lb = max(32 * ii1 + 16, i0 + 33), i1_58855_ub = 32 * ii1 + i0;
              for (register int i1 = i1_60234_lb; i1 <= i1_58855_ub; i1 += 1) {
                const int i2_82356_lb = max(-32 * ii1 + 32 * ii2 - i0 + i1, 32 * ii1 + i0 - i1), i2_78237_ub = -32 * ii1 + 32 * ii2 - i0 + i1 + 31;
                for (register int i2 = i2_82356_lb; i2 <= i2_78237_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
            }
            if (ii1 == 1) {
              const int i1_50110_lb = -i0 + 63, i1_88345_ub = 47;
              for (register int i1 = i1_50110_lb; i1 <= i1_88345_ub; i1 += 1) {
                const int i2_29119_lb = max(i0 + i1 - 63, 32 * ii2 - i0 - i1 + 63), i2_75645_ub = 32 * ii2 - i0 - i1 + 94;
                for (register int i2 = i2_29119_lb; i2 <= i2_75645_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
              const int i1_64730_lb = 48, i1_59294_ub = i0 + 32;
              for (register int i1 = i1_64730_lb; i1 <= i1_59294_ub; i1 += 1) {
                const int i2_72223_lb = max(32 * ii2 - i0 + i1 - 32, i0 - i1 + 32), i2_33825_ub = 32 * ii2 - i0 + i1;
                for (register int i2 = i2_72223_lb; i2 < i2_33825_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
            }
          }
          if (ii1 == 2 && ii2 >= 2) {
            for (register int i0 = 32; i0 <= 46; i0 += 1) {
              const int i1_92636_lb = i0 + 33, i1_71122_ub = -i0 + 126;
              for (register int i1 = i1_92636_lb; i1 <= i1_71122_ub; i1 += 1) {
                const int i2_68994_lb = max(32 * ii2 - i0 + i1 - 64, 32 * ii2 - i0 - i1 + 95), i2_76174_ub = 32 * ii2 - i0 + i1 - 32;
                for (register int i2 = i2_68994_lb; i2 < i2_76174_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
                const int i2_44623_lb = 32 * ii2 - i0 + i1 - 32, i2_7072_ub = 32 * ii2 - i0 - i1 + 126;
                for (register int i2 = i2_44623_lb; i2 <= i2_7072_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
            }
          } else if (ii1 >= 3) {
            const int i0_65121_lb = 32, i0_86456_ub = min(46, 8 * ii2 + 38);
            for (register int i0 = i0_65121_lb; i0 <= i0_86456_ub; i0 += 1) {
              const int i1_93922_lb = 32 * ii1 + i0 - 31, i1_15566_ub = min(32 * ii1 + 15, 32 * ii1 + 16 * ii2 - i0 + 46);
              for (register int i1 = i1_93922_lb; i1 <= i1_15566_ub; i1 += 1) {
                const int i2_42168_lb = max(-32 * ii1 + i0 + i1 - 31, 32 * ii1 + 32 * ii2 - i0 - i1 + 31), i2_58280_ub = 32 * ii1 + 32 * ii2 - i0 - i1 + 62;
                for (register int i2 = i2_42168_lb; i2 <= i2_58280_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
              const int i1_55180_lb = max(32 * ii1 + 16, 32 * ii1 - 16 * ii2 + i0 - 15), i1_13958_ub = 32 * ii1 - i0 + 62;
              for (register int i1 = i1_55180_lb; i1 <= i1_13958_ub; i1 += 1) {
                const int i2_95608_lb = max(-32 * ii1 + 32 * ii2 - i0 + i1, 32 * ii1 + i0 - i1), i2_21738_ub = -32 * ii1 + 32 * ii2 - i0 + i1 + 31;
                for (register int i2 = i2_95608_lb; i2 <= i2_21738_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
            }
          } else if (ii1 == 1 && ii2 == 0) {
            for (register int i0 = 32; i0 <= 38; i0 += 1) {
              const int i1_96945_lb = i0 + 1, i1_31290_ub = -i0 + 78;
              for (register int i1 = i1_96945_lb; i1 <= i1_31290_ub; i1 += 1) {
                const int i2_50431_lb = i0 + i1 - 63, i2_47055_ub = -i0 - i1 + 94;
                for (register int i2 = i2_50431_lb; i2 <= i2_47055_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
              const int i1_19636_lb = i0 + 17, i1_79551_ub = -i0 + 94;
              for (register int i1 = i1_19636_lb; i1 <= i1_79551_ub; i1 += 1) {
                const int i2_39052_lb = i0 - i1 + 32, i2_718_ub = -i0 + i1;
                for (register int i2 = i2_39052_lb; i2 < i2_718_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
            }
          } else if (ii1 == 2 && ii2 <= 1) {
            const int i0_55197_lb = 32, i0_11275_ub = 8 * ii2 + 38;
            for (register int i0 = i0_55197_lb; i0 <= i0_11275_ub; i0 += 1) {
              const int i1_34543_lb = i0 + 33, i1_8561_ub = min(79, 15 * ii2 - i0 + 110);
              for (register int i1 = i1_34543_lb; i1 <= i1_8561_ub; i1 += 1) {
                const int i2_76743_lb = max(i0 + i1 - 95, 28 * ii2 - i0 - i1 + 99), i2_43531_ub = min(i0 - i1 + 63, 29 * ii2 - i0 - i1 + 126);
                for (register int i2 = i2_76743_lb; i2 <= i2_43531_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
                if (ii2 == 1) {
                  const int i2_96035_lb = i0 - i1 + 64, i2_45737_ub = -i0 - i1 + 158;
                  for (register int i2 = i2_96035_lb; i2 <= i2_45737_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
              }
              const int i1_36057_lb = max(80, -16 * ii2 + i0 + 49), i1_40659_ub = -i0 + 126;
              for (register int i1 = i1_36057_lb; i1 <= i1_40659_ub; i1 += 1) {
                const int i2_69161_lb = max(32 * ii2 - i0 + i1 - 64, i0 - i1 + 64), i2_17530_ub = 32 * ii2 - i0 + i1 - 32;
                for (register int i2 = i2_69161_lb; i2 < i2_17530_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
            }
          } else {
            for (register int i0 = 32; i0 <= 46; i0 += 1) {
              const int i1_33097_lb = i0 + 1, i1_69284_ub = -i0 + 94;
              for (register int i1 = i1_33097_lb; i1 <= i1_69284_ub; i1 += 1) {
                if (ii2 == 1) {
                  const int i2_21364_lb = max(i0 + i1 - 63, -i0 - i1 + 95), i2_88277_ub = i0 - i1 + 31;
                  for (register int i2 = i2_21364_lb; i2 <= i2_88277_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
                const int i2_99594_lb = max(max(32 * ii2 - i0 + i1 - 32, i0 - i1 + 32), 32 * ii2 - i0 - i1 + 63), i2_33324_ub = 32 * ii2 - i0 + i1;
                for (register int i2 = i2_99594_lb; i2 < i2_33324_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
                const int i2_26367_lb = 32 * ii2 - i0 + i1, i2_32176_ub = 32 * ii2 - i0 - i1 + 94;
                for (register int i2 = i2_26367_lb; i2 <= i2_32176_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
            }
          }
        } else if (ii1 <= 1 && ii2 == 49) {
          for (register int i0 = 16; i0 <= 46; i0 += 1) {
            if (ii1 == 1) {
              const int i1_63467_lb = max(i0 + 1, -i0 + 63), i1_55950_ub = min(i0 + 32, -i0 + 94);
              for (register int i1 = i1_63467_lb; i1 <= i1_55950_ub; i1 += 1) {
                const int i2_86719_lb = 0, i2_99455_ub = min(i0 + i1 - 64, i0 - i1 + 31);
                for (register int i2 = i2_86719_lb; i2 <= i2_99455_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
                const int i2_35501_lb = max(-i0 + i1 + 1536, -i0 - i1 + 1631), i2_25771_ub = 1599;
                for (register int i2 = i2_35501_lb; i2 <= i2_25771_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
            } else {
              const int i1_16525_lb = max(i0 - 31, -i0 + 31), i1_90699_ub = min(i0, -i0 + 62);
              for (register int i1 = i1_16525_lb; i1 <= i1_90699_ub; i1 += 1) {
                const int i2_53398_lb = 0, i2_51069_ub = min(i0 + i1 - 31, i0 - i1);
                for (register int i2 = i2_53398_lb; i2 < i2_51069_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
                const int i2_15612_lb = max(-i0 + i1 + 1568, -i0 - i1 + 1599), i2_30142_ub = 1599;
                for (register int i2 = i2_15612_lb; i2 <= i2_30142_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
            }
          }
        } else if (ii1 == 0 && ii2 <= 48) {
          for (register int i0 = 16; i0 <= 31; i0 += 1) {
            if (ii2 == 0) {
              const int i1_92231_lb = -i0 + 31, i1_47010_ub = i0 - 15;
              for (register int i1 = i1_92231_lb; i1 < i1_47010_ub; i1 += 1) {
                const int i2_52306_lb = i0 + i1 - 31, i2_77745_ub = -i0 - i1 + 62;
                for (register int i2 = i2_52306_lb; i2 <= i2_77745_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
            }
            const int i1_64540_lb = max(-16 * ii2 + i0 - 15, -i0 + 31), i1_95774_ub = i0;
            for (register int i1 = i1_64540_lb; i1 <= i1_95774_ub; i1 += 1) {
              if (ii2 == 0) {
                const int i2_57181_lb = i0 + i1 - 31, i2_97637_ub = i0 - i1;
                for (register int i2 = i2_57181_lb; i2 < i2_97637_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
              const int i2_81410_lb = max(max(i0 - i1, 32 * ii2 - i0 + i1), 32 * ii2 - i0 - i1 + 31), i2_78545_ub = 32 * ii2 - i0 + i1 + 31;
              for (register int i2 = i2_81410_lb; i2 <= i2_78545_ub; i2 += 1) {
                S1(i0, i1, i2);
              }
              const int i2_2266_lb = 32 * ii2 - i0 + i1 + 32, i2_81004_ub = 32 * ii2 - i0 - i1 + 62;
              for (register int i2 = i2_2266_lb; i2 <= i2_81004_ub; i2 += 1) {
                S1(i0, i1, i2);
              }
            }
          }
          const int i0_28222_lb = 32, i0_28633_ub = min(46, 7 * ii2 + (ii2 + 3) / 4 + 38);
          for (register int i0 = i0_28222_lb; i0 <= i0_28633_ub; i0 += 1) {
            if (ii2 >= 2) {
              const int i1_29532_lb = i0 - 31, i1_33741_ub = 15;
              for (register int i1 = i1_29532_lb; i1 <= i1_33741_ub; i1 += 1) {
                const int i2_68297_lb = 32 * ii2 - i0 - i1 + 31, i2_92999_ub = 32 * ii2 - i0 - i1 + 62;
                for (register int i2 = i2_68297_lb; i2 <= i2_92999_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
            } else {
              const int i1_89691_lb = i0 - 31, i1_71368_ub = min(15, 15 * ii2 - i0 + 46);
              for (register int i1 = i1_89691_lb; i1 <= i1_71368_ub; i1 += 1) {
                const int i2_8806_lb = max(i0 + i1 - 31, 28 * ii2 - i0 - i1 + 35), i2_41545_ub = min(i0 - i1 - 1, 29 * ii2 - i0 - i1 + 62);
                for (register int i2 = i2_8806_lb; i2 <= i2_41545_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
                if (ii2 == 1) {
                  const int i2_13491_lb = i0 - i1, i2_25332_ub = -i0 - i1 + 94;
                  for (register int i2 = i2_13491_lb; i2 <= i2_25332_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
              }
            }
            const int i1_32244_lb = max(16, -16 * ii2 + i0 - 15), i1_83242_ub = -i0 + 62;
            for (register int i1 = i1_32244_lb; i1 <= i1_83242_ub; i1 += 1) {
              const int i2_92753_lb = max(i0 - i1, 32 * ii2 - i0 + i1), i2_47856_ub = 32 * ii2 - i0 + i1 + 31;
              for (register int i2 = i2_92753_lb; i2 <= i2_47856_ub; i2 += 1) {
                S1(i0, i1, i2);
              }
            }
          }
        } else {
          for (register int i0 = 16; i0 <= 46; i0 += 1) {
            if (ii1 >= 3 && i0 >= 32) {
              const int i1_75855_lb = 32 * ii1 + i0 - 31, i1_21967_ub = 32 * ii1 - i0 + 62;
              for (register int i1 = i1_75855_lb; i1 <= i1_21967_ub; i1 += 1) {
                const int i2_50715_lb = 0, i2_28162_ub = min(-32 * ii1 + i0 + i1 - 31, 32 * ii1 + i0 - i1);
                for (register int i2 = i2_50715_lb; i2 < i2_28162_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
                const int i2_99712_lb = max(-32 * ii1 - i0 + i1 + 1568, 32 * ii1 - i0 - i1 + 1599), i2_15256_ub = 1599;
                for (register int i2 = i2_99712_lb; i2 <= i2_15256_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
            } else if (i0 <= 31) {
              const int i1_23936_lb = 32 * ii1 - i0 + 31, i1_56893_ub = 32 * ii1 + i0;
              for (register int i1 = i1_23936_lb; i1 <= i1_56893_ub; i1 += 1) {
                const int i2_29245_lb = 0, i2_5346_ub = min(-32 * ii1 + i0 + i1 - 31, 32 * ii1 + i0 - i1);
                for (register int i2 = i2_29245_lb; i2 < i2_5346_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
                const int i2_51791_lb = max(-32 * ii1 - i0 + i1 + 1568, 32 * ii1 - i0 - i1 + 1599), i2_31512_ub = 1599;
                for (register int i2 = i2_51791_lb; i2 <= i2_31512_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
            } else {
              const int i1_2702_lb = i0 + 33, i1_96365_ub = -i0 + 126;
              for (register int i1 = i1_2702_lb; i1 <= i1_96365_ub; i1 += 1) {
                const int i2_76497_lb = 0, i2_32234_ub = min(i0 + i1 - 96, i0 - i1 + 63);
                for (register int i2 = i2_76497_lb; i2 <= i2_32234_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
                const int i2_30106_lb = max(-i0 + i1 + 1504, -i0 - i1 + 1663), i2_61147_ub = 1599;
                for (register int i2 = i2_30106_lb; i2 <= i2_61147_ub; i2 += 1) {
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
    total = 0;
    for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++) {
        total += A[T % 2][i][j];
      }
    }
    fprintf(stderr, "|sum: %e\t", total);
    for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++) {
        sum_err_sqr +=
            (A[T % 2][i][j] - (total / N)) * (A[T % 2][i][j] - (total / N));
      }
    }
    fprintf(stderr, "|rms(A) = %7.2f\t", sqrt(sum_err_sqr));
    for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++) {
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
