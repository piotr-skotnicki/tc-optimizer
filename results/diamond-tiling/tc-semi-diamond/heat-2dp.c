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
/* ./tc ../examples/pluto+/heat-2dp-nt-t.scop.c --semi-diamond-tiling --omp-for-codegen --iterative-tc --debug -b 32 --use-macros */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define S1_I(t,i,j) A[(t + 1) % 2][i][j] = (((0.125 * ((A[t % 2][i + 1 == N ? 0 : i + 1][j] - (2.0 * A[t % 2][i][j])) + A[t % 2][i == 0 ? N - 1 : i - 1][j])) + (0.125 * ((A[t % 2][i][j + 1 == N ? 0 : j + 1] - (2.0 * A[t % 2][i][j])) + A[t % 2][i][j == 0 ? N - 1 : j - 1]))) + A[t % 2][i][j]);
#define S1(t,i,j) S1_I((t),(i),(j))
#pragma scop
for (register int ii0 = 0; ii0 <= 15; ii0 += 1) {
  for (register int k = 0; k <= 2; k += 1) {
    const int ii1_46960_lb = 0, ii1_53605_ub = k + 47;
    #pragma omp parallel for
    for (register int ii1 = ii1_46960_lb; ii1 <= ii1_53605_ub; ii1 += 1) {
      if (ii1 + 1 >= k) {
        if (ii0 == 15 && k == 1 && ii1 == 0) {
          for (register int ii2 = 0; ii2 <= 48; ii2 += 1) {
            const int i0_55369_lb = 480, i0_95047_ub = min(498, 16 * ii2 + 495);
            for (register int i0 = i0_55369_lb; i0 <= i0_95047_ub; i0 += 1) {
              const int i1_44005_lb = max(i0 - 467, -i0 + 499), i1_80745_ub = min(i0 - 436, -i0 + 530);
              for (register int i1 = i1_44005_lb; i1 <= i1_80745_ub; i1 += 1) {
                if (ii2 == 0) {
                  const int i2_99771_lb = i0 - i1 - 436, i2_95206_ub = i0 - 480;
                  for (register int i2 = i2_99771_lb; i2 < i2_95206_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  const int i2_6395_lb = i0 + i1 - 499, i2_56706_ub = i0 - 480;
                  for (register int i2 = i2_6395_lb; i2 < i2_56706_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
                const int i2_45555_lb = max(max(max(i0 - 480, 32 * ii2 - i0 + 480), 32 * ii2 - i0 + i1 + 436), 32 * ii2 - i0 - i1 + 499), i2_8852_ub = 32 * ii2 - i0 + 511;
                for (register int i2 = i2_45555_lb; i2 <= i2_8852_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
                const int i2_5680_lb = 32 * ii2 - i0 + 512, i2_8655_ub = 32 * ii2 - i0 + i1 + 467;
                for (register int i2 = i2_5680_lb; i2 <= i2_8655_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
                const int i2_93714_lb = 32 * ii2 - i0 + 512, i2_28174_ub = 32 * ii2 - i0 - i1 + 530;
                for (register int i2 = i2_93714_lb; i2 <= i2_28174_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
            }
            const int i0_32483_lb = 480, i0_40877_ub = min(498, 16 * ii2 + 495);
            for (register int i0 = i0_32483_lb; i0 <= i0_40877_ub; i0 += 1) {
              const int i1_539_lb = 0, i1_88466_ub = min(i0 - 468, -i0 + 498);
              for (register int i1 = i1_539_lb; i1 <= i1_88466_ub; i1 += 1) {
                if (ii2 == 0) {
                  const int i2_10515_lb = i0 - i1 - 468, i2_8711_ub = i0 - 480;
                  for (register int i2 = i2_10515_lb; i2 < i2_8711_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
                const int i2_31554_lb = max(max(i0 - 480, 32 * ii2 - i0 + 480), 32 * ii2 - i0 + i1 + 468), i2_6837_ub = 32 * ii2 - i0 + 511;
                for (register int i2 = i2_31554_lb; i2 <= i2_6837_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
                const int i2_62207_lb = 32 * ii2 - i0 + 512, i2_25697_ub = 32 * ii2 - i0 + i1 + 499;
                for (register int i2 = i2_62207_lb; i2 <= i2_25697_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
              const int i1_77804_lb = max(i0 + 1101, -i0 + 2067), i1_98137_ub = 1599;
              for (register int i1 = i1_77804_lb; i1 <= i1_98137_ub; i1 += 1) {
                if (ii2 == 0) {
                  const int i2_53435_lb = i0 + i1 - 2067, i2_41116_ub = i0 - 480;
                  for (register int i2 = i2_53435_lb; i2 < i2_41116_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
                const int i2_51742_lb = max(max(i0 - 480, 32 * ii2 - i0 + 480), 32 * ii2 - i0 - i1 + 2067), i2_73235_ub = 32 * ii2 - i0 + 511;
                for (register int i2 = i2_51742_lb; i2 <= i2_73235_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
                const int i2_11360_lb = 32 * ii2 - i0 + 512, i2_23463_ub = 32 * ii2 - i0 - i1 + 2098;
                for (register int i2 = i2_11360_lb; i2 <= i2_23463_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
            }
          }
        } else if (ii0 <= 13 && k == 1 && ii1 == 0) {
          for (register int ii2 = 0; ii2 <= 48; ii2 += 1) {
            const int i0_20561_lb = 32 * ii0, i0_68053_ub = min(32 * ii0 + 30, 32 * ii0 + 8 * ii2 + 22);
            for (register int i0 = i0_20561_lb; i0 <= i0_68053_ub; i0 += 1) {
              if (32 * ii0 + 16 * ii2 >= i0 + 1) {
                const int i1_66923_lb = max(-32 * ii0 + i0 + 1, 32 * ii0 - i0 + 31), i1_26956_ub = 31;
                for (register int i1 = i1_66923_lb; i1 <= i1_26956_ub; i1 += 1) {
                  const int i2_24760_lb = 32 * ii0 + 32 * ii2 - i0 - i1 + 31, i2_12479_ub = 32 * ii0 + 32 * ii2 - i0 - i1 + 62;
                  for (register int i2 = i2_24760_lb; i2 <= i2_12479_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
              } else {
                const int i1_35808_lb = max(-32 * ii0 + i0 + 1, 32 * ii0 - i0 + 31), i1_46792_ub = min(31, 32 * ii0 + 15 * ii2 - i0 + 46);
                for (register int i1 = i1_35808_lb; i1 <= i1_46792_ub; i1 += 1) {
                  const int i2_37486_lb = max(-32 * ii0 + i0 + i1 - 31, 32 * ii0 + 32 * ii2 - i0 - i1 + 31), i2_45875_ub = min(-32 * ii0 + i0 - i1 + 31, 32 * ii0 + 29 * ii2 - i0 - i1 + 62);
                  for (register int i2 = i2_37486_lb; i2 <= i2_45875_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  const int i2_91319_lb = -32 * ii0 + i0 - i1 + 32, i2_69969_ub = 32 * ii0 + 32 * ii2 - i0 + i1;
                  for (register int i2 = i2_91319_lb; i2 < i2_69969_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  const int i2_86752_lb = max(32 * ii0 + 32 * ii2 - i0 + i1, -32 * ii0 + i0 - i1 + 32), i2_8210_ub = 32 * ii0 + 32 * ii2 - i0 - i1 + 62;
                  for (register int i2 = i2_86752_lb; i2 <= i2_8210_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
              }
              const int i1_58435_lb = max(32, -32 * ii0 - 16 * ii2 + i0 + 17), i1_97267_ub = min(-32 * ii0 + i0 + 32, 32 * ii0 - i0 + 62);
              for (register int i1 = i1_58435_lb; i1 <= i1_97267_ub; i1 += 1) {
                const int i2_16921_lb = max(32 * ii0 + 32 * ii2 - i0 + i1 - 32, -32 * ii0 + i0 - i1 + 32), i2_6341_ub = 32 * ii0 + 32 * ii2 - i0 + i1;
                for (register int i2 = i2_16921_lb; i2 < i2_6341_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
            }
            const int i0_20456_lb = 32 * ii0, i0_79129_ub = min(32 * ii0 + 30, 32 * ii0 + 8 * ii2 + 22);
            for (register int i0 = i0_20456_lb; i0 <= i0_79129_ub; i0 += 1) {
              const int i1_48390_lb = max(0, -32 * ii0 - 16 * ii2 + i0 - 15), i1_14612_ub = min(-32 * ii0 + i0, 32 * ii0 - i0 + 30);
              for (register int i1 = i1_48390_lb; i1 <= i1_14612_ub; i1 += 1) {
                const int i2_77266_lb = max(-32 * ii0 + i0 - i1, 32 * ii0 + 32 * ii2 - i0 + i1), i2_1825_ub = 32 * ii0 + 32 * ii2 - i0 + i1 + 31;
                for (register int i2 = i2_77266_lb; i2 <= i2_1825_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
              const int i1_55728_lb = max(-32 * ii0 + i0 + 1569, 32 * ii0 - i0 + 1599), i1_45360_ub = min(1599, 32 * ii0 + 16 * ii2 - i0 + 1614);
              for (register int i1 = i1_55728_lb; i1 <= i1_45360_ub; i1 += 1) {
                const int i2_75060_lb = max(-32 * ii0 + i0 + i1 - 1599, 32 * ii0 + 32 * ii2 - i0 - i1 + 1599), i2_67088_ub = 32 * ii0 + 32 * ii2 - i0 - i1 + 1630;
                for (register int i2 = i2_75060_lb; i2 <= i2_67088_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
            }
          }
        } else if (ii0 == 14 && ii1 + 1 == k) {
          const int ii2_68823_lb = max(3 * k - 6, -2 * k + 2), ii2_59694_ub = 48;
          for (register int ii2 = ii2_68823_lb; ii2 <= ii2_59694_ub; ii2 += 1) {
            if (k == 2) {
              for (register int i0 = 464; i0 <= 479; i0 += 1) {
                const int i1_27747_lb = -i0 + 511, i1_22081_ub = 47;
                for (register int i1 = i1_27747_lb; i1 <= i1_22081_ub; i1 += 1) {
                  const int i2_16341_lb = max(0, 32 * ii2 - i0 - i1 + 511), i2_68859_ub = 32 * ii2;
                  for (register int i2 = i2_16341_lb; i2 < i2_68859_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  if (i1 + 435 >= i0) {
                    const int i2_34560_lb = max(32 * ii2, i0 + i1 - 511), i2_68501_ub = 32 * ii2 - i0 - i1 + 542;
                    for (register int i2 = i2_34560_lb; i2 <= i2_68501_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  } else {
                    const int i2_15652_lb = max(32 * ii2, i0 + i1 - 511), i2_88398_ub = 32 * ii2 - i0 - i1 + 542;
                    for (register int i2 = i2_15652_lb; i2 <= i2_88398_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                }
                const int i1_14376_lb = 48, i1_6971_ub = i0 - 415;
                for (register int i1 = i1_14376_lb; i1 < i1_6971_ub; i1 += 1) {
                  const int i2_58367_lb = max(0, 32 * ii2 - i0 + i1 + 416), i2_17481_ub = 32 * ii2;
                  for (register int i2 = i2_58367_lb; i2 < i2_17481_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  const int i2_15181_lb = max(32 * ii2, i0 - i1 - 416), i2_33154_ub = 32 * ii2 - i0 + i1 + 447;
                  for (register int i2 = i2_15181_lb; i2 <= i2_33154_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
              }
            } else {
              const int i0_14748_lb = 448, i0_48455_ub = min(478, 8 * ii2 + 470);
              for (register int i0 = i0_14748_lb; i0 <= i0_48455_ub; i0 += 1) {
                if (i0 >= 16 * ii2 + 448) {
                  const int i1_39495_lb = max(i0 - 447, -i0 + 479), i1_51557_ub = min(min(31, -i0 + 498), 15 * ii2 - i0 + 494);
                  for (register int i1 = i1_39495_lb; i1 <= i1_51557_ub; i1 += 1) {
                    const int i2_27584_lb = max(i0 + i1 - 479, 32 * ii2 - i0 - i1 + 479), i2_87885_ub = min(i0 - i1 - 417, 29 * ii2 - i0 - i1 + 510);
                    for (register int i2 = i2_27584_lb; i2 <= i2_87885_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                    const int i2_66169_lb = i0 - i1 - 416, i2_4850_ub = 32 * ii2 - i0 - i1 + 510;
                    for (register int i2 = i2_66169_lb; i2 <= i2_4850_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                  const int i1_6062_lb = max(32, -16 * ii2 + i0 - 431), i1_38250_ub = min(i0 - 416, -i0 + 498);
                  for (register int i1 = i1_6062_lb; i1 <= i1_38250_ub; i1 += 1) {
                    const int i2_50210_lb = max(i0 - i1 - 416, 32 * ii2 - i0 + i1 + 416), i2_81122_ub = 32 * ii2 - i0 + i1 + 447;
                    for (register int i2 = i2_50210_lb; i2 <= i2_81122_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                  if (ii2 == 1) {
                    const int i1_5338_lb = max(i0 - 447, -i0 + 499), i1_35385_ub = 31;
                    for (register int i1 = i1_5338_lb; i1 <= i1_35385_ub; i1 += 1) {
                      const int i2_40816_lb = i0 + i1 - 479, i2_60496_ub = -i0 - i1 + 542;
                      for (register int i2 = i2_40816_lb; i2 <= i2_60496_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                  }
                  const int i1_24770_lb = max(max(32, -16 * ii2 + i0 - 431), -i0 + 499), i1_84915_ub = min(i0 - 416, -i0 + 510);
                  for (register int i1 = i1_24770_lb; i1 <= i1_84915_ub; i1 += 1) {
                    const int i2_82577_lb = max(i0 - i1 - 416, 32 * ii2 - i0 + i1 + 416), i2_57463_ub = 32 * ii2 - i0 + i1 + 447;
                    for (register int i2 = i2_82577_lb; i2 <= i2_57463_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                } else {
                  const int i1_70127_lb = max(i0 - 447, -i0 + 479), i1_33489_ub = 31;
                  for (register int i1 = i1_70127_lb; i1 <= i1_33489_ub; i1 += 1) {
                    const int i2_25964_lb = 32 * ii2 - i0 - i1 + 479, i2_2131_ub = 32 * ii2 - i0 - i1 + 510;
                    for (register int i2 = i2_25964_lb; i2 <= i2_2131_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                  const int i1_21887_lb = 32, i1_56693_ub = min(i0 - 416, -i0 + 510);
                  for (register int i1 = i1_21887_lb; i1 <= i1_56693_ub; i1 += 1) {
                    const int i2_9102_lb = 32 * ii2 - i0 + i1 + 416, i2_96606_ub = 32 * ii2 - i0 + i1 + 447;
                    for (register int i2 = i2_9102_lb; i2 <= i2_96606_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                }
              }
            }
            if (k == 1) {
              const int i0_74174_lb = 448, i0_40635_ub = min(472, 8 * ii2 + 470);
              for (register int i0 = i0_74174_lb; i0 <= i0_40635_ub; i0 += 1) {
                if (16 * ii2 + 463 >= i0) {
                  const int i1_29760_lb = 0, i1_5274_ub = i0 - 467;
                  for (register int i1 = i1_29760_lb; i1 < i1_5274_ub; i1 += 1) {
                    if (ii2 == 1) {
                      const int i2_89090_lb = i0 - i1 - 448, i2_85607_ub = i0 - 448;
                      for (register int i2 = i2_89090_lb; i2 < i2_85607_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                    const int i2_73183_lb = max(i0 - 448, 32 * ii2 - i0 + i1 + 448), i2_33026_ub = 32 * ii2 - i0 + i1 + 467;
                    for (register int i2 = i2_73183_lb; i2 <= i2_33026_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                    const int i2_89844_lb = 32 * ii2 - i0 + i1 + 468, i2_55705_ub = 32 * ii2 - i0 + i1 + 479;
                    for (register int i2 = i2_89844_lb; i2 <= i2_55705_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                }
                const int i1_54228_lb = max(max(0, i0 - 467), -16 * ii2 + i0 - 463), i1_12258_ub = min(i0 - 448, -i0 + 478);
                for (register int i1 = i1_54228_lb; i1 <= i1_12258_ub; i1 += 1) {
                  const int i2_93955_lb = max(i0 - i1 - 448, 32 * ii2 - i0 + i1 + 448), i2_20790_ub = 32 * ii2 - i0 + i1 + 479;
                  for (register int i2 = i2_93955_lb; i2 <= i2_20790_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
                const int i1_93380_lb = max(i0 + 1121, -i0 + 2047), i1_99293_ub = min(min(1599, 16 * ii2 - i0 + 2062), -i0 + 2066);
                for (register int i1 = i1_93380_lb; i1 <= i1_99293_ub; i1 += 1) {
                  const int i2_56176_lb = max(i0 + i1 - 2047, 32 * ii2 - i0 - i1 + 2047), i2_50548_ub = 32 * ii2 - i0 - i1 + 2078;
                  for (register int i2 = i2_56176_lb; i2 <= i2_50548_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
                if (16 * ii2 + 463 >= i0) {
                  const int i1_76141_lb = -i0 + 2067, i1_80946_ub = 1599;
                  for (register int i1 = i1_76141_lb; i1 <= i1_80946_ub; i1 += 1) {
                    if (ii2 == 1) {
                      const int i2_35463_lb = i0 + i1 - 2047, i2_58718_ub = i0 - 448;
                      for (register int i2 = i2_35463_lb; i2 < i2_58718_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                    const int i2_54761_lb = max(i0 - 448, 32 * ii2 - i0 - i1 + 2047), i2_21942_ub = 32 * ii2 - i0 - i1 + 2066;
                    for (register int i2 = i2_54761_lb; i2 <= i2_21942_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                    const int i2_8559_lb = 32 * ii2 - i0 - i1 + 2067, i2_97077_ub = 32 * ii2 - i0 - i1 + 2078;
                    for (register int i2 = i2_8559_lb; i2 <= i2_97077_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                }
              }
              const int i0_24073_lb = 473, i0_30446_ub = min(478, 16 * ii2 + 463);
              for (register int i0 = i0_24073_lb; i0 <= i0_30446_ub; i0 += 1) {
                const int i1_53770_lb = 0, i1_49527_ub = -i0 + 478;
                for (register int i1 = i1_53770_lb; i1 <= i1_49527_ub; i1 += 1) {
                  if (ii2 == 1) {
                    const int i2_27052_lb = i0 - i1 - 448, i2_44296_ub = i0 - 448;
                    for (register int i2 = i2_27052_lb; i2 < i2_44296_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                  const int i2_90163_lb = max(i0 - 448, 32 * ii2 - i0 + i1 + 448), i2_73164_ub = 32 * ii2 - i0 + i1 + 467;
                  for (register int i2 = i2_90163_lb; i2 <= i2_73164_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  const int i2_65923_lb = max(i0 - 448, 32 * ii2 - i0 + i1 + 468), i2_95605_ub = 32 * ii2 - i0 + i1 + 479;
                  for (register int i2 = i2_65923_lb; i2 <= i2_95605_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
                const int i1_75123_lb = i0 + 1121, i1_39106_ub = 1599;
                for (register int i1 = i1_75123_lb; i1 <= i1_39106_ub; i1 += 1) {
                  if (ii2 == 1) {
                    const int i2_44984_lb = i0 + i1 - 2047, i2_81319_ub = i0 - 448;
                    for (register int i2 = i2_44984_lb; i2 < i2_81319_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                  const int i2_94811_lb = max(i0 - 448, 32 * ii2 - i0 - i1 + 2047), i2_99212_ub = 32 * ii2 - i0 - i1 + 2066;
                  for (register int i2 = i2_94811_lb; i2 <= i2_99212_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  const int i2_93577_lb = max(i0 - 448, 32 * ii2 - i0 - i1 + 2067), i2_88766_ub = 32 * ii2 - i0 - i1 + 2078;
                  for (register int i2 = i2_93577_lb; i2 <= i2_88766_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
              }
            }
          }
        } else if (k <= 1 && ii1 >= k) {
          for (register int ii2 = 0; ii2 <= 48; ii2 += 1) {
            if (k == 0) {
              const int i0_4412_lb = 32 * ii0, i0_76179_ub = min(482, 32 * ii0 + 14);
              for (register int i0 = i0_4412_lb; i0 <= i0_76179_ub; i0 += 1) {
                const int i1_53857_lb = max(32 * ii1 + i0 - 435, -32 * ii0 + 32 * ii1 + i0 + 33), i1_80553_ub = min(32 * ii0 + 32 * ii1 - i0 + 62, 32 * ii1 - i0 + 530);
                for (register int i1 = i1_53857_lb; i1 <= i1_80553_ub; i1 += 1) {
                  const int i2_73477_lb = max(0, 32 * ii0 + 32 * ii2 - i0), i2_5673_ub = 32 * ii2;
                  for (register int i2 = i2_73477_lb; i2 < i2_5673_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  const int i2_55624_lb = max(32 * ii2, -32 * ii0 + i0), i2_44590_ub = 32 * ii0 + 32 * ii2 - i0 + 31;
                  for (register int i2 = i2_55624_lb; i2 <= i2_44590_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
              }
              if (ii1 == 0) {
                const int i0_27615_lb = 32 * ii0, i0_64183_ub = min(482, 32 * ii0 + 14);
                for (register int i0 = i0_27615_lb; i0 <= i0_64183_ub; i0 += 1) {
                  const int i1_41667_lb = max(i0 - 467, -32 * ii0 + i0 + 1), i1_51689_ub = min(32 * ii0 - i0 + 30, -i0 + 498);
                  for (register int i1 = i1_41667_lb; i1 <= i1_51689_ub; i1 += 1) {
                    const int i2_10982_lb = max(0, 32 * ii0 + 32 * ii2 - i0), i2_95438_ub = 32 * ii2;
                    for (register int i2 = i2_10982_lb; i2 < i2_95438_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                    const int i2_1216_lb = max(32 * ii2, -32 * ii0 + i0), i2_54386_ub = 32 * ii0 + 32 * ii2 - i0 + 31;
                    for (register int i2 = i2_1216_lb; i2 <= i2_54386_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                }
              }
            } else if (ii0 <= 14) {
              const int i0_56086_lb = 32 * ii0, i0_7731_ub = min(32 * ii0 + 30, 32 * ii0 + 8 * ii2 + 22);
              for (register int i0 = i0_56086_lb; i0 <= i0_7731_ub; i0 += 1) {
                {
                  const int i1_27551_lb = max(-32 * ii0 + 32 * ii1 + i0 + 1, 32 * ii0 + 32 * ii1 - i0 + 31), i1_22009_ub = min(min(min(-32 * ii0 + 32 * ii1 + i0 + 32, 32 * ii0 + 32 * ii1 + 16 * ii2 - i0 + 46), 32 * ii0 + 32 * ii1 - i0 + 62), 32 * ii0 - i0 + 1566);
                  for (register int i1 = i1_27551_lb; i1 <= i1_22009_ub; i1 += 1) {
                    const int i2_19689_lb = max(32 * ii0 - 32 * ii1 + 32 * ii2 - i0 + i1 - 32, -32 * ii0 + 32 * ii1 + i0 - i1 + 32), i2_19026_ub = -32 * ii0 - 32 * ii1 + i0 + i1 - 31;
                    for (register int i2 = i2_19689_lb; i2 < i2_19026_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                    const int i2_61116_lb = max(max(32 * ii0 - 32 * ii1 + 32 * ii2 - i0 + i1 - 32, -32 * ii0 - 32 * ii1 + i0 + i1 - 31), 32 * ii0 + 32 * ii1 + 32 * ii2 - i0 - i1 + 31), i2_81025_ub = 32 * ii0 + 32 * ii1 + 32 * ii2 - i0 - i1 + 62;
                    for (register int i2 = i2_61116_lb; i2 <= i2_81025_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                    const int i2_346_lb = 32 * ii0 + 32 * ii1 + 32 * ii2 - i0 - i1 + 63, i2_72279_ub = 32 * ii0 - 32 * ii1 + 32 * ii2 - i0 + i1;
                    for (register int i2 = i2_346_lb; i2 < i2_72279_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                  if (ii2 == 0 && i0 >= 32 * ii0 + 16) {
                    const int i1_80237_lb = -32 * ii0 + 32 * ii1 + i0 + 17, i1_10275_ub = min(32 * ii0 + 32 * ii1 - i0 + 62, 32 * ii0 - i0 + 1566);
                    for (register int i1 = i1_80237_lb; i1 <= i1_10275_ub; i1 += 1) {
                      const int i2_61046_lb = -32 * ii0 + 32 * ii1 + i0 - i1 + 32, i2_16592_ub = 32 * ii0 - 32 * ii1 - i0 + i1;
                      for (register int i2 = i2_61046_lb; i2 < i2_16592_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                  } else if (ii2 == 0 && 32 * ii0 + 15 >= i0) {
                    const int i1_13585_lb = 32 * ii0 + 32 * ii1 - i0 + 47, i1_81810_ub = min(-32 * ii0 + 32 * ii1 + i0 + 32, 32 * ii0 - i0 + 1566);
                    for (register int i1 = i1_13585_lb; i1 <= i1_81810_ub; i1 += 1) {
                      const int i2_9123_lb = -32 * ii0 + 32 * ii1 + i0 - i1 + 32, i2_67442_ub = 32 * ii0 - 32 * ii1 - i0 + i1;
                      for (register int i2 = i2_9123_lb; i2 < i2_67442_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                  }
                }
                if (ii1 == 48) {
                  if (ii2 == 0) {
                    const int i1_62363_lb = max(-32 * ii0 + i0 + 1537, 32 * ii0 - i0 + 1567), i1_82600_ub = min(-32 * ii0 + i0 + 1552, 32 * ii0 - i0 + 1582);
                    for (register int i1 = i1_62363_lb; i1 <= i1_82600_ub; i1 += 1) {
                      const int i2_89467_lb = -32 * ii0 + i0 + i1 - 1567, i2_34339_ub = 32 * ii0 - i0 - i1 + 1598;
                      for (register int i2 = i2_89467_lb; i2 <= i2_34339_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                  }
                  const int i1_27190_lb = max(max(-32 * ii0 + i0 + 1537, -32 * ii0 - 16 * ii2 + i0 + 1553), 32 * ii0 - i0 + 1567), i1_33435_ub = min(-32 * ii0 + i0 + 1568, 32 * ii0 - i0 + 1598);
                  for (register int i1 = i1_27190_lb; i1 <= i1_33435_ub; i1 += 1) {
                    const int i2_14875_lb = max(-32 * ii0 + i0 + i1 - 1567, 32 * ii0 + 32 * ii2 - i0 - i1 + 1567), i2_85210_ub = -32 * ii0 + i0 - i1 + 1567;
                    for (register int i2 = i2_14875_lb; i2 <= i2_85210_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                    const int i2_1476_lb = max(max(32 * ii0 + 32 * ii2 - i0 + i1 - 1568, 32 * ii0 + 32 * ii2 - i0 - i1 + 1567), -32 * ii0 + i0 - i1 + 1568), i2_25857_ub = 32 * ii0 + 32 * ii2 - i0 + i1 - 1536;
                    for (register int i2 = i2_1476_lb; i2 < i2_25857_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                    const int i2_97000_lb = 32 * ii0 + 32 * ii2 - i0 + i1 - 1536, i2_2692_ub = 32 * ii0 + 32 * ii2 - i0 - i1 + 1598;
                    for (register int i2 = i2_97000_lb; i2 <= i2_2692_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                }
              }
            } else {
              const int i0_80243_lb = 480, i0_53086_ub = min(498, 11 * ii2 - (ii2 + 1) / 3 + 492);
              for (register int i0 = i0_80243_lb; i0 <= i0_53086_ub; i0 += 1) {
                const int i1_26776_lb = max(32 * ii1 + i0 - 467, 32 * ii1 - i0 + 499), i1_24146_ub = 32 * ii1 + 19;
                for (register int i1 = i1_26776_lb; i1 <= i1_24146_ub; i1 += 1) {
                  const int i2_75096_lb = max(i0 - 480, 32 * ii1 + 32 * ii2 - i0 - i1 + 499), i2_62817_ub = 32 * ii2;
                  for (register int i2 = i2_75096_lb; i2 < i2_62817_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  const int i2_43173_lb = max(32 * ii2, -32 * ii1 + i0 + i1 - 499), i2_52564_ub = 32 * ii1 + 32 * ii2 - i0 - i1 + 530;
                  for (register int i2 = i2_43173_lb; i2 <= i2_52564_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
                const int i1_43842_lb = max(32 * ii1 + 20, 32 * ii1 + i0 - 467), i1_59871_ub = min(min(32 * ii1 + 44, 32 * ii1 - i0 + 530), 32 * ii1 + 32 * ii2 - 2 * i0 + 1010);
                for (register int i1 = i1_43842_lb; i1 <= i1_59871_ub; i1 += 1) {
                  const int i2_24843_lb = max(i0 - 480, 32 * ii2 - i0 + 480), i2_40431_ub = 32 * ii2;
                  for (register int i2 = i2_24843_lb; i2 < i2_40431_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  const int i2_70146_lb = max(32 * ii2, i0 - 480), i2_2241_ub = 32 * ii2 - i0 + 511;
                  for (register int i2 = i2_70146_lb; i2 <= i2_2241_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
                if (ii2 == 0) {
                  const int i1_57024_lb = 32 * ii1 - 2 * i0 + 1011, i1_83731_ub = min(32 * ii1 + 44, 32 * ii1 - i0 + 530);
                  for (register int i1 = i1_57024_lb; i1 <= i1_83731_ub; i1 += 1) {
                    const int i2_84051_lb = i0 - 480, i2_66147_ub = -i0 + 511;
                    for (register int i2 = i2_84051_lb; i2 <= i2_66147_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                }
                const int i1_67526_lb = 32 * ii1 + 45, i1_62767_ub = min(32 * ii1 + i0 - 436, 32 * ii1 - i0 + 530);
                for (register int i1 = i1_67526_lb; i1 <= i1_62767_ub; i1 += 1) {
                  const int i2_65100_lb = max(i0 - 480, -32 * ii1 + 32 * ii2 - i0 + i1 + 436), i2_73345_ub = 32 * ii2;
                  for (register int i2 = i2_65100_lb; i2 < i2_73345_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  const int i2_97106_lb = max(32 * ii2, 32 * ii1 + i0 - i1 - 436), i2_8642_ub = 32 * ii2 - i0 + 511;
                  for (register int i2 = i2_97106_lb; i2 <= i2_8642_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  const int i2_6780_lb = 32 * ii2 - i0 + 512, i2_11981_ub = -32 * ii1 + 32 * ii2 - i0 + i1 + 467;
                  for (register int i2 = i2_6780_lb; i2 <= i2_11981_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
              }
              if (ii2 == 0) {
                for (register int i0 = 493; i0 <= 495; i0 += 1) {
                  const int i1_54190_lb = 32 * ii1 + i0 - 467, i1_90852_ub = 32 * ii1 - i0 + 530;
                  for (register int i1 = i1_54190_lb; i1 <= i1_90852_ub; i1 += 1) {
                    const int i2_27301_lb = i0 - 480, i2_50786_ub = -i0 + 511;
                    for (register int i2 = i2_27301_lb; i2 <= i2_50786_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                }
              }
            }
            if (k == 0 && ii1 == 0) {
              const int i0_43939_lb = 32 * ii0, i0_70429_ub = min(482, 32 * ii0 + 14);
              for (register int i0 = i0_43939_lb; i0 <= i0_70429_ub; i0 += 1) {
                const int i1_74932_lb = max(i0 + 1101, -32 * ii0 + i0 + 1569), i1_35387_ub = min(32 * ii0 - i0 + 1598, -i0 + 2066);
                for (register int i1 = i1_74932_lb; i1 <= i1_35387_ub; i1 += 1) {
                  const int i2_33246_lb = max(0, 32 * ii0 + 32 * ii2 - i0), i2_34457_ub = 32 * ii2;
                  for (register int i2 = i2_33246_lb; i2 < i2_34457_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  const int i2_4303_lb = max(32 * ii2, -32 * ii0 + i0), i2_93440_ub = 32 * ii0 + 32 * ii2 - i0 + 31;
                  for (register int i2 = i2_4303_lb; i2 <= i2_93440_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
              }
            }
          }
        }
        if (k <= 1) {
          if (ii0 <= 14 && k == 1 && ii1 == 0) {
            const int i0_94328_lb = 32 * ii0, i0_29146_ub = 32 * ii0 + 30;
            for (register int i0 = i0_94328_lb; i0 <= i0_29146_ub; i0 += 1) {
              const int i1_33871_lb = max(-32 * ii0 + i0 + 1, 32 * ii0 - i0 + 31), i1_80827_ub = min(-32 * ii0 + i0 + 32, 32 * ii0 - i0 + 62);
              for (register int i1 = i1_33871_lb; i1 <= i1_80827_ub; i1 += 1) {
                const int i2_47740_lb = 0, i2_7247_ub = min(-32 * ii0 + i0 + i1 - 32, -32 * ii0 + i0 - i1 + 31);
                for (register int i2 = i2_47740_lb; i2 <= i2_7247_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
                const int i2_80910_lb = max(32 * ii0 - i0 + i1 + 1536, 32 * ii0 - i0 - i1 + 1599), i2_31791_ub = 1599;
                for (register int i2 = i2_80910_lb; i2 <= i2_31791_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
            }
          } else if (ii0 == 15 && k == 1 && ii1 == 0) {
            for (register int i0 = 480; i0 <= 498; i0 += 1) {
              const int i1_94558_lb = max(i0 - 467, -i0 + 499), i1_54847_ub = min(i0 - 436, -i0 + 530);
              for (register int i1 = i1_94558_lb; i1 <= i1_54847_ub; i1 += 1) {
                const int i2_38134_lb = 0, i2_8017_ub = min(min(i0 - 480, i0 + i1 - 499), i0 - i1 - 436);
                for (register int i2 = i2_38134_lb; i2 < i2_8017_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
                const int i2_63489_lb = max(max(-i0 + 2048, -i0 + i1 + 2004), -i0 - i1 + 2067), i2_44914_ub = 1599;
                for (register int i2 = i2_63489_lb; i2 <= i2_44914_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
            }
          } else if (k == 0) {
            const int i0_19998_lb = 32 * ii0, i0_73694_ub = min(482, 32 * ii0 + 14);
            for (register int i0 = i0_19998_lb; i0 <= i0_73694_ub; i0 += 1) {
              const int i1_69523_lb = max(32 * ii1 + i0 - 435, -32 * ii0 + 32 * ii1 + i0 + 33), i1_90541_ub = min(32 * ii0 + 32 * ii1 - i0 + 62, 32 * ii1 - i0 + 530);
              for (register int i1 = i1_69523_lb; i1 <= i1_90541_ub; i1 += 1) {
                const int i2_64546_lb = 0, i2_96824_ub = -32 * ii0 + i0;
                for (register int i2 = i2_64546_lb; i2 < i2_96824_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
                const int i2_41327_lb = 32 * ii0 - i0 + 1568, i2_24837_ub = 1599;
                for (register int i2 = i2_41327_lb; i2 <= i2_24837_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
            }
            if (ii1 == 0) {
              const int i0_67253_lb = 32 * ii0, i0_16259_ub = min(482, 32 * ii0 + 14);
              for (register int i0 = i0_67253_lb; i0 <= i0_16259_ub; i0 += 1) {
                const int i1_76576_lb = max(i0 - 467, -32 * ii0 + i0 + 1), i1_16851_ub = min(32 * ii0 - i0 + 30, -i0 + 498);
                for (register int i1 = i1_76576_lb; i1 <= i1_16851_ub; i1 += 1) {
                  const int i2_67069_lb = 0, i2_80879_ub = -32 * ii0 + i0;
                  for (register int i2 = i2_67069_lb; i2 < i2_80879_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  const int i2_10291_lb = 32 * ii0 - i0 + 1568, i2_77749_ub = 1599;
                  for (register int i2 = i2_10291_lb; i2 <= i2_77749_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
              }
            }
          } else if (ii0 <= 14) {
            const int i0_26378_lb = 32 * ii0, i0_60514_ub = 32 * ii0 + 30;
            for (register int i0 = i0_26378_lb; i0 <= i0_60514_ub; i0 += 1) {
              const int i1_58576_lb = max(-32 * ii0 + 32 * ii1 + i0 + 1, 32 * ii0 + 32 * ii1 - i0 + 31), i1_74118_ub = min(-32 * ii0 + 32 * ii1 + i0 + 32, 32 * ii0 + 32 * ii1 - i0 + 62);
              for (register int i1 = i1_58576_lb; i1 <= i1_74118_ub; i1 += 1) {
                const int i2_67762_lb = 0, i2_39487_ub = min(-32 * ii0 + i0 + i1 - 1568, -32 * ii0 + i0 - i1 + 1567);
                for (register int i2 = i2_67762_lb; i2 <= i2_39487_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
                if (ii1 <= 47) {
                  const int i2_22261_lb = 0, i2_57509_ub = min(-32 * ii0 - 32 * ii1 + i0 + i1 - 32, -32 * ii0 + 32 * ii1 + i0 - i1 + 31);
                  for (register int i2 = i2_22261_lb; i2 <= i2_57509_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
                const int i2_4275_lb = max(32 * ii0 - 32 * ii1 - i0 + i1 + 1536, 32 * ii0 + 32 * ii1 - i0 - i1 + 1599), i2_16820_ub = 1599;
                for (register int i2 = i2_4275_lb; i2 <= i2_16820_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
            }
          } else {
            for (register int i0 = 480; i0 <= 498; i0 += 1) {
              const int i1_24837_lb = max(32 * ii1 + i0 - 467, 32 * ii1 - i0 + 499), i1_8549_ub = min(32 * ii1 + i0 - 436, 32 * ii1 - i0 + 530);
              for (register int i1 = i1_24837_lb; i1 <= i1_8549_ub; i1 += 1) {
                const int i2_87324_lb = 0, i2_61187_ub = min(min(i0 - 480, -32 * ii1 + i0 + i1 - 499), 32 * ii1 + i0 - i1 - 436);
                for (register int i2 = i2_87324_lb; i2 < i2_61187_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
                const int i2_82243_lb = max(max(-i0 + 2048, -32 * ii1 - i0 + i1 + 2004), 32 * ii1 - i0 - i1 + 2067), i2_56847_ub = 1599;
                for (register int i2 = i2_82243_lb; i2 <= i2_56847_ub; i2 += 1) {
                  S1(i0, i1, i2);
                }
              }
            }
          }
          if (ii1 == 0) {
            if (k == 1) {
              if (ii0 <= 14) {
                const int i0_68080_lb = 32 * ii0, i0_63142_ub = min(32 * ii0 + 30, 16 * ii0 + 248);
                for (register int i0 = i0_68080_lb; i0 <= i0_63142_ub; i0 += 1) {
                  if (ii0 == 14) {
                    const int i1_70023_lb = 0, i1_9407_ub = i0 - 467;
                    for (register int i1 = i1_70023_lb; i1 < i1_9407_ub; i1 += 1) {
                      const int i2_4331_lb = 0, i2_53628_ub = i0 - i1 - 448;
                      for (register int i2 = i2_4331_lb; i2 < i2_53628_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                      const int i2_42019_lb = -i0 + i1 + 2016, i2_80908_ub = -i0 + i1 + 2035;
                      for (register int i2 = i2_42019_lb; i2 <= i2_80908_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                      const int i2_70479_lb = -i0 + i1 + 2036, i2_9088_ub = 1599;
                      for (register int i2 = i2_70479_lb; i2 <= i2_9088_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                  }
                  const int i1_78139_lb = max(0, i0 - 467), i1_97122_ub = min(-32 * ii0 + i0, 32 * ii0 - i0 + 30);
                  for (register int i1 = i1_78139_lb; i1 <= i1_97122_ub; i1 += 1) {
                    const int i2_86837_lb = 0, i2_4517_ub = -32 * ii0 + i0 - i1;
                    for (register int i2 = i2_86837_lb; i2 < i2_4517_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                    const int i2_57636_lb = 32 * ii0 - i0 + i1 + 1568, i2_61766_ub = 1599;
                    for (register int i2 = i2_57636_lb; i2 <= i2_61766_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                  const int i1_94987_lb = max(-32 * ii0 + i0 + 1569, 32 * ii0 - i0 + 1599), i1_41750_ub = min(1599, -i0 + 2066);
                  for (register int i1 = i1_94987_lb; i1 <= i1_41750_ub; i1 += 1) {
                    const int i2_1253_lb = 0, i2_17249_ub = -32 * ii0 + i0 + i1 - 1599;
                    for (register int i2 = i2_1253_lb; i2 < i2_17249_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                    const int i2_15611_lb = 32 * ii0 - i0 - i1 + 3167, i2_5528_ub = 1599;
                    for (register int i2 = i2_15611_lb; i2 <= i2_5528_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                  if (ii0 == 14) {
                    const int i1_34069_lb = -i0 + 2067, i1_44319_ub = 1599;
                    for (register int i1 = i1_34069_lb; i1 <= i1_44319_ub; i1 += 1) {
                      const int i2_64290_lb = 0, i2_75258_ub = i0 + i1 - 2047;
                      for (register int i2 = i2_64290_lb; i2 < i2_75258_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                      const int i2_52869_lb = -i0 - i1 + 3615, i2_67966_ub = -i0 - i1 + 3634;
                      for (register int i2 = i2_52869_lb; i2 <= i2_67966_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                      const int i2_52797_lb = -i0 - i1 + 3635, i2_51464_ub = 1599;
                      for (register int i2 = i2_52797_lb; i2 <= i2_51464_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                  }
                }
                if (ii0 == 14) {
                  for (register int i0 = 473; i0 <= 478; i0 += 1) {
                    const int i1_14606_lb = 0, i1_11188_ub = -i0 + 478;
                    for (register int i1 = i1_14606_lb; i1 <= i1_11188_ub; i1 += 1) {
                      const int i2_46637_lb = 0, i2_18938_ub = i0 - i1 - 448;
                      for (register int i2 = i2_46637_lb; i2 < i2_18938_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                      const int i2_64816_lb = -i0 + i1 + 2016, i2_88656_ub = -i0 + i1 + 2035;
                      for (register int i2 = i2_64816_lb; i2 <= i2_88656_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                      const int i2_16198_lb = -i0 + i1 + 2036, i2_51647_ub = 1599;
                      for (register int i2 = i2_16198_lb; i2 <= i2_51647_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                    const int i1_97744_lb = i0 + 1121, i1_94337_ub = 1599;
                    for (register int i1 = i1_97744_lb; i1 <= i1_94337_ub; i1 += 1) {
                      const int i2_48769_lb = 0, i2_934_ub = i0 + i1 - 2047;
                      for (register int i2 = i2_48769_lb; i2 < i2_934_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                      const int i2_15207_lb = -i0 - i1 + 3615, i2_6405_ub = 1567;
                      for (register int i2 = i2_15207_lb; i2 <= i2_6405_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                      for (register int i2 = 1568; i2 <= 1599; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                  }
                }
              } else {
                for (register int i0 = 480; i0 <= 498; i0 += 1) {
                  const int i1_27443_lb = 0, i1_80119_ub = min(i0 - 468, -i0 + 498);
                  for (register int i1 = i1_27443_lb; i1 <= i1_80119_ub; i1 += 1) {
                    const int i2_2185_lb = 0, i2_77864_ub = min(i0 - 480, i0 - i1 - 468);
                    for (register int i2 = i2_2185_lb; i2 < i2_77864_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                    const int i2_24439_lb = max(-i0 + 2048, -i0 + i1 + 2036), i2_66475_ub = 1599;
                    for (register int i2 = i2_24439_lb; i2 <= i2_66475_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                  const int i1_69474_lb = max(i0 + 1101, -i0 + 2067), i1_93660_ub = 1599;
                  for (register int i1 = i1_69474_lb; i1 <= i1_93660_ub; i1 += 1) {
                    const int i2_34441_lb = 0, i2_38624_ub = min(i0 - 480, i0 + i1 - 2067);
                    for (register int i2 = i2_34441_lb; i2 < i2_38624_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                    const int i2_45124_lb = max(-i0 + 2048, -i0 - i1 + 3635), i2_75606_ub = 1599;
                    for (register int i2 = i2_45124_lb; i2 <= i2_75606_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                }
              }
            } else {
              const int i0_59502_lb = 32 * ii0, i0_76083_ub = min(482, 32 * ii0 + 14);
              for (register int i0 = i0_59502_lb; i0 <= i0_76083_ub; i0 += 1) {
                const int i1_86794_lb = max(i0 + 1101, -32 * ii0 + i0 + 1569), i1_6139_ub = min(32 * ii0 - i0 + 1598, -i0 + 2066);
                for (register int i1 = i1_86794_lb; i1 <= i1_6139_ub; i1 += 1) {
                  const int i2_11373_lb = 0, i2_67962_ub = -32 * ii0 + i0;
                  for (register int i2 = i2_11373_lb; i2 < i2_67962_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  const int i2_11148_lb = 32 * ii0 - i0 + 1568, i2_27571_ub = 1599;
                  for (register int i2 = i2_11148_lb; i2 <= i2_27571_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
              }
            }
          }
        } else if (ii0 == 14 && ii1 == 1) {
          for (register int i0 = 464; i0 <= 479; i0 += 1) {
            const int i1_38260_lb = -i0 + 511, i1_68378_ub = i0 - 415;
            for (register int i1 = i1_38260_lb; i1 < i1_68378_ub; i1 += 1) {
              const int i2_26178_lb = 0, i2_53467_ub = min(i0 + i1 - 511, i0 - i1 - 416);
              for (register int i2 = i2_26178_lb; i2 < i2_53467_ub; i2 += 1) {
                S1(i0, i1, i2);
              }
              const int i2_91136_lb = max(-i0 + i1 + 1984, -i0 - i1 + 2079), i2_5230_ub = 1599;
              for (register int i2 = i2_91136_lb; i2 <= i2_5230_ub; i2 += 1) {
                S1(i0, i1, i2);
              }
            }
          }
        }
        if (k == 2 && ii1 >= 2) {
          for (register int ii2 = 0; ii2 <= 49; ii2 += 1) {
            if (ii0 <= 14 && ii1 >= 48) {
              const int i0_85535_lb = 32 * ii0 + 16, i0_7457_ub = 32 * ii0 + 31;
              for (register int i0 = i0_85535_lb; i0 <= i0_7457_ub; i0 += 1) {
                if (ii1 == 49) {
                  if (ii2 == 49 && 16 * ii0 + 249 >= i0) {
                    for (register int i2 = 1568; i2 <= 1599; i2 += 1) {
                      S1(i0, 32 * ii0 - i0 + 1599, i2);
                    }
                  }
                  if (ii0 == 14) {
                    const int i1_1674_lb = -i0 + 2047, i1_92906_ub = i0 + 1100;
                    for (register int i1 = i1_1674_lb; i1 <= i1_92906_ub; i1 += 1) {
                      if (ii2 <= 48) {
                        const int i2_70548_lb = max(0, 32 * ii2 - i0 - i1 + 2047), i2_87500_ub = 32 * ii2;
                        for (register int i2 = i2_70548_lb; i2 < i2_87500_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                        const int i2_2918_lb = max(32 * ii2, i0 + i1 - 2047), i2_4990_ub = 32 * ii2 - i0 - i1 + 2078;
                        for (register int i2 = i2_2918_lb; i2 <= i2_4990_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                      } else {
                        const int i2_26124_lb = 0, i2_48043_ub = i0 + i1 - 2047;
                        for (register int i2 = i2_26124_lb; i2 < i2_48043_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                        const int i2_96948_lb = -i0 - i1 + 3615, i2_1978_ub = 1599;
                        for (register int i2 = i2_96948_lb; i2 <= i2_1978_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                      }
                    }
                  }
                  const int i1_24126_lb = max(max(i0 + 1101, 32 * ii0 + ii2 - i0 + 1551), 32 * ii0 - i0 + 1599), i1_95_ub = min(15 * ii2 + 1583, -32 * ii0 + i0 + 1567);
                  for (register int i1 = i1_24126_lb; i1 <= i1_95_ub; i1 += 1) {
                    if (ii2 <= 48) {
                      const int i2_24470_lb = max(max(0, 32 * ii0 + 32 * ii2 - i0 + i1 - 1568), 32 * ii0 + 32 * ii2 - i0 - i1 + 1599), i2_35499_ub = 32 * ii2;
                      for (register int i2 = i2_24470_lb; i2 < i2_35499_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                      if (ii2 == 0) {
                        const int i2_68057_lb = -32 * ii0 + i0 + i1 - 1599, i2_35618_ub = -32 * ii0 + i0 - i1 + 1567;
                        for (register int i2 = i2_68057_lb; i2 <= i2_35618_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                      }
                      const int i2_79422_lb = max(32 * ii2, -32 * ii0 + i0 - i1 + 1568), i2_4019_ub = 32 * ii0 + 32 * ii2 - i0 + i1 - 1536;
                      for (register int i2 = i2_79422_lb; i2 < i2_4019_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                      const int i2_77214_lb = max(32 * ii0 + 32 * ii2 - i0 + i1 - 1536, -32 * ii0 + i0 - i1 + 1568), i2_17682_ub = 32 * ii0 + 32 * ii2 - i0 - i1 + 1630;
                      for (register int i2 = i2_77214_lb; i2 <= i2_17682_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    } else {
                      const int i2_88749_lb = 0, i2_3393_ub = min(-32 * ii0 + i0 + i1 - 1600, -32 * ii0 + i0 - i1 + 1567);
                      for (register int i2 = i2_88749_lb; i2 <= i2_3393_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                      const int i2_71150_lb = max(32 * ii0 - i0 + i1, 32 * ii0 - i0 - i1 + 3167), i2_96237_ub = 1599;
                      for (register int i2 = i2_71150_lb; i2 <= i2_96237_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                  }
                  if (ii2 == 0) {
                    const int i1_24975_lb = 1584, i1_51164_ub = -32 * ii0 + i0 + 1567;
                    for (register int i1 = i1_24975_lb; i1 <= i1_51164_ub; i1 += 1) {
                      const int i2_68233_lb = -32 * ii0 + i0 - i1 + 1568, i2_26863_ub = 32 * ii0 - i0 + i1 - 1536;
                      for (register int i2 = i2_68233_lb; i2 < i2_26863_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                  }
                  if (ii2 <= 48) {
                    const int i2_74973_lb = 32 * ii2, i2_36701_ub = 32 * ii2 + 31;
                    for (register int i2 = i2_74973_lb; i2 <= i2_36701_ub; i2 += 1) {
                      S1(i0, -32 * ii0 + i0 + 1568, i2);
                    }
                  } else {
                    for (register int i2 = 1568; i2 <= 1599; i2 += 1) {
                      S1(i0, -32 * ii0 + i0 + 1568, i2);
                    }
                  }
                } else {
                  if (ii2 == 0) {
                    const int i1_29607_lb = 32 * ii0 - i0 + 1567, i1_1484_ub = -32 * ii0 + i0 + 1520;
                    for (register int i1 = i1_29607_lb; i1 <= i1_1484_ub; i1 += 1) {
                      const int i2_96852_lb = -32 * ii0 + i0 + i1 - 1567, i2_32526_ub = 32 * ii0 - i0 - i1 + 1598;
                      for (register int i2 = i2_96852_lb; i2 <= i2_32526_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                  }
                  const int i1_22826_lb = max(-32 * ii0 - 16 * ii2 + i0 + 1521, 32 * ii0 - i0 + 1567), i1_22976_ub = -32 * ii0 + i0 + 1536;
                  for (register int i1 = i1_22826_lb; i1 <= i1_22976_ub; i1 += 1) {
                    if (ii2 <= 48) {
                      const int i2_96921_lb = max(max(0, 32 * ii0 + 32 * ii2 - i0 + i1 - 1536), 32 * ii0 + 32 * ii2 - i0 - i1 + 1567), i2_36127_ub = 32 * ii2;
                      for (register int i2 = i2_96921_lb; i2 < i2_36127_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                      if (ii2 == 0) {
                        const int i2_41307_lb = -32 * ii0 + i0 + i1 - 1567, i2_37399_ub = -32 * ii0 + i0 - i1 + 1535;
                        for (register int i2 = i2_41307_lb; i2 <= i2_37399_ub; i2 += 1) {
                          S1(i0, i1, i2);
                        }
                      }
                      const int i2_36222_lb = max(32 * ii2, -32 * ii0 + i0 - i1 + 1536), i2_82129_ub = 32 * ii0 + 32 * ii2 - i0 + i1 - 1504;
                      for (register int i2 = i2_36222_lb; i2 < i2_82129_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                      const int i2_89250_lb = 32 * ii0 + 32 * ii2 - i0 + i1 - 1504, i2_20631_ub = 32 * ii0 + 32 * ii2 - i0 - i1 + 1598;
                      for (register int i2 = i2_89250_lb; i2 <= i2_20631_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    } else {
                      const int i2_17747_lb = 0, i2_68672_ub = min(-32 * ii0 + i0 + i1 - 1568, -32 * ii0 + i0 - i1 + 1535);
                      for (register int i2 = i2_17747_lb; i2 <= i2_68672_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                      const int i2_24650_lb = max(32 * ii0 - i0 + i1 + 32, 32 * ii0 - i0 - i1 + 3135), i2_11313_ub = 1599;
                      for (register int i2 = i2_24650_lb; i2 <= i2_11313_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                  }
                }
              }
            } else if (ii0 <= 14 && ii1 <= 47) {
              const int i0_86354_lb = 32 * ii0 + 16, i0_29752_ub = 32 * ii0 + 31;
              for (register int i0 = i0_86354_lb; i0 <= i0_29752_ub; i0 += 1) {
                if (ii2 <= 48) {
                  const int i1_14706_lb = 32 * ii0 + 32 * ii1 - i0 + 31, i1_73856_ub = 32 * ii1 + 15;
                  for (register int i1 = i1_14706_lb; i1 <= i1_73856_ub; i1 += 1) {
                    const int i2_25989_lb = max(0, 32 * ii0 + 32 * ii1 + 32 * ii2 - i0 - i1 + 31), i2_56034_ub = 32 * ii2;
                    for (register int i2 = i2_25989_lb; i2 < i2_56034_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                    const int i2_41372_lb = max(32 * ii2, -32 * ii0 - 32 * ii1 + i0 + i1 - 31), i2_10575_ub = 32 * ii0 + 32 * ii1 + 32 * ii2 - i0 - i1 + 62;
                    for (register int i2 = i2_41372_lb; i2 <= i2_10575_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                  const int i1_82897_lb = 32 * ii1 + 16, i1_16346_ub = -32 * ii0 + 32 * ii1 + i0;
                  for (register int i1 = i1_82897_lb; i1 <= i1_16346_ub; i1 += 1) {
                    const int i2_47276_lb = max(0, 32 * ii0 - 32 * ii1 + 32 * ii2 - i0 + i1), i2_30185_ub = 32 * ii2;
                    for (register int i2 = i2_47276_lb; i2 < i2_30185_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                    const int i2_25697_lb = max(32 * ii2, -32 * ii0 + 32 * ii1 + i0 - i1), i2_76883_ub = 32 * ii0 - 32 * ii1 + 32 * ii2 - i0 + i1 + 31;
                    for (register int i2 = i2_25697_lb; i2 <= i2_76883_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                } else {
                  const int i1_48021_lb = 32 * ii0 + 32 * ii1 - i0 + 31, i1_22549_ub = -32 * ii0 + 32 * ii1 + i0;
                  for (register int i1 = i1_48021_lb; i1 <= i1_22549_ub; i1 += 1) {
                    const int i2_25761_lb = 0, i2_87200_ub = min(-32 * ii0 - 32 * ii1 + i0 + i1 - 31, -32 * ii0 + 32 * ii1 + i0 - i1);
                    for (register int i2 = i2_25761_lb; i2 < i2_87200_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                    const int i2_61878_lb = max(32 * ii0 - 32 * ii1 - i0 + i1 + 1568, 32 * ii0 + 32 * ii1 - i0 - i1 + 1599), i2_39034_ub = 1599;
                    for (register int i2 = i2_61878_lb; i2 <= i2_39034_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                }
              }
            } else {
              for (register int i0 = 484; i0 <= 499; i0 += 1) {
                const int i1_76433_lb = 32 * ii1 - i0 + 499, i1_59549_ub = min(32 * ii1 + i0 - 468, 32 * ii1 + 16 * ii2 - i0 + 514);
                for (register int i1 = i1_76433_lb; i1 <= i1_59549_ub; i1 += 1) {
                  if (ii2 <= 48) {
                    const int i2_85314_lb = max(max(0, -32 * ii1 + 32 * ii2 - i0 + i1 + 468), 32 * ii1 + 32 * ii2 - i0 - i1 + 499), i2_65683_ub = 32 * ii2;
                    for (register int i2 = i2_85314_lb; i2 < i2_65683_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                    if (ii2 == 0) {
                      const int i2_96532_lb = 32 * ii1 + i0 - i1 - 468, i2_19413_ub = -32 * ii1 + i0 + i1 - 499;
                      for (register int i2 = i2_96532_lb; i2 < i2_19413_ub; i2 += 1) {
                        S1(i0, i1, i2);
                      }
                    }
                    const int i2_34355_lb = max(32 * ii2, -32 * ii1 + i0 + i1 - 499), i2_37535_ub = 32 * ii1 + 32 * ii2 - i0 - i1 + 530;
                    for (register int i2 = i2_34355_lb; i2 <= i2_37535_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                    const int i2_30726_lb = 32 * ii1 + 32 * ii2 - i0 - i1 + 531, i2_37062_ub = -32 * ii1 + 32 * ii2 - i0 + i1 + 499;
                    for (register int i2 = i2_30726_lb; i2 <= i2_37062_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  } else {
                    const int i2_67287_lb = 0, i2_61785_ub = min(-32 * ii1 + i0 + i1 - 499, 32 * ii1 + i0 - i1 - 468);
                    for (register int i2 = i2_67287_lb; i2 < i2_61785_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                    const int i2_27270_lb = max(-32 * ii1 - i0 + i1 + 2036, 32 * ii1 - i0 - i1 + 2067), i2_93276_ub = 1599;
                    for (register int i2 = i2_27270_lb; i2 <= i2_93276_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                }
                if (ii2 == 0) {
                  const int i1_34171_lb = 32 * ii1 - i0 + 515, i1_68643_ub = 32 * ii1 + i0 - 467;
                  for (register int i1 = i1_34171_lb; i1 < i1_68643_ub; i1 += 1) {
                    const int i2_3851_lb = 32 * ii1 + i0 - i1 - 468, i2_33420_ub = -32 * ii1 - i0 + i1 + 499;
                    for (register int i2 = i2_3851_lb; i2 <= i2_33420_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                }
              }
            }
          }
        }
        if (ii0 <= 13 && k == 2 && ii1 == 1) {
          for (register int ii2 = 0; ii2 <= 49; ii2 += 1) {
            const int i0_63605_lb = 32 * ii0 + 16, i0_27038_ub = 32 * ii0 + 31;
            for (register int i0 = i0_63605_lb; i0 <= i0_27038_ub; i0 += 1) {
              if (ii2 <= 48) {
                const int i1_60715_lb = 32 * ii0 - i0 + 63, i1_27978_ub = 47;
                for (register int i1 = i1_60715_lb; i1 <= i1_27978_ub; i1 += 1) {
                  const int i2_65940_lb = max(0, 32 * ii0 + 32 * ii2 - i0 - i1 + 63), i2_86476_ub = 32 * ii2;
                  for (register int i2 = i2_65940_lb; i2 < i2_86476_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  const int i2_15178_lb = max(32 * ii2, -32 * ii0 + i0 + i1 - 63), i2_27818_ub = 32 * ii0 + 32 * ii2 - i0 - i1 + 94;
                  for (register int i2 = i2_15178_lb; i2 <= i2_27818_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
                const int i1_25511_lb = 48, i1_38505_ub = -32 * ii0 + i0 + 32;
                for (register int i1 = i1_25511_lb; i1 <= i1_38505_ub; i1 += 1) {
                  const int i2_47355_lb = max(0, 32 * ii0 + 32 * ii2 - i0 + i1 - 32), i2_18296_ub = 32 * ii2;
                  for (register int i2 = i2_47355_lb; i2 < i2_18296_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  const int i2_14406_lb = max(32 * ii2, -32 * ii0 + i0 - i1 + 32), i2_49021_ub = 32 * ii0 + 32 * ii2 - i0 + i1;
                  for (register int i2 = i2_14406_lb; i2 < i2_49021_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
              } else {
                const int i1_83980_lb = 32 * ii0 - i0 + 63, i1_27291_ub = -32 * ii0 + i0 + 32;
                for (register int i1 = i1_83980_lb; i1 <= i1_27291_ub; i1 += 1) {
                  const int i2_68434_lb = 0, i2_34687_ub = min(-32 * ii0 + i0 + i1 - 64, -32 * ii0 + i0 - i1 + 31);
                  for (register int i2 = i2_68434_lb; i2 <= i2_34687_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  const int i2_64826_lb = max(32 * ii0 - i0 + i1 + 1536, 32 * ii0 - i0 - i1 + 1631), i2_15512_ub = 1599;
                  for (register int i2 = i2_64826_lb; i2 <= i2_15512_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
              }
            }
          }
        } else if (ii0 == 15 && k == 2 && ii1 == 1) {
          for (register int ii2 = 0; ii2 <= 49; ii2 += 1) {
            for (register int i0 = 484; i0 <= 499; i0 += 1) {
              const int i1_41741_lb = -i0 + 531, i1_27820_ub = min(i0 - 436, 16 * ii2 - i0 + 546);
              for (register int i1 = i1_41741_lb; i1 <= i1_27820_ub; i1 += 1) {
                if (ii2 <= 48) {
                  const int i2_84015_lb = max(max(0, 32 * ii2 - i0 + i1 + 436), 32 * ii2 - i0 - i1 + 531), i2_61945_ub = 32 * ii2;
                  for (register int i2 = i2_84015_lb; i2 < i2_61945_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  if (ii2 == 0) {
                    const int i2_61240_lb = i0 - i1 - 436, i2_85356_ub = i0 + i1 - 531;
                    for (register int i2 = i2_61240_lb; i2 < i2_85356_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                  const int i2_29424_lb = max(32 * ii2, i0 + i1 - 531), i2_24845_ub = 32 * ii2 - i0 - i1 + 562;
                  for (register int i2 = i2_29424_lb; i2 <= i2_24845_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  const int i2_28746_lb = 32 * ii2 - i0 - i1 + 563, i2_90139_ub = 32 * ii2 - i0 + i1 + 467;
                  for (register int i2 = i2_28746_lb; i2 <= i2_90139_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                } else {
                  const int i2_52824_lb = 0, i2_11038_ub = min(i0 + i1 - 531, i0 - i1 - 436);
                  for (register int i2 = i2_52824_lb; i2 < i2_11038_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  const int i2_76616_lb = max(-i0 + i1 + 2004, -i0 - i1 + 2099), i2_84354_ub = 1599;
                  for (register int i2 = i2_76616_lb; i2 <= i2_84354_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
              }
              if (ii2 == 0) {
                const int i1_38856_lb = -i0 + 547, i1_18479_ub = i0 - 435;
                for (register int i1 = i1_38856_lb; i1 < i1_18479_ub; i1 += 1) {
                  const int i2_39212_lb = i0 - i1 - 436, i2_2563_ub = -i0 + i1 + 467;
                  for (register int i2 = i2_39212_lb; i2 <= i2_2563_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
              }
            }
          }
        }
      } else {
        for (register int ii2 = 0; ii2 <= 49; ii2 += 1) {
          if (ii0 <= 14) {
            const int i0_51584_lb = 32 * ii0 + 16, i0_37107_ub = 32 * ii0 + 31;
            for (register int i0 = i0_51584_lb; i0 <= i0_37107_ub; i0 += 1) {
              if (ii2 == 49) {
                for (register int i2 = 1568; i2 <= 1599; i2 += 1) {
                  S1(i0, 32 * ii0 - i0 + 31, i2);
                }
              }
              const int i1_88147_lb = max(32 * ii0 + ii2 - i0 - 17, 32 * ii0 - i0 + 31), i1_62087_ub = min(min(15 * ii2 + 15, -32 * ii0 + i0 - 1), -i0 + 498);
              for (register int i1 = i1_88147_lb; i1 <= i1_62087_ub; i1 += 1) {
                if (ii2 <= 48) {
                  const int i2_51883_lb = max(max(0, 32 * ii0 + 32 * ii2 - i0 + i1), 32 * ii0 + 32 * ii2 - i0 - i1 + 31), i2_76248_ub = 32 * ii2;
                  for (register int i2 = i2_51883_lb; i2 < i2_76248_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  if (ii2 == 0) {
                    const int i2_94200_lb = -32 * ii0 + i0 + i1 - 31, i2_45532_ub = -32 * ii0 + i0 - i1;
                    for (register int i2 = i2_94200_lb; i2 < i2_45532_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                  const int i2_91620_lb = max(32 * ii2, -32 * ii0 + i0 - i1), i2_35942_ub = 32 * ii0 + 32 * ii2 - i0 + i1 + 31;
                  for (register int i2 = i2_91620_lb; i2 <= i2_35942_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  const int i2_73353_lb = max(-32 * ii0 + i0 - i1, 32 * ii0 + 32 * ii2 - i0 + i1 + 32), i2_91987_ub = 32 * ii0 + 32 * ii2 - i0 - i1 + 62;
                  for (register int i2 = i2_73353_lb; i2 <= i2_91987_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                } else {
                  const int i2_97887_lb = 0, i2_50945_ub = min(-32 * ii0 + i0 + i1 - 31, -32 * ii0 + i0 - i1);
                  for (register int i2 = i2_97887_lb; i2 < i2_50945_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  const int i2_93695_lb = max(32 * ii0 - i0 + i1 + 1568, 32 * ii0 - i0 - i1 + 1599), i2_27311_ub = 1599;
                  for (register int i2 = i2_93695_lb; i2 <= i2_27311_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
              }
              if (ii0 == 14 && ii2 == 49) {
                const int i1_92143_lb = -i0 + 499, i1_38794_ub = i0 - 448;
                for (register int i1 = i1_92143_lb; i1 < i1_38794_ub; i1 += 1) {
                  const int i2_17451_lb = 0, i2_61319_ub = i0 - i1 - 448;
                  for (register int i2 = i2_17451_lb; i2 < i2_61319_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  const int i2_49832_lb = -i0 + i1 + 2016, i2_10419_ub = 1599;
                  for (register int i2 = i2_49832_lb; i2 <= i2_10419_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
                if (i0 >= 474) {
                  for (register int i2 = 1568; i2 <= 1599; i2 += 1) {
                    S1(i0, i0 - 448, i2);
                  }
                }
              }
              if (ii2 == 49 && 16 * ii0 + 249 >= i0) {
                for (register int i2 = 1568; i2 <= 1599; i2 += 1) {
                  S1(i0, -32 * ii0 + i0, i2);
                }
              } else {
                if (ii0 == 14 && ii2 >= 1 && ii2 <= 48) {
                  const int i1_7604_lb = -i0 + 499, i1_98377_ub = i0 - 448;
                  for (register int i1 = i1_7604_lb; i1 < i1_98377_ub; i1 += 1) {
                    const int i2_54856_lb = 32 * ii2 - i0 + i1 + 448, i2_75541_ub = 32 * ii2 - i0 + i1 + 479;
                    for (register int i2 = i2_54856_lb; i2 <= i2_75541_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                } else if (ii2 == 0) {
                  const int i1_51837_lb = 16, i1_52117_ub = -32 * ii0 + i0;
                  for (register int i1 = i1_51837_lb; i1 < i1_52117_ub; i1 += 1) {
                    const int i2_11911_lb = -32 * ii0 + i0 - i1, i2_39984_ub = 32 * ii0 - i0 + i1 + 31;
                    for (register int i2 = i2_11911_lb; i2 <= i2_39984_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                }
                if (ii2 <= 48) {
                  const int i2_14205_lb = 32 * ii2, i2_80146_ub = 32 * ii2 + 31;
                  for (register int i2 = i2_14205_lb; i2 <= i2_80146_ub; i2 += 1) {
                    S1(i0, -32 * ii0 + i0, i2);
                  }
                }
              }
            }
          } else {
            for (register int i0 = 484; i0 <= 499; i0 += 1) {
              if (ii2 == 0) {
                const int i1_42031_lb = -i0 + 499, i1_40557_ub = i0 - 483;
                for (register int i1 = i1_42031_lb; i1 < i1_40557_ub; i1 += 1) {
                  const int i2_60699_lb = i0 + i1 - 499, i2_31736_ub = -i0 - i1 + 530;
                  for (register int i2 = i2_60699_lb; i2 <= i2_31736_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                }
              }
              const int i1_32544_lb = max(-16 * ii2 + i0 - 483, -i0 + 499), i1_74938_ub = i0 - 467;
              for (register int i1 = i1_32544_lb; i1 < i1_74938_ub; i1 += 1) {
                if (ii2 <= 48) {
                  const int i2_82681_lb = max(max(0, 32 * ii2 - i0 + i1 + 468), 32 * ii2 - i0 - i1 + 499), i2_42592_ub = 32 * ii2;
                  for (register int i2 = i2_82681_lb; i2 < i2_42592_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  if (ii2 == 0) {
                    const int i2_2250_lb = i0 + i1 - 499, i2_74824_ub = i0 - i1 - 468;
                    for (register int i2 = i2_2250_lb; i2 < i2_74824_ub; i2 += 1) {
                      S1(i0, i1, i2);
                    }
                  }
                  const int i2_81386_lb = max(32 * ii2, i0 - i1 - 468), i2_36053_ub = 32 * ii2 - i0 + i1 + 499;
                  for (register int i2 = i2_81386_lb; i2 <= i2_36053_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  const int i2_36143_lb = 32 * ii2 - i0 + i1 + 500, i2_31218_ub = 32 * ii2 - i0 - i1 + 530;
                  for (register int i2 = i2_36143_lb; i2 <= i2_31218_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                } else {
                  const int i2_62824_lb = 0, i2_98169_ub = min(i0 + i1 - 499, i0 - i1 - 468);
                  for (register int i2 = i2_62824_lb; i2 < i2_98169_ub; i2 += 1) {
                    S1(i0, i1, i2);
                  }
                  const int i2_52611_lb = max(-i0 + i1 + 2036, -i0 - i1 + 2067), i2_8074_ub = 1599;
                  for (register int i2 = i2_52611_lb; i2 <= i2_8074_ub; i2 += 1) {
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
