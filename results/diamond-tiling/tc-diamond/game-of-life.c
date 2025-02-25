/*
 * John Conway's Game of Life - 2D, Periodic - B3/S23
 * Adapted from Pochoir test bench
 *
 * Irshad Pananilath: irshad@csa.iisc.ernet.in
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>
#include <time.h>

/*
 * N is the length of one side of the GoL world
 * T is the number of generations to evolve
 */
#ifdef HAS_DECLS
#include "decls.h"
#else
#define N 3200L
#define T 500L
#endif

#ifdef TIME
#define IF_TIME(foo) foo;
#else
#define IF_TIME(foo)
#endif

double t_start, t_end;


#define NUM_FP_OPS 8

/* Define our arrays */
int life[2][N][N];

double rtclock()
{
    struct timezone Tzp;
    struct timeval Tp;
    int stat;
    stat = gettimeofday (&Tp, &Tzp);
    if (stat != 0) printf("Error return from gettimeofday: %d",stat);
    return(Tp.tv_sec + Tp.tv_usec*1.0e-6);
}



/* Subtract the `struct timeval' values X and Y,
 * storing the result in RESULT.
 *
 * Return 1 if the difference is negative, otherwise 0.
 */
int timeval_subtract(struct timeval *result, struct timeval *x, struct timeval *y) {
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

/*
 * Calculate the status of a cell in the next evolution given
 * the neighbors it has
 */
int b2s23(int cell, int neighbors) {
  if((cell == 1 && ((neighbors < 2) || (neighbors > 3)))) {
    return 0;
  }

  if((cell == 1 && (neighbors == 2 || neighbors == 3))) {
    return 1;
  }

  if((cell == 0 && neighbors == 3)) {
    return 1;
  }

  return cell;
}

/*
 * Print the final array
 */
void print_points(int t) {
  int a,b;

  for(a = 0; a < N; a++) {
    for(b = 0; b < N; b++) {
      fprintf(stderr, "%c ", (life[t%2][a][b] ? '1' : '.'));
    }
    fprintf(stderr, "\n");
  }
  fprintf(stderr, "\n");
}

int main(int argc, char * argv[]) {
  long int t, i, j, k;
  long int ip, im, jp, jm;
  int a, b, neighbors, flag;

  // for timekeeping
  int ts_return = -1;
  struct timeval start, end, result;
  double tdiff = 0.0;

  double total = 0.0;

#ifdef DEBUG
  printf("Size of the world = %ldx%ld\nNumber of generations to evolve = %ld\n", N, N, T);
#endif

  /* Initialization */
  srand(42); // seed with a constant value to verify results

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      life[0][i][j] = (rand() & 0x1) ? 1 : 0;
      life[1][i][j] = 0;
    }
  }

  /* // preset patterns */
  /* for (i = 0; i < N; i++) { */
  /*   for (j = 0; j < N; j++) { */
  /*     life[0][i][j] = 0; */
  /*     life[1][i][j] = 0; */
  /*   } */
  /* } */

  /* // toad - oscillator - period 2  */
  /* life[0][2][2] = 1; */
  /* life[0][2][3] = 1; */
  /* life[0][2][4] = 1; */
  /* life[0][3][1] = 1; */
  /* life[0][3][2] = 1; */
  /* life[0][3][3] = 1; */


  /* // pulsar - oscillator - period 3 */
  /* life[0][2][4] = 1; */
  /* life[0][2][5] = 1; */
  /* life[0][2][6] = 1; */
  /* life[0][2][10] = 1; */
  /* life[0][2][11] = 1; */
  /* life[0][2][12] = 1; */
  /* life[0][4][2] = 1; */
  /* life[0][4][7] = 1; */
  /* life[0][4][9] = 1; */
  /* life[0][4][14] = 1; */
  /* life[0][5][2] = 1; */
  /* life[0][5][7] = 1; */
  /* life[0][5][9] = 1; */
  /* life[0][5][14] = 1; */
  /* life[0][6][2] = 1; */
  /* life[0][6][7] = 1; */
  /* life[0][6][9] = 1; */
  /* life[0][6][14] = 1; */
  /* life[0][7][4] = 1; */
  /* life[0][7][5] = 1; */
  /* life[0][7][6] = 1; */
  /* life[0][7][10] = 1; */
  /* life[0][7][11] = 1; */
  /* life[0][7][12] = 1; */
  /* life[0][9][4] = 1; */
  /* life[0][9][5] = 1; */
  /* life[0][9][6] = 1; */
  /* life[0][9][10] = 1; */
  /* life[0][9][11] = 1; */
  /* life[0][9][12] = 1; */
  /* life[0][10][2] = 1; */
  /* life[0][10][7] = 1; */
  /* life[0][10][9] = 1; */
  /* life[0][10][14] = 1; */
  /* life[0][11][2] = 1; */
  /* life[0][11][7] = 1; */
  /* life[0][11][9] = 1; */
  /* life[0][11][14] = 1; */
  /* life[0][12][2] = 1; */
  /* life[0][12][7] = 1; */
  /* life[0][12][9] = 1; */
  /* life[0][12][14] = 1; */
  /* life[0][14][4] = 1; */
  /* life[0][14][5] = 1; */
  /* life[0][14][6] = 1; */
  /* life[0][14][10] = 1; */
  /* life[0][14][11] = 1; */
  /* life[0][14][12] = 1; */

#ifdef DEBUG
  // display the initial setting
  print_points(0);
#endif

    IF_TIME(t_start = rtclock());

/* TC Optimizing Compiler 0.4.1 */
/* ./tc ../examples/pluto/life.scop.c --diamond-tiling --omp-for-codegen --isl-map-tc --inline --debug -b 32 --drop-bounds */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#pragma scop
if (T + N <= 65) {
  for (int ii2 = 0; ii2 <= floord(T + N - 3, 32); ii2 += 1) {
    for (int i0 = 0; i0 < min(T, N - 34); i0 += 1) {
      for (int i1 = i0 + 33; i1 < N - 1; i1 += 1) {
        for (int i2 = max(1, 32 * ii2 - i0); i2 <= min(N - 2, 32 * ii2 - i0 + 31); i2 += 1) {
          life[(i0 + 1) % 2][i1][i2] = b2s23(life[i0 % 2][i1][i2], ((((((life[i0 % 2][i1 - 1][i2 + 1] + life[i0 % 2][i1 - 1][i2]) + life[i0 % 2][i1 - 1][i2 - 1]) + life[i0 % 2][i1][i2 + 1]) + life[i0 % 2][i1][i2 - 1]) + life[i0 % 2][i1 + 1][i2 + 1]) + life[i0 % 2][i1 + 1][i2]) + life[i0 % 2][i1 + 1][i2 - 1]);
        }
      }
    }
  }
  for (int k = 1; k <= min(2, T); k += 1) {
    if (k == 2) {
      for (int ii2 = 0; ii2 <= floord(T + N - 4, 32); ii2 += 1) {
        for (int i0 = max(1, -N + 32 * ii2 + 3); i0 <= min(T - 1, 32 * ii2 + 31); i0 += 1) {
          for (int i1 = 1; i1 <= min(N - 2, i0); i1 += 1) {
            for (int i2 = max(1, 32 * ii2 - i0 + 1); i2 <= min(N - 2, 32 * ii2 - i0 + 32); i2 += 1) {
              life[(i0 + 1) % 2][i1][i2] = b2s23(life[i0 % 2][i1][i2], ((((((life[i0 % 2][i1 - 1][i2 + 1] + life[i0 % 2][i1 - 1][i2]) + life[i0 % 2][i1 - 1][i2 - 1]) + life[i0 % 2][i1][i2 + 1]) + life[i0 % 2][i1][i2 - 1]) + life[i0 % 2][i1 + 1][i2 + 1]) + life[i0 % 2][i1 + 1][i2]) + life[i0 % 2][i1 + 1][i2 - 1]);
            }
          }
        }
      }
    } else {
      for (int ii2 = 0; ii2 <= min(floord(N - 3, 16), floord(T + N - 3, 32)); ii2 += 1) {
        for (int i0 = max(0, -N + 32 * ii2 + 2); i0 < min(T, N - 2); i0 += 1) {
          for (int i1 = i0 + 1; i1 <= min(N - 2, i0 + 32); i1 += 1) {
            for (int i2 = max(1, 32 * ii2 - i0); i2 <= min(N - 2, 32 * ii2 - i0 + 31); i2 += 1) {
              life[(i0 + 1) % 2][i1][i2] = b2s23(life[i0 % 2][i1][i2], ((((((life[i0 % 2][i1 - 1][i2 + 1] + life[i0 % 2][i1 - 1][i2]) + life[i0 % 2][i1 - 1][i2 - 1]) + life[i0 % 2][i1][i2 + 1]) + life[i0 % 2][i1][i2 - 1]) + life[i0 % 2][i1 + 1][i2 + 1]) + life[i0 % 2][i1 + 1][i2]) + life[i0 % 2][i1 + 1][i2 - 1]);
            }
          }
        }
      }
    }
  }
} else if (T >= 1) {
  if (N <= 33) {
    for (int k = 1; k <= 2; k += 1) {
      if (k == 2) {
        for (int ii2 = 0; ii2 <= min(-((-N + 37) / 32) + 2, (T + N - 4) / 32); ii2 += 1) {
          for (int i0 = max(1, -N + 32 * ii2 + 3); i0 <= min(min(61, T - 1), 32 * ii2 + 31); i0 += 1) {
            for (int i1 = 1; i1 <= min(min(N - 2, i0), -i0 + 62); i1 += 1) {
              for (int i2 = max(1, 32 * ii2 - i0 + 1); i2 <= min(N - 2, 32 * ii2 - i0 + 32); i2 += 1) {
                life[(i0 + 1) % 2][i1][i2] = b2s23(life[i0 % 2][i1][i2], ((((((life[i0 % 2][i1 - 1][i2 + 1] + life[i0 % 2][i1 - 1][i2]) + life[i0 % 2][i1 - 1][i2 - 1]) + life[i0 % 2][i1][i2 + 1]) + life[i0 % 2][i1][i2 - 1]) + life[i0 % 2][i1 + 1][i2 + 1]) + life[i0 % 2][i1 + 1][i2]) + life[i0 % 2][i1 + 1][i2 - 1]);
              }
            }
          }
        }
      } else {
        for (int ii2 = 0; ii2 <= floord(N - 3, 16); ii2 += 1) {
          for (int i0 = max(0, -N + 32 * ii2 + 2); i0 < N - 2; i0 += 1) {
            for (int i1 = i0 + 1; i1 < N - 1; i1 += 1) {
              for (int i2 = max(1, 32 * ii2 - i0); i2 <= min(N - 2, 32 * ii2 - i0 + 31); i2 += 1) {
                life[(i0 + 1) % 2][i1][i2] = b2s23(life[i0 % 2][i1][i2], ((((((life[i0 % 2][i1 - 1][i2 + 1] + life[i0 % 2][i1 - 1][i2]) + life[i0 % 2][i1 - 1][i2 - 1]) + life[i0 % 2][i1][i2 + 1]) + life[i0 % 2][i1][i2 - 1]) + life[i0 % 2][i1 + 1][i2 + 1]) + life[i0 % 2][i1 + 1][i2]) + life[i0 % 2][i1 + 1][i2 - 1]);
              }
            }
          }
        }
      }
    }
  }
  for (int ii0 = max(0, -((N + 28) / 31) + 2); ii0 < min(floord(T + 15, 32), (T + N - 2) / 32); ii0 += 1) {
    if (ii0 == 0) {
      #pragma omp parallel for
      for (int ii1 = 0; ii1 < (N - 3) / 32; ii1 += 1) {
        for (int ii2 = 0; ii2 <= (N + 12) / 32; ii2 += 1) {
          for (int i0 = max(0, -N + 32 * ii2 + 2); i0 <= min(14, N - 32 * ii1 - 35); i0 += 1) {
            for (int i1 = 32 * ii1 + i0 + 33; i1 <= min(N - 2, 32 * ii1 - i0 + 62); i1 += 1) {
              for (int i2 = max(1, 32 * ii2 - i0); i2 <= min(N - 2, 32 * ii2 - i0 + 31); i2 += 1) {
                life[(i0 + 1) % 2][i1][i2] = b2s23(life[i0 % 2][i1][i2], ((((((life[i0 % 2][i1 - 1][i2 + 1] + life[i0 % 2][i1 - 1][i2]) + life[i0 % 2][i1 - 1][i2 - 1]) + life[i0 % 2][i1][i2 + 1]) + life[i0 % 2][i1][i2 - 1]) + life[i0 % 2][i1 + 1][i2 + 1]) + life[i0 % 2][i1 + 1][i2]) + life[i0 % 2][i1 + 1][i2 - 1]);
              }
            }
          }
        }
      }
      if (N <= 49) {
        for (int ii2 = 0; ii2 <= min((N - 4) / 32 + 1, (T + N - 3) / 32); ii2 += 1) {
          for (int i0 = max(0, -N + 32 * ii2 + 2); i0 <= min(30, T - 1); i0 += 1) {
            for (int i1 = i0 + 1; i1 <= min(min(N - 2, i0 + 32), -i0 + 62); i1 += 1) {
              for (int i2 = max(1, 32 * ii2 - i0); i2 <= min(N - 2, 32 * ii2 - i0 + 31); i2 += 1) {
                life[(i0 + 1) % 2][i1][i2] = b2s23(life[i0 % 2][i1][i2], ((((((life[i0 % 2][i1 - 1][i2 + 1] + life[i0 % 2][i1 - 1][i2]) + life[i0 % 2][i1 - 1][i2 - 1]) + life[i0 % 2][i1][i2 + 1]) + life[i0 % 2][i1][i2 - 1]) + life[i0 % 2][i1 + 1][i2 + 1]) + life[i0 % 2][i1 + 1][i2]) + life[i0 % 2][i1 + 1][i2 - 1]);
              }
            }
          }
        }
      }
    }
    if (N >= 18 && N + 32 * ii0 >= 50) {
      if (ii0 == 0) {
        for (int ii2 = 0; ii2 <= min((N - 4) / 32 + 1, (T + N - 3) / 32); ii2 += 1) {
          for (int i0 = max(0, -N + 32 * ii2 + 2); i0 <= min(30, T - 1); i0 += 1) {
            for (int i1 = i0 + 1; i1 <= min(i0 + 32, -i0 + 62); i1 += 1) {
              for (int i2 = max(1, 32 * ii2 - i0); i2 <= min(N - 2, 32 * ii2 - i0 + 31); i2 += 1) {
                life[(i0 + 1) % 2][i1][i2] = b2s23(life[i0 % 2][i1][i2], ((((((life[i0 % 2][i1 - 1][i2 + 1] + life[i0 % 2][i1 - 1][i2]) + life[i0 % 2][i1 - 1][i2 - 1]) + life[i0 % 2][i1][i2 + 1]) + life[i0 % 2][i1][i2 - 1]) + life[i0 % 2][i1 + 1][i2 + 1]) + life[i0 % 2][i1 + 1][i2]) + life[i0 % 2][i1 + 1][i2 - 1]);
              }
            }
          }
        }
      }
      #pragma omp parallel for
      for (int ii1 = max(0, -ii0 + 1); ii1 < (N + 14) / 32; ii1 += 1) {
        for (int ii2 = 0; ii2 <= min(min(-ii1 + (N - 3) / 16, (N - 4) / 32 + 1), -ii0 + (T + N - 3) / 32); ii2 += 1) {
          for (int i0 = max(max(32 * ii0, -N + 32 * ii0 + 32 * ii1 + 33), -N + 32 * ii0 + 32 * ii2 + 2); i0 <= min(min(T - 1, 32 * ii0 + 30), N + 32 * ii0 - 32 * ii1 - 3); i0 += 1) {
            for (int i1 = max(-32 * ii0 + 32 * ii1 + i0 + 1, 32 * ii0 + 32 * ii1 - i0 + 31); i1 <= min(min(N - 2, -32 * ii0 + 32 * ii1 + i0 + 32), 32 * ii0 + 32 * ii1 - i0 + 62); i1 += 1) {
              for (int i2 = max(1, 32 * ii0 + 32 * ii2 - i0); i2 <= min(N - 2, 32 * ii0 + 32 * ii2 - i0 + 31); i2 += 1) {
                life[(i0 + 1) % 2][i1][i2] = b2s23(life[i0 % 2][i1][i2], ((((((life[i0 % 2][i1 - 1][i2 + 1] + life[i0 % 2][i1 - 1][i2]) + life[i0 % 2][i1 - 1][i2 - 1]) + life[i0 % 2][i1][i2 + 1]) + life[i0 % 2][i1][i2 - 1]) + life[i0 % 2][i1 + 1][i2 + 1]) + life[i0 % 2][i1 + 1][i2]) + life[i0 % 2][i1 + 1][i2 - 1]);
              }
            }
          }
        }
      }
      if (N >= 34 && 32 * ii0 + 65 >= T + N) {
        for (int ii2 = 0; ii2 < 2; ii2 += 1) {
          for (int i0 = 32 * ii0 + 16; i0 < T; i0 += 1) {
            for (int i1 = 32 * ii0 - i0 + 31; i1 <= -32 * ii0 + i0; i1 += 1) {
              for (int i2 = max(1, 32 * ii0 + 32 * ii2 - i0 + 16); i2 <= min(N - 2, 32 * ii0 + 32 * ii2 - i0 + 47); i2 += 1) {
                life[(i0 + 1) % 2][i1][i2] = b2s23(life[i0 % 2][i1][i2], ((((((life[i0 % 2][i1 - 1][i2 + 1] + life[i0 % 2][i1 - 1][i2]) + life[i0 % 2][i1 - 1][i2 - 1]) + life[i0 % 2][i1][i2 + 1]) + life[i0 % 2][i1][i2 - 1]) + life[i0 % 2][i1 + 1][i2 + 1]) + life[i0 % 2][i1 + 1][i2]) + life[i0 % 2][i1 + 1][i2 - 1]);
              }
            }
          }
        }
      }
    }
    if (N <= 33) {
      for (int ii2 = 0; ii2 <= min((N + 11) / 32 + 1, -ii0 + (T + N + 13) / 32 - 1); ii2 += 1) {
        for (int i0 = max(max(32 * ii0 + 16, -N + 32 * ii0 + 33), -N + 32 * ii0 + 32 * ii2 + 18); i0 <= min(min(T - 1, 32 * ii0 + 61), 32 * ii0 + 32 * ii2 + 46); i0 += 1) {
          for (int i1 = max(1, 32 * ii0 - i0 + 31); i1 <= min(min(N - 2, -32 * ii0 + i0), 32 * ii0 - i0 + 62); i1 += 1) {
            for (int i2 = max(1, 32 * ii0 + 32 * ii2 - i0 + 16); i2 <= min(N - 2, 32 * ii0 + 32 * ii2 - i0 + 47); i2 += 1) {
              life[(i0 + 1) % 2][i1][i2] = b2s23(life[i0 % 2][i1][i2], ((((((life[i0 % 2][i1 - 1][i2 + 1] + life[i0 % 2][i1 - 1][i2]) + life[i0 % 2][i1 - 1][i2 - 1]) + life[i0 % 2][i1][i2 + 1]) + life[i0 % 2][i1][i2 - 1]) + life[i0 % 2][i1 + 1][i2 + 1]) + life[i0 % 2][i1 + 1][i2]) + life[i0 % 2][i1 + 1][i2 - 1]);
            }
          }
        }
      }
    } else if (T + N >= 32 * ii0 + 66) {
      #pragma omp parallel for
      for (int ii1 = 0; ii1 <= min((N - 2) / 32, -ii0 + (T + N - 2) / 32 - 1); ii1 += 1) {
        if (ii1 >= 1) {
          for (int ii2 = 0; ii2 <= min(min(-ii1 + (N + 5) / 16, (N - 4) / 32 + 1), -ii0 + (T + N + 13) / 32 - 1); ii2 += 1) {
            for (int i0 = max(max(32 * ii0 + 16, -N + 32 * ii0 + 32 * ii1 + 33), -N + 32 * ii0 + 32 * ii2 + 18); i0 <= min(min(T - 1, 32 * ii0 + 46), N + 32 * ii0 - 32 * ii1 + 29); i0 += 1) {
              for (int i1 = max(-32 * ii0 + 32 * ii1 + i0 - 31, 32 * ii0 + 32 * ii1 - i0 + 31); i1 <= min(min(N - 2, -32 * ii0 + 32 * ii1 + i0), 32 * ii0 + 32 * ii1 - i0 + 62); i1 += 1) {
                for (int i2 = max(1, 32 * ii0 + 32 * ii2 - i0 + 16); i2 <= min(N - 2, 32 * ii0 + 32 * ii2 - i0 + 47); i2 += 1) {
                  life[(i0 + 1) % 2][i1][i2] = b2s23(life[i0 % 2][i1][i2], ((((((life[i0 % 2][i1 - 1][i2 + 1] + life[i0 % 2][i1 - 1][i2]) + life[i0 % 2][i1 - 1][i2 - 1]) + life[i0 % 2][i1][i2 + 1]) + life[i0 % 2][i1][i2 - 1]) + life[i0 % 2][i1 + 1][i2 + 1]) + life[i0 % 2][i1 + 1][i2]) + life[i0 % 2][i1 + 1][i2 - 1]);
                }
              }
            }
          }
        } else if (ii0 >= 1) {
          for (int ii2 = 0; ii2 <= min((N + 11) / 32 + 1, -ii0 + (T + N + 13) / 32 - 1); ii2 += 1) {
            for (int i0 = max(32 * ii0 + 16, -N + 32 * ii0 + 32 * ii2 + 18); i0 <= min(min(T - 1, 32 * ii0 + 61), 32 * ii0 + 32 * ii2 + 46); i0 += 1) {
              for (int i1 = max(1, 32 * ii0 - i0 + 31); i1 <= min(-32 * ii0 + i0, 32 * ii0 - i0 + 62); i1 += 1) {
                for (int i2 = max(1, 32 * ii0 + 32 * ii2 - i0 + 16); i2 <= min(N - 2, 32 * ii0 + 32 * ii2 - i0 + 47); i2 += 1) {
                  life[(i0 + 1) % 2][i1][i2] = b2s23(life[i0 % 2][i1][i2], ((((((life[i0 % 2][i1 - 1][i2 + 1] + life[i0 % 2][i1 - 1][i2]) + life[i0 % 2][i1 - 1][i2 - 1]) + life[i0 % 2][i1][i2 + 1]) + life[i0 % 2][i1][i2 - 1]) + life[i0 % 2][i1 + 1][i2 + 1]) + life[i0 % 2][i1 + 1][i2]) + life[i0 % 2][i1 + 1][i2 - 1]);
                }
              }
            }
          }
        } else {
          for (int ii2 = 0; ii2 <= min((N - 6) / 32 + 2, (T + N - 4) / 32); ii2 += 1) {
            for (int i0 = max(1, -N + 32 * ii2 + 3); i0 <= min(min(61, T - 1), 32 * ii2 + 31); i0 += 1) {
              for (int i1 = 1; i1 <= min(i0, -i0 + 62); i1 += 1) {
                for (int i2 = max(1, 32 * ii2 - i0 + 1); i2 <= min(N - 2, 32 * ii2 - i0 + 32); i2 += 1) {
                  life[(i0 + 1) % 2][i1][i2] = b2s23(life[i0 % 2][i1][i2], ((((((life[i0 % 2][i1 - 1][i2 + 1] + life[i0 % 2][i1 - 1][i2]) + life[i0 % 2][i1 - 1][i2 - 1]) + life[i0 % 2][i1][i2 + 1]) + life[i0 % 2][i1][i2 - 1]) + life[i0 % 2][i1 + 1][i2 + 1]) + life[i0 % 2][i1 + 1][i2]) + life[i0 % 2][i1 + 1][i2 - 1]);
                }
              }
            }
          }
        }
      }
    }
  }
  if (((T + 15) % 32) + N >= 49 && (T + 15) % 32 >= 16) {
    if (T <= 16) {
      #pragma omp parallel for
      for (int ii1 = 0; ii1 < (N - 3) / 32; ii1 += 1) {
        for (int ii2 = 0; ii2 <= min((T + N - 3) / 32, (N + 12) / 32); ii2 += 1) {
          for (int i0 = max(0, -N + 32 * ii2 + 2); i0 <= min(min(14, T - 1), N - 32 * ii1 - 35); i0 += 1) {
            for (int i1 = 32 * ii1 + i0 + 33; i1 <= min(N - 2, 32 * ii1 - i0 + 62); i1 += 1) {
              for (int i2 = max(1, 32 * ii2 - i0); i2 <= min(N - 2, 32 * ii2 - i0 + 31); i2 += 1) {
                life[(i0 + 1) % 2][i1][i2] = b2s23(life[i0 % 2][i1][i2], ((((((life[i0 % 2][i1 - 1][i2 + 1] + life[i0 % 2][i1 - 1][i2]) + life[i0 % 2][i1 - 1][i2 - 1]) + life[i0 % 2][i1][i2 + 1]) + life[i0 % 2][i1][i2 - 1]) + life[i0 % 2][i1 + 1][i2 + 1]) + life[i0 % 2][i1 + 1][i2]) + life[i0 % 2][i1 + 1][i2 - 1]);
              }
            }
          }
        }
      }
    }
    #pragma omp parallel for
    for (int ii1 = 0; ii1 < -((T - 1) / 32) + (T + N - 2) / 32; ii1 += 1) {
      for (int ii2 = 0; ii2 <= -((T - 1) / 32) + (T + N - 3) / 32; ii2 += 1) {
        for (int i0 = max(max(max(-((T - 1) % 32) + T - 1, -((T - 1) % 32) + T - N + 32 * ii2 + 1), -16 * ii1 + 16 * ((T - 1) / 32) + 16), -((T - 1) % 32) + T - N + 32 * ii1 + 32); i0 < T; i0 += 1) {
          for (int i1 = -((T - 1) % 32) + T + 32 * ii1 - i0 + 30; i1 <= min(N - 2, ((T - 1) % 32) - T + 32 * ii1 + i0 + 33); i1 += 1) {
            for (int i2 = max(1, -((T - 1) % 32) + T + 32 * ii2 - i0 - 1); i2 <= min(N - 2, -((T - 1) % 32) + T + 32 * ii2 - i0 + 30); i2 += 1) {
              life[(i0 + 1) % 2][i1][i2] = b2s23(life[i0 % 2][i1][i2], ((((((life[i0 % 2][i1 - 1][i2 + 1] + life[i0 % 2][i1 - 1][i2]) + life[i0 % 2][i1 - 1][i2 - 1]) + life[i0 % 2][i1][i2 + 1]) + life[i0 % 2][i1][i2 - 1]) + life[i0 % 2][i1 + 1][i2 + 1]) + life[i0 % 2][i1 + 1][i2]) + life[i0 % 2][i1 + 1][i2 - 1]);
            }
          }
        }
        if (T <= 16 && ii1 == 0) {
          for (int i0 = max(0, -N + 32 * ii2 + 2); i0 < T; i0 += 1) {
            for (int i1 = i0 + 1; i1 <= i0 + 32; i1 += 1) {
              for (int i2 = max(1, 32 * ii2 - i0); i2 <= min(N - 2, 32 * ii2 - i0 + 31); i2 += 1) {
                life[(i0 + 1) % 2][i1][i2] = b2s23(life[i0 % 2][i1][i2], ((((((life[i0 % 2][i1 - 1][i2 + 1] + life[i0 % 2][i1 - 1][i2]) + life[i0 % 2][i1 - 1][i2 - 1]) + life[i0 % 2][i1][i2 + 1]) + life[i0 % 2][i1][i2 - 1]) + life[i0 % 2][i1 + 1][i2 + 1]) + life[i0 % 2][i1 + 1][i2]) + life[i0 % 2][i1 + 1][i2 - 1]);
              }
            }
          }
        }
      }
    }
    if (T <= 16) {
      for (int ii2 = 0; ii2 <= (T + N - 4) / 32; ii2 += 1) {
        for (int i0 = max(1, -N + 32 * ii2 + 3); i0 < T; i0 += 1) {
          for (int i1 = 1; i1 <= i0; i1 += 1) {
            for (int i2 = max(1, 32 * ii2 - i0 + 1); i2 <= min(N - 2, 32 * ii2 - i0 + 32); i2 += 1) {
              life[(i0 + 1) % 2][i1][i2] = b2s23(life[i0 % 2][i1][i2], ((((((life[i0 % 2][i1 - 1][i2 + 1] + life[i0 % 2][i1 - 1][i2]) + life[i0 % 2][i1 - 1][i2 - 1]) + life[i0 % 2][i1][i2 + 1]) + life[i0 % 2][i1][i2 - 1]) + life[i0 % 2][i1 + 1][i2 + 1]) + life[i0 % 2][i1 + 1][i2]) + life[i0 % 2][i1 + 1][i2 - 1]);
            }
          }
        }
      }
    }
  }
}
#pragma endscop

    IF_TIME(t_end = rtclock());
    IF_TIME(fprintf(stdout, "%0.6lfs\n", t_end - t_start));


  if (fopen(".test", "r")) {
      print_points(T);
  }

  return 0;
}

// icc -O3 -DTIME -DDEBUG life.c -o op-life
