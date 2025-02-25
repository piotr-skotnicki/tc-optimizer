#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "util.h"

#define N 2000000
#define T 1000

#pragma declarations
double a[N];
double b[N];
#pragma enddeclarations

double t_start, t_end;

void init_array()
{
    int j;

    for (j=0; j<N; j++) {
        a[j] = ((double)j)/N;
    }
}


void print_array()
{
    int j;

    for (j=0; j<N; j++) {
        fprintf(stderr, "%lf ", a[j]);
        if (j%80 == 20) fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
}


int main()
{
    int t, i, j;

    init_array();

    IF_TIME(t_start = rtclock());

/* TC Optimizing Compiler 0.4.1 */
/* ./tc ../examples/pluto-perfect/jacobi-1d-imper.scop.c --diamond-tiling --omp-for-codegen --isl-map-tc --inline --debug -b 32 --drop-bounds */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#pragma scop
if (N >= 4) {
  for (int ii0 = 0; ii0 <= floord(T - 1, 32); ii0 += 1) {
    if (ii0 == 0) {
      #pragma omp parallel for
      for (int ii1 = 0; ii1 <= floord(N - 6, 32); ii1 += 1) {
        for (int i0 = 0; i0 <= min(min(7, T - 1), -16 * ii1 + N / 2 - 3); i0 += 1) {
          for (int i1 = 32 * ii1 + 4 * i0 + 2; i1 <= min(min(2 * N - 32 * ii1 - 10, 32 * ii1 + 31), N + 2 * i0 - 3); i1 += 1) {
            if (N + 2 * i0 >= i1 + 4) {
              b[-2 * i0 + i1 + 2] = (0.33333 * ((a[-2 * i0 + i1 + 1] + a[-2 * i0 + i1 + 2]) + a[-2 * i0 + i1 + 3]));
            }
            if (i1 >= 32 * ii1 + 4 * i0 + 4) {
              a[-2 * i0 + i1 + 1] = b[-2 * i0 + i1 + 1];
            }
          }
        }
      }
    }
    #pragma omp parallel for
    for (int ii1 = 2 * ii0; ii1 <= min(2 * ii0 + (N + 11) / 32, (N + 2 * T - 5) / 32); ii1 += 1) {
      for (int i0 = max(32 * ii0, -N + 16 * ii1 + N / 2 + 2); i0 <= min(min(T - 1, 32 * ii0 + 15), 64 * ii0 - 16 * ii1 + N / 2 + 13); i0 += 1) {
        {
          for (int i1 = max(max(64 * ii0 + 32, 32 * ii1), -128 * ii0 + 32 * ii1 + 4 * i0 - 30); i1 < min(-128 * ii0 + 32 * ii1 + 4 * i0 - 28, N + 2 * i0 - 3); i1 += 1) {
            b[-2 * i0 + i1 + 2] = (0.33333 * ((a[-2 * i0 + i1 + 1] + a[-2 * i0 + i1 + 2]) + a[-2 * i0 + i1 + 3]));
          }
          for (int i1 = max(max(64 * ii0 + 32, 32 * ii1), -128 * ii0 + 32 * ii1 + 4 * i0 - 28); i1 <= min(min(32 * ii1 + 31, N + 2 * i0 - 3), -128 * ii0 + 32 * ii1 + 4 * i0 + 3); i1 += 1) {
            if (N + 2 * i0 >= i1 + 4 && 32 * ii1 + 4 * i0 + 1 >= 128 * ii0 + i1) {
              b[-2 * i0 + i1 + 2] = (0.33333 * ((a[-2 * i0 + i1 + 1] + a[-2 * i0 + i1 + 2]) + a[-2 * i0 + i1 + 3]));
            }
            a[-2 * i0 + i1 + 1] = b[-2 * i0 + i1 + 1];
          }
        }
        if (ii1 == 2 * ii0) {
          b[2] = (0.33333 * ((a[1] + a[2]) + a[3]));
          if (N >= 5 && i0 == 32 * ii0 + 15) {
            b[3] = (0.33333 * ((a[2] + a[3]) + a[4]));
          }
          for (int i1 = max(-64 * ii0 + 4 * i0 - 28, 2 * i0 + 1); i1 <= min(min(64 * ii0 + 31, N + 2 * i0 - 3), -64 * ii0 + 4 * i0 + 3); i1 += 1) {
            if (N + 2 * i0 >= i1 + 4 && 4 * i0 + 1 >= 64 * ii0 + i1) {
              b[-2 * i0 + i1 + 2] = (0.33333 * ((a[-2 * i0 + i1 + 1] + a[-2 * i0 + i1 + 2]) + a[-2 * i0 + i1 + 3]));
            }
            a[-2 * i0 + i1 + 1] = b[-2 * i0 + i1 + 1];
          }
        }
      }
    }
    if (N == 4 && T >= 32 * ii0 + 16) {
      a[2] = b[2];
    }
    if (N + 2 * T >= 64 * ii0 + 37) {
      for (int k = max(2, -N + 7); k <= min(4, -4 * ii0 + (T - 1) / 8 + 1); k += 1) {
        #pragma omp parallel for
        for (int ii1 = max(2 * ii0 + (k + 1) / 2 - 1, 4 * ii0 + k - (T - 4 * ii0 - k + 14) / 14); ii1 <= min((N + 2 * T - 5) / 32, 2 * ii0 + (N + 16 * k - 5) / 32); ii1 += 1) {
          if (k <= 3 && ii1 >= 2 * ii0 + 1) {
            if (k == 3) {
              for (int i0 = max(32 * ii0 + 16, -N + 16 * ii1 + N / 2 + 2); i0 <= min(min(T - 1, 32 * ii0 + 31), 64 * ii0 - 16 * ii1 + N / 2 + 45); i0 += 1) {
                if (ii1 == 2 * ii0 + 1) {
                  b[2] = (0.33333 * ((a[1] + a[2]) + a[3]));
                }
                for (int i1 = max(max(32 * ii1, -128 * ii0 + 32 * ii1 + 4 * i0 - 94), 2 * i0 + 1); i1 < min(-128 * ii0 + 32 * ii1 + 4 * i0 - 92, N + 2 * i0 - 3); i1 += 1) {
                  b[-2 * i0 + i1 + 2] = (0.33333 * ((a[-2 * i0 + i1 + 1] + a[-2 * i0 + i1 + 2]) + a[-2 * i0 + i1 + 3]));
                }
                for (int i1 = max(max(32 * ii1, -128 * ii0 + 32 * ii1 + 4 * i0 - 92), 2 * i0 + 1); i1 <= min(min(32 * ii1 + 31, -128 * ii0 + 32 * ii1 + 4 * i0 - 61), N + 2 * i0 - 3); i1 += 1) {
                  if (32 * ii1 + 4 * i0 >= 128 * ii0 + i1 + 63 && N + 2 * i0 >= i1 + 4) {
                    b[-2 * i0 + i1 + 2] = (0.33333 * ((a[-2 * i0 + i1 + 1] + a[-2 * i0 + i1 + 2]) + a[-2 * i0 + i1 + 3]));
                  }
                  a[-2 * i0 + i1 + 1] = b[-2 * i0 + i1 + 1];
                }
              }
            } else {
              for (int i0 = max(32 * ii0 + 8, -N + 16 * ii1 + N / 2 + 2); i0 <= min(min(T - 1, 32 * ii0 + 23), 64 * ii0 - 16 * ii1 + N / 2 + 29); i0 += 1) {
                for (int i1 = max(32 * ii1, -128 * ii0 + 32 * ii1 + 4 * i0 - 62); i1 < min(-128 * ii0 + 32 * ii1 + 4 * i0 - 60, N + 2 * i0 - 3); i1 += 1) {
                  b[-2 * i0 + i1 + 2] = (0.33333 * ((a[-2 * i0 + i1 + 1] + a[-2 * i0 + i1 + 2]) + a[-2 * i0 + i1 + 3]));
                }
                for (int i1 = max(32 * ii1, -128 * ii0 + 32 * ii1 + 4 * i0 - 60); i1 <= min(min(32 * ii1 + 31, -128 * ii0 + 32 * ii1 + 4 * i0 - 29), N + 2 * i0 - 3); i1 += 1) {
                  if (32 * ii1 + 4 * i0 >= 128 * ii0 + i1 + 31 && N + 2 * i0 >= i1 + 4) {
                    b[-2 * i0 + i1 + 2] = (0.33333 * ((a[-2 * i0 + i1 + 1] + a[-2 * i0 + i1 + 2]) + a[-2 * i0 + i1 + 3]));
                  }
                  a[-2 * i0 + i1 + 1] = b[-2 * i0 + i1 + 1];
                }
              }
            }
          } else if (k == 4) {
            for (int i0 = max(max(32 * ii0 + 24, 64 * ii0 - 16 * ii1 + 47), -N + 16 * ii1 + N / 2 + 2); i0 <= min(min(min(T - 1, 32 * ii0 + 39), 16 * ii1 + 15), 64 * ii0 - 16 * ii1 + N / 2 + 61); i0 += 1) {
              for (int i1 = max(32 * ii1, -128 * ii0 + 32 * ii1 + 4 * i0 - 126); i1 < min(-128 * ii0 + 32 * ii1 + 4 * i0 - 124, N + 2 * i0 - 3); i1 += 1) {
                b[-2 * i0 + i1 + 2] = (0.33333 * ((a[-2 * i0 + i1 + 1] + a[-2 * i0 + i1 + 2]) + a[-2 * i0 + i1 + 3]));
              }
              for (int i1 = max(max(32 * ii1, -128 * ii0 + 32 * ii1 + 4 * i0 - 124), 2 * i0 + 1); i1 <= min(min(32 * ii1 + 31, -128 * ii0 + 32 * ii1 + 4 * i0 - 93), N + 2 * i0 - 3); i1 += 1) {
                if (32 * ii1 + 4 * i0 >= 128 * ii0 + i1 + 95 && N + 2 * i0 >= i1 + 4) {
                  b[-2 * i0 + i1 + 2] = (0.33333 * ((a[-2 * i0 + i1 + 1] + a[-2 * i0 + i1 + 2]) + a[-2 * i0 + i1 + 3]));
                }
                a[-2 * i0 + i1 + 1] = b[-2 * i0 + i1 + 1];
              }
            }
          } else {
            a[2] = b[2];
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
#ifdef MPI
        if (my_rank == 0) {
            print_array();
        }
#else
        print_array();
#endif
    }

    return 0;
}
