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
/* ./tc ../examples/pluto-perfect/jacobi-1d-imper.scop.c --semi-diamond-tiling --omp-for-codegen --isl-map-tc --inline --debug -b 16 --drop-bounds */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#pragma scop
if (N >= 4) {
  for (int ii0 = 0; ii0 <= floord(T - 1, 16); ii0 += 1) {
    for (int k = max(0, -N + (N + 1) / 2 + 3); k <= 1; k += 1) {
      #pragma omp parallel for
      for (int ii1 = 2 * ii0; ii1 <= min((N + 2 * T - 5) / 16, 2 * ii0 + k + floord(N - 7 * k - 6, 16)); ii1 += 1) {
        if (k == 1) {
          for (int i0 = max(16 * ii0, -N + 8 * ii1 + N / 2 + 2); i0 <= min(min(T - 1, 16 * ii0 + 7), 32 * ii0 - 8 * ii1 + N / 2 + 5); i0 += 1) {
            {
              for (int i1 = max(max(32 * ii0 + 16, 16 * ii1), -64 * ii0 + 16 * ii1 + 4 * i0 - 14); i1 < min(-64 * ii0 + 16 * ii1 + 4 * i0 - 12, N + 2 * i0 - 3); i1 += 1) {
                b[-2 * i0 + i1 + 2] = (0.33333 * ((a[-2 * i0 + i1 + 1] + a[-2 * i0 + i1 + 2]) + a[-2 * i0 + i1 + 3]));
              }
              for (int i1 = max(max(32 * ii0 + 16, 16 * ii1), -64 * ii0 + 16 * ii1 + 4 * i0 - 12); i1 <= min(min(16 * ii1 + 15, N + 2 * i0 - 3), -64 * ii0 + 16 * ii1 + 4 * i0 + 3); i1 += 1) {
                if (N + 2 * i0 >= i1 + 4 && 16 * ii1 + 4 * i0 + 1 >= 64 * ii0 + i1) {
                  b[-2 * i0 + i1 + 2] = (0.33333 * ((a[-2 * i0 + i1 + 1] + a[-2 * i0 + i1 + 2]) + a[-2 * i0 + i1 + 3]));
                }
                a[-2 * i0 + i1 + 1] = b[-2 * i0 + i1 + 1];
              }
            }
            if (ii1 == 2 * ii0) {
              b[2] = (0.33333 * ((a[1] + a[2]) + a[3]));
              if (N >= 5 && i0 == 16 * ii0 + 7) {
                b[3] = (0.33333 * ((a[2] + a[3]) + a[4]));
              }
              for (int i1 = max(-32 * ii0 + 4 * i0 - 12, 2 * i0 + 1); i1 <= min(min(32 * ii0 + 15, N + 2 * i0 - 3), -32 * ii0 + 4 * i0 + 3); i1 += 1) {
                if (N + 2 * i0 >= i1 + 4 && 4 * i0 + 1 >= 32 * ii0 + i1) {
                  b[-2 * i0 + i1 + 2] = (0.33333 * ((a[-2 * i0 + i1 + 1] + a[-2 * i0 + i1 + 2]) + a[-2 * i0 + i1 + 3]));
                }
                a[-2 * i0 + i1 + 1] = b[-2 * i0 + i1 + 1];
              }
            }
          }
        } else {
          for (int i0 = 16 * ii0; i0 <= min(min(T - 1, 16 * ii0 + 3), 32 * ii0 - 8 * ii1 + N / 2 - 3); i0 += 1) {
            for (int i1 = -64 * ii0 + 16 * ii1 + 4 * i0 + 2; i1 <= min(min(2 * N + 64 * ii0 - 16 * ii1 - 10, 16 * ii1 + 15), N + 2 * i0 - 3); i1 += 1) {
              if (N + 2 * i0 >= i1 + 4) {
                b[-2 * i0 + i1 + 2] = (0.33333 * ((a[-2 * i0 + i1 + 1] + a[-2 * i0 + i1 + 2]) + a[-2 * i0 + i1 + 3]));
              }
              if (64 * ii0 + i1 >= 16 * ii1 + 4 * i0 + 4) {
                a[-2 * i0 + i1 + 1] = b[-2 * i0 + i1 + 1];
              }
            }
          }
        }
      }
    }
    if (T >= 16 * ii0 + 8) {
      for (int k = 2; k <= min(min(min(4, N / 2 + 1), -4 * ii0 + (T - 1) / 4 + 1), -2 * ii0 + (N + 2 * T - 6) / 16 + 2); k += 1) {
        if (k == 2) {
          a[2] = b[2];
        }
        #pragma omp parallel for
        for (int ii1 = max(2 * ii0 + 1, 4 * ii0 + k - (T + k + 4) / 8); ii1 <= min((N + 2 * T - 5) / 16, 2 * ii0 + (N + 8 * k - 5) / 16); ii1 += 1) {
          if (k <= 3) {
            if (k == 3) {
              for (int i0 = max(16 * ii0 + 8, -N + 8 * ii1 + N / 2 + 2); i0 <= min(min(T - 1, 16 * ii0 + 15), 32 * ii0 - 8 * ii1 + N / 2 + 21); i0 += 1) {
                if (ii1 == 2 * ii0 + 1) {
                  b[2] = (0.33333 * ((a[1] + a[2]) + a[3]));
                }
                for (int i1 = max(max(16 * ii1, -64 * ii0 + 16 * ii1 + 4 * i0 - 46), 2 * i0 + 1); i1 < min(-64 * ii0 + 16 * ii1 + 4 * i0 - 44, N + 2 * i0 - 3); i1 += 1) {
                  b[-2 * i0 + i1 + 2] = (0.33333 * ((a[-2 * i0 + i1 + 1] + a[-2 * i0 + i1 + 2]) + a[-2 * i0 + i1 + 3]));
                }
                for (int i1 = max(max(16 * ii1, -64 * ii0 + 16 * ii1 + 4 * i0 - 44), 2 * i0 + 1); i1 <= min(min(16 * ii1 + 15, -64 * ii0 + 16 * ii1 + 4 * i0 - 29), N + 2 * i0 - 3); i1 += 1) {
                  if (16 * ii1 + 4 * i0 >= 64 * ii0 + i1 + 31 && N + 2 * i0 >= i1 + 4) {
                    b[-2 * i0 + i1 + 2] = (0.33333 * ((a[-2 * i0 + i1 + 1] + a[-2 * i0 + i1 + 2]) + a[-2 * i0 + i1 + 3]));
                  }
                  a[-2 * i0 + i1 + 1] = b[-2 * i0 + i1 + 1];
                }
              }
            } else {
              for (int i0 = max(16 * ii0 + 4, -N + 8 * ii1 + N / 2 + 2); i0 <= min(min(T - 1, 16 * ii0 + 11), 32 * ii0 - 8 * ii1 + N / 2 + 13); i0 += 1) {
                for (int i1 = max(16 * ii1, -64 * ii0 + 16 * ii1 + 4 * i0 - 30); i1 < min(-64 * ii0 + 16 * ii1 + 4 * i0 - 28, N + 2 * i0 - 3); i1 += 1) {
                  b[-2 * i0 + i1 + 2] = (0.33333 * ((a[-2 * i0 + i1 + 1] + a[-2 * i0 + i1 + 2]) + a[-2 * i0 + i1 + 3]));
                }
                for (int i1 = max(16 * ii1, -64 * ii0 + 16 * ii1 + 4 * i0 - 28); i1 <= min(min(16 * ii1 + 15, -64 * ii0 + 16 * ii1 + 4 * i0 - 13), N + 2 * i0 - 3); i1 += 1) {
                  if (16 * ii1 + 4 * i0 >= 64 * ii0 + i1 + 15 && N + 2 * i0 >= i1 + 4) {
                    b[-2 * i0 + i1 + 2] = (0.33333 * ((a[-2 * i0 + i1 + 1] + a[-2 * i0 + i1 + 2]) + a[-2 * i0 + i1 + 3]));
                  }
                  a[-2 * i0 + i1 + 1] = b[-2 * i0 + i1 + 1];
                }
              }
            }
          } else {
            for (int i0 = max(max(16 * ii0 + 12, 32 * ii0 - 8 * ii1 + 23), -N + 8 * ii1 + N / 2 + 2); i0 <= min(T - 1, 16 * ii0 + 15); i0 += 1) {
              for (int i1 = max(16 * ii1, 64 * ii0 - 16 * ii1 + 47); i1 < min(-64 * ii0 + 16 * ii1 + 4 * i0 - 44, N + 2 * i0 - 2); i1 += 1) {
                if (16 * ii1 + 4 * i0 >= 64 * ii0 + i1 + 47 && N + 2 * i0 >= i1 + 4) {
                  b[-2 * i0 + i1 + 2] = (0.33333 * ((a[-2 * i0 + i1 + 1] + a[-2 * i0 + i1 + 2]) + a[-2 * i0 + i1 + 3]));
                }
                a[-2 * i0 + i1 + 1] = b[-2 * i0 + i1 + 1];
              }
            }
          }
        }
      }
      if (N >= 7 && N <= 11 && 32 * ii0 + 37 == N + 2 * T) {
        a[N - 2] = b[N - 2];
      } else if (N <= 5 && T >= 16 * ii0 + 16) {
        #pragma omp parallel for
        for (int ii1 = 2 * ii0 + 1; ii1 < N + 2 * ii0 - 2; ii1 += 1) {
          a[-2 * ii0 + ii1 + 1] = b[-2 * ii0 + ii1 + 1];
        }
      }
    } else {
      #pragma omp parallel for
      for (int ii1 = 6 * ii0 - (T + 3) / 4 + 3; ii1 <= floord(N + 2 * T - 5, 16); ii1 += 1) {
        for (int i0 = max(16 * ii0 + 4, -N + 8 * ii1 + N / 2 + 2); i0 < T; i0 += 1) {
          for (int i1 = 16 * ii1; i1 < min(-64 * ii0 + 16 * ii1 + 4 * i0 - 12, N + 2 * i0 - 2); i1 += 1) {
            if (16 * ii1 + 4 * i0 >= 64 * ii0 + i1 + 15 && N + 2 * i0 >= i1 + 4) {
              b[-2 * i0 + i1 + 2] = (0.33333 * ((a[-2 * i0 + i1 + 1] + a[-2 * i0 + i1 + 2]) + a[-2 * i0 + i1 + 3]));
            }
            a[-2 * i0 + i1 + 1] = b[-2 * i0 + i1 + 1];
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
