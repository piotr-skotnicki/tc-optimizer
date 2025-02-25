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
#define N 1600000L
#define T 1000L
#endif

#define NUM_FP_OPS 4

/* Define our arrays */
double A[2][N + 2];
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
  long int t, i, j, k;
  const int BASE = 1024;

  // for timekeeping
  int ts_return = -1;
  struct timeval start, end, result;
  double tdiff = 0.0;

  printf("Number of points = %ld\t|Number of timesteps = %ld\t", N, T);

  /* Initialization */
  srand(42); // seed with a constant value to verify results

  for (i = 1; i < N + 1; i++) {
    A[0][i] = 1.0 * (rand() % BASE);
  }

#ifdef TIME
  gettimeofday(&start, 0);
#endif

/* TC Optimizing Compiler 0.4.1 */
/* ./tc ../examples/pluto/heat-1d.scop.c --diamond-tiling --omp-for-codegen --isl-map-tc --inline --debug -b 128 --drop-bounds */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#pragma scop
if (T + N >= 256) {
  if (N <= 127) {
    for (register int k = 1; k <= 2; k += 1) {
      if (k == 2) {
        for (register int i0 = 1; i0 <= min(253, T - 1); i0 += 1) {
          for (register int i1 = 1; i1 <= min(min(N, i0), -i0 + 254); i1 += 1) {
            A[(i0 + 1) % 2][i1] = (0.250 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
          }
        }
      } else {
        for (register int i0 = 0; i0 < N; i0 += 1) {
          for (register int i1 = i0 + 1; i1 <= N; i1 += 1) {
            A[(i0 + 1) % 2][i1] = (0.250 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
          }
        }
      }
    }
    if (N <= 63) {
      for (register int ii0 = 1; ii0 < floord(T + N, 128); ii0 += 1) {
        for (register int i0 = -N + 128 * ii0 + 127; i0 <= min(T - 1, 128 * ii0 + 253); i0 += 1) {
          for (register int i1 = max(1, 128 * ii0 - i0 + 127); i1 <= min(N, 128 * ii0 - i0 + 254); i1 += 1) {
            A[(i0 + 1) % 2][i1] = (0.250 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
          }
        }
      }
    }
  } else if (N <= 191) {
    for (register int k = max(0, -N + 129); k <= 2; k += 1) {
      if (k >= 1) {
        if (k == 1) {
          for (register int i0 = 0; i0 <= min(126, T - 1); i0 += 1) {
            for (register int i1 = i0 + 1; i1 <= min(min(N, i0 + 128), -i0 + 254); i1 += 1) {
              A[(i0 + 1) % 2][i1] = (0.250 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
            }
          }
        } else {
          for (register int i0 = 1; i0 <= min(253, T - 1); i0 += 1) {
            for (register int i1 = 1; i1 <= min(i0, -i0 + 254); i1 += 1) {
              A[(i0 + 1) % 2][i1] = (0.250 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
            }
          }
        }
      } else {
        for (register int i0 = 0; i0 < N - 128; i0 += 1) {
          for (register int i1 = i0 + 129; i1 <= N; i1 += 1) {
            A[(i0 + 1) % 2][i1] = (0.250 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
          }
        }
      }
      if (k == 2) {
        for (register int i0 = -N + 255; i0 < min(T, N); i0 += 1) {
          for (register int i1 = max(i0 + 1, -i0 + 255); i1 <= N; i1 += 1) {
            A[(i0 + 1) % 2][i1] = (0.250 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
          }
        }
      }
    }
  }
  if (N >= 64) {
    for (register int ii0 = max(0, -((N + 64) / 128) + 2); ii0 <= min(floord(T - 1, 128), (T + N) / 128 - 1); ii0 += 1) {
      if (ii0 == 0) {
        #pragma omp parallel for
        for (register int ii1 = 0; ii1 < (N - 1) / 128; ii1 += 1) {
          for (register int i0 = 0; i0 <= min(min(62, T - 1), N - 128 * ii1 - 129); i0 += 1) {
            for (register int i1 = 128 * ii1 + i0 + 129; i1 <= min(N, 128 * ii1 - i0 + 254); i1 += 1) {
              A[(i0 + 1) % 2][i1] = (0.250 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
            }
          }
        }
      }
      for (register int k = 1; k <= min(2, -2 * ii0 + (T - 1) / 64 + 1); k += 1) {
        if (ii0 == 0 && k == 1) {
          for (register int i0 = 0; i0 <= min(126, T - 1); i0 += 1) {
            for (register int i1 = i0 + 1; i1 <= min(i0 + 128, -i0 + 254); i1 += 1) {
              A[(i0 + 1) % 2][i1] = (0.250 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
            }
          }
        } else if (ii0 == 0 && k == 2) {
          for (register int i0 = 1; i0 <= min(253, T - 1); i0 += 1) {
            for (register int i1 = 1; i1 <= min(i0, -i0 + 254); i1 += 1) {
              A[(i0 + 1) % 2][i1] = (0.250 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
            }
          }
        }
        #pragma omp parallel for
        for (register int ii1 = max(0, -ii0 + 1); ii1 < min(-ii0 + (T + N) / 128, (N + 64 * k) / 128); ii1 += 1) {
          if (k == 1) {
            for (register int i0 = max(128 * ii0, -N + 128 * ii0 + 128 * ii1 + 127); i0 <= min(min(T - 1, 128 * ii0 + 126), N + 128 * ii0 - 128 * ii1 - 1); i0 += 1) {
              for (register int i1 = max(-128 * ii0 + 128 * ii1 + i0 + 1, 128 * ii0 + 128 * ii1 - i0 + 127); i1 <= min(min(N, -128 * ii0 + 128 * ii1 + i0 + 128), 128 * ii0 + 128 * ii1 - i0 + 254); i1 += 1) {
                A[(i0 + 1) % 2][i1] = (0.250 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
              }
            }
          } else if (ii1 >= 1) {
            for (register int i0 = max(128 * ii0 + 64, -N + 128 * ii0 + 128 * ii1 + 127); i0 <= min(min(T - 1, 128 * ii0 + 190), N + 128 * ii0 - 128 * ii1 + 127); i0 += 1) {
              for (register int i1 = max(-128 * ii0 + 128 * ii1 + i0 - 127, 128 * ii0 + 128 * ii1 - i0 + 127); i1 <= min(min(N, -128 * ii0 + 128 * ii1 + i0), 128 * ii0 + 128 * ii1 - i0 + 254); i1 += 1) {
                A[(i0 + 1) % 2][i1] = (0.250 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
              }
            }
          } else {
            for (register int i0 = 128 * ii0 + 64; i0 <= min(T - 1, 128 * ii0 + 253); i0 += 1) {
              for (register int i1 = max(1, 128 * ii0 - i0 + 127); i1 <= min(min(N, -128 * ii0 + i0), 128 * ii0 - i0 + 254); i1 += 1) {
                A[(i0 + 1) % 2][i1] = (0.250 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
              }
            }
          }
        }
      }
      if (T <= 64 && ii0 == 0) {
        for (register int i0 = 1; i0 < T; i0 += 1) {
          for (register int i1 = 1; i1 <= i0; i1 += 1) {
            A[(i0 + 1) % 2][i1] = (0.250 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
          }
        }
      }
    }
  }
} else {
  for (register int i0 = 0; i0 < min(T, N - 128); i0 += 1) {
    for (register int i1 = i0 + 129; i1 <= N; i1 += 1) {
      A[(i0 + 1) % 2][i1] = (0.250 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
    }
  }
  for (register int k = 1; k <= min(2, T); k += 1) {
    if (k == 2) {
      for (register int i0 = 1; i0 < T; i0 += 1) {
        for (register int i1 = 1; i1 <= min(N, i0); i1 += 1) {
          A[(i0 + 1) % 2][i1] = (0.250 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
        }
      }
    } else {
      for (register int i0 = 0; i0 < min(T, N); i0 += 1) {
        for (register int i1 = i0 + 1; i1 <= min(N, i0 + 128); i1 += 1) {
          A[(i0 + 1) % 2][i1] = (0.250 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
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
  printf("|MFLOPS =  %f\t",
         ((((double)NUM_FP_OPS * N * T) / tdiff) / 1000000L));
#endif

#ifdef VERIFY
  for (i = 1; i < N + 1; i++) {
    total += A[T % 2][i];
  }
  fprintf(stderr, "|sum: %e\t", total);
  for (i = 1; i < N + 1; i++) {
    sum_err_sqr += (A[T % 2][i] - (total / N)) * (A[T % 2][i] - (total / N));
  }
  fprintf(stderr, "|rms(A) = %7.2f\t", sqrt(sum_err_sqr));
  for (i = 1; i < N + 1; i++) {
    chtotal += ((char *)A[T % 2])[i];
  }
  fprintf(stderr, "|sum(rep(A)) = %d\n", chtotal);
#endif
  return 0;
}

// icc -O3 -fp-model precise heat_1d_np.c -o op-heat-1d-np -lm
// /* @ begin PrimeTile (num_tiling_levels=1; first_depth=1; last_depth=-1;
// boundary_tiling_level=-1;) @*/
// /* @ begin PrimeRegTile (scalar_replacement=0; T1t3=8; T1t4=8; ) @*/
// /* @ end @*/
// ,t2,t3,t4,t5,t6)
