/*
 * Calculating the price of American Put Option
 * Adapted from Pochoir test bench
 *
 * Irshad Pananilath: irshad@csa.iisc.ernet.in
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>

#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

#define SIMPLE 0
#define N_RANK 1

/* apop_pochoir -S 100 -E 95 -r 10 -V 30 -T 1 -s 2000000 -t 10000 */

/* #define DEFAULT_S 100.00 */
/* #define DEFAULT_E 100.00 */
/* #define DEFAULT_r   0.10 */
/* #define DEFAULT_V   0.25 */
/* #define DEFAULT_T   1.00 */
/*  */
/* #define DEFAULT_s 100 */
/* #define DEFAULT_t 100 */

#define DEFAULT_S 100.00
#define DEFAULT_E 95.00
#define DEFAULT_r   0.10
#define DEFAULT_V   0.30
#define DEFAULT_T   1.00

#define DEFAULT_s 2000000
#define DEFAULT_t 10000

double C[3][DEFAULT_s + 1]; // TODO::change the space dimension
double F[2][DEFAULT_s + 1]; // TODO::change the space dimension

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

void print_usage( char *prog )
{
  printf( "Usage: %s [ options ]\n\n", prog );

  printf( "Options:\n" );

  printf( "\t-S value : spot price ( default: %0.2lf )\n", DEFAULT_S );
  printf( "\t-E value : exercise price ( default: %0.2lf )\n", DEFAULT_E );
  printf( "\t-r value : interest rate ( default: %0.2lf )\n", DEFAULT_r * 100 );
  printf( "\t-V value : volatility ( default: %0.2lf )\n", DEFAULT_V * 100 );
  printf( "\t-T value : time to mature in years ( default: %0.2lf )\n\n", DEFAULT_T );

  printf( "\t-s value : steps in space dimension ( default: %d )\n", DEFAULT_s );
  printf( "\t-t value : steps in time dimension ( default: %d )\n\n", DEFAULT_t );

  printf( "\t-i               : Run iterative stencil\n\n" );

  printf( "\t-h               : print this help screen\n\n" );
}

void computeCoeffs( double r, double V, double T, int ns, int nt) {
  double V2 = V * V;
  double dt = T / nt;
  double r1 = 1.0 / ( 1.0 + r * dt );
  double r2 = dt / ( 1.0 + r * dt );
  int x;

  for ( x = 0; x <= ns; ++x ) {
    C[0][x] = r2 * 0.5 * x * ( - r + V2 * x );
    C[1][x] = r1 * ( 1 - V2 * x * x * dt );
    C[2][x] = r2 * 0.5 * x * ( r + V2 * x );
  }
}

int main( int argc, char *argv[ ] ) {
  // for timekeeping
  int ts_return = -1;
  struct timeval start, end, result;
  double tdiff = 0.0;

  double S, E, r, V, T;
  int ns, nt;
  int i, t;
  double v, price1;

  S = DEFAULT_S;
  E = DEFAULT_E;
  r = DEFAULT_r;
  V = DEFAULT_V;
  T = DEFAULT_T;

  ns = DEFAULT_s;
  nt = DEFAULT_t;

#ifdef DEBUG
  printf( "\nStencil-based DP for the price of American put option ( Run with option -h for help ).\n\n" );
  printf( "Parameters:\n\n" );

  printf( "\t spot price = %0.2lf\n", S );
  printf( "\t exercise price = %0.2lf\n", E );
  printf( "\t interest rate = %0.2lf\%\n", r * 100 );
  printf( "\t volatility = %0.2lf\%\n", V * 100 );
  printf( "\t time to mature ( in years ) = %0.2lf\n\n", T );

  printf( "\t steps in space dimension = %d\n", ns );
  printf( "\t steps in time dimension = %d\n\n", nt );
#endif

  ns = ns + ( ns & 1 );
  double dS = 2.0 * S / ns;

  /* computeCoeffs( r, V, T, ns, nt ); */

  // initialize
  double V2 = V * V;
  double dt = T / nt;
  double r1 = 1.0 / ( 1.0 + r * dt );
  double r2 = dt / ( 1.0 + r * dt );
  int x;

  for ( x = 0; x <= ns; ++x ) {
    C[0][x] = r2 * 0.5 * x * ( - r + V2 * x );
    C[1][x] = r1 * ( 1 - V2 * x * x * dt );
    C[2][x] = r2 * 0.5 * x * ( r + V2 * x );
  }

  for ( i = 0; i <= ns; ++i )
    F[0][i] = max( 0.0, E - i * dS );

  F[1][0] = E;

  for (t = 0; t < nt; ++t ) {
    F[0][ns] = 0;
    F[1][ns] = 0;
  }

#ifdef TIME
  gettimeofday(&start, 0);
#endif

/* TC Optimizing Compiler 0.4.1 */
/* ./tc ../examples/pluto/apop.scop.c --diamond-tiling --omp-for-codegen --iterative-tc --inline --debug -b 128 --drop-bounds */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#pragma scop
if (nt + ns >= 257) {
  if (ns <= 128) {
    for (int k = 1; k <= 2; k += 1) {
      if (k == 2) {
        for (int i0 = 1; i0 <= min(253, nt - 1); i0 += 1) {
          for (int i1 = 1; i1 <= min(min(ns - 1, i0), -i0 + 254); i1 += 1) {
            F[(i0 + 1) % 2][i1] = max(((C[0][i1] * F[i0 % 2][i1 - 1]) + (C[1][i1] * F[i0 % 2][i1])) + (C[2][i1] * F[i0 % 2][i1 + 1]), E - ((i1) * dS));
          }
        }
      } else {
        for (int i0 = 0; i0 < ns - 1; i0 += 1) {
          for (int i1 = i0 + 1; i1 < ns; i1 += 1) {
            F[(i0 + 1) % 2][i1] = max(((C[0][i1] * F[i0 % 2][i1 - 1]) + (C[1][i1] * F[i0 % 2][i1])) + (C[2][i1] * F[i0 % 2][i1 + 1]), E - ((i1) * dS));
          }
        }
      }
    }
    if (ns <= 64) {
      for (int ii0 = 1; ii0 < floord(nt + ns - 1, 128); ii0 += 1) {
        for (int i0 = -ns + 128 * ii0 + 128; i0 <= min(nt - 1, 128 * ii0 + 253); i0 += 1) {
          for (int i1 = max(1, 128 * ii0 - i0 + 127); i1 <= min(ns - 1, 128 * ii0 - i0 + 254); i1 += 1) {
            F[(i0 + 1) % 2][i1] = max(((C[0][i1] * F[i0 % 2][i1 - 1]) + (C[1][i1] * F[i0 % 2][i1])) + (C[2][i1] * F[i0 % 2][i1 + 1]), E - ((i1) * dS));
          }
        }
      }
    }
  } else if (ns <= 192) {
    for (int k = max(0, -ns + 130); k <= 2; k += 1) {
      if (k >= 1) {
        if (k == 1) {
          for (int i0 = 0; i0 <= min(126, nt - 1); i0 += 1) {
            for (int i1 = i0 + 1; i1 <= min(min(ns - 1, i0 + 128), -i0 + 254); i1 += 1) {
              F[(i0 + 1) % 2][i1] = max(((C[0][i1] * F[i0 % 2][i1 - 1]) + (C[1][i1] * F[i0 % 2][i1])) + (C[2][i1] * F[i0 % 2][i1 + 1]), E - ((i1) * dS));
            }
          }
        } else {
          for (int i0 = 1; i0 <= min(253, nt - 1); i0 += 1) {
            for (int i1 = 1; i1 <= min(i0, -i0 + 254); i1 += 1) {
              F[(i0 + 1) % 2][i1] = max(((C[0][i1] * F[i0 % 2][i1 - 1]) + (C[1][i1] * F[i0 % 2][i1])) + (C[2][i1] * F[i0 % 2][i1 + 1]), E - ((i1) * dS));
            }
          }
        }
      } else {
        for (int i0 = 0; i0 < ns - 129; i0 += 1) {
          for (int i1 = i0 + 129; i1 < ns; i1 += 1) {
            F[(i0 + 1) % 2][i1] = max(((C[0][i1] * F[i0 % 2][i1 - 1]) + (C[1][i1] * F[i0 % 2][i1])) + (C[2][i1] * F[i0 % 2][i1 + 1]), E - ((i1) * dS));
          }
        }
      }
      if (k == 2) {
        for (int i0 = -ns + 256; i0 < min(nt, ns - 1); i0 += 1) {
          for (int i1 = max(i0 + 1, -i0 + 255); i1 < ns; i1 += 1) {
            F[(i0 + 1) % 2][i1] = max(((C[0][i1] * F[i0 % 2][i1 - 1]) + (C[1][i1] * F[i0 % 2][i1])) + (C[2][i1] * F[i0 % 2][i1 + 1]), E - ((i1) * dS));
          }
        }
      }
    }
  }
  if (ns >= 65) {
    for (int ii0 = max(0, -((ns + 63) / 128) + 2); ii0 <= min(floord(nt - 1, 128), (nt + ns - 1) / 128 - 1); ii0 += 1) {
      if (ii0 == 0) {
        #pragma omp parallel for
        for (int ii1 = 0; ii1 < (ns - 2) / 128; ii1 += 1) {
          for (int i0 = 0; i0 <= min(min(62, nt - 1), ns - 128 * ii1 - 130); i0 += 1) {
            for (int i1 = 128 * ii1 + i0 + 129; i1 <= min(ns - 1, 128 * ii1 - i0 + 254); i1 += 1) {
              F[(i0 + 1) % 2][i1] = max(((C[0][i1] * F[i0 % 2][i1 - 1]) + (C[1][i1] * F[i0 % 2][i1])) + (C[2][i1] * F[i0 % 2][i1 + 1]), E - ((i1) * dS));
            }
          }
        }
      }
      for (int k = 1; k <= min(2, -2 * ii0 + (nt - 1) / 64 + 1); k += 1) {
        if (ii0 == 0 && k == 1) {
          for (int i0 = 0; i0 <= min(126, nt - 1); i0 += 1) {
            for (int i1 = i0 + 1; i1 <= min(i0 + 128, -i0 + 254); i1 += 1) {
              F[(i0 + 1) % 2][i1] = max(((C[0][i1] * F[i0 % 2][i1 - 1]) + (C[1][i1] * F[i0 % 2][i1])) + (C[2][i1] * F[i0 % 2][i1 + 1]), E - ((i1) * dS));
            }
          }
        } else if (ii0 == 0 && k == 2) {
          for (int i0 = 1; i0 <= min(253, nt - 1); i0 += 1) {
            for (int i1 = 1; i1 <= min(i0, -i0 + 254); i1 += 1) {
              F[(i0 + 1) % 2][i1] = max(((C[0][i1] * F[i0 % 2][i1 - 1]) + (C[1][i1] * F[i0 % 2][i1])) + (C[2][i1] * F[i0 % 2][i1 + 1]), E - ((i1) * dS));
            }
          }
        }
        #pragma omp parallel for
        for (int ii1 = max(0, -ii0 + 1); ii1 < min(-ii0 + (nt + ns - 1) / 128, (ns + 64 * k - 1) / 128); ii1 += 1) {
          if (k == 1) {
            for (int i0 = max(128 * ii0, -ns + 128 * ii0 + 128 * ii1 + 128); i0 <= min(min(nt - 1, 128 * ii0 + 126), ns + 128 * ii0 - 128 * ii1 - 2); i0 += 1) {
              for (int i1 = max(-128 * ii0 + 128 * ii1 + i0 + 1, 128 * ii0 + 128 * ii1 - i0 + 127); i1 <= min(min(ns - 1, -128 * ii0 + 128 * ii1 + i0 + 128), 128 * ii0 + 128 * ii1 - i0 + 254); i1 += 1) {
                F[(i0 + 1) % 2][i1] = max(((C[0][i1] * F[i0 % 2][i1 - 1]) + (C[1][i1] * F[i0 % 2][i1])) + (C[2][i1] * F[i0 % 2][i1 + 1]), E - ((i1) * dS));
              }
            }
          } else if (ii1 >= 1) {
            for (int i0 = max(128 * ii0 + 64, -ns + 128 * ii0 + 128 * ii1 + 128); i0 <= min(min(nt - 1, 128 * ii0 + 190), ns + 128 * ii0 - 128 * ii1 + 126); i0 += 1) {
              for (int i1 = max(-128 * ii0 + 128 * ii1 + i0 - 127, 128 * ii0 + 128 * ii1 - i0 + 127); i1 <= min(min(ns - 1, -128 * ii0 + 128 * ii1 + i0), 128 * ii0 + 128 * ii1 - i0 + 254); i1 += 1) {
                F[(i0 + 1) % 2][i1] = max(((C[0][i1] * F[i0 % 2][i1 - 1]) + (C[1][i1] * F[i0 % 2][i1])) + (C[2][i1] * F[i0 % 2][i1 + 1]), E - ((i1) * dS));
              }
            }
          } else {
            for (int i0 = 128 * ii0 + 64; i0 <= min(nt - 1, 128 * ii0 + 253); i0 += 1) {
              for (int i1 = max(1, 128 * ii0 - i0 + 127); i1 <= min(min(ns - 1, -128 * ii0 + i0), 128 * ii0 - i0 + 254); i1 += 1) {
                F[(i0 + 1) % 2][i1] = max(((C[0][i1] * F[i0 % 2][i1 - 1]) + (C[1][i1] * F[i0 % 2][i1])) + (C[2][i1] * F[i0 % 2][i1 + 1]), E - ((i1) * dS));
              }
            }
          }
        }
      }
      if (nt <= 64 && ii0 == 0) {
        for (int i0 = 1; i0 < nt; i0 += 1) {
          for (int i1 = 1; i1 <= i0; i1 += 1) {
            F[(i0 + 1) % 2][i1] = max(((C[0][i1] * F[i0 % 2][i1 - 1]) + (C[1][i1] * F[i0 % 2][i1])) + (C[2][i1] * F[i0 % 2][i1 + 1]), E - ((i1) * dS));
          }
        }
      }
    }
  }
} else {
  for (int i0 = 0; i0 < min(nt, ns - 129); i0 += 1) {
    for (int i1 = i0 + 129; i1 < ns; i1 += 1) {
      F[(i0 + 1) % 2][i1] = max(((C[0][i1] * F[i0 % 2][i1 - 1]) + (C[1][i1] * F[i0 % 2][i1])) + (C[2][i1] * F[i0 % 2][i1 + 1]), E - ((i1) * dS));
    }
  }
  for (int k = 1; k <= min(2, nt); k += 1) {
    if (k == 2) {
      for (int i0 = 1; i0 < nt; i0 += 1) {
        for (int i1 = 1; i1 <= min(ns - 1, i0); i1 += 1) {
          F[(i0 + 1) % 2][i1] = max(((C[0][i1] * F[i0 % 2][i1 - 1]) + (C[1][i1] * F[i0 % 2][i1])) + (C[2][i1] * F[i0 % 2][i1 + 1]), E - ((i1) * dS));
        }
      }
    } else {
      for (int i0 = 0; i0 < min(nt, ns - 1); i0 += 1) {
        for (int i1 = i0 + 1; i1 <= min(ns - 1, i0 + 128); i1 += 1) {
          F[(i0 + 1) % 2][i1] = max(((C[0][i1] * F[i0 % 2][i1 - 1]) + (C[1][i1] * F[i0 % 2][i1])) + (C[2][i1] * F[i0 % 2][i1 + 1]), E - ((i1) * dS));
        }
      }
    }
  }
}
#pragma endscop

  price1 = F[nt%2][( ns >> 1 )];

#ifdef TIME
  gettimeofday(&end, 0);

  ts_return = timeval_subtract(&result, &end, &start);
  tdiff = (double) (result.tv_sec + result.tv_usec * 1.0e-6);

  printf( "\n\nIterative Stencil:\n" );
  fprintf(stderr, "\t option price = %.2lf\n", price1);

  printf("Time taken: %7.5lfms\n", tdiff * 1.0e3);
#endif

  return 0;
}
