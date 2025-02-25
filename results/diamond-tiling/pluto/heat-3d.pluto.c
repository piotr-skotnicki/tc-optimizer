#include <omp.h>
#include <math.h>
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

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

  int t1, t2, t3, t4, t5, t6, t7, t8;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
for (t1=-1;t1<=24;t1++) {
  lbp=ceild(t1,2);
  ubp=floord(t1+38,2);
#pragma omp parallel for private(lbv,ubv,t3,t4,t5,t6,t7,t8)
  for (t2=lbp;t2<=ubp;t2++) {
    for (t3=max(0,ceild(t1-1,2));t3<=floord(t1+39,2);t3++) {
      for (t4=max(max(0,ceild(t1-1,2)),t3-19);t4<=min(floord(t1+39,2),t3+19);t4++) {
        for (t5=max(max(max(max(max(0,4*t1),8*t2-150),8*t3-150),8*t4-150),8*t1-8*t2+1);t5<=min(min(min(min(min(98,4*t1+7),8*t2+6),8*t3+6),8*t4+6),8*t1-8*t2+157);t5++) {
          for (t6=max(max(8*t2,t5+1),-8*t1+8*t2+2*t5-7);t6<=min(min(8*t2+7,t5+150),-8*t1+8*t2+2*t5);t6++) {
            for (t7=max(8*t3,t5+1);t7<=min(8*t3+7,t5+150);t7++) {
              lbv=max(8*t4,t5+1);
              ubv=min(8*t4+7,t5+150);
#pragma ivdep
#pragma vector always
              for (t8=lbv;t8<=ubv;t8++) {
                A[(t5 + 1) % 2][(-t5+t6)][(-t5+t7)][(-t5+t8)] = ((((0.125 * ((A[t5 % 2][(-t5+t6) + 1][(-t5+t7)][(-t5+t8)] - (2.0 * A[t5 % 2][(-t5+t6)][(-t5+t7)][(-t5+t8)])) + A[t5 % 2][(-t5+t6) - 1][(-t5+t7)][(-t5+t8)])) + (0.125 * ((A[t5 % 2][(-t5+t6)][(-t5+t7) + 1][(-t5+t8)] - (2.0 * A[t5 % 2][(-t5+t6)][(-t5+t7)][(-t5+t8)])) + A[t5 % 2][(-t5+t6)][(-t5+t7) - 1][(-t5+t8)]))) + (0.125 * ((A[t5 % 2][(-t5+t6)][(-t5+t7)][(-t5+t8) - 1] - (2.0 * A[t5 % 2][(-t5+t6)][(-t5+t7)][(-t5+t8)])) + A[t5 % 2][(-t5+t6)][(-t5+t7)][(-t5+t8) + 1]))) + A[t5 % 2][(-t5+t6)][(-t5+t7)][(-t5+t8)]);;
              }
            }
          }
        }
      }
    }
  }
}

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
