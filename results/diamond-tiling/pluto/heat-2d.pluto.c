#include <omp.h>
#include <math.h>
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

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
#define N 4000L
#define T 1000L
#endif

#define NUM_FP_OPS 10

/* Define our arrays */
double A[2][N + 2][N + 2];
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
  long int i, j;
  const int BASE = 1024;

  // for timekeeping
  int ts_return = -1;
  struct timeval start, end, result;
  double tdiff = 0.0;

  printf("Number of points = %ld\t|Number of timesteps = %ld\t", N * N, T);

  /* Initialization */
  srand(42); // seed with a constant value to verify results

  for (i = 1; i < N + 1; i++) {
    for (j = 1; j < N + 1; j++) {
      A[0][i][j] = 1.0 * (rand() % BASE);
    }
  }

#ifdef TIME
  gettimeofday(&start, 0);
#endif

  int t1, t2, t3, t4, t5, t6;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
for (t1=-1;t1<=124;t1++) {
  lbp=ceild(t1,2);
  ubp=floord(t1+500,2);
#pragma omp parallel for private(lbv,ubv,t3,t4,t5,t6)
  for (t2=lbp;t2<=ubp;t2++) {
    for (t3=max(0,ceild(t1-1,2));t3<=floord(t1+501,2);t3++) {
      for (t4=max(max(max(0,8*t1),16*t3-4000),16*t1-16*t2+1);t4<=min(min(min(999,8*t1+15),16*t2+14),16*t3+14);t4++) {
        for (t5=max(max(16*t2,t4+1),-16*t1+16*t2+2*t4-15);t5<=min(min(16*t2+15,t4+4000),-16*t1+16*t2+2*t4);t5++) {
          lbv=max(16*t3,t4+1);
          ubv=min(16*t3+15,t4+4000);
#pragma ivdep
#pragma vector always
          for (t6=lbv;t6<=ubv;t6++) {
            A[(t4 + 1) % 2][(-t4+t5)][(-t4+t6)] = (((0.125 * ((A[t4 % 2][(-t4+t5) + 1][(-t4+t6)] - (2.0 * A[t4 % 2][(-t4+t5)][(-t4+t6)])) + A[t4 % 2][(-t4+t5) - 1][(-t4+t6)])) + (0.125 * ((A[t4 % 2][(-t4+t5)][(-t4+t6) + 1] - (2.0 * A[t4 % 2][(-t4+t5)][(-t4+t6)])) + A[t4 % 2][(-t4+t5)][(-t4+t6) - 1]))) + A[t4 % 2][(-t4+t5)][(-t4+t6)]);;
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

  printf("|Time taken =  %7.5lfms\t", tdiff * 1.0e3);
  printf("|MFLOPS =  %f\n",
         ((((double)NUM_FP_OPS * N * N * T) / tdiff) / 1000000L));
#endif

  if (fopen(".test", "r")) {
    for (int i = 1; i < N + 1; i++) {
      for (int j = 1; j < N + 1; j++) {
        total += A[T % 2][i][j];
      }
    }
    fprintf(stderr, "|sum: %e\t", total);
    for (int i = 1; i < N + 1; i++) {
      for (int j = 1; j < N + 1; j++) {
        sum_err_sqr +=
            (A[T % 2][i][j] - (total / N)) * (A[T % 2][i][j] - (total / N));
      }
    }
    fprintf(stderr, "|rms(A) = %7.2f\t", sqrt(sum_err_sqr));
    for (int i = 1; i < N + 1; i++) {
      for (int j = 1; j < N + 1; j++) {
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
