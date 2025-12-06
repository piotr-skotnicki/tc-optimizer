#include <omp.h>
#include <math.h>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

/*
 * Discretized 1D heat equation stencil with non periodic boundary conditions
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
double A[2][N];
double total = 0;
double sum_err_sqr = 0;
long int chtotal = 0;
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
  const int BASE = 1024;

  // for timekeeping
  int ts_return = -1;
  struct timeval start, end, result;
  double tdiff = 0.0;
  long count = 0;

  fprintf(stderr, "Number of points = %ld\t|Number of timesteps = %ld\t", N, T);

  /* Initialization */
  srand(42); // seed with a constant value to verify results

  for (int i = 0; i < N; i++) {
    A[0][i] = 1.0 * (rand() % BASE);
  }

#ifdef TIME
  gettimeofday(&start, 0);
#endif

  int _N = N - 1;
  int _T = T;

  int t1, t2, t3, t4, t5;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
if ((_N >= 0) && (_T >= 1)) {
  for (t1=-1;t1<=floord(_T-1,64);t1++) {
    lbp=ceild(t1,2);
    ubp=min(min(floord(2*_T+_N-2,256),floord(128*t1+_N+127,256)),floord(256*t1+_N+254,256));
#pragma omp parallel for private(lbv,ubv,t3,t4,t5)
    for (t2=lbp;t2<=ubp;t2++) {
      if ((_N >= 2) && (t1 <= floord(256*t2-_N,128)) && (t2 >= ceild(_N,256))) {
        if (_N%2 == 0) {
          A[(((256*t2-_N)/2) + 1) % 2][(_N/2)] = (0.125 * ((A[((256*t2-_N)/2) % 2][(_N/2) == _N ? 0 : (_N/2) + 1] - (2.0 * A[((256*t2-_N)/2) % 2][(_N/2)])) + A[((256*t2-_N)/2) % 2][(_N/2) == 0 ? _N : (_N/2) - 1]));;
        }
      }
      if ((_N == 0) && (t1 == 2*t2)) {
        for (t3=64*t1;t3<=min(_T-1,64*t1+127);t3++) {
          if (t1%2 == 0) {
            A[(t3 + 1) % 2][0] = (0.125 * ((A[t3 % 2][0 == _N ? 0 : 0 + 1] - (2.0 * A[t3 % 2][0])) + A[t3 % 2][0 == 0 ? _N : 0 - 1]));;
          }
        }
      }
      if (_N >= 1) {
        for (t3=max(max(0,ceild(256*t2-_N+1,2)),64*t1);t3<=min(min(floord(256*t1-256*t2+_N+253,2),_T-1),64*t1+127);t3++) {
          lbv=max(max(128*t2,t3),-128*t1+128*t2+2*t3-127);
          ubv=min(min(floord(2*t3+_N-1,2),128*t2+127),-128*t1+128*t2+2*t3);
#pragma ivdep
#pragma vector always
          for (t5=lbv;t5<=ubv;t5++) {
            A[(t3 + 1) % 2][(-t3+t5)] = (0.125 * ((A[t3 % 2][(-t3+t5) == _N ? 0 : (-t3+t5) + 1] - (2.0 * A[t3 % 2][(-t3+t5)])) + A[t3 % 2][(-t3+t5) == 0 ? _N : (-t3+t5) - 1]));;
          }
          lbv=max(max(128*t2,t3),-128*t1+128*t2+2*t3-127);
          ubv=min(min(floord(2*t3+_N,2),128*t2+127),-128*t1+128*t2+2*t3);
#pragma ivdep
#pragma vector always
          for (t5=lbv;t5<=ubv;t5++) {
            A[(t3 + 1) % 2][(t3-t5+_N)] = (0.125 * ((A[t3 % 2][(t3-t5+_N) == _N ? 0 : (t3-t5+_N) + 1] - (2.0 * A[t3 % 2][(t3-t5+_N)])) + A[t3 % 2][(t3-t5+_N) == 0 ? _N : (t3-t5+_N) - 1]));;
          }
        }
      }
      if ((_N >= 2) && (t1 <= min(floord(256*t2-_N,128),floord(256*t2+2*_T-_N-256,256)))) {
        if (_N%2 == 0) {
          A[(((256*t1-256*t2+_N+254)/2) + 1) % 2][(_N/2)] = (0.125 * ((A[((256*t1-256*t2+_N+254)/2) % 2][(_N/2) == _N ? 0 : (_N/2) + 1] - (2.0 * A[((256*t1-256*t2+_N+254)/2) % 2][(_N/2)])) + A[((256*t1-256*t2+_N+254)/2) % 2][(_N/2) == 0 ? _N : (_N/2) - 1]));;
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
         ((((double)NUM_FP_OPS * N * T) / tdiff) / 1000000L));
#endif

  if (fopen(".test", "r")) {
    total = 0;
    for (int i = 0; i < N; i++) {
      total += A[T % 2][i];
    }
    fprintf(stderr, "|sum: %e\t", total);
    for (int i = 0; i < N; i++) {
      sum_err_sqr += (A[T % 2][i] - (total / N)) * (A[T % 2][i] - (total / N));
    }
    fprintf(stderr, "|rms(A) = %7.2f\t", sqrt(sum_err_sqr));
    for (int i = 0; i < N; i++) {
      chtotal += ((char *)A[T % 2])[i];
    }
    fprintf(stderr, "|sum(rep(A)) = %ld\n", chtotal);
  }
  return 0;
}

// icc -O3 -fp-model precise heat_1d_np.c -o op-heat-1d-np -lm
// /* @ begin PrimeTile (num_tiling_levels=1; first_depth=1; last_depth=-1;
// boundary_tiling_level=-1;) @*/
// /* @ begin PrimeRegTile (scalar_replacement=0; T1t3=8; T1t4=8; ) @*/
// /* @ end @*/
// ,t2,t3,t4,t5,t6)