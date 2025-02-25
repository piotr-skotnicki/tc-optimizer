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

  int t1, t2, t3, t4;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
for (t1=-1;t1<=15;t1++) {
  lbp=ceild(t1,2);
  ubp=floord(t1+25000,2);
#pragma omp parallel for private(lbv,ubv,t3,t4)
  for (t2=lbp;t2<=ubp;t2++) {
    for (t3=max(max(0,64*t1),128*t1-128*t2+1);t3<=min(min(999,64*t1+127),128*t2+126);t3++) {
      lbv=max(max(128*t2,t3+1),-128*t1+128*t2+2*t3-127);
      ubv=min(min(128*t2+127,t3+1600000),-128*t1+128*t2+2*t3);
#pragma ivdep
#pragma vector always
      for (t4=lbv;t4<=ubv;t4++) {
        A[(t3 + 1) % 2][(-t3+t4)] = (0.250 * ((A[t3 % 2][(-t3+t4) + 1] - (2.0 * A[t3 % 2][(-t3+t4)])) + A[t3 % 2][(-t3+t4) - 1]));;
      }
    }
  }
}

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
