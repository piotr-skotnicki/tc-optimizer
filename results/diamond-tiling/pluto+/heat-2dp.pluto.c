#include <omp.h>
#include <math.h>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
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

  int t1, t2, t3, t4, t5, t6, t7;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
if ((_N >= 0) && (_T >= 1)) {
  for (t1=-1;t1<=floord(_T-1,16);t1++) {
    lbp=ceild(t1,2);
    ubp=min(min(floord(2*_T+_N-2,64),floord(32*t1+_N+31,64)),floord(64*t1+_N+62,64));
#pragma omp parallel for private(lbv,ubv,t3,t4,t5,t6,t7)
    for (t2=lbp;t2<=ubp;t2++) {
      for (t3=max(max(0,ceild(t1-1,2)),ceild(64*t2-_N-62,64));t3<=min(floord(2*_T+_N-2,64),floord(32*t1+_N+62,64));t3++) {
        if ((t1 <= floord(64*t2-_N,32)) && (t2 >= max(ceild(_N,64),t3+1))) {
          if (_N%2 == 0) {
            lbv=max(ceild(64*t2-_N,2),32*t3);
            ubv=32*t3+31;
#pragma ivdep
#pragma vector always
            for (t7=lbv;t7<=ubv;t7++) {
              A[(((64*t2-_N)/2) + 1) % 2][(_N/2)][((-64*t2+2*t7+_N)/2)] = (((0.125 * ((A[((64*t2-_N)/2) % 2][(_N/2) == _N ? 0 : (_N/2) + 1][((-64*t2+2*t7+_N)/2)] - (2.0 * A[((64*t2-_N)/2) % 2][(_N/2)][((-64*t2+2*t7+_N)/2)])) + A[((64*t2-_N)/2) % 2][(_N/2) == 0 ? _N : (_N/2) - 1][((-64*t2+2*t7+_N)/2)])) + (0.125 * ((A[((64*t2-_N)/2) % 2][(_N/2)][((-64*t2+2*t7+_N)/2) == _N ? 0 : ((-64*t2+2*t7+_N)/2) + 1] - (2.0 * A[((64*t2-_N)/2) % 2][(_N/2)][((-64*t2+2*t7+_N)/2)])) + A[((64*t2-_N)/2) % 2][(_N/2)][((-64*t2+2*t7+_N)/2) == 0 ? _N : ((-64*t2+2*t7+_N)/2) - 1]))) + A[((64*t2-_N)/2) % 2][(_N/2)][((-64*t2+2*t7+_N)/2)]);;
            }
            lbv=max(ceild(64*t2-_N,2),32*t3);
            ubv=32*t3+31;
#pragma ivdep
#pragma vector always
            for (t7=lbv;t7<=ubv;t7++) {
              A[(((64*t2-_N)/2) + 1) % 2][(_N/2)][((64*t2-2*t7+_N)/2)] = (((0.125 * ((A[((64*t2-_N)/2) % 2][(_N/2) == _N ? 0 : (_N/2) + 1][((64*t2-2*t7+_N)/2)] - (2.0 * A[((64*t2-_N)/2) % 2][(_N/2)][((64*t2-2*t7+_N)/2)])) + A[((64*t2-_N)/2) % 2][(_N/2) == 0 ? _N : (_N/2) - 1][((64*t2-2*t7+_N)/2)])) + (0.125 * ((A[((64*t2-_N)/2) % 2][(_N/2)][((64*t2-2*t7+_N)/2) == _N ? 0 : ((64*t2-2*t7+_N)/2) + 1] - (2.0 * A[((64*t2-_N)/2) % 2][(_N/2)][((64*t2-2*t7+_N)/2)])) + A[((64*t2-_N)/2) % 2][(_N/2)][((64*t2-2*t7+_N)/2) == 0 ? _N : ((64*t2-2*t7+_N)/2) - 1]))) + A[((64*t2-_N)/2) % 2][(_N/2)][((64*t2-2*t7+_N)/2)]);;
            }
          }
        }
        if ((t1 <= floord(64*t3-_N,32)) && (t2 <= t3-1) && (t3 >= ceild(_N,64))) {
          if (_N%2 == 0) {
            for (t6=max(max(ceild(64*t3-_N,2),32*t2),-32*t1+32*t2+64*t3-_N-31);t6<=min(32*t2+31,-32*t1+32*t2+64*t3-_N);t6++) {
              A[(((64*t3-_N)/2) + 1) % 2][((-64*t3+2*t6+_N)/2)][(_N/2)] = (((0.125 * ((A[((64*t3-_N)/2) % 2][((-64*t3+2*t6+_N)/2) == _N ? 0 : ((-64*t3+2*t6+_N)/2) + 1][(_N/2)] - (2.0 * A[((64*t3-_N)/2) % 2][((-64*t3+2*t6+_N)/2)][(_N/2)])) + A[((64*t3-_N)/2) % 2][((-64*t3+2*t6+_N)/2) == 0 ? _N : ((-64*t3+2*t6+_N)/2) - 1][(_N/2)])) + (0.125 * ((A[((64*t3-_N)/2) % 2][((-64*t3+2*t6+_N)/2)][(_N/2) == _N ? 0 : (_N/2) + 1] - (2.0 * A[((64*t3-_N)/2) % 2][((-64*t3+2*t6+_N)/2)][(_N/2)])) + A[((64*t3-_N)/2) % 2][((-64*t3+2*t6+_N)/2)][(_N/2) == 0 ? _N : (_N/2) - 1]))) + A[((64*t3-_N)/2) % 2][((-64*t3+2*t6+_N)/2)][(_N/2)]);;
            }
            for (t6=max(max(ceild(64*t3-_N,2),32*t2),-32*t1+32*t2+64*t3-_N-31);t6<=min(32*t2+31,-32*t1+32*t2+64*t3-_N);t6++) {
              A[(((64*t3-_N)/2) + 1) % 2][((64*t3-2*t6+_N)/2)][(_N/2)] = (((0.125 * ((A[((64*t3-_N)/2) % 2][((64*t3-2*t6+_N)/2) == _N ? 0 : ((64*t3-2*t6+_N)/2) + 1][(_N/2)] - (2.0 * A[((64*t3-_N)/2) % 2][((64*t3-2*t6+_N)/2)][(_N/2)])) + A[((64*t3-_N)/2) % 2][((64*t3-2*t6+_N)/2) == 0 ? _N : ((64*t3-2*t6+_N)/2) - 1][(_N/2)])) + (0.125 * ((A[((64*t3-_N)/2) % 2][((64*t3-2*t6+_N)/2)][(_N/2) == _N ? 0 : (_N/2) + 1] - (2.0 * A[((64*t3-_N)/2) % 2][((64*t3-2*t6+_N)/2)][(_N/2)])) + A[((64*t3-_N)/2) % 2][((64*t3-2*t6+_N)/2)][(_N/2) == 0 ? _N : (_N/2) - 1]))) + A[((64*t3-_N)/2) % 2][((64*t3-2*t6+_N)/2)][(_N/2)]);;
            }
          }
        }
        if ((_N == 0) && (t1 == 2*t2) && (t1 == 2*t3)) {
          for (t4=16*t1;t4<=min(_T-1,16*t1+31);t4++) {
            if (t1%2 == 0) {
              A[(t4 + 1) % 2][0][0] = (((0.125 * ((A[t4 % 2][0 == _N ? 0 : 0 + 1][0] - (2.0 * A[t4 % 2][0][0])) + A[t4 % 2][0 == 0 ? _N : 0 - 1][0])) + (0.125 * ((A[t4 % 2][0][0 == _N ? 0 : 0 + 1] - (2.0 * A[t4 % 2][0][0])) + A[t4 % 2][0][0 == 0 ? _N : 0 - 1]))) + A[t4 % 2][0][0]);;
            }
          }
        }
        if ((_N >= 2) && (t1 <= floord(64*t2-_N,32)) && (t2 == t3) && (t2 >= ceild(_N,64))) {
          if (_N%2 == 0) {
            A[(((64*t2-_N)/2) + 1) % 2][(_N/2)][(_N/2)] = (((0.125 * ((A[((64*t2-_N)/2) % 2][(_N/2) == _N ? 0 : (_N/2) + 1][(_N/2)] - (2.0 * A[((64*t2-_N)/2) % 2][(_N/2)][(_N/2)])) + A[((64*t2-_N)/2) % 2][(_N/2) == 0 ? _N : (_N/2) - 1][(_N/2)])) + (0.125 * ((A[((64*t2-_N)/2) % 2][(_N/2)][(_N/2) == _N ? 0 : (_N/2) + 1] - (2.0 * A[((64*t2-_N)/2) % 2][(_N/2)][(_N/2)])) + A[((64*t2-_N)/2) % 2][(_N/2)][(_N/2) == 0 ? _N : (_N/2) - 1]))) + A[((64*t2-_N)/2) % 2][(_N/2)][(_N/2)]);;
          }
        }
        if (_N >= 1) {
          for (t4=max(max(max(0,ceild(64*t2-_N+1,2)),ceild(64*t3-_N+1,2)),16*t1);t4<=min(min(min(floord(64*t1-64*t2+_N+61,2),_T-1),16*t1+31),32*t3+31);t4++) {
            for (t6=max(max(32*t2,t4),-32*t1+32*t2+2*t4-31);t6<=min(min(floord(2*t4+_N-1,2),32*t2+31),-32*t1+32*t2+2*t4);t6++) {
              lbv=max(32*t3,t4);
              ubv=min(floord(2*t4+_N-1,2),32*t3+31);
#pragma ivdep
#pragma vector always
              for (t7=lbv;t7<=ubv;t7++) {
                A[(t4 + 1) % 2][(-t4+t6)][(-t4+t7)] = (((0.125 * ((A[t4 % 2][(-t4+t6) == _N ? 0 : (-t4+t6) + 1][(-t4+t7)] - (2.0 * A[t4 % 2][(-t4+t6)][(-t4+t7)])) + A[t4 % 2][(-t4+t6) == 0 ? _N : (-t4+t6) - 1][(-t4+t7)])) + (0.125 * ((A[t4 % 2][(-t4+t6)][(-t4+t7) == _N ? 0 : (-t4+t7) + 1] - (2.0 * A[t4 % 2][(-t4+t6)][(-t4+t7)])) + A[t4 % 2][(-t4+t6)][(-t4+t7) == 0 ? _N : (-t4+t7) - 1]))) + A[t4 % 2][(-t4+t6)][(-t4+t7)]);;
              }
            }
            for (t6=max(max(32*t2,t4),-32*t1+32*t2+2*t4-31);t6<=min(min(floord(2*t4+_N,2),32*t2+31),-32*t1+32*t2+2*t4);t6++) {
              lbv=max(32*t3,t4);
              ubv=min(floord(2*t4+_N-1,2),32*t3+31);
#pragma ivdep
#pragma vector always
              for (t7=lbv;t7<=ubv;t7++) {
                A[(t4 + 1) % 2][(t4-t6+_N)][(-t4+t7)] = (((0.125 * ((A[t4 % 2][(t4-t6+_N) == _N ? 0 : (t4-t6+_N) + 1][(-t4+t7)] - (2.0 * A[t4 % 2][(t4-t6+_N)][(-t4+t7)])) + A[t4 % 2][(t4-t6+_N) == 0 ? _N : (t4-t6+_N) - 1][(-t4+t7)])) + (0.125 * ((A[t4 % 2][(t4-t6+_N)][(-t4+t7) == _N ? 0 : (-t4+t7) + 1] - (2.0 * A[t4 % 2][(t4-t6+_N)][(-t4+t7)])) + A[t4 % 2][(t4-t6+_N)][(-t4+t7) == 0 ? _N : (-t4+t7) - 1]))) + A[t4 % 2][(t4-t6+_N)][(-t4+t7)]);;
              }
            }
            for (t6=max(max(32*t2,t4),-32*t1+32*t2+2*t4-31);t6<=min(min(floord(2*t4+_N-1,2),32*t2+31),-32*t1+32*t2+2*t4);t6++) {
              lbv=max(32*t3,t4);
              ubv=min(floord(2*t4+_N,2),32*t3+31);
#pragma ivdep
#pragma vector always
              for (t7=lbv;t7<=ubv;t7++) {
                A[(t4 + 1) % 2][(-t4+t6)][(t4-t7+_N)] = (((0.125 * ((A[t4 % 2][(-t4+t6) == _N ? 0 : (-t4+t6) + 1][(t4-t7+_N)] - (2.0 * A[t4 % 2][(-t4+t6)][(t4-t7+_N)])) + A[t4 % 2][(-t4+t6) == 0 ? _N : (-t4+t6) - 1][(t4-t7+_N)])) + (0.125 * ((A[t4 % 2][(-t4+t6)][(t4-t7+_N) == _N ? 0 : (t4-t7+_N) + 1] - (2.0 * A[t4 % 2][(-t4+t6)][(t4-t7+_N)])) + A[t4 % 2][(-t4+t6)][(t4-t7+_N) == 0 ? _N : (t4-t7+_N) - 1]))) + A[t4 % 2][(-t4+t6)][(t4-t7+_N)]);;
              }
            }
            for (t6=max(max(32*t2,t4),-32*t1+32*t2+2*t4-31);t6<=min(min(floord(2*t4+_N,2),32*t2+31),-32*t1+32*t2+2*t4);t6++) {
              lbv=max(32*t3,t4);
              ubv=min(floord(2*t4+_N,2),32*t3+31);
#pragma ivdep
#pragma vector always
              for (t7=lbv;t7<=ubv;t7++) {
                A[(t4 + 1) % 2][(t4-t6+_N)][(t4-t7+_N)] = (((0.125 * ((A[t4 % 2][(t4-t6+_N) == _N ? 0 : (t4-t6+_N) + 1][(t4-t7+_N)] - (2.0 * A[t4 % 2][(t4-t6+_N)][(t4-t7+_N)])) + A[t4 % 2][(t4-t6+_N) == 0 ? _N : (t4-t6+_N) - 1][(t4-t7+_N)])) + (0.125 * ((A[t4 % 2][(t4-t6+_N)][(t4-t7+_N) == _N ? 0 : (t4-t7+_N) + 1] - (2.0 * A[t4 % 2][(t4-t6+_N)][(t4-t7+_N)])) + A[t4 % 2][(t4-t6+_N)][(t4-t7+_N) == 0 ? _N : (t4-t7+_N) - 1]))) + A[t4 % 2][(t4-t6+_N)][(t4-t7+_N)]);;
              }
            }
          }
        }
        if ((_N >= 2) && (t1 <= min(min(floord(64*t2-_N,32),floord(64*t2+64*t3-_N,64)),floord(64*t2+2*_T-_N-64,64)))) {
          if (_N%2 == 0) {
            lbv=max(ceild(64*t1-64*t2+_N+62,2),32*t3);
            ubv=min(32*t3+31,32*t1-32*t2+_N+30);
#pragma ivdep
#pragma vector always
            for (t7=lbv;t7<=ubv;t7++) {
              A[(((64*t1-64*t2+_N+62)/2) + 1) % 2][(_N/2)][((-64*t1+64*t2+2*t7-_N-62)/2)] = (((0.125 * ((A[((64*t1-64*t2+_N+62)/2) % 2][(_N/2) == _N ? 0 : (_N/2) + 1][((-64*t1+64*t2+2*t7-_N-62)/2)] - (2.0 * A[((64*t1-64*t2+_N+62)/2) % 2][(_N/2)][((-64*t1+64*t2+2*t7-_N-62)/2)])) + A[((64*t1-64*t2+_N+62)/2) % 2][(_N/2) == 0 ? _N : (_N/2) - 1][((-64*t1+64*t2+2*t7-_N-62)/2)])) + (0.125 * ((A[((64*t1-64*t2+_N+62)/2) % 2][(_N/2)][((-64*t1+64*t2+2*t7-_N-62)/2) == _N ? 0 : ((-64*t1+64*t2+2*t7-_N-62)/2) + 1] - (2.0 * A[((64*t1-64*t2+_N+62)/2) % 2][(_N/2)][((-64*t1+64*t2+2*t7-_N-62)/2)])) + A[((64*t1-64*t2+_N+62)/2) % 2][(_N/2)][((-64*t1+64*t2+2*t7-_N-62)/2) == 0 ? _N : ((-64*t1+64*t2+2*t7-_N-62)/2) - 1]))) + A[((64*t1-64*t2+_N+62)/2) % 2][(_N/2)][((-64*t1+64*t2+2*t7-_N-62)/2)]);;
            }
            lbv=max(ceild(64*t1-64*t2+_N+62,2),32*t3);
            ubv=min(32*t3+31,32*t1-32*t2+_N+31);
#pragma ivdep
#pragma vector always
            for (t7=lbv;t7<=ubv;t7++) {
              A[(((64*t1-64*t2+_N+62)/2) + 1) % 2][(_N/2)][((64*t1-64*t2-2*t7+3*_N+62)/2)] = (((0.125 * ((A[((64*t1-64*t2+_N+62)/2) % 2][(_N/2) == _N ? 0 : (_N/2) + 1][((64*t1-64*t2-2*t7+3*_N+62)/2)] - (2.0 * A[((64*t1-64*t2+_N+62)/2) % 2][(_N/2)][((64*t1-64*t2-2*t7+3*_N+62)/2)])) + A[((64*t1-64*t2+_N+62)/2) % 2][(_N/2) == 0 ? _N : (_N/2) - 1][((64*t1-64*t2-2*t7+3*_N+62)/2)])) + (0.125 * ((A[((64*t1-64*t2+_N+62)/2) % 2][(_N/2)][((64*t1-64*t2-2*t7+3*_N+62)/2) == _N ? 0 : ((64*t1-64*t2-2*t7+3*_N+62)/2) + 1] - (2.0 * A[((64*t1-64*t2+_N+62)/2) % 2][(_N/2)][((64*t1-64*t2-2*t7+3*_N+62)/2)])) + A[((64*t1-64*t2+_N+62)/2) % 2][(_N/2)][((64*t1-64*t2-2*t7+3*_N+62)/2) == 0 ? _N : ((64*t1-64*t2-2*t7+3*_N+62)/2) - 1]))) + A[((64*t1-64*t2+_N+62)/2) % 2][(_N/2)][((64*t1-64*t2-2*t7+3*_N+62)/2)]);;
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
