#include <omp.h>
#include <math.h>
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

/*
 * Order-1, 3D 7 point stencil
 * Adapted from Pochoir test bench
 *
 * Irshad Pananilath: irshad@csa.iisc.ernet.in
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define NUM_FP_OPS 8

#define MAX(a,b) ((a) > (b) ? a : b)
#define MIN(a,b) ((a) < (b) ? a : b)


#define N 256L
#define T 200L

double A[2][N+2][N+2][N+2];

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

int main(int argc, char *argv[])
{
  // for timekeeping
  int ts_return = -1;
  struct timeval start, end, result;
  double tdiff = 0.0;

  int t, i, j, k;

  const int BASE = 1024;
  const double alpha = 0.0876;
  const double beta = 0.0765;
  double total;

    printf("Number of points = %ld\t|Number of timesteps = %ld\t", N, T);

  // initialize variables
  //
  srand(42);
    for (i = 1; i < N+1; i++) {
        for (j = 1; j < N+1; j++) {
            for (k = 1; k < N+1; k++) {
                A[0][i][j][k] = 1.0 * (rand() % BASE);
            }
        }
    }


#ifdef TIME
  gettimeofday(&start, 0);
#endif

  // serial execution - Addition: 6 && Multiplication: 2
  int t1, t2, t3, t4, t5, t6, t7, t8;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
for (t1=-1;t1<=24;t1++) {
  lbp=ceild(t1,2);
  ubp=floord(t1+32,2);
#pragma omp parallel for private(lbv,ubv,t3,t4,t5,t6,t7,t8)
  for (t2=lbp;t2<=ubp;t2++) {
    for (t3=max(0,ceild(t1-1,2));t3<=floord(t1+33,2);t3++) {
      for (t4=max(max(0,ceild(t1-1,2)),t3-16);t4<=min(floord(t1+33,2),t3+16);t4++) {
        for (t5=max(max(max(max(0,8*t1),16*t3-256),16*t4-256),16*t1-16*t2+1);t5<=min(min(min(min(198,8*t1+15),16*t2+14),16*t3+14),16*t4+14);t5++) {
          for (t6=max(max(16*t2,t5+1),-16*t1+16*t2+2*t5-15);t6<=min(min(16*t2+15,t5+256),-16*t1+16*t2+2*t5);t6++) {
            for (t7=max(16*t3,t5+1);t7<=min(16*t3+15,t5+256);t7++) {
              lbv=max(16*t4,t5+1);
              ubv=min(16*t4+15,t5+256);
#pragma ivdep
#pragma vector always
              for (t8=lbv;t8<=ubv;t8++) {
                A[(t5 + 1) % 2][(-t5+t6)][(-t5+t7)][(-t5+t8)] = ((alpha * A[t5 % 2][(-t5+t6)][(-t5+t7)][(-t5+t8)]) + (beta * (((((A[t5 % 2][(-t5+t6) - 1][(-t5+t7)][(-t5+t8)] + A[t5 % 2][(-t5+t6)][(-t5+t7) - 1][(-t5+t8)]) + A[t5 % 2][(-t5+t6)][(-t5+t7)][(-t5+t8) - 1]) + A[t5 % 2][(-t5+t6) + 1][(-t5+t7)][(-t5+t8)]) + A[t5 % 2][(-t5+t6)][(-t5+t7) + 1][(-t5+t8)]) + A[t5 % 2][(-t5+t6)][(-t5+t7)][(-t5+t8) + 1])));;
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
  tdiff = (double) (result.tv_sec + result.tv_usec * 1.0e-6);

  printf("Time taken: %7.5lfms\n", tdiff * 1.0e3);
  printf("|MFLOPS: %f\t", ((((double)NUM_FP_OPS * N *N * N * (T-1)) / tdiff) / 1000000L));
#endif

#ifdef VERIFY
  // display the initial setting
  total = 0.0;
  for (i = 0; i < N; ++i) {
    for (j = 0; j < N; ++j) {
      for (k = 0; k < N; ++k) {
        total += A[T%2][i][j][k];
      }
    }
  }
  fprintf(stderr, "Sum(final): %e\n", total);
#endif

  return 0;
}

// icc -O3 -DTIME -DDEBUG 3d7pt.c -o op-3d7pt
