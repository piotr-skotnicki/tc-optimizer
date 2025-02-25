#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <sys/time.h>
#include <assert.h>

#define N 1000000
#define T 10000
double h[N];
double e[N+1];
#define coeff1 0.5
#define coeff2 0.7

#ifdef TIME
#define IF_TIME(foo) foo;
#else
#define IF_TIME(foo)
#endif

void init_array()
{
    int i, j;

        for (j=0; j<N; j++) {
            h[j] = ((double)j)/N;
            e[j] = ((double)j)/N;
        }
}

void print_array()
{
    int i, j;

    for (j=0; j<N; j++) {
	    fprintf(stderr, "%lf ", h[j]);
	    if (j%80 == 79) fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
}

double rtclock()
{
    struct timezone Tzp;
    struct timeval Tp;
    int stat;
    stat = gettimeofday (&Tp, &Tzp);
    if (stat != 0) printf("Error return from gettimeofday: %d",stat);
    return(Tp.tv_sec + Tp.tv_usec*1.0e-6);
}

int main()
{
    int t, i, j, k, l;

    double t_start, t_end;

    init_array();

	IF_TIME(t_start = rtclock());

/* TC Optimizing Compiler 0.4.1 */
/* ./tc ../examples/pluto-perfect/fdtd-1d.scop.c --semi-diamond-tiling --omp-for-codegen --floyd-warshall-tc --inline --debug -b 16 --drop-bounds */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#pragma scop
if (T >= 1 && N >= 1) {
  for (int ii0 = 0; ii0 <= T / 16; ii0 += 1) {
    for (int k = max(max(0, -N + 2), -N - 8 * ii0 + N / 2 + 2); k <= min(min(2, -2 * ii0 + T / 8 + 1), -ii0 + (T + N - ii0 - 1) / 15 + 1); k += 1) {
      #pragma omp parallel for
      for (int ii1 = max(ii0, ii0 + k - 1); ii1 <= min(min(min(ii0 + (N + 15) / 16, (T + N) / 16), 2 * ii0 + k + floord(N - 6 * k - 3, 16)), ii0 + k + floord(N - 7 * k - 2, 16)); ii1 += 1) {
        if (k <= 1) {
          if (k == 1) {
            for (int i0 = max(max(1, 16 * ii0), -N + 16 * ii1); i0 <= min(min(min(T, 16 * ii0 + 15), 15 * ii0 + ii1 + 14), N + 32 * ii0 - 16 * ii1 + 14); i0 += 1) {
              for (int i1 = max(max(16 * ii0 + 16, 16 * ii1), -32 * ii0 + 16 * ii1 + 2 * i0 - 15); i1 <= min(min(16 * ii1 + 15, N + i0), -32 * ii0 + 16 * ii1 + 2 * i0 + 1); i1 += 1) {
                if (N + i0 >= i1 + 1 && 16 * ii1 + 2 * i0 >= 32 * ii0 + i1) {
                  e[-i0 + i1] = (e[-i0 + i1] - (0.5 * (h[-i0 + i1] - h[-i0 + i1 - 1])));
                }
                if (32 * ii0 + i1 + 14 >= 16 * ii1 + 2 * i0) {
                  h[-i0 + i1 - 1] = (h[-i0 + i1 - 1] - (0.69999999999999996 * (e[-i0 + i1] - e[-i0 + i1 - 1])));
                }
              }
              if (ii1 == ii0) {
                for (int i1 = i0 + 1; i1 <= min(min(16 * ii0 + 15, N + i0), -16 * ii0 + 2 * i0 + 1); i1 += 1) {
                  if (N + i0 >= i1 + 1 && 2 * i0 >= 16 * ii0 + i1) {
                    e[-i0 + i1] = (e[-i0 + i1] - (0.5 * (h[-i0 + i1] - h[-i0 + i1 - 1])));
                  }
                  h[-i0 + i1 - 1] = (h[-i0 + i1 - 1] - (0.69999999999999996 * (e[-i0 + i1] - e[-i0 + i1 - 1])));
                }
              }
            }
          } else {
            for (int i0 = max(1, 16 * ii0); i0 <= min(min(T, 16 * ii0 + 7), N + 32 * ii0 - 16 * ii1 - 2); i0 += 1) {
              for (int i1 = -32 * ii0 + 16 * ii1 + 2 * i0 + 1; i1 <= min(16 * ii1 + 15, N + i0); i1 += 1) {
                if (N + i0 >= i1 + 1) {
                  e[-i0 + i1] = (e[-i0 + i1] - (0.5 * (h[-i0 + i1] - h[-i0 + i1 - 1])));
                }
                if (32 * ii0 + i1 >= 16 * ii1 + 2 * i0 + 2) {
                  h[-i0 + i1 - 1] = (h[-i0 + i1 - 1] - (0.69999999999999996 * (e[-i0 + i1] - e[-i0 + i1 - 1])));
                }
              }
            }
          }
        } else {
          for (int i0 = max(16 * ii0 + 8, -N + 16 * ii1); i0 <= min(T, 16 * ii0 + 15); i0 += 1) {
            for (int i1 = 16 * ii1; i1 <= min(-32 * ii0 + 16 * ii1 + 2 * i0 - 15, N + i0); i1 += 1) {
              if (16 * ii1 + 2 * i0 >= 32 * ii0 + i1 + 16 && N + i0 >= i1 + 1) {
                e[-i0 + i1] = (e[-i0 + i1] - (0.5 * (h[-i0 + i1] - h[-i0 + i1 - 1])));
              }
              h[-i0 + i1 - 1] = (h[-i0 + i1 - 1] - (0.69999999999999996 * (e[-i0 + i1] - e[-i0 + i1 - 1])));
            }
          }
        }
      }
    }
  }
}
#pragma endscop

    IF_TIME(t_end = rtclock());
    IF_TIME(printf("%0.6lfs\n", t_end - t_start));

    if (fopen(".test", "r")) {
        print_array();
    }

    return 0;
}
