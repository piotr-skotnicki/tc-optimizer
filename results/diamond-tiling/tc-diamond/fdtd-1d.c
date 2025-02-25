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
/* ./tc ../examples/pluto-perfect/fdtd-1d.scop.c --diamond-tiling --omp-for-codegen --floyd-warshall-tc --inline --debug -b 32 --drop-bounds */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#pragma scop
{
  if (N >= 1) {
    if (N == 1) {
      for (int i0 = 1; i0 <= min(30, T); i0 += 1) {
        h[0] = (h[0] - (0.69999999999999996 * (e[1] - e[0])));
      }
    }
    for (int k = max(0, -N + 3); k <= 1; k += 1) {
      #pragma omp parallel for
      for (int ii1 = 0; ii1 <= min(floord(T + N, 32), k + floord(N - 14 * k - 3, 32)); ii1 += 1) {
        if (k == 1) {
          for (int i0 = max(1, -N + 32 * ii1); i0 <= min(min(min(31, T), 32 * ii1 + 30), N - 32 * ii1 + 30); i0 += 1) {
            for (int i1 = max(max(32 * ii1, 32 * ii1 + 2 * i0 - 31), i0 + 1); i1 <= min(min(32 * ii1 + 31, N + i0), 32 * ii1 + 2 * i0 + 1); i1 += 1) {
              if (N + i0 >= i1 + 1 && 32 * ii1 + 2 * i0 >= i1) {
                e[-i0 + i1] = (e[-i0 + i1] - (0.5 * (h[-i0 + i1] - h[-i0 + i1 - 1])));
              }
              if (i1 + 30 >= 32 * ii1 + 2 * i0) {
                h[-i0 + i1 - 1] = (h[-i0 + i1 - 1] - (0.69999999999999996 * (e[-i0 + i1] - e[-i0 + i1 - 1])));
              }
            }
          }
        } else {
          for (int i0 = 1; i0 <= min(min(15, T), N - 32 * ii1 - 2); i0 += 1) {
            for (int i1 = 32 * ii1 + 2 * i0 + 1; i1 <= min(32 * ii1 + 31, N + i0); i1 += 1) {
              if (N + i0 >= i1 + 1) {
                e[-i0 + i1] = (e[-i0 + i1] - (0.5 * (h[-i0 + i1] - h[-i0 + i1 - 1])));
              }
              if (i1 >= 32 * ii1 + 2 * i0 + 2) {
                h[-i0 + i1 - 1] = (h[-i0 + i1 - 1] - (0.69999999999999996 * (e[-i0 + i1] - e[-i0 + i1 - 1])));
              }
            }
          }
        }
      }
    }
    #pragma omp parallel for
    for (int ii1 = 1; ii1 <= min((N + 31) / 32, floord(T + N, 32)); ii1 += 1) {
      for (int i0 = max(16, -N + 32 * ii1); i0 <= min(31, T); i0 += 1) {
        for (int i1 = 32 * ii1; i1 <= min(32 * ii1 + 2 * i0 - 31, N + i0); i1 += 1) {
          if (32 * ii1 + 2 * i0 >= i1 + 32 && N + i0 >= i1 + 1) {
            e[-i0 + i1] = (e[-i0 + i1] - (0.5 * (h[-i0 + i1] - h[-i0 + i1 - 1])));
          }
          h[-i0 + i1 - 1] = (h[-i0 + i1 - 1] - (0.69999999999999996 * (e[-i0 + i1] - e[-i0 + i1 - 1])));
        }
      }
      if (ii1 == 1) {
        for (int i0 = 32; i0 <= min(62, T); i0 += 1) {
          for (int i1 = i0 + 1; i1 <= min(63, N + i0); i1 += 1) {
            if (N + i0 >= i1 + 1) {
              e[-i0 + i1] = (e[-i0 + i1] - (0.5 * (h[-i0 + i1] - h[-i0 + i1 - 1])));
            }
            h[-i0 + i1 - 1] = (h[-i0 + i1 - 1] - (0.69999999999999996 * (e[-i0 + i1] - e[-i0 + i1 - 1])));
          }
        }
      } else {
        for (int i0 = 32; i0 <= min(min(47, T), N - 32 * ii1 + 62); i0 += 1) {
          for (int i1 = 32 * ii1 + 2 * i0 - 63; i1 <= min(32 * ii1 + 31, N + i0); i1 += 1) {
            if (N + i0 >= i1 + 1) {
              e[-i0 + i1] = (e[-i0 + i1] - (0.5 * (h[-i0 + i1] - h[-i0 + i1 - 1])));
            }
            if (i1 + 62 >= 32 * ii1 + 2 * i0) {
              h[-i0 + i1 - 1] = (h[-i0 + i1 - 1] - (0.69999999999999996 * (e[-i0 + i1] - e[-i0 + i1 - 1])));
            }
          }
        }
      }
    }
  }
  for (int ii0 = 1; ii0 < min(floord(T + 16, 32), floord(T + N, 32)); ii0 += 1) {
    for (int k = max(1, -((N + 15) / 16) + 3); k <= 2; k += 1) {
      if (k == 1) {
        for (int i0 = max(32 * ii0, -N + 32 * ii0 + 32); i0 <= min(min(T, N + 32 * ii0 - 2), 32 * ii0 + 31); i0 += 1) {
          for (int i1 = max(32 * ii0 + 32, -32 * ii0 + 2 * i0 + 1); i1 <= min(min(32 * ii0 + 63, N + i0), -32 * ii0 + 2 * i0 + 33); i1 += 1) {
            if (N + i0 >= i1 + 1 && 2 * i0 + 32 >= 32 * ii0 + i1) {
              e[-i0 + i1] = (e[-i0 + i1] - (0.5 * (h[-i0 + i1] - h[-i0 + i1 - 1])));
            }
            if (32 * ii0 + i1 >= 2 * i0 + 2) {
              h[-i0 + i1 - 1] = (h[-i0 + i1 - 1] - (0.69999999999999996 * (e[-i0 + i1] - e[-i0 + i1 - 1])));
            }
          }
        }
      }
      #pragma omp parallel for
      for (int ii1 = ii0 - k + 3; ii1 <= min((T + N) / 32, ii0 + (N + 16 * k - 1) / 32); ii1 += 1) {
        if (k == 2) {
          for (int i0 = max(32 * ii0 + 16, -N + 32 * ii1); i0 <= min(min(min(T, 32 * ii0 + 47), 16 * ii0 + 16 * ii1 + 15), N + 64 * ii0 - 32 * ii1 + 62); i0 += 1) {
            for (int i1 = max(32 * ii1, -64 * ii0 + 32 * ii1 + 2 * i0 - 63); i1 <= min(min(32 * ii1 + 31, -64 * ii0 + 32 * ii1 + 2 * i0 - 31), N + i0); i1 += 1) {
              if (32 * ii1 + 2 * i0 >= 64 * ii0 + i1 + 32 && N + i0 >= i1 + 1) {
                e[-i0 + i1] = (e[-i0 + i1] - (0.5 * (h[-i0 + i1] - h[-i0 + i1 - 1])));
              }
              if (32 * ii0 + 31 >= i0) {
                h[-i0 + i1 - 1] = (h[-i0 + i1 - 1] - (0.69999999999999996 * (e[-i0 + i1] - e[-i0 + i1 - 1])));
              } else if (64 * ii0 + i1 + 62 >= 32 * ii1 + 2 * i0) {
                h[-i0 + i1 - 1] = (h[-i0 + i1 - 1] - (0.69999999999999996 * (e[-i0 + i1] - e[-i0 + i1 - 1])));
              }
            }
          }
          if (ii1 == ii0 + 1) {
            for (int i0 = 32 * ii0 + 32; i0 <= min(T, 32 * ii0 + 62); i0 += 1) {
              for (int i1 = i0 + 1; i1 <= min(32 * ii0 + 63, N + i0); i1 += 1) {
                if (N + i0 >= i1 + 1) {
                  e[-i0 + i1] = (e[-i0 + i1] - (0.5 * (h[-i0 + i1] - h[-i0 + i1 - 1])));
                }
                h[-i0 + i1 - 1] = (h[-i0 + i1 - 1] - (0.69999999999999996 * (e[-i0 + i1] - e[-i0 + i1 - 1])));
              }
            }
          }
        } else {
          for (int i0 = max(32 * ii0, -N + 32 * ii1); i0 <= min(min(T, 32 * ii0 + 31), N + 64 * ii0 - 32 * ii1 + 30); i0 += 1) {
            for (int i1 = max(32 * ii1, -64 * ii0 + 32 * ii1 + 2 * i0 - 31); i1 <= min(min(32 * ii1 + 31, N + i0), -64 * ii0 + 32 * ii1 + 2 * i0 + 1); i1 += 1) {
              if (N + i0 >= i1 + 1 && 32 * ii1 + 2 * i0 >= 64 * ii0 + i1) {
                e[-i0 + i1] = (e[-i0 + i1] - (0.5 * (h[-i0 + i1] - h[-i0 + i1 - 1])));
              }
              if (64 * ii0 + i1 + 30 >= 32 * ii1 + 2 * i0) {
                h[-i0 + i1 - 1] = (h[-i0 + i1 - 1] - (0.69999999999999996 * (e[-i0 + i1] - e[-i0 + i1 - 1])));
              }
            }
          }
        }
      }
    }
  }
  if (T >= 32 && (T - 16) % 32 >= 16) {
    for (int i0 = max(-(T % 32) + T, -(T % 32) + T - N + 32); i0 <= T; i0 += 1) {
      for (int i1 = -(T % 32) + T + 32; i1 <= min(N + i0, (T % 32) - T + 2 * i0 + 33); i1 += 1) {
        if (N + i0 >= i1 + 1 && ((T - 16) % 32) + 2 * i0 + 16 >= T + i1) {
          e[-i0 + i1] = (e[-i0 + i1] - (0.5 * (h[-i0 + i1] - h[-i0 + i1 - 1])));
        }
        h[-i0 + i1 - 1] = (h[-i0 + i1 - 1] - (0.69999999999999996 * (e[-i0 + i1] - e[-i0 + i1 - 1])));
      }
    }
    #pragma omp parallel for
    for (int ii1 = T / 32 + 2; ii1 <= (T + N) / 32; ii1 += 1) {
      for (int i0 = max(-N + 32 * ii1, -(T % 32) + T); i0 <= T; i0 += 1) {
        for (int i1 = 32 * ii1; i1 <= min(N + i0, 2 * (T % 32) - 2 * T + 32 * ii1 + 2 * i0 + 1); i1 += 1) {
          if (N + i0 >= i1 + 1 && 2 * ((T - 16) % 32) + 32 * ii1 + 2 * i0 >= 2 * T + i1 + 32) {
            e[-i0 + i1] = (e[-i0 + i1] - (0.5 * (h[-i0 + i1] - h[-i0 + i1 - 1])));
          }
          h[-i0 + i1 - 1] = (h[-i0 + i1 - 1] - (0.69999999999999996 * (e[-i0 + i1] - e[-i0 + i1 - 1])));
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
