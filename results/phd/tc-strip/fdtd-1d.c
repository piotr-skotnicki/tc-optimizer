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

/* TC Optimizing Compiler 0.4.2 */
/* ./tc ../examples/pluto-perfect/fdtd-1d.scop.c --strip-tiling --omp-for-codegen --floyd-warshall-tc --inline --debug -b 32 --drop-bounds */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#pragma scop
if (N >= 1) {
  for (int ii0 = 0; ii0 <= floord(T - 1, 32); ii0 += 1) {
    if (N >= 2) {
      #pragma omp parallel for
      for (int ii1 = ii0; ii1 <= ii0 + floord(N - 2, 32); ii1 += 1) {
        if (ii1 >= ii0 + 1) {
          for (int i1 = 32 * ii1 + 2; i1 <= min(N + 32 * ii0, 32 * ii1 + 33); i1 += 1) {
            e[-32 * ii0 + i1 - 1] = (e[-32 * ii0 + i1 - 1] - (0.5 * (h[-32 * ii0 + i1 - 1] - h[-32 * ii0 + i1 - 2])));
          }
        } else {
          for (int i1 = 32 * ii0 + 2; i1 <= min(N + 32 * ii0, 32 * ii0 + 33); i1 += 1) {
            e[-32 * ii0 + i1 - 1] = (e[-32 * ii0 + i1 - 1] - (0.5 * (h[-32 * ii0 + i1 - 1] - h[-32 * ii0 + i1 - 2])));
          }
        }
      }
      #pragma omp parallel for
      for (int ii1 = ii0; ii1 <= ii0 + (N - 1) / 32; ii1 += 1) {
        for (int i1 = 32 * ii1 + 2; i1 <= min(N + 32 * ii0 + 1, 32 * ii1 + 33); i1 += 1) {
          h[-32 * ii0 + i1 - 2] = (h[-32 * ii0 + i1 - 2] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 1] - e[-32 * ii0 + i1 - 2])));
        }
      }
      if (T >= 32 * ii0 + 2) {
        #pragma omp parallel for
        for (int ii1 = ii0; ii1 <= ii0 + floord(N - 1, 32); ii1 += 1) {
          for (int i1 = max(32 * ii0 + 3, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 1, 32 * ii1 + 33); i1 += 1) {
            e[-32 * ii0 + i1 - 2] = (e[-32 * ii0 + i1 - 2] - (0.5 * (h[-32 * ii0 + i1 - 2] - h[-32 * ii0 + i1 - 3])));
          }
        }
        #pragma omp parallel for
        for (int ii1 = ii0; ii1 <= ii0 + N / 32; ii1 += 1) {
          for (int i1 = max(32 * ii0 + 3, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 2, 32 * ii1 + 33); i1 += 1) {
            h[-32 * ii0 + i1 - 3] = (h[-32 * ii0 + i1 - 3] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 2] - e[-32 * ii0 + i1 - 3])));
          }
        }
        if (T >= 32 * ii0 + 3) {
          #pragma omp parallel for
          for (int ii1 = ii0; ii1 <= ii0 + N / 32; ii1 += 1) {
            for (int i1 = max(32 * ii0 + 4, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 2, 32 * ii1 + 33); i1 += 1) {
              e[-32 * ii0 + i1 - 3] = (e[-32 * ii0 + i1 - 3] - (0.5 * (h[-32 * ii0 + i1 - 3] - h[-32 * ii0 + i1 - 4])));
            }
          }
          if (T >= 32 * ii0 + 32) {
            #pragma omp parallel for
            for (int ii1 = ii0; ii1 <= min((T + N - 1) / 32 - 1, ii0 + (N + 1) / 32); ii1 += 1) {
              for (int i1 = max(32 * ii0 + 4, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 3, 32 * ii1 + 33); i1 += 1) {
                h[-32 * ii0 + i1 - 4] = (h[-32 * ii0 + i1 - 4] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 3] - e[-32 * ii0 + i1 - 4])));
              }
            }
            for (int i1 = -((T + N + 31) % 32) + T + N + 1; i1 <= N + 32 * ii0 + 3; i1 += 1) {
              h[-32 * ii0 + i1 - 4] = (h[-32 * ii0 + i1 - 4] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 3] - e[-32 * ii0 + i1 - 4])));
            }
          } else {
            #pragma omp parallel for
            for (int ii1 = ii0; ii1 <= ii0 + (N + 1) / 32; ii1 += 1) {
              for (int i1 = max(32 * ii0 + 4, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 3, 32 * ii1 + 33); i1 += 1) {
                h[-32 * ii0 + i1 - 4] = (h[-32 * ii0 + i1 - 4] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 3] - e[-32 * ii0 + i1 - 4])));
              }
            }
          }
          if (T >= 32 * ii0 + 4) {
            #pragma omp parallel for
            for (int ii1 = ii0; ii1 <= ii0 + (N + 1) / 32; ii1 += 1) {
              for (int i1 = max(32 * ii0 + 5, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 3, 32 * ii1 + 33); i1 += 1) {
                e[-32 * ii0 + i1 - 4] = (e[-32 * ii0 + i1 - 4] - (0.5 * (h[-32 * ii0 + i1 - 4] - h[-32 * ii0 + i1 - 5])));
              }
            }
            #pragma omp parallel for
            for (int ii1 = ii0; ii1 <= ii0 + (N + 2) / 32; ii1 += 1) {
              if (T >= 32 * ii0 + 34 && T + N >= 32 * ii1 + 33) {
                for (int i1 = max(32 * ii0 + 5, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 4, 32 * ii1 + 33); i1 += 1) {
                  h[-32 * ii0 + i1 - 5] = (h[-32 * ii0 + i1 - 5] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 4] - e[-32 * ii0 + i1 - 5])));
                }
              } else if (T >= 32 * ii0 + 32 && 32 * ii1 + 32 >= T + N) {
                for (int i1 = 32 * ii1 + 2; i1 <= N + 32 * ii0 + 4; i1 += 1) {
                  h[-32 * ii0 + i1 - 5] = (h[-32 * ii0 + i1 - 5] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 4] - e[-32 * ii0 + i1 - 5])));
                }
              } else if (32 * ii0 + 31 >= T) {
                for (int i1 = max(32 * ii0 + 5, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 4, 32 * ii1 + 33); i1 += 1) {
                  h[-32 * ii0 + i1 - 5] = (h[-32 * ii0 + i1 - 5] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 4] - e[-32 * ii0 + i1 - 5])));
                }
              } else {
                for (int i1 = max(32 * ii0 + 5, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 4, 32 * ii1 + 33); i1 += 1) {
                  h[-32 * ii0 + i1 - 5] = (h[-32 * ii0 + i1 - 5] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 4] - e[-32 * ii0 + i1 - 5])));
                }
              }
            }
            if (T >= 32 * ii0 + 5) {
              #pragma omp parallel for
              for (int ii1 = ii0; ii1 <= ii0 + (N + 2) / 32; ii1 += 1) {
                for (int i1 = max(32 * ii0 + 6, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 4, 32 * ii1 + 33); i1 += 1) {
                  e[-32 * ii0 + i1 - 5] = (e[-32 * ii0 + i1 - 5] - (0.5 * (h[-32 * ii0 + i1 - 5] - h[-32 * ii0 + i1 - 6])));
                }
              }
              #pragma omp parallel for
              for (int ii1 = ii0; ii1 <= ii0 + (N + 3) / 32; ii1 += 1) {
                if (T >= 32 * ii0 + 35 && T + N >= 32 * ii1 + 33) {
                  for (int i1 = max(32 * ii0 + 6, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 5, 32 * ii1 + 33); i1 += 1) {
                    h[-32 * ii0 + i1 - 6] = (h[-32 * ii0 + i1 - 6] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 5] - e[-32 * ii0 + i1 - 6])));
                  }
                } else if (T >= 32 * ii0 + 32 && 32 * ii1 + 32 >= T + N) {
                  for (int i1 = 32 * ii1 + 2; i1 <= N + 32 * ii0 + 5; i1 += 1) {
                    h[-32 * ii0 + i1 - 6] = (h[-32 * ii0 + i1 - 6] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 5] - e[-32 * ii0 + i1 - 6])));
                  }
                } else if (32 * ii0 + 31 >= T) {
                  for (int i1 = max(32 * ii0 + 6, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 5, 32 * ii1 + 33); i1 += 1) {
                    h[-32 * ii0 + i1 - 6] = (h[-32 * ii0 + i1 - 6] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 5] - e[-32 * ii0 + i1 - 6])));
                  }
                } else {
                  for (int i1 = max(32 * ii0 + 6, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 5, 32 * ii1 + 33); i1 += 1) {
                    h[-32 * ii0 + i1 - 6] = (h[-32 * ii0 + i1 - 6] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 5] - e[-32 * ii0 + i1 - 6])));
                  }
                }
              }
              if (T >= 32 * ii0 + 6) {
                #pragma omp parallel for
                for (int ii1 = ii0; ii1 <= ii0 + (N + 3) / 32; ii1 += 1) {
                  for (int i1 = max(32 * ii0 + 7, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 5, 32 * ii1 + 33); i1 += 1) {
                    e[-32 * ii0 + i1 - 6] = (e[-32 * ii0 + i1 - 6] - (0.5 * (h[-32 * ii0 + i1 - 6] - h[-32 * ii0 + i1 - 7])));
                  }
                }
                #pragma omp parallel for
                for (int ii1 = ii0; ii1 <= ii0 + (N + 4) / 32; ii1 += 1) {
                  if (T >= 32 * ii0 + 36 && T + N >= 32 * ii1 + 33) {
                    for (int i1 = max(32 * ii0 + 7, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 6, 32 * ii1 + 33); i1 += 1) {
                      h[-32 * ii0 + i1 - 7] = (h[-32 * ii0 + i1 - 7] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 6] - e[-32 * ii0 + i1 - 7])));
                    }
                  } else if (T >= 32 * ii0 + 32 && 32 * ii1 + 32 >= T + N) {
                    for (int i1 = 32 * ii1 + 2; i1 <= N + 32 * ii0 + 6; i1 += 1) {
                      h[-32 * ii0 + i1 - 7] = (h[-32 * ii0 + i1 - 7] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 6] - e[-32 * ii0 + i1 - 7])));
                    }
                  } else if (32 * ii0 + 31 >= T) {
                    for (int i1 = max(32 * ii0 + 7, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 6, 32 * ii1 + 33); i1 += 1) {
                      h[-32 * ii0 + i1 - 7] = (h[-32 * ii0 + i1 - 7] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 6] - e[-32 * ii0 + i1 - 7])));
                    }
                  } else {
                    for (int i1 = max(32 * ii0 + 7, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 6, 32 * ii1 + 33); i1 += 1) {
                      h[-32 * ii0 + i1 - 7] = (h[-32 * ii0 + i1 - 7] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 6] - e[-32 * ii0 + i1 - 7])));
                    }
                  }
                }
                if (T >= 32 * ii0 + 7) {
                  #pragma omp parallel for
                  for (int ii1 = ii0; ii1 <= ii0 + (N + 4) / 32; ii1 += 1) {
                    for (int i1 = max(32 * ii0 + 8, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 6, 32 * ii1 + 33); i1 += 1) {
                      e[-32 * ii0 + i1 - 7] = (e[-32 * ii0 + i1 - 7] - (0.5 * (h[-32 * ii0 + i1 - 7] - h[-32 * ii0 + i1 - 8])));
                    }
                  }
                  #pragma omp parallel for
                  for (int ii1 = ii0; ii1 <= ii0 + (N + 5) / 32; ii1 += 1) {
                    if (T >= 32 * ii0 + 37 && T + N >= 32 * ii1 + 33) {
                      for (int i1 = max(32 * ii0 + 8, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 7, 32 * ii1 + 33); i1 += 1) {
                        h[-32 * ii0 + i1 - 8] = (h[-32 * ii0 + i1 - 8] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 7] - e[-32 * ii0 + i1 - 8])));
                      }
                    } else if (T >= 32 * ii0 + 32 && 32 * ii1 + 32 >= T + N) {
                      for (int i1 = 32 * ii1 + 2; i1 <= N + 32 * ii0 + 7; i1 += 1) {
                        h[-32 * ii0 + i1 - 8] = (h[-32 * ii0 + i1 - 8] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 7] - e[-32 * ii0 + i1 - 8])));
                      }
                    } else if (32 * ii0 + 31 >= T) {
                      for (int i1 = max(32 * ii0 + 8, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 7, 32 * ii1 + 33); i1 += 1) {
                        h[-32 * ii0 + i1 - 8] = (h[-32 * ii0 + i1 - 8] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 7] - e[-32 * ii0 + i1 - 8])));
                      }
                    } else {
                      for (int i1 = max(32 * ii0 + 8, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 7, 32 * ii1 + 33); i1 += 1) {
                        h[-32 * ii0 + i1 - 8] = (h[-32 * ii0 + i1 - 8] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 7] - e[-32 * ii0 + i1 - 8])));
                      }
                    }
                  }
                  if (T >= 32 * ii0 + 8) {
                    #pragma omp parallel for
                    for (int ii1 = ii0; ii1 <= ii0 + (N + 5) / 32; ii1 += 1) {
                      for (int i1 = max(32 * ii0 + 9, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 7, 32 * ii1 + 33); i1 += 1) {
                        e[-32 * ii0 + i1 - 8] = (e[-32 * ii0 + i1 - 8] - (0.5 * (h[-32 * ii0 + i1 - 8] - h[-32 * ii0 + i1 - 9])));
                      }
                    }
                    #pragma omp parallel for
                    for (int ii1 = ii0; ii1 <= ii0 + (N + 6) / 32; ii1 += 1) {
                      if (T >= 32 * ii0 + 38 && T + N >= 32 * ii1 + 33) {
                        for (int i1 = max(32 * ii0 + 9, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 8, 32 * ii1 + 33); i1 += 1) {
                          h[-32 * ii0 + i1 - 9] = (h[-32 * ii0 + i1 - 9] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 8] - e[-32 * ii0 + i1 - 9])));
                        }
                      } else if (T >= 32 * ii0 + 32 && 32 * ii1 + 32 >= T + N) {
                        for (int i1 = 32 * ii1 + 2; i1 <= N + 32 * ii0 + 8; i1 += 1) {
                          h[-32 * ii0 + i1 - 9] = (h[-32 * ii0 + i1 - 9] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 8] - e[-32 * ii0 + i1 - 9])));
                        }
                      } else if (32 * ii0 + 31 >= T) {
                        for (int i1 = max(32 * ii0 + 9, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 8, 32 * ii1 + 33); i1 += 1) {
                          h[-32 * ii0 + i1 - 9] = (h[-32 * ii0 + i1 - 9] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 8] - e[-32 * ii0 + i1 - 9])));
                        }
                      } else {
                        for (int i1 = max(32 * ii0 + 9, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 8, 32 * ii1 + 33); i1 += 1) {
                          h[-32 * ii0 + i1 - 9] = (h[-32 * ii0 + i1 - 9] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 8] - e[-32 * ii0 + i1 - 9])));
                        }
                      }
                    }
                    if (T >= 32 * ii0 + 9) {
                      #pragma omp parallel for
                      for (int ii1 = ii0; ii1 <= ii0 + (N + 6) / 32; ii1 += 1) {
                        for (int i1 = max(32 * ii0 + 10, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 8, 32 * ii1 + 33); i1 += 1) {
                          e[-32 * ii0 + i1 - 9] = (e[-32 * ii0 + i1 - 9] - (0.5 * (h[-32 * ii0 + i1 - 9] - h[-32 * ii0 + i1 - 10])));
                        }
                      }
                      #pragma omp parallel for
                      for (int ii1 = ii0; ii1 <= ii0 + (N + 7) / 32; ii1 += 1) {
                        if (T >= 32 * ii0 + 39 && T + N >= 32 * ii1 + 33) {
                          for (int i1 = max(32 * ii0 + 10, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 9, 32 * ii1 + 33); i1 += 1) {
                            h[-32 * ii0 + i1 - 10] = (h[-32 * ii0 + i1 - 10] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 9] - e[-32 * ii0 + i1 - 10])));
                          }
                        } else if (T >= 32 * ii0 + 32 && 32 * ii1 + 32 >= T + N) {
                          for (int i1 = 32 * ii1 + 2; i1 <= N + 32 * ii0 + 9; i1 += 1) {
                            h[-32 * ii0 + i1 - 10] = (h[-32 * ii0 + i1 - 10] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 9] - e[-32 * ii0 + i1 - 10])));
                          }
                        } else if (32 * ii0 + 31 >= T) {
                          for (int i1 = max(32 * ii0 + 10, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 9, 32 * ii1 + 33); i1 += 1) {
                            h[-32 * ii0 + i1 - 10] = (h[-32 * ii0 + i1 - 10] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 9] - e[-32 * ii0 + i1 - 10])));
                          }
                        } else {
                          for (int i1 = max(32 * ii0 + 10, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 9, 32 * ii1 + 33); i1 += 1) {
                            h[-32 * ii0 + i1 - 10] = (h[-32 * ii0 + i1 - 10] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 9] - e[-32 * ii0 + i1 - 10])));
                          }
                        }
                      }
                      if (T >= 32 * ii0 + 10) {
                        #pragma omp parallel for
                        for (int ii1 = ii0; ii1 <= ii0 + (N + 7) / 32; ii1 += 1) {
                          for (int i1 = max(32 * ii0 + 11, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 9, 32 * ii1 + 33); i1 += 1) {
                            e[-32 * ii0 + i1 - 10] = (e[-32 * ii0 + i1 - 10] - (0.5 * (h[-32 * ii0 + i1 - 10] - h[-32 * ii0 + i1 - 11])));
                          }
                        }
                        #pragma omp parallel for
                        for (int ii1 = ii0; ii1 <= ii0 + (N + 8) / 32; ii1 += 1) {
                          if (T >= 32 * ii0 + 40 && T + N >= 32 * ii1 + 33) {
                            for (int i1 = max(32 * ii0 + 11, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 10, 32 * ii1 + 33); i1 += 1) {
                              h[-32 * ii0 + i1 - 11] = (h[-32 * ii0 + i1 - 11] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 10] - e[-32 * ii0 + i1 - 11])));
                            }
                          } else if (T >= 32 * ii0 + 32 && 32 * ii1 + 32 >= T + N) {
                            for (int i1 = 32 * ii1 + 2; i1 <= N + 32 * ii0 + 10; i1 += 1) {
                              h[-32 * ii0 + i1 - 11] = (h[-32 * ii0 + i1 - 11] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 10] - e[-32 * ii0 + i1 - 11])));
                            }
                          } else if (32 * ii0 + 31 >= T) {
                            for (int i1 = max(32 * ii0 + 11, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 10, 32 * ii1 + 33); i1 += 1) {
                              h[-32 * ii0 + i1 - 11] = (h[-32 * ii0 + i1 - 11] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 10] - e[-32 * ii0 + i1 - 11])));
                            }
                          } else {
                            for (int i1 = max(32 * ii0 + 11, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 10, 32 * ii1 + 33); i1 += 1) {
                              h[-32 * ii0 + i1 - 11] = (h[-32 * ii0 + i1 - 11] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 10] - e[-32 * ii0 + i1 - 11])));
                            }
                          }
                        }
                        if (T >= 32 * ii0 + 11) {
                          #pragma omp parallel for
                          for (int ii1 = ii0; ii1 <= ii0 + (N + 8) / 32; ii1 += 1) {
                            for (int i1 = max(32 * ii0 + 12, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 10, 32 * ii1 + 33); i1 += 1) {
                              e[-32 * ii0 + i1 - 11] = (e[-32 * ii0 + i1 - 11] - (0.5 * (h[-32 * ii0 + i1 - 11] - h[-32 * ii0 + i1 - 12])));
                            }
                          }
                          #pragma omp parallel for
                          for (int ii1 = ii0; ii1 <= ii0 + (N + 9) / 32; ii1 += 1) {
                            if (T >= 32 * ii0 + 41 && T + N >= 32 * ii1 + 33) {
                              for (int i1 = max(32 * ii0 + 12, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 11, 32 * ii1 + 33); i1 += 1) {
                                h[-32 * ii0 + i1 - 12] = (h[-32 * ii0 + i1 - 12] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 11] - e[-32 * ii0 + i1 - 12])));
                              }
                            } else if (T >= 32 * ii0 + 32 && 32 * ii1 + 32 >= T + N) {
                              for (int i1 = 32 * ii1 + 2; i1 <= N + 32 * ii0 + 11; i1 += 1) {
                                h[-32 * ii0 + i1 - 12] = (h[-32 * ii0 + i1 - 12] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 11] - e[-32 * ii0 + i1 - 12])));
                              }
                            } else if (32 * ii0 + 31 >= T) {
                              for (int i1 = max(32 * ii0 + 12, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 11, 32 * ii1 + 33); i1 += 1) {
                                h[-32 * ii0 + i1 - 12] = (h[-32 * ii0 + i1 - 12] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 11] - e[-32 * ii0 + i1 - 12])));
                              }
                            } else {
                              for (int i1 = max(32 * ii0 + 12, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 11, 32 * ii1 + 33); i1 += 1) {
                                h[-32 * ii0 + i1 - 12] = (h[-32 * ii0 + i1 - 12] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 11] - e[-32 * ii0 + i1 - 12])));
                              }
                            }
                          }
                          if (T >= 32 * ii0 + 12) {
                            #pragma omp parallel for
                            for (int ii1 = ii0; ii1 <= ii0 + (N + 9) / 32; ii1 += 1) {
                              for (int i1 = max(32 * ii0 + 13, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 11, 32 * ii1 + 33); i1 += 1) {
                                e[-32 * ii0 + i1 - 12] = (e[-32 * ii0 + i1 - 12] - (0.5 * (h[-32 * ii0 + i1 - 12] - h[-32 * ii0 + i1 - 13])));
                              }
                            }
                            #pragma omp parallel for
                            for (int ii1 = ii0; ii1 <= ii0 + (N + 10) / 32; ii1 += 1) {
                              if (T >= 32 * ii0 + 42 && T + N >= 32 * ii1 + 33) {
                                for (int i1 = max(32 * ii0 + 13, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 12, 32 * ii1 + 33); i1 += 1) {
                                  h[-32 * ii0 + i1 - 13] = (h[-32 * ii0 + i1 - 13] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 12] - e[-32 * ii0 + i1 - 13])));
                                }
                              } else if (T >= 32 * ii0 + 32 && 32 * ii1 + 32 >= T + N) {
                                for (int i1 = 32 * ii1 + 2; i1 <= N + 32 * ii0 + 12; i1 += 1) {
                                  h[-32 * ii0 + i1 - 13] = (h[-32 * ii0 + i1 - 13] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 12] - e[-32 * ii0 + i1 - 13])));
                                }
                              } else if (32 * ii0 + 31 >= T) {
                                for (int i1 = max(32 * ii0 + 13, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 12, 32 * ii1 + 33); i1 += 1) {
                                  h[-32 * ii0 + i1 - 13] = (h[-32 * ii0 + i1 - 13] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 12] - e[-32 * ii0 + i1 - 13])));
                                }
                              } else {
                                for (int i1 = max(32 * ii0 + 13, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 12, 32 * ii1 + 33); i1 += 1) {
                                  h[-32 * ii0 + i1 - 13] = (h[-32 * ii0 + i1 - 13] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 12] - e[-32 * ii0 + i1 - 13])));
                                }
                              }
                            }
                            if (T >= 32 * ii0 + 13) {
                              #pragma omp parallel for
                              for (int ii1 = ii0; ii1 <= ii0 + (N + 10) / 32; ii1 += 1) {
                                for (int i1 = max(32 * ii0 + 14, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 12, 32 * ii1 + 33); i1 += 1) {
                                  e[-32 * ii0 + i1 - 13] = (e[-32 * ii0 + i1 - 13] - (0.5 * (h[-32 * ii0 + i1 - 13] - h[-32 * ii0 + i1 - 14])));
                                }
                              }
                              #pragma omp parallel for
                              for (int ii1 = ii0; ii1 <= ii0 + (N + 11) / 32; ii1 += 1) {
                                if (T >= 32 * ii0 + 43 && T + N >= 32 * ii1 + 33) {
                                  for (int i1 = max(32 * ii0 + 14, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 13, 32 * ii1 + 33); i1 += 1) {
                                    h[-32 * ii0 + i1 - 14] = (h[-32 * ii0 + i1 - 14] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 13] - e[-32 * ii0 + i1 - 14])));
                                  }
                                } else if (T >= 32 * ii0 + 32 && 32 * ii1 + 32 >= T + N) {
                                  for (int i1 = 32 * ii1 + 2; i1 <= N + 32 * ii0 + 13; i1 += 1) {
                                    h[-32 * ii0 + i1 - 14] = (h[-32 * ii0 + i1 - 14] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 13] - e[-32 * ii0 + i1 - 14])));
                                  }
                                } else if (32 * ii0 + 31 >= T) {
                                  for (int i1 = max(32 * ii0 + 14, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 13, 32 * ii1 + 33); i1 += 1) {
                                    h[-32 * ii0 + i1 - 14] = (h[-32 * ii0 + i1 - 14] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 13] - e[-32 * ii0 + i1 - 14])));
                                  }
                                } else {
                                  for (int i1 = max(32 * ii0 + 14, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 13, 32 * ii1 + 33); i1 += 1) {
                                    h[-32 * ii0 + i1 - 14] = (h[-32 * ii0 + i1 - 14] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 13] - e[-32 * ii0 + i1 - 14])));
                                  }
                                }
                              }
                              if (T >= 32 * ii0 + 14) {
                                #pragma omp parallel for
                                for (int ii1 = ii0; ii1 <= ii0 + (N + 11) / 32; ii1 += 1) {
                                  for (int i1 = max(32 * ii0 + 15, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 13, 32 * ii1 + 33); i1 += 1) {
                                    e[-32 * ii0 + i1 - 14] = (e[-32 * ii0 + i1 - 14] - (0.5 * (h[-32 * ii0 + i1 - 14] - h[-32 * ii0 + i1 - 15])));
                                  }
                                }
                                #pragma omp parallel for
                                for (int ii1 = ii0; ii1 <= ii0 + (N + 12) / 32; ii1 += 1) {
                                  if (T >= 32 * ii0 + 44 && T + N >= 32 * ii1 + 33) {
                                    for (int i1 = max(32 * ii0 + 15, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 14, 32 * ii1 + 33); i1 += 1) {
                                      h[-32 * ii0 + i1 - 15] = (h[-32 * ii0 + i1 - 15] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 14] - e[-32 * ii0 + i1 - 15])));
                                    }
                                  } else if (T >= 32 * ii0 + 32 && 32 * ii1 + 32 >= T + N) {
                                    for (int i1 = 32 * ii1 + 2; i1 <= N + 32 * ii0 + 14; i1 += 1) {
                                      h[-32 * ii0 + i1 - 15] = (h[-32 * ii0 + i1 - 15] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 14] - e[-32 * ii0 + i1 - 15])));
                                    }
                                  } else if (32 * ii0 + 31 >= T) {
                                    for (int i1 = max(32 * ii0 + 15, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 14, 32 * ii1 + 33); i1 += 1) {
                                      h[-32 * ii0 + i1 - 15] = (h[-32 * ii0 + i1 - 15] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 14] - e[-32 * ii0 + i1 - 15])));
                                    }
                                  } else {
                                    for (int i1 = max(32 * ii0 + 15, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 14, 32 * ii1 + 33); i1 += 1) {
                                      h[-32 * ii0 + i1 - 15] = (h[-32 * ii0 + i1 - 15] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 14] - e[-32 * ii0 + i1 - 15])));
                                    }
                                  }
                                }
                                if (T >= 32 * ii0 + 15) {
                                  #pragma omp parallel for
                                  for (int ii1 = ii0; ii1 <= ii0 + (N + 12) / 32; ii1 += 1) {
                                    for (int i1 = max(32 * ii0 + 16, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 14, 32 * ii1 + 33); i1 += 1) {
                                      e[-32 * ii0 + i1 - 15] = (e[-32 * ii0 + i1 - 15] - (0.5 * (h[-32 * ii0 + i1 - 15] - h[-32 * ii0 + i1 - 16])));
                                    }
                                  }
                                  #pragma omp parallel for
                                  for (int ii1 = ii0; ii1 <= ii0 + (N + 13) / 32; ii1 += 1) {
                                    if (T >= 32 * ii0 + 45 && T + N >= 32 * ii1 + 33) {
                                      for (int i1 = max(32 * ii0 + 16, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 15, 32 * ii1 + 33); i1 += 1) {
                                        h[-32 * ii0 + i1 - 16] = (h[-32 * ii0 + i1 - 16] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 15] - e[-32 * ii0 + i1 - 16])));
                                      }
                                    } else if (T >= 32 * ii0 + 32 && 32 * ii1 + 32 >= T + N) {
                                      for (int i1 = 32 * ii1 + 2; i1 <= N + 32 * ii0 + 15; i1 += 1) {
                                        h[-32 * ii0 + i1 - 16] = (h[-32 * ii0 + i1 - 16] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 15] - e[-32 * ii0 + i1 - 16])));
                                      }
                                    } else if (32 * ii0 + 31 >= T) {
                                      for (int i1 = max(32 * ii0 + 16, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 15, 32 * ii1 + 33); i1 += 1) {
                                        h[-32 * ii0 + i1 - 16] = (h[-32 * ii0 + i1 - 16] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 15] - e[-32 * ii0 + i1 - 16])));
                                      }
                                    } else {
                                      for (int i1 = max(32 * ii0 + 16, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 15, 32 * ii1 + 33); i1 += 1) {
                                        h[-32 * ii0 + i1 - 16] = (h[-32 * ii0 + i1 - 16] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 15] - e[-32 * ii0 + i1 - 16])));
                                      }
                                    }
                                  }
                                  if (T >= 32 * ii0 + 16) {
                                    #pragma omp parallel for
                                    for (int ii1 = ii0; ii1 <= ii0 + (N + 13) / 32; ii1 += 1) {
                                      for (int i1 = max(32 * ii0 + 17, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 15, 32 * ii1 + 33); i1 += 1) {
                                        e[-32 * ii0 + i1 - 16] = (e[-32 * ii0 + i1 - 16] - (0.5 * (h[-32 * ii0 + i1 - 16] - h[-32 * ii0 + i1 - 17])));
                                      }
                                    }
                                    #pragma omp parallel for
                                    for (int ii1 = ii0; ii1 <= ii0 + (N + 14) / 32; ii1 += 1) {
                                      if (T >= 32 * ii0 + 46 && T + N >= 32 * ii1 + 33) {
                                        for (int i1 = max(32 * ii0 + 17, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 16, 32 * ii1 + 33); i1 += 1) {
                                          h[-32 * ii0 + i1 - 17] = (h[-32 * ii0 + i1 - 17] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 16] - e[-32 * ii0 + i1 - 17])));
                                        }
                                      } else if (T >= 32 * ii0 + 32 && 32 * ii1 + 32 >= T + N) {
                                        for (int i1 = 32 * ii1 + 2; i1 <= N + 32 * ii0 + 16; i1 += 1) {
                                          h[-32 * ii0 + i1 - 17] = (h[-32 * ii0 + i1 - 17] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 16] - e[-32 * ii0 + i1 - 17])));
                                        }
                                      } else if (32 * ii0 + 31 >= T) {
                                        for (int i1 = max(32 * ii0 + 17, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 16, 32 * ii1 + 33); i1 += 1) {
                                          h[-32 * ii0 + i1 - 17] = (h[-32 * ii0 + i1 - 17] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 16] - e[-32 * ii0 + i1 - 17])));
                                        }
                                      } else {
                                        for (int i1 = max(32 * ii0 + 17, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 16, 32 * ii1 + 33); i1 += 1) {
                                          h[-32 * ii0 + i1 - 17] = (h[-32 * ii0 + i1 - 17] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 16] - e[-32 * ii0 + i1 - 17])));
                                        }
                                      }
                                    }
                                    if (T >= 32 * ii0 + 17) {
                                      #pragma omp parallel for
                                      for (int ii1 = ii0; ii1 <= ii0 + (N + 14) / 32; ii1 += 1) {
                                        for (int i1 = max(32 * ii0 + 18, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 16, 32 * ii1 + 33); i1 += 1) {
                                          e[-32 * ii0 + i1 - 17] = (e[-32 * ii0 + i1 - 17] - (0.5 * (h[-32 * ii0 + i1 - 17] - h[-32 * ii0 + i1 - 18])));
                                        }
                                      }
                                      #pragma omp parallel for
                                      for (int ii1 = ii0; ii1 <= ii0 + (N + 15) / 32; ii1 += 1) {
                                        if (T >= 32 * ii0 + 47 && T + N >= 32 * ii1 + 33) {
                                          for (int i1 = max(32 * ii0 + 18, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 17, 32 * ii1 + 33); i1 += 1) {
                                            h[-32 * ii0 + i1 - 18] = (h[-32 * ii0 + i1 - 18] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 17] - e[-32 * ii0 + i1 - 18])));
                                          }
                                        } else if (T >= 32 * ii0 + 32 && 32 * ii1 + 32 >= T + N) {
                                          for (int i1 = 32 * ii1 + 2; i1 <= N + 32 * ii0 + 17; i1 += 1) {
                                            h[-32 * ii0 + i1 - 18] = (h[-32 * ii0 + i1 - 18] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 17] - e[-32 * ii0 + i1 - 18])));
                                          }
                                        } else if (32 * ii0 + 31 >= T) {
                                          for (int i1 = max(32 * ii0 + 18, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 17, 32 * ii1 + 33); i1 += 1) {
                                            h[-32 * ii0 + i1 - 18] = (h[-32 * ii0 + i1 - 18] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 17] - e[-32 * ii0 + i1 - 18])));
                                          }
                                        } else {
                                          for (int i1 = max(32 * ii0 + 18, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 17, 32 * ii1 + 33); i1 += 1) {
                                            h[-32 * ii0 + i1 - 18] = (h[-32 * ii0 + i1 - 18] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 17] - e[-32 * ii0 + i1 - 18])));
                                          }
                                        }
                                      }
                                      if (T >= 32 * ii0 + 18) {
                                        #pragma omp parallel for
                                        for (int ii1 = ii0; ii1 <= ii0 + (N + 15) / 32; ii1 += 1) {
                                          for (int i1 = max(32 * ii0 + 19, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 17, 32 * ii1 + 33); i1 += 1) {
                                            e[-32 * ii0 + i1 - 18] = (e[-32 * ii0 + i1 - 18] - (0.5 * (h[-32 * ii0 + i1 - 18] - h[-32 * ii0 + i1 - 19])));
                                          }
                                        }
                                        #pragma omp parallel for
                                        for (int ii1 = ii0; ii1 <= ii0 + (N + 16) / 32; ii1 += 1) {
                                          if (T >= 32 * ii0 + 48 && T + N >= 32 * ii1 + 33) {
                                            for (int i1 = max(32 * ii0 + 19, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 18, 32 * ii1 + 33); i1 += 1) {
                                              h[-32 * ii0 + i1 - 19] = (h[-32 * ii0 + i1 - 19] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 18] - e[-32 * ii0 + i1 - 19])));
                                            }
                                          } else if (T >= 32 * ii0 + 32 && 32 * ii1 + 32 >= T + N) {
                                            for (int i1 = 32 * ii1 + 2; i1 <= N + 32 * ii0 + 18; i1 += 1) {
                                              h[-32 * ii0 + i1 - 19] = (h[-32 * ii0 + i1 - 19] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 18] - e[-32 * ii0 + i1 - 19])));
                                            }
                                          } else if (32 * ii0 + 31 >= T) {
                                            for (int i1 = max(32 * ii0 + 19, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 18, 32 * ii1 + 33); i1 += 1) {
                                              h[-32 * ii0 + i1 - 19] = (h[-32 * ii0 + i1 - 19] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 18] - e[-32 * ii0 + i1 - 19])));
                                            }
                                          } else {
                                            for (int i1 = max(32 * ii0 + 19, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 18, 32 * ii1 + 33); i1 += 1) {
                                              h[-32 * ii0 + i1 - 19] = (h[-32 * ii0 + i1 - 19] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 18] - e[-32 * ii0 + i1 - 19])));
                                            }
                                          }
                                        }
                                        if (T >= 32 * ii0 + 19) {
                                          #pragma omp parallel for
                                          for (int ii1 = ii0; ii1 <= ii0 + (N + 16) / 32; ii1 += 1) {
                                            for (int i1 = max(32 * ii0 + 20, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 18, 32 * ii1 + 33); i1 += 1) {
                                              e[-32 * ii0 + i1 - 19] = (e[-32 * ii0 + i1 - 19] - (0.5 * (h[-32 * ii0 + i1 - 19] - h[-32 * ii0 + i1 - 20])));
                                            }
                                          }
                                          #pragma omp parallel for
                                          for (int ii1 = ii0; ii1 <= ii0 + (N + 17) / 32; ii1 += 1) {
                                            if (T >= 32 * ii0 + 49 && T + N >= 32 * ii1 + 33) {
                                              for (int i1 = max(32 * ii0 + 20, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 19, 32 * ii1 + 33); i1 += 1) {
                                                h[-32 * ii0 + i1 - 20] = (h[-32 * ii0 + i1 - 20] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 19] - e[-32 * ii0 + i1 - 20])));
                                              }
                                            } else if (T >= 32 * ii0 + 32 && 32 * ii1 + 32 >= T + N) {
                                              for (int i1 = 32 * ii1 + 2; i1 <= N + 32 * ii0 + 19; i1 += 1) {
                                                h[-32 * ii0 + i1 - 20] = (h[-32 * ii0 + i1 - 20] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 19] - e[-32 * ii0 + i1 - 20])));
                                              }
                                            } else if (32 * ii0 + 31 >= T) {
                                              for (int i1 = max(32 * ii0 + 20, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 19, 32 * ii1 + 33); i1 += 1) {
                                                h[-32 * ii0 + i1 - 20] = (h[-32 * ii0 + i1 - 20] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 19] - e[-32 * ii0 + i1 - 20])));
                                              }
                                            } else {
                                              for (int i1 = max(32 * ii0 + 20, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 19, 32 * ii1 + 33); i1 += 1) {
                                                h[-32 * ii0 + i1 - 20] = (h[-32 * ii0 + i1 - 20] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 19] - e[-32 * ii0 + i1 - 20])));
                                              }
                                            }
                                          }
                                          if (T >= 32 * ii0 + 20) {
                                            #pragma omp parallel for
                                            for (int ii1 = ii0; ii1 <= ii0 + (N + 17) / 32; ii1 += 1) {
                                              for (int i1 = max(32 * ii0 + 21, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 19, 32 * ii1 + 33); i1 += 1) {
                                                e[-32 * ii0 + i1 - 20] = (e[-32 * ii0 + i1 - 20] - (0.5 * (h[-32 * ii0 + i1 - 20] - h[-32 * ii0 + i1 - 21])));
                                              }
                                            }
                                            #pragma omp parallel for
                                            for (int ii1 = ii0; ii1 <= ii0 + (N + 18) / 32; ii1 += 1) {
                                              if (T >= 32 * ii0 + 50 && T + N >= 32 * ii1 + 33) {
                                                for (int i1 = max(32 * ii0 + 21, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 20, 32 * ii1 + 33); i1 += 1) {
                                                  h[-32 * ii0 + i1 - 21] = (h[-32 * ii0 + i1 - 21] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 20] - e[-32 * ii0 + i1 - 21])));
                                                }
                                              } else if (T >= 32 * ii0 + 32 && 32 * ii1 + 32 >= T + N) {
                                                for (int i1 = 32 * ii1 + 2; i1 <= N + 32 * ii0 + 20; i1 += 1) {
                                                  h[-32 * ii0 + i1 - 21] = (h[-32 * ii0 + i1 - 21] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 20] - e[-32 * ii0 + i1 - 21])));
                                                }
                                              } else if (32 * ii0 + 31 >= T) {
                                                for (int i1 = max(32 * ii0 + 21, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 20, 32 * ii1 + 33); i1 += 1) {
                                                  h[-32 * ii0 + i1 - 21] = (h[-32 * ii0 + i1 - 21] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 20] - e[-32 * ii0 + i1 - 21])));
                                                }
                                              } else {
                                                for (int i1 = max(32 * ii0 + 21, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 20, 32 * ii1 + 33); i1 += 1) {
                                                  h[-32 * ii0 + i1 - 21] = (h[-32 * ii0 + i1 - 21] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 20] - e[-32 * ii0 + i1 - 21])));
                                                }
                                              }
                                            }
                                            if (T >= 32 * ii0 + 21) {
                                              #pragma omp parallel for
                                              for (int ii1 = ii0; ii1 <= ii0 + (N + 18) / 32; ii1 += 1) {
                                                for (int i1 = max(32 * ii0 + 22, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 20, 32 * ii1 + 33); i1 += 1) {
                                                  e[-32 * ii0 + i1 - 21] = (e[-32 * ii0 + i1 - 21] - (0.5 * (h[-32 * ii0 + i1 - 21] - h[-32 * ii0 + i1 - 22])));
                                                }
                                              }
                                              #pragma omp parallel for
                                              for (int ii1 = ii0; ii1 <= ii0 + (N + 19) / 32; ii1 += 1) {
                                                if (T >= 32 * ii0 + 51 && T + N >= 32 * ii1 + 33) {
                                                  for (int i1 = max(32 * ii0 + 22, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 21, 32 * ii1 + 33); i1 += 1) {
                                                    h[-32 * ii0 + i1 - 22] = (h[-32 * ii0 + i1 - 22] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 21] - e[-32 * ii0 + i1 - 22])));
                                                  }
                                                } else if (T >= 32 * ii0 + 32 && 32 * ii1 + 32 >= T + N) {
                                                  for (int i1 = 32 * ii1 + 2; i1 <= N + 32 * ii0 + 21; i1 += 1) {
                                                    h[-32 * ii0 + i1 - 22] = (h[-32 * ii0 + i1 - 22] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 21] - e[-32 * ii0 + i1 - 22])));
                                                  }
                                                } else if (32 * ii0 + 31 >= T) {
                                                  for (int i1 = max(32 * ii0 + 22, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 21, 32 * ii1 + 33); i1 += 1) {
                                                    h[-32 * ii0 + i1 - 22] = (h[-32 * ii0 + i1 - 22] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 21] - e[-32 * ii0 + i1 - 22])));
                                                  }
                                                } else {
                                                  for (int i1 = max(32 * ii0 + 22, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 21, 32 * ii1 + 33); i1 += 1) {
                                                    h[-32 * ii0 + i1 - 22] = (h[-32 * ii0 + i1 - 22] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 21] - e[-32 * ii0 + i1 - 22])));
                                                  }
                                                }
                                              }
                                              if (T >= 32 * ii0 + 22) {
                                                #pragma omp parallel for
                                                for (int ii1 = ii0; ii1 <= ii0 + (N + 19) / 32; ii1 += 1) {
                                                  for (int i1 = max(32 * ii0 + 23, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 21, 32 * ii1 + 33); i1 += 1) {
                                                    e[-32 * ii0 + i1 - 22] = (e[-32 * ii0 + i1 - 22] - (0.5 * (h[-32 * ii0 + i1 - 22] - h[-32 * ii0 + i1 - 23])));
                                                  }
                                                }
                                                #pragma omp parallel for
                                                for (int ii1 = ii0; ii1 <= ii0 + (N + 20) / 32; ii1 += 1) {
                                                  if (T >= 32 * ii0 + 52 && T + N >= 32 * ii1 + 33) {
                                                    for (int i1 = max(32 * ii0 + 23, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 22, 32 * ii1 + 33); i1 += 1) {
                                                      h[-32 * ii0 + i1 - 23] = (h[-32 * ii0 + i1 - 23] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 22] - e[-32 * ii0 + i1 - 23])));
                                                    }
                                                  } else if (T >= 32 * ii0 + 32 && 32 * ii1 + 32 >= T + N) {
                                                    for (int i1 = 32 * ii1 + 2; i1 <= N + 32 * ii0 + 22; i1 += 1) {
                                                      h[-32 * ii0 + i1 - 23] = (h[-32 * ii0 + i1 - 23] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 22] - e[-32 * ii0 + i1 - 23])));
                                                    }
                                                  } else if (32 * ii0 + 31 >= T) {
                                                    for (int i1 = max(32 * ii0 + 23, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 22, 32 * ii1 + 33); i1 += 1) {
                                                      h[-32 * ii0 + i1 - 23] = (h[-32 * ii0 + i1 - 23] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 22] - e[-32 * ii0 + i1 - 23])));
                                                    }
                                                  } else {
                                                    for (int i1 = max(32 * ii0 + 23, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 22, 32 * ii1 + 33); i1 += 1) {
                                                      h[-32 * ii0 + i1 - 23] = (h[-32 * ii0 + i1 - 23] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 22] - e[-32 * ii0 + i1 - 23])));
                                                    }
                                                  }
                                                }
                                                if (T >= 32 * ii0 + 23) {
                                                  #pragma omp parallel for
                                                  for (int ii1 = ii0; ii1 <= ii0 + (N + 20) / 32; ii1 += 1) {
                                                    for (int i1 = max(32 * ii0 + 24, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 22, 32 * ii1 + 33); i1 += 1) {
                                                      e[-32 * ii0 + i1 - 23] = (e[-32 * ii0 + i1 - 23] - (0.5 * (h[-32 * ii0 + i1 - 23] - h[-32 * ii0 + i1 - 24])));
                                                    }
                                                  }
                                                  #pragma omp parallel for
                                                  for (int ii1 = ii0; ii1 <= ii0 + (N + 21) / 32; ii1 += 1) {
                                                    if (T >= 32 * ii0 + 53 && T + N >= 32 * ii1 + 33) {
                                                      for (int i1 = max(32 * ii0 + 24, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 23, 32 * ii1 + 33); i1 += 1) {
                                                        h[-32 * ii0 + i1 - 24] = (h[-32 * ii0 + i1 - 24] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 23] - e[-32 * ii0 + i1 - 24])));
                                                      }
                                                    } else if (T >= 32 * ii0 + 32 && 32 * ii1 + 32 >= T + N) {
                                                      for (int i1 = 32 * ii1 + 2; i1 <= N + 32 * ii0 + 23; i1 += 1) {
                                                        h[-32 * ii0 + i1 - 24] = (h[-32 * ii0 + i1 - 24] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 23] - e[-32 * ii0 + i1 - 24])));
                                                      }
                                                    } else if (32 * ii0 + 31 >= T) {
                                                      for (int i1 = max(32 * ii0 + 24, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 23, 32 * ii1 + 33); i1 += 1) {
                                                        h[-32 * ii0 + i1 - 24] = (h[-32 * ii0 + i1 - 24] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 23] - e[-32 * ii0 + i1 - 24])));
                                                      }
                                                    } else {
                                                      for (int i1 = max(32 * ii0 + 24, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 23, 32 * ii1 + 33); i1 += 1) {
                                                        h[-32 * ii0 + i1 - 24] = (h[-32 * ii0 + i1 - 24] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 23] - e[-32 * ii0 + i1 - 24])));
                                                      }
                                                    }
                                                  }
                                                  if (T >= 32 * ii0 + 24) {
                                                    #pragma omp parallel for
                                                    for (int ii1 = ii0; ii1 <= ii0 + (N + 21) / 32; ii1 += 1) {
                                                      for (int i1 = max(32 * ii0 + 25, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 23, 32 * ii1 + 33); i1 += 1) {
                                                        e[-32 * ii0 + i1 - 24] = (e[-32 * ii0 + i1 - 24] - (0.5 * (h[-32 * ii0 + i1 - 24] - h[-32 * ii0 + i1 - 25])));
                                                      }
                                                    }
                                                    #pragma omp parallel for
                                                    for (int ii1 = ii0; ii1 <= ii0 + (N + 22) / 32; ii1 += 1) {
                                                      if (T >= 32 * ii0 + 54 && T + N >= 32 * ii1 + 33) {
                                                        for (int i1 = max(32 * ii0 + 25, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 24, 32 * ii1 + 33); i1 += 1) {
                                                          h[-32 * ii0 + i1 - 25] = (h[-32 * ii0 + i1 - 25] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 24] - e[-32 * ii0 + i1 - 25])));
                                                        }
                                                      } else if (T >= 32 * ii0 + 32 && 32 * ii1 + 32 >= T + N) {
                                                        for (int i1 = 32 * ii1 + 2; i1 <= N + 32 * ii0 + 24; i1 += 1) {
                                                          h[-32 * ii0 + i1 - 25] = (h[-32 * ii0 + i1 - 25] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 24] - e[-32 * ii0 + i1 - 25])));
                                                        }
                                                      } else if (32 * ii0 + 31 >= T) {
                                                        for (int i1 = max(32 * ii0 + 25, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 24, 32 * ii1 + 33); i1 += 1) {
                                                          h[-32 * ii0 + i1 - 25] = (h[-32 * ii0 + i1 - 25] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 24] - e[-32 * ii0 + i1 - 25])));
                                                        }
                                                      } else {
                                                        for (int i1 = max(32 * ii0 + 25, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 24, 32 * ii1 + 33); i1 += 1) {
                                                          h[-32 * ii0 + i1 - 25] = (h[-32 * ii0 + i1 - 25] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 24] - e[-32 * ii0 + i1 - 25])));
                                                        }
                                                      }
                                                    }
                                                    if (T >= 32 * ii0 + 25) {
                                                      #pragma omp parallel for
                                                      for (int ii1 = ii0; ii1 <= ii0 + (N + 22) / 32; ii1 += 1) {
                                                        for (int i1 = max(32 * ii0 + 26, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 24, 32 * ii1 + 33); i1 += 1) {
                                                          e[-32 * ii0 + i1 - 25] = (e[-32 * ii0 + i1 - 25] - (0.5 * (h[-32 * ii0 + i1 - 25] - h[-32 * ii0 + i1 - 26])));
                                                        }
                                                      }
                                                      #pragma omp parallel for
                                                      for (int ii1 = ii0; ii1 <= ii0 + (N + 23) / 32; ii1 += 1) {
                                                        if (T >= 32 * ii0 + 55 && T + N >= 32 * ii1 + 33) {
                                                          for (int i1 = max(32 * ii0 + 26, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 25, 32 * ii1 + 33); i1 += 1) {
                                                            h[-32 * ii0 + i1 - 26] = (h[-32 * ii0 + i1 - 26] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 25] - e[-32 * ii0 + i1 - 26])));
                                                          }
                                                        } else if (T >= 32 * ii0 + 32 && 32 * ii1 + 32 >= T + N) {
                                                          for (int i1 = 32 * ii1 + 2; i1 <= N + 32 * ii0 + 25; i1 += 1) {
                                                            h[-32 * ii0 + i1 - 26] = (h[-32 * ii0 + i1 - 26] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 25] - e[-32 * ii0 + i1 - 26])));
                                                          }
                                                        } else if (32 * ii0 + 31 >= T) {
                                                          for (int i1 = max(32 * ii0 + 26, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 25, 32 * ii1 + 33); i1 += 1) {
                                                            h[-32 * ii0 + i1 - 26] = (h[-32 * ii0 + i1 - 26] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 25] - e[-32 * ii0 + i1 - 26])));
                                                          }
                                                        } else {
                                                          for (int i1 = max(32 * ii0 + 26, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 25, 32 * ii1 + 33); i1 += 1) {
                                                            h[-32 * ii0 + i1 - 26] = (h[-32 * ii0 + i1 - 26] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 25] - e[-32 * ii0 + i1 - 26])));
                                                          }
                                                        }
                                                      }
                                                      if (T >= 32 * ii0 + 26) {
                                                        #pragma omp parallel for
                                                        for (int ii1 = ii0; ii1 <= ii0 + (N + 23) / 32; ii1 += 1) {
                                                          for (int i1 = max(32 * ii0 + 27, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 25, 32 * ii1 + 33); i1 += 1) {
                                                            e[-32 * ii0 + i1 - 26] = (e[-32 * ii0 + i1 - 26] - (0.5 * (h[-32 * ii0 + i1 - 26] - h[-32 * ii0 + i1 - 27])));
                                                          }
                                                        }
                                                        #pragma omp parallel for
                                                        for (int ii1 = ii0; ii1 <= ii0 + (N + 24) / 32; ii1 += 1) {
                                                          if (T >= 32 * ii0 + 56 && T + N >= 32 * ii1 + 33) {
                                                            for (int i1 = max(32 * ii0 + 27, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 26, 32 * ii1 + 33); i1 += 1) {
                                                              h[-32 * ii0 + i1 - 27] = (h[-32 * ii0 + i1 - 27] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 26] - e[-32 * ii0 + i1 - 27])));
                                                            }
                                                          } else if (T >= 32 * ii0 + 32 && 32 * ii1 + 32 >= T + N) {
                                                            for (int i1 = 32 * ii1 + 2; i1 <= N + 32 * ii0 + 26; i1 += 1) {
                                                              h[-32 * ii0 + i1 - 27] = (h[-32 * ii0 + i1 - 27] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 26] - e[-32 * ii0 + i1 - 27])));
                                                            }
                                                          } else if (32 * ii0 + 31 >= T) {
                                                            for (int i1 = max(32 * ii0 + 27, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 26, 32 * ii1 + 33); i1 += 1) {
                                                              h[-32 * ii0 + i1 - 27] = (h[-32 * ii0 + i1 - 27] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 26] - e[-32 * ii0 + i1 - 27])));
                                                            }
                                                          } else {
                                                            for (int i1 = max(32 * ii0 + 27, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 26, 32 * ii1 + 33); i1 += 1) {
                                                              h[-32 * ii0 + i1 - 27] = (h[-32 * ii0 + i1 - 27] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 26] - e[-32 * ii0 + i1 - 27])));
                                                            }
                                                          }
                                                        }
                                                        if (T >= 32 * ii0 + 27) {
                                                          #pragma omp parallel for
                                                          for (int ii1 = ii0; ii1 <= ii0 + (N + 24) / 32; ii1 += 1) {
                                                            for (int i1 = max(32 * ii0 + 28, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 26, 32 * ii1 + 33); i1 += 1) {
                                                              e[-32 * ii0 + i1 - 27] = (e[-32 * ii0 + i1 - 27] - (0.5 * (h[-32 * ii0 + i1 - 27] - h[-32 * ii0 + i1 - 28])));
                                                            }
                                                          }
                                                          #pragma omp parallel for
                                                          for (int ii1 = ii0; ii1 <= ii0 + (N + 25) / 32; ii1 += 1) {
                                                            if (T >= 32 * ii0 + 57 && T + N >= 32 * ii1 + 33) {
                                                              for (int i1 = max(32 * ii0 + 28, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 27, 32 * ii1 + 33); i1 += 1) {
                                                                h[-32 * ii0 + i1 - 28] = (h[-32 * ii0 + i1 - 28] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 27] - e[-32 * ii0 + i1 - 28])));
                                                              }
                                                            } else if (T >= 32 * ii0 + 32 && 32 * ii1 + 32 >= T + N) {
                                                              for (int i1 = 32 * ii1 + 2; i1 <= N + 32 * ii0 + 27; i1 += 1) {
                                                                h[-32 * ii0 + i1 - 28] = (h[-32 * ii0 + i1 - 28] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 27] - e[-32 * ii0 + i1 - 28])));
                                                              }
                                                            } else if (32 * ii0 + 31 >= T) {
                                                              for (int i1 = max(32 * ii0 + 28, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 27, 32 * ii1 + 33); i1 += 1) {
                                                                h[-32 * ii0 + i1 - 28] = (h[-32 * ii0 + i1 - 28] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 27] - e[-32 * ii0 + i1 - 28])));
                                                              }
                                                            } else {
                                                              for (int i1 = max(32 * ii0 + 28, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 27, 32 * ii1 + 33); i1 += 1) {
                                                                h[-32 * ii0 + i1 - 28] = (h[-32 * ii0 + i1 - 28] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 27] - e[-32 * ii0 + i1 - 28])));
                                                              }
                                                            }
                                                          }
                                                          if (T >= 32 * ii0 + 28) {
                                                            #pragma omp parallel for
                                                            for (int ii1 = ii0; ii1 <= ii0 + (N + 25) / 32; ii1 += 1) {
                                                              for (int i1 = max(32 * ii0 + 29, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 27, 32 * ii1 + 33); i1 += 1) {
                                                                e[-32 * ii0 + i1 - 28] = (e[-32 * ii0 + i1 - 28] - (0.5 * (h[-32 * ii0 + i1 - 28] - h[-32 * ii0 + i1 - 29])));
                                                              }
                                                            }
                                                            #pragma omp parallel for
                                                            for (int ii1 = ii0; ii1 <= ii0 + (N + 26) / 32; ii1 += 1) {
                                                              if (T >= 32 * ii0 + 58 && T + N >= 32 * ii1 + 33) {
                                                                for (int i1 = max(32 * ii0 + 29, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 28, 32 * ii1 + 33); i1 += 1) {
                                                                  h[-32 * ii0 + i1 - 29] = (h[-32 * ii0 + i1 - 29] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 28] - e[-32 * ii0 + i1 - 29])));
                                                                }
                                                              } else if (T >= 32 * ii0 + 32 && 32 * ii1 + 32 >= T + N) {
                                                                for (int i1 = 32 * ii1 + 2; i1 <= N + 32 * ii0 + 28; i1 += 1) {
                                                                  h[-32 * ii0 + i1 - 29] = (h[-32 * ii0 + i1 - 29] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 28] - e[-32 * ii0 + i1 - 29])));
                                                                }
                                                              } else if (32 * ii0 + 31 >= T) {
                                                                for (int i1 = max(32 * ii0 + 29, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 28, 32 * ii1 + 33); i1 += 1) {
                                                                  h[-32 * ii0 + i1 - 29] = (h[-32 * ii0 + i1 - 29] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 28] - e[-32 * ii0 + i1 - 29])));
                                                                }
                                                              } else {
                                                                for (int i1 = max(32 * ii0 + 29, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 28, 32 * ii1 + 33); i1 += 1) {
                                                                  h[-32 * ii0 + i1 - 29] = (h[-32 * ii0 + i1 - 29] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 28] - e[-32 * ii0 + i1 - 29])));
                                                                }
                                                              }
                                                            }
                                                            if (T >= 32 * ii0 + 29) {
                                                              #pragma omp parallel for
                                                              for (int ii1 = ii0; ii1 <= ii0 + (N + 26) / 32; ii1 += 1) {
                                                                for (int i1 = max(32 * ii0 + 30, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 28, 32 * ii1 + 33); i1 += 1) {
                                                                  e[-32 * ii0 + i1 - 29] = (e[-32 * ii0 + i1 - 29] - (0.5 * (h[-32 * ii0 + i1 - 29] - h[-32 * ii0 + i1 - 30])));
                                                                }
                                                              }
                                                              #pragma omp parallel for
                                                              for (int ii1 = ii0; ii1 <= ii0 + (N + 27) / 32; ii1 += 1) {
                                                                if (T >= 32 * ii0 + 59 && T + N >= 32 * ii1 + 33) {
                                                                  for (int i1 = max(32 * ii0 + 30, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 29, 32 * ii1 + 33); i1 += 1) {
                                                                    h[-32 * ii0 + i1 - 30] = (h[-32 * ii0 + i1 - 30] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 29] - e[-32 * ii0 + i1 - 30])));
                                                                  }
                                                                } else if (T >= 32 * ii0 + 32 && 32 * ii1 + 32 >= T + N) {
                                                                  for (int i1 = 32 * ii1 + 2; i1 <= N + 32 * ii0 + 29; i1 += 1) {
                                                                    h[-32 * ii0 + i1 - 30] = (h[-32 * ii0 + i1 - 30] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 29] - e[-32 * ii0 + i1 - 30])));
                                                                  }
                                                                } else if (32 * ii0 + 31 >= T) {
                                                                  for (int i1 = max(32 * ii0 + 30, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 29, 32 * ii1 + 33); i1 += 1) {
                                                                    h[-32 * ii0 + i1 - 30] = (h[-32 * ii0 + i1 - 30] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 29] - e[-32 * ii0 + i1 - 30])));
                                                                  }
                                                                } else {
                                                                  for (int i1 = max(32 * ii0 + 30, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 29, 32 * ii1 + 33); i1 += 1) {
                                                                    h[-32 * ii0 + i1 - 30] = (h[-32 * ii0 + i1 - 30] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 29] - e[-32 * ii0 + i1 - 30])));
                                                                  }
                                                                }
                                                              }
                                                              if (T >= 32 * ii0 + 30) {
                                                                #pragma omp parallel for
                                                                for (int ii1 = ii0; ii1 <= ii0 + (N + 27) / 32; ii1 += 1) {
                                                                  for (int i1 = max(32 * ii0 + 31, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 29, 32 * ii1 + 33); i1 += 1) {
                                                                    e[-32 * ii0 + i1 - 30] = (e[-32 * ii0 + i1 - 30] - (0.5 * (h[-32 * ii0 + i1 - 30] - h[-32 * ii0 + i1 - 31])));
                                                                  }
                                                                }
                                                                #pragma omp parallel for
                                                                for (int ii1 = ii0; ii1 <= ii0 + (N + 28) / 32; ii1 += 1) {
                                                                  if (T >= 32 * ii0 + 60 && T + N >= 32 * ii1 + 33) {
                                                                    for (int i1 = max(32 * ii0 + 31, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 30, 32 * ii1 + 33); i1 += 1) {
                                                                      h[-32 * ii0 + i1 - 31] = (h[-32 * ii0 + i1 - 31] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 30] - e[-32 * ii0 + i1 - 31])));
                                                                    }
                                                                  } else if (T >= 32 * ii0 + 32 && 32 * ii1 + 32 >= T + N) {
                                                                    for (int i1 = 32 * ii1 + 2; i1 <= N + 32 * ii0 + 30; i1 += 1) {
                                                                      h[-32 * ii0 + i1 - 31] = (h[-32 * ii0 + i1 - 31] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 30] - e[-32 * ii0 + i1 - 31])));
                                                                    }
                                                                  } else if (32 * ii0 + 31 >= T) {
                                                                    for (int i1 = max(32 * ii0 + 31, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 30, 32 * ii1 + 33); i1 += 1) {
                                                                      h[-32 * ii0 + i1 - 31] = (h[-32 * ii0 + i1 - 31] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 30] - e[-32 * ii0 + i1 - 31])));
                                                                    }
                                                                  } else {
                                                                    for (int i1 = max(32 * ii0 + 31, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 30, 32 * ii1 + 33); i1 += 1) {
                                                                      h[-32 * ii0 + i1 - 31] = (h[-32 * ii0 + i1 - 31] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 30] - e[-32 * ii0 + i1 - 31])));
                                                                    }
                                                                  }
                                                                }
                                                                if (T >= 32 * ii0 + 31) {
                                                                  #pragma omp parallel for
                                                                  for (int ii1 = ii0; ii1 <= ii0 + (N + 28) / 32; ii1 += 1) {
                                                                    for (int i1 = max(32 * ii0 + 32, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 30, 32 * ii1 + 33); i1 += 1) {
                                                                      e[-32 * ii0 + i1 - 31] = (e[-32 * ii0 + i1 - 31] - (0.5 * (h[-32 * ii0 + i1 - 31] - h[-32 * ii0 + i1 - 32])));
                                                                    }
                                                                  }
                                                                  #pragma omp parallel for
                                                                  for (int ii1 = ii0; ii1 <= ii0 + (N + 29) / 32; ii1 += 1) {
                                                                    if (T >= 32 * ii0 + 61 && T + N >= 32 * ii1 + 33) {
                                                                      for (int i1 = max(32 * ii0 + 32, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 31, 32 * ii1 + 33); i1 += 1) {
                                                                        h[-32 * ii0 + i1 - 32] = (h[-32 * ii0 + i1 - 32] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 31] - e[-32 * ii0 + i1 - 32])));
                                                                      }
                                                                    } else if (T >= 32 * ii0 + 32 && 32 * ii1 + 32 >= T + N) {
                                                                      for (int i1 = 32 * ii1 + 2; i1 <= N + 32 * ii0 + 31; i1 += 1) {
                                                                        h[-32 * ii0 + i1 - 32] = (h[-32 * ii0 + i1 - 32] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 31] - e[-32 * ii0 + i1 - 32])));
                                                                      }
                                                                    } else if (T >= 32 * ii0 + 32 && T + N >= 32 * ii1 + 33) {
                                                                      for (int i1 = max(32 * ii0 + 32, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 31, 32 * ii1 + 33); i1 += 1) {
                                                                        h[-32 * ii0 + i1 - 32] = (h[-32 * ii0 + i1 - 32] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 31] - e[-32 * ii0 + i1 - 32])));
                                                                      }
                                                                    } else {
                                                                      for (int i1 = max(T + 1, 32 * ii1 + 2); i1 <= min(T + N, 32 * ii1 + 33); i1 += 1) {
                                                                        h[-T + i1 - 1] = (h[-T + i1 - 1] - (0.69999999999999996 * (e[-T + i1] - e[-T + i1 - 1])));
                                                                      }
                                                                    }
                                                                  }
                                                                  if (T >= 32 * ii0 + 32) {
                                                                    #pragma omp parallel for
                                                                    for (int ii1 = ii0; ii1 <= ii0 + (N + 29) / 32; ii1 += 1) {
                                                                      for (int i1 = max(32 * ii0 + 33, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 31, 32 * ii1 + 33); i1 += 1) {
                                                                        e[-32 * ii0 + i1 - 32] = (e[-32 * ii0 + i1 - 32] - (0.5 * (h[-32 * ii0 + i1 - 32] - h[-32 * ii0 + i1 - 33])));
                                                                      }
                                                                    }
                                                                    #pragma omp parallel for
                                                                    for (int ii1 = ii0; ii1 <= ii0 + (N - 2) / 32 + 1; ii1 += 1) {
                                                                      for (int i1 = max(32 * ii0 + 33, 32 * ii1 + 2); i1 <= min(N + 32 * ii0 + 32, 32 * ii1 + 33); i1 += 1) {
                                                                        h[-32 * ii0 + i1 - 33] = (h[-32 * ii0 + i1 - 33] - (0.69999999999999996 * (e[-32 * ii0 + i1 - 32] - e[-32 * ii0 + i1 - 33])));
                                                                      }
                                                                    }
                                                                  }
                                                                }
                                                              }
                                                            }
                                                          }
                                                        }
                                                      }
                                                    }
                                                  }
                                                }
                                              }
                                            }
                                          }
                                        }
                                      }
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    } else {
      for (int i0 = 32 * ii0 + 1; i0 <= min(T, 32 * ii0 + 32); i0 += 1) {
        h[0] = (h[0] - (0.69999999999999996 * (e[1] - e[0])));
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
