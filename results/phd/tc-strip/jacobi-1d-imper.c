#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "util.h"

#define N 2000000
#define T 1000

#pragma declarations
double a[N];
double b[N];
#pragma enddeclarations

double t_start, t_end;

void init_array()
{
    int j;

    for (j=0; j<N; j++) {
        a[j] = ((double)j)/N;
    }
}


void print_array()
{
    int j;

    for (j=0; j<N; j++) {
        fprintf(stderr, "%lf ", a[j]);
        if (j%80 == 20) fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
}


int main()
{
    int t, i, j;

    init_array();

    IF_TIME(t_start = rtclock());

/* TC Optimizing Compiler 0.4.2 */
/* ./tc ../examples/pluto-perfect/jacobi-1d-imper.scop.c --strip-tiling --omp-for-codegen --isl-map-tc --inline --debug -b 32 --drop-bounds */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#pragma scop
if (N >= 4) {
  for (int ii0 = 0; ii0 <= floord(T - 1, 32); ii0 += 1) {
    {
      #pragma omp parallel for
      for (int ii1 = 2 * ii0; ii1 <= 2 * ii0 + floord(N - 4, 32); ii1 += 1) {
        if (ii1 >= 2 * ii0 + 1) {
          for (int i1 = 32 * ii1; i1 <= min(N + 64 * ii0 - 4, 32 * ii1 + 31); i1 += 1) {
            b[-64 * ii0 + i1 + 2] = (0.33333 * ((a[-64 * ii0 + i1 + 1] + a[-64 * ii0 + i1 + 2]) + a[-64 * ii0 + i1 + 3]));
          }
        } else {
          for (int i1 = 64 * ii0; i1 <= min(N + 64 * ii0 - 4, 64 * ii0 + 31); i1 += 1) {
            b[-64 * ii0 + i1 + 2] = (0.33333 * ((a[-64 * ii0 + i1 + 1] + a[-64 * ii0 + i1 + 2]) + a[-64 * ii0 + i1 + 3]));
          }
        }
      }
      #pragma omp parallel for
      for (int ii1 = 2 * ii0; ii1 <= 2 * ii0 + floord(N - 3, 32); ii1 += 1) {
        for (int i1 = max(64 * ii0 + 1, 32 * ii1); i1 <= min(N + 64 * ii0 - 3, 32 * ii1 + 31); i1 += 1) {
          a[-64 * ii0 + i1 + 1] = b[-64 * ii0 + i1 + 1];
        }
      }
    }
    if (T >= 32 * ii0 + 2) {
      #pragma omp parallel for
      for (int ii1 = 2 * ii0; ii1 <= 2 * ii0 + (N - 2) / 32; ii1 += 1) {
        if (ii1 >= 2 * ii0 + 1) {
          for (int i1 = 32 * ii1; i1 <= min(N + 64 * ii0 - 2, 32 * ii1 + 31); i1 += 1) {
            b[-64 * ii0 + i1] = (0.33333 * ((a[-64 * ii0 + i1 - 1] + a[-64 * ii0 + i1]) + a[-64 * ii0 + i1 + 1]));
          }
        } else {
          for (int i1 = 64 * ii0 + 2; i1 <= min(N + 64 * ii0 - 2, 64 * ii0 + 31); i1 += 1) {
            b[-64 * ii0 + i1] = (0.33333 * ((a[-64 * ii0 + i1 - 1] + a[-64 * ii0 + i1]) + a[-64 * ii0 + i1 + 1]));
          }
        }
      }
      #pragma omp parallel for
      for (int ii1 = 2 * ii0; ii1 <= 2 * ii0 + floord(N - 1, 32); ii1 += 1) {
        for (int i1 = max(64 * ii0 + 3, 32 * ii1); i1 <= min(N + 64 * ii0 - 1, 32 * ii1 + 31); i1 += 1) {
          a[-64 * ii0 + i1 - 1] = b[-64 * ii0 + i1 - 1];
        }
      }
      if (T >= 32 * ii0 + 3) {
        #pragma omp parallel for
        for (int ii1 = 2 * ii0; ii1 <= 2 * ii0 + N / 32; ii1 += 1) {
          if (ii1 >= 2 * ii0 + 1) {
            for (int i1 = 32 * ii1; i1 <= min(N + 64 * ii0, 32 * ii1 + 31); i1 += 1) {
              b[-64 * ii0 + i1 - 2] = (0.33333 * ((a[-64 * ii0 + i1 - 3] + a[-64 * ii0 + i1 - 2]) + a[-64 * ii0 + i1 - 1]));
            }
          } else {
            for (int i1 = 64 * ii0 + 4; i1 <= min(N + 64 * ii0, 64 * ii0 + 31); i1 += 1) {
              b[-64 * ii0 + i1 - 2] = (0.33333 * ((a[-64 * ii0 + i1 - 3] + a[-64 * ii0 + i1 - 2]) + a[-64 * ii0 + i1 - 1]));
            }
          }
        }
        #pragma omp parallel for
        for (int ii1 = 2 * ii0; ii1 <= 2 * ii0 + (N + 1) / 32; ii1 += 1) {
          for (int i1 = max(64 * ii0 + 5, 32 * ii1); i1 <= min(N + 64 * ii0 + 1, 32 * ii1 + 31); i1 += 1) {
            a[-64 * ii0 + i1 - 3] = b[-64 * ii0 + i1 - 3];
          }
        }
        if (T >= 32 * ii0 + 4) {
          #pragma omp parallel for
          for (int ii1 = 2 * ii0; ii1 <= 2 * ii0 + (N + 2) / 32; ii1 += 1) {
            if (ii1 >= 2 * ii0 + 1) {
              for (int i1 = 32 * ii1; i1 <= min(N + 64 * ii0 + 2, 32 * ii1 + 31); i1 += 1) {
                b[-64 * ii0 + i1 - 4] = (0.33333 * ((a[-64 * ii0 + i1 - 5] + a[-64 * ii0 + i1 - 4]) + a[-64 * ii0 + i1 - 3]));
              }
            } else {
              for (int i1 = 64 * ii0 + 6; i1 <= min(N + 64 * ii0 + 2, 64 * ii0 + 31); i1 += 1) {
                b[-64 * ii0 + i1 - 4] = (0.33333 * ((a[-64 * ii0 + i1 - 5] + a[-64 * ii0 + i1 - 4]) + a[-64 * ii0 + i1 - 3]));
              }
            }
          }
          #pragma omp parallel for
          for (int ii1 = 2 * ii0; ii1 <= 2 * ii0 + (N + 3) / 32; ii1 += 1) {
            for (int i1 = max(64 * ii0 + 7, 32 * ii1); i1 <= min(N + 64 * ii0 + 3, 32 * ii1 + 31); i1 += 1) {
              a[-64 * ii0 + i1 - 5] = b[-64 * ii0 + i1 - 5];
            }
          }
          if (T >= 32 * ii0 + 5) {
            #pragma omp parallel for
            for (int ii1 = 2 * ii0; ii1 <= 2 * ii0 + (N + 4) / 32; ii1 += 1) {
              if (ii1 >= 2 * ii0 + 1) {
                for (int i1 = 32 * ii1; i1 <= min(N + 64 * ii0 + 4, 32 * ii1 + 31); i1 += 1) {
                  b[-64 * ii0 + i1 - 6] = (0.33333 * ((a[-64 * ii0 + i1 - 7] + a[-64 * ii0 + i1 - 6]) + a[-64 * ii0 + i1 - 5]));
                }
              } else {
                for (int i1 = 64 * ii0 + 8; i1 <= min(N + 64 * ii0 + 4, 64 * ii0 + 31); i1 += 1) {
                  b[-64 * ii0 + i1 - 6] = (0.33333 * ((a[-64 * ii0 + i1 - 7] + a[-64 * ii0 + i1 - 6]) + a[-64 * ii0 + i1 - 5]));
                }
              }
            }
            #pragma omp parallel for
            for (int ii1 = 2 * ii0; ii1 <= 2 * ii0 + (N + 5) / 32; ii1 += 1) {
              for (int i1 = max(64 * ii0 + 9, 32 * ii1); i1 <= min(N + 64 * ii0 + 5, 32 * ii1 + 31); i1 += 1) {
                a[-64 * ii0 + i1 - 7] = b[-64 * ii0 + i1 - 7];
              }
            }
            if (T >= 32 * ii0 + 6) {
              #pragma omp parallel for
              for (int ii1 = 2 * ii0; ii1 <= 2 * ii0 + (N + 6) / 32; ii1 += 1) {
                if (ii1 >= 2 * ii0 + 1) {
                  for (int i1 = 32 * ii1; i1 <= min(N + 64 * ii0 + 6, 32 * ii1 + 31); i1 += 1) {
                    b[-64 * ii0 + i1 - 8] = (0.33333 * ((a[-64 * ii0 + i1 - 9] + a[-64 * ii0 + i1 - 8]) + a[-64 * ii0 + i1 - 7]));
                  }
                } else {
                  for (int i1 = 64 * ii0 + 10; i1 <= min(N + 64 * ii0 + 6, 64 * ii0 + 31); i1 += 1) {
                    b[-64 * ii0 + i1 - 8] = (0.33333 * ((a[-64 * ii0 + i1 - 9] + a[-64 * ii0 + i1 - 8]) + a[-64 * ii0 + i1 - 7]));
                  }
                }
              }
              #pragma omp parallel for
              for (int ii1 = 2 * ii0; ii1 <= 2 * ii0 + (N + 7) / 32; ii1 += 1) {
                for (int i1 = max(64 * ii0 + 11, 32 * ii1); i1 <= min(N + 64 * ii0 + 7, 32 * ii1 + 31); i1 += 1) {
                  a[-64 * ii0 + i1 - 9] = b[-64 * ii0 + i1 - 9];
                }
              }
              if (T >= 32 * ii0 + 7) {
                #pragma omp parallel for
                for (int ii1 = 2 * ii0; ii1 <= 2 * ii0 + (N + 8) / 32; ii1 += 1) {
                  if (ii1 >= 2 * ii0 + 1) {
                    for (int i1 = 32 * ii1; i1 <= min(N + 64 * ii0 + 8, 32 * ii1 + 31); i1 += 1) {
                      b[-64 * ii0 + i1 - 10] = (0.33333 * ((a[-64 * ii0 + i1 - 11] + a[-64 * ii0 + i1 - 10]) + a[-64 * ii0 + i1 - 9]));
                    }
                  } else {
                    for (int i1 = 64 * ii0 + 12; i1 <= min(N + 64 * ii0 + 8, 64 * ii0 + 31); i1 += 1) {
                      b[-64 * ii0 + i1 - 10] = (0.33333 * ((a[-64 * ii0 + i1 - 11] + a[-64 * ii0 + i1 - 10]) + a[-64 * ii0 + i1 - 9]));
                    }
                  }
                }
                #pragma omp parallel for
                for (int ii1 = 2 * ii0; ii1 <= 2 * ii0 + (N + 9) / 32; ii1 += 1) {
                  for (int i1 = max(64 * ii0 + 13, 32 * ii1); i1 <= min(N + 64 * ii0 + 9, 32 * ii1 + 31); i1 += 1) {
                    a[-64 * ii0 + i1 - 11] = b[-64 * ii0 + i1 - 11];
                  }
                }
                if (T >= 32 * ii0 + 8) {
                  #pragma omp parallel for
                  for (int ii1 = 2 * ii0; ii1 <= 2 * ii0 + (N + 10) / 32; ii1 += 1) {
                    if (ii1 >= 2 * ii0 + 1) {
                      for (int i1 = 32 * ii1; i1 <= min(N + 64 * ii0 + 10, 32 * ii1 + 31); i1 += 1) {
                        b[-64 * ii0 + i1 - 12] = (0.33333 * ((a[-64 * ii0 + i1 - 13] + a[-64 * ii0 + i1 - 12]) + a[-64 * ii0 + i1 - 11]));
                      }
                    } else {
                      for (int i1 = 64 * ii0 + 14; i1 <= min(N + 64 * ii0 + 10, 64 * ii0 + 31); i1 += 1) {
                        b[-64 * ii0 + i1 - 12] = (0.33333 * ((a[-64 * ii0 + i1 - 13] + a[-64 * ii0 + i1 - 12]) + a[-64 * ii0 + i1 - 11]));
                      }
                    }
                  }
                  #pragma omp parallel for
                  for (int ii1 = 2 * ii0; ii1 <= 2 * ii0 + (N + 11) / 32; ii1 += 1) {
                    for (int i1 = max(64 * ii0 + 15, 32 * ii1); i1 <= min(N + 64 * ii0 + 11, 32 * ii1 + 31); i1 += 1) {
                      a[-64 * ii0 + i1 - 13] = b[-64 * ii0 + i1 - 13];
                    }
                  }
                  if (T >= 32 * ii0 + 9) {
                    #pragma omp parallel for
                    for (int ii1 = 2 * ii0; ii1 <= 2 * ii0 + (N + 12) / 32; ii1 += 1) {
                      for (int i1 = max(64 * ii0 + 16, 32 * ii1); i1 <= min(N + 64 * ii0 + 12, 32 * ii1 + 31); i1 += 1) {
                        b[-64 * ii0 + i1 - 14] = (0.33333 * ((a[-64 * ii0 + i1 - 15] + a[-64 * ii0 + i1 - 14]) + a[-64 * ii0 + i1 - 13]));
                      }
                    }
                    #pragma omp parallel for
                    for (int ii1 = 2 * ii0; ii1 <= 2 * ii0 + (N + 13) / 32; ii1 += 1) {
                      for (int i1 = max(64 * ii0 + 17, 32 * ii1); i1 <= min(N + 64 * ii0 + 13, 32 * ii1 + 31); i1 += 1) {
                        a[-64 * ii0 + i1 - 15] = b[-64 * ii0 + i1 - 15];
                      }
                    }
                    if (T >= 32 * ii0 + 10) {
                      #pragma omp parallel for
                      for (int ii1 = 2 * ii0; ii1 <= 2 * ii0 + (N + 14) / 32; ii1 += 1) {
                        for (int i1 = max(64 * ii0 + 18, 32 * ii1); i1 <= min(N + 64 * ii0 + 14, 32 * ii1 + 31); i1 += 1) {
                          b[-64 * ii0 + i1 - 16] = (0.33333 * ((a[-64 * ii0 + i1 - 17] + a[-64 * ii0 + i1 - 16]) + a[-64 * ii0 + i1 - 15]));
                        }
                      }
                      #pragma omp parallel for
                      for (int ii1 = 2 * ii0; ii1 <= 2 * ii0 + (N + 15) / 32; ii1 += 1) {
                        for (int i1 = max(64 * ii0 + 19, 32 * ii1); i1 <= min(N + 64 * ii0 + 15, 32 * ii1 + 31); i1 += 1) {
                          a[-64 * ii0 + i1 - 17] = b[-64 * ii0 + i1 - 17];
                        }
                      }
                      if (T >= 32 * ii0 + 11) {
                        #pragma omp parallel for
                        for (int ii1 = 2 * ii0; ii1 <= 2 * ii0 + (N + 16) / 32; ii1 += 1) {
                          for (int i1 = max(64 * ii0 + 20, 32 * ii1); i1 <= min(N + 64 * ii0 + 16, 32 * ii1 + 31); i1 += 1) {
                            b[-64 * ii0 + i1 - 18] = (0.33333 * ((a[-64 * ii0 + i1 - 19] + a[-64 * ii0 + i1 - 18]) + a[-64 * ii0 + i1 - 17]));
                          }
                        }
                        #pragma omp parallel for
                        for (int ii1 = 2 * ii0; ii1 <= 2 * ii0 + (N + 17) / 32; ii1 += 1) {
                          for (int i1 = max(64 * ii0 + 21, 32 * ii1); i1 <= min(N + 64 * ii0 + 17, 32 * ii1 + 31); i1 += 1) {
                            a[-64 * ii0 + i1 - 19] = b[-64 * ii0 + i1 - 19];
                          }
                        }
                        if (T >= 32 * ii0 + 12) {
                          #pragma omp parallel for
                          for (int ii1 = 2 * ii0; ii1 <= 2 * ii0 + (N + 18) / 32; ii1 += 1) {
                            for (int i1 = max(64 * ii0 + 22, 32 * ii1); i1 <= min(N + 64 * ii0 + 18, 32 * ii1 + 31); i1 += 1) {
                              b[-64 * ii0 + i1 - 20] = (0.33333 * ((a[-64 * ii0 + i1 - 21] + a[-64 * ii0 + i1 - 20]) + a[-64 * ii0 + i1 - 19]));
                            }
                          }
                          #pragma omp parallel for
                          for (int ii1 = 2 * ii0; ii1 <= 2 * ii0 + (N + 19) / 32; ii1 += 1) {
                            for (int i1 = max(64 * ii0 + 23, 32 * ii1); i1 <= min(N + 64 * ii0 + 19, 32 * ii1 + 31); i1 += 1) {
                              a[-64 * ii0 + i1 - 21] = b[-64 * ii0 + i1 - 21];
                            }
                          }
                          if (T >= 32 * ii0 + 13) {
                            #pragma omp parallel for
                            for (int ii1 = 2 * ii0; ii1 <= 2 * ii0 + (N + 20) / 32; ii1 += 1) {
                              for (int i1 = max(64 * ii0 + 24, 32 * ii1); i1 <= min(N + 64 * ii0 + 20, 32 * ii1 + 31); i1 += 1) {
                                b[-64 * ii0 + i1 - 22] = (0.33333 * ((a[-64 * ii0 + i1 - 23] + a[-64 * ii0 + i1 - 22]) + a[-64 * ii0 + i1 - 21]));
                              }
                            }
                            #pragma omp parallel for
                            for (int ii1 = 2 * ii0; ii1 <= 2 * ii0 + (N + 21) / 32; ii1 += 1) {
                              for (int i1 = max(64 * ii0 + 25, 32 * ii1); i1 <= min(N + 64 * ii0 + 21, 32 * ii1 + 31); i1 += 1) {
                                a[-64 * ii0 + i1 - 23] = b[-64 * ii0 + i1 - 23];
                              }
                            }
                            if (T >= 32 * ii0 + 14) {
                              #pragma omp parallel for
                              for (int ii1 = 2 * ii0; ii1 <= 2 * ii0 + (N + 22) / 32; ii1 += 1) {
                                for (int i1 = max(64 * ii0 + 26, 32 * ii1); i1 <= min(N + 64 * ii0 + 22, 32 * ii1 + 31); i1 += 1) {
                                  b[-64 * ii0 + i1 - 24] = (0.33333 * ((a[-64 * ii0 + i1 - 25] + a[-64 * ii0 + i1 - 24]) + a[-64 * ii0 + i1 - 23]));
                                }
                              }
                              #pragma omp parallel for
                              for (int ii1 = 2 * ii0; ii1 <= 2 * ii0 + (N + 23) / 32; ii1 += 1) {
                                for (int i1 = max(64 * ii0 + 27, 32 * ii1); i1 <= min(N + 64 * ii0 + 23, 32 * ii1 + 31); i1 += 1) {
                                  a[-64 * ii0 + i1 - 25] = b[-64 * ii0 + i1 - 25];
                                }
                              }
                              if (T >= 32 * ii0 + 15) {
                                #pragma omp parallel for
                                for (int ii1 = 2 * ii0; ii1 <= 2 * ii0 + (N + 24) / 32; ii1 += 1) {
                                  for (int i1 = max(64 * ii0 + 28, 32 * ii1); i1 <= min(N + 64 * ii0 + 24, 32 * ii1 + 31); i1 += 1) {
                                    b[-64 * ii0 + i1 - 26] = (0.33333 * ((a[-64 * ii0 + i1 - 27] + a[-64 * ii0 + i1 - 26]) + a[-64 * ii0 + i1 - 25]));
                                  }
                                }
                                #pragma omp parallel for
                                for (int ii1 = 2 * ii0; ii1 <= 2 * ii0 + (N + 25) / 32; ii1 += 1) {
                                  for (int i1 = max(64 * ii0 + 29, 32 * ii1); i1 <= min(N + 64 * ii0 + 25, 32 * ii1 + 31); i1 += 1) {
                                    a[-64 * ii0 + i1 - 27] = b[-64 * ii0 + i1 - 27];
                                  }
                                }
                                if (T >= 32 * ii0 + 16) {
                                  #pragma omp parallel for
                                  for (int ii1 = 2 * ii0; ii1 <= 2 * ii0 + (N + 26) / 32; ii1 += 1) {
                                    for (int i1 = max(64 * ii0 + 30, 32 * ii1); i1 <= min(N + 64 * ii0 + 26, 32 * ii1 + 31); i1 += 1) {
                                      b[-64 * ii0 + i1 - 28] = (0.33333 * ((a[-64 * ii0 + i1 - 29] + a[-64 * ii0 + i1 - 28]) + a[-64 * ii0 + i1 - 27]));
                                    }
                                  }
                                  #pragma omp parallel for
                                  for (int ii1 = 2 * ii0; ii1 <= 2 * ii0 + (N + 27) / 32; ii1 += 1) {
                                    for (int i1 = max(64 * ii0 + 31, 32 * ii1); i1 <= min(N + 64 * ii0 + 27, 32 * ii1 + 31); i1 += 1) {
                                      a[-64 * ii0 + i1 - 29] = b[-64 * ii0 + i1 - 29];
                                    }
                                  }
                                  if (T >= 32 * ii0 + 17) {
                                    #pragma omp parallel for
                                    for (int ii1 = 2 * ii0 + 1; ii1 <= 2 * ii0 + (N + 28) / 32; ii1 += 1) {
                                      if (ii1 >= 2 * ii0 + 2) {
                                        for (int i1 = 32 * ii1; i1 <= min(N + 64 * ii0 + 28, 32 * ii1 + 31); i1 += 1) {
                                          b[-64 * ii0 + i1 - 30] = (0.33333 * ((a[-64 * ii0 + i1 - 31] + a[-64 * ii0 + i1 - 30]) + a[-64 * ii0 + i1 - 29]));
                                        }
                                      } else {
                                        for (int i1 = 64 * ii0 + 32; i1 <= min(N + 64 * ii0 + 28, 64 * ii0 + 63); i1 += 1) {
                                          b[-64 * ii0 + i1 - 30] = (0.33333 * ((a[-64 * ii0 + i1 - 31] + a[-64 * ii0 + i1 - 30]) + a[-64 * ii0 + i1 - 29]));
                                        }
                                      }
                                    }
                                    #pragma omp parallel for
                                    for (int ii1 = 2 * ii0 + 1; ii1 <= 2 * ii0 + (N + 29) / 32; ii1 += 1) {
                                      for (int i1 = max(64 * ii0 + 33, 32 * ii1); i1 <= min(N + 64 * ii0 + 29, 32 * ii1 + 31); i1 += 1) {
                                        a[-64 * ii0 + i1 - 31] = b[-64 * ii0 + i1 - 31];
                                      }
                                    }
                                    if (T >= 32 * ii0 + 18) {
                                      #pragma omp parallel for
                                      for (int ii1 = 2 * ii0 + 1; ii1 <= 2 * ii0 + (N - 2) / 32 + 1; ii1 += 1) {
                                        if (ii1 >= 2 * ii0 + 2) {
                                          for (int i1 = 32 * ii1; i1 <= min(N + 64 * ii0 + 30, 32 * ii1 + 31); i1 += 1) {
                                            b[-64 * ii0 + i1 - 32] = (0.33333 * ((a[-64 * ii0 + i1 - 33] + a[-64 * ii0 + i1 - 32]) + a[-64 * ii0 + i1 - 31]));
                                          }
                                        } else {
                                          for (int i1 = 64 * ii0 + 34; i1 <= min(N + 64 * ii0 + 30, 64 * ii0 + 63); i1 += 1) {
                                            b[-64 * ii0 + i1 - 32] = (0.33333 * ((a[-64 * ii0 + i1 - 33] + a[-64 * ii0 + i1 - 32]) + a[-64 * ii0 + i1 - 31]));
                                          }
                                        }
                                      }
                                      #pragma omp parallel for
                                      for (int ii1 = 2 * ii0 + 1; ii1 <= 2 * ii0 + (N + 31) / 32; ii1 += 1) {
                                        for (int i1 = max(64 * ii0 + 35, 32 * ii1); i1 <= min(N + 64 * ii0 + 31, 32 * ii1 + 31); i1 += 1) {
                                          a[-64 * ii0 + i1 - 33] = b[-64 * ii0 + i1 - 33];
                                        }
                                      }
                                      if (T >= 32 * ii0 + 19) {
                                        #pragma omp parallel for
                                        for (int ii1 = 2 * ii0 + 1; ii1 <= 2 * ii0 + N / 32 + 1; ii1 += 1) {
                                          if (ii1 >= 2 * ii0 + 2) {
                                            for (int i1 = 32 * ii1; i1 <= min(N + 64 * ii0 + 32, 32 * ii1 + 31); i1 += 1) {
                                              b[-64 * ii0 + i1 - 34] = (0.33333 * ((a[-64 * ii0 + i1 - 35] + a[-64 * ii0 + i1 - 34]) + a[-64 * ii0 + i1 - 33]));
                                            }
                                          } else {
                                            for (int i1 = 64 * ii0 + 36; i1 <= min(N + 64 * ii0 + 32, 64 * ii0 + 63); i1 += 1) {
                                              b[-64 * ii0 + i1 - 34] = (0.33333 * ((a[-64 * ii0 + i1 - 35] + a[-64 * ii0 + i1 - 34]) + a[-64 * ii0 + i1 - 33]));
                                            }
                                          }
                                        }
                                        #pragma omp parallel for
                                        for (int ii1 = 2 * ii0 + 1; ii1 <= 2 * ii0 + (N + 1) / 32 + 1; ii1 += 1) {
                                          for (int i1 = max(64 * ii0 + 37, 32 * ii1); i1 <= min(N + 64 * ii0 + 33, 32 * ii1 + 31); i1 += 1) {
                                            a[-64 * ii0 + i1 - 35] = b[-64 * ii0 + i1 - 35];
                                          }
                                        }
                                        if (T >= 32 * ii0 + 20) {
                                          #pragma omp parallel for
                                          for (int ii1 = 2 * ii0 + 1; ii1 <= 2 * ii0 + (N + 2) / 32 + 1; ii1 += 1) {
                                            if (ii1 >= 2 * ii0 + 2) {
                                              for (int i1 = 32 * ii1; i1 <= min(N + 64 * ii0 + 34, 32 * ii1 + 31); i1 += 1) {
                                                b[-64 * ii0 + i1 - 36] = (0.33333 * ((a[-64 * ii0 + i1 - 37] + a[-64 * ii0 + i1 - 36]) + a[-64 * ii0 + i1 - 35]));
                                              }
                                            } else {
                                              for (int i1 = 64 * ii0 + 38; i1 <= min(N + 64 * ii0 + 34, 64 * ii0 + 63); i1 += 1) {
                                                b[-64 * ii0 + i1 - 36] = (0.33333 * ((a[-64 * ii0 + i1 - 37] + a[-64 * ii0 + i1 - 36]) + a[-64 * ii0 + i1 - 35]));
                                              }
                                            }
                                          }
                                          #pragma omp parallel for
                                          for (int ii1 = 2 * ii0 + 1; ii1 <= 2 * ii0 + (N + 3) / 32 + 1; ii1 += 1) {
                                            for (int i1 = max(64 * ii0 + 39, 32 * ii1); i1 <= min(N + 64 * ii0 + 35, 32 * ii1 + 31); i1 += 1) {
                                              a[-64 * ii0 + i1 - 37] = b[-64 * ii0 + i1 - 37];
                                            }
                                          }
                                          if (T >= 32 * ii0 + 21) {
                                            #pragma omp parallel for
                                            for (int ii1 = 2 * ii0 + 1; ii1 <= 2 * ii0 + (N + 4) / 32 + 1; ii1 += 1) {
                                              if (ii1 >= 2 * ii0 + 2) {
                                                for (int i1 = 32 * ii1; i1 <= min(N + 64 * ii0 + 36, 32 * ii1 + 31); i1 += 1) {
                                                  b[-64 * ii0 + i1 - 38] = (0.33333 * ((a[-64 * ii0 + i1 - 39] + a[-64 * ii0 + i1 - 38]) + a[-64 * ii0 + i1 - 37]));
                                                }
                                              } else {
                                                for (int i1 = 64 * ii0 + 40; i1 <= min(N + 64 * ii0 + 36, 64 * ii0 + 63); i1 += 1) {
                                                  b[-64 * ii0 + i1 - 38] = (0.33333 * ((a[-64 * ii0 + i1 - 39] + a[-64 * ii0 + i1 - 38]) + a[-64 * ii0 + i1 - 37]));
                                                }
                                              }
                                            }
                                            #pragma omp parallel for
                                            for (int ii1 = 2 * ii0 + 1; ii1 <= 2 * ii0 + (N + 5) / 32 + 1; ii1 += 1) {
                                              for (int i1 = max(64 * ii0 + 41, 32 * ii1); i1 <= min(N + 64 * ii0 + 37, 32 * ii1 + 31); i1 += 1) {
                                                a[-64 * ii0 + i1 - 39] = b[-64 * ii0 + i1 - 39];
                                              }
                                            }
                                            if (T >= 32 * ii0 + 22) {
                                              #pragma omp parallel for
                                              for (int ii1 = 2 * ii0 + 1; ii1 <= 2 * ii0 + (N + 6) / 32 + 1; ii1 += 1) {
                                                if (ii1 >= 2 * ii0 + 2) {
                                                  for (int i1 = 32 * ii1; i1 <= min(N + 64 * ii0 + 38, 32 * ii1 + 31); i1 += 1) {
                                                    b[-64 * ii0 + i1 - 40] = (0.33333 * ((a[-64 * ii0 + i1 - 41] + a[-64 * ii0 + i1 - 40]) + a[-64 * ii0 + i1 - 39]));
                                                  }
                                                } else {
                                                  for (int i1 = 64 * ii0 + 42; i1 <= min(N + 64 * ii0 + 38, 64 * ii0 + 63); i1 += 1) {
                                                    b[-64 * ii0 + i1 - 40] = (0.33333 * ((a[-64 * ii0 + i1 - 41] + a[-64 * ii0 + i1 - 40]) + a[-64 * ii0 + i1 - 39]));
                                                  }
                                                }
                                              }
                                              #pragma omp parallel for
                                              for (int ii1 = 2 * ii0 + 1; ii1 <= 2 * ii0 + (N + 7) / 32 + 1; ii1 += 1) {
                                                for (int i1 = max(64 * ii0 + 43, 32 * ii1); i1 <= min(N + 64 * ii0 + 39, 32 * ii1 + 31); i1 += 1) {
                                                  a[-64 * ii0 + i1 - 41] = b[-64 * ii0 + i1 - 41];
                                                }
                                              }
                                              if (T >= 32 * ii0 + 23) {
                                                #pragma omp parallel for
                                                for (int ii1 = 2 * ii0 + 1; ii1 <= 2 * ii0 + (N + 8) / 32 + 1; ii1 += 1) {
                                                  if (ii1 >= 2 * ii0 + 2) {
                                                    for (int i1 = 32 * ii1; i1 <= min(N + 64 * ii0 + 40, 32 * ii1 + 31); i1 += 1) {
                                                      b[-64 * ii0 + i1 - 42] = (0.33333 * ((a[-64 * ii0 + i1 - 43] + a[-64 * ii0 + i1 - 42]) + a[-64 * ii0 + i1 - 41]));
                                                    }
                                                  } else {
                                                    for (int i1 = 64 * ii0 + 44; i1 <= min(N + 64 * ii0 + 40, 64 * ii0 + 63); i1 += 1) {
                                                      b[-64 * ii0 + i1 - 42] = (0.33333 * ((a[-64 * ii0 + i1 - 43] + a[-64 * ii0 + i1 - 42]) + a[-64 * ii0 + i1 - 41]));
                                                    }
                                                  }
                                                }
                                                #pragma omp parallel for
                                                for (int ii1 = 2 * ii0 + 1; ii1 <= 2 * ii0 + (N + 9) / 32 + 1; ii1 += 1) {
                                                  for (int i1 = max(64 * ii0 + 45, 32 * ii1); i1 <= min(N + 64 * ii0 + 41, 32 * ii1 + 31); i1 += 1) {
                                                    a[-64 * ii0 + i1 - 43] = b[-64 * ii0 + i1 - 43];
                                                  }
                                                }
                                                if (T >= 32 * ii0 + 24) {
                                                  #pragma omp parallel for
                                                  for (int ii1 = 2 * ii0 + 1; ii1 <= 2 * ii0 + (N + 10) / 32 + 1; ii1 += 1) {
                                                    if (ii1 >= 2 * ii0 + 2) {
                                                      for (int i1 = 32 * ii1; i1 <= min(N + 64 * ii0 + 42, 32 * ii1 + 31); i1 += 1) {
                                                        b[-64 * ii0 + i1 - 44] = (0.33333 * ((a[-64 * ii0 + i1 - 45] + a[-64 * ii0 + i1 - 44]) + a[-64 * ii0 + i1 - 43]));
                                                      }
                                                    } else {
                                                      for (int i1 = 64 * ii0 + 46; i1 <= min(N + 64 * ii0 + 42, 64 * ii0 + 63); i1 += 1) {
                                                        b[-64 * ii0 + i1 - 44] = (0.33333 * ((a[-64 * ii0 + i1 - 45] + a[-64 * ii0 + i1 - 44]) + a[-64 * ii0 + i1 - 43]));
                                                      }
                                                    }
                                                  }
                                                  #pragma omp parallel for
                                                  for (int ii1 = 2 * ii0 + 1; ii1 <= 2 * ii0 + (N + 11) / 32 + 1; ii1 += 1) {
                                                    for (int i1 = max(64 * ii0 + 47, 32 * ii1); i1 <= min(N + 64 * ii0 + 43, 32 * ii1 + 31); i1 += 1) {
                                                      a[-64 * ii0 + i1 - 45] = b[-64 * ii0 + i1 - 45];
                                                    }
                                                  }
                                                  if (T >= 32 * ii0 + 25) {
                                                    #pragma omp parallel for
                                                    for (int ii1 = 2 * ii0 + 1; ii1 <= 2 * ii0 + (N + 12) / 32 + 1; ii1 += 1) {
                                                      for (int i1 = max(64 * ii0 + 48, 32 * ii1); i1 <= min(N + 64 * ii0 + 44, 32 * ii1 + 31); i1 += 1) {
                                                        b[-64 * ii0 + i1 - 46] = (0.33333 * ((a[-64 * ii0 + i1 - 47] + a[-64 * ii0 + i1 - 46]) + a[-64 * ii0 + i1 - 45]));
                                                      }
                                                    }
                                                    #pragma omp parallel for
                                                    for (int ii1 = 2 * ii0 + 1; ii1 <= 2 * ii0 + (N + 13) / 32 + 1; ii1 += 1) {
                                                      for (int i1 = max(64 * ii0 + 49, 32 * ii1); i1 <= min(N + 64 * ii0 + 45, 32 * ii1 + 31); i1 += 1) {
                                                        a[-64 * ii0 + i1 - 47] = b[-64 * ii0 + i1 - 47];
                                                      }
                                                    }
                                                    if (T >= 32 * ii0 + 26) {
                                                      #pragma omp parallel for
                                                      for (int ii1 = 2 * ii0 + 1; ii1 <= 2 * ii0 + (N + 14) / 32 + 1; ii1 += 1) {
                                                        for (int i1 = max(64 * ii0 + 50, 32 * ii1); i1 <= min(N + 64 * ii0 + 46, 32 * ii1 + 31); i1 += 1) {
                                                          b[-64 * ii0 + i1 - 48] = (0.33333 * ((a[-64 * ii0 + i1 - 49] + a[-64 * ii0 + i1 - 48]) + a[-64 * ii0 + i1 - 47]));
                                                        }
                                                      }
                                                      #pragma omp parallel for
                                                      for (int ii1 = 2 * ii0 + 1; ii1 <= 2 * ii0 + (N + 15) / 32 + 1; ii1 += 1) {
                                                        for (int i1 = max(64 * ii0 + 51, 32 * ii1); i1 <= min(N + 64 * ii0 + 47, 32 * ii1 + 31); i1 += 1) {
                                                          a[-64 * ii0 + i1 - 49] = b[-64 * ii0 + i1 - 49];
                                                        }
                                                      }
                                                      if (T >= 32 * ii0 + 27) {
                                                        #pragma omp parallel for
                                                        for (int ii1 = 2 * ii0 + 1; ii1 <= 2 * ii0 + (N + 16) / 32 + 1; ii1 += 1) {
                                                          for (int i1 = max(64 * ii0 + 52, 32 * ii1); i1 <= min(N + 64 * ii0 + 48, 32 * ii1 + 31); i1 += 1) {
                                                            b[-64 * ii0 + i1 - 50] = (0.33333 * ((a[-64 * ii0 + i1 - 51] + a[-64 * ii0 + i1 - 50]) + a[-64 * ii0 + i1 - 49]));
                                                          }
                                                        }
                                                        #pragma omp parallel for
                                                        for (int ii1 = 2 * ii0 + 1; ii1 <= 2 * ii0 + (N + 17) / 32 + 1; ii1 += 1) {
                                                          for (int i1 = max(64 * ii0 + 53, 32 * ii1); i1 <= min(N + 64 * ii0 + 49, 32 * ii1 + 31); i1 += 1) {
                                                            a[-64 * ii0 + i1 - 51] = b[-64 * ii0 + i1 - 51];
                                                          }
                                                        }
                                                        if (T >= 32 * ii0 + 28) {
                                                          #pragma omp parallel for
                                                          for (int ii1 = 2 * ii0 + 1; ii1 <= 2 * ii0 + (N + 18) / 32 + 1; ii1 += 1) {
                                                            for (int i1 = max(64 * ii0 + 54, 32 * ii1); i1 <= min(N + 64 * ii0 + 50, 32 * ii1 + 31); i1 += 1) {
                                                              b[-64 * ii0 + i1 - 52] = (0.33333 * ((a[-64 * ii0 + i1 - 53] + a[-64 * ii0 + i1 - 52]) + a[-64 * ii0 + i1 - 51]));
                                                            }
                                                          }
                                                          #pragma omp parallel for
                                                          for (int ii1 = 2 * ii0 + 1; ii1 <= 2 * ii0 + (N + 19) / 32 + 1; ii1 += 1) {
                                                            for (int i1 = max(64 * ii0 + 55, 32 * ii1); i1 <= min(N + 64 * ii0 + 51, 32 * ii1 + 31); i1 += 1) {
                                                              a[-64 * ii0 + i1 - 53] = b[-64 * ii0 + i1 - 53];
                                                            }
                                                          }
                                                          if (T >= 32 * ii0 + 29) {
                                                            #pragma omp parallel for
                                                            for (int ii1 = 2 * ii0 + 1; ii1 <= 2 * ii0 + (N + 20) / 32 + 1; ii1 += 1) {
                                                              for (int i1 = max(64 * ii0 + 56, 32 * ii1); i1 <= min(N + 64 * ii0 + 52, 32 * ii1 + 31); i1 += 1) {
                                                                b[-64 * ii0 + i1 - 54] = (0.33333 * ((a[-64 * ii0 + i1 - 55] + a[-64 * ii0 + i1 - 54]) + a[-64 * ii0 + i1 - 53]));
                                                              }
                                                            }
                                                            #pragma omp parallel for
                                                            for (int ii1 = 2 * ii0 + 1; ii1 <= 2 * ii0 + (N + 21) / 32 + 1; ii1 += 1) {
                                                              for (int i1 = max(64 * ii0 + 57, 32 * ii1); i1 <= min(N + 64 * ii0 + 53, 32 * ii1 + 31); i1 += 1) {
                                                                a[-64 * ii0 + i1 - 55] = b[-64 * ii0 + i1 - 55];
                                                              }
                                                            }
                                                            if (T >= 32 * ii0 + 30) {
                                                              #pragma omp parallel for
                                                              for (int ii1 = 2 * ii0 + 1; ii1 <= 2 * ii0 + (N + 22) / 32 + 1; ii1 += 1) {
                                                                for (int i1 = max(64 * ii0 + 58, 32 * ii1); i1 <= min(N + 64 * ii0 + 54, 32 * ii1 + 31); i1 += 1) {
                                                                  b[-64 * ii0 + i1 - 56] = (0.33333 * ((a[-64 * ii0 + i1 - 57] + a[-64 * ii0 + i1 - 56]) + a[-64 * ii0 + i1 - 55]));
                                                                }
                                                              }
                                                              #pragma omp parallel for
                                                              for (int ii1 = 2 * ii0 + 1; ii1 <= 2 * ii0 + (N + 23) / 32 + 1; ii1 += 1) {
                                                                for (int i1 = max(64 * ii0 + 59, 32 * ii1); i1 <= min(N + 64 * ii0 + 55, 32 * ii1 + 31); i1 += 1) {
                                                                  a[-64 * ii0 + i1 - 57] = b[-64 * ii0 + i1 - 57];
                                                                }
                                                              }
                                                              if (T >= 32 * ii0 + 31) {
                                                                #pragma omp parallel for
                                                                for (int ii1 = 2 * ii0 + 1; ii1 <= 2 * ii0 + (N + 24) / 32 + 1; ii1 += 1) {
                                                                  for (int i1 = max(64 * ii0 + 60, 32 * ii1); i1 <= min(N + 64 * ii0 + 56, 32 * ii1 + 31); i1 += 1) {
                                                                    b[-64 * ii0 + i1 - 58] = (0.33333 * ((a[-64 * ii0 + i1 - 59] + a[-64 * ii0 + i1 - 58]) + a[-64 * ii0 + i1 - 57]));
                                                                  }
                                                                }
                                                                #pragma omp parallel for
                                                                for (int ii1 = 2 * ii0 + 1; ii1 <= 2 * ii0 + (N + 25) / 32 + 1; ii1 += 1) {
                                                                  for (int i1 = max(64 * ii0 + 61, 32 * ii1); i1 <= min(N + 64 * ii0 + 57, 32 * ii1 + 31); i1 += 1) {
                                                                    a[-64 * ii0 + i1 - 59] = b[-64 * ii0 + i1 - 59];
                                                                  }
                                                                }
                                                                if (T >= 32 * ii0 + 32) {
                                                                  #pragma omp parallel for
                                                                  for (int ii1 = 2 * ii0 + 1; ii1 <= 2 * ii0 + (N + 26) / 32 + 1; ii1 += 1) {
                                                                    for (int i1 = max(64 * ii0 + 62, 32 * ii1); i1 <= min(N + 64 * ii0 + 58, 32 * ii1 + 31); i1 += 1) {
                                                                      b[-64 * ii0 + i1 - 60] = (0.33333 * ((a[-64 * ii0 + i1 - 61] + a[-64 * ii0 + i1 - 60]) + a[-64 * ii0 + i1 - 59]));
                                                                    }
                                                                  }
                                                                  #pragma omp parallel for
                                                                  for (int ii1 = 2 * ii0 + 1; ii1 <= 2 * ii0 + (N + 27) / 32 + 1; ii1 += 1) {
                                                                    for (int i1 = max(64 * ii0 + 63, 32 * ii1); i1 <= min(N + 64 * ii0 + 59, 32 * ii1 + 31); i1 += 1) {
                                                                      a[-64 * ii0 + i1 - 61] = b[-64 * ii0 + i1 - 61];
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
  }
}
#pragma endscop

    IF_TIME(t_end = rtclock());
    IF_TIME(fprintf(stdout, "%0.6lfs\n", t_end - t_start));

    if (fopen(".test", "r")) {
#ifdef MPI
        if (my_rank == 0) {
            print_array();
        }
#else
        print_array();
#endif
    }

    return 0;
}
