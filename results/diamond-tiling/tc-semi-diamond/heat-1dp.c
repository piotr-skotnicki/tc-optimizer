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
 
 /* TC Optimizing Compiler 0.4.1 */
 /* ./tc ../examples/pluto+/heat-1dp-nt.scop.c --semi-diamond-tiling --omp-for-codegen --iterative-tc --debug --inline -b 128 */
 #define min(x,y)    ((x) < (y) ? (x) : (y))
 #define max(x,y)    ((x) > (y) ? (x) : (y))
 #define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
 #pragma scop
 if (N >= 257) {
   for (register int ii0 = 0; ii0 <= floord(T - 1, 128); ii0 += 1) {
     if (N >= 385) {
       for (register int i0 = 128 * ii0; i0 <= min(T - N + (N + 1) / 2 + 64 * ((N - 1) / 128) - 66, -N + 128 * ii0 + (N + 1) / 2 + 64 * ((N - 1) / 128) + 62); i0 += 1) {
         for (register int i1 = max(-128 * ii0 + i0 + 1, -T + i0 + 129); i1 <= min(-((N - 1) % 128) + T - i0 - 3, -((N - 1) % 128) + 128 * ii0 - i0 + 125); i1 += 1) {
           A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
         }
         for (register int i1 = max(N - 128 * ii0 + i0 - 127, -T + N + i0 + 1); i1 <= min(T - i0 + 254, 128 * ii0 - i0 + 382); i1 += 1) {
           A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
         }
       }
       if (N >= 513) {
         for (register int i0 = 128 * ii0; i0 <= min(T - N + (N + 1) / 2 + 64 * ((N - 1) / 128) - 66, -N + 128 * ii0 + (N + 1) / 2 + 64 * ((N - 1) / 128) + 62); i0 += 1) {
           for (register int i1 = max(N - 128 * ii0 + i0 - 127, -T + N + i0 + 1); i1 < min(-((N - 1) % 128) + T + N - i0 - 130, -((N - 1) % 128) + N + 128 * ii0 - i0 - 2); i1 += 1) {
             A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
           }
         }
       }
     }
     if (2 * T >= ((N - 1) % 128) + 256 * ii0 + 4) {
       for (register int k = 1; k <= min(min(2, ((N - 1) % 128) + 1), -ii0 + (T - ii0 - 2) / 127 + 1); k += 1) {
         {
           if (N >= 385 && k == 2) {
             for (register int i0 = 128 * ii0; i0 <= 128 * ii0 + 126; i0 += 1) {
               for (register int i1 = max(-128 * ii0 + i0 + 1, 128 * ii0 - i0 + 127); i1 <= min(-128 * ii0 + i0 + 128, 128 * ii0 - i0 + 254); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
             for (register int i0 = -N + 128 * ii0 + N / 2 + 64 * ((N - 2) / 128) + 64; i0 <= -((N - 2) % 128) + 128 * ii0 + 125; i0 += 1) {
               for (register int i1 = -((N - 2) % 128) + 128 * ii0 - i0 + 125; i1 <= min(-128 * ii0 + i0, 128 * ii0 - i0 + 126); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 && (N + i0 + 1) % 128 == 0 ? N - 1 : i1 - 1]));
               }
             }
             for (register int i0 = -((N - 2) % 128) + 128 * ii0 + 126; i0 <= 128 * ii0 + 126; i0 += 1) {
               for (register int i1 = 0; i1 <= min(-128 * ii0 + i0, 128 * ii0 - i0 + 126); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
               }
               for (register int i1 = max(N - 128 * ii0 + i0 - 127, -((N - 2) % 128) + N + 128 * ii0 - i0 + 125); i1 < N; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           }
           if (k == 2) {
             #pragma omp parallel for
             for (register int ii1 = 1; ii1 < floord(N - 1, 128) - 1; ii1 += 1) {
               for (register int i0 = 128 * ii0; i0 <= 128 * ii0 + 126; i0 += 1) {
                 for (register int i1 = max(-128 * ii0 + 128 * ii1 + i0 + 1, 128 * ii0 + 128 * ii1 - i0 + 127); i1 <= min(-128 * ii0 + 128 * ii1 + i0 + 128, 128 * ii0 + 128 * ii1 - i0 + 254); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             }
           }
           if (N <= 384 && 2 * T + 253 >= N + 256 * ii0) {
             if (N >= 258 && k == 2) {
               if (T >= N + 128 * ii0) {
                 for (register int i0 = 128 * ii0; i0 <= 128 * ii0 + 126; i0 += 1) {
                   for (register int i1 = max(-128 * ii0 + i0 + 1, 128 * ii0 - i0 + 127); i1 <= min(-128 * ii0 + i0 + 128, 128 * ii0 - i0 + 254); i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                 }
               }
               for (register int ii1_2 = 1; ii1_2 <= ii0 + floord(-T + N - 1, 128) + 1; ii1_2 += 1) {
                 if (ii1_2 == 1) {
                   for (register int i0 = 128 * ii0; i0 <= min(128 * ii0 + 126, T - N + N / 2 - 1); i0 += 1) {
                     for (register int i1 = max(-128 * ii0 + i0 + 1, 128 * ii0 - i0 + 127); i1 <= min(-128 * ii0 + i0 + 128, 128 * ii0 - i0 + 254); i1 += 1) {
                       A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                     }
                   }
                   for (register int i0 = max(128 * ii0, T - N + N / 2); i0 < 64 * ii0 + T / 2 - 64; i0 += 1) {
                     for (register int i1 = 128 * ii0 - i0 + 127; i1 <= -128 * ii0 + i0 + 128; i1 += 1) {
                       A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                     }
                   }
                 } else {
                   for (register int i0 = -N + 128 * ii0 + N / 2 + 192; i0 <= -N + 64 * ii0 + (T + N) / 2 + 127; i0 += 1) {
                     for (register int i1 = max(0, -N + 128 * ii0 - i0 + 383); i1 <= -128 * ii0 + i0; i1 += 1) {
                       A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
                     }
                     for (register int i1 = 128 * ii0 - i0 + 383; i1 < N; i1 += 1) {
                       A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                     }
                   }
                 }
                 if (128 * ii0 + 127 * ii1_2 + 128 >= T) {
                   for (register int i0 = max(max(T - N + N / 2, -N + 64 * ii0 + 63 * ii1_2 + (T + N + ii1_2) / 2 + 1), 64 * ii0 - 64 * ii1_2 + (T + ii1_2 + 1) / 2 - 1); i0 <= min(64 * ii0 + 63 * ii1_2 + (T + ii1_2) / 2 - 65, 64 * ii0 - 64 * ii1_2 + (T + ii1_2 + 1) / 2 + 125); i0 += 1) {
                     for (register int i1 = max(-128 * ii0 + i0 + 1, 128 * ii0 - i0 + 127); i1 <= T - N - i0 + 254; i1 += 1) {
                       A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                     }
                     if (ii1_2 == 2) {
                       for (register int i1 = max(0, -N + 128 * ii0 - i0 + 383); i1 <= T - N - i0 + 254; i1 += 1) {
                         A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
                       }
                     }
                     for (register int i1 = max(max(max(0, -128 * ii0 - 127 * ii1_2 + i0 + 128), 128 * ii0 - 127 * ii1_2 - i0 + 254), T - N - i0 + 255); i1 <= min(min(-128 * ii0 - 127 * ii1_2 + i0 + 254, -T + i0 + 256), 128 * ii0 - 127 * ii1_2 - i0 + 380); i1 += 1) {
                       A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][ii1_2 == 2 && i1 == 0 ? N - 1 : i1 - 1]));
                     }
                     if (ii1_2 == 1) {
                       for (register int i1 = -T + i0 + 257; i1 <= min(-128 * ii0 + i0 + 128, 128 * ii0 - i0 + 254); i1 += 1) {
                         A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                       }
                     } else {
                       for (register int i1 = max(N - 128 * ii0 + i0 - 127, 128 * ii0 - i0 + 383); i1 < N; i1 += 1) {
                         A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                       }
                     }
                   }
                 }
                 if (ii1_2 == 1) {
                   for (register int i0 = 64 * ii0 + (T + 1) / 2 - 1; i0 <= 128 * ii0 + 126; i0 += 1) {
                     for (register int i1 = -128 * ii0 + i0 + 1; i1 <= 128 * ii0 - i0 + 254; i1 += 1) {
                       A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                     }
                   }
                   if (T >= 128 * ii0 + 256) {
                     for (register int i0 = T - N + N / 2; i0 <= 128 * ii0 + 126; i0 += 1) {
                       for (register int i1 = -128 * ii0 + i0 + 1; i1 <= 128 * ii0 - i0 + 254; i1 += 1) {
                         A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                       }
                     }
                   }
                 } else {
                   for (register int i0 = 64 * ii0 + (T + 1) / 2 - 1; i0 <= T - N + 254; i0 += 1) {
                     for (register int i1 = max(0, -N + 128 * ii0 - i0 + 383); i1 <= 128 * ii0 - i0 + 126; i1 += 1) {
                       A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
                     }
                     for (register int i1 = max(N - 128 * ii0 + i0 - 127, 128 * ii0 - i0 + 383); i1 < N; i1 += 1) {
                       A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                     }
                   }
                   for (register int i0 = max(T - N + 255, 64 * ii0 + (T + 1) / 2 - 1); i0 <= 128 * ii0 + 126; i0 += 1) {
                     for (register int i1 = 0; i1 <= 128 * ii0 - i0 + 126; i1 += 1) {
                       A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
                     }
                     for (register int i1 = max(N - 128 * ii0 + i0 - 127, 128 * ii0 - i0 + 383); i1 < N; i1 += 1) {
                       A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                     }
                   }
                 }
               }
               if (T + 128 >= N + 128 * ii0) {
                 for (register int i0 = -N + 128 * ii0 + N / 2 + 192; i0 <= 128 * ii0 + 126; i0 += 1) {
                   for (register int i1 = max(0, -N + 128 * ii0 - i0 + 383); i1 <= min(-128 * ii0 + i0, 128 * ii0 - i0 + 126); i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
                   }
                   for (register int i1 = max(N - 128 * ii0 + i0 - 127, 128 * ii0 - i0 + 383); i1 < N; i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                 }
               }
             } else if (N + 128 * ii0 >= T + 255) {
               for (register int i0 = 128 * ii0; i0 <= min(T - N + (N + 1) / 2 + 126, -N + 128 * ii0 + (N + 1) / 2 + 254); i0 += 1) {
                 for (register int i1 = max(-128 * ii0 + i0 + 1, -T + i0 + 129); i1 <= min(T - i0 - 2, 128 * ii0 - i0 + 126); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
                 for (register int i1 = max(-128 * ii0 + i0 + 129, -T + i0 + 257); i1 <= min(T - i0 + 126, 128 * ii0 - i0 + 254); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
                 for (register int i1 = max(N - 128 * ii0 + i0 - 127, -T + N + i0 + 1); i1 <= min(T - i0 + 254, 128 * ii0 - i0 + 382); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][N + 128 * ii0 == T + 255 && N + i0 == T + 255 && i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             } else if (N >= 383) {
               for (register int i0 = 128 * ii0; i0 <= -N + 128 * ii0 + 446; i0 += 1) {
                 for (register int i1 = -128 * ii0 + i0 + 1; i1 <= 128 * ii0 - i0 + 126; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
                 for (register int i1 = -128 * ii0 + i0 + 129; i1 <= 128 * ii0 - i0 + 254; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
                 for (register int i1 = N - 128 * ii0 + i0 - 127; i1 <= 128 * ii0 - i0 + 382; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][N == 383 && i0 == 128 * ii0 && i1 == 382 ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             } else {
               for (register int i0 = 128 * ii0; i0 <= min(T - N + (N + 1) / 2 + 126, -N + 128 * ii0 + (N + 1) / 2 + 254); i0 += 1) {
                 for (register int i1 = 0; i1 <= min(T - N - i0 + 254, -N + 128 * ii0 - i0 + 382); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
                 }
                 for (register int i1 = max(-128 * ii0 + i0 + 1, -N + 128 * ii0 - i0 + 383); i1 <= min(128 * ii0 - i0 + 126, T - N - i0 + 254); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
                 for (register int i1 = max(max(-128 * ii0 + i0 + 1, -T + i0 + 129), T - N - i0 + 255); i1 <= min(T - i0 - 2, 128 * ii0 - i0 + 126); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
                 for (register int i1 = max(-128 * ii0 + i0 + 129, -T + i0 + 257); i1 <= min(T - i0 + 126, 128 * ii0 - i0 + 254); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
                 for (register int i1 = max(N - 128 * ii0 + i0 - 127, 128 * ii0 - i0 + 255); i1 <= min(min(N - 1, T - i0 + 126), 128 * ii0 - i0 + 382); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
                 for (register int i1 = max(max(N - 128 * ii0 + i0 - 127, -T + N + i0 + 1), T - i0 + 127); i1 <= min(min(N - 1, T - i0 + 254), 128 * ii0 - i0 + 382); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             }
           }
           if (N >= 258 && k == 2) {
             if (N >= 385) {
               for (register int i0 = -N + 128 * ii0 + N / 2 + 64 * ((N - 2) / 128) + 64; i0 <= 128 * ii0 + 126; i0 += 1) {
                 for (register int i1 = max(-((N - 2) % 128) + N - 128 * ii0 + i0 - 129, -((N - 2) % 128) + N + 128 * ii0 - i0 - 3); i1 <= min(N - 128 * ii0 + i0 - 128, -((N - 2) % 128) + N + 128 * ii0 - i0 + 124); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             } else {
               for (register int i0 = -N + 128 * ii0 + N / 2 + 192; i0 <= min(128 * ii0 + 126, -N + 64 * ii0 + (T + N + 1) / 2 + 190); i0 += 1) {
                 for (register int i1 = max(-128 * ii0 + i0 + 129, 128 * ii0 - i0 + 255); i1 <= min(N - 128 * ii0 + i0 - 128, 128 * ii0 - i0 + 382); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
               for (register int i0 = -N + 64 * ii0 + (T + N + 1) / 2 + 191; i0 <= 128 * ii0 + 126; i0 += 1) {
                 for (register int i1 = -128 * ii0 + i0 + 129; i1 <= 128 * ii0 - i0 + 382; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             }
           }
         }
         if (N >= 385 && T >= 128 * ii0 + 66 && k == 1) {
           if (T + 128 >= N + 128 * ii0) {
             for (register int i0 = 128 * ii0; i0 <= 128 * ii0 + 62; i0 += 1) {
               for (register int i1 = -128 * ii0 + i0 + 129; i1 <= 128 * ii0 - i0 + 254; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
             for (register int i0 = 128 * ii0; i0 <= 128 * ii0 + 62; i0 += 1) {
               for (register int i1 = max(-128 * ii0 + i0 + 1, -((N - 1) % 128) + 128 * ii0 - i0 + 126); i1 <= 128 * ii0 - i0 + 126; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
             for (register int i0 = 128 * ii0; i0 <= -N + 128 * ii0 + (N + 1) / 2 + 64 * ((N - 1) / 128) + 126; i0 += 1) {
               for (register int i1 = 0; i1 <= min(-128 * ii0 + i0, -((N - 1) % 128) + 128 * ii0 - i0 + 125); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
               }
               for (register int i1 = max(N - 128 * ii0 + i0 - 127, -((N - 1) % 128) + N + 128 * ii0 - i0 - 2); i1 <= min(N - 1, -((N - 1) % 128) + N + 128 * ii0 - i0 + 125); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           }
           #pragma omp parallel for
           for (register int ii1 = 0; ii1 < min(ii0 + floord(-T + N - 1, 128), (N - 1) / 128 - 1); ii1 += 1) {
             {
               for (register int i0 = 128 * ii0; i0 <= min(T - 66, 128 * ii0 + 62); i0 += 1) {
                 for (register int i1 = max(-128 * ii0 + 128 * ii1 + i0 + 129, -T + 128 * ii1 + i0 + 257); i1 <= min(min(-T + N + i0, T + 128 * ii1 - i0 + 126), 128 * ii0 + 128 * ii1 - i0 + 254); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
                 for (register int i1 = -T + N + i0 + 1; i1 <= min(N - 128 * ii0 + i0 - 128, 128 * ii0 + 128 * ii1 - i0 + 254); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
               if (ii1 == 0) {
                 for (register int i0 = 128 * ii0; i0 <= min(T - 66, 128 * ii0 + 62); i0 += 1) {
                   if (T >= 128 * ii0 + 129) {
                     for (register int i1 = max(-128 * ii0 + i0 + 1, -((N - 1) % 128) + 128 * ii0 - i0 + 126); i1 <= 128 * ii0 - i0 + 126; i1 += 1) {
                       A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                     }
                   } else {
                     for (register int i1 = max(-T + i0 + 129, -((N - 1) % 128) + T - i0 - 2); i1 < T - i0 - 1; i1 += 1) {
                       A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                     }
                   }
                 }
               }
             }
             if (ii1 == 0) {
               for (register int i0 = 128 * ii0; i0 <= min(T - N + (N + 1) / 2 + 64 * ((N - 1) / 128) - 2, -N + 128 * ii0 + (N + 1) / 2 + 64 * ((N - 1) / 128) + 126); i0 += 1) {
                 for (register int i1 = 0; i1 <= min(min(-T + i0 + 128, -((N - 1) % 128) + T - i0 - 3), -((N - 1) % 128) + 128 * ii0 - i0 + 125); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
                 }
                 for (register int i1 = max(0, -T + i0 + 129); i1 <= min(-128 * ii0 + i0, -((N - 1) % 128) + 128 * ii0 - i0 + 125); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
                 }
                 if (128 * ii0 + 128 >= T) {
                   for (register int i1 = max(-T + N + i0 + 1, -((N - 1) % 128) + T + N - i0 - 130); i1 < min(N, -((N - 1) % 128) + T + N - i0 - 2); i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                 } else {
                   for (register int i1 = max(N - 128 * ii0 + i0 - 127, -((N - 1) % 128) + N + 128 * ii0 - i0 - 2); i1 <= min(N - 1, -((N - 1) % 128) + N + 128 * ii0 - i0 + 125); i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                 }
               }
             }
           }
           #pragma omp parallel for
           for (register int ii1 = max(1, ii0 + floord(-T + N - 1, 128)); ii1 < floord(N - 1, 128) - 1; ii1 += 1) {
             for (register int i0 = 128 * ii0; i0 <= 128 * ii0 + 62; i0 += 1) {
               for (register int i1 = -128 * ii0 + 128 * ii1 + i0 + 129; i1 <= min(N - 128 * ii0 + i0 - 128, 128 * ii0 + 128 * ii1 - i0 + 254); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           }
         } else if (N >= 385 && T >= 128 * ii0 + 3 && 128 * ii0 + 65 >= T && k == 1) {
           for (register int i0 = 128 * ii0; i0 < T - N + 128 * ii0 + (N + 1) / 2 - 64 * floord(2 * T - N - 3, 128) - 1; i0 += 1) {
             for (register int i1 = 0; i1 < T - N + 256 * ii0 - i0 - 128 * floord(2 * T - N - 3, 128) - 1; i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
             }
             for (register int i1 = -T + N + i0 + 1; i1 < min(N, T + 256 * ii0 - i0 - 128 * floord(2 * T - N - 3, 128) - 1); i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
           }
         } else if (N >= 385 && 128 * ii0 + 2 == T && k == 1 && (N - 1) % 128 == 0) {
           A[1][N - 1] = (0.125 * ((A[0][0] - (2.0 * A[0][N - 1])) + A[0][N - 2]));
         }
       }
       if (N >= 385 && T >= 128 * ii0 + 129 && (N - 1) % 128 == 0) {
         #pragma omp parallel for
         for (register int ii1 = 0; ii1 < (N - 129) / 128; ii1 += 1) {
           for (register int i0 = 128 * ii0; i0 <= 128 * ii0 + 126; i0 += 1) {
             for (register int i1 = max(max(-128 * ii0 + 128 * ii1 + i0 + 1, 128 * ii0 + 128 * ii1 - i0 + 127), 128 * ii0 - i0 + 255); i1 <= min(-128 * ii0 + 128 * ii1 + i0 + 128, 128 * ii0 + 128 * ii1 - i0 + 254); i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
             if (ii1 == 0) {
               for (register int i1 = max(-128 * ii0 + i0 + 1, 128 * ii0 - i0 + 127); i1 <= min(-128 * ii0 + i0 + 128, 128 * ii0 - i0 + 254); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           }
           if (ii1 == 0) {
             for (register int i0 = 128 * ii0 + 63; i0 <= 128 * ii0 + 126; i0 += 1) {
               A[(i0 + 1) % 2][128 * ii0 - i0 + 126] = (0.125 * ((A[i0 % 2][128 * ii0 - i0 + 127] - (2.0 * A[i0 % 2][128 * ii0 - i0 + 126])) + A[i0 % 2][i0 == 128 * ii0 + 126 ? N - 1 : 128 * ii0 - i0 + 125]));
             }
           }
         }
         for (register int i0 = 128 * ii0 + 63; i0 <= 128 * ii0 + 126; i0 += 1) {
           A[(i0 + 1) % 2][N - 128 * ii0 + i0 - 128] = (0.125 * ((A[i0 % 2][N - 128 * ii0 + i0 - 127] - (2.0 * A[i0 % 2][N - 128 * ii0 + i0 - 128])) + A[i0 % 2][N - 128 * ii0 + i0 - 129]));
         }
       } else if (N >= 385 && 128 * ii0 + 128 >= T) {
         #pragma omp parallel for
         for (register int ii1 = 0; ii1 < (N - 1) / 128 - 1; ii1 += 1) {
           for (register int i0 = 128 * ii0; i0 < T - 1; i0 += 1) {
             for (register int i1 = max(max(T + 128 * ii1 - i0 - 1, T - i0 + 127), -T + 128 * ii1 + i0 + 129); i1 <= min(T + 128 * ii1 - i0 + 126, -T + 128 * ii1 + i0 + 256); i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
             if (ii1 == 0) {
               for (register int i1 = max(T - i0 - 1, -T + i0 + 129); i1 <= min(T - i0 + 126, -T + i0 + 256); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           }
           if (ii1 == 0) {
             for (register int i0 = max(128 * ii0, T - N + N / 2 + 64 * ((N - 1) / 128) - 64); i0 < -((N - 1) % 128) + T - 1; i0 += 1) {
               for (register int i1 = -((N - 1) % 128) + T - i0 - 2; i1 <= min(T - i0 - 2, -T + i0 + 128); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 && (-T + N + i0 + 1) % 128 == 0 ? N - 1 : i1 - 1]));
               }
             }
             for (register int i0 = max(128 * ii0, -((N - 1) % 128) + T - 1); i0 < T - 1; i0 += 1) {
               for (register int i1 = 0; i1 <= min(T - i0 - 2, -T + i0 + 128); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
               }
               for (register int i1 = max(-T + N + i0 + 1, -((N - 1) % 128) + T + N - i0 - 2); i1 < N; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           }
         }
         for (register int i0 = max(128 * ii0, T - N + N / 2 + 64 * ((N - 1) / 128) - 64); i0 < T - 1; i0 += 1) {
           for (register int i1 = max(-((N - 1) % 128) + T + N - i0 - 130, -((N - 1) % 128) - T + N + i0); i1 <= min(-T + N + i0, -((N - 1) % 128) + T + N - i0 - 3); i1 += 1) {
             A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
           }
         }
       }
     } else if (N >= 385) {
       for (register int i0 = 128 * ii0; i0 < T - 1; i0 += 1) {
         for (register int i1 = -T + i0 + 129; i1 <= T - i0 + 126; i1 += 1) {
           A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
         }
       }
       for (register int i0 = 128 * ii0; i0 < T - 1; i0 += 1) {
         for (register int i1 = 0; i1 < T - i0 - 1; i1 += 1) {
           A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
         }
         for (register int i1 = -T + N + i0 + 1; i1 < N; i1 += 1) {
           A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
         }
       }
       #pragma omp parallel for
       for (register int ii1 = 1; ii1 < floord(N - 1, 128) - 1; ii1 += 1) {
         for (register int i0 = 128 * ii0; i0 < T - 1; i0 += 1) {
           for (register int i1 = -T + 128 * ii1 + i0 + 129; i1 <= T + 128 * ii1 - i0 + 126; i1 += 1) {
             A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
           }
         }
       }
       for (register int ii1_2 = N <= 512 ? 3 : N - 127 * N / 128 - 1; ii1_2 <= (N <= 512 ? 3 : N - 127 * N / 128 - 1); ii1_2 += 1) {
         for (register int i0 = 128 * ii0; i0 < T - 1; i0 += 1) {
           for (register int i1 = -T + 128 * ii1_2 + i0 + 1; i1 < T + 128 * ii1_2 - i0 - 1; i1 += 1) {
             A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
           }
         }
       }
     }
     if (N >= 385) {
       #pragma omp parallel for
       for (register int ii1 = 0; ii1 <= (N - 1) / 128; ii1 += 1) {
         if (128 * ii0 + 128 >= T) {
           for (register int i0 = max(max(max(T - 64, 128 * ii0), T - 64 * ii1 + 64), T - N + 64 * ii1 + N / 2); i0 < T; i0 += 1) {
             if (N >= 128 * ii1 + 257) {
               for (register int i1 = T + 128 * ii1 - i0 - 1; i1 <= -T + 128 * ii1 + i0 + 128; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             } else {
               for (register int i1 = T + 128 * ii1 - i0 - 1; i1 <= min(-T + N + i0, -T + 128 * ii1 + i0 + 128); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i0 + 1 == T && i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           }
           if (ii1 <= 1) {
             for (register int i0 = max(T - 64, 128 * ii0); i0 < T; i0 += 1) {
               if (ii1 == 1) {
                 for (register int i1 = T - i0 + 127; i1 <= -T + i0 + 256; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               } else {
                 for (register int i1 = T - i0 - 1; i1 <= -T + i0 + 128; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i0 + 1 == T && i1 == 0 ? N - 1 : i1 - 1]));
                 }
               }
             }
           }
         } else if (N >= 128 * ii1 + 257) {
           for (register int i0 = 128 * ii0 + 64; i0 <= 128 * ii0 + 127; i0 += 1) {
             if (ii1 >= 1) {
               for (register int i1 = 128 * ii0 + 128 * ii1 - i0 + 127; i1 <= -128 * ii0 + 128 * ii1 + i0; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             } else {
               for (register int i1 = 128 * ii0 - i0 + 127; i1 <= -128 * ii0 + i0; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i0 == 128 * ii0 + 127 && i1 == 0 ? N - 1 : i1 - 1]));
               }
             }
           }
         } else {
           for (register int i0 = max(128 * ii0 + 64, -N + 128 * ii0 + 64 * ii1 + N / 2 + 128); i0 <= 128 * ii0 + 127; i0 += 1) {
             for (register int i1 = 128 * ii0 + 128 * ii1 - i0 + 127; i1 <= min(N - 128 * ii0 + i0 - 128, -128 * ii0 + 128 * ii1 + i0); i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i0 == 128 * ii0 + 127 && i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
           }
         }
       }
     } else {
       if (N >= 258 && 128 * ii0 + 128 >= T && 2 * T + 253 >= N + 256 * ii0) {
         #pragma omp parallel for
         for (register int ii1 = 0; ii1 <= 1; ii1 += 1) {
           for (register int ii1_2 = ii1 + 1; ii1_2 <= 2; ii1_2 += 1) {
             if (ii1_2 == 2) {
               for (register int i0 = max(128 * ii0, T - N + N / 2 + 64); i0 < T - 1; i0 += 1) {
                 if (ii1 == 1) {
                   for (register int i1 = max(T - i0 + 127, -T + i0 + 257); i1 <= min(-T + N + i0, T - i0 + 254); i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                 } else {
                   for (register int i1 = max(0, T - N - i0 + 255); i1 <= min(T - i0 - 2, -T + i0 + 128); i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
                   }
                   for (register int i1 = max(-T + N + i0 + 1, T - i0 + 255); i1 < N; i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                 }
               }
             } else {
               for (register int i0 = 128 * ii0; i0 < T - 1; i0 += 1) {
                 for (register int i1 = max(T - i0 - 1, -T + i0 + 129); i1 <= min(T - i0 + 126, -T + i0 + 256); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             }
           }
         }
       }
       if (N >= 259 && 2 * T + 253 >= N + 256 * ii0) {
         for (register int i0 = 128 * ii0 + 64; i0 <= min(T - 65, 128 * ii0 + 127); i0 += 1) {
           for (register int i1 = 128 * ii0 - i0 + 127; i1 <= -128 * ii0 + i0; i1 += 1) {
             A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i0 == 128 * ii0 + 127 && i1 == 0 ? N - 1 : i1 - 1]));
           }
         }
         for (register int i0 = max(T - 64, 128 * ii0); i0 <= min(T - 1, 128 * ii0 + 127); i0 += 1) {
           for (register int i1 = 128 * ii0 - i0 + 127; i1 < T - i0 - 1; i1 += 1) {
             A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i0 == 128 * ii0 + 127 && i1 == 0 ? N - 1 : i1 - 1]));
           }
           for (register int i1 = T - i0 - 1; i1 <= -T + i0 + 128; i1 += 1) {
             A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i0 + 1 == T && i1 == 0 ? N - 1 : i1 - 1]));
           }
           for (register int i1 = -T + i0 + 129; i1 <= -128 * ii0 + i0; i1 += 1) {
             A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
           }
         }
         if (N + 128 * ii0 >= T + 1) {
           if (T >= 128 * ii0 + 128) {
             for (register int i0 = 128 * ii0 + 64; i0 <= 128 * ii0 + 127; i0 += 1) {
               for (register int i1 = 128 * ii0 - i0 + 255; i1 <= -128 * ii0 + i0 + 128; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           } else {
             for (register int i0 = max(T - 64, 128 * ii0); i0 < T; i0 += 1) {
               for (register int i1 = T - i0 + 127; i1 <= -T + i0 + 256; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           }
           if (128 * ii0 + 128 >= T) {
             for (register int i0 = T - N + N / 2 + 128; i0 < T; i0 += 1) {
               for (register int i1 = T - i0 + 255; i1 <= -T + N + i0; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i0 + 1 == T && i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           }
         } else {
           for (register int i0 = 128 * ii0 + 64; i0 <= 128 * ii0 + 127; i0 += 1) {
             for (register int i1 = 128 * ii0 - i0 + 255; i1 <= -128 * ii0 + i0 + 128; i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
           }
         }
         if (T >= 128 * ii0 + 129) {
           for (register int i0 = -N + 128 * ii0 + N / 2 + 256; i0 <= 128 * ii0 + 127; i0 += 1) {
             for (register int i1 = 128 * ii0 - i0 + 383; i1 < N - 128 * ii0 + i0 - 127; i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i0 == 128 * ii0 + 127 && i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
           }
         }
       } else if (N + 256 * ii0 >= 2 * T + 254) {
         for (register int k = 2 * ii0 - (T + 62) / 64 + 3; k <= min(3, N - 256); k += 1) {
           #pragma omp parallel for
           for (register int ii1 = 0; ii1 < k; ii1 += 1) {
             if (T >= 128 * ii0 + 64 * ii1 + 1 && 64 * k >= 127 * ii1 + 65) {
               for (register int ii1_2 = -k + 3; ii1_2 <= -2 * k + 6; ii1_2 += 1) {
                 if (k + ii1_2 == 3) {
                   for (register int i0 = max(128 * ii0, T + 64 * k - 256); i0 < min(T, T + 64 * k - 129); i0 += 1) {
                     if (ii1 == 0) {
                       for (register int i1 = max(T - i0 - 1, -T - 128 * k + i0 + 385); i1 <= min(T - i0 + 126, -T - 128 * k + i0 + 512); i1 += 1) {
                         A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][k == 3 && i0 + 1 == T && i1 == 0 ? N - 1 : i1 - 1]));
                       }
                     } else {
                       for (register int i1 = T - i0 + 127; i1 <= -T + i0 + 256; i1 += 1) {
                         A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                       }
                     }
                   }
                 } else {
                   for (register int i0 = 128 * ii0; i0 < T - 1; i0 += 1) {
                     for (register int i1 = 0; i1 < T - i0 - 1; i1 += 1) {
                       A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
                     }
                     for (register int i1 = -T + N + i0 + 1; i1 < N; i1 += 1) {
                       A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                     }
                   }
                 }
               }
             } else if (k == 3 && 128 * ii0 + ii1 + 63 >= T) {
               for (register int i0 = max(128 * ii0, T - N + (N + ii1) / 2 + 127); i0 < T; i0 += 1) {
                 if (ii1 == 2) {
                   for (register int i1 = T - i0 + 255; i1 <= -T + N + i0; i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i0 + 1 == T && i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                 } else {
                   for (register int i1 = T - i0 + 127; i1 <= -T + i0 + 256; i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                 }
               }
             } else {
               for (register int i0 = 128 * ii0; i0 < T - 1; i0 += 1) {
                 for (register int i1 = -T + i0 + 257; i1 <= T - i0 + 254; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             }
           }
         }
       }
     }
     for (register int k = max(N - 255, 2 * ii0 - (T + 62) / 64 + 3); k <= 3; k += 1) {
       #pragma omp parallel for
       for (register int ii1 = 0; ii1 <= -k + 3; ii1 += 1) {
         if (ii1 == 0) {
           if (T >= 128 * ii0 + 129) {
             for (register int i0 = 128 * ii0 + 64 * k - 128; i0 <= min(min(128 * ii0 + 127, 128 * ii0 + 64 * k - 2), 64 * ii0 + 64 * k + (T + 1) / 2 - 130); i0 += 1) {
               for (register int i1 = max(128 * ii0 - i0 + 127, -128 * ii0 - 128 * k + i0 + 257); i1 <= min(min(128 * ii0 - i0 + 254, T - N - i0 + 254), -128 * ii0 - 128 * k + i0 + 384); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][k == 3 && i0 == 128 * ii0 + 127 && i1 == 0 ? N - 1 : i1 - 1]));
               }
               for (register int i1 = max(128 * ii0 - i0 + 127, T - N - i0 + 255); i1 <= min(T - i0 - 2, -T + i0 + 128); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i0 == 128 * ii0 + 127 && i1 == 0 ? N - 1 : i1 - 1]));
               }
               if (N == 257 && k == 2 && i0 + 129 >= T) {
                 A[(i0 + 1) % 2][T - i0 - 2] = (0.125 * ((A[i0 % 2][T - i0 - 1] - (2.0 * A[i0 % 2][T - i0 - 2])) + A[i0 % 2][T - i0 - 3]));
               }
               for (register int i1 = T - i0 - 1; i1 <= -T - 128 * k + i0 + 512; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
               for (register int i1 = max(T - N - i0 + 255, -T - 128 * k + i0 + 513); i1 <= min(128 * ii0 - i0 + 254, -128 * ii0 - 128 * k + i0 + 384); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
             if (N == 257 && k == 2) {
               for (register int i0 = 64 * ii0 + (T + 1) / 2 - 1; i0 <= 128 * ii0 + 126; i0 += 1) {
                 for (register int i1 = -128 * ii0 + i0 + 1; i1 <= 128 * ii0 - i0 + 254; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             }
           } else {
             for (register int i0 = max(128 * ii0, T + 64 * k - 256); i0 < min(T, T + 64 * k - 129); i0 += 1) {
               for (register int i1 = max(T - i0 - 1, -T - 128 * k + i0 + 385); i1 <= min(T - i0 + 126, -T - 128 * k + i0 + 512); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][k == 3 && i0 + 1 == T && i1 == 0 ? N - 1 : i1 - 1]));
               }
             }
           }
         }
         if (N == 257 && k == 2) {
           if (T >= 128 * ii0 + 129) {
             for (register int i0 = 128 * ii0 + 63; i0 <= 128 * ii0 + 126; i0 += 1) {
               if (ii1 == 1) {
                 A[(i0 + 1) % 2][-128 * ii0 + i0 + 129] = (0.125 * ((A[i0 % 2][-128 * ii0 + i0 + 130] - (2.0 * A[i0 % 2][-128 * ii0 + i0 + 129])) + A[i0 % 2][-128 * ii0 + i0 + 128]));
               } else {
                 A[(i0 + 1) % 2][128 * ii0 - i0 + 126] = (0.125 * ((A[i0 % 2][128 * ii0 - i0 + 127] - (2.0 * A[i0 % 2][128 * ii0 - i0 + 126])) + A[i0 % 2][i0 == 128 * ii0 + 126 ? 256 : 128 * ii0 - i0 + 125]));
               }
             }
           } else {
             for (register int i0 = max(T - 65, 128 * ii0); i0 < T - 1; i0 += 1) {
               if (ii1 == 1) {
                 A[(i0 + 1) % 2][-T + i0 + 257] = (0.125 * ((A[i0 % 2][-T + i0 + 258] - (2.0 * A[i0 % 2][-T + i0 + 257])) + A[i0 % 2][-T + i0 + 256]));
               } else {
                 A[(i0 + 1) % 2][T - i0 - 2] = (0.125 * ((A[i0 % 2][T - i0 - 1] - (2.0 * A[i0 % 2][T - i0 - 2])) + A[i0 % 2][i0 + 2 == T ? 256 : T - i0 - 3]));
               }
             }
           }
         }
       }
       if (N + 128 * ii0 >= T + 128 && T >= 128 * ii0 + 128 && k == 3) {
         for (register int i0 = 128 * ii0 + 64; i0 <= 128 * ii0 + 127; i0 += 1) {
           for (register int i1 = 128 * ii0 - i0 + 255; i1 <= -128 * ii0 + i0 + 128; i1 += 1) {
             A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
           }
         }
       }
       if (k == 3) {
         #pragma omp parallel for
         for (register int ii1 = max(1, T - 128 * ii0 - 126); ii1 <= 2; ii1 += 1) {
           for (register int i0 = max(128 * ii0, T + 63 * ii1 - 127); i0 < T; i0 += 1) {
             if (ii1 == 1) {
               for (register int i1 = T - i0 + 127; i1 <= -T + i0 + 256; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             } else {
               for (register int i1 = 256; i1 < N; i1 += 1) {
                 A[T % 2][i1] = (0.125 * ((A[-(T % 2) + 1][i1 + 1 == N ? 0 : 257] - (2.0 * A[-(T % 2) + 1][i1])) + A[-(T % 2) + 1][i1 - 1]));
               }
             }
           }
         }
         if (T >= 128 * ii0 + 129) {
           #pragma omp parallel for
           for (register int ii1 = max(1, -T + 64 * ii0 + (T + N) / 2 - 62); ii1 <= 2; ii1 += 1) {
             for (register int i0 = 128 * ii0 + 63 * ii1 + 1; i0 <= 128 * ii0 + 127; i0 += 1) {
               if (ii1 == 1) {
                 for (register int i1 = 128 * ii0 - i0 + 255; i1 <= -128 * ii0 + i0 + 128; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               } else {
                 for (register int i1 = 256; i1 < N; i1 += 1) {
                   A[0][i1] = (0.125 * ((A[1][i1 + 1 == N ? 0 : 257] - (2.0 * A[1][i1])) + A[1][i1 - 1]));
                 }
               }
             }
           }
         }
       }
     }
   }
 } else if (T >= 1 && N >= 1) {
   for (register int ii0 = 0; ii0 <= min(floord(T - 1, 128), floord(T - 64 * ((N + 127) / 128) - 2, 128) + 1); ii0 += 1) {
     if (N <= 128) {
       for (register int i0 = 128 * ii0; i0 <= min(T - 1, 128 * ii0 + 127); i0 += 1) {
         for (register int i1 = 0; i1 < N; i1 += 1) {
           A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
         }
       }
     } else {
       for (register int i0 = 128 * ii0; i0 <= min(T - 2, 128 * ii0 + 126); i0 += 1) {
         for (register int i1 = 0; i1 <= min(T - i0 - 2, 128 * ii0 - i0 + 126); i1 += 1) {
           A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
         }
         for (register int i1 = max(-128 * ii0 + i0 + 1, 128 * ii0 - i0 + 127); i1 <= min(min(N - 128 * ii0 + i0 - 128, T - i0 - 2), 128 * ii0 - i0 + 254); i1 += 1) {
           A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
         }
         for (register int i1 = max(max(T - i0 - 1, -128 * ii0 + i0 + 1), -T + i0 + 129); i1 <= min(min(-T + N + i0, T - i0 + 126), 128 * ii0 - i0 + 254); i1 += 1) {
           A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
         }
         for (register int i1 = max(max(T - i0 - 1, -T + N + i0 + 1), -128 * ii0 + i0 + 1); i1 <= min(N - 128 * ii0 + i0 - 128, 128 * ii0 - i0 + 254); i1 += 1) {
           A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
         }
         for (register int i1 = max(N - 128 * ii0 + i0 - 127, 128 * ii0 - i0 + 127); i1 < min(N, T - i0 - 1); i1 += 1) {
           A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
         }
         for (register int i1 = max(max(N - 128 * ii0 + i0 - 127, T - i0 - 1), -T + N + i0 + 1); i1 < N; i1 += 1) {
           A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
         }
       }
       if (N == 256 && 128 * ii0 + 128 == T) {
         #pragma omp parallel for
         for (register int ii1 = 0; ii1 <= 1; ii1 += 1) {
           for (register int i0 = T - 64; i0 < T; i0 += 1) {
             if (ii1 == 1) {
               for (register int i1 = T - i0 + 127; i1 <= -T + i0 + 256; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i0 + 1 == T && i1 == 255 ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             } else {
               for (register int i1 = T - i0 - 1; i1 <= -T + i0 + 128; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i0 + 1 == T && i1 == 0 ? 255 : i1 - 1]));
               }
             }
           }
         }
       } else if (T >= 128 * ii0 + 192 && 128 * ii0 + 383 >= T + N) {
         #pragma omp parallel for
         for (register int ii1 = 0; ii1 <= 1; ii1 += 1) {
           for (register int i0 = max(128 * ii0 + 64, -N + 128 * ii0 + 63 * ii1 + N / 2 + 129); i0 <= 128 * ii0 + 127; i0 += 1) {
             if (ii1 == 1) {
               for (register int i1 = 128 * ii0 - i0 + 255; i1 < N - 128 * ii0 + i0 - 127; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i0 == 128 * ii0 + 127 && i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             } else {
               for (register int i1 = 128 * ii0 - i0 + 127; i1 <= -128 * ii0 + i0; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i0 == 128 * ii0 + 127 && i1 == 0 ? N - 1 : i1 - 1]));
               }
             }
           }
         }
       } else if (T + N >= 128 * ii0 + 384 && T >= 128 * ii0 + 129) {
         #pragma omp parallel for
         for (register int ii1 = 0; ii1 <= 1; ii1 += 1) {
           if (ii1 == 1) {
             for (register int i0 = -N + 128 * ii0 + N / 2 + 192; i0 <= 128 * ii0 + 127; i0 += 1) {
               for (register int i1 = 128 * ii0 - i0 + 255; i1 < N - 128 * ii0 + i0 - 127; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i0 == 128 * ii0 + 127 && i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           } else if (N <= 255) {
             for (register int i0 = 128 * ii0 + 64; i0 <= 128 * ii0 + 127; i0 += 1) {
               for (register int i1 = 128 * ii0 - i0 + 127; i1 <= -128 * ii0 + i0; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i0 == 128 * ii0 + 127 && i1 == 0 ? N - 1 : i1 - 1]));
               }
             }
           } else {
             for (register int i0 = 128 * ii0 + 64; i0 <= 128 * ii0 + 127; i0 += 1) {
               for (register int i1 = 128 * ii0 - i0 + 127; i1 <= -128 * ii0 + i0; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i0 == 128 * ii0 + 127 && i1 == 0 ? 255 : i1 - 1]));
               }
             }
           }
         }
       } else if (128 * ii0 + 128 >= T) {
         for (register int i0 = max(T - 64, 128 * ii0); i0 <= T - N + N / 2 + 63; i0 += 1) {
           for (register int i1 = T - i0 - 1; i1 <= -T + i0 + 128; i1 += 1) {
             A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
           }
         }
         for (register int i0 = max(128 * ii0, T - N + N / 2 + 64); i0 < T; i0 += 1) {
           for (register int i1 = T - i0 - 1; i1 <= -T + i0 + 128; i1 += 1) {
             A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i0 + 1 == T && i1 == 0 ? N - 1 : i1 - 1]));
           }
         }
         if (128 * ii0 + 128 == T) {
           for (register int i0 = T - N + N / 2 + 64; i0 < T; i0 += 1) {
             for (register int i1 = T - i0 + 127; i1 <= -T + N + i0; i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i0 + 1 == T && i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
           }
         } else {
           for (register int i0 = max(128 * ii0, T - N + N / 2 + 64); i0 < T; i0 += 1) {
             for (register int i1 = T - i0 + 127; i1 <= -T + N + i0; i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i0 + 1 == T && i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
           }
         }
       } else {
         #pragma omp parallel for
         for (register int ii1 = 0; ii1 <= 1; ii1 += 1) {
           if (ii1 == 1) {
             for (register int i0 = -N + 128 * ii0 + N / 2 + 192; i0 <= min(128 * ii0 + 127, -N + 64 * ii0 + (T + N) / 2 + 127); i0 += 1) {
               for (register int i1 = 128 * ii0 - i0 + 255; i1 < N - 128 * ii0 + i0 - 127; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i0 == 128 * ii0 + 127 && i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           }
           for (register int i0 = max(128 * ii0 + 64, -N + 64 * ii0 + 94 * ii1 + (T + N + ii1 + 1) / 2 + 33); i0 <= 128 * ii0 + 127; i0 += 1) {
             if (ii1 == 1) {
               for (register int i1 = 128 * ii0 - i0 + 255; i1 < N - 128 * ii0 + i0 - 127; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i0 == 128 * ii0 + 127 && i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             } else {
               for (register int i1 = 128 * ii0 - i0 + 127; i1 <= -128 * ii0 + i0; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i0 == 128 * ii0 + 127 && i1 == 0 ? N - 1 : i1 - 1]));
               }
             }
           }
         }
       }
     }
   }
   if (N >= 129 && (T - 1) % 128 == 0) {
     #pragma omp parallel for
     for (register int ii1 = 0; ii1 <= 1; ii1 += 1) {
       if (ii1 == 1) {
         for (register int i1 = 128; i1 < N; i1 += 1) {
           A[1][i1] = (0.125 * ((A[0][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[0][i1])) + A[0][i1 - 1]));
         }
       } else {
         for (register int i1 = 0; i1 <= 127; i1 += 1) {
           A[1][i1] = (0.125 * ((A[0][i1 + 1] - (2.0 * A[0][i1])) + A[0][i1 == 0 ? N - 1 : i1 - 1]));
         }
       }
     }
   }
 }
 #pragma endscop
 
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
 