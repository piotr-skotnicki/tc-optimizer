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
 
 /* TC Optimizing Compiler 0.4.2 */
 /* ./tc ../examples/pluto+/heat-1dp-nt.scop.c --diamond-tiling --omp-for-codegen --iterative-tc --debug --inline -b 128 */
 #define min(x,y)    ((x) < (y) ? (x) : (y))
 #define max(x,y)    ((x) > (y) ? (x) : (y))
 #define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
 #pragma scop
 if (N <= 128) {
   if (T >= 129) {
     for (register int i0 = 0; i0 <= min(255, T - 1); i0 += 1) {
       for (register int i1 = 0; i1 < N; i1 += 1) {
         A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
       }
     }
   } else {
     for (register int i0 = 0; i0 < T; i0 += 1) {
       for (register int i1 = 0; i1 < N; i1 += 1) {
         A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
       }
     }
   }
   for (register int ii0 = 1; ii0 < floord(T - 1, 128); ii0 += 1) {
     for (register int i0 = 128 * ii0 + 128; i0 <= min(T - 1, 128 * ii0 + 255); i0 += 1) {
       for (register int i1 = 0; i1 < N; i1 += 1) {
         A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
       }
     }
   }
 } else if (N <= 384) {
   if (T + N <= 256 && 193 * N + 16381 >= 510 * T) {
     for (register int k = max(2, -T + 4); k <= 3; k += 1) {
       #pragma omp parallel for
       for (register int ii1 = 0; ii1 < k - 1; ii1 += 1) {
         if (k == 3 && ii1 == 1) {
           for (register int i0 = max(0, T - N + N / 2 + 64); i0 < T; i0 += 1) {
             for (register int i1 = T - i0 + 127; i1 <= -T + N + i0; i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i0 + 1 == T && i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
           }
         } else if (k == 3) {
           for (register int i0 = max(0, T - 64); i0 < T; i0 += 1) {
             for (register int i1 = T - i0 - 1; i1 <= -T + i0 + 128; i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i0 + 1 == T && i1 == 0 ? N - 1 : i1 - 1]));
             }
           }
         } else {
           for (register int i0 = 0; i0 < T - 1; i0 += 1) {
             for (register int i1 = 0; i1 < T - i0 - 1; i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
             }
             for (register int i1 = max(T - i0 - 1, -T + i0 + 129); i1 <= min(-T + N + i0, T - i0 + 126); i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
             for (register int i1 = max(T - i0 - 1, -T + N + i0 + 1); i1 < N; i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
           }
         }
       }
     }
   } else if (T >= N && N <= 256) {
     for (register int k = 2; k <= 3; k += 1) {
       #pragma omp parallel for
       for (register int ii1 = 0; ii1 < k - 1; ii1 += 1) {
         if (k == 3 && ii1 == 1) {
           for (register int i0 = -N + N / 2 + 192; i0 <= 127; i0 += 1) {
             for (register int i1 = -i0 + 255; i1 < N + i0 - 127; i1 += 1) {
               A[(-i0 + 127) % 2][i1] = (0.125 * ((A[(-i0 + 128) % 2][i0 == 127 && i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[(-i0 + 128) % 2][i1])) + A[(-i0 + 128) % 2][i1 - 1]));
             }
           }
           for (register int i0 = 128; i0 <= (N + 1) / 2 + 62; i0 += 1) {
             for (register int i1 = i0 + 1; i1 <= N - i0 + 126; i1 += 1) {
               A[-(i0 % 2) + 1][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
           }
         } else if (k == 3) {
           if (N == 256) {
             for (register int i0 = 64; i0 <= 127; i0 += 1) {
               for (register int i1 = -i0 + 127; i1 <= i0; i1 += 1) {
                 A[(-i0 + 127) % 2][i1] = (0.125 * ((A[(-i0 + 128) % 2][i1 + 1] - (2.0 * A[(-i0 + 128) % 2][i1])) + A[(-i0 + 128) % 2][i0 == 127 && i1 == 0 ? 255 : i1 - 1]));
               }
             }
           } else {
             for (register int i0 = 64; i0 <= 127; i0 += 1) {
               for (register int i1 = -i0 + 127; i1 <= i0; i1 += 1) {
                 A[(-i0 + 127) % 2][i1] = (0.125 * ((A[(-i0 + 128) % 2][i1 + 1] - (2.0 * A[(-i0 + 128) % 2][i1])) + A[(-i0 + 128) % 2][i0 == 127 && i1 == 0 ? N - 1 : i1 - 1]));
               }
             }
           }
           for (register int i0 = 128; i0 <= min(190, T - 1); i0 += 1) {
             for (register int i1 = i0 - 127; i1 <= -i0 + 254; i1 += 1) {
               A[-(i0 % 2) + 1][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
           }
         } else {
           for (register int i0 = 0; i0 <= 126; i0 += 1) {
             for (register int i1 = 0; i1 <= -i0 + 126; i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
             }
             for (register int i1 = max(i0 + 1, -i0 + 127); i1 <= min(N + i0 - 128, -i0 + 254); i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
             for (register int i1 = max(N + i0 - 127, -i0 + 127); i1 < N; i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
           }
         }
       }
     }
   }
   if (N <= 256) {
     for (register int ii0 = 1; ii0 <= min(floord(T - N, 128), floord(T + N, 128) - 3); ii0 += 1) {
       for (register int k = 2; k <= 3; k += 1) {
         #pragma omp parallel for
         for (register int ii1 = 0; ii1 < k - 1; ii1 += 1) {
           if (k == 3) {
             if (ii1 == 1) {
               for (register int i0 = -N + 128 * ii0 + N / 2 + 192; i0 <= 128 * ii0 + 127; i0 += 1) {
                 for (register int i1 = 128 * ii0 - i0 + 255; i1 < N - 128 * ii0 + i0 - 127; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i0 == 128 * ii0 + 127 && i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             } else {
               for (register int i0 = 128 * ii0 + 64; i0 <= -N + 128 * ii0 + N / 2 + 191; i0 += 1) {
                 for (register int i1 = 128 * ii0 - i0 + 127; i1 <= -128 * ii0 + i0; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             }
             for (register int i0 = max(128 * ii0 + 64 * ii1 + 64, -N + 128 * ii0 - 63 * ii1 + N / 2 + 192); i0 <= min(min(128 * ii0 + 190, 128 * ii0 - 64 * ii1 + (N + 1) / 2 + 126), -63 * N + 128 * ii0 + 63 * ii1 + 16255); i0 += 1) {
               for (register int i1 = max(-128 * ii0 + 128 * ii1 + i0 - 127, 128 * ii0 - i0 + 127); i1 <= min(min(-128 * ii0 + 126 * ii1 + i0, N + 128 * ii0 - i0 + 126), 128 * ii0 + 253 * ii1 - i0 + 254); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][ii1 == 0 && i0 == 128 * ii0 + 127 && i1 == 0 ? N - 1 : i1 - 1]));
               }
             }
             if (N == 256 && ii1 == 0) {
               for (register int i0 = 128 * ii0 + 128; i0 <= 128 * ii0 + 190; i0 += 1) {
                 for (register int i1 = -128 * ii0 + i0 - 127; i1 <= 128 * ii0 - i0 + 254; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             }
           } else {
             for (register int i0 = 128 * ii0; i0 <= 128 * ii0 + 126; i0 += 1) {
               for (register int i1 = 0; i1 <= min(-128 * ii0 + i0, 128 * ii0 - i0 + 126); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
               }
               if (N >= 253 && N <= 254) {
                 for (register int i1 = max(-128 * ii0 + i0 + 1, 128 * ii0 - i0 + 127); i1 <= min(N + 128 * ii0 - i0 - 2, -384 * ii0 + i0 + 128 * ((-2 * ii0 + 2 * i0 + 128) / 127)); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               } else if (2 * ((43 * N - 43) / 170) + 66 >= N) {
                 for (register int i1 = max(-128 * ii0 + i0 + 1, 128 * ii0 - i0 + 127); i1 <= min(N + 128 * ii0 - i0 - 2, -384 * ii0 + i0 + 128 * ((-2 * ii0 + 2 * i0 + 128) / 127)); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               } else if (N <= 252 && i0 >= 128 * ii0 + 63) {
                 for (register int i1 = -128 * ii0 + i0 + 1; i1 < N + 128 * ii0 - i0 - 1; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               } else if (N <= 252 && 128 * ii0 + 62 >= i0) {
                 for (register int i1 = 128 * ii0 - i0 + 127; i1 <= min(N + 128 * ii0 - i0 - 2, -128 * ii0 + i0 + 128); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               } else {
                 for (register int i1 = max(-128 * ii0 + i0 + 1, 128 * ii0 - i0 + 127); i1 <= min(N + 128 * ii0 - i0 - 2, -128 * ii0 + i0 + 128); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
               for (register int i1 = max(N + 128 * ii0 - i0 - 1, -128 * ii0 + i0 + 1); i1 <= min(N - 128 * ii0 + i0 - 128, 128 * ii0 - i0 + 254); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
               for (register int i1 = max(N - 128 * ii0 + i0 - 127, N + 128 * ii0 - i0 - 1); i1 < N; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           }
         }
       }
     }
     if (T >= 320 && ((T + N + 128) % 128) + 64 >= N) {
       for (register int k = 2; k <= 3; k += 1) {
         #pragma omp parallel for
         for (register int ii1 = 0; ii1 < k - 1; ii1 += 1) {
           if (k == 3 && ii1 == 0) {
             for (register int i0 = -((T + N + 128) % 128) + T + N - 192; i0 < T / 2 + 64 * ((T + N + 128) / 128) - 192; i0 += 1) {
               for (register int i1 = -((T + N + 128) % 128) + T + N - i0 - 129; i1 <= ((T + N + 128) % 128) - T - N + i0 + 256; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           }
           for (register int i0 = max(max(-((T + N + 128) % 128) + T + N - 256, 126 * k + 64 * ((T + N + 128) / 128) + (T + k + 1) / 2 - 572), -((T + N + 128) % 128) + T + 95 * k + 31 * ii1 + (N + k + ii1) / 2 - 382); i0 < min(-((T + N + 128) % 128) + T + N + 64 * k - 257, -((T + N + 128) % 128) + T + N + 63 * k - 63 * ii1 - 254); i0 += 1) {
             {
               if ((ii1 == 0 && i0 <= 255) || (T + N >= 512 && ii1 == 0 && T + N >= ((T + N) % 128) + i0 + 129 && ((T + N) % 128) + 64 >= N)) {
                 for (register int i1 = max(0, -((T + N + 128) % 128) + T + N + 127 * k - i0 - 510); i1 <= min(min(T - i0 - 2, -((T + N + 128) % 128) + T + N + 127 * k - i0 - 384), ((T + N + 128) % 128) - T - N + i0 + 256); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
                 }
                 if (k == 2 && ((T + N + 128) % 128) + 64 >= N) {
                   for (register int i1 = max(-((T + N + 128) % 128) + T + N - i0 - 129, ((T + N + 128) % 128) - T - N + i0 + 257); i1 <= min(-((T + N + 128) % 128) + T + 2 * N - i0 - 258, 3 * ((T + N + 128) % 128) - 3 * T - 3 * N + i0 + 128 * ((2 * i0 - 2 * ((T + N + 128) / 128) + 134) / 127) + 768); i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                 }
               }
               if (ii1 + 2 == k && ((T + N + 128) % 128) + 64 >= N) {
                 for (register int i1 = max(max(-((T + N + 128) % 128) + T + N + 127 * k - i0 - 382, -((T + N + 128) % 128) + T + 2 * N - 128 * k - i0 - 1), ((T + N + 128) % 128) - T - N + i0 + 257); i1 <= min(-((T + N + 128) % 128) + T + N + 128 * k - i0 - 258, ((T + N + 128) % 128) - T + i0 + 128); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][k == 3 && i1 + 1 == N && (i0 + 1) % 128 == 0 ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
                 if (k == 2) {
                   for (register int i1 = max(-((T + N + 128) % 128) + T + 2 * N - i0 - 257, ((T + N + 128) % 128) - T + i0 + 129); i1 < N; i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                 }
               } else if ((k == 3 && ii1 == 0 && i0 <= 255) || (T + N >= 512 && k == 3 && ii1 == 0 && T + N >= ((T + N) % 128) + i0 + 129 && ((T + N) % 128) + 64 >= N)) {
                 for (register int i1 = T - i0 - 1; i1 <= ((T + N + 128) % 128) - T - N + i0 + 256; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             }
             if (k == 3 && ii1 == 0 && ((T + N + 128) % 128) + 64 >= N && ((T + N + 128) % 128) + i0 + 128 >= T + N) {
               for (register int i1 = ((T + N + 128) % 128) - T - N + i0 + 129; i1 < -((T + N + 128) % 128) + T + N - i0 - 1; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           }
           if (k == 3 && ii1 == 1) {
             for (register int i0 = -((T + N) % 128) + T + N - 128; i0 < -((T + N) % 128) + T + N + (N + 1) / 2 - 193; i0 += 1) {
               for (register int i1 = ((T + N + 128) % 128) - T - N + i0 + 257; i1 < -((T + N + 128) % 128) + T + 2 * N - i0 - 129; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           }
         }
       }
     } else if (T + N >= 512 && T + N + 128 * floord(-T + N - 1, 128) >= 384) {
       if (N == 256) {
         if (T % 128 >= 2) {
           for (register int i0 = -(T % 128) + T - 128; i0 < -(T % 128) + T - 1; i0 += 1) {
             for (register int i1 = 0; i1 <= min(-(T % 128) + T - i0 - 2, (T % 128) - T + i0 + 128); i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? 255 : i1 - 1]));
             }
             for (register int i1 = max(-(T % 128) + T - i0 - 1, (T % 128) - T + i0 + 129); i1 <= min(-(T % 128) + T - i0 + 126, (T % 128) - T + i0 + 256); i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
             for (register int i1 = max(-(T % 128) + T - i0 + 127, (T % 128) - T + i0 + 257); i1 <= 255; i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 == 255 ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
           }
         } else {
           for (register int i0 = -(T % 128) + T - 128; i0 < -(T % 128) + T - 1; i0 += 1) {
             for (register int i1 = 0; i1 <= min(-(T % 128) + T - i0 - 2, (T % 128) - T + i0 + 128); i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? 255 : i1 - 1]));
             }
             if ((T - 1) % 128 == 0) {
               for (register int i1 = max(T - i0 - 2, -T + i0 + 130); i1 <= min(T - i0 + 125, -3 * T + i0 + 128 * ((-T + 128 * i0 + 8384) / 8128) + 387); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             } else {
               for (register int i1 = max(T - i0 - 1, -T + i0 + 129); i1 <= min(T - i0 + 126, -3 * T + i0 + 128 * ((-T + 128 * i0 + 8383) / 8128) + 384); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
             for (register int i1 = max(-(T % 128) + T - i0 + 127, (T % 128) - T + i0 + 257); i1 <= 255; i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 == 255 ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
           }
         }
       } else if (N <= 254) {
         for (register int i0 = -((T - N + 128) % 128) + T - N + 128; i0 <= -((T - N + 128) % 128) + T - N + 254; i0 += 1) {
           for (register int i1 = 0; i1 <= min(((T - N + 128) % 128) - T + N + i0 - 128, -((T - N + 128) % 128) + T - N - i0 + 254); i1 += 1) {
             A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
           }
           if (N <= 252 && ((T - N + 128) % 128) + 2 * N >= 512 && ((T - N + 128) % 128) + N + i0 >= T + 191) {
             for (register int i1 = ((T - N + 128) % 128) - T + N + i0 - 127; i1 < -((T - N + 128) % 128) + T - i0 + 127; i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
           } else if (N >= 253 && ((T - N + 128) % 128) + 2 * N >= 512) {
             for (register int i1 = max(((T - N + 128) % 128) - T + N + i0 - 127, -((T - N + 128) % 128) + T - N - i0 + 255); i1 <= min(-((T - N + 128) % 128) + T - i0 + 126, 3 * ((T - N + 128) % 128) - 3 * T + 3 * N + i0 + 128 * ((2 * i0 - 2 * ((T - N + 128) / 128) + 128) / 127) - 384); i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
           } else {
             for (register int i1 = -((T - N + 128) % 128) + T - N - i0 + 255; i1 <= min(-((T - N + 128) % 128) + T - i0 + 126, ((T - N + 128) % 128) - T + N + i0); i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
           }
           for (register int i1 = max(-((T - N + 128) % 128) + T - i0 + 127, ((T - N + 128) % 128) - T + N + i0 - 127); i1 <= min(((T - N + 128) % 128) - T + 2 * N + i0 - 256, -((T - N + 128) % 128) + T - N - i0 + 382); i1 += 1) {
             A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
           }
           for (register int i1 = max(((T - N + 128) % 128) - T + 2 * N + i0 - 255, -((T - N + 128) % 128) + T - i0 + 127); i1 < N; i1 += 1) {
             A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
           }
         }
       } else {
         for (register int i0 = -((T - 127) % 128) + T - 127; i0 < -((T - 127) % 128) + T; i0 += 1) {
           for (register int i1 = 0; i1 <= min(-((T - 127) % 128) + T - i0 - 1, ((T - 127) % 128) - T + i0 + 127); i1 += 1) {
             A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? 254 : i1 - 1]));
           }
           if ((T - 127) % 128 >= 3) {
             for (register int i1 = max(-((T - 127) % 128) + T - i0, ((T - 127) % 128) - T + i0 + 128); i1 <= min(T - i0 - 1, ((T - 127) % 128) - T + i0 + 255); i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
             for (register int i1 = max(T - i0, ((T - 127) % 128) - T + i0 + 128); i1 <= 127; i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
           } else {
             if (T >= i0 + 66 && (T - 1) % 128 == 0) {
               for (register int i1 = T - i0 - 2; i1 < T - i0; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
             for (register int i1 = max(T - i0, -T + i0 + 130); i1 <= 127; i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
           }
           for (register int i1 = max(128, T - i0); i1 <= min(min(190, T - i0 + 124), -((T - 127) % 128) + T - i0 + 126); i1 += 1) {
             if (((T - 127) % 128 >= 3 && ((T - 127) % 128) + i0 + 255 >= T + i1) || (i0 - i1) % 128 == 0 || (i0 + 257 >= T + i1 && (T - 1) % 128 == 0)) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
           }
           if ((T - 127) % 128 >= 2 && ((T - 127) % 128) + i0 + 63 >= T) {
             A[(i0 + 1) % 2][-((T - 127) % 128) + T - i0 + 127] = (0.125 * ((A[i0 % 2][-((T - 127) % 128) + T - i0 + 128] - (2.0 * A[i0 % 2][-((T - 127) % 128) + T - i0 + 127])) + A[i0 % 2][-((T - 127) % 128) + T - i0 + 126]));
           }
           for (register int i1 = max(-((T - 127) % 128) + T - i0 + 127, ((T - 127) % 128) - T + i0 + 255); i1 <= 254; i1 += 1) {
             A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 == 254 ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
           }
         }
       }
       #pragma omp parallel for
       for (register int ii1 = 0; ii1 <= 1; ii1 += 1) {
         if (ii1 == 1 && ((T - N + 128) % 128) + 2 * N >= 512) {
           if ((T - N + 128) % 128 >= 1) {
             for (register int i0 = -((T - N) % 128) + T - 2 * N + N / 2 + 320; i0 <= -((T - N) % 128) + T - N + 255; i0 += 1) {
               for (register int i1 = -((T - N) % 128) + T - N - i0 + 383; i1 < ((T - N) % 128) - T + 2 * N + i0 - 255; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N && (i0 + 1) % 128 == 0 ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
             for (register int i0 = -((T - N + 128) % 128) + T - N + 256; i0 <= min(T - 1, -((T - N + 128) % 128) + T - N + (N + 1) / 2 + 190); i0 += 1) {
               for (register int i1 = ((T - N + 128) % 128) - T + N + i0 - 127; i1 <= -((T - N + 128) % 128) + T - i0 + 254; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           } else {
             for (register int i0 = T - 64; i0 < T; i0 += 1) {
               for (register int i1 = T - i0 + 127; i1 <= -T + i0 + 256; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i0 + 1 == T && i1 == 255 ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           }
         } else if (2 * ((T - N + 128) % 128) + N >= 384) {
           for (register int i0 = -((T - N + 128) % 128) + T - N + 192; i0 <= min(T - 1, -((T - N) % 128) + T - N + 318); i0 += 1) {
             if (N <= 255 && ((T - N + 128) % 128) + 2 * N >= 512 && ((T - N + 128) % 128) + N + i0 >= T + 256 && 2 * ((T - N + 128) % 128) + N >= 384) {
               for (register int i1 = ((T - N + 128) % 128) - T + N + i0 - 255; i1 <= -((T - N + 128) % 128) + T - N - i0 + 382; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             } else if (N == 256 && T % 128 >= 64 && (T % 128) + i0 >= T) {
               for (register int i1 = (T % 128) - T + i0 + 1; i1 < -(T % 128) + T - i0 + 127; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             } else {
               for (register int i1 = -((T - N + 128) % 128) + T - N - i0 + 255; i1 <= ((T - N + 128) % 128) - T + N + i0 - 128; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 && (i0 + 1) % 128 == 0 ? N - 1 : i1 - 1]));
               }
             }
           }
         } else {
           for (register int i0 = -((T - N + 128) % 128) + T - N + 192; i0 <= T - N + N / 2 + 63; i0 += 1) {
             for (register int i1 = -((T - N + 128) % 128) + T - N - i0 + 255; i1 <= ((T - N + 128) % 128) - T + N + i0 - 128; i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
           }
           for (register int i0 = T - N + N / 2 + 64; i0 <= -((T - N + 128) % 128) + T - N + 255; i0 += 1) {
             for (register int i1 = -((T - N + 128) % 128) + T - N - i0 + 255; i1 <= ((T - N + 128) % 128) - T + N + i0 - 128; i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 && (i0 + 1) % 128 == 0 ? N - 1 : i1 - 1]));
             }
           }
           for (register int i0 = -((T - N + 128) % 128) + T - N + 256; i0 < T; i0 += 1) {
             for (register int i1 = ((T - N + 128) % 128) - T + N + i0 - 255; i1 <= -((T - N + 128) % 128) + T - N - i0 + 382; i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
           }
         }
       }
     } else if (T >= 129 && N >= T + 1) {
       for (register int k = 2; k <= 3; k += 1) {
         #pragma omp parallel for
         for (register int ii1 = 0; ii1 < k - 1; ii1 += 1) {
           if (k == 3 && ii1 == 1) {
             for (register int i0 = -N + N / 2 + 192; i0 <= min(T - 1, (N + 1) / 2 + 62); i0 += 1) {
               for (register int i1 = max(i0 + 1, -i0 + 255); i1 <= min(N + i0 - 128, N - i0 + 126); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i0 == 127 && i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           } else if (k == 3) {
             for (register int i0 = 64; i0 <= -N + N / 2 + 191; i0 += 1) {
               for (register int i1 = -i0 + 127; i1 <= i0; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
             for (register int i0 = -N + N / 2 + 192; i0 <= 127; i0 += 1) {
               for (register int i1 = -i0 + 127; i1 <= i0; i1 += 1) {
                 A[(-i0 + 127) % 2][i1] = (0.125 * ((A[(-i0 + 128) % 2][i1 + 1] - (2.0 * A[(-i0 + 128) % 2][i1])) + A[(-i0 + 128) % 2][i0 == 127 && i1 == 0 ? N - 1 : i1 - 1]));
               }
             }
             for (register int i0 = 128; i0 <= min(190, T - 1); i0 += 1) {
               if (2 * T + 125 * N <= 32383) {
                 for (register int i1 = i0 - 127; i1 <= -i0 + 254; i1 += 1) {
                   A[-(i0 % 2) + 1][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               } else {
                 for (register int i1 = i0 - 127; i1 < T - i0 - 1; i1 += 1) {
                   A[-(i0 % 2) + 1][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
                 for (register int i1 = max(i0 - 127, T - i0 - 1); i1 <= -i0 + 254; i1 += 1) {
                   A[-(i0 % 2) + 1][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             }
           } else {
             for (register int i0 = 0; i0 <= 126; i0 += 1) {
               for (register int i1 = 0; i1 <= -i0 + 126; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
               }
               for (register int i1 = max(i0 + 1, -i0 + 127); i1 <= min(N + i0 - 128, -i0 + 254); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
               for (register int i1 = max(N + i0 - 127, -i0 + 127); i1 < N; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           }
         }
       }
     }
     for (register int ii0 = max(max(1, floord(T + 64, 128) - 1), floord(T + N, 128) - 2); ii0 <= floord(T - 1, 128); ii0 += 1) {
       {
         if (N >= 255 && T >= 128 * ii0 + 65 && 128 * ii0 + 127 >= T) {
           for (register int i0 = 128 * ii0; i0 < T - 1; i0 += 1) {
             for (register int i1 = 0; i1 <= min(T - i0 - 2, -128 * ii0 + i0); i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
             }
             if (T + 62 * N >= 128 * ii0 + 15937) {
               for (register int i1 = max(128 * ii0 - i0 + 127, -T + i0 + 129); i1 <= min(min(T - i0 + 126, -384 * ii0 + i0 + 128 * ((-2 * ii0 + 2 * i0 + 128) / 127)), T - N - 389 * ii0 + 2 * i0 + (-5 * ii0 + 5 * i0 + 5) / 127 + 257); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             } else {
               for (register int i1 = max(128 * ii0 - i0 + 127, -T + i0 + 129); i1 <= min(min(T - i0 + 126, -128 * ii0 + i0 + 128), T - 389 * ii0 + 2 * i0 + (-5 * ii0 + 5 * i0 + 5) / 127 + 2); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
             for (register int i1 = max(128 * ii0 - i0 + 127, T - N - 389 * ii0 + 2 * i0 + (-5 * ii0 + 5 * i0 + 5) / 127 + 258); i1 <= min(T - i0 + 126, -128 * ii0 + i0 + 128); i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
             for (register int i1 = max(N + 128 * ii0 - i0 - 1, -T + N + i0 + 1); i1 < N; i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
           }
         } else if (N >= 255 && N + 128 * ii0 >= T + 193) {
           for (register int i0 = 128 * ii0; i0 < T - 1; i0 += 1) {
             for (register int i1 = 0; i1 <= min(T - i0 - 2, -128 * ii0 + i0); i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
             }
             for (register int i1 = max(128 * ii0 - i0 + 127, -T + i0 + 129); i1 <= min(T - i0 + 126, T - N - 389 * ii0 + 2 * i0 + (-5 * ii0 + 5 * i0 + 5) / 127 + 257); i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
             if (N == 256 && 128 * ii0 + 63 == T) {
               for (register int i1 = max(T - i0 + 64, -2 * T + 2 * i0 - (5 * T - 5 * i0 + 187) / 127 + 194); i1 <= min(T - i0 + 126, -T + i0 + 191); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             } else {
               for (register int i1 = max(max(128 * ii0 - i0 + 127, -T + i0 + 129), T - N - 389 * ii0 + 2 * i0 + (-5 * ii0 + 5 * i0 + 5) / 127 + 258); i1 <= min(T - i0 + 126, -128 * ii0 + i0 + 128); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
             for (register int i1 = max(N + 128 * ii0 - i0 - 1, -T + N + i0 + 1); i1 < N; i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
           }
         } else if (T + N >= 128 * ii0 + 256 && N + 128 * ii0 >= T + 130 && 128 * ii0 + 64 >= T && 128 * ii0 + 318 >= T + N && T + 15492 >= 61 * N + 128 * ii0) {
           for (register int i0 = 128 * ii0; i0 < T - 1; i0 += 1) {
             for (register int i1 = 0; i1 <= min(T - i0 - 2, -128 * ii0 + i0); i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
             }
             if (N == 255 && 128 * ii0 + 63 == T) {
               for (register int i1 = max(T - i0 + 64, -T + i0 + 129); i1 <= min(T - i0 + 126, -2 * T + 2 * i0 - (5 * T - 5 * i0 + 187) / 127 + 194); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             } else {
               for (register int i1 = max(128 * ii0 - i0 + 127, -T + i0 + 129); i1 <= min(min(T - i0 + 126, -128 * ii0 + i0 + 128), T - N - 389 * ii0 + 2 * i0 + (-5 * ii0 + 5 * i0 + 5) / 127 + 257); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
             for (register int i1 = max(max(128 * ii0 - i0 + 127, -T + i0 + 129), T - N - 389 * ii0 + 2 * i0 + (-5 * ii0 + 5 * i0 + 5) / 127 + 258); i1 <= min(T - i0 + 126, -128 * ii0 + i0 + 128); i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
             for (register int i1 = max(N + 128 * ii0 - i0 - 1, -T + N + i0 + 1); i1 < N; i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
           }
         } else if (N >= 255 && 128 * ii0 + 64 == T) {
           for (register int i0 = T - 64; i0 < T - 1; i0 += 1) {
             for (register int i1 = 0; i1 <= min(T - i0 - 2, -T + i0 + 64); i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
             }
             for (register int i1 = max(T - i0 + 63, -T + i0 + 129); i1 <= min(T - i0 + 126, -T + i0 + 192); i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
             for (register int i1 = max(T + N - i0 - 65, -T + N + i0 + 1); i1 < N; i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
           }
         }
         if (N >= 255 && T >= 128 * ii0 + 43 && 128 * ii0 + 127 >= T) {
           if (T + 5312 >= 21 * N + 128 * ii0) {
             #pragma omp parallel for
             for (register int ii1 = 0; ii1 <= min(min(-N + 256, -6 * ii0 + (T + N - 2 * ii0 - 4) / 21 - 14), 2 * ii0 - (T + 19 * N - 2 * ii0 + 48) / 63 + 79); ii1 += 1) {
               if (N == 255 && ii1 == 0 && (T + 1) % 2 == 0) {
                 A[((T + 1) / 2) + 64 * ii0 - 2 * ((T + 128 * ii0 + 1) / 4)][((T - 1) / 2) - 64 * ii0] = (0.125 * ((A[((T - 1) / 2) + 64 * ii0 - 2 * ((T + 128 * ii0 - 1) / 4)][((T + 1) / 2) - 64 * ii0] - (2.0 * A[((T - 1) / 2) + 64 * ii0 - 2 * ((T + 128 * ii0 - 1) / 4)][((T - 1) / 2) - 64 * ii0])) + A[((T - 1) / 2) + 64 * ii0 - 2 * ((T + 128 * ii0 - 1) / 4)][((T - 3) / 2) - 64 * ii0]));
                 if (T >= 128 * ii0 + 65) {
                   A[((T + 1) / 2) + 64 * ii0 - 2 * ((T + 128 * ii0 + 1) / 4)][((-T + 255) / 2) + 64 * ii0] = (0.125 * ((A[((T - 1) / 2) + 64 * ii0 - 2 * ((T + 128 * ii0 - 1) / 4)][((-T + 257) / 2) + 64 * ii0] - (2.0 * A[((T - 1) / 2) + 64 * ii0 - 2 * ((T + 128 * ii0 - 1) / 4)][((-T + 255) / 2) + 64 * ii0])) + A[((T - 1) / 2) + 64 * ii0 - 2 * ((T + 128 * ii0 - 1) / 4)][((-T + 253) / 2) + 64 * ii0]));
                 } else {
                   A[((T + 1) / 2) + 64 * ii0 - 2 * ((T + 128 * ii0 + 1) / 4)][((-T + 255) / 2) + 64 * ii0] = (0.125 * ((A[((T - 1) / 2) + 64 * ii0 - 2 * ((T + 128 * ii0 - 1) / 4)][((-T + 257) / 2) + 64 * ii0] - (2.0 * A[((T - 1) / 2) + 64 * ii0 - 2 * ((T + 128 * ii0 - 1) / 4)][((-T + 255) / 2) + 64 * ii0])) + A[((T - 1) / 2) + 64 * ii0 - 2 * ((T + 128 * ii0 - 1) / 4)][((-T + 253) / 2) + 64 * ii0]));
                 }
               }
               if (128 * ii0 + 64 >= T && ii1 == 0) {
                 for (register int i0 = -N + 64 * ii0 + (T + N) / 2 + 128; i0 < T; i0 += 1) {
                   for (register int i1 = T - i0 - 1; i1 <= -128 * ii0 + i0; i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i0 + 1 == T && i1 == 0 ? N - 1 : i1 - 1]));
                   }
                   for (register int i1 = 128 * ii0 - i0 + 127; i1 <= min(-T + i0 + 128, T - N - 389 * ii0 + 2 * i0 + (-5 * ii0 + 5 * i0 + 5) / 127 + 257); i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                   if (N == 255) {
                     for (register int i1 = max(128 * ii0 - i0 + 127, T - 389 * ii0 + 2 * i0 + (-5 * ii0 + 5 * i0 + 5) / 127 + 3); i1 <= -T + i0 + 128; i1 += 1) {
                       A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                     }
                   }
                 }
               } else {
                 for (register int i0 = -N + 64 * ii0 - ii1 + (T + N + ii1) / 2 + 128; i0 < T; i0 += 1) {
                   if (ii1 == 0) {
                     for (register int i1 = T - i0 - 1; i1 <= -128 * ii0 + i0; i1 += 1) {
                       A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i0 + 1 == T && i1 == 0 ? N - 1 : i1 - 1]));
                     }
                     for (register int i1 = max(-128 * ii0 + i0 + 1, 128 * ii0 - i0 + 127); i1 <= -T + i0 + 128; i1 += 1) {
                       A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                     }
                   } else {
                     for (register int i1 = T - i0 + 127; i1 <= min(-128 * ii0 + i0 + 128, T - 389 * ii0 + 2 * i0 + (-5 * ii0 + 5 * i0 + 5) / 127 + 2); i1 += 1) {
                       A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                     }
                     for (register int i1 = max(T - i0 + 127, T - 389 * ii0 + 2 * i0 + (-5 * ii0 + 5 * i0 + 5) / 127 + 3); i1 <= -128 * ii0 + i0 + 128; i1 += 1) {
                       A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                     }
                     for (register int i1 = max(-128 * ii0 + i0 + 129, 128 * ii0 - i0 + 254); i1 <= -T + i0 + 255; i1 += 1) {
                       A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i0 + 1 == T && i1 == 254 ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                     }
                   }
                 }
               }
             }
             if (N == 255 && 128 * ii0 + 127 == T) {
               for (register int i0 = T - 64; i0 < T; i0 += 1) {
                 for (register int i1 = T - i0 + 127; i1 <= -T + i0 + 255; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i0 + 1 == T && i1 == 254 ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             }
             if (N == 255 && T >= 128 * ii0 + 84 && 128 * ii0 + 126 >= T) {
               for (register int i0 = 64 * ii0 + T / 2; i0 < T; i0 += 1) {
                 for (register int i1 = T - i0 + 127; i1 <= -128 * ii0 + i0 + 128; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
                 for (register int i1 = max(-128 * ii0 + i0 + 129, 128 * ii0 - i0 + 254); i1 <= -T + i0 + 255; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i0 + 1 == T && i1 == 254 ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             }
           } else {
             for (register int i0 = 64 * ii0 + T / 2; i0 < T; i0 += 1) {
               for (register int i1 = T - i0 - 1; i1 <= -128 * ii0 + i0; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i0 + 1 == T && i1 == 0 ? 255 : i1 - 1]));
               }
               for (register int i1 = 128 * ii0 - i0 + 127; i1 <= -T + i0 + 128; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           }
           if (N == 256 && T >= 128 * ii0 + 63) {
             for (register int i0 = 64 * ii0 + T / 2; i0 < T; i0 += 1) {
               for (register int i1 = T - i0 + 127; i1 <= min(-128 * ii0 + i0 + 128, T - 389 * ii0 + 2 * i0 + (-5 * ii0 + 5 * i0 + 5) / 127 + 1); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
               for (register int i1 = max(T - i0 + 127, T - 389 * ii0 + 2 * i0 + (-5 * ii0 + 5 * i0 + 5) / 127 + 2); i1 <= -128 * ii0 + i0 + 128; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
               for (register int i1 = max(-128 * ii0 + i0 + 129, 128 * ii0 - i0 + 255); i1 <= -T + i0 + 256; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i0 + 1 == T && i1 == 255 ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           } else if (128 * ii0 + 318 >= T + N) {
             for (register int i0 = 64 * ii0 + T / 2; i0 < T; i0 += 1) {
               for (register int i1 = T - i0 + 127; i1 <= min(-128 * ii0 + i0 + 128, T - N - 389 * ii0 + 2 * i0 + (-5 * ii0 + 5 * i0 + 5) / 127 + 257); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
               if (N == 255 && 128 * ii0 + 63 == T) {
                 for (register int i1 = max(T - i0 + 127, -2 * T + 2 * i0 - (5 * T - 5 * i0 + 187) / 127 + 195); i1 <= -T + i0 + 191; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               } else {
                 for (register int i1 = max(T - i0 + 127, T - N - 389 * ii0 + 2 * i0 + (-5 * ii0 + 5 * i0 + 5) / 127 + 258); i1 <= -128 * ii0 + i0 + 128; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
               for (register int i1 = N + 128 * ii0 - i0 - 1; i1 <= -T + N + i0; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i0 + 1 == T && i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           }
         } else if (N == 256 && T >= 128 * ii0 + 2 && 128 * ii0 + 42 >= T) {
           #pragma omp parallel for
           for (register int ii1 = 0; ii1 <= 1; ii1 += 1) {
             for (register int i0 = 64 * ii0 + T / 2; i0 < T; i0 += 1) {
               if (ii1 == 1) {
                 for (register int i1 = T - i0 + 127; i1 <= -128 * ii0 + i0 + 128; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
                 for (register int i1 = 128 * ii0 - i0 + 255; i1 <= -T + i0 + 256; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i0 + 1 == T && i1 == 255 ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               } else {
                 for (register int i1 = T - i0 - 1; i1 <= -128 * ii0 + i0; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i0 + 1 == T && i1 == 0 ? 255 : i1 - 1]));
                 }
                 for (register int i1 = 128 * ii0 - i0 + 127; i1 <= -T + i0 + 128; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             }
           }
         } else if (N == 255 && 128 * ii0 + 42 >= T) {
           #pragma omp parallel for
           for (register int ii1 = 0; ii1 <= 1; ii1 += 1) {
             for (register int i0 = 64 * ii0 + T / 2; i0 < T; i0 += 1) {
               if (ii1 == 1) {
                 for (register int i1 = T - i0 + 127; i1 <= -128 * ii0 + i0 + 128; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
                 for (register int i1 = 128 * ii0 - i0 + 254; i1 <= -T + i0 + 255; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i0 + 1 == T && i1 == 254 ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               } else {
                 for (register int i1 = T - i0 - 1; i1 <= -128 * ii0 + i0; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i0 + 1 == T && i1 == 0 ? 254 : i1 - 1]));
                 }
                 for (register int i1 = 128 * ii0 - i0 + 127; i1 <= -T + i0 + 128; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             }
           }
         } else if (N <= 254 && N + 128 * ii0 >= T + 130 && T >= 128 * ii0 + 65) {
           for (register int i0 = 128 * ii0; i0 < T - 1; i0 += 1) {
             for (register int i1 = 0; i1 <= min(T - i0 - 2, -128 * ii0 + i0); i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
             }
             for (register int i1 = max(128 * ii0 - i0 + 127, -T + i0 + 129); i1 <= min(min(T - i0 + 126, -128 * ii0 + i0 + 128), T - N - 389 * ii0 + 2 * i0 + (-5 * ii0 + 5 * i0 + 5) / 127 + 257); i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
             for (register int i1 = max(128 * ii0 - i0 + 127, T - N - 389 * ii0 + 2 * i0 + (-5 * ii0 + 5 * i0 + 5) / 127 + 258); i1 <= min(T - i0 + 126, -128 * ii0 + i0 + 128); i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
             for (register int i1 = max(N + 128 * ii0 - i0 - 1, -T + N + i0 + 1); i1 < N; i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
           }
         } else if (N <= 254 && N + 128 * ii0 >= T + 128 && T >= 128 * ii0 + 64 && T + 129 >= N + 128 * ii0) {
           for (register int i0 = 128 * ii0; i0 < T - 1; i0 += 1) {
             for (register int i1 = 0; i1 <= min(T - i0 - 2, -128 * ii0 + i0); i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
             }
             if (N + 128 * ii0 == T + 128) {
               for (register int i1 = max(-T + i0 + 129, T - N - i0 + 255); i1 <= min(T - i0 + 126, -3 * T + 3 * N + i0 + 128 * ((-T + N + 128 * i0 + 8127) / 8128) - 384); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             } else {
               for (register int i1 = max(-T + i0 + 129, T - N - i0 + 256); i1 <= min(-T + N + i0 - 1, T - i0 + 126); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
             for (register int i1 = max(N + 128 * ii0 - i0 - 1, -T + N + i0 + 1); i1 < N; i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
           }
         } else if (N + 128 * ii0 >= T + 128 && 128 * ii0 + 255 >= T + N) {
           for (register int i0 = 128 * ii0; i0 < T - 1; i0 += 1) {
             for (register int i1 = 0; i1 <= min(T - i0 - 2, -128 * ii0 + i0); i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
             }
             for (register int i1 = max(128 * ii0 - i0 + 127, -T + i0 + 129); i1 <= min(T - i0 + 126, -128 * ii0 + i0 + 128); i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
             for (register int i1 = max(N + 128 * ii0 - i0 - 1, -T + N + i0 + 1); i1 < N; i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
           }
         }
         for (register int k = ii0 - (T - N + 127 * ii0 + 127) / 255 + 2; k <= -ii0 + (T + 126 * ii0 + 125) / 254 + 2; k += 1) {
           #pragma omp parallel for
           for (register int ii1 = 0; ii1 <= min(k - 2, 2 * ii0 - (2 * T - N + 130) / 128 + 2); ii1 += 1) {
             if (128 * ii0 + 63 * k + 2 >= T && ii1 == 0) {
               if (k == 3) {
                 for (register int i0 = 128 * ii0 + 64; i0 < T; i0 += 1) {
                   for (register int i1 = max(-128 * ii0 + i0 - 127, 128 * ii0 - i0 + 127); i1 <= min(-128 * ii0 + i0, 128 * ii0 - i0 + 254); i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i0 == 128 * ii0 + 127 && i1 == 0 ? N - 1 : i1 - 1]));
                   }
                 }
               } else if (N <= 254) {
                 for (register int i0 = 128 * ii0; i0 < 64 * ii0 + (T + N + 1) / 2 - 65; i0 += 1) {
                   for (register int i1 = 0; i1 <= min(T - i0 - 2, -128 * ii0 + i0); i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
                   }
                   for (register int i1 = max(128 * ii0 - i0 + 127, -T + i0 + 129); i1 <= min(N + 128 * ii0 - i0 - 2, -T + N + i0); i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                   for (register int i1 = max(-T + N + i0 + 1, 128 * ii0 - i0 + 127); i1 <= min(N + 128 * ii0 - i0 - 2, -128 * ii0 + i0 + 128); i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                   for (register int i1 = N + 128 * ii0 - i0 - 1; i1 <= min(N - 1, T - i0 + 126); i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                   for (register int i1 = max(-T + N + i0 + 1, T - i0 + 127); i1 < N; i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                 }
                 for (register int i0 = 64 * ii0 + (T + N + 1) / 2 - 65; i0 < T - 1; i0 += 1) {
                   for (register int i1 = 0; i1 < T - i0 - 1; i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
                   }
                   for (register int i1 = -T + i0 + 129; i1 <= min(-T + N + i0, T - i0 + 126); i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                   for (register int i1 = -T + N + i0 + 1; i1 < N; i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                 }
               } else {
                 for (register int i0 = T - 128; i0 < T - 1; i0 += 1) {
                   for (register int i1 = 0; i1 <= min(T - i0 - 2, -T + i0 + 128); i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? 254 : i1 - 1]));
                   }
                   for (register int i1 = max(T - i0 - 1, -T + i0 + 129); i1 <= min(T - i0 + 125, -T + i0 + 255); i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                   if (T >= i0 + 66) {
                     A[(i0 + 1) % 2][-T + i0 + 256] = (0.125 * ((A[i0 % 2][-T + i0 + 257] - (2.0 * A[i0 % 2][-T + i0 + 256])) + A[i0 % 2][-T + i0 + 255]));
                   }
                   A[(i0 + 1) % 2][T - i0 + 126] = (0.125 * ((A[i0 % 2][i0 + 128 == T ? 0 : T - i0 + 127] - (2.0 * A[i0 % 2][T - i0 + 126])) + A[i0 % 2][T - i0 + 125]));
                   for (register int i1 = max(T - i0 + 127, -T + i0 + 256); i1 <= 254; i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 == 254 ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                 }
               }
             } else if (k == 3 && ii1 == 1) {
               for (register int i0 = -N + 128 * ii0 + N / 2 + 192; i0 < T; i0 += 1) {
                 for (register int i1 = max(-128 * ii0 + i0 + 1, 128 * ii0 - i0 + 255); i1 <= min(N - 128 * ii0 + i0 - 128, N + 128 * ii0 - i0 + 126); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i0 == 128 * ii0 + 127 && i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             } else {
               for (register int i0 = 128 * ii0; i0 <= 128 * ii0 + 126; i0 += 1) {
                 for (register int i1 = 0; i1 <= min(-128 * ii0 + i0, 128 * ii0 - i0 + 126); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
                 }
                 if (N == 238 && 128 * ii0 + 129 == T) {
                   for (register int i1 = max(T - i0 - 2, -T + i0 + 130); i1 <= min(T - i0 + 107, -3 * T + i0 + 128 * ((-T + 128 * i0 + 8384) / 8128) + 387); i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                 } else if (N == 254 && 128 * ii0 + 129 == T) {
                   for (register int i1 = max(T - i0 - 2, -T + i0 + 130); i1 <= min(T - i0 + 123, -T + i0 + 255); i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                   if (T >= i0 + 67) {
                     for (register int i1 = -T + i0 + 256; i1 <= -T + i0 + 257; i1 += 1) {
                       A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                     }
                   }
                 } else if (N >= 240 && 128 * ii0 + 129 == T) {
                   for (register int i1 = max(T - i0 - 2, -T + i0 + 130); i1 <= min(T + N - i0 - 131, -3 * T + i0 + 128 * ((-T + 128 * i0 + 8384) / 8128) + 387); i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                 } else if (128 * ii0 + 129 == T && 2 * ((43 * N - 43) / 170) + 119 >= N) {
                   for (register int i1 = max(T - i0 - 2, -T + i0 + 130); i1 <= min(T + N - i0 - 131, -3 * T + i0 + 128 * ((-T + 128 * i0 + 8384) / 8128) + 387); i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                 } else {
                   for (register int i1 = max(-128 * ii0 + i0 + 1, 128 * ii0 - i0 + 127); i1 <= min(N + 128 * ii0 - i0 - 2, -384 * ii0 + i0 + 128 * ((-2 * ii0 + 2 * i0 + 128) / 127)); i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                 }
                 for (register int i1 = max(N + 128 * ii0 - i0 - 1, -128 * ii0 + i0 + 1); i1 <= min(N - 128 * ii0 + i0 - 128, 128 * ii0 - i0 + 254); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
                 for (register int i1 = max(N - 128 * ii0 + i0 - 127, N + 128 * ii0 - i0 - 1); i1 < N; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             }
           }
           if (2 * T >= N + 256 * ii0 + 126 && k == 3) {
             for (register int i0 = -N + 128 * ii0 + N / 2 + 192; i0 <= 128 * ii0 + (N + 1) / 2 + 62; i0 += 1) {
               for (register int i1 = max(-128 * ii0 + i0 + 1, 128 * ii0 - i0 + 255); i1 <= min(N - 128 * ii0 + i0 - 128, N + 128 * ii0 - i0 + 126); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i0 == 128 * ii0 + 127 && i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           }
         }
         if (N <= 254 && 128 * ii0 + 128 >= T) {
           if (128 * ii0 + 256 >= T + N) {
             for (register int i0 = 64 * ii0 + T / 2; i0 < T; i0 += 1) {
               for (register int i1 = T - i0 - 1; i1 <= -128 * ii0 + i0; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i0 + 1 == T && i1 == 0 ? N - 1 : i1 - 1]));
               }
               for (register int i1 = max(-128 * ii0 + i0 + 1, 128 * ii0 - i0 + 127); i1 <= min(min(N + 128 * ii0 - i0 - 2, -T + i0 + 128), T - N - i0 + 254); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
               for (register int i1 = max(N + 128 * ii0 - i0 - 1, -128 * ii0 + i0 + 1); i1 <= min(-T + i0 + 128, T - N - i0 + 254); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
               if (128 * ii0 + 64 >= T) {
                 for (register int i1 = max(128 * ii0 - i0 + 127, T - N - i0 + 255); i1 <= min(min(N + 128 * ii0 - i0 - 2, -T + i0 + 128), T - N - 389 * ii0 + 2 * i0 + (-5 * ii0 + 5 * i0 + 5) / 127 + 257); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
                 for (register int i1 = max(N + 128 * ii0 - i0 - 1, T - N - i0 + 255); i1 <= -T + i0 + 128; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               } else {
                 for (register int i1 = T - N - i0 + 255; i1 <= -T + i0 + 128; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
               for (register int i1 = max(128 * ii0 - i0 + 127, T - N - 389 * ii0 + 2 * i0 + (-5 * ii0 + 5 * i0 + 5) / 127 + 258); i1 <= -T + i0 + 128; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           } else if (128 * ii0 + 64 >= T) {
             for (register int i0 = 64 * ii0 + T / 2; i0 < T; i0 += 1) {
               for (register int i1 = T - i0 - 1; i1 <= -128 * ii0 + i0; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i0 + 1 == T && i1 == 0 ? N - 1 : i1 - 1]));
               }
               for (register int i1 = 128 * ii0 - i0 + 127; i1 <= min(-T + i0 + 128, T - N - 389 * ii0 + 2 * i0 + (-5 * ii0 + 5 * i0 + 5) / 127 + 257); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
               for (register int i1 = max(128 * ii0 - i0 + 127, T - N - 389 * ii0 + 2 * i0 + (-5 * ii0 + 5 * i0 + 5) / 127 + 258); i1 <= -T + i0 + 128; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           }
           #pragma omp parallel for
           for (register int ii1 = max(2 * ii0 - (T + 63) / 64 + 2, ii0 - (T + N - 2 * ii0 + 121) / 126 + 3); ii1 <= ii0 - (T - N + 128) / 128; ii1 += 1) {
             if (ii1 == 0) {
               for (register int i0 = 64 * ii0 + T / 2; i0 <= -N + 64 * ii0 + (T + N) / 2 + 127; i0 += 1) {
                 for (register int i1 = T - i0 - 1; i1 <= -128 * ii0 + i0; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
                 for (register int i1 = max(-128 * ii0 + i0 + 1, 128 * ii0 - i0 + 127); i1 <= -T + i0 + 128; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             }
             for (register int i0 = max(64 * ii0 + T / 2, -N + 64 * ii0 - 63 * ii1 + (T + N) / 2 + 128); i0 < T; i0 += 1) {
               if (ii1 == 1) {
                 for (register int i1 = T - i0 + 127; i1 <= min(min(N + 128 * ii0 - i0 - 2, -128 * ii0 + i0 + 128), T - N - 389 * ii0 + 2 * i0 + (-5 * ii0 + 5 * i0 + 5) / 127 + 257); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
                 for (register int i1 = max(T - i0 + 127, T - N - 389 * ii0 + 2 * i0 + (-5 * ii0 + 5 * i0 + 5) / 127 + 258); i1 <= -128 * ii0 + i0 + 128; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
                 for (register int i1 = N + 128 * ii0 - i0 - 1; i1 <= -T + N + i0; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i0 + 1 == T && i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               } else {
                 for (register int i1 = T - i0 - 1; i1 <= -128 * ii0 + i0; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i0 + 1 == T && i1 == 0 ? N - 1 : i1 - 1]));
                 }
                 for (register int i1 = max(-128 * ii0 + i0 + 1, 128 * ii0 - i0 + 127); i1 <= -T + i0 + 128; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             }
           }
           if (T + 128 >= N + 128 * ii0) {
             for (register int i0 = T - N + N / 2 + 64; i0 < T; i0 += 1) {
               for (register int i1 = T - i0 + 127; i1 <= -T + N + i0; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i0 + 1 == T && i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           }
         } else if (N == 255 && 128 * ii0 + 128 == T) {
           #pragma omp parallel for
           for (register int ii1 = 0; ii1 <= 1; ii1 += 1) {
             for (register int i0 = T - 64; i0 < T; i0 += 1) {
               if (ii1 == 1) {
                 for (register int i1 = T - i0 + 127; i1 <= -T + i0 + 255; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i0 + 1 == T && i1 == 254 ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               } else {
                 for (register int i1 = T - i0 - 1; i1 <= -T + i0 + 128; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i0 + 1 == T && i1 == 0 ? 254 : i1 - 1]));
                 }
               }
             }
           }
         }
       }
       if (N == 256 && 128 * ii0 + 1 == T) {
         #pragma omp parallel for
         for (register int ii1 = 0; ii1 <= 1; ii1 += 1) {
           if (ii1 == 0) {
             A[1][0] = (0.125 * ((A[0][1] - (2.0 * A[0][0])) + A[0][255]));
           }
           A[1][ii1 + 127] = (0.125 * ((A[0][ii1 + 128] - (2.0 * A[0][ii1 + 127])) + A[0][ii1 + 126]));
           if (ii1 == 1) {
             A[1][255] = (0.125 * ((A[0][0] - (2.0 * A[0][255])) + A[0][254]));
           }
         }
       }
     }
     if (T + N <= 256 && 510 * T >= 193 * N + 16382) {
       for (register int k = 2; k <= 3; k += 1) {
         #pragma omp parallel for
         for (register int ii1 = 0; ii1 < k - 1; ii1 += 1) {
           if (k == 3) {
             if (ii1 == 0) {
               for (register int i0 = T - 64; i0 <= T - N + N / 2 + 63; i0 += 1) {
                 for (register int i1 = T - i0 - 1; i1 <= -T + i0 + 128; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             }
             for (register int i0 = T - N + N / 2 + 64; i0 < T; i0 += 1) {
               if (ii1 == 1) {
                 for (register int i1 = T - i0 + 127; i1 <= -T + N + i0; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i0 + 1 == T && i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               } else {
                 for (register int i1 = T - i0 - 1; i1 <= -T + i0 + 128; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i0 + 1 == T && i1 == 0 ? N - 1 : i1 - 1]));
                 }
               }
             }
           } else {
             for (register int i0 = 0; i0 < T - 1; i0 += 1) {
               for (register int i1 = 0; i1 < T - i0 - 1; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
               }
               for (register int i1 = max(T - i0 - 1, -T + i0 + 129); i1 <= min(-T + N + i0, T - i0 + 126); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
               for (register int i1 = max(T - i0 - 1, -T + N + i0 + 1); i1 < N; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           }
         }
       }
     }
   }
   if (T >= 65 && T <= 128 && T + N >= 257) {
     if (N >= 257) {
       for (register int i0 = 0; i0 <= T - N + (N + 1) / 2 + 126; i0 += 1) {
         for (register int i1 = 0; i1 <= T - N - i0 + 254; i1 += 1) {
           A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
         }
         for (register int i1 = max(-T + i0 + 129, T - N - i0 + 255); i1 < T - i0 - 1; i1 += 1) {
           A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
         }
         for (register int i1 = -T + i0 + 257; i1 <= T - i0 + 126; i1 += 1) {
           A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
         }
         for (register int i1 = max(-T + N + i0 + 1, T - i0 + 127); i1 <= min(N - 1, T - i0 + 254); i1 += 1) {
           A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
         }
       }
     }
     for (register int k = 2; k <= 3; k += 1) {
       if (k == 3) {
         #pragma omp parallel for
         for (register int ii1 = 0; ii1 <= 1; ii1 += 1) {
           if (ii1 == 1) {
             for (register int i0 = max(T - 64, T - N + N / 2 + 64); i0 < T; i0 += 1) {
               for (register int i1 = T - i0 + 127; i1 <= min(-T + N + i0, -T + i0 + 256); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i0 + 1 == T && i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           } else if (N <= 255) {
             for (register int i0 = T - 64; i0 < T; i0 += 1) {
               for (register int i1 = T - i0 - 1; i1 <= -T + i0 + 128; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i0 + 1 == T && i1 == 0 ? N - 1 : i1 - 1]));
               }
             }
           } else {
             for (register int i0 = T - 64; i0 < T; i0 += 1) {
               for (register int i1 = T - i0 - 1; i1 <= -T + i0 + 128; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i0 + 1 == T && i1 == 0 ? N - 1 : i1 - 1]));
               }
             }
           }
         }
         for (register int i0 = T - N + N / 2 + 128; i0 < T; i0 += 1) {
           for (register int i1 = T - i0 + 255; i1 <= -T + N + i0; i1 += 1) {
             A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i0 + 1 == T && i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
           }
         }
       } else if (N >= 257) {
         #pragma omp parallel for
         for (register int ii1 = 0; ii1 <= 1; ii1 += 1) {
           for (register int ii1_2 = ii1 + 1; ii1_2 <= 2; ii1_2 += 1) {
             if (ii1_2 == 2) {
               for (register int i0 = max(0, T - N + N / 2 + 64); i0 < T - 1; i0 += 1) {
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
               for (register int i0 = 0; i0 < T - 1; i0 += 1) {
                 for (register int i1 = max(T - i0 - 1, -T + i0 + 129); i1 <= min(T - i0 + 126, -T + i0 + 256); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             }
           }
         }
       } else {
         for (register int i0 = 0; i0 < T - 1; i0 += 1) {
           for (register int i1 = 0; i1 < T - i0 - 1; i1 += 1) {
             A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
           }
           for (register int i1 = max(T - i0 - 1, -T + i0 + 129); i1 <= min(-T + N + i0, T - i0 + 126); i1 += 1) {
             A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
           }
           for (register int i1 = max(T - i0 - 1, -T + N + i0 + 1); i1 < N; i1 += 1) {
             A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
           }
         }
       }
     }
   } else if (T >= 129 && N >= 257) {
     for (register int ii0 = 0; ii0 <= (81727 * T - 16256 * N + 1999041) / 10461056; ii0 += 1) {
       if (ii0 >= 1) {
         for (register int i0 = 128 * ii0; i0 <= min(T - N + (N + 1) / 2 + 126, -N + 128 * ii0 + (N + 1) / 2 + 254); i0 += 1) {
           for (register int i1 = 0; i1 <= min(min(-128 * ii0 + i0, T - N - i0 + 254), -N + 128 * ii0 - i0 + 382); i1 += 1) {
             A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
           }
           for (register int i1 = max(max(N - 128 * ii0 + i0 - 127, -T + N + i0 + 1), 128 * ii0 - i0 + 255); i1 <= min(min(N + 128 * ii0 - i0 - 2, -382 * ii0 + 3 * i0 - (-2 * ii0 + 2 * i0 + 129) / 127 + 130), T - N - 389 * ii0 + 2 * i0 + (-5 * ii0 + 5 * i0 + 5) / 127 + 257); i1 += 1) {
             A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
           }
           for (register int i1 = max(max(max(N - 128 * ii0 + i0 - 127, -T + N + i0 + 1), 128 * ii0 - i0 + 255), T - N - 389 * ii0 + 2 * i0 + (-5 * ii0 + 5 * i0 + 5) / 127 + 258); i1 <= min(min(N + 128 * ii0 - i0 - 2, -128 * ii0 + i0 + 256), -2 * T + 2 * N - 369 * ii0 + 5 * i0 + floord(8 * T - 8 * N + 19 * ii0 - 16 * i0 + 2, 67) - 109); i1 += 1) {
             A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
           }
           for (register int i1 = max(max(max(N - 128 * ii0 + i0 - 127, 128 * ii0 - i0 + 255), -2 * T + 2 * N - 369 * ii0 + 5 * i0 + floord(8 * T - 8 * N + 19 * ii0 - 16 * i0 + 2, 67) - 108), -382 * ii0 + 3 * i0 - (-2 * ii0 + 2 * i0 + 129) / 127 + 131); i1 <= min(N + 128 * ii0 - i0 - 2, -128 * ii0 + i0 + 256); i1 += 1) {
             A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
           }
           for (register int i1 = max(max(N - 128 * ii0 + i0 - 127, N + 128 * ii0 - i0 - 1), -T + N + i0 + 1); i1 <= min(min(N - 1, T - i0 + 254), 128 * ii0 - i0 + 382); i1 += 1) {
             A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
           }
         }
       } else {
         for (register int i0 = 0; i0 <= -N + (N + 1) / 2 + 254; i0 += 1) {
           for (register int i1 = 0; i1 <= -N - i0 + 382; i1 += 1) {
             A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
           }
           for (register int i1 = max(i0 + 1, -N - i0 + 383); i1 <= -i0 + 126; i1 += 1) {
             A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
           }
           for (register int i1 = i0 + 129; i1 <= -i0 + 254; i1 += 1) {
             A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
           }
           for (register int i1 = max(N + i0 - 127, -i0 + 255); i1 <= min(N - 1, -i0 + 382); i1 += 1) {
             A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
           }
         }
       }
       {
         if (127 * N + 49408 * ii0 >= 386 * T + 15998) {
           for (register int ii1_2 = 1; ii1_2 <= 2; ii1_2 += 1) {
             if (ii1_2 == 2) {
               for (register int i0 = max(128 * ii0, -N + 64 * ii0 + (T + N) / 2 + 128); i0 < T - 1; i0 += 1) {
                 for (register int i1 = max(0, T - N - i0 + 255); i1 <= min(T - i0 - 2, -128 * ii0 + i0); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
                 }
                 for (register int i1 = max(max(N + 128 * ii0 - i0 - 1, -T + N + i0 + 1), T - i0 + 255); i1 < N; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             } else {
               for (register int i0 = 128 * ii0; i0 < T - 1; i0 += 1) {
                 for (register int i1 = max(128 * ii0 - i0 + 127, -T + i0 + 129); i1 <= min(T - i0 + 126, -128 * ii0 + i0 + 128); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             }
           }
         } else if (ii0 == 0) {
           for (register int ii1_2 = 1; ii1_2 <= 2; ii1_2 += 1) {
             if (ii1_2 == 2) {
               if (N >= T + 129) {
                 for (register int i0 = -N + N / 2 + 192; i0 <= -N + (T + N) / 2 + 127; i0 += 1) {
                   for (register int i1 = max(0, -N - i0 + 383); i1 <= i0; i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
                   }
                   for (register int i1 = -i0 + 383; i1 < N; i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                 }
                 for (register int i0 = -N + (T + N) / 2 + 128; i0 < (T + 1) / 2 - 1; i0 += 1) {
                   for (register int i1 = max(0, -N - i0 + 383); i1 <= min(i0, -i0 + 126); i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
                   }
                   for (register int i1 = max(N + i0 - 127, -i0 + 383); i1 < N; i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                 }
               } else {
                 for (register int i0 = -N + N / 2 + 192; i0 <= min(126, floord(T + 1, 2) - 2); i0 += 1) {
                   for (register int i1 = max(0, -N - i0 + 383); i1 <= min(i0, -i0 + 126); i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
                   }
                   for (register int i1 = max(N + i0 - 127, -i0 + 383); i1 < N; i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                 }
               }
               for (register int i0 = (T + 1) / 2 - 1; i0 <= min(126, T - N + 254); i0 += 1) {
                 for (register int i1 = max(0, -N - i0 + 383); i1 <= -i0 + 126; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
                 }
                 for (register int i1 = max(N + i0 - 127, -i0 + 383); i1 < N; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
               for (register int i0 = max(T - N + 255, (T + 1) / 2 - 1); i0 <= 126; i0 += 1) {
                 for (register int i1 = 0; i1 <= -i0 + 126; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
                 }
                 for (register int i1 = max(N + i0 - 127, -i0 + 383); i1 < N; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             } else {
               for (register int i0 = 0; i0 <= min(126, T - 129); i0 += 1) {
                 for (register int i1 = max(i0 + 1, -i0 + 127); i1 <= min(i0 + 128, -i0 + 254); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
               for (register int i0 = T - 128; i0 <= 126; i0 += 1) {
                 for (register int i1 = max(i0 + 1, -i0 + 127); i1 <= min(i0 + 128, -i0 + 254); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             }
           }
         }
         #pragma omp parallel for
         for (register int ii1 = max(max(0, -ii0 + 1), 2 * ii0 - (386 * T - 127 * N + 15230 * ii0 + 48316) / 32319 + 1); ii1 <= min(1, -ii0 + (T - 1) / 128); ii1 += 1) {
           {
             if (ii1 == 0) {
               for (register int i0 = 128 * ii0; i0 < -2 * T + N + 374 * ii0 + (10 * T - 5 * N + 40 * ii0 - 26) / 132 - 122; i0 += 1) {
                 for (register int i1 = max(128 * ii0 - i0 + 127, -T + i0 + 129); i1 <= min(T - i0 + 126, -128 * ii0 + i0 + 128); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
               for (register int i0 = max(128 * ii0, -2 * T + N + 374 * ii0 + (10 * T - 5 * N + 40 * ii0 - 26) / 132 - 122); i0 < 170 * ii0 + floord(-127 * T + 127 * N + 44 * ii0 + 82, 386) - 42; i0 += 1) {
                 for (register int i1 = 128 * ii0 - i0 + 127; i1 <= -128 * ii0 + i0 + 128; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             }
             for (register int i0 = max(max(max(128 * ii0, -N + 128 * ii0 + 64 * ii1 + N / 2 + 128), -2 * T + N + 374 * ii0 - 245 * ii1 + floord(10 * T - 5 * N + 40 * ii0 + 21 * ii1 - 26, 132) - 122), 170 * ii0 - 42 * ii1 + floord(-127 * T + 127 * N + 44 * ii0 + 149 * ii1 + 82, 386) - 42); i0 <= min(min(T - 2, 128 * ii0 + 126), 64 * ii0 - 64 * ii1 + (T + 1) / 2 + 62); i0 += 1) {
               for (register int i1 = max(max(-128 * ii0 + 128 * ii1 + i0 + 1, 128 * ii0 + 128 * ii1 - i0 + 127), -T + i0 + 129); i1 <= min(min(min(min(N - 128 * ii0 + i0 - 128, T - i0 + 126), -128 * ii0 + 128 * ii1 + i0 + 128), 128 * ii0 + 128 * ii1 - i0 + 254), T - N - 389 * ii0 + 254 * ii1 + 2 * i0 + (-5 * ii0 + 61 * ii1 + 5 * i0 + 5) / 127 + 257); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
               if (ii1 == 0) {
                 for (register int i1 = T - N - 389 * ii0 + 2 * i0 + (-5 * ii0 + 5 * i0 + 5) / 127 + 258; i1 <= min(min(T - i0 + 126, -128 * ii0 + i0 + 128), 128 * ii0 - i0 + 254); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               } else {
                 if (ii0 == 0) {
                   for (register int i1 = T - i0 + 127; i1 <= -T + N + i0; i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                 } else {
                   for (register int i1 = T - i0 + 127; i1 <= -T + N + i0; i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                 }
                 for (register int i1 = max(-T + N + i0 + 1, T - i0 + 127); i1 <= min(N - 128 * ii0 + i0 - 128, 128 * ii0 - i0 + 382); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             }
             if (N == 384 && ii0 >= 1 && ii1 == 1) {
               for (register int i0 = 64 * ii0 + (T + 1) / 2 - 1; i0 <= 170 * ii0 - (127 * T - 44 * ii0 + 210) / 386 + 167; i0 += 1) {
                 for (register int i1 = -128 * ii0 + i0 + 129; i1 <= 128 * ii0 - i0 + 382; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
               for (register int i0 = max(64 * ii0 + (T + 1) / 2 - 1, 170 * ii0 - (127 * T - 44 * ii0 + 210) / 386 + 168); i0 <= 128 * ii0 + 126; i0 += 1) {
                 for (register int i1 = -128 * ii0 + i0 + 129; i1 <= 128 * ii0 - i0 + 382; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             }
             if (ii1 == 1) {
               for (register int i0 = 64 * ii0 + (T + 1) / 2 - 1; i0 <= min(128 * ii0 + 126, -N + 64 * ii0 + (T + N + 1) / 2 + 190); i0 += 1) {
                 if (ii0 == 0) {
                   for (register int i1 = i0 + 129; i1 <= -T + N + i0; i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                 } else {
                   for (register int i1 = -128 * ii0 + i0 + 129; i1 <= -T + N + i0; i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                 }
                 for (register int i1 = max(-T + N + i0 + 1, -128 * ii0 + i0 + 129); i1 <= min(N - 128 * ii0 + i0 - 128, 128 * ii0 - i0 + 382); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
               if (N <= 383 && ii0 >= 1) {
                 for (register int i0 = -N + 64 * ii0 + (T + N + 1) / 2 + 191; i0 <= 128 * ii0 + 126; i0 += 1) {
                   for (register int i1 = -128 * ii0 + i0 + 129; i1 <= 128 * ii0 - i0 + 382; i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                 }
               } else if (ii0 == 0) {
                 for (register int i0 = -N + (T + N + 1) / 2 + 191; i0 <= 126; i0 += 1) {
                   for (register int i1 = i0 + 129; i1 <= -i0 + 382; i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                 }
               }
             }
           }
           if (128 * ii0 + 128 >= T && ii1 == 0) {
             for (register int i0 = max(128 * ii0, -N + 64 * ii0 + (T + N) / 2 + 128); i0 < T - 1; i0 += 1) {
               for (register int i1 = max(0, T - N - i0 + 255); i1 <= min(T - i0 - 2, -128 * ii0 + i0); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
               }
               for (register int i1 = max(max(N + 128 * ii0 - i0 - 1, -T + N + i0 + 1), T - i0 + 255); i1 < N; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           } else if (ii1 == 0) {
             if (N + 128 * ii0 >= T + 129) {
               for (register int i0 = -N + 128 * ii0 + N / 2 + 192; i0 <= -N + 64 * ii0 + (T + N) / 2 + 127; i0 += 1) {
                 for (register int i1 = max(0, -N + 128 * ii0 - i0 + 383); i1 <= -128 * ii0 + i0; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
                 }
                 for (register int i1 = 128 * ii0 - i0 + 383; i1 < N; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
               for (register int i0 = -N + 64 * ii0 + (T + N) / 2 + 128; i0 < 64 * ii0 + (T + 1) / 2 - 1; i0 += 1) {
                 for (register int i1 = max(0, -N + 128 * ii0 - i0 + 383); i1 <= min(-128 * ii0 + i0, 128 * ii0 - i0 + 126); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
                 }
                 for (register int i1 = max(N - 128 * ii0 + i0 - 127, 128 * ii0 - i0 + 383); i1 < N; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             } else {
               for (register int i0 = -N + 128 * ii0 + N / 2 + 192; i0 <= min(128 * ii0 + 126, 64 * ii0 + (T + 1) / 2 - 2); i0 += 1) {
                 for (register int i1 = max(0, -N + 128 * ii0 - i0 + 383); i1 <= min(-128 * ii0 + i0, 128 * ii0 - i0 + 126); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
                 }
                 for (register int i1 = max(N - 128 * ii0 + i0 - 127, 128 * ii0 - i0 + 383); i1 < N; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             }
             for (register int i0 = 64 * ii0 + (T + 1) / 2 - 1; i0 <= min(T - N + 254, 128 * ii0 + 126); i0 += 1) {
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
         if (N + 128 * ii0 >= T + 256) {
           if (N + 128 * ii0 >= T + 257) {
             for (register int i0 = 128 * ii0; i0 < T - 1; i0 += 1) {
               for (register int i1 = max(128 * ii0 - i0 + 255, -T + i0 + 257); i1 <= min(T - i0 + 254, -128 * ii0 + i0 + 256); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           } else {
             for (register int i0 = T - N + 256; i0 <= min(T - 2, T - 2 * N + (5 * N - 58) / 132 + 749); i0 += 1) {
               for (register int i1 = max(-T + i0 + 257, T - N - i0 + 511); i1 <= min(-T + N + i0, T - i0 + 254); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
             for (register int i0 = T - 2 * N + (5 * N - 58) / 132 + 750; i0 < T - 1; i0 += 1) {
               for (register int i1 = -T + i0 + 257; i1 <= T - i0 + 254; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           }
         } else if (128 * ii0 + 128 >= T) {
           for (register int i0 = -N + 64 * ii0 + (T + N) / 2 + 128; i0 <= min(64 * ii0 + (T + N + 1) / 2 - 130, -2 * T + N + 374 * ii0 + (10 * T - 5 * N + 40 * ii0 - 6) / 132); i0 += 1) {
             for (register int i1 = max(128 * ii0 - i0 + 255, -T + i0 + 257); i1 <= min(-T + N + i0, T - i0 + 254); i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
           }
           for (register int i0 = max(-N + 64 * ii0 + (T + N) / 2 + 128, -2 * T + N + 374 * ii0 + floord(10 * T - 5 * N + 40 * ii0 - 6, 132) + 1); i0 < 170 * ii0 + floord(-127 * T + 127 * N + 44 * ii0 + 126, 386); i0 += 1) {
             for (register int i1 = 128 * ii0 - i0 + 255; i1 <= -T + N + i0; i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
           }
           for (register int i0 = max(max(-N + 64 * ii0 + (T + N) / 2 + 128, -2 * T + N + 374 * ii0 + (10 * T - 5 * N + 40 * ii0 - 6) / 132 + 1), 170 * ii0 + floord(-127 * T + 127 * N + 44 * ii0 + 126, 386)); i0 < 64 * ii0 + (T + N + 1) / 2 - 129; i0 += 1) {
             for (register int i1 = max(128 * ii0 - i0 + 255, -T + i0 + 257); i1 <= min(min(N + 128 * ii0 - i0 - 2, -T + N + i0), T - N - 389 * ii0 + 2 * i0 + (-5 * ii0 + 5 * i0 + 5) / 127 + 257); i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
             for (register int i1 = T - N - 389 * ii0 + 2 * i0 + (-5 * ii0 + 5 * i0 + 5) / 127 + 258; i1 <= min(N + 128 * ii0 - i0 - 2, -T + N + i0); i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
             for (register int i1 = N + 128 * ii0 - i0 - 1; i1 <= min(-T + N + i0, T - i0 + 254); i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
           }
           for (register int i0 = 64 * ii0 + (T + N + 1) / 2 - 129; i0 < T - 1; i0 += 1) {
             for (register int i1 = -T + i0 + 257; i1 <= min(-T + N + i0, T - i0 + 254); i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
           }
         }
       }
       if (T >= N + 128 * ii0) {
         #pragma omp parallel for
         for (register int ii1 = 0; ii1 <= 2; ii1 += 1) {
           if (ii1 == 0) {
             for (register int i0 = 128 * ii0 + 64; i0 <= 128 * ii0 + 127; i0 += 1) {
               for (register int i1 = 128 * ii0 - i0 + 127; i1 <= -128 * ii0 + i0; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i0 == 128 * ii0 + 127 && i1 == 0 ? N - 1 : i1 - 1]));
               }
             }
           }
           for (register int i0 = max(128 * ii0 - 64 * ii1 + 128, -N + 128 * ii0 + 63 * ii1 + N / 2 + 130); i0 <= min(128 * ii0 + 190, 128 * ii0 - 63 * ii1 + (N + 1) / 2 + 124); i0 += 1) {
             if (ii1 >= 1) {
               for (register int i1 = max(-128 * ii0 + 128 * ii1 + i0 - 127, 128 * ii0 + 128 * ii1 - i0 + 127); i1 <= min(min(min(N - 128 * ii0 + i0 - 128, -128 * ii0 + 128 * ii1 + i0), N + 128 * ii0 - i0 + 126), 128 * ii0 + 128 * ii1 - i0 + 254); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][ii1 == 2 && i0 == 128 * ii0 + 127 && i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             } else {
               for (register int i1 = -128 * ii0 + i0 - 127; i1 <= 128 * ii0 - i0 + 254; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           }
         }
       } else if (N + 128 * ii0 >= T + 256 && 128 * ii0 + 64 >= T) {
         #pragma omp parallel for
         for (register int ii1 = 0; ii1 <= 2; ii1 += 1) {
           for (register int i0 = 64 * ii0 + T / 2; i0 < T; i0 += 1) {
             if (ii1 == 2) {
               for (register int i1 = T - i0 + 255; i1 <= -128 * ii0 + i0 + 256; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][N + 128 * ii0 == T + 256 && i0 + 1 == T && i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
               for (register int i1 = max(N + 128 * ii0 - i0 - 1, -128 * ii0 + i0 + 257); i1 <= -T + N + i0; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i0 + 1 == T && i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             } else if (ii1 == 1) {
               if (128 * ii0 + 64 == T && i0 + 1 == T) {
                 for (register int i1 = 128; i1 <= -N + 449; i1 += 1) {
                   A[0][i1] = (0.125 * ((A[1][i1 + 1] - (2.0 * A[1][i1])) + A[1][i1 - 1]));
                 }
               }
               for (register int i1 = max(T - i0 + 127, T - N - 389 * ii0 + 2 * i0 + (-5 * ii0 + 5 * i0 + 5) / 127 + 258); i1 <= -128 * ii0 + i0 + 128; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
               for (register int i1 = 128 * ii0 - i0 + 255; i1 <= -T + i0 + 256; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             } else {
               for (register int i1 = T - i0 - 1; i1 <= -128 * ii0 + i0; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i0 + 1 == T && i1 == 0 ? N - 1 : i1 - 1]));
               }
               for (register int i1 = 128 * ii0 - i0 + 127; i1 <= -T + i0 + 128; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           }
         }
       } else if (N + 128 * ii0 >= T + 256 && T >= 128 * ii0 + 65 && 128 * ii0 + 127 >= T) {
         #pragma omp parallel for
         for (register int ii1 = 0; ii1 <= 2; ii1 += 1) {
           for (register int i0 = 64 * ii0 + T / 2; i0 < T; i0 += 1) {
             if (ii1 == 2) {
               for (register int i1 = T - i0 + 255; i1 <= -128 * ii0 + i0 + 256; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][N + 128 * ii0 == T + 256 && i0 + 1 == T && i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
               for (register int i1 = max(N + 128 * ii0 - i0 - 1, -128 * ii0 + i0 + 257); i1 <= -T + N + i0; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i0 + 1 == T && i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             } else if (ii1 == 1) {
               for (register int i1 = T - i0 + 127; i1 <= -128 * ii0 + i0 + 128; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
               for (register int i1 = max(-128 * ii0 + i0 + 129, 128 * ii0 - i0 + 255); i1 <= -T + i0 + 256; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             } else {
               for (register int i1 = T - i0 - 1; i1 <= -128 * ii0 + i0; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i0 + 1 == T && i1 == 0 ? N - 1 : i1 - 1]));
               }
               for (register int i1 = max(-128 * ii0 + i0 + 1, 128 * ii0 - i0 + 127); i1 <= -T + i0 + 128; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           }
         }
       } else if (T >= 128 * ii0 + 128) {
         #pragma omp parallel for
         for (register int ii1 = 0; ii1 <= 2; ii1 += 1) {
           if (ii1 <= 1 && N + 256 * ii0 + 255 * ii1 + 127 >= 2 * T) {
             if (ii1 == 0) {
               for (register int i0 = 128 * ii0 + 64; i0 <= 128 * ii0 + 127; i0 += 1) {
                 for (register int i1 = 128 * ii0 - i0 + 127; i1 <= min(-128 * ii0 + i0, T - N - i0 + 254); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i0 == 128 * ii0 + 127 && i1 == 0 ? N - 1 : i1 - 1]));
                 }
                 for (register int i1 = max(128 * ii0 - i0 + 127, T - N - i0 + 255); i1 <= min(T - i0 - 2, -T + i0 + 128); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i0 == 128 * ii0 + 127 && i1 == 0 ? N - 1 : i1 - 1]));
                 }
                 if (ii0 == 0) {
                   for (register int i1 = T - i0 - 1; i1 <= -T + i0 + 128; i1 += 1) {
                     A[(-i0 + 127) % 2][i1] = (0.125 * ((A[(-i0 + 128) % 2][i1 + 1] - (2.0 * A[(-i0 + 128) % 2][i1])) + A[(-i0 + 128) % 2][i1 - 1]));
                   }
                 } else {
                   for (register int i1 = T - i0 - 1; i1 <= -T + i0 + 128; i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][128 * ii0 + 128 == T && i0 + 1 == T && i1 == 0 ? N - 1 : i1 - 1]));
                   }
                 }
                 for (register int i1 = max(max(128 * ii0 - i0 + 127, -T + i0 + 129), T - N - i0 + 255); i1 <= -128 * ii0 + i0; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             }
             for (register int i0 = 128 * ii0 - 64 * ii1 + 128; i0 <= min(T - 1, 128 * ii0 - 63 * ii1 + 190); i0 += 1) {
               for (register int i1 = max(-128 * ii0 + i0 - 127, 128 * ii0 + 126 * ii1 - i0 + 129); i1 <= min(-128 * ii0 + i0 + 128, 128 * ii0 + 128 * ii1 - i0 + 254); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
             if (ii1 == 1) {
               for (register int i0 = 128 * ii0 + 128; i0 <= min(T - 1, 128 * ii0 + 190); i0 += 1) {
                 for (register int i1 = -128 * ii0 + i0 + 1; i1 <= 128 * ii0 - i0 + 382; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             }
           } else if (T >= 128 * ii0 + 129 && ii1 == 2) {
             for (register int i0 = -N + 128 * ii0 + N / 2 + 256; i0 < min(T, 128 * ii0 + (N + 1) / 2 - 1); i0 += 1) {
               if (i0 >= 128 * ii0 + 128) {
                 for (register int i1 = -128 * ii0 + i0 + 129; i1 <= N + 128 * ii0 - i0 + 126; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               } else {
                 for (register int i1 = 128 * ii0 - i0 + 383; i1 <= min(-T + N + i0, T - i0 + 254); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
                 if (ii0 == 0) {
                   for (register int i1 = T - i0 + 255; i1 <= -T + N + i0; i1 += 1) {
                     A[(-i0 + 127) % 2][i1] = (0.125 * ((A[(-i0 + 128) % 2][i1 + 1] - (2.0 * A[(-i0 + 128) % 2][i1])) + A[(-i0 + 128) % 2][i1 - 1]));
                   }
                 } else {
                   for (register int i1 = T - i0 + 255; i1 <= -T + N + i0; i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                 }
                 for (register int i1 = max(-T + N + i0 + 1, 128 * ii0 - i0 + 383); i1 < N - 128 * ii0 + i0 - 127; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i0 == 128 * ii0 + 127 && i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             }
           } else if (2 * T >= N + 256 * ii0 + 128 && ii1 == 0) {
             for (register int i0 = 128 * ii0 + 64; i0 <= 128 * ii0 + 190; i0 += 1) {
               for (register int i1 = max(-128 * ii0 + i0 - 127, 128 * ii0 - i0 + 127); i1 <= min(-128 * ii0 + i0, 128 * ii0 - i0 + 254); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i0 == 128 * ii0 + 127 && i1 == 0 ? N - 1 : i1 - 1]));
               }
             }
           } else {
             for (register int i0 = T - N + N / 2 + 128; i0 < T; i0 += 1) {
               for (register int i1 = T - i0 + 255; i1 <= -T + N + i0; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i0 + 1 == T && i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           }
         }
       } else {
         #pragma omp parallel for
         for (register int ii1 = 0; ii1 <= 2; ii1 += 1) {
           if (ii1 >= 1) {
             if (ii1 == 1) {
               for (register int i0 = 64 * ii0 + T / 2; i0 < min(64 * ii0 + (T + N) / 2 - 128, -2 * T + N + 374 * ii0 + (10 * T - 5 * N + 40 * ii0 - 1) / 132); i0 += 1) {
                 for (register int i1 = T - i0 + 127; i1 <= min(min(-128 * ii0 + i0 + 128, 128 * ii0 - i0 + 254), T - N - 389 * ii0 + 2 * i0 + (-5 * ii0 + 5 * i0 + 5) / 127 + 257); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
                 for (register int i1 = max(T - i0 + 127, T - N - 389 * ii0 + 2 * i0 + (-5 * ii0 + 5 * i0 + 5) / 127 + 258); i1 <= min(-128 * ii0 + i0 + 128, 128 * ii0 - i0 + 254); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
                 for (register int i1 = 128 * ii0 - i0 + 255; i1 <= -128 * ii0 + i0 + 128; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
                 for (register int i1 = max(-128 * ii0 + i0 + 129, 128 * ii0 - i0 + 255); i1 <= T - N - 389 * ii0 + 2 * i0 + (-5 * ii0 + 5 * i0 + 5) / 127 + 257; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
                 for (register int i1 = max(max(-128 * ii0 + i0 + 129, 128 * ii0 - i0 + 255), T - N - 389 * ii0 + 2 * i0 + (-5 * ii0 + 5 * i0 + 5) / 127 + 258); i1 <= -T + i0 + 256; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
               for (register int i0 = max(64 * ii0 + T / 2, -2 * T + N + 374 * ii0 + (10 * T - 5 * N + 40 * ii0 - 1) / 132); i0 < 64 * ii0 + (T + N) / 2 - 128; i0 += 1) {
                 for (register int i1 = T - i0 + 127; i1 <= min(-128 * ii0 + i0 + 128, 128 * ii0 - i0 + 254); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
                 for (register int i1 = 128 * ii0 - i0 + 255; i1 <= -T + i0 + 256; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             }
             for (register int i0 = max(T - N + 62 * ii1 + (N + ii1) / 2 + 3, 64 * ii0 - 63 * ii1 + (T + N + ii1 + 1) / 2 - 66); i0 < T; i0 += 1) {
               if (ii1 == 2) {
                 for (register int i1 = T - i0 + 255; i1 <= -T + N + i0; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i0 + 1 == T && i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               } else {
                 for (register int i1 = T - i0 + 127; i1 <= min(min(-128 * ii0 + i0 + 128, 128 * ii0 - i0 + 254), T - N - 389 * ii0 + 2 * i0 + (-5 * ii0 + 5 * i0 + 5) / 127 + 257); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
                 for (register int i1 = max(T - i0 + 127, T - N - 389 * ii0 + 2 * i0 + (-5 * ii0 + 5 * i0 + 5) / 127 + 258); i1 <= min(-128 * ii0 + i0 + 128, 128 * ii0 - i0 + 254); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
                 for (register int i1 = 128 * ii0 - i0 + 255; i1 <= -128 * ii0 + i0 + 128; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
                 for (register int i1 = max(-128 * ii0 + i0 + 129, 128 * ii0 - i0 + 255); i1 <= min(N + 128 * ii0 - i0 - 2, T - N - 389 * ii0 + 2 * i0 + (-5 * ii0 + 5 * i0 + 5) / 127 + 257); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
                 for (register int i1 = max(max(-128 * ii0 + i0 + 129, 128 * ii0 - i0 + 255), T - N - 389 * ii0 + 2 * i0 + (-5 * ii0 + 5 * i0 + 5) / 127 + 258); i1 < N + 128 * ii0 - i0 - 1; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
                 for (register int i1 = max(N + 128 * ii0 - i0 - 1, -128 * ii0 + i0 + 129); i1 <= -T + i0 + 256; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             }
           } else {
             for (register int i0 = 64 * ii0 + T / 2; i0 < T; i0 += 1) {
               for (register int i1 = T - i0 - 1; i1 <= -128 * ii0 + i0; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i0 + 1 == T && i1 == 0 ? N - 1 : i1 - 1]));
               }
               for (register int i1 = max(-128 * ii0 + i0 + 1, 128 * ii0 - i0 + 127); i1 <= min(-T + i0 + 128, T - N - 389 * ii0 + 2 * i0 + (-5 * ii0 + 5 * i0 + 5) / 127 + 257); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
               for (register int i1 = max(128 * ii0 - i0 + 127, T - N - 389 * ii0 + 2 * i0 + (-5 * ii0 + 5 * i0 + 5) / 127 + 258); i1 <= -T + i0 + 128; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           }
         }
       }
     }
     if (N >= 259 && (T - 1) % 128 == 0) {
       #pragma omp parallel for
       for (register int ii1 = 0; ii1 <= 2; ii1 += 1) {
         if (ii1 == 0) {
           A[1][0] = (0.125 * ((A[0][1] - (2.0 * A[0][0])) + A[0][N - 1]));
         }
         if (ii1 <= 1) {
           A[1][ii1 + 127] = (0.125 * ((A[0][ii1 + 128] - (2.0 * A[0][ii1 + 127])) + A[0][ii1 + 126]));
         }
         if (ii1 >= 1) {
           A[1][ii1 + 254] = (0.125 * ((A[0][ii1 + 255] - (2.0 * A[0][ii1 + 254])) + A[0][ii1 + 253]));
           if (ii1 == 2) {
             A[1][N - 1] = (0.125 * ((A[0][0] - (2.0 * A[0][N - 1])) + A[0][N - 2]));
           }
         }
       }
     } else if (T >= 128 * floord(T - 27, 128) + 129 && 128 * floord(T - 27, 128) + 32895 >= T + 127 * N) {
       if (N == 257 && (T - 27) % 128 >= 103) {
         if ((T - 27) % 128 >= 104) {
           for (register int i0 = -((T - 27) % 128) + T + 101; i0 <= (T + 1) / 2 + 64 * ((T - 27) / 128) + 62; i0 += 1) {
             for (register int i1 = 0; i1 <= ((T - 1) % 128) - T + i0 + 1; i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? 256 : i1 - 1]));
             }
             for (register int i1 = max(T - i0 + 229, -T + i0 + 258); i1 <= 256; i1 += 1) {
               if (((T - 27) % 128) + i0 + i1 >= T + 357 || (i0 + i1 + 1) % 128 == 0) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 == 256 ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           }
           for (register int i0 = (T + 1) / 2 + 64 * ((T - 1) / 128) - 1; i0 < T - 1; i0 += 1) {
             for (register int i1 = 0; i1 < T - i0 - 2; i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? 256 : i1 - 1]));
             }
             for (register int i1 = -T + i0 + 258; i1 <= 256; i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 == 256 ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
           }
         } else {
           A[1][256] = (0.125 * ((A[0][0] - (2.0 * A[0][256])) + A[0][255]));
         }
       }
       for (register int k = max(-((T - 1) % 128) + 3, N - (T - 1) / 128 + (T - (T - 1) / 128 - 2) / 127 - 255); k <= 3; k += 1) {
         if ((T + 1) % 2 == 0) {
           #pragma omp parallel for
           for (register int ii1 = 0; ii1 < k; ii1 += 1) {
             for (register int ii1_2 = max(-k + 3, -3 * k + ii1 + 7); ii1_2 <= -2 * k + 6; ii1_2 += 1) {
               if (N == 257 && k == 2) {
                 for (register int i0 = max(max(-((T - 1) % 128) + T - 1, ((T - 49) / 2) - 12 * ii1 + 12 * ii1_2 + 64 * ((T - 1) / 128)), ((T - 25) / 2) + 12 * ii1 + 64 * ((T - 1) / 128)); i0 < T - 1; i0 += 1) {
                   if (ii1 == 1 && ii1_2 == 2) {
                     A[(i0 + 1) % 2][-T + i0 + 257] = (0.125 * ((A[i0 % 2][-T + i0 + 258] - (2.0 * A[i0 % 2][-T + i0 + 257])) + A[i0 % 2][-T + i0 + 256]));
                   } else if (ii1_2 == 2) {
                     A[(i0 + 1) % 2][T - i0 - 2] = (0.125 * ((A[i0 % 2][T - i0 - 1] - (2.0 * A[i0 % 2][T - i0 - 2])) + A[i0 % 2][i0 + 2 == T ? 256 : T - i0 - 3]));
                   } else {
                     for (register int i1 = max(-T + i0 + 129, -((T - 1) % 128) + T - i0 + 126); i1 <= min(T - i0 + 126, ((T - 1) % 128) - T + i0 + 129); i1 += 1) {
                       A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                     }
                   }
                 }
               } else {
                 for (register int i0 = max(T + 12 * ii1 - 25, ((T - 1) / 2) + 64 * ((T - 1) / 128)); i0 < T; i0 += 1) {
                   if (ii1 == 1) {
                     for (register int i1 = T - i0 + 127; i1 <= ((T - 1) % 128) - T + i0 + 129; i1 += 1) {
                       A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                     }
                     for (register int i1 = max(T + 127 * N - i0 - 32512, T - i0 + 230); i1 <= -T + i0 + 256; i1 += 1) {
                       if ((i0 + i1 + 1) % 128 == 0 || ((T - 1) % 128) + i0 + i1 >= T + 255) {
                         A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                       }
                     }
                   } else if (ii1 == 0) {
                     for (register int i1 = T - i0 - 1; i1 <= ((T - 1) % 128) - T + i0 + 1; i1 += 1) {
                       A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i0 + 1 == T && i1 == 0 ? N - 1 : i1 - 1]));
                     }
                     for (register int i1 = -((T - 1) % 128) + T - i0 + 126; i1 <= -T + i0 + 128; i1 += 1) {
                       A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                     }
                   } else {
                     for (register int i1 = 256; i1 < N; i1 += 1) {
                       A[1][i1] = (0.125 * ((A[0][i1 + 1 == N ? 0 : 257] - (2.0 * A[0][i1])) + A[0][i1 - 1]));
                     }
                   }
                 }
               }
             }
           }
         } else {
           for (register int i0 = max(-((T - 1) % 128) + T - 1, (T / 2) + 13 * k + 64 * ((T - 1) / 128) - 39); i0 < T + k - 3; i0 += 1) {
             if (k == 3) {
               for (register int i1 = T - i0 - 1; i1 <= ((T - 1) % 128) - T + i0 + 1; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i0 + 1 == T && i1 == 0 ? 256 : i1 - 1]));
               }
               for (register int i1 = -((T - 1) % 128) + T - i0 + 126; i1 <= -T + i0 + 128; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             } else {
               for (register int i1 = max(-T + i0 + 129, -((T - 1) % 128) + T - i0 + 126); i1 <= min(T - i0 + 126, ((T - 1) % 128) - T + i0 + 129); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           }
           if (k == 2) {
             for (register int i0 = (T / 2) + 64 * ((T - 1) / 128) - 1; i0 < T - 1; i0 += 1) {
               A[(i0 + 1) % 2][T - i0 - 2] = (0.125 * ((A[i0 % 2][T - i0 - 1] - (2.0 * A[i0 % 2][T - i0 - 2])) + A[i0 % 2][i0 + 2 == T ? 256 : T - i0 - 3]));
             }
           }
           if (k == 3) {
             #pragma omp parallel for
             for (register int ii1 = 1; ii1 <= 2; ii1 += 1) {
               if (ii1 == 2) {
                 A[0][256] = (0.125 * ((A[1][0] - (2.0 * A[1][256])) + A[1][255]));
               } else {
                 for (register int i0 = (T / 2) + 64 * ((T - 1) / 128); i0 < T; i0 += 1) {
                   for (register int i1 = T - i0 + 127; i1 <= ((T - 1) % 128) - T + i0 + 129; i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                   for (register int i1 = T - i0 + 229; i1 <= -T + i0 + 256; i1 += 1) {
                     if (((T - 1) % 128) + i0 + i1 >= T + 255 || (i0 + i1 + 1) % 128 == 0) {
                       A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                     }
                   }
                 }
               }
             }
           } else {
             for (register int i0 = T - 14; i0 < T - 1; i0 += 1) {
               if (((T - 27) % 128) + 2 * i0 >= 2 * T + 101) {
                 A[(i0 + 1) % 2][-T + i0 + 257] = (0.125 * ((A[i0 % 2][-T + i0 + 258] - (2.0 * A[i0 % 2][-T + i0 + 257])) + A[i0 % 2][-T + i0 + 256]));
               } else if (T % 2 == 0 && (T - 2 * i0 - 2) % 128 == 0) {
                 A[(i0 + 1) % 2][-T + i0 + 257] = (0.125 * ((A[i0 % 2][-T + i0 + 258] - (2.0 * A[i0 % 2][-T + i0 + 257])) + A[i0 % 2][-T + i0 + 256]));
               }
             }
           }
         }
       }
     } else if (N >= 258 && T + 128 * floord(-81727 * T + 16256 * N - 1999042, 10461056) >= 2) {
       for (register int i0 = 128 * ((81727 * T - 16256 * N + 12460097) / 10461056); i0 <= T - N + (N + 1) / 2 + 126; i0 += 1) {
         for (register int i1 = 0; i1 <= min(T - N - i0 + 254, i0 - 128 * ((81727 * T - 16256 * N + 12460097) / 10461056)); i1 += 1) {
           A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
         }
         for (register int i1 = max(-T + N + i0 + 1, -i0 + 128 * ((81727 * T - 16256 * N + 12460097) / 10461056) + 255); i1 <= min(N - i0 + 128 * ((81727 * T - 16256 * N + 12460097) / 10461056) - 2, i0 - 128 * ((81727 * T - 16256 * N + 12460097) / 10461056) + 256); i1 += 1) {
           A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
         }
         for (register int i1 = max(-T + N + i0 + 1, N - i0 + 128 * ((81727 * T - 16256 * N + 12460097) / 10461056) - 1); i1 <= min(N - 1, T - i0 + 254); i1 += 1) {
           A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
         }
       }
       #pragma omp parallel for
       for (register int ii1 = 0; ii1 <= 1; ii1 += 1) {
         for (register int ii1_2 = ii1 + 1; ii1_2 <= 2; ii1_2 += 1) {
           if (ii1_2 == 2) {
             for (register int i0 = max(128 * ((81727 * T - 16256 * N + 12460097) / 10461056), -N + (T + N) / 2 + 64 * ((81727 * T - 16256 * N + 12460097) / 10461056) + 128); i0 < T - 1; i0 += 1) {
               if (ii1 == 1) {
                 for (register int i1 = max(-T + i0 + 257, -i0 + 128 * ((81727 * T - 16256 * N + 12460097) / 10461056) + 255); i1 <= min(min(-T + N + i0, T - i0 + 254), i0 - 128 * ((81727 * T - 16256 * N + 12460097) / 10461056) + 256); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               } else {
                 for (register int i1 = max(0, T - N - i0 + 255); i1 <= min(T - i0 - 2, i0 - 128 * ((81727 * T - 16256 * N + 12460097) / 10461056)); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
                 }
                 for (register int i1 = max(max(-T + N + i0 + 1, T - i0 + 255), -((T - 2) % 128) + T + N - i0 - 3); i1 < N; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             }
           } else {
             for (register int i0 = 128 * ((81727 * T - 16256 * N + 12460097) / 10461056); i0 < T - 1; i0 += 1) {
               for (register int i1 = max(-T + i0 + 129, -i0 + 128 * ((81727 * T - 16256 * N + 12460097) / 10461056) + 127); i1 <= min(T - i0 + 126, i0 - 128 * ((81727 * T - 16256 * N + 12460097) / 10461056) + 128); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           }
         }
       }
       #pragma omp parallel for
       for (register int ii1 = 0; ii1 <= 2; ii1 += 1) {
         if (ii1 <= 1) {
           for (register int i0 = T / 2 + 64 * ((81727 * T - 16256 * N + 12460097) / 10461056); i0 < T; i0 += 1) {
             if (ii1 == 1) {
               for (register int i1 = T - i0 + 127; i1 <= i0 - 128 * ((81727 * T - 16256 * N + 12460097) / 10461056) + 128; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
               for (register int i1 = -i0 + 128 * ((81727 * T - 16256 * N + 12460097) / 10461056) + 255; i1 <= -T + i0 + 256; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             } else {
               for (register int i1 = T - i0 - 1; i1 <= i0 - 128 * ((81727 * T - 16256 * N + 12460097) / 10461056); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i0 + 1 == T && i1 == 0 ? N - 1 : i1 - 1]));
               }
               for (register int i1 = -i0 + 128 * ((81727 * T - 16256 * N + 12460097) / 10461056) + 127; i1 <= -T + i0 + 128; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           }
         } else {
           for (register int i0 = max(T - N + N / 2 + 128, T / 2 + 64 * ((81727 * T - 16256 * N + 12460097) / 10461056)); i0 < T; i0 += 1) {
             for (register int i1 = T - i0 + 255; i1 <= min(-T + N + i0, i0 - 128 * ((81727 * T - 16256 * N + 12460097) / 10461056) + 256); i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i0 + 1 == T && i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
             for (register int i1 = max(N - i0 + 128 * ((81727 * T - 16256 * N + 12460097) / 10461056) - 1, i0 - 128 * ((81727 * T - 16256 * N + 12460097) / 10461056) + 257); i1 <= -T + N + i0; i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i0 + 1 == T && i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
           }
         }
       }
     }
   }
   if (T >= 2 && T <= 64 && N >= 258) {
     for (register int i0 = 0; i0 <= T - N + (N + 1) / 2 + 126; i0 += 1) {
       for (register int i1 = 0; i1 <= T - N - i0 + 254; i1 += 1) {
         A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
       }
       for (register int i1 = -T + N + i0 + 1; i1 <= min(N - 1, T - i0 + 254); i1 += 1) {
         A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
       }
     }
     for (register int k = 2; k <= 3; k += 1) {
       #pragma omp parallel for
       for (register int ii1 = 0; ii1 < k; ii1 += 1) {
         if (k == 2) {
           for (register int ii1_2 = ii1 + 1; ii1_2 <= 2; ii1_2 += 1) {
             for (register int i0 = 0; i0 < T - 1; i0 += 1) {
               if (ii1 == 1 && ii1_2 == 2) {
                 for (register int i1 = -T + i0 + 257; i1 <= min(-T + N + i0, T - i0 + 254); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               } else if (ii1_2 == 2) {
                 for (register int i1 = max(0, T - N - i0 + 255); i1 < T - i0 - 1; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
                 }
                 for (register int i1 = max(-T + N + i0 + 1, T - i0 + 255); i1 < N; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               } else {
                 for (register int i1 = -T + i0 + 129; i1 <= T - i0 + 126; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             }
           }
         } else if (ii1 <= 1) {
           for (register int i0 = 0; i0 < T; i0 += 1) {
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
         } else {
           for (register int i0 = max(0, T - N + N / 2 + 128); i0 < T; i0 += 1) {
             for (register int i1 = T - i0 + 255; i1 <= -T + N + i0; i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i0 + 1 == T && i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
           }
         }
       }
     }
   } else if (T >= 2 && T <= 64 && T + N >= 257 && N <= 257) {
     if (N <= 256) {
       for (register int k = 2; k <= 3; k += 1) {
         #pragma omp parallel for
         for (register int ii1 = 0; ii1 < k - 1; ii1 += 1) {
           if (ii1 == 0) {
             for (register int i0 = 0; i0 < T + k - 3; i0 += 1) {
               if (k == 2) {
                 for (register int i1 = 0; i1 < T - i0 - 1; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
                 }
                 for (register int i1 = -T + i0 + 129; i1 <= min(-T + N + i0, T - i0 + 126); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
                 for (register int i1 = -T + N + i0 + 1; i1 < N; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               } else {
                 for (register int i1 = T - i0 - 1; i1 <= -T + i0 + 128; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i0 + 1 == T && i1 == 0 ? N - 1 : i1 - 1]));
                 }
               }
             }
           } else {
             for (register int i0 = max(0, T - N + N / 2 + 64); i0 < T; i0 += 1) {
               for (register int i1 = T - i0 + 127; i1 <= -T + N + i0; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i0 + 1 == T && i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           }
         }
       }
     } else {
       for (register int k = 1; k <= 3; k += 1) {
         #pragma omp parallel for
         for (register int ii1 = 0; ii1 < k; ii1 += 1) {
           {
             if (k >= ii1 + 2) {
               for (register int i0 = 0; i0 < T + k - 3; i0 += 1) {
                 if (k == 2 && ii1 == 0) {
                   for (register int i1 = -T + i0 + 129; i1 <= T - i0 + 126; i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                 } else if (ii1 == 1) {
                   for (register int i1 = T - i0 + 127; i1 <= -T + i0 + 256; i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                 } else {
                   for (register int i1 = T - i0 - 1; i1 <= -T + i0 + 128; i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i0 + 1 == T && i1 == 0 ? 256 : i1 - 1]));
                   }
                 }
               }
             }
             if (k <= 2) {
               if (k == 2) {
                 for (register int i0 = 0; i0 < T - 1; i0 += 1) {
                   if (ii1 == 1) {
                     A[(i0 + 1) % 2][-T + i0 + 257] = (0.125 * ((A[i0 % 2][-T + i0 + 258] - (2.0 * A[i0 % 2][-T + i0 + 257])) + A[i0 % 2][-T + i0 + 256]));
                   } else {
                     A[(i0 + 1) % 2][T - i0 - 2] = (0.125 * ((A[i0 % 2][T - i0 - 1] - (2.0 * A[i0 % 2][T - i0 - 2])) + A[i0 % 2][i0 + 2 == T ? 256 : T - i0 - 3]));
                   }
                 }
               } else {
                 for (register int i0 = 0; i0 < T - 1; i0 += 1) {
                   for (register int i1 = 0; i1 < T - i0 - 2; i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? 256 : i1 - 1]));
                   }
                   for (register int i1 = -T + i0 + 258; i1 <= 256; i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 == 256 ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                 }
               }
             }
           }
           if (k == 3 && ii1 == 2) {
             A[T % 2][256] = (0.125 * ((A[-(T % 2) + 1][0] - (2.0 * A[-(T % 2) + 1][256])) + A[-(T % 2) + 1][255]));
           }
         }
       }
     }
   }
   if (T == 1 && N >= 256) {
     #pragma omp parallel for
     for (register int ii1 = 0; ii1 <= min(2, N - 255); ii1 += 1) {
       if (ii1 == 2) {
         for (register int i1 = 256; i1 < N; i1 += 1) {
           A[1][i1] = (0.125 * ((A[0][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[0][i1])) + A[0][i1 - 1]));
         }
       } else if (ii1 == 1) {
         for (register int i1 = 128; i1 <= 255; i1 += 1) {
           A[1][i1] = (0.125 * ((A[0][N == 256 && i1 == 255 ? 0 : i1 + 1] - (2.0 * A[0][i1])) + A[0][i1 - 1]));
         }
       } else {
         for (register int i1 = 0; i1 <= 127; i1 += 1) {
           A[1][i1] = (0.125 * ((A[0][i1 + 1] - (2.0 * A[0][i1])) + A[0][i1 == 0 ? N - 1 : i1 - 1]));
         }
       }
     }
   }
 } else {
   if (T >= 129) {
     for (register int ii0 = 0; ii0 <= (T - 2) / 128; ii0 += 1) {
       if (ii0 >= 1) {
         for (register int k = max(max(1, -(N % 128) + 2), ii0 - N / 128 + floord(-T + N + 2 * ii0 - 2 * (N / 128), 126) + 2); k <= 2; k += 1) {
           if (k == 2) {
             for (register int i0 = 128 * ii0; i0 <= min(T - 2, 128 * ii0 + 126); i0 += 1) {
               for (register int i1 = max(max(-128 * ii0 + i0 + 1, 128 * ii0 - i0 + 127), -T + i0 + 129); i1 <= min(min(T - i0 + 126, -128 * ii0 + i0 + 128), 128 * ii0 - i0 + 254); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
             if (T >= 128 * ii0 + 129) {
               for (register int i0 = -N + 128 * ii0 + N / 2 + 64 * ((N - 1) / 128) + 64; i0 <= -((N - 1) % 128) + 128 * ii0 + 126; i0 += 1) {
                 for (register int i1 = -((N - 1) % 128) + 128 * ii0 - i0 + 126; i1 <= min(-128 * ii0 + i0, 128 * ii0 - i0 + 126); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 && (N + i0 + 1) % 128 == 0 ? N - 1 : i1 - 1]));
                 }
               }
               for (register int i0 = -((N - 1) % 128) + 128 * ii0 + 127; i0 <= 128 * ii0 + 126; i0 += 1) {
                 for (register int i1 = 0; i1 <= min(-128 * ii0 + i0, 128 * ii0 - i0 + 126); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
                 }
                 for (register int i1 = max(N - 128 * ii0 + i0 - 127, -((N - 1) % 128) + N + 128 * ii0 - i0 + 126); i1 < N; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             } else {
               for (register int i0 = max(128 * ii0, -N + 64 * ii0 + (T + N) / 2 + 64 * ((N - 1) / 128)); i0 < T - 1; i0 += 1) {
                 for (register int i1 = max(0, -((N - 1) % 128) + T - i0 - 2); i1 <= min(T - i0 - 2, -128 * ii0 + i0); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
                 }
                 for (register int i1 = max(max(N + 128 * ii0 - i0 - 1, -T + N + i0 + 1), -((N - 1) % 128) + T + N - i0 - 2); i1 < N; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             }
             if (N + 128 * ii0 >= T + 384 && 128 * ii0 + 128 >= T) {
               for (register int i0 = 128 * ii0; i0 < T - 1; i0 += 1) {
                 for (register int i1 = max(128 * ii0 - i0 + 255, -T + i0 + 257); i1 <= min(T - i0 + 254, -128 * ii0 + i0 + 256); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             }
             if (128 * ii0 + 128 >= T) {
               #pragma omp parallel for
               for (register int ii1 = 2; ii1 < floord(N - 1, 128) - 1; ii1 += 1) {
                 if (T + 128 * ii1 + 255 >= N + 128 * ii0) {
                   if (ii1 == 2) {
                     for (register int i0 = 128 * ii0; i0 < 128 * ii0 + N / 2 - 256; i0 += 1) {
                       for (register int i1 = 128 * ii0 - i0 + 383; i1 <= -128 * ii0 + i0 + 384; i1 += 1) {
                         A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                       }
                     }
                   }
                   for (register int i0 = max(128 * ii0, 128 * ii0 - 64 * ii1 + 192); i0 < 128 * ii0 - 64 * ii1 + N / 2 - 128; i0 += 1) {
                     for (register int i1 = 128 * ii0 + 128 * ii1 - i0 + 127; i1 <= -128 * ii0 + 128 * ii1 + i0 + 128; i1 += 1) {
                       A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                     }
                   }
                   for (register int i0 = 128 * ii0 - 64 * ii1 + N / 2 - 128; i0 < T - 1; i0 += 1) {
                     if (ii1 == 2) {
                       for (register int i1 = max(128 * ii0 - i0 + 383, -T + i0 + 385); i1 < N + 128 * ii0 - i0 - 129; i1 += 1) {
                         A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                       }
                     }
                     for (register int i1 = max(max(128 * ii0 + 128 * ii1 - i0 + 127, -T + 128 * ii1 + i0 + 129), 128 * ii0 - i0 + 511); i1 < N + 128 * ii0 - i0 - 129; i1 += 1) {
                       A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                     }
                     for (register int i1 = max(N + 128 * ii0 - i0 - 129, -T + 128 * ii1 + i0 + 129); i1 <= min(T + 128 * ii1 - i0 + 126, -128 * ii0 + 128 * ii1 + i0 + 128); i1 += 1) {
                       A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                     }
                   }
                 } else {
                   for (register int i0 = 128 * ii0; i0 < T - 1; i0 += 1) {
                     if (N >= 128 * ii1 + 385) {
                       for (register int i1 = max(max(128 * ii0 + 128 * ii1 - i0 + 127, -T + 128 * ii1 + i0 + 129), 128 * ii0 - i0 + 511); i1 <= min(T + 128 * ii1 - i0 + 126, -128 * ii0 + 128 * ii1 + i0 + 128); i1 += 1) {
                         A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                       }
                     } else {
                       for (register int i1 = max(max(128 * ii0 + 128 * ii1 - i0 + 127, -T + 128 * ii1 + i0 + 129), 128 * ii0 - i0 + 511); i1 <= min(T + 128 * ii1 - i0 + 126, -128 * ii0 + 128 * ii1 + i0 + 128); i1 += 1) {
                         A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                       }
                     }
                     if (ii1 == 2) {
                       for (register int i1 = max(128 * ii0 - i0 + 383, -T + i0 + 385); i1 <= min(T - i0 + 382, -128 * ii0 + i0 + 384); i1 += 1) {
                         A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                       }
                     }
                   }
                 }
               }
               if (T + 383 >= N + 128 * ii0) {
                 for (register int i0 = 128 * ii0; i0 < T - 1; i0 += 1) {
                   for (register int i1 = max(128 * ii0 - i0 + 255, -T + i0 + 257); i1 <= min(T - i0 + 254, -128 * ii0 + i0 + 256); i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                 }
               }
               for (register int i0 = max(128 * ii0, -N + 64 * ii0 + (T + N) / 2 + 64 * ((N - 1) / 128)); i0 < T - 1; i0 += 1) {
                 for (register int i1 = max(-((N - 1) % 128) + N + 128 * ii0 - i0 - 2, -((N - 1) % 128) - T + N + i0); i1 <= min(min(-T + N + i0, -((N - 1) % 128) + T + N - i0 - 3), -((N - 1) % 128) + N - 128 * ii0 + i0 - 1); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             } else {
               #pragma omp parallel for
               for (register int ii1 = 1; ii1 < N / 128 - 2; ii1 += 1) {
                 for (register int i0 = 128 * ii0; i0 <= 128 * ii0 + 126; i0 += 1) {
                   for (register int i1 = max(max(-128 * ii0 + 128 * ii1 + i0 + 1, 128 * ii0 + 128 * ii1 - i0 + 127), 128 * ii0 - i0 + 511); i1 <= min(-128 * ii0 + 128 * ii1 + i0 + 128, 128 * ii0 + 128 * ii1 - i0 + 254); i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                   if (ii1 <= 2) {
                     for (register int i1 = max(-128 * ii0 + 128 * ii1 + i0 + 1, 128 * ii0 + 128 * ii1 - i0 + 127); i1 <= min(-128 * ii0 + 128 * ii1 + i0 + 128, 128 * ii0 + 128 * ii1 - i0 + 254); i1 += 1) {
                       A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                     }
                   }
                 }
               }
               if (N % 128 >= 1) {
                 for (register int i0 = 128 * ii0; i0 <= 128 * ii0 + 126; i0 += 1) {
                   if (N <= 640) {
                     for (register int i1 = max(-((N - 1) % 128) + N - 128 * ii0 + i0 - 256, -((N - 1) % 128) + N + 128 * ii0 - i0 - 130); i1 < min(N + 128 * ii0 - i0 - 129, -((N - 1) % 128) + N - 128 * ii0 + i0 - 128); i1 += 1) {
                       A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                     }
                   }
                   for (register int i1 = max(max(128 * ii0 - i0 + 511, -((N + 127) % 128) + N - 128 * ii0 + i0 - 256), -((N + 127) % 128) + N + 128 * ii0 - i0 - 130); i1 < min(N + 128 * ii0 - i0 - 129, -((N + 127) % 128) + N - 128 * ii0 + i0 - 128); i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                   for (register int i1 = max(N + 128 * ii0 - i0 - 129, -((N - 1) % 128) + N - 128 * ii0 + i0 - 256); i1 < min(-((N - 1) % 128) + N - 128 * ii0 + i0 - 128, -((N - 1) % 128) + N + 128 * ii0 - i0 - 2); i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                 }
               }
               for (register int i0 = -N + 128 * ii0 + N / 2 + 64 * ((N - 1) / 128) + 64; i0 <= 128 * ii0 + 126; i0 += 1) {
                 for (register int i1 = max(-((N - 1) % 128) + N - 128 * ii0 + i0 - 128, -((N - 1) % 128) + N + 128 * ii0 - i0 - 2); i1 <= min(N - 128 * ii0 + i0 - 128, -((N - 1) % 128) + N + 128 * ii0 - i0 + 125); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             }
           } else if (T >= 128 * ii0 + 3 && 128 * ii0 + 253 >= T) {
             for (register int i0 = 128 * ii0; i0 <= min(T - N + (N + 1) / 2 + 64 * ((N - 1) / 128) - 2, -N + 128 * ii0 + (N + 1) / 2 + 64 * ((N - 1) / 128) + 126); i0 += 1) {
               for (register int i1 = 0; i1 <= min(min(-128 * ii0 + i0, -((N - 1) % 128) + T - i0 - 3), -((N - 1) % 128) + 128 * ii0 - i0 + 125); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
               }
               if (128 * ii0 + 128 >= T) {
                 for (register int i1 = max(-T + N + i0 + 1, -((N - 1) % 128) + N + 128 * ii0 - i0 - 2); i1 <= min(N + 128 * ii0 - i0 - 2, -((N - 1) % 128) + N - 128 * ii0 + i0 - 1); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               } else {
                 for (register int i1 = max(N - 128 * ii0 + i0 - 127, -((N - 1) % 128) + N + 128 * ii0 - i0 - 2); i1 <= min(N + 128 * ii0 - i0 - 2, -((N - 1) % 128) + N - 128 * ii0 + i0 - 1); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
               for (register int i1 = max(max(N - 128 * ii0 + i0 - 127, N + 128 * ii0 - i0 - 1), -T + N + i0 + 1); i1 <= min(min(N - 1, -((N - 1) % 128) + T + N - i0 - 3), -((N - 1) % 128) + N + 128 * ii0 - i0 + 125); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           } else if (T >= 128 * ii0 + 254 && (N - 127) % 128 >= 2) {
             for (register int i0 = 128 * ii0; i0 <= -N + 128 * ii0 + (N + 1) / 2 + 64 * ((N - 1) / 128) + 126; i0 += 1) {
               for (register int i1 = 0; i1 <= min(-128 * ii0 + i0, -((N - 1) % 128) + 128 * ii0 - i0 + 125); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
               }
               for (register int i1 = max(N - 128 * ii0 + i0 - 127, -((N - 1) % 128) + N + 128 * ii0 - i0 - 2); i1 <= min(N + 128 * ii0 - i0 - 2, -((N - 1) % 128) + N - 128 * ii0 + i0 - 1); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
               for (register int i1 = max(N - 128 * ii0 + i0 - 127, N + 128 * ii0 - i0 - 1); i1 <= min(N - 1, -((N - 1) % 128) + N + 128 * ii0 - i0 + 125); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           } else if (T >= 128 * ii0 + 254 && (N + 1) % 128 == 0) {
             for (register int i0 = 128 * ii0; i0 <= 128 * ii0 + 63; i0 += 1) {
               if (128 * ii0 + 62 >= i0) {
                 A[(i0 + 1) % 2][N - 128 * ii0 + i0 - 127] = (0.125 * ((A[i0 % 2][N - 128 * ii0 + i0 - 126] - (2.0 * A[i0 % 2][N - 128 * ii0 + i0 - 127])) + A[i0 % 2][N - 128 * ii0 + i0 - 128]));
               }
               A[(i0 + 1) % 2][N + 128 * ii0 - i0 - 1] = (0.125 * ((A[i0 % 2][i0 == 128 * ii0 ? 0 : N + 128 * ii0 - i0] - (2.0 * A[i0 % 2][N + 128 * ii0 - i0 - 1])) + A[i0 % 2][N + 128 * ii0 - i0 - 2]));
             }
           } else {
             A[1][N - 1] = (0.125 * ((A[0][0] - (2.0 * A[0][N - 1])) + A[0][N - 2]));
           }
         }
         if (128 * ii0 + 128 >= T) {
           for (register int i0 = 64 * ii0 + T / 2; i0 < T; i0 += 1) {
             for (register int i1 = T - i0 - 1; i1 <= -128 * ii0 + i0; i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i0 + 1 == T && i1 == 0 ? N - 1 : i1 - 1]));
             }
             for (register int i1 = max(-128 * ii0 + i0 + 1, 128 * ii0 - i0 + 127); i1 <= -T + i0 + 128; i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
           }
         }
         #pragma omp parallel for
         for (register int ii1 = max(0, ii0 - (T - ii0 + 125) / 127 + 2); ii1 <= min(2, (N - 1) / 128 - 2); ii1 += 1) {
           if (T + 125 >= 128 * ii0 + 127 * ii1) {
             if (T >= 128 * ii0 + 129) {
               for (register int i0 = 128 * ii0 + 64; i0 <= min(min(T - 1, 128 * ii0 + 190), 128 * ii0 - 126 * ii1 + 316); i0 += 1) {
                 if (ii1 == 1 && i0 >= 128 * ii0 + 128) {
                   for (register int i1 = -128 * ii0 + i0 + 1; i1 <= 128 * ii0 - i0 + 382; i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                 } else if (ii1 == 0) {
                   for (register int i1 = max(-128 * ii0 + i0 - 127, 128 * ii0 - i0 + 127); i1 <= min(-128 * ii0 + i0, 128 * ii0 - i0 + 254); i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i0 == 128 * ii0 + 127 && i1 == 0 ? N - 1 : i1 - 1]));
                   }
                 } else if (ii1 == 1 && i0 >= 128 * ii0 + 65) {
                   for (register int i1 = 128 * ii0 - i0 + 255; i1 <= -128 * ii0 + i0 + 128; i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                 } else if (ii1 == 2) {
                   for (register int i1 = 319; i1 <= 320; i1 += 1) {
                     A[1][i1] = (0.125 * ((A[0][i1 + 1] - (2.0 * A[0][i1])) + A[0][i1 - 1]));
                   }
                 } else {
                   for (register int i1 = 191; i1 <= 192; i1 += 1) {
                     A[1][i1] = (0.125 * ((A[0][i1 + 1] - (2.0 * A[0][i1])) + A[0][i1 - 1]));
                   }
                 }
               }
               if (ii1 == 2) {
                 for (register int i0 = 128 * ii0 + 65; i0 <= min(T - 1, 128 * ii0 + 190); i0 += 1) {
                   for (register int i1 = max(-128 * ii0 + i0 + 129, 128 * ii0 - i0 + 383); i1 <= min(-128 * ii0 + i0 + 256, 128 * ii0 - i0 + 510); i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                 }
               }
             } else {
               for (register int i0 = 64 * ii0 + T / 2; i0 < T; i0 += 1) {
                 if (i0 + 63 >= T) {
                   for (register int i1 = T - i0 + 127; i1 <= min(N + 128 * ii0 - i0 - 130, -128 * ii0 + i0 + 128); i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                 } else {
                   for (register int i1 = 191; i1 <= T - 128 * ii0 + 64; i1 += 1) {
                     A[T - 128 * ii0 - 127][i1] = (0.125 * ((A[-T + 128 * ii0 + 128][i1 + 1] - (2.0 * A[-T + 128 * ii0 + 128][i1])) + A[-T + 128 * ii0 + 128][i1 - 1]));
                   }
                 }
                 for (register int i1 = max(-128 * ii0 + i0 + 129, 128 * ii0 - i0 + 255); i1 <= min(N + 128 * ii0 - i0 - 130, -T + i0 + 256); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
                 if (i0 + 63 >= T) {
                   for (register int i1 = N + 128 * ii0 - i0 - 129; i1 <= -T + i0 + 256; i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                 }
               }
             }
           } else {
             for (register int i0 = 64 * ii0 + T / 2; i0 < T; i0 += 1) {
               for (register int i1 = T - i0 + 255; i1 <= min(N + 128 * ii0 - i0 - 130, -128 * ii0 + i0 + 256); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
               for (register int i1 = max(-128 * ii0 + i0 + 257, 128 * ii0 - i0 + 383); i1 <= min(N + 128 * ii0 - i0 - 130, -T + i0 + 384); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
               for (register int i1 = N + 128 * ii0 - i0 - 129; i1 <= -T + i0 + 384; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           }
         }
         if (T >= 128 * ii0 + 129) {
           #pragma omp parallel for
           for (register int ii1 = 3; ii1 < floord(N - 1, 128) - 1; ii1 += 1) {
             for (register int i0 = 128 * ii0 + 64; i0 <= min(T - 1, 128 * ii0 + 190); i0 += 1) {
               if (i0 >= 128 * ii0 + 65) {
                 for (register int i1 = max(-128 * ii0 + 128 * ii1 + i0 - 127, 128 * ii0 + 128 * ii1 - i0 + 127); i1 <= min(-128 * ii0 + 128 * ii1 + i0, 128 * ii0 + 128 * ii1 - i0 + 254); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               } else {
                 for (register int i1 = 128 * ii1 + 63; i1 <= 128 * ii1 + 64; i1 += 1) {
                   A[1][i1] = (0.125 * ((A[0][i1 + 1] - (2.0 * A[0][i1])) + A[0][i1 - 1]));
                 }
               }
             }
           }
           if (((N - 1) % 128) + 128 * ii0 + 128 >= T) {
             for (register int i0 = 128 * ii0 + 64; i0 <= min(T - 1, 128 * ii0 + 190); i0 += 1) {
               if (i0 >= 128 * ii0 + 128) {
                 for (register int i1 = -((N - 1) % 128) + N - 128 * ii0 + i0 - 256; i1 < -((N - 1) % 128) + N + 128 * ii0 - i0 + 126; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               } else {
                 for (register int i1 = -((N - 1) % 128) + N + 128 * ii0 - i0 - 2; i1 < -((N - 1) % 128) + N - 128 * ii0 + i0 - 128; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             }
           }
         } else {
           if (N >= 769 && 128 * ii0 + 127 >= T) {
             for (register int i0 = 64 * ii0 + T / 2; i0 < T; i0 += 1) {
               for (register int i1 = T - i0 + 383; i1 <= -128 * ii0 + i0 + 384; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
               for (register int i1 = max(-128 * ii0 + i0 + 385, 128 * ii0 - i0 + 511); i1 <= -T + i0 + 512; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           } else if (N >= 641 && 128 * ii0 + 128 == T) {
             for (register int i0 = T - 64; i0 < T; i0 += 1) {
               for (register int i1 = T - i0 + 383; i1 <= -T + i0 + 512; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           }
           #pragma omp parallel for
           for (register int ii1 = 4; ii1 < floord(N - 1, 128) - 2; ii1 += 1) {
             for (register int i0 = max(T - 64, 128 * ii0); i0 < T; i0 += 1) {
               for (register int i1 = T + 128 * ii1 - i0 - 1; i1 <= -T + 128 * ii1 + i0 + 128; i1 += 1) {
                 if (2 * i0 + 1 >= ((i0 + i1 + 129) % 128) + 256 * ii0) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             }
           }
           if (T + 384 >= N + 128 * ii0) {
             for (register int i0 = max(T - 64, 128 * ii0); i0 < T; i0 += 1) {
               for (register int i1 = T - i0 + 255; i1 <= -T + i0 + 384; i1 += 1) {
                 if (i0 + i1 + 1 >= N + 128 * ii0 || 2 * i0 + 1 >= ((i0 + i1 + 129) % 128) + 256 * ii0) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             }
           } else {
             if (N + 128 * ii0 >= ((N - 1) % 128) + T + 514) {
               if (128 * ii0 + 127 >= T) {
                 for (register int i0 = 64 * ii0 + T / 2; i0 < T; i0 += 1) {
                   for (register int i1 = max(128 * ii0 - i0 + 511, -((N + 127) % 128) + T + N - i0 - 258); i1 < min(-((N + 127) % 128) + N - 128 * ii0 + i0 - 256, -((N + 127) % 128) + N + 128 * ii0 - i0 - 130); i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                   if (N <= 768) {
                     for (register int i1 = T - i0 + 383; i1 <= min(N + 128 * ii0 - i0 - 130, -128 * ii0 + i0 + 384); i1 += 1) {
                       A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                     }
                   }
                   for (register int i1 = max(-128 * ii0 + i0 + 385, -((N + 127) % 128) + N + 128 * ii0 - i0 - 130); i1 < min(N + 128 * ii0 - i0 - 129, -((N + 127) % 128) - T + N + i0 - 128); i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                   for (register int i1 = N + 128 * ii0 - i0 - 129; i1 < -((N - 1) % 128) - T + N + i0 - 128; i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                 }
               } else {
                 for (register int i0 = T - 64; i0 < T; i0 += 1) {
                   for (register int i1 = -((N - 1) % 128) + T + N - i0 - 258; i1 < -((N - 1) % 128) - T + N + i0 - 128; i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                 }
               }
             }
             if (128 * ii0 + 512 >= T + N || ((N - 1) % 128) + 128 * ii0 >= T) {
               for (register int i0 = max(128 * ii0, 64 * ii0 + (T + N) / 2 - 64 * ((N + 127) / 128)); i0 < T; i0 += 1) {
                 for (register int i1 = -((N + 127) % 128) + T + N - i0 - 130; i1 < min(N + 128 * ii0 - i0 - 129, -((N + 127) % 128) + N - 128 * ii0 + i0 - 128); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
                 for (register int i1 = N + 128 * ii0 - i0 - 129; i1 < -((N + 127) % 128) - T + N + i0; i1 += 1) {
                   if (2 * i0 + 1 >= ((i0 + i1 + 129) % 128) + 256 * ii0) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                 }
               }
             }
             if (256 * ii0 + 128 * floord(-T + N - 1, 128) >= T + 257) {
               if (N >= 128 * ii0 + 128 * floord(-T + N - 1, 128) + 129) {
                 for (register int i0 = max(T - 64, 128 * ii0); i0 < (T + N) / 2 - 64 * floord(-T + N - 1, 128) - 64; i0 += 1) {
                   for (register int i1 = T + 128 * ii0 - i0 + 128 * floord(-T + N - 1, 128) - 1; i1 <= -T + 128 * ii0 + i0 + 128 * floord(-T + N - 1, 128) + 128; i1 += 1) {
                     if (2 * i0 + 1 >= ((i0 + i1 + 129) % 128) + 256 * ii0) {
                       A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                     }
                   }
                 }
               }
               for (register int i0 = max(64 * ii0 + T / 2, (T + N) / 2 - 64 * floord(-T + N - 1, 128) - 64); i0 < T; i0 += 1) {
                 if (N >= 128 * ii0 + 128 * floord(-T + N - 1, 128) + 129) {
                   for (register int i1 = T + 128 * ii0 - i0 + 128 * floord(-T + N - 1, 128) - 1; i1 < N + 128 * ii0 - i0 - 1; i1 += 1) {
                     if (2 * i0 + 1 >= ((i0 + i1 + 129) % 128) + 256 * ii0) {
                       A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                     }
                   }
                 } else {
                   for (register int i1 = T + 128 * ii0 - i0 + 128 * floord(-T + N - 1, 128) - 1; i1 <= min(N + 128 * ii0 - i0 - 2, i0 + 128 * floord(-T + N - 1, 128)); i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                 }
                 for (register int i1 = N + 128 * ii0 - i0 - 1; i1 <= min(-T + N + i0, -T + 128 * ii0 + i0 + 128 * floord(-T + N - 1, 128) + 128); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i0 + 1 == T && i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             } else {
               for (register int i0 = T - 64; i0 < T; i0 += 1) {
                 for (register int i1 = T - i0 + 383; i1 <= -T + i0 + 512; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             }
           }
         }
         #pragma omp parallel for
         for (register int ii1 = max(ii0 + floord(-T + N - 1, 128) + 1, (N - 1) / 128 - 1); ii1 <= (N - 1) / 128; ii1 += 1) {
           if (T >= 128 * ii0 + 129) {
             for (register int i0 = max(128 * ii0 + 64, -N + 128 * ii0 + 64 * ii1 + N / 2 + 128); i0 <= min(min(T - 1, 128 * ii0 + 190), 128 * ii0 - 64 * ii1 + (N + 1) / 2 + 126); i0 += 1) {
               if (i0 >= 128 * ii0 + 128) {
                 for (register int i1 = -128 * ii0 + 128 * ii1 + i0 - 127; i1 <= min(N + 128 * ii0 - i0 + 126, 128 * ii0 + 128 * ii1 - i0 + 254); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               } else {
                 for (register int i1 = 128 * ii0 + 128 * ii1 - i0 + 127; i1 <= min(N - 128 * ii0 + i0 - 128, -128 * ii0 + 128 * ii1 + i0); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i0 == 128 * ii0 + 127 && i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             }
           } else {
             for (register int i0 = T - N + 64 * ii1 + N / 2; i0 < T; i0 += 1) {
               for (register int i1 = T + 128 * ii1 - i0 - 1; i1 <= -T + N + i0; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i0 + 1 == T && i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           }
         }
       } else {
         for (register int k = max(0, 43 * ((N + 2) / 128) - (N + (N + 2) / 128 + 2) / 3 + 1); k <= 2; k += 1) {
           if (k >= 1) {
             #pragma omp parallel for
             for (register int ii1 = 0; ii1 < k + (N - 1) / 128 - 2; ii1 += 1) {
               if (k == 2) {
                 for (register int i0 = max(0, -N + 64 * ii1 + N / 2 + 128); i0 <= 126; i0 += 1) {
                   if (N >= 128 * ii1 + 257) {
                     for (register int i1 = max(max(128 * ii1 + i0 + 1, 128 * ii1 - i0 + 127), i0 + 129); i1 <= min(128 * ii1 + i0 + 128, 128 * ii1 - i0 + 254); i1 += 1) {
                       A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                     }
                     if (ii1 == 0) {
                       for (register int i1 = max(i0 + 1, -i0 + 127); i1 <= min(i0 + 128, -i0 + 254); i1 += 1) {
                         A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                       }
                     }
                   } else {
                     for (register int i1 = max(128 * ii1 + i0 + 1, 128 * ii1 - i0 + 127); i1 <= min(N + i0 - 128, 128 * ii1 - i0 + 254); i1 += 1) {
                       A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                     }
                   }
                 }
               } else {
                 for (register int i0 = 0; i0 <= 62; i0 += 1) {
                   for (register int i1 = 128 * ii1 + i0 + 129; i1 <= min(N + i0 - 128, 128 * ii1 - i0 + 254); i1 += 1) {
                     A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                 }
                 if (ii1 == 0) {
                   for (register int i0 = 0; i0 <= 62; i0 += 1) {
                     for (register int i1 = max(i0 + 1, -((N - 1) % 128) - i0 + 126); i1 <= -i0 + 126; i1 += 1) {
                       A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                     }
                   }
                 }
               }
               if (ii1 == 0) {
                 if (k == 2) {
                   for (register int i0 = -N + N / 2 + 64 * ((N - 1) / 128) + 64; i0 <= -((N - 1) % 128) + 126; i0 += 1) {
                     for (register int i1 = -((N - 1) % 128) - i0 + 126; i1 <= min(i0, -i0 + 126); i1 += 1) {
                       A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 && (N + i0 + 1) % 128 == 0 ? N - 1 : i1 - 1]));
                     }
                   }
                   for (register int i0 = -((N - 1) % 128) + 127; i0 <= 126; i0 += 1) {
                     for (register int i1 = 0; i1 <= min(i0, -i0 + 126); i1 += 1) {
                       A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
                     }
                     for (register int i1 = max(N + i0 - 127, -((N - 1) % 128) + N - i0 + 126); i1 < N; i1 += 1) {
                       A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                     }
                   }
                 } else {
                   for (register int i0 = 0; i0 <= -N + (N + 1) / 2 + 64 * ((N - 1) / 128) + 126; i0 += 1) {
                     for (register int i1 = 0; i1 <= min(i0, -((N - 1) % 128) - i0 + 125); i1 += 1) {
                       A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
                     }
                     for (register int i1 = max(N + i0 - 127, -((N - 1) % 128) + N - i0 - 2); i1 <= min(N - 1, -((N - 1) % 128) + N - i0 + 125); i1 += 1) {
                       A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                     }
                   }
                 }
               }
             }
           } else {
             for (register int i0 = 0; i0 <= -N + (N + 1) / 2 + 64 * ((N - 1) / 128) + 62; i0 += 1) {
               for (register int i1 = i0 + 1; i1 <= -((N - 1) % 128) - i0 + 125; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
               for (register int i1 = N + i0 - 127; i1 <= -i0 + 382; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
             if (N >= 513) {
               for (register int i0 = 0; i0 <= -N + (N + 1) / 2 + 64 * ((N - 1) / 128) + 62; i0 += 1) {
                 for (register int i1 = N + i0 - 127; i1 < -((N - 1) % 128) + N - i0 - 2; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             }
           }
         }
         #pragma omp parallel for
         for (register int ii1 = 0; ii1 <= (N - 1) / 128; ii1 += 1) {
           if (ii1 >= 2) {
             for (register int i0 = max(64, -N + 64 * ii1 + N / 2 + 128); i0 <= 127; i0 += 1) {
               if (128 * ii1 + 256 >= N) {
                 for (register int i1 = 128 * ii1 - i0 + 127; i1 <= min(N + i0 - 128, 128 * ii1 + i0); i1 += 1) {
                   A[(-i0 + 127) % 2][i1] = (0.125 * ((A[(-i0 + 128) % 2][i0 == 127 && i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[(-i0 + 128) % 2][i1])) + A[(-i0 + 128) % 2][i1 - 1]));
                 }
               } else if (i0 >= 65) {
                 for (register int i1 = 128 * ii1 - i0 + 127; i1 <= 128 * ii1 + i0; i1 += 1) {
                   A[-(i0 % 2) + 1][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               } else {
                 for (register int i1 = 128 * ii1 + 63; i1 <= 128 * ii1 + 64; i1 += 1) {
                   A[1][i1] = (0.125 * ((A[0][i1 + 1] - (2.0 * A[0][i1])) + A[0][i1 - 1]));
                 }
               }
             }
             if (N >= 128 * ii1 + 257) {
               for (register int i0 = 128; i0 <= min(190, T - 1); i0 += 1) {
                 if (ii1 >= 3) {
                   for (register int i1 = 128 * ii1 + i0 - 127; i1 <= 128 * ii1 - i0 + 254; i1 += 1) {
                     A[-(i0 % 2) + 1][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                 } else {
                   for (register int i1 = i0 + 129; i1 <= -i0 + 510; i1 += 1) {
                     A[-(i0 % 2) + 1][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                   }
                 }
               }
             } else {
               for (register int i0 = 128; i0 <= min(min(190, T - 1), -64 * ii1 + (N + 1) / 2 + 126); i0 += 1) {
                 for (register int i1 = 128 * ii1 + i0 - 127; i1 <= min(N - i0 + 126, 128 * ii1 - i0 + 254); i1 += 1) {
                   A[-(i0 % 2) + 1][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             }
           } else {
             for (register int i0 = 64; i0 <= min(190, T - 1); i0 += 1) {
               if (ii1 == 1 && i0 >= 128) {
                 for (register int i1 = i0 + 1; i1 <= -i0 + 382; i1 += 1) {
                   A[-(i0 % 2) + 1][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               } else if (ii1 == 0) {
                 for (register int i1 = max(i0 - 127, -i0 + 127); i1 <= min(i0, -i0 + 254); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i0 == 127 && i1 == 0 ? N - 1 : i1 - 1]));
                 }
               } else if (i0 >= 65) {
                 for (register int i1 = -i0 + 255; i1 <= i0 + 128; i1 += 1) {
                   A[-(i0 % 2) + 1][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               } else {
                 for (register int i1 = 191; i1 <= 192; i1 += 1) {
                   A[1][i1] = (0.125 * ((A[0][i1 + 1] - (2.0 * A[0][i1])) + A[0][i1 - 1]));
                 }
               }
             }
           }
         }
       }
     }
     if (T >= 257 && (T - 1) % 128 == 0) {
       #pragma omp parallel for
       for (register int ii1 = 0; ii1 <= min(2, (N - 1) / 128 - 2); ii1 += 1) {
         if (ii1 <= 1) {
           if (ii1 == 0) {
             A[1][0] = (0.125 * ((A[0][1] - (2.0 * A[0][0])) + A[0][N - 1]));
           }
           A[1][ii1 + 127] = (0.125 * ((A[0][ii1 + 128] - (2.0 * A[0][ii1 + 127])) + A[0][ii1 + 126]));
           if (ii1 == 1) {
             A[1][255] = (0.125 * ((A[0][256] - (2.0 * A[0][255])) + A[0][254]));
           }
         } else {
           A[1][256] = (0.125 * ((A[0][257] - (2.0 * A[0][256])) + A[0][255]));
           A[1][383] = (0.125 * ((A[0][384] - (2.0 * A[0][383])) + A[0][382]));
         }
       }
       #pragma omp parallel for
       for (register int ii1 = 3; ii1 < floord(N - 1, 128) - 2; ii1 += 1) {
         if (ii1 == 3) {
           A[1][384] = (0.125 * ((A[0][385] - (2.0 * A[0][384])) + A[0][383]));
         }
         for (register int i1 = max(128 * ii1, ((ii1 - 3) % 127) + 511); i1 <= 128 * ii1 + 127; i1 += 127) {
           A[1][i1] = (0.125 * ((A[0][i1 + 1] - (2.0 * A[0][i1])) + A[0][i1 - 1]));
         }
       }
       if (N >= 641) {
         if (N <= 768) {
           A[1][384] = (0.125 * ((A[0][385] - (2.0 * A[0][384])) + A[0][383]));
         }
         for (register int i1 = max(511, N - 384); i1 < N - 129; i1 += 1) {
           if ((i1 + 257 >= N && (i1 + 1) % 128 == 0) || (N >= i1 + 257 && i1 % 128 == 0)) {
             A[1][i1] = (0.125 * ((A[0][i1 + 1] - (2.0 * A[0][i1])) + A[0][i1 - 1]));
           }
         }
       }
       #pragma omp parallel for
       for (register int ii1 = (N - 1) / 128 - 1; ii1 <= (N - 1) / 128; ii1 += 1) {
         if (N >= 128 * ii1 + 129) {
           if (N >= 128 * ii1 + 130) {
             A[1][128 * ii1] = (0.125 * ((A[0][128 * ii1 + 1] - (2.0 * A[0][128 * ii1])) + A[0][128 * ii1 - 1]));
           }
           for (register int i1 = -((N - 128 * ii1 - 3) % 127) + N - 3; i1 <= 128 * ii1 + 127; i1 += 127) {
             A[1][i1] = (0.125 * ((A[0][i1 + 1] - (2.0 * A[0][i1])) + A[0][i1 - 1]));
           }
         } else {
           if (N >= 128 * ii1 + 2) {
             A[1][128 * ii1] = (0.125 * ((A[0][128 * ii1 + 1] - (2.0 * A[0][128 * ii1])) + A[0][128 * ii1 - 1]));
           }
           A[1][N - 1] = (0.125 * ((A[0][0] - (2.0 * A[0][N - 1])) + A[0][N - 2]));
         }
       }
     } else if (T == 129) {
       #pragma omp parallel for
       for (register int ii1 = 0; ii1 <= min(2, floord(N - 1, 128) - 3); ii1 += 1) {
         if (ii1 == 2) {
           A[1][256] = (0.125 * ((A[0][257] - (2.0 * A[0][256])) + A[0][255]));
           A[1][383] = (0.125 * ((A[0][384] - (2.0 * A[0][383])) + A[0][382]));
         } else if (ii1 == 1) {
           A[1][128] = (0.125 * ((A[0][129] - (2.0 * A[0][128])) + A[0][127]));
           A[1][255] = (0.125 * ((A[0][256] - (2.0 * A[0][255])) + A[0][254]));
         } else {
           A[1][0] = (0.125 * ((A[0][1] - (2.0 * A[0][0])) + A[0][N - 1]));
           A[1][127] = (0.125 * ((A[0][128] - (2.0 * A[0][127])) + A[0][126]));
         }
       }
       if (N <= 640) {
         if (N >= 513) {
           A[1][256] = (0.125 * ((A[0][257] - (2.0 * A[0][256])) + A[0][255]));
           A[1][383] = (0.125 * ((A[0][384] - (2.0 * A[0][383])) + A[0][382]));
         } else {
           A[1][128] = (0.125 * ((A[0][129] - (2.0 * A[0][128])) + A[0][127]));
           A[1][255] = (0.125 * ((A[0][256] - (2.0 * A[0][255])) + A[0][254]));
         }
       }
       #pragma omp parallel for
       for (register int ii1 = 3; ii1 < floord(N - 1, 128) - 2; ii1 += 1) {
         if (ii1 == 3) {
           A[1][384] = (0.125 * ((A[0][385] - (2.0 * A[0][384])) + A[0][383]));
         }
         for (register int i1 = max(128 * ii1, ((ii1 - 3) % 127) + 511); i1 <= 128 * ii1 + 127; i1 += 127) {
           A[1][i1] = (0.125 * ((A[0][i1 + 1] - (2.0 * A[0][i1])) + A[0][i1 - 1]));
         }
       }
       if (N >= 641) {
         if (N <= 768) {
           A[1][384] = (0.125 * ((A[0][385] - (2.0 * A[0][384])) + A[0][383]));
         }
         for (register int i1 = max(511, N - 384); i1 < N - 129; i1 += 1) {
           if ((N >= i1 + 257 && i1 % 128 == 0) || (i1 + 257 >= N && (i1 + 1) % 128 == 0)) {
             A[1][i1] = (0.125 * ((A[0][i1 + 1] - (2.0 * A[0][i1])) + A[0][i1 - 1]));
           }
         }
       }
       #pragma omp parallel for
       for (register int ii1 = (N - 1) / 128 - 1; ii1 <= (N - 1) / 128; ii1 += 1) {
         if (N >= 128 * ii1 + 129) {
           if (N >= 128 * ii1 + 130) {
             A[1][128 * ii1] = (0.125 * ((A[0][128 * ii1 + 1] - (2.0 * A[0][128 * ii1])) + A[0][128 * ii1 - 1]));
           }
           for (register int i1 = -((N - 128 * ii1 - 3) % 127) + N - 3; i1 <= 128 * ii1 + 127; i1 += 127) {
             A[1][i1] = (0.125 * ((A[0][i1 + 1] - (2.0 * A[0][i1])) + A[0][i1 - 1]));
           }
         } else {
           if (N >= 128 * ii1 + 2) {
             A[1][128 * ii1] = (0.125 * ((A[0][128 * ii1 + 1] - (2.0 * A[0][128 * ii1])) + A[0][128 * ii1 - 1]));
           }
           A[1][N - 1] = (0.125 * ((A[0][0] - (2.0 * A[0][N - 1])) + A[0][N - 2]));
         }
       }
     }
   }
   if (T >= 2 && T <= 128 && (N + 127) % 128 >= 1) {
     if (2 * T >= ((N - 1) % 128) + 132) {
       for (register int i0 = 0; i0 < T - N + (N + 1) / 2 + 64 * ((N - 2) / 128) - 65; i0 += 1) {
         for (register int i1 = -T + i0 + 129; i1 < -((N - 2) % 128) + T - i0 - 3; i1 += 1) {
           A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
         }
         for (register int i1 = -T + N + i0 + 1; i1 <= T - i0 + 254; i1 += 1) {
           A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
         }
       }
       if (N >= 514) {
         for (register int ii1_3 = 2 * T + 381 >= N ? 3 : (N - 2) / 128 - 1; ii1_3 <= (2 * T + 381 >= N ? 3 : (N - 2) / 128 - 1); ii1_3 += 1) {
           if (N >= 2 * T + 382 || 1) {
             for (register int i0 = 0; i0 < T - N + 64 * ii1_3 + (N + 1) / 2 - 1; i0 += 1) {
               for (register int i1 = -T + N + i0 + 1; i1 < T + 128 * ii1_3 - i0 - 1; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           }
         }
       }
     }
     for (register int k = max(1, -((N - 1) / 128) + (-2 * T + N - (N - 1) / 128 + 2) / 127 + 2); k <= 2; k += 1) {
       if (k == 2) {
         for (register int i0 = 0; i0 < T - 1; i0 += 1) {
           for (register int i1 = max(T - i0 - 1, -T + i0 + 129); i1 <= min(T - i0 + 126, -T + i0 + 256); i1 += 1) {
             A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
           }
         }
         for (register int i0 = max(0, T - N + N / 2 + 64 * ((N - 2) / 128) - 64); i0 < -((N - 2) % 128) + T - 2; i0 += 1) {
           for (register int i1 = -((N - 2) % 128) + T - i0 - 3; i1 <= min(T - i0 - 2, -T + i0 + 128); i1 += 1) {
             A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 && (-T + N + i0 + 1) % 128 == 0 ? N - 1 : i1 - 1]));
           }
         }
         for (register int i0 = max(0, -((N - 2) % 128) + T - 2); i0 < T - 1; i0 += 1) {
           for (register int i1 = 0; i1 <= min(T - i0 - 2, -T + i0 + 128); i1 += 1) {
             A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
           }
           for (register int i1 = max(-T + N + i0 + 1, -((N - 2) % 128) + T + N - i0 - 3); i1 < N; i1 += 1) {
             A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
           }
         }
         #pragma omp parallel for
         for (register int ii1 = 1; ii1 < floord(N - 1, 128) - 1; ii1 += 1) {
           for (register int i0 = 0; i0 < T - 1; i0 += 1) {
             for (register int i1 = max(T + 128 * ii1 - i0 - 1, -T + 128 * ii1 + i0 + 129); i1 <= min(T + 128 * ii1 - i0 + 126, -T + 128 * ii1 + i0 + 256); i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
           }
         }
         for (register int i0 = max(0, T - N + N / 2 + 64 * ((N - 2) / 128) - 64); i0 < T - 1; i0 += 1) {
           for (register int i1 = max(-((N - 2) % 128) + T + N - i0 - 131, -((N - 2) % 128) - T + N + i0 - 1); i1 <= min(-T + N + i0, -((N - 2) % 128) + T + N - i0 - 4); i1 += 1) {
             A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
           }
         }
       } else {
         for (register int i0 = 0; i0 < T - 65; i0 += 1) {
           for (register int i1 = -T + i0 + 257; i1 <= T - i0 + 126; i1 += 1) {
             A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
           }
         }
         for (register int i0 = 0; i0 < T - 65; i0 += 1) {
           for (register int i1 = max(-T + i0 + 129, -((N - 2) % 128) + T - i0 - 3); i1 < T - i0 - 1; i1 += 1) {
             A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
           }
         }
         for (register int i0 = 0; i0 < T - N + (N + 1) / 2 + 64 * ((N - 2) / 128) - 1; i0 += 1) {
           for (register int i1 = 0; i1 <= min(-T + i0 + 128, -((N - 2) % 128) + T - i0 - 4); i1 += 1) {
             A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
           }
           for (register int i1 = max(-T + N + i0 + 1, -((N - 2) % 128) + T + N - i0 - 131); i1 < min(N, -((N - 2) % 128) + T + N - i0 - 3); i1 += 1) {
             A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
           }
         }
         #pragma omp parallel for
         for (register int ii1 = 1; ii1 < floord(N - 1, 128) - 1; ii1 += 1) {
           for (register int i0 = 0; i0 < T - 65; i0 += 1) {
             for (register int i1 = -T + 128 * ii1 + i0 + 257; i1 <= min(-T + N + i0, T + 128 * ii1 - i0 + 126); i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
           }
         }
       }
     }
     for (register int i0 = max(0, T - 64); i0 < T; i0 += 1) {
       for (register int i1 = T - i0 - 1; i1 <= -T + i0 + 128; i1 += 1) {
         A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i0 + 1 == T && i1 == 0 ? N - 1 : i1 - 1]));
       }
     }
     #pragma omp parallel for
     for (register int ii1 = 1; ii1 < floord(N - 1, 128) - 1; ii1 += 1) {
       for (register int i0 = max(max(0, T - 64), T - 64 * ii1 + 64); i0 < T; i0 += 1) {
         for (register int i1 = T + 128 * ii1 - i0 - 1; i1 <= -T + 128 * ii1 + i0 + 128; i1 += 1) {
           A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
         }
       }
       if (ii1 == 1) {
         for (register int i1 = T + 127; i1 <= -T + 256; i1 += 1) {
           A[1][i1] = (0.125 * ((A[0][i1 + 1] - (2.0 * A[0][i1])) + A[0][i1 - 1]));
         }
         if (T >= 65) {
           for (register int i1 = 191; i1 <= 192; i1 += 1) {
             A[-(T % 2) + 1][i1] = (0.125 * ((A[T % 2][i1 + 1] - (2.0 * A[T % 2][i1])) + A[T % 2][i1 - 1]));
           }
         }
         for (register int i0 = max(1, T - 63); i0 < T; i0 += 1) {
           for (register int i1 = T - i0 + 127; i1 <= -T + i0 + 256; i1 += 1) {
             A[-(i0 % 2) + 1][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
           }
         }
       }
     }
     #pragma omp parallel for
     for (register int ii1 = (N - 1) / 128 - 1; ii1 <= (N - 1) / 128; ii1 += 1) {
       for (register int i0 = max(max(0, T - 64), T - N + 64 * ii1 + N / 2); i0 < T; i0 += 1) {
         for (register int i1 = T + 128 * ii1 - i0 - 1; i1 <= min(-T + N + i0, -T + 128 * ii1 + i0 + 128); i1 += 1) {
           A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i0 + 1 == T && i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
         }
       }
     }
   } else if (T == 1) {
     #pragma omp parallel for
     for (register int ii1 = 0; ii1 <= (N - 1) / 128; ii1 += 1) {
       if (N >= 128 * ii1 + 257) {
         for (register int i1 = max(256, 128 * ii1); i1 <= 128 * ii1 + 127; i1 += 1) {
           A[1][i1] = (0.125 * ((A[0][i1 + 1] - (2.0 * A[0][i1])) + A[0][i1 - 1]));
         }
         if (ii1 == 1) {
           for (register int i1 = 128; i1 <= 255; i1 += 1) {
             A[1][i1] = (0.125 * ((A[0][i1 + 1] - (2.0 * A[0][i1])) + A[0][i1 - 1]));
           }
         } else if (ii1 == 0) {
           for (register int i1 = 0; i1 <= 127; i1 += 1) {
             A[1][i1] = (0.125 * ((A[0][i1 + 1] - (2.0 * A[0][i1])) + A[0][i1 == 0 ? N - 1 : i1 - 1]));
           }
         }
       } else {
         for (register int i1 = 128 * ii1; i1 <= min(N - 1, 128 * ii1 + 127); i1 += 1) {
           A[1][i1] = (0.125 * ((A[0][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[0][i1])) + A[0][i1 - 1]));
         }
       }
     }
   }
   if (T >= 66 && T <= 128 && (N - 1) % 128 == 0) {
     for (register int k = 0; k <= 3; k += 1) {
       if (k >= 1) {
         #pragma omp parallel for
         for (register int ii1 = 0; ii1 < k - 1; ii1 += 1) {
           if (k == 3) {
             for (register int i0 = T - 64; i0 < T; i0 += 1) {
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
           } else {
             for (register int i0 = 0; i0 < T - 1; i0 += 1) {
               for (register int i1 = max(T - i0 - 1, -T + i0 + 129); i1 <= min(T - i0 + 126, -T + i0 + 256); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           }
           if (k == 2 && ii1 == 0) {
             for (register int i0 = T - 65; i0 < T - 1; i0 += 1) {
               A[(i0 + 1) % 2][T - i0 - 2] = (0.125 * ((A[i0 % 2][T - i0 - 1] - (2.0 * A[i0 % 2][T - i0 - 2])) + A[i0 % 2][i0 + 2 == T ? N - 1 : T - i0 - 3]));
             }
           }
         }
         if (k == 3) {
           #pragma omp parallel for
           for (register int ii1 = 2; ii1 <= (N - 1) / 128; ii1 += 1) {
             for (register int i0 = max(T - 64, ((-N - 1) / 2) + T + 64 * ii1); i0 < T; i0 += 1) {
               for (register int i1 = T + 128 * ii1 - i0 - 1; i1 <= min(-T + N + i0, -T + 128 * ii1 + i0 + 128); i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][128 * ii1 + 1 == N && i0 + 1 == T && i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           }
         } else {
           #pragma omp parallel for
           for (register int ii1 = k - 1; ii1 < (N - 129) / 128; ii1 += 1) {
             if (k == 2) {
               for (register int i0 = 0; i0 < T - 1; i0 += 1) {
                 for (register int i1 = max(T + 128 * ii1 - i0 - 1, -T + 128 * ii1 + i0 + 129); i1 <= min(T + 128 * ii1 - i0 + 126, -T + 128 * ii1 + i0 + 256); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             } else {
               for (register int i0 = 0; i0 < T - 65; i0 += 1) {
                 for (register int i1 = -T + 128 * ii1 + i0 + 257; i1 <= min(-T + N + i0, T + 128 * ii1 - i0 + 126); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
               if (ii1 == 0) {
                 for (register int i0 = 0; i0 < T - 65; i0 += 1) {
                   A[(i0 + 1) % 2][T - i0 - 2] = (0.125 * ((A[i0 % 2][T - i0 - 1] - (2.0 * A[i0 % 2][T - i0 - 2])) + A[i0 % 2][T - i0 - 3]));
                 }
               }
             }
             if (k == 1 && ii1 == 0) {
               for (register int i0 = 0; i0 < T - 1; i0 += 1) {
                 for (register int i1 = 0; i1 <= min(T - i0 - 3, -T + i0 + 128); i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
                 }
                 for (register int i1 = max(T + N - i0 - 130, -T + N + i0 + 1); i1 < N; i1 += 1) {
                   A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
                 }
               }
             }
           }
           if (k == 2) {
             for (register int i0 = T - 65; i0 < T - 1; i0 += 1) {
               A[(i0 + 1) % 2][-T + N + i0] = (0.125 * ((A[i0 % 2][-T + N + i0 + 1] - (2.0 * A[i0 % 2][-T + N + i0])) + A[i0 % 2][-T + N + i0 - 1]));
             }
           }
         }
       } else {
         for (register int i0 = 0; i0 < T - 65; i0 += 1) {
           for (register int i1 = -T + i0 + 129; i1 < T - i0 - 2; i1 += 1) {
             A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
           }
           if (N == 385) {
             for (register int i1 = -T + i0 + 386; i1 <= T - i0 + 254; i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
           }
         }
         if (N >= 513) {
           for (register int i0 = 0; i0 < T - 65; i0 += 1) {
             for (register int i1 = -T + N + i0 + 1; i1 < T + N - i0 - 130; i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
           }
         }
       }
     }
   } else if (T >= 2 && T <= 65 && (N - 1) % 128 == 0) {
     for (register int k = 1; k <= 3; k += 1) {
       if (k >= 2) {
         for (register int i0 = max(0, T + k - 67); i0 < T + k - 3; i0 += 1) {
           if (k == 2) {
             for (register int i1 = -T + i0 + 129; i1 <= T - i0 + 126; i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
           } else {
             for (register int i1 = T - i0 - 1; i1 <= -T + i0 + 128; i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i0 + 1 == T && i1 == 0 ? N - 1 : i1 - 1]));
             }
           }
         }
         if (k == 2) {
           for (register int i0 = 0; i0 < T - 1; i0 += 1) {
             A[(i0 + 1) % 2][T - i0 - 2] = (0.125 * ((A[i0 % 2][T - i0 - 1] - (2.0 * A[i0 % 2][T - i0 - 2])) + A[i0 % 2][i0 + 2 == T ? N - 1 : T - i0 - 3]));
           }
         }
       } else {
         for (register int i0 = 0; i0 < T - 1; i0 += 1) {
           for (register int i1 = 0; i1 < T - i0 - 2; i1 += 1) {
             A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 == 0 ? N - 1 : i1 - 1]));
           }
           for (register int i1 = -T + N + i0 + 1; i1 < N; i1 += 1) {
             A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
           }
         }
       }
       if (k == 3) {
         #pragma omp parallel for
         for (register int ii1 = 1; ii1 <= (N - 1) / 128; ii1 += 1) {
           for (register int i0 = max(max(max(0, T - 64), T - 64 * ii1 + 64), ((-N - 1) / 2) + T + 64 * ii1); i0 < T; i0 += 1) {
             for (register int i1 = T + 128 * ii1 - i0 - 1; i1 <= min(-T + N + i0, -T + 128 * ii1 + i0 + 128); i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][128 * ii1 + 1 == N && i0 + 1 == T && i1 + 1 == N ? 0 : i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
           }
           if (ii1 == 1) {
             for (register int i0 = max(0, T - 64); i0 < T; i0 += 1) {
               for (register int i1 = T - i0 + 127; i1 <= -T + i0 + 256; i1 += 1) {
                 A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
               }
             }
           }
         }
       } else if (k == 2) {
         #pragma omp parallel for
         for (register int ii1 = 1; ii1 < (N - 129) / 128; ii1 += 1) {
           for (register int i0 = 0; i0 < T - 1; i0 += 1) {
             for (register int i1 = -T + 128 * ii1 + i0 + 129; i1 <= T + 128 * ii1 - i0 + 126; i1 += 1) {
               A[(i0 + 1) % 2][i1] = (0.125 * ((A[i0 % 2][i1 + 1] - (2.0 * A[i0 % 2][i1])) + A[i0 % 2][i1 - 1]));
             }
           }
         }
         for (register int i0 = 0; i0 < T - 1; i0 += 1) {
           A[(i0 + 1) % 2][-T + N + i0] = (0.125 * ((A[i0 % 2][-T + N + i0 + 1] - (2.0 * A[i0 % 2][-T + N + i0])) + A[i0 % 2][-T + N + i0 - 1]));
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
 