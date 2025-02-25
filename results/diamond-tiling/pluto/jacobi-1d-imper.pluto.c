#include <omp.h>
#include <math.h>
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

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

  int t1, t2, t3, t4;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
for (t1=-1;t1<=124;t1++) {
  lbp=ceild(t1,2);
  ubp=floord(t1+125000,2);
#pragma omp parallel for private(lbv,ubv,t3,t4)
  for (t2=lbp;t2<=ubp;t2++) {
    if (t1 == 2*t2) {
      if (t1%2 == 0) {
        b[2] = (0.33333 * ((a[2 - 1] + a[2]) + a[2 + 1]));;
      }
      for (t4=16*t1+5;t4<=16*t1+6;t4++) {
        if (t1%2 == 0) {
          a[(-16*t1+t4-3)] = b[(-16*t1+t4-3)];;
        }
      }
    }
    for (t3=max(0,8*t1);t3<=min(min(8*t1+7,16*t2-2),16*t1-16*t2+999999);t3++) {
      for (t4=32*t2;t4<=-32*t1+32*t2+4*t3;t4++) {
        b[(-2*t3+t4)] = (0.33333 * ((a[(-2*t3+t4) - 1] + a[(-2*t3+t4)]) + a[(-2*t3+t4) + 1]));;
        a[(-2*t3+t4-1)] = b[(-2*t3+t4-1)];;
      }
      for (t4=-32*t1+32*t2+4*t3+1;t4<=min(2*t3+1999999,-32*t1+32*t2+4*t3+2);t4++) {
        a[(-2*t3+t4-1)] = b[(-2*t3+t4-1)];;
      }
    }
    if (t1 == 2*t2-125000) {
      for (t3=8*t1+1;t3<=8*t1+7;t3++) {
        for (t4=16*t1+2000000;t4<=2*t3+1999998;t4++) {
          if (t1%2 == 0) {
            b[(-2*t3+t4)] = (0.33333 * ((a[(-2*t3+t4) - 1] + a[(-2*t3+t4)]) + a[(-2*t3+t4) + 1]));;
          }
          if (t1%2 == 0) {
            a[(-2*t3+t4-1)] = b[(-2*t3+t4-1)];;
          }
        }
        if (t1%2 == 0) {
          a[1999998] = b[1999998];;
        }
      }
    }
    for (t3=max(max(0,16*t2-1),16*t1-16*t2+2);t3<=8*t1+7;t3++) {
      b[2] = (0.33333 * ((a[2 - 1] + a[2]) + a[2 + 1]));;
      for (t4=2*t3+3;t4<=-32*t1+32*t2+4*t3;t4++) {
        b[(-2*t3+t4)] = (0.33333 * ((a[(-2*t3+t4) - 1] + a[(-2*t3+t4)]) + a[(-2*t3+t4) + 1]));;
        a[(-2*t3+t4-1)] = b[(-2*t3+t4-1)];;
      }
      for (t4=-32*t1+32*t2+4*t3+1;t4<=-32*t1+32*t2+4*t3+2;t4++) {
        a[(-2*t3+t4-1)] = b[(-2*t3+t4-1)];;
      }
    }
    for (t3=8*t1+8;t3<=min(min(999,16*t2+14),16*t1-16*t2+16);t3++) {
      b[2] = (0.33333 * ((a[2 - 1] + a[2]) + a[2 + 1]));;
      for (t4=2*t3+3;t4<=32*t2+31;t4++) {
        b[(-2*t3+t4)] = (0.33333 * ((a[(-2*t3+t4) - 1] + a[(-2*t3+t4)]) + a[(-2*t3+t4) + 1]));;
        a[(-2*t3+t4-1)] = b[(-2*t3+t4-1)];;
      }
    }
    for (t3=8*t1+8;t3<=min(min(999,16*t2-999984),16*t1-16*t2+1000013);t3++) {
      for (t4=-32*t1+32*t2+4*t3-31;t4<=-32*t1+32*t2+4*t3-30;t4++) {
        b[(-2*t3+t4)] = (0.33333 * ((a[(-2*t3+t4) - 1] + a[(-2*t3+t4)]) + a[(-2*t3+t4) + 1]));;
      }
      for (t4=-32*t1+32*t2+4*t3-29;t4<=2*t3+1999998;t4++) {
        b[(-2*t3+t4)] = (0.33333 * ((a[(-2*t3+t4) - 1] + a[(-2*t3+t4)]) + a[(-2*t3+t4) + 1]));;
        a[(-2*t3+t4-1)] = b[(-2*t3+t4-1)];;
      }
      a[1999998] = b[1999998];;
    }
    for (t3=max(max(8*t1+8,16*t2-999983),16*t1-16*t2+17);t3<=min(999,8*t1+15);t3++) {
      for (t4=-32*t1+32*t2+4*t3-31;t4<=-32*t1+32*t2+4*t3-30;t4++) {
        b[(-2*t3+t4)] = (0.33333 * ((a[(-2*t3+t4) - 1] + a[(-2*t3+t4)]) + a[(-2*t3+t4) + 1]));;
      }
      for (t4=-32*t1+32*t2+4*t3-29;t4<=32*t2+31;t4++) {
        b[(-2*t3+t4)] = (0.33333 * ((a[(-2*t3+t4) - 1] + a[(-2*t3+t4)]) + a[(-2*t3+t4) + 1]));;
        a[(-2*t3+t4-1)] = b[(-2*t3+t4-1)];;
      }
    }
    if ((t1 == 2*t2-125000) && (t1 <= 122)) {
      for (t4=16*t1+2000025;t4<=16*t1+2000026;t4++) {
        if (t1%2 == 0) {
          b[(-16*t1+t4-28)] = (0.33333 * ((a[(-16*t1+t4-28) - 1] + a[(-16*t1+t4-28)]) + a[(-16*t1+t4-28) + 1]));;
        }
      }
      if (t1%2 == 0) {
        a[1999998] = b[1999998];;
      }
    }
  }
}

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
