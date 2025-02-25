#include <omp.h>
#include <math.h>
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

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

  int t1, t2, t3, t4;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
for (t1=-1;t1<=625;t1++) {
  lbp=ceild(t1,2);
  ubp=floord(t1+62500,2);
#pragma omp parallel for private(lbv,ubv,t3,t4)
  for (t2=lbp;t2<=ubp;t2++) {
    if ((t1 >= 2) && (t1 == 2*t2)) {
      if (t1%2 == 0) {
        h[0] = (h[0] - (0.69999999999999996 * (e[0 + 1] - e[0])));;
      }
    }
    if ((t1 == 2*t2-62500) && (t1 >= 2)) {
      if (t1%2 == 0) {
        h[999999] = (h[999999] - (0.69999999999999996 * (e[999999 + 1] - e[999999])));;
      }
    }
    for (t3=max(max(1,16*t1),32*t1-32*t2+1);t3<=min(min(10000,16*t1+15),32*t1-32*t2+31);t3++) {
      for (t4=max(32*t2,t3+1);t4<=-32*t1+32*t2+2*t3;t4++) {
        e[(-t3+t4)] = (e[(-t3+t4)] - (0.5 * (h[(-t3+t4)] - h[(-t3+t4) - 1])));;
        h[(-t3+t4-1)] = (h[(-t3+t4-1)] - (0.69999999999999996 * (e[(-t3+t4-1) + 1] - e[(-t3+t4-1)])));;
      }
      h[(-32*t1+32*t2+t3)] = (h[(-32*t1+32*t2+t3)] - (0.69999999999999996 * (e[(-32*t1+32*t2+t3) + 1] - e[(-32*t1+32*t2+t3)])));;
    }
    if (t1 == 2*t2) {
      for (t3=16*t1+16;t3<=min(10000,16*t1+30);t3++) {
        for (t4=t3+1;t4<=16*t1+31;t4++) {
          if (t1%2 == 0) {
            e[(-t3+t4)] = (e[(-t3+t4)] - (0.5 * (h[(-t3+t4)] - h[(-t3+t4) - 1])));;
          }
          if (t1%2 == 0) {
            h[(-t3+t4-1)] = (h[(-t3+t4-1)] - (0.69999999999999996 * (e[(-t3+t4-1) + 1] - e[(-t3+t4-1)])));;
          }
        }
      }
    }
    for (t3=max(max(1,16*t1),32*t1-32*t2+32);t3<=min(min(10000,16*t1+15),32*t1-32*t2+999999);t3++) {
      for (t4=32*t2;t4<=-32*t1+32*t2+2*t3;t4++) {
        e[(-t3+t4)] = (e[(-t3+t4)] - (0.5 * (h[(-t3+t4)] - h[(-t3+t4) - 1])));;
        h[(-t3+t4-1)] = (h[(-t3+t4-1)] - (0.69999999999999996 * (e[(-t3+t4-1) + 1] - e[(-t3+t4-1)])));;
      }
      h[(-32*t1+32*t2+t3)] = (h[(-32*t1+32*t2+t3)] - (0.69999999999999996 * (e[(-32*t1+32*t2+t3) + 1] - e[(-32*t1+32*t2+t3)])));;
    }
    if (t1 == 2*t2-62500) {
      for (t3=16*t1+1;t3<=16*t1+15;t3++) {
        for (t4=16*t1+1000000;t4<=t3+999999;t4++) {
          if (t1%2 == 0) {
            e[(-t3+t4)] = (e[(-t3+t4)] - (0.5 * (h[(-t3+t4)] - h[(-t3+t4) - 1])));;
          }
          if (t1%2 == 0) {
            h[(-t3+t4-1)] = (h[(-t3+t4-1)] - (0.69999999999999996 * (e[(-t3+t4-1) + 1] - e[(-t3+t4-1)])));;
          }
        }
        if (t1%2 == 0) {
          h[999999] = (h[999999] - (0.69999999999999996 * (e[999999 + 1] - e[999999])));;
        }
      }
    }
    if (t1 == 2*t2-62500) {
      for (t3=16*t1+16;t3<=min(10000,16*t1+29);t3++) {
        if (t1%2 == 0) {
          e[(-16*t1+t3+999969)] = (e[(-16*t1+t3+999969)] - (0.5 * (h[(-16*t1+t3+999969)] - h[(-16*t1+t3+999969) - 1])));;
        }
        for (t4=-16*t1+2*t3+999970;t4<=t3+999999;t4++) {
          if (t1%2 == 0) {
            e[(-t3+t4)] = (e[(-t3+t4)] - (0.5 * (h[(-t3+t4)] - h[(-t3+t4) - 1])));;
          }
          if (t1%2 == 0) {
            h[(-t3+t4-1)] = (h[(-t3+t4-1)] - (0.69999999999999996 * (e[(-t3+t4-1) + 1] - e[(-t3+t4-1)])));;
          }
        }
        if (t1%2 == 0) {
          h[999999] = (h[999999] - (0.69999999999999996 * (e[999999 + 1] - e[999999])));;
        }
      }
    }
    if (t1 >= 2*t2-62499) {
      for (t3=max(max(1,16*t1+16),32*t1-32*t2+32);t3<=min(10000,16*t1+30);t3++) {
        e[(-32*t1+32*t2+t3-31)] = (e[(-32*t1+32*t2+t3-31)] - (0.5 * (h[(-32*t1+32*t2+t3-31)] - h[(-32*t1+32*t2+t3-31) - 1])));;
        for (t4=-32*t1+32*t2+2*t3-30;t4<=32*t2+31;t4++) {
          e[(-t3+t4)] = (e[(-t3+t4)] - (0.5 * (h[(-t3+t4)] - h[(-t3+t4) - 1])));;
          h[(-t3+t4-1)] = (h[(-t3+t4-1)] - (0.69999999999999996 * (e[(-t3+t4-1) + 1] - e[(-t3+t4-1)])));;
        }
      }
    }
    if ((t1 == 2*t2-62500) && (t1 <= 622)) {
      if (t1%2 == 0) {
        e[999999] = (e[999999] - (0.5 * (h[999999] - h[999999 - 1])));;
      }
      if (t1%2 == 0) {
        h[999999] = (h[999999] - (0.69999999999999996 * (e[999999 + 1] - e[999999])));;
      }
    }
    if ((t1 <= min(623,2*t2-1)) && (t1 >= 2*t2-62499)) {
      e[(-16*t1+32*t2)] = (e[(-16*t1+32*t2)] - (0.5 * (h[(-16*t1+32*t2)] - h[(-16*t1+32*t2) - 1])));;
    }
  }
}

    IF_TIME(t_end = rtclock());
    IF_TIME(printf("%0.6lfs\n", t_end - t_start));

    if (fopen(".test", "r")) {
        print_array();
    }

    return 0;
}
