int main()
{
#if 0
# define NMAX 1000
#else
  int NMAX;
#endif

  int i,j,k;
  
  double a[NMAX][NMAX];
  double b[NMAX][NMAX];

#pragma scop
  for (i=1; i<NMAX; i++) {
    for (j=0; j<NMAX; j++) {
      for (k=0; k<i; k++) {
S1:     b[j][k] += a[i][k] * b[j][i];
      }
    }
  }
#pragma endscop
}

