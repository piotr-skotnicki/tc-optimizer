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
  double c[NMAX][NMAX];

#pragma scop
  for (i=0; i<NMAX; i++) {
    for (j=0; j<NMAX; j++) {
      for (k=0; k<j-1; k++) {
S1:     c[i][k] += a[j][k] * b[i][j];
S2:     c[i][j] += a[j][j] * b[i][j];
      }
S3:   c[i][j] += a[j][j] * b[i][j];
    }
  }
#pragma endscop
}

