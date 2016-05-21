int main()
{
#if 0
# define N 1000
#else
  int N;
#endif

  int i,j,k;
  
  double a[N][N];
  double b[N][N];

#pragma scop
  for (i=0; i<N; i++) {
    for (j=0; j<N; j++) {
      for (k=i+1; k<N; k++) {
        if (k == i+1)
S1:       b[j][i] /= a[i][i];
S2:     b[j][k] -= a[i][k] * b[j][i];
      }
    }
  }
#pragma endscop
}

