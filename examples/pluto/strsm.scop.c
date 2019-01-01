int main()
{
#if 0
# define N 1000
#else
  int N;
#endif
  
  double a[N][N];
  double b[N][N];

#pragma scop
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) {
      for (int k=i+1; k<N; k++) {
        if (k == i+1)
S1:       b[j][i] /= a[i][i];
S2:     b[j][k] -= a[i][k] * b[j][i];
      }
    }
  }
#pragma endscop
}
