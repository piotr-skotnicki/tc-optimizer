int main()
{
#if 0
# define N 2000
# define T 1000
#else
  short N;
  short T;
#endif

  double a[N][N];
  double b[N][N];

#pragma scop
  for (int t = 0; t < T; t++) {
    for (int i = 2; i < N-1; i++) {
      for (int j = 2; j < N-1; j++) {
S1:     b[i][j] = 0.2 * (a[i][j] + a[i][j-1] + a[i][1+j] + a[1+i][j] + a[i-1][j]);
      }
    }
    for (int i = 2; i < N-1; i++) {
      for (int j = 2; j < N-1; j++) {
S2:     a[i][j] = b[i][j];
      }
    }
  }
#pragma endscop
}
