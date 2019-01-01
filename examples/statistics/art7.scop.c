int main()
{
#if 1
#define N 1024
#else
  int N;
#endif

  int a[N][N];
  int b[N][N];
  int c[N][N];
  int d[N][N];
  
#pragma scop
  for (int k = 0; k < N; ++k) {
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < N; ++j) {
S1:     a[i][j] = a[i][j];
S2:     b[i][j] = b[i][j] + a[i][j];
S3:     c[i][j] = c[i][j];
S4:     d[i][j] = d[i][j] - c[i][j];
      }
    }
  }
#pragma endscop
}
