int main()
{
#if 1
# define N 1600
# define T 500
#else
  int N;
  int T;
#endif

  double A[T+1][N][N];
int _N = N - 1;
int _T = T;

#pragma scop
  for (int t = 0; t < _T; t++) {
    for (int i = 0; i <= _N; i++) {
      for (int j = 0; j <= _N; j++) {
S1:     A[t + 1][i][j] =
            0.125 * (A[t][i == _N ? 0 : i + 1][j] - 2.0 * A[t][i][j] +
                     A[t][(i == 0) ? (_N) : (i - 1)][j]) +
            0.125 * (A[t][i][j == _N ? 0 : j + 1] - 2.0 * A[t][i][j] +
                     A[t][i][(j == 0) ? (_N) : (j - 1)]) +
            A[t][i][j];
      }
    }
  }
#pragma endscop
}
