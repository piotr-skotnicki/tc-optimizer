int main()
{
#if 0
# define N 300
# define T 200
#else
  int N;
  int T;
#endif

  double A[2][N][N][N];
  short _N = N - 1;
  short _T = T;

#pragma scop
  for (int t = 0; t < _T; t++) {
    for (int i = 0; i <= _N; i++) {
      for (int j = 0; j <= _N; j++) {
        for (int k = 0; k <= _N; k++) {
S1:       A[(t + 1) % 2][i][j][k] =
              0.125 * (A[t % 2][i == _N ? 0 : i + 1][j][k] -
                       2.0 * A[t % 2][i][j][k] +
                       A[t % 2][(i == 0) ? (_N) : (i - 1)][j][k]) +
              0.125 * (A[t % 2][i][j == _N ? 0 : j + 1][k] -
                       2.0 * A[t % 2][i][j][k] +
                       A[t % 2][i][(j == 0) ? (_N) : (j - 1)][k]) +
              0.125 * (A[t % 2][i][j][(k == 0) ? (_N) : (k - 1)] -
                       2.0 * A[t % 2][i][j][k] +
                       A[t % 2][i][j][k == _N ? 0 : k + 1]) +
              A[t % 2][i][j][k];
        }
      }
    }
  }
#pragma endscop
}
