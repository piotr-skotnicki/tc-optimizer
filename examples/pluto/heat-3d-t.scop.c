int main()
{
#if 0
# define N 150
# define T 100
#else
  int N;
  int T;
#endif

  double A[T][N+2][N+2][N+2];

#pragma scop
  for (int t = 0; t < T - 1; t++) {
    for (int i = 1; i < N + 1; i++) {
      for (int j = 1; j < N + 1; j++) {
        for (int k = 1; k < N + 1; k++) {
S1:       A[t + 1][i][j][k] =
              0.125 * (A[t][i + 1][j][k] - 2.0 * A[t][i][j][k] +
                       A[t][i - 1][j][k]) +
              0.125 * (A[t][i][j + 1][k] - 2.0 * A[t][i][j][k] +
                       A[t][i][j - 1][k]) +
              0.125 * (A[t][i][j][k - 1] - 2.0 * A[t][i][j][k] +
                       A[t][i][j][k + 1]) +
              A[t][i][j][k];
        }
      }
    }
  }
#pragma endscop
}
