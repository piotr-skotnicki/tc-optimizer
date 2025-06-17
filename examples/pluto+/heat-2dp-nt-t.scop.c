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

#pragma scop
  for (int t = 0; t < T; t++) {
    for (int i = 0; i <= (N - 1); i++) {
      for (int j = 0; j <= (N - 1); j++) {
S1:     A[t + 1][i][j] =
            0.125 * (A[t][i == (N - 1) ? 0 : i + 1][j] - 2.0 * A[t][i][j] +
                     A[t][(i == 0) ? (N - 1) : (i - 1)][j]) +
            0.125 * (A[t][i][j == (N - 1) ? 0 : j + 1] - 2.0 * A[t][i][j] +
                     A[t][i][(j == 0) ? (N - 1) : (j - 1)]) +
            A[t][i][j];
      }
    }
  }
#pragma endscop
}
