int main()
{
#if 0
# define N 150
# define T 100
#else
  int T;
  int N;
#endif

  int A[T+1][N][N][N];

#pragma scop
  for (int t = 0; t < T; ++t) {
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < N; ++j) {
        for (int k = 0; k < N; ++k) {
S1:       A[t+1][i][j][k] = 0.125 * (A[t][i == N-1 ? 0 : i+1][j][k] - 2.0 * A[t][i][j][k] + A[t][i == 0 ? N-1 : i-1][j][k])
                          + 0.125 * (A[t][i][j == N-1 ? 0 : j+1][k] - 2.0 * A[t][i][j][k] + A[t][i][j == 0 ? N-1 : j-1][k])
                          + 0.125 * (A[t][i][j][k == N-1 ? 0 : k+1] - 2.0 * A[t][i][j][k] + A[t][i][j][k == 0 ? N-1 : k-1])
                          + A[t][i][j][k];
        }
      }
    }
  }
#pragma endscop
}
