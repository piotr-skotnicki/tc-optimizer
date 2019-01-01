int main()
{
#if 0
# define N 4000
# define T 1000
#else
  int T;
  int N;
#endif

  int A[2][N+2][N+2];

#pragma scop
  for (int t = 0; t < T; ++t) {
    for (int i = 1; i < N+1; ++i) {
      for (int j = 1; j < N+1; ++j) {
S1:     A[(t+1)%2][i][j] = 0.125 * (A[t%2][i+1][j] - 2.0 * A[t%2][i][j] + A[t%2][i-1][j])
                         + 0.125 * (A[t%2][i][j+1] - 2.0 * A[t%2][i][j] + A[t%2][i][j-1])
                         + A[t%2][i][j];
      }
    }
  }
#pragma endscop
}
