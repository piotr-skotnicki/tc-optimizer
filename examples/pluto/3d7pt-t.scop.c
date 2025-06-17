int main()
{
#if 0
# define N 256
# define T 200
#else
  int N;
  int T;
#endif

  double A[T][N+2][N+2][N+2];
  const double alpha = 0.0876;
  const double beta = 0.0765;

#pragma scop
  for (int t = 0; t < T-1; t++) {
      for (int i = 1; i < N+1; i++) {
          for (int j = 1; j < N+1; j++) {
              for (int k = 1; k < N+1; k++) {
S1:               A[t+1][i][j][k] = alpha * (A[t][i][j][k])
                      + beta * (A[t][i-1][j][k] + A[t][i][j-1][k] + A[t][i][j][k-1] +
                                A[t][i+1][j][k] + A[t][i][j+1][k] + A[t][i][j][k+1]);
              }
          }
      }
  }
#pragma endscop
}
