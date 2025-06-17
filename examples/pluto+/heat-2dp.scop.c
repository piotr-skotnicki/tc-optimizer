int main()
{
#if 0
# define N 1600
# define T 500
#else
  int N;
  int T;
#endif

  double A[2][N][N];
  
short _N = N - 1;
short _T = T;

#pragma scop
  for (int t = 0; t < _T; t++) {
    for (int i = 0; i <= _N; i++) {
      for (int j = 0; j <= _N; j++) {
S1:     A[(t + 1) % 2][i][j] =
            0.125 * (A[t % 2][i == _N ? 0 : i + 1][j] - 2.0 * A[t % 2][i][j] +
                     A[t % 2][(i == 0) ? _N : (i - 1)][j]) +
            0.125 * (A[t % 2][i][j == _N ? 0 : j + 1] - 2.0 * A[t % 2][i][j] +
                     A[t % 2][i][(j == 0) ? _N : (j - 1)]) +
            A[t % 2][i][j];
      }
    }
  }
#pragma endscop

}
