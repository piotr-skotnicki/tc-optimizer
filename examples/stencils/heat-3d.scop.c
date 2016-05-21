int main()
{
#if 0
# define N 150
# define T 100
#else
  int T;
  int N;
#endif
  
  int t,i,j,k;
  
  int A[2][N+2][N+2][N+2];
  
#pragma scop
  for (t = 0; t < T; ++t) {
    for (i = 1; i < N+1; ++i) {
      for (j = 1; j < N+1; ++j) {
        for (k = 1; k < N+1; ++k) {
S1:       A[(t+1)%2][i][j][k] = 0.125 * (A[t%2][i+1][j][k] - 2.0 * A[t%2][i][j][k] + A[t%2][i-1][j][k])
                              + 0.125 * (A[t%2][i][j+1][k] - 2.0 * A[t%2][i][j][k] + A[t%2][i][j-1][k])
                              + 0.125 * (A[t%2][i][j][k-1] - 2.0 * A[t%2][i][j][k] + A[t%2][i][j][k+1])
                              + A[t%2][i][j][k];
        }
      }
    }
  }
#pragma endscop
}

