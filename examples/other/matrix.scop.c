int main()
{
#if 0
# define M 4000
# define N 4000
# define K 4000
#else
  int M;
  int N;
  int K;
#endif
  
  double A[M][K];
  double B[K][N];
  double C[M][N];

#pragma scop
  for (int i = 0; i < M ; ++i) {
    for (int j = 0; j < N ; ++j) {
S1:   C[i][j] = 0;
      for (int k = 0; k < K; ++k) {
S2:     C[i][j] = C[i][j] + (A[i][k] * B[k][j]);
      }
    }
  }
#pragma endscop
}
