int main()
{
#if 0
# define M 1024
# define N 1024
# define K 1024
#else
  int M;
  int N;
  int K;
#endif
  
  double alpha;
  double beta;
  double A[M][K];
  double B[K][N];
  double C[M][N];

#pragma scop
  for (int i = 0; i < M; i++)
    for (int j = 0; j < N; j++)  
      for (int k = 0; k < K; k++)
        C[i][j] = beta * C[i][j] + alpha * A[i][k] * B[k][j];
#pragma endscop
}
