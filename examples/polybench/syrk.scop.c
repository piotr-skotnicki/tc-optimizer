int main()
{
#if 0
# define M 1000
# define N 1000
#else
  int M;
  int N;
#endif

  int i,j,k;
  
  double alpha;
  double beta;
  double C[N][N];
  double A[N][M];

#pragma scop
  for (i = 0; i < N; i++) {
    for (j = 0; j <= i; j++) {
S1:   C[i][j] *= beta;
    }
    for (k = 0; k < M; k++) {
      for (j = 0; j <= i; j++) {
S2:     C[i][j] += alpha * A[i][k] * A[j][k];
      }
    }
  }
#pragma endscop
}

