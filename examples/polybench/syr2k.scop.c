int main()
{
#if 0
# define M 2000
# define N 2600
#else
  int M;
  int N;
#endif
  
  int i,j,k;
  
  double A[N][M];
  double B[N][M];
  double C[N][N];
  double alpha,beta;

#pragma scop
  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
S1:   C[i][j] *= beta;
      
  for (i = 0; i < N; i++) {
    for (k = 0; k < M; k++) {
      for (j = 0; j < N; j++) {
S2:     C[i][j] += A[j][k] * alpha*B[i][k] + B[j][k] * alpha*A[i][k];
	  }
    }
  }
#pragma endscop
}

