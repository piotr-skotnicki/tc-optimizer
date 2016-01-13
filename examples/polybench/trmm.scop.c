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
  
  double temp;
  double alpha;
  double A[M][M];
  double B[M][N];

#pragma scop
  for (i = 0; i < M; i++) {
    for (j = 0; j < N; j++) {
      for (k = i+1; k < M; k++) {
S1:     B[i][j] += A[k][i] * B[k][j];
      }
S2:   B[i][j] = alpha * B[i][j];
    }
  }
#pragma endscop
}

