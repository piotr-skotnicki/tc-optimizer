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
  
  double temp2;
  double alpha;
  double beta;
  double C[M][N];
  double A[M][M];
  double B[M][N];

#pragma scop
   for (i = 0; i < M; i++) {
     for (j = 0; j < N; j++) {
S1:    temp2 = 0;
       for (k = 0; k < i; k++) {
S2:      C[k][j] += alpha * B[i][j] * A[i][k];
S3:      temp2 += B[k][j] * A[i][k];
       }
S4:    C[i][j] = beta * C[i][j] + alpha*B[i][j] * A[i][i] + alpha * temp2;
     }
   } 
#pragma endscop
}

