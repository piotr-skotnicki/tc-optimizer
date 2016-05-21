int main()
{
#if 0
# define M 1024
#else
  int M;
#endif

  int i,j,k;
  
  double cos1[M][M+13];
  double temp2d[M][M+23];
  double block[M][M+43];
  double sum2[M][M];

#pragma scop
  for (i = 0; i < M; i++) {
    for (j = 0; j < M; j++) {
S1:   temp2d[i][j] = 0.0;
      for (k = 0; k < M; k++) {
S2:     temp2d[i][j] = temp2d[i][j] + block[i][k] * cos1[j][k];
      }
    }
  }

  for (i = 0; i < M; i++) {
    for (j = 0; j < M; j++) {
S3:   sum2[i][j] = 0.0;
      for (k = 0; k < M; k++) {
S4:     sum2[i][j] = sum2[i][j] + cos1[i][k] * temp2d[k][j];
      }
S5:   block[i][j] = (sum2[i][j]);
    }
  }
#pragma endscop
}

