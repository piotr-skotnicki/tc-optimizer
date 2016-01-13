double SCALAR_VAL(double);

int main()
{
#if 0
# define N 1000
# define M 1000
#else
  int N;
  int M;
#endif

  int i,j,k;
  
  double float_n;
  double data[N][M];
  double cov[M][M];
  double mean[M];

#pragma scop
  for (j = 0; j < M; j++) {
S1: mean[j] = SCALAR_VAL(0.0);
    for (i = 0; i < N; i++) {
S2:   mean[j] += data[i][j];
    }
S3: mean[j] /= float_n;
  }

  for (i = 0; i < N; i++) {
    for (j = 0; j < M; j++) {
S4:   data[i][j] -= mean[j];
    }
  }

  for (i = 0; i < M; i++) {
    for (j = i; j < M; j++) {
S5:   cov[i][j] = SCALAR_VAL(0.0);
      for (k = 0; k < N; k++) {
S6:     cov[i][j] += data[k][i] * data[k][j];
	    }
S7:   cov[i][j] /= (float_n - SCALAR_VAL(1.0));
S8:   cov[j][i] = cov[i][j];
    }
  }
#pragma endscop
}

