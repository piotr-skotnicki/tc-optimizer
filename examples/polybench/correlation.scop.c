double SCALAR_VAL(double);
double SQRT_FUN(double);

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
  
  double eps;
  double float_n;
  double data[N][M];
  double corr[M][M];
  double mean[M];
  double stddev[M];

#pragma scop
  for (j = 0; j < M; j++) {
S1: mean[j] = SCALAR_VAL(0.0);
    for (i = 0; i < N; i++) {
S2:   mean[j] += data[i][j];
    }
S3: mean[j] /= float_n;
  }

  for (j = 0; j < M; j++) {
S4: stddev[j] = SCALAR_VAL(0.0);
    for (i = 0; i < N; i++) {
S5:   stddev[j] += (data[i][j] - mean[j]) * (data[i][j] - mean[j]);
    }
S6: stddev[j] /= float_n;
S7: stddev[j] = SQRT_FUN(stddev[j]);
S8: stddev[j] = stddev[j] <= eps ? SCALAR_VAL(1.0) : stddev[j];
  }

  for (i = 0; i < N; i++) {
    for (j = 0; j < M; j++) {
S9:   data[i][j] -= mean[j];
S10:  data[i][j] /= SQRT_FUN(float_n) * stddev[j];
    }
  }

  for (i = 0; i < M-1; i++) {
S11:corr[i][i] = SCALAR_VAL(1.0);
    for (j = i+1; j < M; j++) {
S12:  corr[i][j] = SCALAR_VAL(0.0);
      for (k = 0; k < N; k++) {
S13:    corr[i][j] += (data[k][i] * data[k][j]);
      }
S14:  corr[j][i] = corr[i][j];
    }
  }
  
S15:corr[M-1][M-1] = SCALAR_VAL(1.0);
#pragma endscop
}

