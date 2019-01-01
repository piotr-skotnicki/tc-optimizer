double SCALAR_VAL(double);

int main()
{
#if 0
# define _PB_N 1000
# define _PB_M 1000
#else
  int _PB_N;
  int _PB_M;
#endif
  
  double float_n;
  double data[_PB_N][_PB_M];
  double cov[_PB_M][_PB_M];
  double mean[_PB_M];

#pragma scop
  for (int j = 0; j < _PB_M; j++) {
S1: mean[j] = SCALAR_VAL(0.0);
    for (int i = 0; i < _PB_N; i++) {
S2:   mean[j] += data[i][j];
    }
S3: mean[j] /= float_n;
  }

  for (int i = 0; i < _PB_N; i++) {
    for (int j = 0; j < _PB_M; j++) {
S4:   data[i][j] -= mean[j];
    }
  }

  for (int i = 0; i < _PB_M; i++) {
    for (int j = i; j < _PB_M; j++) {
S5:   cov[i][j] = SCALAR_VAL(0.0);
      for (int k = 0; k < _PB_N; k++) {
S6:     cov[i][j] += data[k][i] * data[k][j];
	    }
S7:   cov[i][j] /= (float_n - SCALAR_VAL(1.0));
S8:   cov[j][i] = cov[i][j];
    }
  }
#pragma endscop
}
