int main()
{
#if 0
# define _PB_M 2000
# define _PB_N 2600
#else
  int _PB_M;
  int _PB_N;
#endif
  
  int i,j,k;
  
  double A[_PB_N][_PB_M];
  double B[_PB_N][_PB_M];
  double C[_PB_N][_PB_N];
  double alpha,beta;

#pragma scop
  for (i = 0; i < _PB_N; i++) {
    for (j = 0; j < _PB_N; j++) {
S1:   C[i][j] *= beta;
    }
  }
      
  for (i = 0; i < _PB_N; i++) {
    for (k = 0; k < _PB_M; k++) {
      for (j = 0; j < _PB_N; j++) {
S2:     C[i][j] += A[j][k] * alpha * B[i][k] + B[j][k] * alpha * A[i][k];
      }
    }
  }
#pragma endscop
}

