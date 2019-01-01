int main()
{
#if 0
# define _PB_M 2000
# define _PB_N 2600
#else
  int _PB_M;
  int _PB_N;
#endif
  
  double A[_PB_N][_PB_M];
  double B[_PB_N][_PB_M];
  double C[_PB_N][_PB_N];
  double alpha,beta;

#pragma scop
  for (int i = 0; i < _PB_N; i++) {
    for (int j = 0; j < _PB_N; j++) {
S1:   C[i][j] *= beta;
    }
  }
      
  for (int i = 0; i < _PB_N; i++) {
    for (int k = 0; k < _PB_M; k++) {
      for (int j = 0; j < _PB_N; j++) {
S2:     C[i][j] += A[j][k] * alpha * B[i][k] + B[j][k] * alpha * A[i][k];
      }
    }
  }
#pragma endscop
}
