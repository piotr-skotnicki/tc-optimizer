double SCALAR_VAL(double);

int main()
{
#if 0
# define _PB_TSTEPS 20
# define _PB_N 10
#else
  int _PB_TSTEPS;
  int _PB_N;
#endif
  
  double A[_PB_N][_PB_N][_PB_N];
  double B[_PB_N][_PB_N][_PB_N];

#pragma scop
  for (int t = 1; t <= _PB_TSTEPS; t++) {
    for (int i = 1; i < _PB_N-1; i++) {
      for (int j = 1; j < _PB_N-1; j++) {
        for (int k = 1; k < _PB_N-1; k++) {
S1:         B[i][j][k] =   SCALAR_VAL(0.125) * (A[i+1][j][k] - SCALAR_VAL(2.0) * A[i][j][k] + A[i-1][j][k])
                         + SCALAR_VAL(0.125) * (A[i][j+1][k] - SCALAR_VAL(2.0) * A[i][j][k] + A[i][j-1][k])
                         + SCALAR_VAL(0.125) * (A[i][j][k+1] - SCALAR_VAL(2.0) * A[i][j][k] + A[i][j][k-1])
                         + A[i][j][k];
        }
      }
    }
    for (int i = 1; i < _PB_N-1; i++) {
      for (int j = 1; j < _PB_N-1; j++) {
        for (int k = 1; k < _PB_N-1; k++) {
S2:       A[i][j][k] =   SCALAR_VAL(0.125) * (B[i+1][j][k] - SCALAR_VAL(2.0) * B[i][j][k] + B[i-1][j][k])
                       + SCALAR_VAL(0.125) * (B[i][j+1][k] - SCALAR_VAL(2.0) * B[i][j][k] + B[i][j-1][k])
                       + SCALAR_VAL(0.125) * (B[i][j][k+1] - SCALAR_VAL(2.0) * B[i][j][k] + B[i][j][k-1])
                       + B[i][j][k];
        }
      }
    }
  }
#pragma endscop
}
