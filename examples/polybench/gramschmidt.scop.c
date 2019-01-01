double SQRT_FUN(double);
double SCALAR_VAL(double);

int main()
{
#if 0
# define _PB_M 2000
# define _PB_N 2600
#else
  int _PB_N;
  int _PB_M;
#endif
  
  double A[_PB_M][_PB_N];
  double R[_PB_N][_PB_N];
  double Q[_PB_M][_PB_N];
  double nrm;

#pragma scop
  for (int k = 0; k < _PB_N; k++) {
S1: nrm = SCALAR_VAL(0.0);
    for (int i = 0; i < _PB_M; i++) {
S2:   nrm += A[i][k] * A[i][k];
    }
      
S3: R[k][k] = SQRT_FUN(nrm);
    
    for (int i = 0; i < _PB_M; i++) {
S4:   Q[i][k] = A[i][k] / R[k][k];
    }
      
    for (int j = k + 1; j < _PB_N; j++) {
S5:   R[k][j] = SCALAR_VAL(0.0);
      for (int i = 0; i < _PB_M; i++) {
S6:     R[k][j] += Q[i][k] * A[i][j];
      }
      for (int i = 0; i < _PB_M; i++) {
S7:     A[i][j] = A[i][j] - Q[i][k] * R[k][j];
      }
    }
  }
#pragma endscop
}
