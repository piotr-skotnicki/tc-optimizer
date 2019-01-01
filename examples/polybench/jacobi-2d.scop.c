double SCALAR_VAL(double);

int main()
{
#if 0
# define _PB_TSTEPS 1000
# define _PB_N 1000
#else
  int _PB_TSTEPS;
  int _PB_N;
#endif

  double A[_PB_N][_PB_N];
  double B[_PB_N][_PB_N];

#pragma scop
  for (int t = 0; t < _PB_TSTEPS; t++) {
    for (int i = 1; i < _PB_N - 1; i++)
      for (int j = 1; j < _PB_N - 1; j++)
S1:     B[i][j] = SCALAR_VAL(0.2) * (A[i][j] + A[i][j-1] + A[i][1+j] + A[1+i][j] + A[i-1][j]);
        
    for (int i = 1; i < _PB_N - 1; i++)
      for (int j = 1; j < _PB_N - 1; j++)
S2:     A[i][j] = SCALAR_VAL(0.2) * (B[i][j] + B[i][j-1] + B[i][1+j] + B[1+i][j] + B[i-1][j]);
  }
#pragma endscop
}
