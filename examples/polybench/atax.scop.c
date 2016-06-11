double SCALAR_VAL(double);

int main()
{
#if 0
# define _PB_M 4
# define _PB_N 4
#else
  int _PB_M;
  int _PB_N;
#endif

  int i,j,k;
  
  double A[_PB_M][_PB_N];
  double x[_PB_N];
  double y[_PB_N];
  double tmp[_PB_M];

#pragma scop
  for (i = 0; i < _PB_N; i++) {
S1: y[i] = 0;
  }
  
  for (i = 0; i < _PB_M; i++) {
S2: tmp[i] = SCALAR_VAL(0.0);
    for (j = 0; j < _PB_N; j++) {
S3:   tmp[i] = tmp[i] + A[i][j] * x[j];
    }    
    for (j = 0; j < _PB_N; j++) {
S4:   y[j] = y[j] + A[i][j] * tmp[i];
    }
  }
#pragma endscop
}

