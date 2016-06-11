double SCALAR_VAL(double);

int main()
{
#if 0
# define _PB_N 1000
#else
  int _PB_N;
#endif

  int i,j;
  
  double alpha;
  double beta;
  double A[_PB_N][_PB_N];
  double B[_PB_N][_PB_N];
  double tmp[_PB_N];
  double x[_PB_N];
  double y[_PB_N];

#pragma scop
  for (i = 0; i < _PB_N; i++) {
S1: tmp[i] = SCALAR_VAL(0.0);
S2: y[i] = SCALAR_VAL(0.0);
    for (j = 0; j < _PB_N; j++) {
S3:   tmp[i] = A[i][j] * x[j] + tmp[i];
S4:   y[i] = B[i][j] * x[j] + y[i];
    }
S5: y[i] = alpha * tmp[i] + beta * y[i];
  }
#pragma endscop
}

