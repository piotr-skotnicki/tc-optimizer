double SCALAR_VAL(double);

int main()
{
#if 0
# define N 1000
#else
  int N;
#endif

  int i,j;
  
  double alpha;
  double beta;
  double A[N][N];
  double B[N][N];
  double tmp[N];
  double x[N];
  double y[N];

#pragma scop
  for (i = 0; i < N; i++) {
S1: tmp[i] = SCALAR_VAL(0.0);
S2: y[i] = SCALAR_VAL(0.0);
    for (j = 0; j < N; j++) {
S3:   tmp[i] = A[i][j] * x[j] + tmp[i];
S4:   y[i] = B[i][j] * x[j] + y[i];
    }
S5: y[i] = alpha * tmp[i] + beta * y[i];
  }
#pragma endscop
}

