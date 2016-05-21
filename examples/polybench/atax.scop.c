double SCALAR_VAL(double);

int main()
{
#if 0
# define M 4
# define N 4
#else
  int M;
  int N;
#endif

  int i,j,k;
  
  double A[M][N];
  double x[N];
  double y[N];
  double tmp[M];

#pragma scop
  for (i = 0; i < N; i++) {
S1: y[i] = 0;
  }
  
  for (i = 0; i < M; i++) {
S2: tmp[i] = SCALAR_VAL(0.0);
    for (j = 0; j < N; j++) {
S3:   tmp[i] = tmp[i] + A[i][j] * x[j];
    }    
    for (j = 0; j < N; j++) {
S4:   y[j] = y[j] + A[i][j] * tmp[i];
    }
  }
#pragma endscop
}

