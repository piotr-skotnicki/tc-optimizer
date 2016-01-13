double SCALAR_VAL(double);

int main()
{
#if 0
# define TSTEPS 1000
# define N 1000
#else
  int TSTEPS;
  int N;
#endif

  int t,i,j,k;
  
  double A[N][N][N];
  double B[N][N][N];

#pragma scop
  for (t = 1; t <= TSTEPS; t++) {
    for (i = 1; i < N-1; i++) {
      for (j = 1; j < N-1; j++) {
        for (k = 1; k < N-1; k++) {
S1:         B[i][j][k] =   SCALAR_VAL(0.125) * (A[i+1][j][k] - SCALAR_VAL(2.0) * A[i][j][k] + A[i-1][j][k])
                         + SCALAR_VAL(0.125) * (A[i][j+1][k] - SCALAR_VAL(2.0) * A[i][j][k] + A[i][j-1][k])
                         + SCALAR_VAL(0.125) * (A[i][j][k+1] - SCALAR_VAL(2.0) * A[i][j][k] + A[i][j][k-1])
                         + A[i][j][k];
        }
      }
    }
    for (i = 1; i < N-1; i++) {
      for (j = 1; j < N-1; j++) {
        for (k = 1; k < N-1; k++) {
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

