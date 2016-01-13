double SQRT_FUN(double);
double SCALAR_VAL(double);

int main()
{
#if 0
# define N 1000
# define M 1000
#else
  int N;
  int M;
#endif

  int i,j,k;
  
  double A[M][N];
  double R[N][N];
  double Q[M][N];
  double nrm;

#pragma scop
  for (k = 0; k < N; k++) {
S1: nrm = SCALAR_VAL(0.0);
    for (i = 0; i < M; i++) {
S2:   nrm += A[i][k] * A[i][k];
    }
      
S3: R[k][k] = SQRT_FUN(nrm);
    
    for (i = 0; i < M; i++) {
S4:   Q[i][k] = A[i][k] / R[k][k];
    }
      
    for (j = k + 1; j < N; j++) {
S5:   R[k][j] = SCALAR_VAL(0.0);
	  for (i = 0; i < M; i++) {
S6:     R[k][j] += Q[i][k] * A[i][j];
      }
      for (i = 0; i < M; i++) {
S7:     A[i][j] = A[i][j] - Q[i][k] * R[k][j];
      }
    }
  }
#pragma endscop
}

