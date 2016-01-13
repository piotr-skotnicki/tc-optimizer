double SQRT_FUN(double);

int main()
{
#if 0
# define N 1000
#else
  int N;
#endif

  int i,j,k;

  double A[N][N];

#pragma scop
  for (i = 0; i < N; i++) {
    for (j = 0; j < i; j++) {
      for (k = 0; k < j; k++) {
S1:     A[i][j] -= A[i][k] * A[j][k];
      }
S2:   A[i][j] /= A[j][j];
    }
    for (k = 0; k < i; k++) {
S3:   A[i][i] -= A[i][k] * A[i][k];
    }
S4: A[i][i] = SQRT_FUN(A[i][i]);
  }
#pragma endscop
}

