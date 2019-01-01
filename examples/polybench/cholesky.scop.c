double SQRT_FUN(double);

int main()
{
#if 0
# define _PB_N 4000
#else
  int _PB_N;
#endif

  double A[_PB_N][_PB_N];

#pragma scop
  for (int i = 0; i < _PB_N; i++) {
    for (int j = 0; j < i; j++) {
      for (int k = 0; k < j; k++) {
S1:     A[i][j] -= A[i][k] * A[j][k];
      }
S2:   A[i][j] /= A[j][j];
    }
    for (int k = 0; k < i; k++) {
S3:   A[i][i] -= A[i][k] * A[i][k];
    }
S4: A[i][i] = SQRT_FUN(A[i][i]);
  }
#pragma endscop
}
