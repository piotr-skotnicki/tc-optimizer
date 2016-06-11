int main()
{
#if 0
# define _PB_N 4000
#else
  int _PB_N;
#endif

  int i,j,k;
  
  double A[_PB_N][_PB_N];

#pragma scop
  for (i = 0; i < _PB_N; i++) {
    for (j = 0; j < i; j++) {
      for (k = 0; k < j; k++) {
S1:     A[i][j] -= A[i][k] * A[k][j];
      }
S2:   A[i][j] /= A[j][j];
    }
    for (j = i; j < _PB_N; j++) {
      for (k = 0; k < i; k++) {
S3:     A[i][j] -= A[i][k] * A[k][j];
      }
    }
  }
#pragma endscop
}

