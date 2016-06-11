int main()
{
#if 0
# define _PB_M 1000
# define _PB_N 1000
#else
  int _PB_M;
  int _PB_N;
#endif

  int i,j,k;
  
  double alpha;
  double beta;
  double C[_PB_N][_PB_N];
  double A[_PB_N][_PB_M];

#pragma scop
  for (i = 0; i < _PB_N; i++) {
    for (j = 0; j <= i; j++) {
S1:   C[i][j] *= beta;
    }
    for (k = 0; k < _PB_M; k++) {
      for (j = 0; j <= i; j++) {
S2:     C[i][j] += alpha * A[i][k] * A[j][k];
      }
    }
  }
#pragma endscop
}

