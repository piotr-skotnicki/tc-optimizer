int main()
{
#if 0
# define _PB_M 1000
# define _PB_N 1000
#else
  int _PB_M;
  int _PB_N;
#endif
  
  double temp;
  double alpha;
  double A[_PB_M][_PB_M];
  double B[_PB_M][_PB_N];

#pragma scop
  for (int i = 0; i < _PB_M; i++) {
    for (int j = 0; j < _PB_N; j++) {
      for (int k = i+1; k < _PB_M; k++) {
S1:     B[i][j] += A[k][i] * B[k][j];
      }
S2:   B[i][j] = alpha * B[i][j];
    }
  }
#pragma endscop
}
