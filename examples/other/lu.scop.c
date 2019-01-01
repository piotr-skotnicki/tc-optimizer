int main()
{
#if 0
# define _PB_N 4000
#else
  int _PB_N;
#endif
  
  double A[_PB_N][_PB_N];

#pragma scop
  for (int k = 0; k < _PB_N; k++) {
    for (int i = 0; i < 1; i++) {
      for (int j = k+1; j < _PB_N; j++) {
S1:     A[k][j] = A[k][j] / A[k][k];
      }
    }
    for (int i = k+1; i < _PB_N; i++) {
      for (int j = k+1; j < _PB_N; j++) {
S2:     A[i][j] = A[i][j] - A[i][k]*A[k][j];
      }
    }
  }
#pragma endscop
}
