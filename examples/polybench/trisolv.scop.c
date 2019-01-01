int main()
{
#if 0
# define _PB_N 1000
#else
  int _PB_N;
#endif
  
  double L[_PB_N][_PB_N];
  double x[_PB_N];
  double b[_PB_N];

#pragma scop
  for (int i = 0; i < _PB_N; i++) {
S1: x[i] = b[i];
    for (int j = 0; j < i; j++) {
S2:   x[i] -= L[i][j] * x[j];
    }  
S3: x[i] = x[i] / L[i][i];
  }
#pragma endscop
}
