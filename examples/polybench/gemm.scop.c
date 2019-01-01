int main()
{
#if 0
# define _PB_NI 1000
# define _PB_NJ 1000
# define _PB_NK 1000
#else
  int _PB_NI;
  int _PB_NJ;
  int _PB_NK;
#endif
  
  double alpha;
  double beta;
  double C[_PB_NI][_PB_NJ];
  double A[_PB_NI][_PB_NK];
  double B[_PB_NK][_PB_NJ];

#pragma scop
  for (int i = 0; i < _PB_NI; i++) {
    for (int j = 0; j < _PB_NJ; j++) {
S1:   C[i][j] *= beta;
    }
    for (int k = 0; k < _PB_NK; k++) {
      for (int j = 0; j < _PB_NJ; j++) {
S2:     C[i][j] += alpha * A[i][k] * B[k][j];
      }
    }
  }
#pragma endscop
}
