double SCALAR_VAL(double);

int main()
{
#if 0
# define _PB_NI 1000
# define _PB_NJ 1000
# define _PB_NK 1000
# define _PB_NL 1000
#else
  int _PB_NI;
  int _PB_NJ;
  int _PB_NK;
  int _PB_NL;
#endif

  int i,j,k;
  
  double alpha;
  double beta;
  double tmp[_PB_NI][_PB_NJ];
  double A[_PB_NI][_PB_NK];
  double B[_PB_NK][_PB_NJ];
  double C[_PB_NJ][_PB_NL];
  double D[_PB_NI][_PB_NL];

#pragma scop
  for (i = 0; i < _PB_NI; i++) {
    for (j = 0; j < _PB_NJ; j++) {
S1:   tmp[i][j] = SCALAR_VAL(0.0);
      for (k = 0; k < _PB_NK; ++k) {
S2:     tmp[i][j] += alpha * A[i][k] * B[k][j];
      }
    }
  }
  for (i = 0; i < _PB_NI; i++) {
    for (j = 0; j < _PB_NL; j++) {
S3:   D[i][j] *= beta;
      for (k = 0; k < _PB_NJ; ++k) {
S4:     D[i][j] += tmp[i][k] * C[k][j];
      }
    }
  }
#pragma endscop
}

