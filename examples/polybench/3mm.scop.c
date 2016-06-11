double SCALAR_VAL(double);

int main()
{
#if 0
# define _PB_NI 1000
# define _PB_NJ 1000
# define _PB_NK 1000
# define _PB_NL 1000
# define _PB_NM 1000
#else
  int _PB_NI;
  int _PB_NJ;
  int _PB_NK;
  int _PB_NL;
  int _PB_NM;
#endif

  int i,j,k;
  
  double E[_PB_NI][_PB_NJ];
  double A[_PB_NI][_PB_NK];
  double B[_PB_NK][_PB_NJ];
  double F[_PB_NJ][_PB_NL];
  double C[_PB_NJ][_PB_NM];
  double D[_PB_NM][_PB_NL];
  double G[_PB_NI][_PB_NL];

#pragma scop
  for (i = 0; i < _PB_NI; i++) {
    for (j = 0; j < _PB_NJ; j++) {
S1:   E[i][j] = SCALAR_VAL(0.0);
      for (k = 0; k < _PB_NK; ++k) {
S2:     E[i][j] += A[i][k] * B[k][j];
      }
    }
  }
    
  for (i = 0; i < _PB_NJ; i++) {
    for (j = 0; j < _PB_NL; j++) {
S3:   F[i][j] = SCALAR_VAL(0.0);
      for (k = 0; k < _PB_NM; ++k) {
S4:     F[i][j] += C[i][k] * D[k][j];
      }
    }
  }

  for (i = 0; i < _PB_NI; i++) {
    for (j = 0; j < _PB_NL; j++) {
S5:   G[i][j] = SCALAR_VAL(0.0);
      for (k = 0; k < _PB_NJ; ++k) {
S6:     G[i][j] += E[i][k] * F[k][j];
      }
    }
  }
#pragma endscop
}

