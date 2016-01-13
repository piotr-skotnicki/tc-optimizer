double SCALAR_VAL(double);

int main()
{
#if 0
# define NI 1000
# define NJ 1000
# define NK 1000
# define NL 1000
#else
  int NI;
  int NJ;
  int NK;
  int NL;
#endif

  int i,j,k;
  
  double alpha;
  double beta;
  double tmp[NI][NJ];
  double A[NI][NK];
  double B[NK][NJ];
  double C[NJ][NL];
  double D[NI][NL];

#pragma scop
  for (i = 0; i < NI; i++) {
    for (j = 0; j < NJ; j++) {
S1:   tmp[i][j] = SCALAR_VAL(0.0);
      for (k = 0; k < NK; ++k) {
S2:     tmp[i][j] += alpha * A[i][k] * B[k][j];
      }
    }
  }
  for (i = 0; i < NI; i++) {
    for (j = 0; j < NL; j++) {
S3:   D[i][j] *= beta;
      for (k = 0; k < NJ; ++k) {
S4:     D[i][j] += tmp[i][k] * C[k][j];
      }
    }
  }
#pragma endscop
}

