double SCALAR_VAL(double);

int main()
{
#if 0
# define NI 1000
# define NJ 1000
# define NK 1000
# define NL 1000
# define NM 1000
#else
  int NI;
  int NJ;
  int NK;
  int NL;
  int NM;
#endif

  int i,j,k;
  
  double E[NI][NJ];
  double A[NI][NK];
  double B[NK][NJ];
  double F[NJ][NL];
  double C[NJ][NM];
  double D[NM][NL];
  double G[NI][NL];

#pragma scop
  for (i = 0; i < NI; i++) {
    for (j = 0; j < NJ; j++) {
S1:   E[i][j] = SCALAR_VAL(0.0);
      for (k = 0; k < NK; ++k) {
S2:     E[i][j] += A[i][k] * B[k][j];
      }
    }
  }
    
  for (i = 0; i < NJ; i++) {
    for (j = 0; j < NL; j++) {
S3:   F[i][j] = SCALAR_VAL(0.0);
      for (k = 0; k < NM; ++k) {
S4:     F[i][j] += C[i][k] * D[k][j];
      }
    }
  }

  for (i = 0; i < NI; i++) {
    for (j = 0; j < NL; j++) {
S5:   G[i][j] = SCALAR_VAL(0.0);
      for (k = 0; k < NJ; ++k) {
S6:     G[i][j] += E[i][k] * F[k][j];
      }
    }
  }
#pragma endscop
}

