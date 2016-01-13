int main()
{
#if 0
# define NI 1000
# define NJ 1000
# define NK 1000
#else
  int NI;
  int NJ;
  int NK;
#endif

  int i,j,k;
  
  double alpha;
  double beta;
  double C[NI][NJ];
  double A[NI][NK];
  double B[NK][NJ];

#pragma scop
  for (i = 0; i < NI; i++) {
    for (j = 0; j < NJ; j++) {
S1:   C[i][j] *= beta;
    }
    for (k = 0; k < NK; k++) {
      for (j = 0; j < NJ; j++) {
S2:     C[i][j] += alpha * A[i][k] * B[k][j];
      }
    }
  }
#pragma endscop
}

