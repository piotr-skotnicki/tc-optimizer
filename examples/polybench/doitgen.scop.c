double SCALAR_VAL(double);

int main()
{
#if 0
# define NR 1000
# define NQ 1000
# define NP 1000
#else
  int NR;
  int NQ;
  int NP;
#endif
  
  int r,q,p,s;
  
  double A[NR][NQ][NP];
  double C4[NP][NP];
  double sum[NP];

#pragma scop
  for (r = 0; r < NR; r++) {
    for (q = 0; q < NQ; q++) {
      for (p = 0; p < NP; p++) {
S1:     sum[p] = SCALAR_VAL(0.0);
        for (s = 0; s < NP; s++)
S2:       sum[p] += A[r][q][s] * C4[s][p];
      }
      for (p = 0; p < NP; p++)
S3:     A[r][q][p] = sum[p];
    }
  }
#pragma endscop
}

