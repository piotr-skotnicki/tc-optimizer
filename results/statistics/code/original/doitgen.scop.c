double SCALAR_VAL(double);

int main()
{
#if 0
# define _PB_NR 1000
# define _PB_NQ 1000
# define _PB_NP 1000
#else
  int _PB_NR;
  int _PB_NQ;
  int _PB_NP;
#endif
  
  int r,q,p,s;
  
  double A[_PB_NR][_PB_NQ][_PB_NP];
  double C4[_PB_NP][_PB_NP];
  double sum[_PB_NP];

#pragma scop
  for (r = 0; r < _PB_NR; r++) {
    for (q = 0; q < _PB_NQ; q++) {
      for (p = 0; p < _PB_NP; p++) {
S1:     sum[p] = SCALAR_VAL(0.0);
        for (s = 0; s < _PB_NP; s++)
S2:       sum[p] += A[r][q][s] * C4[s][p];
      }
      for (p = 0; p < _PB_NP; p++)
S3:     A[r][q][p] = sum[p];
    }
  }
#pragma endscop
}

