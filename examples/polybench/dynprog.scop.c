int main()
{
#if 0
# define _PB_TSTEPS 4000
# define _PB_LENGTH 4000
#else
  int _PB_TSTEPS;
  int _PB_LENGTH;
#endif

  double c[_PB_LENGTH][_PB_LENGTH];
  double W[_PB_LENGTH][_PB_LENGTH];
  double sum_c[_PB_LENGTH][_PB_LENGTH][_PB_LENGTH];
  double out_l;

#pragma scop
  for (int iter = 0; iter < _PB_TSTEPS; iter++) {
    for (int i = 0; i <= _PB_LENGTH - 1; i++) {
      for (int j = 0; j <= _PB_LENGTH - 1; j++) {
S1:     c[i][j] = 0;
      }
    }

    for (int i = 0; i <= _PB_LENGTH - 2; i++) {
      for (int j = i + 1; j <= _PB_LENGTH - 1; j++) {
S2:     sum_c[i][j][i] = 0;
        for (int k = i + 1; k <= j-1; k++) {
S3:       sum_c[i][j][k] = sum_c[i][j][k - 1] + c[i][k] + c[k][j];
        }
S4:     c[i][j] = sum_c[i][j][j-1] + W[i][j];
      }
    }
S5: out_l += c[0][_PB_LENGTH - 1];
  }
#pragma endscop

}
