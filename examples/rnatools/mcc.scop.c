int paired(char, char);

int main()
{
#if 0
# define _PB_N 100
#else
  int _PB_N;
#endif

  int l;
  double ERT;
  char RNA[_PB_N];
  double Q[_PB_N+2][_PB_N+2];
  double Qbp[_PB_N+2][_PB_N+2];

#pragma scop
  if (_PB_N >= 1 && l >= 0 && l <= 5) {
    for (int i = _PB_N; i >= 1; i--) {
      for (int j = i+1; j <= _PB_N; j++) {
S1:     Q[i][j] = Q[i][j-1];
        for (int k = i; k < j-l; k++) {
S2:       Qbp[i][j] = paired(RNA[i-1], RNA[j-1]) * Q[i+1][j-1] * ERT;
S3:       Q[i][j] += Q[i][k-1] * Qbp[k][j];
        }
      }
    }
  }
#pragma endscop

}
