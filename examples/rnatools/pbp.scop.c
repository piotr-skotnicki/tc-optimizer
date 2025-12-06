int paired(char, char);

int main()
{
#if 0
# define _PB_N 100
#else
  int _PB_N;
#endif

  double ERT;
  char RNA[_PB_N];
  double Q[_PB_N+2][_PB_N+2];
  double Pbp[_PB_N+2][_PB_N+2];
  double Qbp[_PB_N+2][_PB_N+2];

#pragma scop
  if (_PB_N >= 1) {
    for (int i = 1; i <= _PB_N; i++) {
      for (int j = 1; j <= _PB_N; j++) {
S1:     Pbp[i][j] = (Q[1][i-1] * Qbp[i][j] * Q[j+1][_PB_N]) / Q[1][_PB_N];
        for (int p = 1; p < i; p++) {
          for (int q = j+1; q <= _PB_N; q++) {
S2:         Pbp[i][j] += paired(RNA[p-1], RNA[q-1]) * ((Pbp[p][q] * ERT * Q[p+1][i-1] * Qbp[i][j] * Q[j+1][q-1]) / (Qbp[p][q] == 0.0 ? 1.0 : Qbp[p][q]));
          }
        }
      }
    }
  }
#pragma endscop

}
