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

  if (_PB_N >= 1) {
    for (int k1 = 1; k1 <= _PB_N; k1 += 1) {
      for (int k2 = 1; k2 <= _PB_N; k2 += 1) {
S1:     Pbp[k1][k2] = (((Q[1][k1 - 1] * Qbp[k1][k2]) * Q[k2 + 1][_PB_N]) / Q[1][_PB_N]);
      }
    }
#pragma scop
    for (int k1 = 2; k1 <= _PB_N; k1 += 1) {
      for (int k2 = -_PB_N + 1; k2 < 0; k2 += 1) {
        for (int k3 = 1; k3 < k1; k3 += 1) {
          for (int k4 = -k2 + 1; k4 <= _PB_N; k4 += 1) {
S2:         Pbp[k1][-k2] += (paired(RNA[k3 - 1], RNA[k4 - 1]) * (((((Pbp[k3][k4] * ERT) * Q[k3 + 1][k1 - 1]) * Qbp[k1][-k2]) * Q[-k2 + 1][k4 - 1]) / ((Qbp[k3][k4] == 0.0) ? 1.0 : Qbp[k3][k4])));
          }
        }
      }
    }
#pragma endscop
  }

}
