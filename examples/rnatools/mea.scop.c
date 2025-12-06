int paired(char, char);
double max_score(double, double);

int main()
{
#if 0
# define _PB_N 100
#else
  int _PB_N;
#endif

  int l;
  double gamma;
  char RNA[_PB_N];
  double Pu[_PB_N+1][_PB_N+1];
  double Pbp[_PB_N+2][_PB_N+2];
  double M[_PB_N+1][_PB_N+1];

#pragma scop
  if (_PB_N >= 1 && l >= 0 && l <= 5) {
    for (int i = _PB_N; i >= 1; i--) {
      for (int j = i; j <= _PB_N; j++) {
S1:     M[i][j] = M[i][j-1] + Pu[j][j];
        for (int k = i; k < j-l; k++) {
S2:       M[i][j] = max_score(M[i][j], paired(RNA[k-1], RNA[j-1]) * (M[i][k-1] + M[k+1][j-1] + gamma * Pbp[k][j]));
        }
      }
    }
  }
#pragma endscop

}
