int paired(char, char);

int main()
{
#if 0
# define _PB_N 100
#else
  int _PB_N;
#endif

  int l;
  char RNA[_PB_N];
  int C[_PB_N+1][_PB_N+1];

#pragma scop
  if (_PB_N >= 1 && l >= 0 && l <= 5) {
    for (int i = _PB_N-1; i >= 1; i--) {
      for (int j = i+1; j <= _PB_N; j++) {
S1:     C[i][j] = C[i][j-1];
        for (int k = i; k < j-l; k++) {
S2:       C[i][j] += paired(RNA[k-1], RNA[j-1]) * C[i][k-1] * C[k+1][j-1];
        }
      }
    }
  }
#pragma endscop

}
