int main()
{
#if 0
# define _PB_N 256
#else
  int _PB_N;
#endif
  
  double path[_PB_N][_PB_N];

#pragma scop
  for (int k = 0; k < _PB_N; k++) {
    for (int i = 0; i < _PB_N; i++) {
      for (int j = 0; j < _PB_N; j++) {
S1:     path[i][j] = path[i][j] < path[i][k] + path[k][j] ? path[i][j] : path[i][k] + path[k][j];
      }
    }
  }
#pragma endscop
}
