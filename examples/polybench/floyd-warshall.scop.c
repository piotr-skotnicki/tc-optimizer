int main()
{
#if 0
# define N 1000
#else
  int N;
#endif

  int i,j,k;
  
  double path[N][N];

#pragma scop
  for (k = 0; k < N; k++) {
    for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++) {
S1:     path[i][j] = path[i][j] < path[i][k] + path[k][j] ? path[i][j] : path[i][k] + path[k][j];
      }
    }
  }
#pragma endscop
}

