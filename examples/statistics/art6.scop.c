int min(int, int);

int main()
{
#if 0
#define N 128
#else
  int N;
#endif

  int path[N][N];

#pragma scop
  for (int k = 0; k < N; ++k) {
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < N; ++j) {
        path[i][j] = path[i][j] < path[i][k] + path[k][j] ?
                     path[i][j] : path[i][k] + path[k][j];
      }
    }
  }
#pragma endscop
}
