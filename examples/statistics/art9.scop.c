int min(int, int);

int main()
{
#if 1
#define N 128
#else
  int N;
#endif

  int path[N][N];

#pragma scop
  for (int k = 0; k < N; ++k) {
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < N; ++j) {
        if (i < k && j < k)
S1:       path[i][j] = min(path[i][j], path[i][k] + path[k][j]);
        if (i > k && j < k)
S2:       path[i][j] = min(path[i][j], path[i][k] + path[k][j]);
        if (i < k && j > k)
S3:       path[i][j] = min(path[i][j], path[i][k] + path[k][j]);
        if (i > k && j > k)
S4:       path[i][j] = min(path[i][j], path[i][k] + path[k][j]);
      }
    }
  }
#pragma endscop
}
