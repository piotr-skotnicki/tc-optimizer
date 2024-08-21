int MIN(int, int);

int main()
{
#if 0
# define N 100
#else
  int N;
#endif

  int c[N+10][N+10];
  int w[N+10][N+10];

#pragma scop
  for (int i = N-1; i >= 1; i--) {
    for (int j = i+1; j <= N; j++) {
      for (int k = i+1; k < j; k++) {
S1:     c[i][j] = MIN(c[i][j], w[i][j] + c[i][k] + c[k][j]);
      }
    }
  }
#pragma endscop
}
