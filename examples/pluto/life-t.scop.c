int b2s23(int cell, int neighbors);

int main()
{
#if 0
# define N 16000
# define T 500
#else
  int N;
  int T;
#endif

  int life[T+1][N][N];

#pragma scop
  for (int t = 0; t < T; t++) {
    for (int i = 1; i < N-1; i++) {
      for (int j = 1; j < N-1; j++) {
S1:     life[t+1][i][j] = b2s23(life[t][i][j], life[t][i-1][j+1] + life[t][i-1][j] + life[t][i-1][j-1]
                    + life[t][i][j+1] + life[t][i][j-1]
                    + life[t][i+1][j+1] + life[t][i+1][j] + life[t][i+1][j-1]);
      }
    }
  }
#pragma endscop
}
