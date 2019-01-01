int b2s23(int cell, int neighbors);

int main()
{
#if 0
# define T 500
# define N 2000
#else
  int T;
  int N;
#endif
  
  int life[2][N][N];

#pragma scop
  for (int t = 0; t < T; t++) {
    for (int i = 1; i < N-1; i++) {
      for (int j = 1; j < N-1; j++) {
S1:     life[(t+1)%2][i][j] = b2s23(life[t%2][i][j], life[t%2][i-1][j+1] + life[t%2][i-1][j] + life[t%2][i-1][j-1] 
                    + life[t%2][i][j+1] + life[t%2][i][j-1]                     
                    + life[t%2][i+1][j+1] + life[t%2][i+1][j] + life[t%2][i+1][j-1]);
      }
    }
  }
#pragma endscop
}
