int MIN(int, int);

int main()
{
#if 0
# define N 100
#else
  int N;
#endif

  int W[N+2][N+2];
  int V[N+2][N+2];
  int EHF[N+2][N+2];
  int EFL[N+2][N+2];

#pragma scop
  for (int i = N-1; i >= 0; i--) {
    for (int j = i+1; j < N; j++) {
      for (int k = i+1; k < j; k++) {
        for (int m = k+1; m < j; m++) {
          if (k - i + j - m > 2 && k - i + j - m < 30) {
S1:         V[i][j] = MIN(V[k][m] + EFL[i][j], V[i][j]);
          }
        }
S2:     W[i][j] += MIN(MIN(W[i][k], W[k+1][j]), W[i][j]);
        if (k < j-1) {
S3:       V[i][j] = MIN(W[i+1][k] + W[k+1][j-1], V[i][j]);
        }
      }
S4:   V[i][j] = MIN(MIN(V[i+1][j-1], EHF[i][j]), V[i][j]);
S5:   W[i][j] = MIN(MIN(MIN(W[i+1][j], W[i][j-1]), V[i][j]), W[i][j]);
    }
  }
#pragma endscop
}
