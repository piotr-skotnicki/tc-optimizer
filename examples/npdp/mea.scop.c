int main()
{
#if 0
# define N 100
#else
  int N;
#endif

  double Q[2*N+2][2*N+2];
  double Qbp[2*N+2][2*N+2];
  double Pbp[2*N+2][2*N+2];
  float ERT;

#pragma scop
  for (int i = 0; i < N; i++) {
    for (int j = i+1; j < N; j++) {
S1:   Pbp[i][j] = (Q[0][i] * Q[j][N-1] * Qbp[i][j]) / Q[0][N-1];
      for (int p = 0; p < i; p++) {
        for (int q = j+1; q < N; q++) {
S2:       Pbp[i][j] += (Pbp[p][q] * ERT * Q[p+1][i] * Qbp[i][j] * Q[j+1][q-1]) / (Qbp[p][q] == 0 ? 1 : Qbp[p][q]);
        }
      }
    }
  }
#pragma endscop
}
