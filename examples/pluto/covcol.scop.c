int main()
{
#if 0
# define M 1000
# define N 1000
#else
  int M;
  int N;
#endif
  
  double data[M+1][N+1];
  double symmat[M+1][M+1];

#pragma scop
  for (int j1 = 1; j1 <= M; j1++) {
    for (int j2 = j1; j2 <= M; j2++) {
      for (int i = 1; i <= N; i++) {
S1:     symmat[j1][j2] += data[i][j1] * data[i][j2];
      }
S2:   symmat[j2][j1] = symmat[j1][j2];
    }
  }
#pragma endscop
}
