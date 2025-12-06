int main()
{
#if 1
#define N 6
#define M 9
#else
  int N;
  int M;
#endif

  int A[N+2][M+2];

#pragma scop
  for (int i = 1; i <= N; ++i) {
    for (int j = 1; j <= M; ++j) {
S1:   A[i][j] = A[i][j+1] + A[i+1][j];
    }
  }
#pragma endscop
}
