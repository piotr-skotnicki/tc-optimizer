int main()
{
#if 1
#define N 9
#define M 9
#else
  int N;
  int M;
#endif

  int A[N+1][M+1];

#pragma scop
  for (int i = 1; i <= N; ++i) {
    for (int j = 1; j <= M; ++j) {
S1:   A[i][j] = A[i-1][j];
    }
  }
#pragma endscop
}
