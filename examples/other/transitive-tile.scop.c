int main()
{
#if 0
#define M 4000
#define N 4000
#else
  int M;
  int N;
#endif

  int A[M+1][N+1];

#pragma scop
  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < N; ++j) {
S1:   A[i][j] = A[i+1][j == N-1 ? 0 : j+1];
    }
  }
#pragma endscop
}
