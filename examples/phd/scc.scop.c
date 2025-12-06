int main()
{
#if 0
#define N 6
#define M 6
#else
  int N;
  int M;
#endif

  int A[N+1][M+1];
  int B[N+1];

#pragma scop
  for (int i = 1; i <= N; ++i) {
    for (int j = 1; j <= M; ++j) {
S1:   A[i][j] = A[i][j-1] + A[i+1][j-1] + B[i-1];
    }
S2: B[i] = B[i-1];
  }
#pragma endscop
}
