int main()
{
#if 0
#define M 4000
#define N 4000
#else
  int M;
  int N;
#endif

  int A[M+2][N+2];

#pragma scop
  for (int i = 1; i <= M; ++i) {
    for (int j = 1; j <= N; ++j) {
S1:   A[i][j] = A[i][j+1] + A[i+1][j] + A[i+1][j-1];
    }
  }
#pragma endscop
}
