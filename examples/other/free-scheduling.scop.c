int main()
{
#if 0
#define N 4000
#else
  int N;
#endif

  int A[N+2][N+2];
  int B[N+2];

#pragma scop
  for (int i = 1; i <= N; ++i) {
S1: B[i] = A[i+1][N] + B[i+1];
    for (int j = 1; j <= N; ++j) {
S2:   A[i][j] = A[i-1][j];
    }
  }
#pragma endscop
}
