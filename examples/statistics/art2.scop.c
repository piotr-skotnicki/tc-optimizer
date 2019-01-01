int main()
{
#if 0
#define N 3
#else
  int N;
#endif

  int A[N+1][N+1], B[N+1];

#pragma scop
  for (int i = 0; i <= N; ++i) {
    for (int j = 0; j <= N; ++j) {
S1:   A[i][j] = A[i][j+1];
    }
S2: B[i] = A[i][N] + A[i+1][0];
  }
#pragma endscop
}
