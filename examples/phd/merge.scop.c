int main()
{
#if 0
#define N 9
#else
  int N;
#endif

  int A[N+2][N+2];
  int B[N+2];

#pragma scop
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j <= i; ++j) {
S1:   A[i][j] = A[i][j+1];
    }
S2: B[i] = B[i+1] + A[i][i] + A[i+1][i+1];
  }
#pragma endscop

}
