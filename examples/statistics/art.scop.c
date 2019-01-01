int main()
{
#if 1
#define N 4
#else
  int N;
#endif

  int A[N+2][N+2];
  
#pragma scop
  for (int i = 1; i <= N; ++i) {
    for (int j = 1; j <= N; ++j) {
S1:   A[i][j] = A[i][j+1] + A[i+1][j] + A[i+1][j-1];
    }
  }
#pragma endscop
}
