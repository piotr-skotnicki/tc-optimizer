int main()
{
#if 0
#define N 4000
#else
  int N;
#endif

  int A[N+1][N+1];
  int B[N+1];

#pragma scop
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
S1:   A[i][j+1] = A[i][j] + A[i+1][j];
    }
S2: B[i+1] = B[i] * A[i+1][N]; 
  }
#pragma endscop
}
