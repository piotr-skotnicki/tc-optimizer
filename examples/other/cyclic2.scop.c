int main()
{
#if 0
#define N 4000
#else
  int N;
#endif

  int A[N][N];
  int B[N];

#pragma scop
  for (int i = 0; i < N-1; ++i) {
    for (int j = 0; j < N-1; ++j) {
S1:   A[i][j+1] = (A[i][j] + A[i+1][j]) % 1000007;
    }
S2: B[i+1] = (B[i] + A[i+1][N-1]) % 1000007; 
  }
#pragma endscop
}
