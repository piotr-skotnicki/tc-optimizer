int main()
{
#if 0
# define N 3600000
# define T 1000
#else
  int T;
  int N;
#endif

  int A[T+1][N+2];

#pragma scop
  for (int t = 1; t <= T; ++t) {
    for (int i = 0; i < N; ++i) {
S1:   A[t][i] = A[t-1][(N+i-1)%N] + A[t-1][(i+1)%N];
    }
  }
#pragma endscop
}
