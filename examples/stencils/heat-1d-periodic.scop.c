int main()
{
#if 0
# define N 3600000
# define T 1000
#else
  int T;
  int N;
#endif
  
  int t,i;
  
  int A[2][N];
  
#pragma scop
  for (t = 0; t < T; ++t) {
    for (i = 0; i < N; ++i) {
S1:   A[(t+1)%2][i] = 0.125 * (A[t%2][(N+i+1)%N] - 2.0 * A[t%2][i] + A[t%2][(N+i-1)%N]);
    }
  }
#pragma endscop
}

