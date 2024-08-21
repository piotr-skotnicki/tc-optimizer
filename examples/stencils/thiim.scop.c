int main()
{
#if 0
# define N 150
# define T 100
#else
  int T;
  int N;
#endif

  int E[N];
  int H[N];
  int k1;
  int k2;

#pragma scop
  for (int t = 0; t < T; ++t) {
    for (int i = 1; i < N-1; ++i) {
S1:   E[i] = k1 * E[i] + k2 * (H[i] - H[i-1]);
    }
    for (int i = 1; i < N-1; ++i) {
S2:   H[i] += E[i] - E[i+1];
    }
  }
#pragma endscop
}
