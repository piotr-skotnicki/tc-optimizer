int main()
{
#if 0
# define N 1600000
# define T 1000
#else
  int N;
  int T;
#endif
  
  double A[T+1][N];
  short _N = N - 1;
  short _T = T;
  
#pragma scop
  for (int t = 0; t < _T; t++) {
    for (int i = 0; i <= _N; i++) {
S1:   A[t + 1][i] =
          0.125 * (A[t][i == _N ? 0 : i + 1] - 2.0 * A[t][i] +
                   A[t][((i == 0) ? (_N) : (i - 1))]);
    }
  }
#pragma endscop
}
