int main()
{
#if 0
# define N 1000000
# define T 10000
#else
  int N;
  int T;
#endif

  double h[N];
  double e[N+1];
# define coeff1 0.5
# define coeff2 0.7

  int _T = T;
  int _N = N;
#pragma scop
  for (int t = 1; t <= _T; t++) {
    for (int i = 1; i <= _N-1; i++)
S1:   e[i] = e[i] - coeff1 * (h[i] - h[i-1]);
    for (int i = 0; i <= _N-1; i++)
S2:   h[i] = h[i] - coeff2 * (e[i+1] - e[i]);
  }
#pragma endscop
}
