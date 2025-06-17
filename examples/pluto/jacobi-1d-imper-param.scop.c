int main()
{
#if 0
# define N 2000000
# define T 1000
#else
  int N;
  int T;
#endif

  double a[N];
  double b[N];

  int _T = T;
  int _N = N;
#pragma scop
  for (int t = 0; t < _T; t++) {
    for (int i = 2; i < _N - 1; i++) {
S1:   b[i] = 0.33333 * (a[i-1] + a[i] + a[i + 1]);
    }
    for (int j = 2; j < _N - 1; j++) {
S2:   a[j] = b[j];
    }
  }
#pragma endscop
}
