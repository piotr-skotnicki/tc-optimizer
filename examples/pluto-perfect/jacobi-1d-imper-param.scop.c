int main()
{
#if 0
# define N 2000000
# define T 1000
#else
  short N;
  short T;
#endif

  double a[N];
  double b[N];

  short _N = N;
  short _T = T;
#pragma scop
  if (_N >= 4) {
    for (int k0 = 0; k0 < _T; k0 += 1) {
      for (int k1 = 2 * k0; k1 < _N + 2 * k0 - 2; k1 += 1) {
        if (_N + 2 * k0 >= k1 + 4) {
S1:       b[-2 * k0 + k1 + 2] = (0.33333 * ((a[-2 * k0 + k1 + 1] + a[-2 * k0 + k1 + 2]) + a[-2 * k0 + k1 + 3]));
        }
        if (k1 >= 2 * k0 + 1) {
S2:       a[-2 * k0 + k1 + 1] = b[-2 * k0 + k1 + 1];
        }
      }
    }
  }
#pragma endscop
}
