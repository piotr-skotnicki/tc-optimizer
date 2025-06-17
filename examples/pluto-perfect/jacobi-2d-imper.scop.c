int main()
{
#if 0
# define N 2000
# define T 1000
#else
  short N;
  short T;
#endif

  double a[N][N];
  double b[N][N];

#pragma scop
  for (int k0 = 0; k0 < T; k0 += 1) {
    for (int k1 = 2 * k0; k1 < N + 2 * k0 - 2; k1 += 1) {
      for (int k2 = max(2 * k0, -N + k1 + 4); k2 < min(N + 2 * k0 - 2, N + k1 - 3); k2 += 1) {
        if (N + 2 * k0 >= k1 + 4 && N + 2 * k0 >= k2 + 4) {
S1:       b[-2 * k0 + k1 + 2][-2 * k0 + k2 + 2] = (0.2 * ((((a[-2 * k0 + k1 + 2][-2 * k0 + k2 + 2] + a[-2 * k0 + k1 + 2][-2 * k0 + k2 + 1]) + a[-2 * k0 + k1 + 2][-2 * k0 + k2 + 3]) + a[-2 * k0 + k1 + 3][-2 * k0 + k2 + 2]) + a[-2 * k0 + k1 + 1][-2 * k0 + k2 + 2]));
        }
        if (k1 >= 2 * k0 + 1 && k2 >= 2 * k0 + 1) {
S2:       a[-2 * k0 + k1 + 1][-2 * k0 + k2 + 1] = b[-2 * k0 + k1 + 1][-2 * k0 + k2 + 1];
        }
      }
    }
  }
#pragma endscop
}
