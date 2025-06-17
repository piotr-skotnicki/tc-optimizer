int main()
{
#if 0
# define N 1000000
# define T 10000
#else
  short N;
  short T;
#endif

  double h[N];
  double e[N+1];
# define coeff1 0.5
# define coeff2 0.7

#pragma scop
  for (int c1 = 1; c1 <= T; c1 += 1) {
    for (int c2 = c1 + 1; c2 <= N + c1; c2 += 1) {
      if (N + c1 >= c2 + 1) {
S1:     e[-c1 + c2] = e[-c1 + c2] - coeff1 * (h[-c1 + c2] - h[-c1 + c2 - 1]);
      }
S2:   h[-c1 + c2 - 1] = h[-c1 + c2 - 1] - coeff2 * (e[-c1 + c2 - 1 + 1] - e[-c1 + c2 - 1]);
    }
  }
#pragma endscop

}
