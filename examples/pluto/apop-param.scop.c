double max(double, double);

int main()
{
#if 0
# define _ns 2000000
# define _nt 10000
#else
  int _ns;
  int _nt;
#endif

  double F[2][_ns+1];
  double C[3][_ns+1];
  double E;
  double dS;

#pragma scop
  for (int t = 0; t < _nt; ++t) {
    for (int i = 1; i < _ns; ++i) {
S1:   F[(t+1)%2][i] = max(C[0][i] * F[t%2][i-1] + C[1][i] * F[t%2][i] + C[2][i] * F[t%2][i+1], E - i * dS);
    }
  }
#pragma endscop
}
