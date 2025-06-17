double max(double, double);

int main()
{
#if 0
# define ns 2000000
# define nt 10000
#else
  int ns;
  int nt;
#endif

  double F[nt+1][ns+1];
  double C[3][ns+1];
  double E;
  double dS;

#pragma scop
  for (int t = 0; t < nt; ++t) {
    for (int i = 1; i < ns; ++i) {
S1:   F[t+1][i] = max(C[0][i] * F[t][i-1] + C[1][i] * F[t][i] + C[2][i] * F[t][i+1], E - i * dS);
    }
  }
#pragma endscop
}
