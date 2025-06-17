int main()
{
#if 0
# define _PB_TSTEPS 1000
# define _PB_N 1000
#else
  short _PB_TSTEPS;
  short _PB_N;
#endif
  
  double A[_PB_N];
  double B[_PB_N];

#pragma scop
  for (int c0 = 0; c0 < _PB_TSTEPS; c0 += 1) {
    for (int c1 = 2 * c0 + 1; c1 < _PB_N + 2 * c0; c1 += 1) {
      if (_PB_N + 2 * c0 >= c1 + 2) {
S1:     B[-2 * c0 + c1] = (0.33333 * ((A[-2 * c0 + c1 - 1] + A[-2 * c0 + c1]) + A[-2 * c0 + c1 + 1]));
      }
      if (c1 >= 2 * c0 + 2) {
S2:     A[-2 * c0 + c1 - 1] = (0.33333 * ((B[-2 * c0 + c1 - 2] + B[-2 * c0 + c1 - 1]) + B[-2 * c0 + c1]));
      }
    }
  }
#pragma endscop
}
