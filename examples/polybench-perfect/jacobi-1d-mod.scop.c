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
  for (int t = 0; t < _PB_TSTEPS * 2; t++) {
    for (int i = 1; i < _PB_N-1; i++) {
      if (t % 2 == 0) {
S1:     B[i] = 0.33333 * (A[i-1] + A[i] + A[i+1]);
      }
      if (t % 2 == 1) {
S2:     A[i] = 0.33333 * (B[i-1] + B[i] + B[i+1]);
      }
    }
  }
#pragma endscop
}
