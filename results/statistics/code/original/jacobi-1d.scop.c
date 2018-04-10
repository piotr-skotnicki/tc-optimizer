int main()
{
#if 0
# define _PB_TSTEPS 1000
# define _PB_N 1000
#else
  int _PB_TSTEPS;
  int _PB_N;
#endif

  int t,i;
  
  double A[_PB_N];
  double B[_PB_N];

#pragma scop
  for (t = 0; t < _PB_TSTEPS; t++) {
    for (i = 1; i < _PB_N - 1; i++) {
S1:   B[i] = 0.33333 * (A[i-1] + A[i] + A[i + 1]);
    }
    for (i = 1; i < _PB_N - 1; i++) {
S2:   A[i] = 0.33333 * (B[i-1] + B[i] + B[i + 1]);
    }
  }
#pragma endscop
}

