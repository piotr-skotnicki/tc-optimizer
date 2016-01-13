int main()
{
#if 0
# define TSTEPS 1000
# define N 1000
#else
  int TSTEPS;
  int N;
#endif

  int t,i;
  
  double A[N];
  double B[N];

#pragma scop
  for (t = 0; t < TSTEPS; t++) {
    for (i = 1; i < N - 1; i++) {
S1:   B[i] = 0.33333 * (A[i-1] + A[i] + A[i + 1]);
    }
    for (i = 1; i < N - 1; i++) {
S2:   A[i] = 0.33333 * (B[i-1] + B[i] + B[i + 1]);
    }
  }
#pragma endscop
}

