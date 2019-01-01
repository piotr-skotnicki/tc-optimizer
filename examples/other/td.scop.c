int main()
{
#if 0
# define TSTEPS 1000
# define N 1000
#else
  int TSTEPS;
  int N;
#endif
  
  double A[TSTEPS][N];
  double B[TSTEPS][N];

#pragma scop
  for (int t = 0; t < TSTEPS; t++) {
    for (int i = 1; i < N - 1; i++) {
S1:   B[t][i] = 0.33333 * (A[t][i-1] + A[t][i] + A[t][i + 1]);
    }
    for (int i = 1; i < N - 1; i++) {
S2:   A[t][i] = 0.33333 * (B[t][i-1] + B[t][i] + B[t][i + 1]);
    }
  }
#pragma endscop
}
