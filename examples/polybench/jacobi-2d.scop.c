double SCALAR_VAL(double);

int main()
{
#if 0
# define TSTEPS 1000
# define N 1000
#else
  int TSTEPS;
  int N;
#endif

  int t,i,j;

  double A[N][N];
  double B[N][N];

#pragma scop
  for (t = 0; t < TSTEPS; t++) {
    for (i = 1; i < N - 1; i++)
      for (j = 1; j < N - 1; j++)
S1:     B[i][j] = SCALAR_VAL(0.2) * (A[i][j] + A[i][j-1] + A[i][1+j] + A[1+i][j] + A[i-1][j]);
        
    for (i = 1; i < N - 1; i++)
      for (j = 1; j < N - 1; j++)
S2:     A[i][j] = SCALAR_VAL(0.2) * (B[i][j] + B[i][j-1] + B[i][1+j] + B[1+i][j] + B[i-1][j]);
  }
#pragma endscop
}

