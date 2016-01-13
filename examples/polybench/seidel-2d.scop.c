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

#pragma scop
  for (t = 0; t <= TSTEPS - 1; t++)
    for (i = 1; i<= N - 2; i++)
      for (j = 1; j <= N - 2; j++)
S1:     A[i][j] = (A[i-1][j-1] + A[i-1][j] + A[i-1][j+1]
                + A[i][j-1] + A[i][j] + A[i][j+1]
                + A[i+1][j-1] + A[i+1][j] + A[i+1][j+1])/SCALAR_VAL(9.0);
#pragma endscop
}

