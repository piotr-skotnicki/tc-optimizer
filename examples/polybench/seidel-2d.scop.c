double SCALAR_VAL(double);

int main()
{
#if 0
# define _PB_TSTEPS 20
# define _PB_N 40
#else
  int _PB_TSTEPS;
  int _PB_N;
#endif

  int t,i,j;
  
  double A[_PB_N][_PB_N];

#pragma scop
  for (t = 0; t <= _PB_TSTEPS - 1; t++)
    for (i = 1; i <= _PB_N - 2; i++)
      for (j = 1; j <= _PB_N - 2; j++)
S1:     A[i][j] = (A[i-1][j-1] + A[i-1][j] + A[i-1][j+1]
                + A[i][j-1] + A[i][j] + A[i][j+1]
                + A[i+1][j-1] + A[i+1][j] + A[i+1][j+1])/SCALAR_VAL(9.0);
#pragma endscop
}

