int main()
{
#if 0
# define N 1000
#else
  int N;
#endif

  int i,j;
  
  double x1[N];
  double x2[N];
  double y_1[N];
  double y_2[N];
  double A[N][N];

#pragma scop
  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
S1:   x1[i] = x1[i] + A[i][j] * y_1[j];
      
  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
S2:   x2[i] = x2[i] + A[j][i] * y_2[j];
#pragma endscop
}

