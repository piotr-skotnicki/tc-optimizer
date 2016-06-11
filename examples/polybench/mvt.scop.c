int main()
{
#if 0
# define _PB_N 1000
#else
  int _PB_N;
#endif

  int i,j;
  
  double x1[_PB_N];
  double x2[_PB_N];
  double y_1[_PB_N];
  double y_2[_PB_N];
  double A[_PB_N][_PB_N];

#pragma scop
  for (i = 0; i < _PB_N; i++)
    for (j = 0; j < _PB_N; j++)
S1:   x1[i] = x1[i] + A[i][j] * y_1[j];
      
  for (i = 0; i < _PB_N; i++)
    for (j = 0; j < _PB_N; j++)
S2:   x2[i] = x2[i] + A[j][i] * y_2[j];
#pragma endscop
}

