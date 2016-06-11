int main()
{
#if 0
# define _PB_N 1000
#else
  int _PB_N;
#endif
  
  int i,j;
  int alpha;
  int beta;
  int A[_PB_N][_PB_N];
  int u1[_PB_N];
  int v1[_PB_N];
  int u2[_PB_N];
  int v2[_PB_N];
  int w[_PB_N];
  int x[_PB_N];
  int y[_PB_N];
  int z[_PB_N];

#pragma scop

  for (i = 0; i < _PB_N; i++)
    for (j = 0; j < _PB_N; j++)
S1:   A[i][j] = A[i][j] + u1[i] * v1[j] + u2[i] * v2[j];

  for (i = 0; i < _PB_N; i++)
    for (j = 0; j < _PB_N; j++)
S2:   x[i] = x[i] + beta * A[j][i] * y[j];

  for (i = 0; i < _PB_N; i++)
S3: x[i] = x[i] + z[i];

  for (i = 0; i < _PB_N; i++)
    for (j = 0; j < _PB_N; j++)
S4:   w[i] = w[i] + alpha * A[i][j] * x[j];

#pragma endscop
}

