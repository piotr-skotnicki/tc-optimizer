int main()
{
#if 0
# define N 1000
#else
  int N;
#endif
  
  int i,j;
  int alpha;
  int beta;
  int A[N][N];
  int u1[N];
  int v1[N];
  int u2[N];
  int v2[N];
  int w[N];
  int x[N];
  int y[N];
  int z[N];

#pragma scop

  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
S1:   A[i][j] = A[i][j] + u1[i] * v1[j] + u2[i] * v2[j];

  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
S2:   x[i] = x[i] + beta * A[j][i] * y[j];

  for (i = 0; i < N; i++)
S3: x[i] = x[i] + z[i];

  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
S4:   w[i] = w[i] + alpha * A[i][j] * x[j];

#pragma endscop
}

