int main()
{
#if 0
# define N 1000
#else
  int N;
#endif
  
  int i,j;
  
  double L[N][N];
  double x[N];
  double b[N];

#pragma scop
  for (i = 0; i < N; i++) {
S1: x[i] = b[i];
    for (j = 0; j < i; j++) {
S2:   x[i] -= L[i][j] * x[j];
    }  
S3: x[i] = x[i] / L[i][i];
  }
#pragma endscop
}

