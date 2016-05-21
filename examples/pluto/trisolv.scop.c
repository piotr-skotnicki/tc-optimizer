int main()
{
#if 0
# define N 1000
#else
  int N;
#endif

  int i,j,k;
  
  double B[N][N];
  double L[N][N];

#pragma scop
  for (i=0; i <= N-1; i++) {
    for (j=0; j <= N-1; j++) {
      for (k=0; k <= j-1; k++) {
S1:     B[j][i] = B[j][i] - L[j][k] * B[k][i];
      }
S2:   B[j][i] = B[j][i] / L[j][j];
    }
  }
#pragma endscop
}

