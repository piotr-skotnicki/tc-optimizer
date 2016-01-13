int main()
{
#if 0
# define N 1000
#else
  int N;
#endif

  int i,j,k;
  
  double A[N][N];
  double b[N];
  double x[N];
  double y[N];
  double w;

#pragma scop
  for (i = 0; i < N; i++) {
    for (j = 0; j < i; j++) {
S1:   w = A[i][j];
      for (k = 0; k < j; k++) {
S2:     w -= A[i][k] * A[k][j];
      }
S3:   A[i][j] = w / A[j][j];
    }
    for (j = i; j < N; j++) {
S4:   w = A[i][j];
      for (k = 0; k < i; k++) {
S5:     w -= A[i][k] * A[k][j];
      }
S6:   A[i][j] = w;
    }
  }

  for (i = 0; i < N; i++) {
S7: w = b[i];
    for (j = 0; j < i; j++)
S8:   w -= A[i][j] * y[j];
S9: y[i] = w;
  }

  for (i = N-1; i >= 0; i--) {
S10:w = y[i];
    for (j = i+1; j < N; j++) {
S11:  w -= A[i][j] * x[j];
    }
S12:x[i] = w / A[i][i];
  }
#pragma endscop
}

