int main()
{
#if 0
#define T 500
#define N 500
#else
  int T;
  int N;
#endif
  
  double X[N][N+13];
  double A[N][N+23];
  double B[N][N+37];

#pragma scop
  for (int t = 0; t < T; t++) {
    for (int i1 = 0; i1 < N; i1++) {
      for (int i2 = 1; i2 < N; i2++) {
S1:     X[i1][i2] = X[i1][i2] - X[i1][i2-1] * A[i1][i2] / B[i1][i2-1];
S2:     B[i1][i2] = B[i1][i2] - A[i1][i2] * A[i1][i2] / B[i1][i2-1];
      }
    }

    for (int i1 = 1; i1 < N; i1++) {
      for (int i2 = 0; i2 < N; i2++) {
S3:     X[i1][i2] = X[i1][i2] - X[i1-1][i2] * A[i1][i2] / B[i1-1][i2];
S4:     B[i1][i2] = B[i1][i2] - A[i1][i2] * A[i1][i2] / B[i1-1][i2];
      }
    }
  }
#pragma endscop
}
