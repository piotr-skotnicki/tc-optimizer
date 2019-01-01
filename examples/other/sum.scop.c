int main()
{
#if 0
#define n 4000
#else
  int n;
#endif

  int sum;
  int a[n+1][n+1],b[n+1][n+1],c[n+1];

#pragma scop
  for (int i = 0; i <= n; i++) {
S1: sum = 0;
    for (int j = 0; j <= n; j++) {
S2:   sum = sum + a[i][j] * b[i][j];
    }
S3: c[i] = sum;
  }
#pragma endscop
}
