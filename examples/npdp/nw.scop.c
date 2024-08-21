int MAX(int, int);
int s(unsigned char, unsigned char);

int main()
{
#if 0
# define N 100
#else
  int N;
#endif

  int m1[2*N+2][2*N+2];
  int m2[2*N+2][2*N+2];
  int F[2*N+2][2*N+2];
  int W[2*N+2];
  unsigned char a[2*N+2];
  unsigned char b[2*N+2];
  int INT_MIN;

#pragma scop
  for (int i = 1; i <= N; i++) {
    for (int j = 1; j <= N; j++) {
S1:   m1[i][j] = INT_MIN;
      for (int k = 1; k <= i; k++) {
S2:     m1[i][j] = MAX(m1[i][j], F[i-k][j] - W[k]);
      }
S3:   m2[i][j] = INT_MIN;
      for (int k = 1; k <= j; k++) {
S4:     m2[i][j] = MAX(m2[i][j], F[i][j-k] - W[k]);
      }
S5:   F[i][j] = MAX(0, MAX(F[i-1][j-1] + s(a[i], b[i]), MAX(m1[i][j], m2[i][j])));
    }
  }
#pragma endscop
}
