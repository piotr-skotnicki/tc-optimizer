int MAX(int, int);
int MIN(int, int);
int s(unsigned char, unsigned char);

#define min(a,b) (((a)<(b))?(a):(b))

int main()
{
#if 0
# define N 100
#else
  int N;
#endif

  int m1[N+2][N+2][N+2];
  int m2[N+2][N+2][N+2];
  int m3[N+2][N+2][N+2];
  int m4[N+2][N+2][N+2];
  int m5[N+2][N+2][N+2];
  int m6[N+2][N+2][N+2];
  int H[N+2][N+2][N+2];
  int W[N+2];
  unsigned char a[N+2];
  unsigned char b[N+2];
  unsigned char c[N+2];
  int INT_MIN;

#pragma scop
  for (int i = 1; i <= N; i++) {
    for (int j = 1; j <= N; j++) {
      for (int l = 1; l <= N; l++) {
S1:     m1[i][j][l] = INT_MIN;
        for (int k = 1; k <= i; k++) {
S2:       m1[i][j][l] = MAX(m1[i][j][l] ,H[i-k][j][l] - 2*W[k]);
        }
S3:     m2[i][j][l] = INT_MIN;
        for (int k = 1; k <= j; k++) {
S4:       m2[i][j][l] = MAX(m2[i][j][l], H[i][j-k][l] - 2*W[k]);
        }
S5:     m3[i][j][l] = INT_MIN;
        for (int k = 1; k <= l; k++) {
S6:       m3[i][j][l] = MAX(m3[i][j][l], H[i][j][l-k] - 2*W[k]);
        }
S7:     m4[i][j][l] = INT_MIN;
        for (int k = 1; k <= min(i, j); k++) {
S8:       m4[i][j][l] = MAX(m4[i][j][l], H[i-k][j-k][l] - W[k] + s(a[i], b[j]));
        }
S9:     m5[i][j][l] = INT_MIN;
        for (int k = 1; k <= min(j, l); k++) {
S10:      m5[i][j][l] = MAX(m5[i][j][l], H[i][j-k][l-k] - W[k] + s(b[j], c[l]));
        }
S11:    m6[i][j][l] = INT_MIN;
        for (int k = 1; k <= min(i, l); k++) {
S12:      m6[i][j][l] = MAX(m6[i][j][l], H[i-k][j][l-k] - W[k] + s(a[i], c[l]));
        }
S13:    H[i][j][l] = MAX(0, MAX(H[i-1][j-1][l-1] + s(a[i], b[j]) + s(a[i], c[l]) + s(b[j], c[l]), MAX(m1[i][j][l], MAX(m2[i][j][l], MAX(m3[i][j][l], MAX(m4[i][j][l], MAX(m5[i][j][l], m6[i][j][l])))))));
      }
    }
  }
#pragma endscop
}
