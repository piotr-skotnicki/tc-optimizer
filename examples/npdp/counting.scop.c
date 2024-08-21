int paired(int, int);

int main()
{
#if 0
# define N 100
#else
  int N;
#endif

  int ck[N+10][N+10];
  int l;

#pragma scop
  for (int i = N-2; i >= 1; i--) {
    for (int j = i+2; j <= N; j++) {
      for (int k = i; k <= j-l; k++) {
S1:     ck[i][j] +=  ck[i][j-1] + paired(k, j) ? ck[i][k-1] + ck[k+1][j-1] : 0;
      }
    }
  }
#pragma endscop
}
