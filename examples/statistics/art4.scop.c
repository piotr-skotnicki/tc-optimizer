int main()
{
#if 1
#define N 5
#else
  int N;
#endif

  int B[N+2][N+2];

#pragma scop
  for (int i = 0; i <= 5; ++i) {
    for (int j = 0; j <= 5; ++j) {
S1:   B[i][j] = B[i+1][5-j];
    }
  }
#pragma endscop
}
