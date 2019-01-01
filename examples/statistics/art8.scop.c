int main()
{
#if 0
#define N 5
#else
  int N;
#endif

  int a[N+1][N+1];

#pragma scop
  for (int i = 1; i <= N; ++i) {
    for (int j = 1; j <= N; ++j) {
      if (j == N)
S1:     a[i+1][1] = a[i][N] + a[i][j-1];
      else 
S2:     a[i][j] = a[i][j-1];
    }
  }
#pragma endscop
}

//S1[1, 6]
//[N] -> { S2[1, 6] : N = 6; S2[2, 2] : N = 6; S2[2, 1] : N = 6 }
