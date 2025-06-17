int main()
{
#if 0
# define T 500
# define N 2000
#else
  int T;
  int N;
#endif

  double a[N][N];

#pragma scop
  for (int t = 0; t <= T-1; t++)  {
    for (int i = 1; i <= N-2; i++)  {
      for (int j = 1; j <= N-2; j++)  {
        a[i][j] = (a[i-1][j-1] + a[i-1][j] + a[i-1][j+1] 
                 + a[i][j-1] + a[i][j] + a[i][j+1]
                 + a[i+1][j-1] + a[i+1][j] + a[i+1][j+1])/9.0;
      }
    }
  }
#pragma endscop
}
