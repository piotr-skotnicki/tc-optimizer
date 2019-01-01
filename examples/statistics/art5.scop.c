int main()
{
#if 1
#define loop 256
#define n 257
#else
  int loop;
  int n;
#endif

  int w[n];
  int b[n][n];

#pragma scop
  for (int l = 1; l <= loop; ++l) {
    for (int i = 1; i < n; ++i) {
      for (int k = 0; k < i; ++k) {
S1:     w[i] += b[k][i] * w[(i-k)-1];
      }
    }
  }
#pragma endscop
}
