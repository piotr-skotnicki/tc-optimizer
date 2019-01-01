int main()
{
#if 0
#define N 4000
#else
  int N;
#endif

  int S[N][N];

#pragma scop
  for (int i = -N+1; i <= 0; i++) {
    for (int j = -i+1; j < N; j++) {
      for (int k = 0; k < j+i; k++) {
S1:     S[-i][j] = (S[-i][k-i] + S[k-i+1][j] + S[-i][j]) % 100000;
      }
S2:   S[-i][j] = (S[-i][j]+ S[-i+1][j-1]) % 100000;
    }
  }  
#pragma endscop
}
