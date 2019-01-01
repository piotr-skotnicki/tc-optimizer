int main()
{
#if 0
#define N 4000
#else
  int N;
#endif

  int S[N][N];

#pragma scop
  for (int inew = 1; inew <= N; inew++) {
    for (int j = N-inew+1; j < N; j++) {
      for (int k = 0; k < j-(N-inew); k++) {
S1:     S[N-inew][j] = (S[N-inew][k+N-inew] + S[k+N-inew+1][j] + S[N-inew][j]) % 100000; 
      }
S2:   S[N-inew][j] = (S[N-inew][j]+ S[N-inew+1][j-1]) % 100000;
    }
  }
#pragma endscop
}
