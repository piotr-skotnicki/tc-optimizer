int MAX(int, int);

int main()
{
#if 0
#define N 4000
#else
  int N;
#endif

  int H[N+1][N+1];

#pragma scop
  for (int i=N-3; i>=0; i--){
    for (int j=i+2; j<=N-1; j++) {
S1:   H[i][j] = H[i+1][j-1];
    }
    for (int j=i; j<=N-2; j++) {
      for (int k=j+1; k<=N-1; k++) {
S2:     H[i][k] = MAX(H[i][k], H[k][j] + H[j+1][k]);
      }
    }  
  }
#pragma endscop 
}
