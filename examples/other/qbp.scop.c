int paired(int, int);

int main()
{
#if 0
#define N 4000
#else
  int N;
#endif

  int l,ERT;
  int Q[N][N];
  int Qbp[N][N];

#pragma scop
  if(N>=1 && l>=0 && l<=5) {
    for(int i=N-1; i>=0; i--) {
      for(int j=i+1; j<N; j++) {
S1:     Q[i][j] = Q[i][j-1];
        for(int k=0; k<j-i-l; k++) {
S2:       Qbp[k+i][j] = Q[k+i+1][j-1] * ERT * paired(k+i,j-1);
S3:       Q[i][j] += Q[i][k+i] * Qbp[k+i][j];
        }
      }
    }
  }
#pragma endscop
}
