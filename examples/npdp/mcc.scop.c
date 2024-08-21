int paired(int, int);

int main()
{
#if 0
# define N 100
#else
  int N;
#endif

  double Q1[N+10][N+10];
  double Qbp1[N+10][N+10];
  float ERT;
  int l;

#pragma scop
  if (N >= 1 && l >= 0 && l <= 5) {
    for (int i = N-1; i >= 0; i--) {
      for (int j = i+1; j < N; j++) {
S1:     Q1[i][j] = Q1[i][j-1];
        for(int k = 0; k < j-i-l; k++) {
S2:       Qbp1[k+i][j] = Q1[k+i+1][j-1] * ERT * paired(k+i, j-1);
S3:       Q1[i][j] += Q1[i][k+i] * Qbp1[k+i][j];
        }
      }
    }
  }
#pragma endscop
}
