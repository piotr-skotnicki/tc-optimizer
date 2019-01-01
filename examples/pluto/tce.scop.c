int main()
{
#if 0
# define N 1024
#else
  int N;
#endif
  
  double A[N][N][N][N];
  double T1[N][N][N][N];
  double T2[N][N][N][N];
  double T3[N][N][N][N];
  double B[N][N][N][N];
  double C1[N][N];
  double C2[N][N];
  double C3[N][N];
  double C4[N][N];

#pragma scop
  for(int a=0; a<N; a++)
    for(int q=0; q<N; q++)
      for(int r=0; r<N; r++)
        for(int s=0; s<N; s++)
          for(int p=0; p<N; p++)
S1:         T1[a][q][r][s] = T1[a][q][r][s] + A[p][q][r][s]*C4[p][a];

  for(int a=0; a<N; a++)
    for(int b=0; b<N; b++)
      for(int r=0; r<N; r++)
        for(int s=0; s<N; s++)
          for(int q=0; q<N; q++)
S2:         T2[a][b][r][s] = T2[a][b][r][s] + T1[a][q][r][s]*C3[q][b];

  for(int a=0; a<N; a++)
    for(int b=0; b<N; b++)
      for(int c=0; c<N; c++)
        for(int s=0; s<N; s++)
          for(int r=0; r<N; r++)
S3:         T3[a][b][c][s] = T3[a][b][c][s] + T2[a][b][r][s]*C2[r][c];

  for(int a=0; a<N; a++)
    for(int b=0; b<N; b++)
      for(int c=0; c<N; c++)
        for(int d=0; d<N; d++)
          for(int s=0; s<N; s++)
S4:         B[a][b][c][d] = B[a][b][c][d] + T3[a][b][c][s]*C1[s][d];
#pragma endscop
}
