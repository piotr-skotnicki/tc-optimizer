int main()
{
#if 0
# define _PB_M 2000
# define _PB_N 2600
#else
  int _PB_M;
  int _PB_N;
#endif
  
  double temp2;
  double alpha;
  double beta;
  double C[_PB_M][_PB_N];
  double A[_PB_M][_PB_M];
  double B[_PB_M][_PB_N];

#pragma scop
   for (int i = 0; i < _PB_M; i++) {
     for (int j = 0; j < _PB_N; j++) {
S1:    temp2 = 0;
       for (int k = 0; k < i; k++) {
S2:      C[k][j] += alpha * B[i][j] * A[i][k];
S3:      temp2 += B[k][j] * A[i][k];
       }
S4:    C[i][j] = beta * C[i][j] + alpha*B[i][j] * A[i][i] + alpha * temp2;
     }
   } 
#pragma endscop
}
