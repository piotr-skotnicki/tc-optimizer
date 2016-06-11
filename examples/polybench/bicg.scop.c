double SCALAR_VAL(double);

int main()
{
#if 0
# define _PB_N 1000
# define _PB_M 1000
#else
  int _PB_N;
  int _PB_M;
#endif

  int i,j;
  
  double A[_PB_N][_PB_M];
  double s[_PB_M];
  double q[_PB_N];
  double p[_PB_M];
  double r[_PB_N];

#pragma scop
  for (i = 0; i < _PB_M; i++) {
S1: s[i] = 0;
  }
    
  for (i = 0; i < _PB_N; i++) {
S2: q[i] = SCALAR_VAL(0.0);
    for (j = 0; j < _PB_M; j++) {
S3:   s[j] = s[j] + r[i] * A[i][j];
S4:   q[i] = q[i] + A[i][j] * p[j];
    }
  }
#pragma endscop
}

