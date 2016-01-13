double SCALAR_VAL(double);

int main()
{
#if 0
# define N 1000
# define M 1000
#else
  int N;
  int M;
#endif

  int i,j;
  
  double A[N][M];
  double s[M];
  double q[N];
  double p[M];
  double r[N];

#pragma scop
  for (i = 0; i < M; i++) {
S1: s[i] = 0;
  }
    
  for (i = 0; i < N; i++) {
S2: q[i] = SCALAR_VAL(0.0);
    for (j = 0; j < M; j++) {
S3:   s[j] = s[j] + r[i] * A[i][j];
S4:   q[i] = q[i] + A[i][j] * p[j];
    }
  }
#pragma endscop
}

