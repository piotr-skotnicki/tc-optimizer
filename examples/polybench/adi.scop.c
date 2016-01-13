double SCALAR_VAL(double);

int main()
{
#if 0
# define TSTEPS 1000
# define N 1000
#else
  int TSTEPS;
  int N;
#endif

  int t,i,j;
  
  double u[N][N];
  double v[N][N];
  double p[N][N];
  double q[N][N];
  double DX, DY, DT;
  double B1, B2;
  double mul1, mul2;
  double a, b, c, d, e, f;

#pragma scop
  for (t = 1; t <= TSTEPS; t++) {
    for (i = 1; i < N-1; i++) {
S1:   v[0][i] = SCALAR_VAL(1.0);
S2:   p[i][0] = SCALAR_VAL(0.0);
S3:   q[i][0] = v[0][i];
      for (j = 1; j < N-1; j++) {
S4:     p[i][j] = -c / (a*p[i][j-1]+b);
S5:     q[i][j] = (-d*u[j][i-1]+(SCALAR_VAL(1.0)+SCALAR_VAL(2.0)*d)*u[j][i] - f*u[j][i+1]-a*q[i][j-1])/(a*p[i][j-1]+b);
      }
S6:   v[N-1][i] = SCALAR_VAL(1.0);
      for (j = N-2; j >= 1; j--) {
S7:     v[j][i] = p[i][j] * v[j+1][i] + q[i][j];
      }
    }
    
    for (i = 1; i < N-1; i++) {
S8:   u[i][0] = SCALAR_VAL(1.0);
S9:   p[i][0] = SCALAR_VAL(0.0);
S10:  q[i][0] = u[i][0];
      for (j = 1; j < N-1; j++) {
S11:    p[i][j] = -f / (d*p[i][j-1]+e);
S12:    q[i][j] = (-a*v[i-1][j]+(SCALAR_VAL(1.0)+SCALAR_VAL(2.0)*a)*v[i][j] - c*v[i+1][j]-d*q[i][j-1])/(d*p[i][j-1]+e);
      }
S13:  u[i][N-1] = SCALAR_VAL(1.0);
      for (j = N-2; j >= 1; j--) {
S14:    u[i][j] = p[i][j] * u[i][j+1] + q[i][j];
      }
    }
  }
#pragma endscop
}

