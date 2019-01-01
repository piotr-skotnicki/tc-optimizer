double SCALAR_VAL(double);

int main()
{
#if 0
# define _PB_TSTEPS 1000
# define _PB_N 1000
#else
  int _PB_TSTEPS;
  int _PB_N;
#endif
  
  double u[_PB_N][_PB_N];
  double v[_PB_N][_PB_N];
  double p[_PB_N][_PB_N];
  double q[_PB_N][_PB_N];
  double DX, DY, DT;
  double B1, B2;
  double mul1, mul2;
  double a, b, c, d, e, f;

#pragma scop
  for (int t = 1; t <= _PB_TSTEPS; t++) {
    for (int i = 1; i < _PB_N-1; i++) {
S1:   v[0][i] = SCALAR_VAL(1.0);
S2:   p[i][0] = SCALAR_VAL(0.0);
S3:   q[i][0] = v[0][i];
      for (int j = 1; j < _PB_N-1; j++) {
S4:     p[i][j] = -c / (a*p[i][j-1]+b);
S5:     q[i][j] = (-d*u[j][i-1]+(SCALAR_VAL(1.0)+SCALAR_VAL(2.0)*d)*u[j][i] - f*u[j][i+1]-a*q[i][j-1])/(a*p[i][j-1]+b);
      }
S6:   v[_PB_N-1][i] = SCALAR_VAL(1.0);
      for (int j = _PB_N-2; j >= 1; j--) {
S7:     v[j][i] = p[i][j] * v[j+1][i] + q[i][j];
      }
    }
    
    for (int i = 1; i < _PB_N-1; i++) {
S8:   u[i][0] = SCALAR_VAL(1.0);
S9:   p[i][0] = SCALAR_VAL(0.0);
S10:  q[i][0] = u[i][0];
      for (int j = 1; j < _PB_N-1; j++) {
S11:    p[i][j] = -f / (d*p[i][j-1]+e);
S12:    q[i][j] = (-a*v[i-1][j]+(SCALAR_VAL(1.0)+SCALAR_VAL(2.0)*a)*v[i][j] - c*v[i+1][j]-d*q[i][j-1])/(d*p[i][j-1]+e);
      }
S13:  u[i][_PB_N-1] = SCALAR_VAL(1.0);
      for (int j = _PB_N-2; j >= 1; j--) {
S14:    u[i][j] = p[i][j] * u[i][j+1] + q[i][j];
      }
    }
  }
#pragma endscop
}
