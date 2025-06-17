double SCALAR_VAL(double);

int main()
{
#if 0
# define _PB_TMAX 1000
# define _PB_NX 1000
# define _PB_NY 1000
#else
  short _PB_TMAX;
  short _PB_NX;
  short _PB_NY;
#endif
  
  double ex[_PB_NX][_PB_NY];
  double ey[_PB_NX][_PB_NY];
  double hz[_PB_NX][_PB_NY];
  double _fict_[_PB_TMAX];

#pragma scop
  for (int t = 0; t < _PB_TMAX; t++) {
    for (int ii = 0; ii < 4; ii++) {
      for (int i = (ii==1 ? 1 : 0); i < (ii==0 ? 1 : (ii==3 ? _PB_NX - 1 : _PB_NX)); i++) {
        for (int j = (ii==2 ? 1 : 0 ); j < (ii==3 ? (_PB_NY - 1) : _PB_NY ); j++) {
          if (ii == 0) {
S1:         ey[0][j] = _fict_[t];
          }
          if (ii == 1) {
S2:         ey[i][j] = ey[i][j] - SCALAR_VAL(0.5) * (hz[i][j]-hz[i-1][j]);
          }
          if (ii == 2) {
S3:         ex[i][j] = ex[i][j] - SCALAR_VAL(0.5) * (hz[i][j]-hz[i][j-1]);
          }
          if (ii == 3) {
S4:         hz[i][j] = hz[i][j] - SCALAR_VAL(0.7) * (ex[i][j+1] - ex[i][j] + ey[i+1][j] - ey[i][j]);
          }
        }
      }
    }
  }
#pragma endscop
}
