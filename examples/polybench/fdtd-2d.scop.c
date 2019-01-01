double SCALAR_VAL(double);

int main()
{
#if 0
# define _PB_TMAX 1000
# define _PB_NX 1000
# define _PB_NY 1000
#else
  int _PB_TMAX;
  int _PB_NX;
  int _PB_NY;
#endif
  
  double ex[_PB_NX][_PB_NY];
  double ey[_PB_NX][_PB_NY];
  double hz[_PB_NX][_PB_NY];
  double _fict_[_PB_TMAX];

#pragma scop
  for (int t = 0; t < _PB_TMAX; t++) {
    for (int j = 0; j < _PB_NY; j++)
S1:   ey[0][j] = _fict_[t];
      
    for (int i = 1; i < _PB_NX; i++)
      for (int j = 0; j < _PB_NY; j++)
S2:     ey[i][j] = ey[i][j] - SCALAR_VAL(0.5) * (hz[i][j]-hz[i-1][j]);
        
    for (int i = 0; i < _PB_NX; i++)
      for (int j = 1; j < _PB_NY; j++)
S3:     ex[i][j] = ex[i][j] - SCALAR_VAL(0.5) * (hz[i][j]-hz[i][j-1]);
    
    for (int i = 0; i < _PB_NX - 1; i++)
      for (int j = 0; j < _PB_NY - 1; j++)
S4:     hz[i][j] = hz[i][j] - SCALAR_VAL(0.7) * (ex[i][j+1] - ex[i][j] + ey[i+1][j] - ey[i][j]);
  }
#pragma endscop
}
