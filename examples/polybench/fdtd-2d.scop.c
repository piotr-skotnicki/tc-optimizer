double SCALAR_VAL(double);

int main()
{
#if 0
# define TMAX 1000
# define NX 1000
# define NY 1000
#else
  int TMAX;
  int NX;
  int NY;
#endif

  int t,i,j;
  
  double ex[NX][NY];
  double ey[NX][NY];
  double hz[NX][NY];
  double _fict_[TMAX];

#pragma scop
  for(t = 0; t < TMAX; t++) {
    for (j = 0; j < NY; j++)
S1:   ey[0][j] = _fict_[t];
      
    for (i = 1; i < NX; i++)
      for (j = 0; j < NY; j++)
S2:     ey[i][j] = ey[i][j] - SCALAR_VAL(0.5) * (hz[i][j]-hz[i-1][j]);
        
    for (i = 0; i < NX; i++)
      for (j = 1; j < NY; j++)
S3:     ex[i][j] = ex[i][j] - SCALAR_VAL(0.5) * (hz[i][j]-hz[i][j-1]);
    
    for (i = 0; i < NX - 1; i++)
      for (j = 0; j < NY - 1; j++)
S4:     hz[i][j] = hz[i][j] - SCALAR_VAL(0.7) * (ex[i][j+1] - ex[i][j] + ey[i+1][j] - ey[i][j]);
  }
#pragma endscop
}

