int main()
{
#if 0
# define nx 2048
# define ny 2048
# define tmax 128
#else
  int nx;
  int ny;
  int tmax;
#endif

  double ex[nx][ny+1];
  double ey[nx+1][ny];
  double hz[nx][ny];

#pragma scop
  for (int t = 0; t < tmax; t++)  {
    for (int j = 0; j < ny; j++)
S1:   ey[0][j] = t;
    for (int i = 1; i < nx; i++)
      for (int j = 0; j < ny; j++)
S2:         ey[i][j] = ey[i][j] - 0.5 * (hz[i][j] - hz[i-1][j]);
    for (int i = 0; i < nx; i++)
      for (int j = 1; j < ny; j++)
S3:         ex[i][j] = ex[i][j] - 0.5 * (hz[i][j] - hz[i][j-1]);
    for (int i = 0; i < nx; i++)
      for (int j = 0; j < ny; j++)
S4:     hz[i][j] = hz[i][j] - 0.7 * (ex[i][j+1] - ex[i][j] + ey[i+1][j] - ey[i][j]);
  }
#pragma endscop
}
