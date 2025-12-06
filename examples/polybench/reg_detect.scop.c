int main()
{
#if 0
# define _PB_LENGTH 4000
# define _PB_MAXGRID 4000
# define _PB_NITER 4000
#else
  int _PB_LENGTH;
  int _PB_MAXGRID;
  int _PB_NITER;
#endif

  int sum_tang[_PB_MAXGRID][_PB_MAXGRID];
  int mean[_PB_MAXGRID][_PB_MAXGRID];
  int diff[_PB_MAXGRID][_PB_MAXGRID][_PB_LENGTH];
  int sum_diff[_PB_MAXGRID][_PB_MAXGRID][_PB_LENGTH];
  int path[_PB_MAXGRID][_PB_MAXGRID];

#pragma scop
  for (int t = 0; t < _PB_NITER; t++) {
    for (int j = 0; j <= _PB_MAXGRID - 1; j++) {
      for (int i = j; i <= _PB_MAXGRID - 1; i++) {
        for (int cnt = 0; cnt <= _PB_LENGTH - 1; cnt++) {
          diff[j][i][cnt] = sum_tang[j][i];
        }
      }
    }

    for (int j = 0; j <= _PB_MAXGRID - 1; j++) {
      for (int i = j; i <= _PB_MAXGRID - 1; i++) {
        sum_diff[j][i][0] = diff[j][i][0];
        for (int cnt = 1; cnt <= _PB_LENGTH - 1; cnt++) {
          sum_diff[j][i][cnt] = sum_diff[j][i][cnt - 1] + diff[j][i][cnt];
        }
        mean[j][i] = sum_diff[j][i][_PB_LENGTH - 1];
      }
    }

    for (int i = 0; i <= _PB_MAXGRID - 1; i++) {
      path[0][i] = mean[0][i];
    }

    for (int j = 1; j <= _PB_MAXGRID - 1; j++) {
      for (int i = j; i <= _PB_MAXGRID - 1; i++) {
        path[j][i] = path[j - 1][i - 1] + mean[j][i];
      }
    }
  }
#pragma endscop

}
