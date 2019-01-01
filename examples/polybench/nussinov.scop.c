int max_score(int, int);
int match(int, int);

int main()
{
#if 0
#define _PB_N 4000
#else
  int _PB_N;
#endif

  int table[_PB_N][_PB_N];
  int seq[_PB_N];

#pragma scop
  for (int i = _PB_N-1; i >= 0; i--) {
    for (int j = i+1; j < _PB_N; j++) {
      if (j-1 >= 0) 
S1:     table[i][j] = max_score(table[i][j], table[i][j-1]);
      if (i+1 < _PB_N) 
S2:     table[i][j] = max_score(table[i][j], table[i+1][j]);

      if (j-1 >= 0 && i+1 < _PB_N) {
        if (i < j-1) 
S3:       table[i][j] = max_score(table[i][j], table[i+1][j-1]+match(seq[i], seq[j]));
        else 
S4:       table[i][j] = max_score(table[i][j], table[i+1][j-1]);
      }

      for (int k = i+1; k < j; k++) {
S5:     table[i][j] = max_score(table[i][j], table[i][k] + table[k+1][j]);
      }
    }
  }
#pragma endscop
}
