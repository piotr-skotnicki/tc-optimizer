int max_score(int, int);
int match(int, int);

int main()
{
  int i,j,k;
  int N;
  int table[N][N];
  int seq[N];

#pragma scop
  for (i = N-1; i >= 0; i--) {
    for (j = i+1; j < N; j++) {
      if (j-1 >= 0) 
S1:     table[i][j] = max_score(table[i][j], table[i][j-1]);
      if (i+1 < N) 
S2:     table[i][j] = max_score(table[i][j], table[i+1][j]);

      if (j-1 >= 0 && i+1 < N) {
        if (i < j-1) 
S3:       table[i][j] = max_score(table[i][j], table[i+1][j-1]+match(seq[i], seq[j]));
        else 
S4:       table[i][j] = max_score(table[i][j], table[i+1][j-1]);
      }

      for (k = i+1; k < j; k++) {
S5:     table[i][j] = max_score(table[i][j], table[i][k] + table[k+1][j]);
      }
    }
  }
#pragma endscop
}

