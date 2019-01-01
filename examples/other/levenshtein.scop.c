#if 1
int MIN(int, int);
#else
#define MIN(a,b) (((a)<(b))?(a):(b))
#endif

int main()
{
#if 0
#define M 4000
#define N 4000
#else
  int M;
  int N;
#endif

  int tab[M+1][N+1];
  char a[M], b[N];

#pragma scop
  for (int i = 1; i <= M; ++i) {
    for (int j = 1; j <= N; ++j) {
S1:   tab[i][j] = MIN(MIN(tab[i - 1][j] + 1, tab[i][j - 1] + 1), tab[i - 1][j - 1] + ((a[i - 1] == b[j - 1]) ? 0 : 1));
    }
  }
#pragma endscop
}
