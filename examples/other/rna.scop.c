#if 1
int MAX(int, int);
int can_pair(const int*, int, int);
#else
#define MAX(s1, s2) (((s1) >= (s2)) ? (s1) : (s2))
#define can_pair(RNA, i, j) (((RNA[i] + RNA[j]) == 3 && (i) < (j) - 1) ? 1 : 0)
#endif

int main()
{
#if 0
#define N 4000
#else
  int N;
#endif

  int i,j,k;
  int S[N][N];
  int RNA[N];

#pragma scop
  for (i = N-1; i >= 0; i--) {
    for (j = i+1; j < N; j++) {
      for (k = 0; k < j-i; k++) {
S1:     S[i][j] = MAX(S[i][k+i] + S[k+i+1][j], S[i][j]);
      }
S2:   S[i][j] = MAX(S[i][j], S[i+1][j-1] + can_pair(RNA, i, j));
    }
  }
#pragma endscop
}

