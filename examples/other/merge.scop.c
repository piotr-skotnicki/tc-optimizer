int main()
{
#if 0
#define N 4000
#else
  int N;
#endif

  int A[N][N];

#pragma scop
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
S1:   A[i] = A[i] + A[j] * 5;
    }
S2: A[i] = A[i] + 1;
  }
#pragma endscop
}
