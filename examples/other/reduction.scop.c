int main()
{
#if 0
#define N 4000
#else
  int N;
#endif

  int A[N];

#pragma scop
  for (int i = 0; i < N; i++){
    for (int j = 0; j < i+1; j++) {
S1:   A[i] = A[i] + A[i-j];
    }
  }
#pragma endscop
}
