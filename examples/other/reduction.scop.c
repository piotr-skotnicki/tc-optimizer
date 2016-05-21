int main()
{
#if 0
#define N 4000
#else
  int N;
#endif

  int i,j;
  int A[N];

#pragma scop
  for (i = 0; i < N; i++){
    for (j = 0; j < i+1; j++) {
S1:   A[i] = A[i] + A[i-j];
    }
  }
#pragma endscop
}

