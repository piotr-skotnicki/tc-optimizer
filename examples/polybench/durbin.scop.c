double SCALAR_VAL(double);

int main()
{
#if 0
# define N 1000
#else
  int N;
#endif

  int i,k;
 
  double z[N];
  double alpha;
  double beta;
  double sum;
  double r[N];
  double y[N];

#pragma scop
  for (k = 1; k < N; k++) {
S1: beta = (1-alpha*alpha)*beta;
S2: sum = SCALAR_VAL(0.0);
    for (i = 0; i < k; i++) {
S3:   sum += r[k-i-1]*y[i];
    }
S4: alpha = -(r[k] + sum)/beta;

    for (i = 0; i < k; i++) {
S5:   z[i] = y[i] + alpha*y[k-i-1];
    }
    for (i = 0; i < k; i++) {
S6:   y[i] = z[i];
    }
S7: y[k] = alpha;
  }
#pragma endscop
}

