double SCALAR_VAL(double);

int main()
{
#if 0
# define _PB_N 4000
#else
  int _PB_N;
#endif
 
  double z[_PB_N];
  double alpha;
  double beta;
  double sum;
  double r[_PB_N];
  double y[_PB_N];

#pragma scop
  for (int k = 1; k < _PB_N; k++) {
S1: beta = (1-alpha*alpha)*beta;
S2: sum = SCALAR_VAL(0.0);
    for (int i = 0; i < k; i++) {
S3:   sum += r[k-i-1]*y[i];
    }
S4: alpha = -(r[k] + sum)/beta;

    for (int i = 0; i < k; i++) {
S5:   z[i] = y[i] + alpha*y[k-i-1];
    }
    for (int i = 0; i < k; i++) {
S6:   y[i] = z[i];
    }
S7: y[k] = alpha;
  }
#pragma endscop
}
