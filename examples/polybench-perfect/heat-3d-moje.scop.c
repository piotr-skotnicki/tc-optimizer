double SCALAR_VAL(double);
int max(int, int);

int main()
{
#if 0
# define _PB_TSTEPS 20
# define _PB_N 10
#else
  short _PB_TSTEPS;
  short _PB_N;
#endif

  double A[_PB_N][_PB_N][_PB_N];
  double B[_PB_N][_PB_N][_PB_N];

#pragma scop
  for (int k0 = 0; k0 < _PB_TSTEPS; k0 += 1) {
    for (int k1 = 2 * k0; k1 < _PB_N + 2 * k0 - 1; k1 += 1) {
      for (int k2 = max(2 * k0, -_PB_N + k1 + 3); k2 < min(_PB_N + 2 * k0 - 1, _PB_N + k1 - 2); k2 += 1) {
        for (int k3 = max(max(2 * k0, -_PB_N + k1 + 3), -_PB_N + k2 + 3); k3 < min(min(_PB_N + 2 * k0 - 1, _PB_N + k1 - 2), _PB_N + k2 - 2); k3 += 1) {
          if (_PB_N + 2 * k0 >= k1 + 3 && _PB_N + 2 * k0 >= k2 + 3 && _PB_N + 2 * k0 >= k3 + 3) {
            B[-2 * k0 + k1 + 1][-2 * k0 + k2 + 1][-2 * k0 + k3 + 1] = ((((SCALAR_VAL(0.125) * ((A[-2 * k0 + k1 + 2][-2 * k0 + k2 + 1][-2 * k0 + k3 + 1] - (SCALAR_VAL(2.0) * A[-2 * k0 + k1 + 1][-2 * k0 + k2 + 1][-2 * k0 + k3 + 1])) + A[-2 * k0 + k1][-2 * k0 + k2 + 1][-2 * k0 + k3 + 1])) + (SCALAR_VAL(0.125) * ((A[-2 * k0 + k1 + 1][-2 * k0 + k2 + 2][-2 * k0 + k3 + 1] - (SCALAR_VAL(2.0) * A[-2 * k0 + k1 + 1][-2 * k0 + k2 + 1][-2 * k0 + k3 + 1])) + A[-2 * k0 + k1 + 1][-2 * k0 + k2][-2 * k0 + k3 + 1]))) + (SCALAR_VAL(0.125) * ((A[-2 * k0 + k1 + 1][-2 * k0 + k2 + 1][-2 * k0 + k3 + 2] - (SCALAR_VAL(2.0) * A[-2 * k0 + k1 + 1][-2 * k0 + k2 + 1][-2 * k0 + k3 + 1])) + A[-2 * k0 + k1 + 1][-2 * k0 + k2 + 1][-2 * k0 + k3]))) + A[-2 * k0 + k1 + 1][-2 * k0 + k2 + 1][-2 * k0 + k3 + 1]);
          }
          if (k1 >= 2 * k0 + 1 && k2 >= 2 * k0 + 1 && k3 >= 2 * k0 + 1) {
            A[-2 * k0 + k1][-2 * k0 + k2][-2 * k0 + k3] = ((((SCALAR_VAL(0.125) * ((B[-2 * k0 + k1 + 1][-2 * k0 + k2][-2 * k0 + k3] - (SCALAR_VAL(2.0) * B[-2 * k0 + k1][-2 * k0 + k2][-2 * k0 + k3])) + B[-2 * k0 + k1 - 1][-2 * k0 + k2][-2 * k0 + k3])) + (SCALAR_VAL(0.125) * ((B[-2 * k0 + k1][-2 * k0 + k2 + 1][-2 * k0 + k3] - (SCALAR_VAL(2.0) * B[-2 * k0 + k1][-2 * k0 + k2][-2 * k0 + k3])) + B[-2 * k0 + k1][-2 * k0 + k2 - 1][-2 * k0 + k3]))) + (SCALAR_VAL(0.125) * ((B[-2 * k0 + k1][-2 * k0 + k2][-2 * k0 + k3 + 1] - (SCALAR_VAL(2.0) * B[-2 * k0 + k1][-2 * k0 + k2][-2 * k0 + k3])) + B[-2 * k0 + k1][-2 * k0 + k2][-2 * k0 + k3 - 1]))) + B[-2 * k0 + k1][-2 * k0 + k2][-2 * k0 + k3]);
          }
        }
      }
    }
  }
#pragma endscop
}
