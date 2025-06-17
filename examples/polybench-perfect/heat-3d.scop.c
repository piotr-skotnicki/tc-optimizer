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
  for (int c0 = 1; c0 <= _PB_TSTEPS; c0 += 1) {
    for (int c1 = 2 * c0 + 1; c1 < _PB_N + 2 * c0; c1 += 1) {
      for (int c2 = max(2 * c0 + 1, -_PB_N + c1 + 3); c2 < min(_PB_N + 2 * c0, _PB_N + c1 - 2); c2 += 1) {
        for (int c3 = max(max(2 * c0 + 1, -_PB_N + c1 + 3), -_PB_N + c2 + 3); c3 < min(min(_PB_N + 2 * c0, _PB_N + c1 - 2), _PB_N + c2 - 2); c3 += 1) {
          if (_PB_N + 2 * c0 >= c1 + 2 && _PB_N + 2 * c0 >= c2 + 2 && _PB_N + 2 * c0 >= c3 + 2) {
S1:         B[-2 * c0 + c1][-2 * c0 + c2][-2 * c0 + c3] = ((((SCALAR_VAL(0.125) * ((A[-2 * c0 + c1 + 1][-2 * c0 + c2][-2 * c0 + c3] - (SCALAR_VAL(2.0) * A[-2 * c0 + c1][-2 * c0 + c2][-2 * c0 + c3])) + A[-2 * c0 + c1 - 1][-2 * c0 + c2][-2 * c0 + c3])) + (SCALAR_VAL(0.125) * ((A[-2 * c0 + c1][-2 * c0 + c2 + 1][-2 * c0 + c3] - (SCALAR_VAL(2.0) * A[-2 * c0 + c1][-2 * c0 + c2][-2 * c0 + c3])) + A[-2 * c0 + c1][-2 * c0 + c2 - 1][-2 * c0 + c3]))) + (SCALAR_VAL(0.125) * ((A[-2 * c0 + c1][-2 * c0 + c2][-2 * c0 + c3 + 1] - (SCALAR_VAL(2.0) * A[-2 * c0 + c1][-2 * c0 + c2][-2 * c0 + c3])) + A[-2 * c0 + c1][-2 * c0 + c2][-2 * c0 + c3 - 1]))) + A[-2 * c0 + c1][-2 * c0 + c2][-2 * c0 + c3]);
          }
          if (c1 >= 2 * c0 + 2 && c2 >= 2 * c0 + 2 && c3 >= 2 * c0 + 2) {
S2:         A[-2 * c0 + c1 - 1][-2 * c0 + c2 - 1][-2 * c0 + c3 - 1] = ((((SCALAR_VAL(0.125) * ((B[-2 * c0 + c1][-2 * c0 + c2 - 1][-2 * c0 + c3 - 1] - (SCALAR_VAL(2.0) * B[-2 * c0 + c1 - 1][-2 * c0 + c2 - 1][-2 * c0 + c3 - 1])) + B[-2 * c0 + c1 - 2][-2 * c0 + c2 - 1][-2 * c0 + c3 - 1])) + (SCALAR_VAL(0.125) * ((B[-2 * c0 + c1 - 1][-2 * c0 + c2][-2 * c0 + c3 - 1] - (SCALAR_VAL(2.0) * B[-2 * c0 + c1 - 1][-2 * c0 + c2 - 1][-2 * c0 + c3 - 1])) + B[-2 * c0 + c1 - 1][-2 * c0 + c2 - 2][-2 * c0 + c3 - 1]))) + (SCALAR_VAL(0.125) * ((B[-2 * c0 + c1 - 1][-2 * c0 + c2 - 1][-2 * c0 + c3] - (SCALAR_VAL(2.0) * B[-2 * c0 + c1 - 1][-2 * c0 + c2 - 1][-2 * c0 + c3 - 1])) + B[-2 * c0 + c1 - 1][-2 * c0 + c2 - 1][-2 * c0 + c3 - 2]))) + B[-2 * c0 + c1 - 1][-2 * c0 + c2 - 1][-2 * c0 + c3 - 1]);
          }
        }
      }
    }
  }
#pragma endscop
}
