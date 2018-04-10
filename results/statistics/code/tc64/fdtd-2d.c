/**
 * This version is stamped on Apr. 14, 2015
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* fdtd-2d.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "fdtd-2d.h"


/* Array initialization. */
static
void init_array (int tmax,
		 int nx,
		 int ny,
		 DATA_TYPE POLYBENCH_2D(ex,NX,NY,nx,ny),
		 DATA_TYPE POLYBENCH_2D(ey,NX,NY,nx,ny),
		 DATA_TYPE POLYBENCH_2D(hz,NX,NY,nx,ny),
		 DATA_TYPE POLYBENCH_1D(_fict_,TMAX,tmax))
{
  int i, j;

  for (i = 0; i < tmax; i++)
    _fict_[i] = (DATA_TYPE) i;
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
      {
	ex[i][j] = ((DATA_TYPE) i*(j+1)) / nx;
	ey[i][j] = ((DATA_TYPE) i*(j+2)) / ny;
	hz[i][j] = ((DATA_TYPE) i*(j+3)) / nx;
      }
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int nx,
		 int ny,
		 DATA_TYPE POLYBENCH_2D(ex,NX,NY,nx,ny),
		 DATA_TYPE POLYBENCH_2D(ey,NX,NY,nx,ny),
		 DATA_TYPE POLYBENCH_2D(hz,NX,NY,nx,ny))
{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("ex");
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++) {
      if ((i * nx + j) % 20 == 0) fprintf(POLYBENCH_DUMP_TARGET, "\n");
      fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, ex[i][j]);
    }
  POLYBENCH_DUMP_END("ex");
  POLYBENCH_DUMP_FINISH;

  POLYBENCH_DUMP_BEGIN("ey");
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++) {
      if ((i * nx + j) % 20 == 0) fprintf(POLYBENCH_DUMP_TARGET, "\n");
      fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, ey[i][j]);
    }
  POLYBENCH_DUMP_END("ey");

  POLYBENCH_DUMP_BEGIN("hz");
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++) {
      if ((i * nx + j) % 20 == 0) fprintf(POLYBENCH_DUMP_TARGET, "\n");
      fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, hz[i][j]);
    }
  POLYBENCH_DUMP_END("hz");
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_fdtd_2d(int tmax,
		    int nx,
		    int ny,
		    DATA_TYPE POLYBENCH_2D(ex,NX,NY,nx,ny),
		    DATA_TYPE POLYBENCH_2D(ey,NX,NY,nx,ny),
		    DATA_TYPE POLYBENCH_2D(hz,NX,NY,nx,ny),
		    DATA_TYPE POLYBENCH_1D(_fict_,TMAX,tmax))
{
  int t, i, j;

/* TC Optimizing Compiler 0.2.26 */
/* ./tc ../examples/polybench/fdtd-2d.scop.c --correction-tiling --lex-scheduling --serial-codegen --debug -b 64 --report -R _PB_TMAX=100 -R _PB_NX=200 -R _PB_NY=240 */
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#pragma scop
if (_PB_NY >= 1) {
  for (int ii0 = 0; ii0 <= floord(_PB_TMAX - 1, 64); ii0 += 1) {
    if (_PB_NX >= 2) {
      for (int ii1 = 0; ii1 <= 1; ii1 += 1) {
        if (ii1 == 1) {
          for (int ii2 = 0; ii2 <= floord(_PB_NX - 1, 64); ii2 += 1) {
            for (int ii3 = 0; ii3 <= floord(_PB_NY - 1, 64); ii3 += 1) {
              for (int i2 = max(1, 64 * ii2); i2 <= min(_PB_NX - 1, 64 * ii2 + 63); i2 += 1) {
                for (int i3 = 64 * ii3; i3 <= min(_PB_NY - 1, 64 * ii3 + 63); i3 += 1) {
                  ey[i2][i3] = (ey[i2][i3] - (SCALAR_VAL(0.5) * (hz[i2][i3] - hz[i2 - 1][i3])));
                }
              }
              if (64 * ii3 + 64 >= _PB_NY) {
                for (int i0 = 64 * ii0 + 1; i0 <= min(_PB_TMAX - 1, 64 * ii0 + 63); i0 += 1) {
                  for (int i2 = max(1, 64 * ii2); i2 <= min(_PB_NX - 1, 64 * ii2 + 63); i2 += 1) {
                    ey[i2][_PB_NY - 1] = (ey[i2][_PB_NY - 1] - (SCALAR_VAL(0.5) * (hz[i2][_PB_NY - 1] - hz[i2 - 1][_PB_NY - 1])));
                  }
                }
              }
            }
          }
        } else {
          for (int ii2 = 0; ii2 <= floord(_PB_NY - 1, 64); ii2 += 1) {
            for (int i2 = 64 * ii2; i2 <= min(_PB_NY - 1, 64 * ii2 + 63); i2 += 1) {
              ey[0][i2] = _fict_[64 * ii0];
            }
            if (64 * ii2 + 64 >= _PB_NY) {
              for (int i0 = 64 * ii0 + 1; i0 <= min(_PB_TMAX - 1, 64 * ii0 + 63); i0 += 1) {
                ey[0][_PB_NY - 1] = _fict_[i0];
              }
            }
          }
        }
      }
    } else {
      for (int ii2 = 0; ii2 <= floord(_PB_NY - 1, 64); ii2 += 1) {
        for (int i0 = 64 * ii0; i0 <= min(_PB_TMAX - 1, 64 * ii0 + 63); i0 += 1) {
          for (int i2 = 64 * ii2; i2 <= min(_PB_NY - 1, 64 * ii2 + 63); i2 += 1) {
            ey[0][i2] = _fict_[i0];
          }
        }
      }
    }
    if (_PB_NY >= 2) {
      for (int ii1 = 2; ii1 <= min(3, _PB_NX + 1); ii1 += 1) {
        for (int ii2 = 0; ii2 <= (_PB_NX - ii1 + 1) / 64; ii2 += 1) {
          for (int ii3 = 0; ii3 <= (_PB_NY - ii1 + 1) / 64; ii3 += 1) {
            for (int i0 = 64 * ii0; i0 < min(_PB_TMAX, 64 * ii0 + 63 * ii1 - 125); i0 += 1) {
              if (ii1 == 3 && ii2 == 0 && i0 >= 64 * ii0 + 1) {
                if (_PB_NY >= 64 * ii3 + 66) {
                  for (int i2 = max(0, 64 * ii0 + 64 * ii3 - i0 + 1); i2 <= 64 * ii0 + 64 * ii3 - i0 + 64; i2 += 1) {
                    ey[0][i2] = _fict_[i0];
                  }
                } else {
                  for (int i2 = max(0, 64 * ii0 + 64 * ii3 - i0 + 1); i2 < _PB_NY - 1; i2 += 1) {
                    ey[0][i2] = _fict_[i0];
                  }
                }
              }
              if (ii1 == 3 && _PB_NX >= 64 * ii2 + 66 && i0 >= 64 * ii0 + 1) {
                for (int i2 = max(1, 64 * ii0 + 64 * ii2 - i0 + 1); i2 <= 64 * ii0 + 64 * ii2 - i0 + 64; i2 += 1) {
                  if (_PB_NY >= 64 * ii3 + 66 && 64 * ii2 >= i2 + 1) {
                    for (int i3 = max(0, 64 * ii0 + 64 * ii2 + 64 * ii3 - i0 - i2 + 1); i3 <= 64 * ii0 + 64 * ii2 + 64 * ii3 - i0 - i2 + 64; i3 += 1) {
                      ey[i2][i3] = (ey[i2][i3] - (SCALAR_VAL(0.5) * (hz[i2][i3] - hz[i2 - 1][i3])));
                    }
                  } else if (64 * ii3 + 65 >= _PB_NY) {
                    for (int i3 = max(max(0, 64 * ii0 + 64 * ii3 - i0 + 1), 64 * ii0 + 64 * ii2 + 64 * ii3 - i0 - i2 + 1); i3 < _PB_NY - 1; i3 += 1) {
                      ey[i2][i3] = (ey[i2][i3] - (SCALAR_VAL(0.5) * (hz[i2][i3] - hz[i2 - 1][i3])));
                    }
                  } else {
                    for (int i3 = max(0, 64 * ii0 + 64 * ii3 - i0 + 1); i3 <= 64 * ii0 + 64 * ii3 - i0 + 64; i3 += 1) {
                      ey[i2][i3] = (ey[i2][i3] - (SCALAR_VAL(0.5) * (hz[i2][i3] - hz[i2 - 1][i3])));
                    }
                  }
                }
              } else if (ii1 == 3 && 64 * ii2 + 65 >= _PB_NX && i0 >= 64 * ii0 + 1) {
                for (int i2 = max(1, 64 * ii0 + 64 * ii2 - i0 + 1); i2 < _PB_NX; i2 += 1) {
                  if (_PB_NY >= 64 * ii3 + 66 && 64 * ii2 >= i2 + 1) {
                    for (int i3 = max(0, 64 * ii0 + 64 * ii2 + 64 * ii3 - i0 - i2 + 1); i3 <= 64 * ii0 + 64 * ii2 + 64 * ii3 - i0 - i2 + 64; i3 += 1) {
                      ey[i2][i3] = (ey[i2][i3] - (SCALAR_VAL(0.5) * (hz[i2][i3] - hz[i2 - 1][i3])));
                    }
                  } else if (64 * ii3 + 65 >= _PB_NY) {
                    for (int i3 = max(max(0, 64 * ii0 + 64 * ii3 - i0 + 1), 64 * ii0 + 64 * ii2 + 64 * ii3 - i0 - i2 + 1); i3 < _PB_NY - 1; i3 += 1) {
                      ey[i2][i3] = (ey[i2][i3] - (SCALAR_VAL(0.5) * (hz[i2][i3] - hz[i2 - 1][i3])));
                    }
                  } else {
                    for (int i3 = max(0, 64 * ii0 + 64 * ii3 - i0 + 1); i3 <= 64 * ii0 + 64 * ii3 - i0 + 64; i3 += 1) {
                      ey[i2][i3] = (ey[i2][i3] - (SCALAR_VAL(0.5) * (hz[i2][i3] - hz[i2 - 1][i3])));
                    }
                  }
                }
              }
              if (ii1 == 3 && _PB_NY >= 64 * ii3 + 66 && i0 >= 64 * ii0 + 1) {
                for (int i2 = max(0, 64 * ii0 + 64 * ii2 - i0 + 1); i2 < 64 * ii2; i2 += 1) {
                  for (int i3 = max(1, 64 * ii0 + 64 * ii2 + 64 * ii3 - i0 - i2 + 1); i3 <= 64 * ii0 + 64 * ii2 + 64 * ii3 - i0 - i2 + 64; i3 += 1) {
                    ex[i2][i3] = (ex[i2][i3] - (SCALAR_VAL(0.5) * (hz[i2][i3] - hz[i2][i3 - 1])));
                  }
                }
                if (64 * ii2 + 65 >= _PB_NX) {
                  for (int i2 = 64 * ii2; i2 < _PB_NX - 1; i2 += 1) {
                    for (int i3 = max(1, 64 * ii0 + 64 * ii3 - i0 + 1); i3 <= 64 * ii0 + 64 * ii3 - i0 + 64; i3 += 1) {
                      ex[i2][i3] = (ex[i2][i3] - (SCALAR_VAL(0.5) * (hz[i2][i3] - hz[i2][i3 - 1])));
                    }
                  }
                } else {
                  for (int i2 = 64 * ii2; i2 <= 64 * ii0 + 64 * ii2 - i0 + 64; i2 += 1) {
                    for (int i3 = max(1, 64 * ii0 + 64 * ii3 - i0 + 1); i3 <= 64 * ii0 + 64 * ii3 - i0 + 64; i3 += 1) {
                      ex[i2][i3] = (ex[i2][i3] - (SCALAR_VAL(0.5) * (hz[i2][i3] - hz[i2][i3 - 1])));
                    }
                  }
                }
              } else if (ii1 == 3 && _PB_NX >= 64 * ii2 + 66 && 64 * ii3 + 65 >= _PB_NY && i0 >= 64 * ii0 + 1) {
                for (int i2 = max(0, 64 * ii0 + 64 * ii2 - i0 + 1); i2 <= 64 * ii0 + 64 * ii2 - i0 + 64; i2 += 1) {
                  for (int i3 = max(max(1, 64 * ii0 + 64 * ii3 - i0 + 1), 64 * ii0 + 64 * ii2 + 64 * ii3 - i0 - i2 + 1); i3 < _PB_NY; i3 += 1) {
                    ex[i2][i3] = (ex[i2][i3] - (SCALAR_VAL(0.5) * (hz[i2][i3] - hz[i2][i3 - 1])));
                  }
                }
              } else if (ii1 == 3 && 64 * ii2 + 65 >= _PB_NX && 64 * ii3 + 65 >= _PB_NY && i0 >= 64 * ii0 + 1) {
                for (int i2 = max(0, 64 * ii0 + 64 * ii2 - i0 + 1); i2 < _PB_NX - 1; i2 += 1) {
                  for (int i3 = max(max(1, 64 * ii0 + 64 * ii3 - i0 + 1), 64 * ii0 + 64 * ii2 + 64 * ii3 - i0 - i2 + 1); i3 < _PB_NY; i3 += 1) {
                    ex[i2][i3] = (ex[i2][i3] - (SCALAR_VAL(0.5) * (hz[i2][i3] - hz[i2][i3 - 1])));
                  }
                }
              }
              if (ii1 == 3) {
                for (int i2 = max(0, 64 * ii0 + 64 * ii2 - i0); i2 < 64 * ii2; i2 += 1) {
                  if (_PB_NY >= 64 * ii3 + 66) {
                    for (int i3 = max(0, 64 * ii0 + 64 * ii2 + 64 * ii3 - i0 - i2); i3 <= 64 * ii0 + 64 * ii2 + 64 * ii3 - i0 - i2 + 63; i3 += 1) {
                      hz[i2][i3] = (hz[i2][i3] - (SCALAR_VAL(0.7) * (((ex[i2][i3 + 1] - ex[i2][i3]) + ey[i2 + 1][i3]) - ey[i2][i3])));
                    }
                  } else {
                    for (int i3 = max(0, 64 * ii0 + 64 * ii2 + 64 * ii3 - i0 - i2); i3 < _PB_NY - 1; i3 += 1) {
                      hz[i2][i3] = (hz[i2][i3] - (SCALAR_VAL(0.7) * (((ex[i2][i3 + 1] - ex[i2][i3]) + ey[i2 + 1][i3]) - ey[i2][i3])));
                    }
                  }
                }
                if (64 * ii2 + 65 >= _PB_NX) {
                  for (int i2 = 64 * ii2; i2 < _PB_NX - 1; i2 += 1) {
                    for (int i3 = max(0, 64 * ii0 + 64 * ii3 - i0); i3 < 64 * ii3; i3 += 1) {
                      hz[i2][i3] = (hz[i2][i3] - (SCALAR_VAL(0.7) * (((ex[i2][i3 + 1] - ex[i2][i3]) + ey[i2 + 1][i3]) - ey[i2][i3])));
                    }
                    if (64 * ii3 + 65 >= _PB_NY) {
                      for (int i3 = 64 * ii3; i3 < _PB_NY - 1; i3 += 1) {
                        hz[i2][i3] = (hz[i2][i3] - (SCALAR_VAL(0.7) * (((ex[i2][i3 + 1] - ex[i2][i3]) + ey[i2 + 1][i3]) - ey[i2][i3])));
                      }
                    } else {
                      for (int i3 = 64 * ii3; i3 <= 64 * ii0 + 64 * ii3 - i0 + 63; i3 += 1) {
                        hz[i2][i3] = (hz[i2][i3] - (SCALAR_VAL(0.7) * (((ex[i2][i3 + 1] - ex[i2][i3]) + ey[i2 + 1][i3]) - ey[i2][i3])));
                      }
                    }
                  }
                } else {
                  for (int i2 = 64 * ii2; i2 <= 64 * ii0 + 64 * ii2 - i0 + 63; i2 += 1) {
                    for (int i3 = max(0, 64 * ii0 + 64 * ii3 - i0); i3 < 64 * ii3; i3 += 1) {
                      hz[i2][i3] = (hz[i2][i3] - (SCALAR_VAL(0.7) * (((ex[i2][i3 + 1] - ex[i2][i3]) + ey[i2 + 1][i3]) - ey[i2][i3])));
                    }
                    if (64 * ii3 + 65 >= _PB_NY) {
                      for (int i3 = 64 * ii3; i3 < _PB_NY - 1; i3 += 1) {
                        hz[i2][i3] = (hz[i2][i3] - (SCALAR_VAL(0.7) * (((ex[i2][i3 + 1] - ex[i2][i3]) + ey[i2 + 1][i3]) - ey[i2][i3])));
                      }
                    } else {
                      for (int i3 = 64 * ii3; i3 <= 64 * ii0 + 64 * ii3 - i0 + 63; i3 += 1) {
                        hz[i2][i3] = (hz[i2][i3] - (SCALAR_VAL(0.7) * (((ex[i2][i3 + 1] - ex[i2][i3]) + ey[i2 + 1][i3]) - ey[i2][i3])));
                      }
                    }
                  }
                }
              } else {
                for (int i2 = 64 * ii2; i2 <= min(_PB_NX - 1, 64 * ii2 + 63); i2 += 1) {
                  for (int i3 = max(1, 64 * ii3); i3 <= min(_PB_NY - 1, 64 * ii3 + 63); i3 += 1) {
                    ex[i2][i3] = (ex[i2][i3] - (SCALAR_VAL(0.5) * (hz[i2][i3] - hz[i2][i3 - 1])));
                  }
                }
              }
            }
            if (ii1 == 2 && 64 * ii2 + 64 >= _PB_NX) {
              for (int i0 = 64 * ii0 + 1; i0 <= min(_PB_TMAX - 1, 64 * ii0 + 63); i0 += 1) {
                for (int i3 = max(1, 64 * ii3); i3 <= min(_PB_NY - 1, 64 * ii3 + 63); i3 += 1) {
                  ex[_PB_NX - 1][i3] = (ex[_PB_NX - 1][i3] - (SCALAR_VAL(0.5) * (hz[_PB_NX - 1][i3] - hz[_PB_NX - 1][i3 - 1])));
                }
              }
            }
          }
        }
      }
    }
  }
}
#pragma endscop
}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int tmax = TMAX;
  int nx = NX;
  int ny = NY;

  /* Variable declaration/allocation. */
  POLYBENCH_2D_ARRAY_DECL(ex,DATA_TYPE,NX,NY,nx,ny);
  POLYBENCH_2D_ARRAY_DECL(ey,DATA_TYPE,NX,NY,nx,ny);
  POLYBENCH_2D_ARRAY_DECL(hz,DATA_TYPE,NX,NY,nx,ny);
  POLYBENCH_1D_ARRAY_DECL(_fict_,DATA_TYPE,TMAX,tmax);

  /* Initialize array(s). */
  init_array (tmax, nx, ny,
	      POLYBENCH_ARRAY(ex),
	      POLYBENCH_ARRAY(ey),
	      POLYBENCH_ARRAY(hz),
	      POLYBENCH_ARRAY(_fict_));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_fdtd_2d (tmax, nx, ny,
		  POLYBENCH_ARRAY(ex),
		  POLYBENCH_ARRAY(ey),
		  POLYBENCH_ARRAY(hz),
		  POLYBENCH_ARRAY(_fict_));


  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(nx, ny, POLYBENCH_ARRAY(ex),
				    POLYBENCH_ARRAY(ey),
				    POLYBENCH_ARRAY(hz)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(ex);
  POLYBENCH_FREE_ARRAY(ey);
  POLYBENCH_FREE_ARRAY(hz);
  POLYBENCH_FREE_ARRAY(_fict_);

  return 0;
}
