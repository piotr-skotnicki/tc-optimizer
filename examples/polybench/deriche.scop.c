double SCALAR_VAL(double);

int main()
{
#if 0
# define _PB_W 1000
# define _PB_H 1000
#else
  int _PB_W;
  int _PB_H;
#endif

  int i,j;
  
  double alpha;
  double imgIn[_PB_W][_PB_H];
  double imgOut[_PB_W][_PB_H];
  double y1[_PB_W][_PB_H];
  double y2[_PB_W][_PB_H];
  double xm1, tm1, ym1, ym2;
  double xp1, xp2;
  double tp1, tp2;
  double yp1, yp2;
  double k;
  double a1, a2, a3, a4, a5, a6, a7, a8;
  double b1, b2, c1, c2;

#pragma scop
  for (i=0; i<_PB_W; i++) {
S1: ym1 = SCALAR_VAL(0.0);
S2: ym2 = SCALAR_VAL(0.0);
S3: xm1 = SCALAR_VAL(0.0);
    for (j=0; j<_PB_H; j++) {
S4:   y1[i][j] = a1*imgIn[i][j] + a2*xm1 + b1*ym1 + b2*ym2;
S5:   xm1 = imgIn[i][j];
S6:   ym2 = ym1;
S7:   ym1 = y1[i][j];
    }
  }

  for (i=0; i<_PB_W; i++) {
S8: yp1 = SCALAR_VAL(0.0);
S9: yp2 = SCALAR_VAL(0.0);
S10:xp1 = SCALAR_VAL(0.0);
S11:xp2 = SCALAR_VAL(0.0);
    for (j=_PB_H-1; j>=0; j--) {
S12:  y2[i][j] = a3*xp1 + a4*xp2 + b1*yp1 + b2*yp2;
S13:  xp2 = xp1;
S14:  xp1 = imgIn[i][j];
S15:  yp2 = yp1;
S16:  yp1 = y2[i][j];
    }
  }

  for (i=0; i<_PB_W; i++)
    for (j=0; j<_PB_H; j++) {
S17:  imgOut[i][j] = c1 * (y1[i][j] + y2[i][j]);
    }

  for (j=0; j<_PB_H; j++) {
S18:tm1 = SCALAR_VAL(0.0);
S19:ym1 = SCALAR_VAL(0.0);
S20:ym2 = SCALAR_VAL(0.0);
    for (i=0; i<_PB_W; i++) {
S21:  y1[i][j] = a5*imgOut[i][j] + a6*tm1 + b1*ym1 + b2*ym2;
S22:  tm1 = imgOut[i][j];
S23:  ym2 = ym1;
S24:  ym1 = y1 [i][j];
    }
  }
  
  for (j=0; j<_PB_H; j++) {
S25:tp1 = SCALAR_VAL(0.0);
S26:tp2 = SCALAR_VAL(0.0);
S27:yp1 = SCALAR_VAL(0.0);
S28:yp2 = SCALAR_VAL(0.0);
    for (i=_PB_W-1; i>=0; i--) {
S29:  y2[i][j] = a7*tp1 + a8*tp2 + b1*yp1 + b2*yp2;
S30:  tp2 = tp1;
S31:  tp1 = imgOut[i][j];
S32:  yp2 = yp1;
S33:  yp1 = y2[i][j];
    }
  }

  for (i=0; i<_PB_W; i++)
    for (j=0; j<_PB_H; j++)
S34:  imgOut[i][j] = c2*(y1[i][j] + y2[i][j]);
#pragma endscop
}

