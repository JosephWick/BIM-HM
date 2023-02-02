#include <stdio.h>
#include "stressVertShear_s12.h"

int main(){

  /*
  x1: 200000.000000, x2: 50423.931519, x3: 51745.059718
  y1: 222222.222222, y2: 26392.960684, y3: 37124.434120
  L: 222222.222222, T: 26392.960684, W: 37124.434120
  lambda: 30000.000000
  epsvkk: 0.000000
  */

  double x1 = 200000.0;
  double x2 = 50423.931519;
  double x3 = 51745.059718;

  double y1 = 222222.222222;
  double y2 = 26392.960684;
  double y3 = 37124.434120;

  double L = 222222.222222;
  double T = 26392.960684;
  double W = 37124.434120;

  double lambda = 30000.0;

  double nu = 0.25;

  // r1 and r2
  double r1_out = s12::r1(x1,x2,x3, y1,y2,y3);
  printf("r1: %f\n", r1_out);

  double r2_out = s12::r2(x1,x2,x3, y1,y2,y3);
  printf("r2: %f\n", r2_out);




}
