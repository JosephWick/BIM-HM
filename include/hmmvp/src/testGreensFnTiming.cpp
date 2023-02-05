#include <stdio.h>
#include "stressVertShear_s12.h"

int main(){

  /*
  x1: 200000.000000, x2: 50423.931519, x3: 40000.000000
y1: 200000.000000, y2: 26392.960684, y3: 42754.167781
L: 200000.000000, T: 26392.960684, W: 42754.167781
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

  double G = 30000.0;

  double nu = 0.25;

  // test
  s12::stressVertShear_s12(x1,x1,x3, y1,y2,y3, L,T,W, 0, 0,1,0,0,0,0, G,nu);




}
