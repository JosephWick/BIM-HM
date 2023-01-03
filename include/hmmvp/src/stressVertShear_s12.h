/* stressVertShear_s12.h
 * Joseph Wick, adapted from a mablab code by Sylvain Barbot
 * Header file for s12 component stress calculation
 */

namespace s12 {

  // main function
  double s12::stressVertShear_s12(double x1, double x2, double x3,
                             double q1, double q2, double q3,
                             double L, double T, double W, double theta,
                             double epsv11p, double epsv12p, double epsv13p,
                             double epsv22p, double epsv23p, double epsv33p,
                             double G, double nu);

  // helpers

  // Green's functions
  double r1(double x1, double x2, double x3, double y1, double y2, double y3);
  double r2(double x1, double x2, double x3, double y1, double y2, double y3);

  double IU1d2(double y1,double y2,double y3, double lambda,double epsvkk,double G,
    double epsv11p,double epsv12p,double epsv13p,double epsv22p,double epsv23p,double epsv33p);
  double IU2d1(double y1,double y2, double y3, double lambda, double epsvkk, double G,
    double epsv11p,double epsv12p,double epsv13p,double epsv22p,double epsv23p,double epsv33p);

  // remove inelastic eigenstrain
  int heaviside(double x);
  double Omega(double x);
  double S(double x);

  //
  double J1112d2(double y1,double y2,double y3, double nu,double G,
        double x1,double x2,doube x3);
  double J1113d2(double y1,double y2,double y3, double nu,double G,
      double x1,double x2,double x3);
}
