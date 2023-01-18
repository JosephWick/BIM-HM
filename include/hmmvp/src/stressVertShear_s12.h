/* stressVertShear_s12.h
 * Joseph Wick, adapted from a mablab code by Sylvain Barbot
 * Header file for s12 component stress calculation
 */

namespace s12 {

  // main function
  double stressVertShear_s12(double x1, double x2, double x3,
                             double q1, double q2, double q3,
                             double L, double T, double W, double theta,
                             double epsv11p, double epsv12p, double epsv13p,
                             double epsv22p, double epsv23p, double epsv33p,
                             double G, double nu);

  // helpers

  // Green's functions
  double r1(double x1,double x2,double x3, double y1,double y2,double y3);
  double r2(double x1,double x2,double x3, double y1,double y2,double y3);

  double IU1d2(double y1,double y2,double y3, double lambda,double epsvkk, double nu,double G,
    double epsv11p,double epsv12p,double epsv13p,double epsv22p,double epsv23p,double epsv33p);
  double IU2d1(double y1,double y2, double y3, double lambda,double epsvkk, double nu,double G,
    double epsv11p,double epsv12p,double epsv13p,double epsv22p,double epsv23p,double epsv33p);

  // remove inelastic eigenstrain
  int heaviside(double x);
  double Omega(double x);
  double S(double x);

  //
  double J1112d2(double y1,double y2,double y3, double nu,double G,
        double x1,double x2,double x3);
  double J1113d2(double y1,double y2,double y3, double nu,double G,
      double x1,double x2,double x3);
  double J2112d1(double y1,double y2,double y3, double nu,double G,
      double x1, double x2, double x3)
  double J2113d1(double y1,double y2,double y3, double nu,double G,
      double x1,double x2, double x3);
  double J2123d1(double y1,double y2,double y3, double nu,double G,
      double x1,double x2,double x3);
  double J1212d2(double y1,double y2,double y3, double nu,double G,
      double x1,double x2,double x3);
  double J1213d2(double y1,double y2,double y3, double nu,double G,
      double x1,double x2,double x3);
  double J1213d2(double y1,double y2,double y3, double nu,double G,
      double x1,double x2,double x3);
  double J1223d2(double y1,double y2,double y3, double nu,double G,
      double x1,double x2,double x3);
  double J2212d1(double y1,double y2,double y3, double nu,double G,
      double x1,double x2,double x3);
  double J2213d1(double y1,double y2,double y3, double nu,double G,
      double x1,double x2,double x3);
  double J2223d1(double y1,double y2,double y3, double nu,double G,
      double x1,double x2,double x3);
  double J1312d2(double y1,double y2,double y3, double nu,double G,
      double x1,double x2,double x3);
  double J1313d2(double y1,double y2,double y3, double nu,double G,
      double x1,double x2,double x3);
  double J1323d2(double y1,double y2,double y3, double nu,double G,
      double x1,double x2,double x3);
  double J2312d1(double y1,double y2,double y3, double nu,double G,
      double x1,double x2,double x3);
  double J2313d1(double y1,double y2,double y3, double nu,double G,
      double x1,double x2,double x3);
  double J2323d1(double y1,double y2,double y3, double nu,double G,
      double x1,double x2, double x3);
}
