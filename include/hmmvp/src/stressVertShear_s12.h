/* stressVertShear_s12.h
 * Joseph Wick, adapted from a mablab code by Sylvain Barbot
 * Header file for s12 component stress calculation
 */

// main function
double s12::stressVertShear_s12(double x1, double x2, double x3,
                           double q1, double q2, double q3,
                           double L, double T, double W, double theta,
                           double epsv11p, double epsv12p, double epsv13p,
                           double epsv22p, double epsv23p, double epsv33p,
                           double G, double nu);
