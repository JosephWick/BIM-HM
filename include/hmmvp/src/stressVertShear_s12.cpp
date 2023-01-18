#define _USE_MATH_DEFINES
#include <cmath>

#include "stressVertShear_s12.h"

// main function
double stressVertShear_s12(double x1, double x2, double x3,
                           double y1, double y2, double y3,
                           double L, double T, double W, double theta,
                           double epsv11p, double epsv12p, double epsv13p,
                           double epsv22p, double epsv23p, double epsv33p,
                           double G, double nu)
{

  // x is receiver, y is sender

  // Lame paramter
  double lambda = G*2*nu/(1-2*nu);
  // isotropic strain
  double epsvkk = epsv11p + epsv22p + epsv33p;

  // rotate observation points to the shear-zone-centric system of coords
  double t1 =  (x1-q1)*cos(theta * M_PI/180) + (x2-q2)*sin(theta * M_PI/180);
  double x2 = -(x1-q1)*sin(theta * M_PI/180) + (x2-q2)*cos(theta * M_PI/180);
  double x1 = t1;


  // Displacement gradient
  double u12 = s12::IU1d2(L,T/2,q3+W, lambda,epsvkk,G,epsv11p,epsv12p,epsv13p,epsv22p,epsv23p,epsv33p)
              -s12::IU1d2(L,-T/2,q3+W, lambda,epsvkk,G,epsv11p,epsv12p,epsv13p,epsv22p,epsv23p,epsv33p)
              +s12::IU1d2(L,-T/2,q3, lambda,epsvkk,G,epsv11p,epsv12p,epsv13p,epsv22p,epsv23p,epsv33p)
              -s12::IU1d2(L,T/2,q3, lambda,epsvkk,G,epsv11p,epsv12p,epsv13p,epsv22p,epsv23p,epsv33p)
              -s12::IU1d2(0,T/2,q3+W, lambda,epsvkk,G,epsv11p,epsv12p,epsv13p,epsv22p,epsv23p,epsv33p)
              +s12::IU1d2(0,-T/2,q3+W, lambda,epsvkk,G,epsv11p,epsv12p,epsv13p,epsv22p,epsv23p,epsv33p)
              -s12::IU1d2(0,-T/2,q3, lambda,epsvkk,G,epsv11p,epsv12p,epsv13p,epsv22p,epsv23p,epsv33p)
              +s12::IU1d2(0,T/2,q3, lambda,epsvkk,G,epsv11p,epsv12p,epsv13p,epsv22p,epsv23p,epsv33p);
  double u21 = s12::IU2d1(L,T/2,q3+W, lambda,epsvkk,G,epsv11p,epsv12p,epsv13p,epsv22p,epsv23p,epsv33p)
              -s12::IU2d1(L,-T/2,q3+W, lambda,epsvkk,G,epsv11p,epsv12p,epsv13p,epsv22p,epsv23p,epsv33p)
              +s12::IU2d1(L,-T/2,q3, lambda,epsvkk,G,epsv11p,epsv12p,epsv13p,epsv22p,epsv23p,epsv33p)
              -s12::IU2d1(L,T/2,q3, lambda,epsvkk,G,epsv11p,epsv12p,epsv13p,epsv22p,epsv23p,epsv33p)
              -s12::IU2d1(0,T/2,q3+W, lambda,epsvkk,G,epsv11p,epsv12p,epsv13p,epsv22p,epsv23p,epsv33p)
              +s12::Iu2d1(0,-T/2,q3+W, lambda,epsvkk,G,epsv11p,epsv12p,epsv13p,epsv22p,epsv23p,epsv33p)
              -s12::IU2d1(0,-T/2,q3, lambda,epsvkk,G,epsv11p,epsv12p,epsv13p,epsv22p,epsv23p,epsv33p)
              +s12::IU2d1(0,T/2,q3, lambda,epsvkk,G,epsv11p,epsv12p,epsv13p,epsv22p,epsv23p,epsv33p);

  // strain
  double e12p = (u12+u21)/2;

  // rotate stress field to reference (unprimed) system of coordinates
  double e12 = e12p;

  double epsv12 = epsv12p;

  e12 - e12 - epsv12*S(x1/L)*Omega(x2/T)*S((x3-q3)/W);

  // stress components
  double s12 = 2*G*e12;

  return s12;
}

// Green's functions
double s12::r1(double x1, double x2, double x3, double y1, double y2, double y3){
  return sqrt( pow(x1-y1, 2) + pow(x2-y2,2) + pow(x3-y3,2) );
}
double s12::r2(double x1, double x2, double x3, double y1, double y2, double y3){
  return sqrt( pow(x1-y1, 2) + pow(x2-y2,2) + pow(x3+y3,2) );
}

double s12::IU2d1(double y1,double y2, double y3, double lambda,double epsvkk, double nu,double G,
  double epsv11p,double epsv12p,double epsv13p,double epsv22p,double epsv23p,double epsv33p){
  double s1 = (lambda*epsvkk + 2*G*epsv11p)*s12::J1123d2(y1,y2,y3, nu,G, x1,x2,x3);
  double s2 = 2*G*epsv12p*(s12::J1223d2(y1,y2,y3, nu,G, x1,x2,x3)
              + s12::J1113d2(y1,y2,y3, nu,G, x1,x2,x3));
  double s3 = 2*G*epsv13p*(s12::J1323d2(y1,y2,y3, nu,G, x1,x2,x3)
              + s12::J1112d2(y1,y2,y3, nu,G, x1,x2,x3));

  double s4 = (lambda*epsvkk + 2*G*epsv22p)*s12::J1213d2(y1,y2,y3, nu,G, x1,x2,x3);
  double s5 = 2*G*epsv23p*(s12::J1212d2(y1,y2,y3, nu,G, x1,x2,x3)
              + s12::J1313d2(y1,y2,y3, nu,G, x1,x2,x3));

  double s6 = (lambda*epsvkk + 2*G*epsv33p)*s12::J1312d2(y1,y2,y3, nu,G, x1,x2,x3);

  return s1+s2+s3+s4+s5+s6;
}
double s12::IU1d2(double y1,double y2,double y3, double lambda,double epsvkk, double nu,double G,
  double epsv11p,double epsv12p,double epsv13p,double epsv22p,double epsv23p,epsv33p){
  doubld s1 = (lambda*epsvkk + 2*G*epsv11p)*s12::J2123d1(y1,y2,y3, nu,G, x1,x2,x3);
  double s2 = 2*G*epsv12p*(s12::J2223d1(y1,y2,y3, nu,G, x1,x2,x3)
              + s12::J2113d1(y1,y2,y3, nu,G, x1,x2,x3));
  double s3 = 2*G*epsv13p*(s12::J2323d1(y1,y2,y3, nu,G, x1,x2,x3)
              + s12::J2112d1(y1,y2,y3, nu,G, x1,x2,x3));

  double s4 = (lambda*epsvkk + 2*G*epsv22p)*J2213d1(y1,y2,y3, nu,G, x1,x2,x3);
  double s5 = 2*G*epsv23p*(s12::J2212d1(y1,y2,y3, nu,G, x1,x2,x3)
              + s12::J2313d1(y1,y2,y3, nu,G, x1,x2,x3));

  double s6 = (lambda*epsvkk + 2*G*epsv33p)*(s12::J2312d1(y1,y2,y3, nu,G, x1,x2,x3));

  return s1+s2+s3+s4+s5+s6;
}

// remove inelastic eigenstrain
int s12::heaviside = [](x){
  return x>0;
};
double s12::Omega = [](x){
  return s12::heaviside(x+0.5) - s12::heaviside(x-0.5);
};
double s12::S = [](x){
  return s12::Omega(x-0.5);
};

// displacement gradient
double s12::J1112d2(double y1,double y2,double y3, double nu,double G,
    double x1,double x2,doube x3){
  double lr1 = s12::r1(x1,x2,x3, y1,y2,y3);
  double lr2 = s12::r2(x1,x2,x3, y1,y2,y3);

  double p1 = (1/16)*pow(1-nu,-1)*pow(M_PI,-1)*pow(G,-1);

  double s1 = -4*(nu-1)*pow(lr1,-1)*pow(lr1+x1-y1,-1)*pow(x2-y2,2);
  double s2 = -4*(nu-1)*pow(lr2,-1)*pow(lr2+x1-y1,-1)*pow(x2-y2,2);
  double s3 = -1*(4*nu-3)*(x1-y1)*pow(lr1+x2-y2,-1);
  double s4 = -1*(4*nu-3)*pow(lr1,-1)*(x1-y1)*(x2-y2)*pow(lr1+x2-y2,-1);
  double s5 = (5+4*nu*(2*nu-3))*(x1-y1)*pow(lr2+x2-y2,-1);
  double s6 = (5+4*nu*(2*nu-3))*pow(lr2,-1)*(x1-y1)*(x2-y2)*pow(lr2+x2-y2,-1);
  double s7 = 4*(nu-1)*lr1*(x1-y1)*pow(pow(x1-y1,2)+pow(x3-y3,2),-1)*
                pow(pow(x2-y2,2)+pow(x3-y3,2),-1)*pow(x3-y3,2);
  double s8 = -4*(nu-1)*pow(lr1,-1)*(x1-y1)*pow(x2-y2,2)*pow(pow(x1-y1,2)+
                pow(x3-y3,2),-1)*pow(pow(x2-y2,2)+pow(x3-y3,2),-1)*pow(x3-y3,2);
  double s9 = -4*(nu-1)*(2*nu-1)*(x1-y1)*pow(pow(x1-y1,2)+pow(x2-y2,2),-1)*(x3+y3);
  double s10= 2*pow(lr2,-1)*x3*(x1-y1)*y3*pow(pow(x1-y1,2)+pow(x3+y3,2),-1);
  double s11= 4*(nu-1)*(2*nu-1)*lr2*(x1-y1)*pow(pow(x1-y1,2)+pow(x2-y2,2),-1)*
                pow(x3+y3,2)*pow(pow(x1-y1,2)+pow(x3+y3,2),-1);
  double s12= -4*(nu-1)*(2*nu-1)*pow(lr2,-1)*(x1-y1)*pow(pow(x1-y1,2)+pow(x2-y2,2),-1)*
                pow(x2-y2,2)*pow(x3+y3,2)*pow(pow(x1-y1,2)+pow(x3+y3,2),-1);
  double s13= 4*(nu-1)*lr2*(x1-y1)*pow(x3+y3,2)*pow(pow(x1-y1,2)+pow(x3+y3,2),-1)*
                pow(pow(x2-y2,2)+pow(x3+y3,2),-1);
  double s14= -4*(nu-1)*pow(lr2,-1)*(x1-y1)*pow(x2-y2,2)*pow(x3+y3,2)*
                pow(pow(x1-y1,2)+pow(x3+y3,2),-1)*pow(pow(x2-y2,2)+pow(x3+y3,2),-1);
  double s15= -2*x3*(x1-y1)*pow(x2-y2,2)*y3*pow(pow(x1-y1,2)+pow(x3+y3,2),-1)*
                pow(pow(x1-y1,2)+pow(x2-y2,2)+pow(x3+y3,2),-1.5);
  double s16= (4-4*nu)*log(lr1+x1-y1);
  double s17= (4-4*nu)*log(lr2+x1-y1);

  return p1*(s1+s2+s3+s4+s5+s6+s7+s8+s9+s10+s11+s12+s13+s14+s15+s16+s17);
}

double s12::J1113d2(double y1,double y2,double y3, double nu,double G,
    double x1,double x2,double x3){
  double lr1 = s12::r1(x1,x2,x3, y1,y2,y3);
  double lr2 = s12::r2(x1,x2,x3, y1,y2,y3);

  double p1 = (1/16)*pow(1-nu,-1)*pow(M_PI,-1)*pow(G,-1);

  double s1 = x2*pow(pow(x2,2)+pow(x1-y1,2),-1)*(y1-x1);
  double s2 = 9*x2*pow(9*pow(x2,2)+pow(x1-y1,2),-1)*(y1-x1);
  double s3 = 4*pow(nu,2)*x2*pow(pow(nu*x2,2)+pow(x1-y1,2),-1)*(y1-x1);
  double s4 = -4*(nu-1)*pow(lr1,-1)*pow(lr1+x1-y1,-1)*(x2-y2)*(x3-y3);
  double s5 = -4(nu-1)*lr1*(x1-y1)*pow(pow(x1-y1,2)+pow(x2-y2,2),-1)*(x2-y2)*
                pow(pow(x2-y2,2)+pow(x3-y3,2),-1)*(x3-y3);
  double s6 = -4*(nu-1)*pow(lr1,-1)*(x1-y1)*pow(pow(x1-y1,2)+pow(x2-y2,2),-1)*
                *pow(x2-y2,3)*pow(pow(x2-y2,2)+pow(x3-y3,2),-1)*(x3-y3);
  double s7 = (3-4*nu)*pow(lr1,-1)*(x1-y1)*(x2-y2)*pow(lr1+x3-y3,-1);
  double s8 = 4(nu-1)*pow(lr2,-1)*pow(lr2+x1-y1,-1)*(x2-y2)*(x3+y3);
  double s9 = -(3-6*nu+4*nu*nu)*pow(lr2,-1)*(x1-y1)*(x2-y2)*pow(lr2+x3+y3,-1);
  double s10= -4*(nu-1)*(2*nu-1)*(x1-y1)*pow(pow(x1-y1,2)+pow(x2-y2,2),-2)*(x2-y2)*
                y3*(2*x3+y3);
  double s11= 2*(nu-1)*(2*nu-1)*pow(lr2,-2)*(x1-y1)*pow(pow(x1-y1,2)+pow(x2-y2,2),-1)*
                (x2-y2)*y3*(2*x3+y3);
  double s12= -4*pow(lr2,-1)*(x1-y1)*pow(pow(x1-y1,2)+pow(x2-y2,2),-1)*(x2-y2)*
                (y3-3*nu*(x3+y3)+2*nu*nu*(x3+y3));
  double s13= 4*(nu-1)*lr2*(x1-y1)*pow(pow(x1-y1,2)+pow(x2-y2,2),-1)*(x2-y2)*
                (x3+y3)*pow(pow(x2-y2,2)+pow(x3+y3,2),-1);
  double s14= 4*(nu-1)*pow(lr2,-1)*(x1-y1)*pow(pow(x1-y1,2)+pow(x2-y2,2),-1)*
                pow(x2-y2,3)*(x3+y3)*pow(pow(x2-y2,2)+pow(x3+y3,2),-1);
  double s15= 2*(x1-y1)*pow(pow(x1-y1,2)+pow(x2-y2,2),-1)*(x2-y2)*
                pow(pow(x1-y1,2)+pow(x2-y2,2)+pow(x3+y3,2),-1.5)*(
                -3*nu*(x3*(x3*x3+pow(x1-y1,2)+pow(x2-y2,2))+(x3*(-2*lr2+3*x3)+
                pow(x1-y1,2)+pow(x2-y2,2))*y3-(lr2-3*x3)*y3*y3 +y3*y3 +y3*y3*y3)+
                2*nu*nu*(x3*(x3*x3+pow(x1-y1,2)+(x2-y2,2))+(x3*(-2*lr2+3*x3)+
                pow(x1-y1,2)+pow(x2-y2,2))-(lr2-3*x3)*y3*y3 +y3*y3*y3)+
                y3*(pow(x1-y1,2)+pow(x2-y2,2)-(lr2-x3-y3)*(2*x2+y3)) );
  double s16= 4*pow(lr2,-1)*(x1-y1)*pow(pow(x1-y1,2)+pow(x2-y2,2),-2)*(x2-y2)*(
                -3*nu*(x3+y3)*(pow(x1-y1,2)+pow(x2-y2,2)+pow(x3+y3,2)) +2*nu*nu*
                (x3+y3)*(pow(x1-y1,2)+pow(x2-y2,2)+pow(x3+y3,2)) +y3*(pow(x1-y1,2)+
                pow(x2-y2,2)+(x3+y3)(2*x3+y3))
                );
  double s17 = atan2(x1-y1,-x2);
  double s18 = -3*atan2(3*x2,x1-y1);
  double s19 = 4*nu*atan2(-nu*x2,x1-y1);
  double s20 = (4-4*nu)*atan2(lr1*(x2-y2), (x1-y1)*(x3-y3));
  double s21 = 4*(nu-1)*atan2(lr2*(x2-y2), (x1-y1)*(x3+y3));

  return p1*(s1+s2+s3+s4+s5+s6+s7+s8+s9+s10+s11+s12+s13+s14+s15+s16+s17+s18+s19+s20+s21);
}

double s12::J1123d2(double y1,double y2,double y3, double nu,double G,
    double x1,double x2,double x3) {
  double lr1 = s12::r1(x1,x2,x3, y1,y2,y3);
  double lr2 = s12::r2(x1,x2,x3, y1,y2,y3);

  double p1 = (1/16)*pow(1-nu,-1)*pow(M_PI,-1)*pow(G,-1);

  double s1 = pow(x1,2)*pow(pow(x1,2)*pow(x2-y2,2),-1);
  double s2 = 9*pow(x1,2)*pow(9*pow(x1,2)+pow(x-y2,2),-1);
  double s3 = 4*pow(nu,2)*pow(x1,2)*pow(pow(nu*x,2)+pow(x2-y2,2),-1);
  double s4 = -2*(nu-1)*(2*nu-1)*pow(x1-y1,2)*pow(pow(x1-y1,2)+pow(x2-y2,2),-1);
  double s5 = pow(y1,2)*pow(pow(y1,2)+pow(x2-y2,2),-1);
  double s6 = 9*pow(y1,2)*pow(9*pow(y1,2)+pow(x2-y2,2),-1);
  double s7 = 4*pow(nu*y1,2)*pow(pow(nu*y1,2)+pow(x2-y2,2),-1);
  double s8 = 2*x3*pow(lr2-x2+y2,-1);
  double s9 = 2*pow(lr2,-1)*x3*(y2-x2)*pow(lr2-x2+y2,-1);
  double s10= -(4*nu-3)*pow(lr1+x2-y2,-1)*(x3-y3);
  double s11= -(4*nu-3)*pow(lr1,-1)*(x2-y2)*pow(lr1+x2-y2,-1)*(x3-y3);
  double s12= 2*(2*nu-1)*lr1*pow(x1-y1,2)*pow(pow(x1-y1,2)+pow(x2-y2,2),-1)*
                pow(pow(x1-y1,2)+pow(x3-y3,2),-1)*(x3-y3);
  double s13= -2*(2*nu-1)*pow(lr1,-1)*pow(x1-y1,2)*pow(pow(x1-y1,2)+pow(x2-y2,2),-1)*
                pow(x2-y2,2)*pow(pow(x1-y1,2)+pow(x3-y3,2),-1)*(x3-y3);
  double s14= -(4*nu-3)*pow(lr1,-1)*pow(x2-y2,2)*pow(lr1+x3-y3,-1);
  double s15= -(5+4*nu*(2*nu-3))*pow(lr2+x2-y2,-1)*(x3+y3);
  double s16= -(5+4*nu*(2*nu-3))*pow(lr2,-1)*(x2-y2)*pow(lr2+x2-y2,-1)*(x3+y3);
  double s17= -(3-6*nu+4*pow(nu,2))*pow(lr2,-1)*pow(x2-y2,2)*pow(lr2+x3+y3,-1);
  double s18= 4*(nu-1)*(2*nu-1)*pow(pow(x1-y1,2)+pow(x2-y2,2),-2)*pow(x2-y2,2)*
                y3*(2*x3+y3);
  double s19= -2*(nu-1)*(2*nu-1)*pow(lr2,-2)*pow(pow(x1-y1,2)+pow(x2-y2,2),-1)*
                pow(x2-y2,2)*y3*(2*x3+y3);
  double s20= 2*(2*nu-1)*pow(pow(x1-y1,2)+pow(x2-y2,2),-1)*pow(pow(x1-y1,2)+pow(x3+y3,2),-1)*
                ((2*nu-1)*lr2*pow(x1-y1,2)*(x3+y3)-(nu-1)*y3*(2*x3+y3)*
                (pow(x1-y1,2)+pow(x3+y3,2)));
  double s21= -4*pow(lr2,-1)*pow(pow(x1-y1,2)+pow(x2-y2,2),-2)*pow(x2-y2,2)*
                pow(pow(x1-y1,2)+pow(x3-y3,2),-1)*(pow(x1,4)*y3-4*pow(x1,3))*y1*y3-
                3*nu*(x3+y3)*(pow(x1-y1,2)+pow(x3+y3,2))*(pow(x1-y1,2)+pow(x2-y2,2)+pow(x3+y3,2))+
                2*pow(nu,2)*(x3+y3)*(pow(x1-y1,2)+pow(x3+y3,2))*(pow(x1-y1,2)+pow(x2-y2,2)+
                pow(x3+y3,2))-2*x1*y1*y3*(pow(x2-y2,2)+2*(pow(y1,2)+(x3+y3)*(2*x3+y3)))+
                y3*(2*pow(x3,4)+7*pow(x3,3)*y3+(pow(y1,2)+pow(y3,2))*(pow(y1,2)+
                pow(y2,2)+pow(y3,2)+x3+y3*(6*pow(y1,2)+3*pow(y2,2)+5*pow(y3,2))+
                pow(x3,2)*(4*pow(y1,2)+2*pow(y2,2)+9*pow(y3,2))+pow(x2,2)*(pow(y1,2)+
                (x3+y3)*(2*x3+y3))-2*x2*y2*(pow(y1,2)+(x3+y3)*(2*x3+y3)))+
                pow(x1,2)*y3*(pow(x2-y2,2)+2*(3*pow(y1,2)+(x3+y3)*(2*x3+y3))));
  double s22= -2*pow(pow(x1-y1,2)+pow(x2-y2,2),-1)*pow(x2-y2,2)*
                pow(pow(x1-y1,2)_+pow(x3+y3,2),-1)*pow(pow(x1-y1,2)+pow(x2-y2,2)
                pow(x3+y3,2),-1.5)*(pow(x1,4)*(y3-3*nu*(x3+y3)+2*pow(nu,2)*(x3+y3))
                -4*pow(x1,3)*y1*(y3-3*nu*(x3+y3)+2*pow(nu,2)*(x3+y3))
                -3*nu*(pow(y1,2)+pow(x3+y3,2))*(pow(x3,3)+3*pow(x3,2)*y3+y3*(
                pow(y1,2)+pow(y2,2)-lr2*y3+pow(y3,2))+x3*(pow(y1,2)+pow(y2,2)-
                2*lr2*y3+3*pow(y3,2)))+2*pow(nu,2)*(pow(y1,2)+pow(x3+y3,2))*
                (pow(x3,3)+3*pow(x3,2)*y3+y3*(pow(y1,2)+pow(y2,2)-lr2*y3+
                pow(y3,2))+x3*(pow(y1,2)+pow(y2,2)-2*lr2*y3+3*pow(y3,2)))+
                y3*(2*pow(x3,4)+7*pow(x3,3)*y3+(pow(y1,2)+pow(y3,2))*(pow(y1,2)+
                pow(y2,2)+pow(y3,2))+x3*y3*(6*pow(y1,2)+3*pow(y2,2)+5*pow(y3,2))
                +pow(x3,2)*(4*pow(y1,2)+2*pow(y2,2)+9*pow(y3,2))-lr2*(2*x3+y3)*
                (pow(y1,2)+pow(x3+y3,2)))+2*x2*y2*(3*nu*(x3+y3)*(pow(y1,2)+
                pow(x3+y3,2))-2*pow(nu,2)*(x3+y3)*(pow(y1,2)+pow(x3+y3,2))-
                y3*(pow(y1,2)+(x3+y3)*(2*x3+y3)))+2*x1*y1*(((-1)+nu)
                *((-1)+2*nu)*lr2*y3*(2*x3+y3)+pow(x2,2)*(-y3+3*nu*(x3+y3)
                -2*pow(nu,2)*(x3+y3))+2*x2*y2*(y3-3*nu*(x3+y3)+2*pow(nu,2)*
                (x3+y3))+3*nu*(x3+y3)*(2*pow(y1,2)+pow(y2,2)+2*pow(x3+y3,2))-2*
                pow(nu,2)*(x3+y3)*(2*pow(y1,2)+pow(y2,2)+2*pow(x3+y3,2))-y3*(2*
                pow(y1,2)+pow(y2,2)+2*(x3+y3)*(2*x3+y3)))+pow(x1,2)(-((-1)+nu)*((
                -1)+2*nu)*lr2*y3*(2*x3+y3)+2*x2*y2*(-y3+3*nu*(x3+y3)
                -2*pow(nu,2)*(x3+y3))+pow(x2,2)*(y3-3*nu*(x3+y3)+2*pow(nu,2)*(x3+
                y3))-3*nu*(x3+y3)*(6*pow(y1,2)+pow(y2,2)+2*pow(x3+y3,2)+2*pow(nu,2)*(
                x3+y3)*(6*pow(y1,2)+pow(y2,2)+2*pow(x3+y3,2))+y3*(6*pow(y1,2)+pow(y2,2)
                +2*(x3+y3)*(2*x3+y3))));
  double s23= 2*pow(lr2,-1)*pow(pow(x1-y1,2)+pow(x2-y2,2),-1*
  pow(pow(x1-y1,2)+pow(x3+y3,2),-1)*(-x3*pow(y1,2)*pow(y2,2)+2*pow(x3,4)*
  y3+4*pow(x3,2)*pow(y1,2)*y3+pow(y1,4)*y3+6*pow(x3,2)*pow(y2,2)*y3+2*
  pow(y1,2)*pow(y2,2)*y3+7*pow(x3,3)*pow(y3,2)+6*x3*pow(y1,2)*pow(y3,2)+
  9*x3*pow(y2,2)*pow(y3,2)+9*pow(x3,2)*pow(y3,3)+2*pow(y1,2)*pow(y3,3)+
  3*pow(y2,2)*pow(y3,3)+5*x3*pow(y3,4)+pow(y3,5)-nu*(x3+y3)*(3*pow(x3,4)+
  6*pow(x3,2)*pow(y1,2)+3*pow(y1,4)+9*pow(x3,2)*pow(y2,2)+5*pow(y1,2)*
  pow(y2,2)+6*x3*(2*(pow(x3,2)+pow(y1,2))+3*pow(y2,2))*y3+3*(6*pow(x3,2)+
  2*pow(y1,2)+3*pow(y2,2))*pow(y3,2)+12*x3*pow(y3,3)+3*pow(y3,4))+pow(x1,4)*
  (y3-3*nu*(x3+y3)+2*pow(nu,2)*(x3+y3))-4*pow(x1,3)*y1*(y3-3*nu*(x3+y3)+2
  *pow(nu,2)*(x3+y3))+2*pow(nu,2)*(x3+y3)*(pow(x3,4)+pow(y1,4)+4*pow(x3,3)
  *y3+3*pow(y2,2)*pow(y3,2)+pow(y3,4)+pow(y,2)*(pow(y2,2)+2*pow(y3,2))+pow(x3,2)
  *(2*pow(y1,2)+3*pow(y2,2)+6*pow(y3,2))+x3*(4*pow(y1,2)*y3+6*pow(y2,2)*y3+4*
  pow(y3,3)))+pow(x2,2)*(-x3*pow(y1,2)+6*pow(x3,2)*y3+2*pow(y1,2)*y3+9*x3*
  pow(y3,2)+3*pow(y3,3)+2*pow(nu,2)*(x3+y3)*(pow(y1,2)+3*pow(x3+y3,2))-nu*
  (x3+y3)*(5*pow(y1,2)+9*pow(x3+y3,2)))+2*x2*y2*(x3*pow(y1,2)-6*pow(x3,2)*y3
  -2*pow(y1,2)*y3-9*x3*pow(y3,2)-3*pow(y3,3)-2*pow(nu,2)*(x3+y3)*(pow(y1,2)+3
  *pow(x3+y3,2))+nu*(x3+y3)*(5*pow(y1,2)+9*pow(x3+y3,2)))+2*x1*y1*(pow(x2,2)*x3-2*
  x2*x3*y2+x3*pow(y2,2)-2*pow(x2,2)*y3-4*pow(x3,2)*y3-2*pow(y1,2)*
  y3+4*x2*y2*y3-2*pow(y2,2)*y3-6*x3*pow(y3,2)-2*pow(y3,3)
  -2*pow(nu,2)*(x3+y3)*(pow(x2-y2,2)+2*(pow(y1,2)+pow(x3+y3,2)))+nu*(
  x3+y3)*(5*pow(x2-y2,2)+6*(pow(y1,2)+pow(x3+y3,2))))+pow(x1,2)*((-1)
  *pow(x2,2)*x3+2*x2*x3*y2-x3*pow(y,2)+2*pow(x2,2)*y3+4*pow(x3,2)*
  y3+6*pow(y1,2)*y3-4*x2*y2*y3+2*pow(y2,2)*y3+6*x3*pow(y3,2)+2*
  pow(y3,3)+2*pow(nu,2)*(x3+y3)*(pow(x2-y2,2)+2*(3*pow(y1,2)+pow(x3+y3,2)
))-nu*(x3+y3)*(5*pow(x2-y2,2)+6*(3*pow(y1,2)+pow(x3+y3,2)))));
double s24 = (3-4*nu)*log(lr1+x3-y3);
double s25 = (-3+6*nu-4*pow(nu,2))*log(lr2+x3+y3);

return p1*(s1+s2+s3+s4+s5+s6+s7+s8+s9+s10+s11+s12+s13+s14+s15+s16+s17+s18+s19
        +s20+s21+s22+s23+s24+s25);
}

double s12::J2112d1(double y1,double y2,double y3, double nu,double G,
    double x1, double x2, double x3){
  double lr1 = s12::r1(x1,x2,x3, y1,y2,y3);
  double lr2 = s12::r2(x1,x2,x3, y1,y2,y3);

  return (1/16)*pow(1-nu,-1)*pow(M_PI,-1)*((1+8*((-1)+nu)*nu)*pow(lr2,-1)
    *(x1-y1)+pow(lr1,-1)*(-x1+y1)-4*((-1)+nu)*((-1)+
      2*nu)*pow(lr2,-1)*(x1-y1)*(x3+y3)*pow(lr2+x3+y3,-1)+2*x3*
      (x1-y1)*y3*pow(pow(x1-y1,2)+pow(x2-y2,2)+pow(x3+y3,2),-3/2))*pow(G,-1);

}

double s12::J2113d1(double y1,double y2,double y3, double nu,double G,
    double x1,double x2, double x3){
  double lr1=s12::r1(x1,x2,x3, y1,y2,y3);
  double lr2=s12::r2(x1,x2,x3, y1,y2,y3);

  return (1/16)*pow(1-nu,-1)*pow(M_PI,-1)*(-pow(lr1,-1)*(x1-
      y1)*(x2-y2)*pow(lr1+x3-y3,-1)-((-1)-2*nu+
      4*pow(nu,2))*pow(lr2,-1)*(x1-y1)*(x2-y2)*pow(lr2+x3+y3,-1)
      -4*((-1)+nu)*((-1)+2*nu)*(x1-y1)*pow(pow(x1-y1,2)
      +pow(x2-y2,2),-2)*(x2-y2)*y3*(2*x3+y3)+2*((
      -1)+nu)*((-1)+2*nu)*pow(lr2,-2)*(x1-y1)*pow(pow(x1-y1,2)+
      pow(x2-y2,2),-1)*(x2-y2)*y3*(2*x3+y3)-4*pow(lr2,-1)
      *(x1-y1)*pow(pow(x1-y1,2)+pow(x2-y2,2),-1)*(
      x2-y2)*(y3-3*nu*(x3+y3)+2*pow(nu,2)*(x3+y3))+2*(x1+(-1)
      *y1)*pow(pow(x1-y1,2)+pow(x2-y2,2),-1)*(x2-y2)*pow(
      pow(x1-y1,2)+pow(x2-y2,2)+pow(x3+y3,2),-3/2)*((-3)*nu*
      (x3*(pow(x3,2)+pow(x1-y1,2)+pow(x2-y2,2))+(x3*((-2)*lr2+3*
      x3)+pow(x1-y1,2)+pow(x2-y2,2))*y3-(lr2-3*x3)*
      pow(y3,2)+pow(y3,3))+2*pow(nu,2)*(x3*(pow(x3,2)+pow(x1-y1,2)+pow(x2-y2,2))
      +(x3*((-2)*lr2+3*x3)+pow(x1-y1,2)+pow(x2-y2,2))*y3+
      (-1)*(lr2-3*x3)*pow(y3,2)+pow(y3,3))+y3*(pow(x1-y1,2)+pow(x2-
      y2,2)-(lr2-x3-y3)*(2*x3+y3)))+4*pow(lr2,-1)*(
      x1-y1)*pow(pow(x1-y1,2)+pow(x2-y2,2),-2)*(x2+(-1)
      *y2)*((-3)*nu*(x3+y3)*(pow(x1-y1,2)+pow(x2-y2,2)+pow(x3+
      y3,2))+2*pow(nu,2)*(x3+y3)*(pow(x1-y1,2)+pow(x2-y2,2)+pow(x3+
      y3,2))+y3*(pow(x1-y1,2)+pow(x2-y2,2)+(x3+y3)*(2*x3+y3)
    )))*pow(G,-1);
}

double s12::J2123d1(double y1,double y2,double y3, double nu,double G,
    double x1,double x2,double x3){
  double lr1=s12::r1(x1,x2,x3, y1,y2,y3);
  double lr2=s12::r2(x1,x2,x3, y1,y2,y3);
  return (1/16)*pow(1-nu,-1)*pow(M_PI,-1)*pow(G,-1)*(pow(-lr1,-1)*pow(
      x1-y1,2)*pow(lr1+x3-y3,-1)-((-1)-2*nu+4*
      pow(nu,2))*pow(lr2,-1)*pow(x1-y1,2)*pow(lr2+x3+y3,-1)-4*((-1)+
      nu)*((-1)+2*nu)*pow(x1-y1,2)*pow(pow(x1-y1,2)+pow(x2-
      y2,2),-2)*y3*(2*x3+y3)+2*((-1)+nu)*((-1)+2*nu)*pow(pow(x1+(
      -1)*y1,2)+pow(x2-y2,2),-1)*y3*(2*x3+y3)+2*((-1)+nu)*
      ((-1)+2*nu)*pow(lr2,-2)*pow(x1-y1,2)*pow(pow(x1-y1,2)+pow(x2+(
      -1)*y2,2),-1)*y3*(2*x3+y3)+2*pow(x1-y1,2)*pow(pow(x1+(-1)
      *y1,2)+pow(x2-y2,2),-1)*pow(pow(x1-y1,2)+pow(x2-y2,2)
      +pow(x3+y3,2),-3/2)*((-3)*nu*(x3*(pow(x3,2)+pow(x1-y1,2)+pow(
      x2-y2,2))+(x3*((-2)*lr2+3*x3)+pow(x1-y1,2)+pow(x2-
      y2,2))*y3-(lr2-3*x3)*pow(y3,2)+pow(y3,3))+2*pow(nu,2)*(x3*(
      pow(x3,2)+pow(x1-y1,2)+pow(x2-y2,2))+(x3*((-2)*lr2+3*x3)+pow(
      x1-y1,2)+pow(x2-y2,2))*y3-(lr2-3*x3)*pow(y3,2)+
      pow(y3,3))+y3*(pow(x1-y1,2)+pow(x2-y2,2)-(lr2-x3+(
      -1)*y3)*(2*x3+y3)))+4*pow(lr2,-1)*pow(x1-y1,2)*pow(pow(x1-
      y1,2)+pow(x2-y2,2),-2)*((-3)*nu*(x3+y3)*(pow(x1-y1,
      2)+pow(x2-y2,2)+pow(x3+y3,2))+2*pow(nu,2)*(x3+y3)*(pow(x1-y1,
      2)+pow(x2-y2,2)+pow(x3+y3,2))+y3*(pow(x1-y1,2)+pow(x2-
      y2,2)+(x3+y3)*(2*x3+y3)))+pow(lr2,-1)*pow(pow(x1-y1,2)+pow(x2+(-1)
      *y2,2),-1)*(6*nu*(x3+y3)*(3*pow(x1-y1,2)+pow(x2-
      y2,2)+pow(x3+y3,2))-4*pow(nu,2)*(x3+y3)*(3*pow(x1-y1,2)+pow(x2+
      (-1)*y2,2)+pow(x3+y3,2))-2*y3*(3*pow(x1-y1,2)+pow(x2-
      y2,2)+(x3+y3)*(2*x3+y3)))+log(lr2+x3+y3)-log(lr1+x3-
      y3)+2*(1-2*nu)*nu*log(lr2+x3+y3));
}

double s12::J1212d2(double y1,double y2,double y3, double nu, double G,
    double x1,double x2,double x3){
  double lr1=s12::r1(x1,x2,x3, y1,y2,y3);
  double lr2=s12::r2(x1,x2,x3, y1,y2,y3);

  return (1/16)*pow(1-nu,-1)*pow(M_PI,-1)*((1+8*((-1)+nu)*nu)*pow(lr2,
      -1)*(x2-y2)+pow(lr1,-1)*(-x2+y2)-4*((-1)+nu)*((-1)+
      2*nu)*pow(lr2,-1)*(x2-y2)*(x3+y3)*pow(lr2+x3+y3,-1)+2*x3*
      (x2-y2)*y3*pow(pow(x1-y1,2)+pow(x2-y2,2)+pow(x3+y3,2),
      -3/2))*pow(G,-1);
}

double s12::J1213d2(double y1,double y2,double y3, double nu, double G,
    double x1,double x2,double x3){
  double lr1=s12::r1(x1,x2,x3, y1,y2,y3);
  double lr2=s12::r2(x1,x2,x3, y1,y2,y3);

  return (1/16)*pow(1-nu,-1)*pow(M_PI,-1)*pow(G,-1)*(-pow(lr1,-1)*pow(
      x2-y2,2)*pow(lr1+x3-y3,-1)-((-1)-2*nu+4*
      pow(nu,2))*pow(lr2,-1)*pow(x2-y2,2)*pow(lr2+x3+y3,-1)+2*((-1)+nu)
      *((-1)+2*nu)*pow(pow(x1-y1,2)+pow(x2-y2,2),-1)*y3*(
      2*x3+y3)-4*((-1)+nu)*((-1)+2*nu)*pow(pow(x1-y1,2)+pow(x2+(
      -1)*y2,2),-2)*pow(x2-y2,2)*y3*(2*x3+y3)+2*((-1)+nu)
      *((-1)+2*nu)*pow(lr2,-2)*pow(pow(x1-y1,2)+pow(x2-y2,2),-1)
      *pow(x2-y2,2)*y3*(2*x3+y3)+2*pow(pow(x1-y1,2)+pow(x2+(
      -1)*y2,2),-1)*pow(x2-y2,2)*pow(pow(x1-y1,2)+pow(x2-
      y2,2)+pow(x3+y3,2),-3/2)*((-3)*nu*(x3*(pow(x3,2)+pow(x1-y1,
      2)+pow(x2-y2,2))+(x3*((-2)*lr2+3*x3)+pow(x1-y1,2)+pow(x2+(
      -1)*y2,2))*y3-(lr2-3*x3)*pow(y3,2)+pow(y3,3))+2*pow(nu,2)*(x3*
      (pow(x3,2)+pow(x1-y1,2)+pow(x2-y2,2))+(x3*((-2)*lr2+3*x3)+pow(
      x1-y1,2)+pow(x2-y2,2))*y3-(lr2-3*x3)*pow(y3,2)+
      pow(y3,3))+y3*(pow(x1-y1,2)+pow(x2-y2,2)-(lr2-x3+(
      -1)*y3)*(2*x3+y3)))+4*pow(lr2,-1)*pow(pow(x1-y1,2)+(x2-
      y2,2),-2)*(x2-y2,2)*((-3)*nu*(x3+y3)*(pow(x1-
      y1,2)+pow(x2-y2,2)+pow(x3+y3,2))+2*pow(nu,2)*(x3+y3)*(pow(x1-
      y1,2)+pow(x2-y2,2)+pow(x3+y3,2))+y3*(pow(x1-y1,2)+pow(x2+(-1)
      *y2,2)+(x3+y3)*(2*x3+y3)))+pow(lr2,-1)*pow(pow(x1-y1,2)+pow(x2+(
      -1)*y2,2),-1)*(6*nu*(x3+y3)*(pow(x1-y1,2)+3*pow(x2+(-1)
      *y2,2)+pow(x3+y3,2))-4*pow(nu,2)*(x3+y3)*(pow(x1-y1,2)+3*pow(
      x2-y2,2)+pow(x3+y3,2))-2*y3*(pow(x1-y1,2)+3*pow(x2+(
      -1)*y2,2)+(x3+y3)*(2*x3+y3)))+log(lr2+x3+y3)-log(lr1+x3+(
      -1)*y3)+2*(1-2*nu)*nu*log(lr2+x3+y3));
}

double s12::J1223d2(double y1,double y2,double y3, double nu,double G,
    double x1,double x2,double x3){
  double lr1=s12::r1(x1,x2,x3, y1,y2,y3);
  double lr2=s12::r2(x1,x2,x3, y1,y2,y3);

  return (1/16)*pow(1-nu,-1)*pow(M_PI,-1)*(-pow(lr1,-1)*(x1-
        y1)*(x2-y2)*pow(lr1+x3-y3,-1)-((-1)-2*nu+
        4*pow(nu,2))*pow(lr2,-1)*(x1-y1)*(x2-y2)*pow(lr2+x3+y3,
        -1)-4*((-1)+nu)*((-1)+2*nu)*(x1-y1)*pow(pow(x1-y1,
        2)+pow(x2-y2,2),-2)*(x2-y2)*y3*(2*x3+y3)+2*((
        -1)+nu)*((-1)+2*nu)*pow(lr2,-2)*(x1-y1)*pow(pow(x1-y1,2)+
        pow(x2-y2,2),-1)*(x2-y2)*y3*(2*x3+y3)-4*pow(lr2,
        -1)*(x1-y1)*pow(pow(x1-y1,2)+pow(x2-y2,2),-1)*(
        x2-y2)*(y3-3*nu*(x3+y3)+2*pow(nu,2)*(x3+y3))+2*(x1+(-1)
        *y1)*pow(pow(x1-y1,2)+pow(x2-y2,2),-1)*(x2-y2)*pow(
        pow(x1-y1,2)+pow(x2-y2,2)+pow(x3+y3,2),-3/2)*((-3)*nu*
        (x3*(pow(x3,2)+pow(x1-y1,2)+pow(x2-y2,2))+(x3*((-2)*lr2+3*
        x3)+pow(x1-y1,2)+pow(x2-y2,2))*y3-(lr2-3*x3)*
        pow(y3,2)+pow(y3,3))+2*pow(nu,2)*(x3*(pow(x3,2)+pow(x1-y1,2)+pow(x2-y2,
        2))+(x3*((-2)*lr2+3*x3)+pow(x1-y1,2)+pow(x2-y2,2))*y3+
        (-1)*(lr2-3*x3)*pow(y3,2)+pow(y3,3))+y3*(pow(x1-y1,2)+pow(x2-
        y2,2)-(lr2-x3-y3)*(2*x3+y3)))+4*pow(lr2,-1)*(
        x1-y1)*pow(pow(x1-y1,2)+pow(x2-y2,2),-2)*(x2+(-1)
        *y2)*((-3)*nu*(x3+y3)*(pow(x1-y1,2)+pow(x2-y2,2)+pow(x3+
        y3,2))+2*pow(nu,2)*(x3+y3)*(pow(x1-y1,2)+pow(x2-y2,2)+pow(x3+
        y3,2))+y3*(pow(x1-y1,2)+pow(x2-y2,2)+(x3+y3)*(2*x3+y3)
      )))*pow(G,-1);
}

double s12::J2212d1(double y1,double y2,double y3, double nu,double G,
    double x1,double x2,double x3){
  double lr1=s12::r1(x1,x2,x3, y1,y2,y3);
  double lr2=s12::r2(x1,x2,x3, y1,y2,y3);

  return (1/16)*pow(1-nu,-1)*pow(M_PI,-1)*pow(G,-1)*(pow(x3,2)*pow(pow(x3,2)+pow(x1+
      (-1)*y1,2),-1)+9*pow(x3,2)*pow(9*pow(x3,2)+pow(x1-y1,2),-1)+
      4*pow(nu,2)*pow(x3,2)*pow(pow(nu,2)*pow(x3,2)+pow(x1-y1,2),-1)-((-3)
      +4*nu)*pow(lr1+x1-y1,-1)*(x2-y2)-((-3)+4*nu)
      *pow(lr1,-1)*(x1-y1)*pow(lr1+x1-y1,-1)*(x2-y2)+(
      5+4*nu*((-3)+2*nu))*pow(lr2+x1-y1,-1)*(x2-y2)+(5+
      4*nu*((-3)+2*nu))*pow(lr2,-1)*(x1-y1)*pow(lr2+x1-y1,
      -1)*(x2-y2)-4*((-1)+nu)*pow(lr1,-1)*pow(x1-y1,2)*pow(
      lr1+x2-y2,-1)-4*((-1)+nu)*pow(lr2,-1)*pow(x1-y1,2)
      *pow(lr2+x2-y2,-1)+4*((-1)+nu)*lr1*(x2-y2)*pow(pow(x1+
      (-1)*y1,2)+pow(x3-y3,2),-1)*pow(pow(x2-y2,2)+pow(x3-
      y3,2),-1)*pow(x3-y3,2)-4*((-1)+nu)*pow(lr1,-1)*pow(x1+(
      -1)*y1,2)*(x2-y2)*pow(pow(x1-y1,2)+pow(x3-y3,2),-1)
      *pow(pow(x2-y2,2)+pow(x3-y3,2),-1)*pow(x3-y3,2)+(
      -4)*((-1)+nu)*((-1)+2*nu)*pow(pow(x1-y1,2)+pow(x2-y2,2),
      -1)*(x2-y2)*(x3+y3)+pow(y3,2)*pow(pow(x1-y1,2)+pow(y3,2),
      -1)+9*pow(y3,2)*pow(pow(x1-y1,2)+9*pow(y3,2),-1)+4*pow(nu,2)*pow(y3,2)*
      pow(pow(x1-y1,2)+pow(nu,2)*pow(y3,2),-1)+2*pow(lr2,-1)*x3*(x2-
      y2)*y3*pow(pow(x2-y2,2)+pow(x3+y3,2),-1)+4*((-1)+nu)*((-1)+
      2*nu)*lr2*pow(pow(x1-y1,2)+pow(x2-y2,2),-1)*(x2-
      y2)*pow(x3+y3,2)*pow(pow(x2-y2,2)+pow(x3+y3,2),-1)-4*((-1)+
      nu)*((-1)+2*nu)*pow(lr2,-1)*pow(x1-y1,2)*pow(pow(x1-y1,2)+pow(
      x2-y2,2),-1)*(x2-y2)*pow(x3+y3,2)*pow(pow(x2-y2,2)
      +pow(x3+y3,2),-1)+4*((-1)+nu)*lr2*(x2-y2)*pow(x3+y3,2)
      *pow(pow(x1-y1,2)+pow(x3+y3,2),-1)*pow(pow(x2-y2,2)+pow(x3+
      y3,2),-1)-4*((-1)+nu)*pow(lr2,-1)*pow(x1-y1,2)*(x2+(
      -1)*y2)*pow(x3+y3,2)*pow(pow(x1-y1,2)+pow(x3+y3,2),-1)*pow(pow(x2+(
      -1)*y2,2)+pow(x3+y3,2),-1)-2*x3*pow(x1-y1,2)*(x2+(-1)
      *y2)*y3*pow(pow(x2-y2,2)+pow(x3+y3,2),-1)*pow(pow(x1-y1,2)+
      pow(x2-y2,2)+pow(x3+y3,2),-3/2)+4-4*nu*log(lr1+x2-y2)
      +4-4*nu*log(lr2+x2-y2));
}

double s12::J2213d1(double y1,double y2,double y3, double nu,double G,
    double x1,double x2,double x3){

  double lr1=s12::r1(x1,x2,x3, y1,y2,y3);
  double lr2=s12::r2(x1,x2,x3, y1,y2,y3);

  return (1/16)*pow(1-nu,-1)*pow(M_PI,-1)*pow(G,-1)*(pow(x2,2)*pow(pow(x2,2)+pow(x1+
      (-1)*y1,2),-1)+9*pow(x2,2)*pow(9*pow(x2,2)+pow(x1-y1,2),-1)+
      4*pow(nu,2)*pow(x2,2)*pow(pow(nu,2)*pow(x2,2)+pow(x1-y1,2),-1)+2*x3*pow(lr2+
      (-1)*x1+y1,-1)+2*pow(lr2,-1)*x3*(-x1+y1)*pow(lr2-x1+
      y1,-1)-2*((-1)+nu)*((-1)+2*nu)*pow(pow(x1-y1,2)+pow(x2+(
      -1)*y2,2),-1)*pow(x2-y2,2)+pow(y2,2)*pow(pow(x1-y1,2)+
      pow(y2,2),-1)+9*pow(y2,2)*pow(pow(x1-y1,2)+9*pow(y2,2),-1)+4*
      pow(nu,2)*pow(y2,2)*pow(pow(x1-y1,2)+pow(nu,2)*pow(y2,2),-1)-((-3)+
      4*nu)*pow(lr1+x1-y1,-1)*(x3-y3)-((-3)+4*nu)*
      pow(lr1,-1)*(x1-y1)*pow(lr1+x1-y1,-1)*(x3-y3)+2*
      ((-1)+2*nu)*lr1*pow(pow(x1-y1,2)+pow(x2-y2,2),-1)*pow(x2+(
      -1)*y2,2)*pow(pow(x2-y2,2)+pow(x3-y3,2),-1)*(x3-
      y3)-2*((-1)+2*nu)*pow(lr1,-1)*pow(x1-y1,2)*pow(pow(x1-y1,
      2)+pow(x2-y2,2),-1)*pow(x2-y2,2)*pow(pow(x2-y2,2)+pow(
      x3-y3,2),-1)*(x3-y3)-((-3)+4*nu)*pow(lr1,-1)
      *pow(x1-y1,2)*pow(lr1+x3-y3,-1)-(5+4*nu*((-3)+
      2*nu))*pow(lr2+x1-y1,-1)*(x3+y3)-(5+4*nu*((-3)+2*
      nu))*pow(lr2,-1)*(x1-y1)*pow(lr2+x1-y1,-1)*(x3+y3)+(
      -1)*(3-6*nu+4*pow(nu,2))*pow(lr2,-1)*pow(x1-y1,2)*pow(lr2+x3+y3
      ,-1)+4*((-1)+nu)*((-1)+2*nu)*pow(x1-y1,2)*pow(pow(x1-
      y1,2)+pow(x2-y2,2),-2)*y3*(2*x3+y3)-2*((-1)+nu)*((
      -1)+2*nu)*pow(lr2,-2)*pow(x1-y1,2)*pow(pow(x1-y1,2)+pow(x2+(-1)
      *y2,2),-1)*y3*(2*x3+y3)+2*((-1)+2*nu)*pow(pow(x1-y1,
      2)+pow(x2-y2,2),-1)*pow(pow(x2-y2,2)+pow(x3+y3,2),-1)*
      (((-1)+2*nu)*lr2*pow(x2-y2,2)*(x3+y3)-((-1)+nu)*y3*
      (2*x3+y3)*(pow(x2-y2,2)+pow(x3+y3,2)))-4*pow(lr2,-1)*pow(x1+(
      -1)*y1,2)*pow(pow(x1-y1,2)+pow(x2-y2,2),-2)*pow(pow(x2+(-1)
      *y2,2)+pow(x3+y3,2),-1)*(pow(x2,4)*y3+4*pow(x2,2)*pow(x3,2)*y3+2*
      pow(x3,4)*y3+pow(x2,2)*pow(y1,2)*y3+2*pow(x3,2)*pow(y1,2)*y3-4*pow(x2,3)*y2*
      y3-8*x2*pow(x3,2)*y2*y3-2*x2*pow(y1,2)*y2*y3+6*pow(x2,2)*
      pow(y2,2)*y3+4*pow(x3,2)*pow(y2,2)*y3+pow(y1,2)*pow(y2,2)*y3-4*x2*pow(y2,3)*
      y3+pow(y2,4)*y3+6*pow(x2,2)*x3*pow(y3,2)+7*pow(x3,3)*pow(y3,2)+3*x3*pow(y1,2)*
      pow(y3,2)-12*x2*x3*y2*pow(y3,2)+6*x3*pow(y2,2)*pow(y3,2)+2*pow(x2,2)*
      pow(y3,3)+9*pow(x3,2)*pow(y3,3)+pow(y1,2)*pow(y3,3)-4*x2*y2*pow(y3,3)+2*pow(y2,2)*
      pow(y3,3)+5*x3*pow(y3,4)+pow(y3,5)-3*nu*(x3+y3)*(pow(x2-y2,2)+pow(x3+
      y3,2))*(pow(x1-y1,2)+pow(x2-y2,2)+pow(x3+y3,2))+2*pow(nu,2)*
      (x3+y3)*(pow(x2-y2,2)+pow(x3+y3,2))*(pow(x1-y1,2)+pow(x2+(-1)
      *y2,2)+pow(x3+y3,2))+pow(x1,2)*y3*(pow(x2-y2,2)+(x3+y3)*(2*x3+
      y3))-2*x1*y1*y3*(pow(x2-y2,2)+(x3+y3)*(2*x3+y3)))+(
      -2)*pow(x1-y1,2)*pow(pow(x1-y1,2)+pow(x2-y2,2),-1)*
      pow(pow(x2-y2,2)+pow(x3+y3,2),-1)*pow(pow(x1-y1,2)+pow(x2-
      y2,2)+pow(x3+y3,2),-3/2)*(pow(x2,4)*y3-2*lr2*pow(x2,2)*x3*y3+4*
      pow(x2,2)*pow(x3,2)*y3-2*lr2*pow(x3,3)*y3+2*pow(x3,4)*y3+pow(x2,2)*pow(y1,2)*
      y3+2*pow(x3,2)*pow(y1,2)*y3-4*pow(x2,3)*y2*y3+4*lr2*x2*x3*y2*y3+(
      -8)*x2*pow(x3,2)*y2*y3-2*x2*pow(y1,2)*y2*y3+6*pow(x2,2)*pow(y2,2)*
      y3-2*lr2*x3*pow(y2,2)*y3+4*pow(x3,2)*pow(y2,2)*y3+pow(y1,2)*pow(y2,2)*y3+(
      -4)*x2*pow(y2,3)*y3+pow(y2,4)*y3-lr2*pow(x2,2)*pow(y3,2)+6*pow(x2,2)*x3*
      pow(y3,2)-5*lr2*pow(x3,2)*pow(y3,2)+7*pow(x3,3)*pow(y3,2)+3*x3*pow(y1,2)*pow(y3,2)+
      2*lr2*x2*y2*pow(y3,2)-12*x2*x3*y2*pow(y3,2)-lr2*pow(y2,2)*
      pow(y3,2)+6*x3*pow(y2,2)*pow(y3,2)+2*pow(x2,2)*pow(y3,3)-4*lr2*x3*pow(y3,3)+9*
      pow(x3,2)*pow(y3,3)+pow(y1,2)*pow(y3,3)-4*x2*y2*pow(y3,3)+2*pow(y2,2)*pow(y3,3)+(-1)
      *lr2*pow(y3,4)+5*x3*pow(y3,4)+pow(y3,5)-3*nu*(x3*(pow(x3,2)+pow(y1,2)+pow(x2+(
      -1)*y2,2))+((-2)*lr2*x3+3*pow(x3,2)+pow(y1,2)+pow(x2-y2,2))*y3+(
      -1)*(lr2-3*x3)*pow(y3,2)+pow(y3,3))*(pow(x2-y2,2)+pow(x3+y3,2))+
      2*pow(nu,2)*(x3*(pow(x3,2)+pow(y1,2)+pow(x2-y2,2))+((-2)*lr2*x3+3*
      pow(x3,2)+pow(y1,2)+pow(x2-y2,2))*y3-(lr2-3*x3)*pow(y3,2)+
      pow(y3,3))*(pow(x2-y2,2)+pow(x3+y3,2))+2*x1*y1*(3*nu*(x3+y3)*
      (pow(x2-y2,2)+pow(x3+y3,2))-2*pow(nu,2)*(x3+y3)*(pow(x2-y2,
      2)+pow(x3+y3,2))-y3*(pow(x2-y2,2)+(x3+y3)*(2*x3+y3)))+
      pow(x1,2)*((-3)*nu*(x3+y3)*(pow(x2-y2,2)+pow(x3+y3,2))+2*
      pow(nu,2)*(x3+y3)*(pow(x2-y2,2)+pow(x3+y3,2))+y3*(pow(x2-y2,
      2)+(x3+y3)*(2*x3+y3))))+2*pow(lr2,-1)*pow(pow(x1-y1,2)+pow(x2+(-1)
      *y2,2),-1)*pow(pow(x2-y2,2)+pow(x3+y3,2),-1)*(-x3*
      pow(y1,2)*pow(y2,2)+2*pow(x3,4)*y3+6*pow(x3,2)*pow(y1,2)*y3+4*pow(x3,2)*pow(y2,2)*y3+
      2*pow(y1,2)*pow(y2,2)*y3+pow(y2,4)*y3+7*pow(x3,3)*pow(y3,2)+9*x3*pow(y1,2)*pow(y3,2)+
      6*x3*pow(y2,2)*pow(y3,2)+9*pow(x3,2)*pow(y3,3)+3*pow(y1,2)*pow(y3,3)+2*pow(y2,2)*
      pow(y3,3)+5*x3*pow(y3,4)+pow(y3,5)+pow(x2,4)*(y3-3*nu*(x3+y3)+2*pow(nu,2)*(
      x3+y3))-4*pow(x2,3)*y2*(y3-3*nu*(x3+y3)+2*pow(nu,2)*(x3+y3))+
      2*x2*y2*(x3*pow(y1,2)-4*pow(x3,2)*y3-2*pow(y1,2)*y3-2*
      pow(y2,2)*y3-6*x3*pow(y3,2)-2*pow(y3,3)-2*pow(nu,2)*(x3+y3)*(
      pow(y1,2)+2*pow(y2,2)+2*pow(x3+y3,2))+nu*(x3+y3)*(5*pow(y1,2)+6*pow(y2,2)+6*pow(
      x3+y3,2)))+pow(x2,2)*(-x3*pow(y1,2)+4*pow(x3,2)*y3+2*pow(y1,2)*y3+6*
      pow(y2,2)*y3+6*x3*pow(y3,2)+2*pow(y3,3)+2*pow(nu,2)*(x3+y3)*(pow(y1,2)+6*
      pow(y2,2)+2*pow(x3+y3,2))-nu*(x3+y3)*(5*pow(y1,2)+18*pow(y2,2)+6*pow(
      x3+y3,2)))+pow(x1,2)*(-pow(x2,2)*x3+2*x2*x3*y2-x3*pow(y2,2)+
      2*pow(x2,2)*y3+6*pow(x3,2)*y3-4*x2*y2*y3+2*pow(y2,2)*y3+9*x3*
      pow(y3,2)+3*pow(y3,3)+2*pow(nu,2)*(x3+y3)*(pow(x2-y2,2)+3*pow(x3+y3,2))
      -nu*(x3+y3)*(5*pow(x2-y2,2)+9*pow(x3+y3,2)))+2*x1*
      y1*(pow(x2,2)*x3-2*x2*x3*y2+x3*pow(y2,2)-2*pow(x2,2)*y3-6*
      pow(x3,2)*y3+4*x2*y2*y3-2*pow(y2,2)*y3-9*x3*pow(y3,2)-3*
      pow(y3,3)-2*pow(nu,2)*(x3+y3)*(pow(x2-y2,2)+3*pow(x3+y3,2))+nu*(
      x3+y3)*(5*pow(x2-y2,2)+9*pow(x3+y3,2)))-nu*(x3+y3)*(
      3*pow(x3,4)+12*pow(x3,3)*y3+3*pow(pow(y2,2)+pow(y3,2),2)+3*pow(x3,2)*(3*(y1,2)+2*
      pow(y2,2)+6*pow(y3,2))+pow(y1,2)*(5*pow(y2,2)+9*pow(y3,2))+6*x3*y3*(3*pow(y1,2)+
      2*(pow(y2,2)+pow(y3,2))))+2*pow(nu,2)*(x3+y3)*(pow(x3,4)+4*pow(x3,3)*y3+pow(pow(y2,2)+
      pow(y3,2),2)+pow(y1,2)*(pow(y2,2)+3*pow(y3,2))+pow(x3,2)*(3*pow(y1,2)+2*pow(y2,2)+6*
      pow(y3,2))+x3*(6*pow(y1,2)*y3+4*y3*(pow(y2,2)+pow(y3,2)))))+(3-4*nu)*log(
      lr1+x3-y3)+((-3)+6*nu-4*pow(nu,2))*log(lr2+x3+y3));
}

double s12::J2223d1(double y1,double y2,double y3, double nu,double G,
    double x1,double x2,double x3) {
  double lr1=s12::r1(x1,x2,x3, y1,y2,y3);
  double lr2=s12::r2(x1,x2,x3, y1,y2,y3);

  return (1/16)*pow(1-nu,-1)*pow(M_PI,-1)*pow(G,-1)*(x1*pow(pow(x1,2)+pow(x2+(
      -1)*y2,2),-1)*(-x2+y2)+9*x1*pow(9*pow(x1,2)+pow(x2-y2,
      2),-1)*(-x2+y2)+4*pow(nu,2)*x1*pow(pow(nu,2)*pow(x1,2)+pow(x2-
      y2,2),-1)*(-x2+y2)-4*((-1)+nu)*pow(lr1,-1)*(x1-
      y1)*pow(lr1+x2-y2,-1)*(x3-y3)-4*((-1)+nu)*lr1*(
      x1-y1)*pow(pow(x1-y1,2)+pow(x2-y2,2),-1)*(x2+(-1)
      *y2)*pow(pow(x1-y1,2)+pow(x3-y3,2),-1)*(x3-y3)+(
      -4)*((-1)+nu)*pow(lr1,-1)*pow(x1-y1,3)*pow(pow(x1-y1,2)+pow(x2+
      (-1)*y2,2),-1)*(x2-y2)*pow(pow(x1-y1,2)+pow(x3-
      y3,2),-1)*(x3-y3)-((-3)+4*nu)*pow(lr1,-1)*(x1+(
      -1)*y1)*(x2-y2)*pow(lr1+x3-y3,-1)+4*((-1)+nu)*
      pow(lr2,-1)*(x1-y1)*pow(lr2+x2-y2,-1)*(x3+y3)-(3+
      (-6)*nu+4*pow(nu,2))*pow(lr2,-1)*(x1-y1)*(x2-y2)*pow(lr2+
      x3+y3,-1)-4*((-1)+nu)*((-1)+2*nu)*(x1-y1)*pow(pow(x1+(
      -1)*y1,2)+pow(x2-y2,2),-2)*(x2-y2)*y3*(2*x3+y3)
      +2*((-1)+nu)*((-1)+2*nu)*pow(lr2,-2)*(x1-y1)*pow(pow(x1-
      y1,2)+pow(x2-y2,2),-1)*(x2-y2)*y3*(2*x3+y3)
      -4*pow(lr2,-1)*(x1-y1)*pow(pow(x1-y1,2)+pow(x2-y2,2),
      -1)*(x2-y2)*(y3-3*nu*(x3+y3)+2*pow(nu,2)*(x3+y3))+4*((
      -1)+nu)*lr2*(x1-y1)*pow(pow(x1-y1,2)+pow(x2-y2,2),
      -1)*(x2-y2)*(x3+y3)*pow(pow(x1-y1,2)+pow(x3+y3,2),-1)+
      4*((-1)+nu)*pow(lr2,-1)*pow(x1-y1,3)*pow(pow(x1-y1,2)+pow(x2+(
      -1)*y2,2),-1)*(x2-y2)*(x3+y3)*pow(pow(x1-y1,2)+pow(x3+
      y3,2),-1)+2*(x1-y1)*pow(pow(x1-y1,2)+pow(x2-y2,
      2),-1)*(x2-y2)*pow(pow(x1-y1,2)+pow(x2-y2,2)+pow(x3+
      y3,2),-3/2)*((-3)*nu*(x3*(pow(x3,2)+pow(x1-y1,2)+pow(x2+(-1)
      *y2,2))+(x3*((-2)*lr2+3*x3)+pow(x1-y1,2)+pow(x2-y2,2))
      *y3-(lr2-3*x3)*pow(y3,2)+pow(y3,3))+2*pow(nu,2)*(x3*(pow(x3,2)+pow(x1+(
      -1)*y1,2)+pow(x2-y2,2))+(x3*((-2)*lr2+3*x3)+pow(x1-y1,
      2)+pow(x2-y2,2))*y3-(lr2-3*x3)*pow(y3,2)+pow(y3,3))+y3*(pow(
      x1-y1,2)+pow(x2-y2,2)-(lr2-x3-y3)*(
      2*x3+y3)))+4*pow(lr2,-1)*(x1-y1)*pow(pow(x1-y1,2)+pow(x2+(-1)
      *y2,2),-2)*(x2-y2)*((-3)*nu*(x3+y3)*(pow(x1-y1,
      2)+pow(x2-y2,2)+pow(x3+y3,2))+2*pow(nu,2)*(x3+y3)*(pow(x1-y1,
      2)+pow(x2-y2,2)+pow(x3+y3,2))+y3*(pow(x1-y1,2)+pow(x2-
      y2,2)+(x3+y3)*(2*x3+y3)))+atan2(x2-y2,(-1)*x1)-3*
      atan2(3*x1,x2-y2)+4*nu*atan2(-nu*x1,x2-y2)+(
      4-4*nu)*atan2(lr1*(x1-y1),(x2-y2)*(x3-y3))
      +4*((-1)+nu)*atan2(lr2*(x1-y1),(x2-y2)*(x3+y3)));
}

double s12::J1312d2(double y1,double y2,double y3, double nu,double G,
    double x1,double x2,double x3){
  double lr1=s12::r1(x1,x2,x3, y1,y2,y3);
  double lr2=s12::r2(x1,x2,x3, y1,y2,y3);
  return (-1/16)*pow(1-nu,-1)*pow(M_PI,-1)*pow(G,-1)*(4*((-1)+nu)*((
      -1)+2*nu)*pow(x1-y1,2)*pow(pow(x1-y1,2)+pow(x2-y2,2),
      -1)+pow(lr1+x2-y2,-1)*(x3-y3)+pow(lr1,-1)*(x2-
      y2)*pow(lr1+x2-y2,-1)*(x3-y3)+4*((-1)+nu)*((-1)+2*
      nu)*pow(lr2,-1)*pow(x2-y2,2)*pow(lr2+x3+y3,-1)+pow(lr2+x2-y2,
      -1)*((7+8*((-2)+nu)*nu)*x3+y3+8*((-1)+nu)*nu*y3)+pow(lr2,
      -1)*(x2-y2)*pow(lr2+x2-y2,-1)*((7+8*((-2)+nu)*nu)
      *x3+y3+8*((-1)+nu)*nu*y3)-4*((-1)+nu)*((-1)+2*nu)*lr2*pow(
      x1-y1,2)*pow(pow(x1-y1,2)+pow(x2-y2,2),-1)*(x3+
      y3)*pow(pow(x1-y1,2)+pow(x3+y3,2),-1)+4*((-1)+nu)*((-1)+2*
      nu)*pow(lr2,-1)*pow(x1-y1,2)*pow(pow(x1-y1,2)+pow(x2-y2,
      2),-1)*pow(x2-y2,2)*(x3+y3)*pow(pow(x1-y1,2)+pow(x3+y3,
      2),-1)+2*pow(lr2,-1)*x3*y3*(x3+y3)*pow(pow(x1-y1,2)+pow(x3+y3,
      2),-1)-2*x3*pow(x2-y2,2)*y3*(x3+y3)*pow(pow(x1-y1,
      2)+pow(x3+y3,2),-1)*pow(pow(x1-y1,2)+pow(x2-y2,2)+pow(x3+y3,
      2),-3/2)+4*((-1)+nu)*((-1)+2*nu)*log(lr2+x3+y3));
}

double s12::J1313d2(double y1,double y2,double y3, double nu,double G,
    double x1,double x2,double x3){
  double lr1=s12::r1(x1,x2,x3, y1,y2,y3);
  double lr2=s12::r2(x2,x2,x3, y1,y2,y3);
  return (-1/16)*pow(1-nu,-1)*pow(M_PI,-1)*(pow(lr1,-1)*(x2-y2)+2*
      (7+8*((-2)+nu)*nu)*pow(lr2,-1)*(x2-y2)-2*pow(lr2,-1)*(x2+
      (-1)*y2)*pow(lr2+x3+y3,-1)*(3*x3+2*y3-6*nu*(x3+y3)+4*
      pow(nu,2)*(x3+y3))-2*((-3)+4*nu)*x3*pow(pow(x1-y1,2)+pow(x2+(-1)
      *y2,2),-1)*(x2-y2)*(x3+y3)*pow(pow(x1-y1,2)+pow(x2+(
      -1)*y2,2)+pow(x3+y3,2),-1/2)+(x2-y2)*pow(pow(x1-y1,2)+pow(
      x2-y2,2)+pow(x3+y3,2),-3/2)*(((-7)-8*((-2)+nu)*nu)*
      pow(x1,2)+2*(7+8*((-2)+nu)*nu)*x1*y1+((-7)-8*((-2)+nu)*nu)*
      pow(y1,2)-7*(pow(x3,2)+pow(x2-y2,2))-16*x3*y3-7*pow(y3,2)+(
      -8)*((-2)+nu)*nu*(pow(x2-y2,2)+pow(x3+y3,2))))*pow(G,-1);
}

double s12::J1323d2(double y1,double y2,double y3, double nu,double G,
    double x1,double x2,double x3){
  double lr1=s12::r1(x1,x2,x3, y1,y2,y3);
  double lr2=s12::r2(x1,x2,x3, y1,y2,y3);
  return (-1/16)*pow(1-nu,-1)*pow(M_PI,-1)*((x1-y1)*pow(lr1+x2+(-1)
      *y2,-1)+pow(lr1,-1)*(x1-y1)*(x2-y2)*pow(lr1+x2-
      y2,-1)+(7+8*((-2)+nu)*nu)*(x1-y1)*pow(lr2+x2-y2,
      -1)+(7+8*((-2)+nu)*nu)*pow(lr2,-1)*(x1-y1)*(x2-y2)*
      pow(lr2+x2-y2,-1)-4*((-1)+nu)*((-1)+2*nu)*(x1-
      y1)*pow(pow(x1-y1,2)+pow(x2-y2,2),-1)*(x3+y3)+2*pow(lr2,
      -1)*x3*(-x1+y1)*y3*pow(pow(x1-y1,2)+pow(x3+y3,2),-1)+
      4*((-1)+nu)*lr2*(x1-y1)*pow(pow(x1-y1,2)+pow(x2-y2,
      2),-1)*(x3+y3)*((-3)*x3-y3+2*nu*(x3+y3))*pow(pow(x1+(-1)
      *y1,2)+pow(x3+y3,2),-1)-4*((-1)+nu)*pow(lr2,-1)*(x1-
      y1)*pow(pow(x1-y1,2)+pow(x2-y2,2),-1)*pow(x2-y2,2)*
      (x3+y3)*((-3)*x3-y3+2*nu*(x3+y3))*pow(pow(x1-y1,2)+pow(
      x3+y3,2),-1)+2*x3*(x1-y1)*pow(x2-y2,2)*y3*pow(pow(x1+
      (-1)*y1,2)+pow(x3+y3,2),-1)*pow(pow(x1-y1,2)+pow(x2-y2,
      2)+pow(x3+y3,2),-3/2))*pow(G,-1);
}

double s12::J2312d1(double y1,double y2,double y3, double nu,double G,
    double x1,double x2,double x3){
  double lr1=s12::r1(x1,x2,x3, y1,y2,y3);
  double lr2=s12::r2(x1,x2,x3, y1,y2,y3);
  return (-1/16)*pow(1-nu,-1)*pow(M_PI,-1)*pow(G,-1)*(4*((-1)+nu)*((
      -1)+2*nu)*pow(pow(x1-y1,2)+pow(x2-y2,2),-1)*pow(x2-
      y2,2)+pow(lr1+x1-y1,-1)*(x3-y3)+pow(lr1,-1)*(x1-
      y1)*pow(lr1+x1-y1,-1)*(x3-y3)+4*((-1)+nu)*((-1)+2*
      nu)*pow(lr2,-1)*pow(x1-y1,2)*pow(lr2+x3+y3,-1)+pow(lr2+x1-y1,
      -1)*((7+8*((-2)+nu)*nu)*x3+y3+8*((-1)+nu)*nu*y3)+pow(lr2,
      -1)*(x1-y1)*pow(lr2+x1-y1,-1)*((7+8*((-2)+nu)*nu)
      *x3+y3+8*((-1)+nu)*nu*y3)-4*((-1)+nu)*((-1)+2*nu)*lr2*pow(
      pow(x1-y1,2)+pow(x2-y2,2),-1)*pow(x2-y2,2)*(x3+
      y3)*pow(pow(x2-y2,2)+pow(x3+y3,2),-1)+4*((-1)+nu)*((-1)+2*
      nu)*pow(lr2,-1)*pow(x1-y1,2)*pow(pow(x1-y1,2)+pow(x2-y2,
      2),-1)*pow(x2-y2,2)*(x3+y3)*pow(pow(x2-y2,2)+pow(x3+y3,
      2),-1)+2*pow(lr2,-1)*x3*y3*(x3+y3)*pow(pow(x2-y2,2)+pow(x3+y3,
      2),-1)-2*x3*pow(x1-y1,2)*y3*(x3+y3)*pow(pow(x2-y2,
      2)+(x3+y3,2),-1)*pow(pow(x1-y1,2)+pow(x2-y2,2)+pow(x3+y3,
      2),-3/2)+(4*((-1)+nu)*((-1)+2*nu))*log(lr2+x3+y3));
}

double s12::J2313d1(double y1,double y2,double y3, double nu,double G,
    double x1,double x2,double x3){
  double lr1=s12::r1(x1,x2,x3, y1,y2,y3);
  double lr2=s12::r2(x1,x2,x3, y1,y2,y3);
  return (-1/16)*pow(1-nu,-1)*pow(M_PI,-1)*(pow(lr1+x1-y1,-1)*(
      x2-y2)+pow(lr1,-1)*(x1-y1)*pow(lr1+x1-y1,-1)*(x2+
      (-1)*y2)+(7+8*((-2)+nu)*nu)*pow(lr2+x1-y1,-1)*(x2-
      y2)+(7+8*((-2)+nu)*nu)*pow(lr2,-1)*(x1-y1)*pow(lr2+x1-
      y1,-1)*(x2-y2)-4*((-1)+nu)*((-1)+2*nu)*pow(pow(x1+(-1)
      *y1,2)+pow(x2-y2,2),-1)*(x2-y2)*(x3+y3)+2*pow(lr2,
      -1)*x3*(-x2+y2)*y3*pow(pow(x2-y2,2)+pow(x3+y3,2),-1)+
      4*((-1)+nu)*lr2*pow(pow(x1-y1,2)+pow(x2-y2,2),-1)*(x2+(
      -1)*y2)*(x3+y3)*((-3)*x3-y3+2*nu*(x3+y3))*pow(pow(x2-
      y2,2)+pow(x3+y3,2),-1)-4*((-1)+nu)*pow(lr2,-1)*pow(x1-y1,
      2)*pow(pow(x1-y1,2)+pow(x2-y2,2),-1)*(x2-y2)*(
      x3+y3)*((-3)*x3-y3+2*nu*(x3+y3))*pow(pow(x2-y2,2)+pow(x3+
      y3,2),-1)+2*x3*pow(x1-y1,2)*(x2-y2)*y3*pow(pow(x2+(
      -1)*y2,2)+pow(x3+y3,2),-1)*pow(pow(x1-y1,2)+pow(x2-y2,2)+
      pow(x3+y3,2),-3/2))*pow(G,-1);
}

double s12::J2323d1(double y1,double y2,double y3, double nu,double G,
    double x1,double x2, double x3){
  double lr1=s12::r1(x1,x2,x3, y1,y2,y3);
  double lr2=s12::r2(x1,x2,x3, y1,y2,y3);
  return (-1/16)*pow(1-nu,-1)*pow(M_PI,-1)*(pow(lr1,-1)*(x1-y1)+2*
      (7+8*((-2)+nu)*nu)*pow(lr2,-1)*(x1-y1)-2*pow(lr2,-1)*(x1+
      (-1)*y1)*pow(lr2+x3+y3,-1)*(3*x3+2*y3-6*nu*(x3+y3)+4*
      pow(nu,2)*(x3+y3))-2*((-3)+4*nu)*x3*(x1-y1)*pow(pow(x1-
      y1,2)+pow(x2-y2,2),-1)*(x3+y3)*pow(pow(x1-y1,2)+pow(x2+(
      -1)*y2,2)+pow(x3+y3,2),-1/2)+(x1-y1)*pow(pow(x1-y1,2)+pow(
      x2-y2,2)+pow(x3+y3,2),-3/2)*(((-7)-8*((-2)+nu)*nu)*
      pow(x1,2)+2*(7+8*((-2)+nu)*nu)*x1*y1+((-7)-8*((-2)+nu)*nu)*
      pow(y1,2)-7*(pow(x3,2)+pow(x2-y2,2))-16*x3*y3-7*pow(y3,2)+(
      -8)*((-2)+nu)*nu*(pow(x2-y2,2)+pow(x3+y3,2))))*pow(G,-1);
}
