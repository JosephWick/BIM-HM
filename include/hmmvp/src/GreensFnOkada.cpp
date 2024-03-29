/* GreensFnOkada
 * Joseph Wick, 7/30/2021
 *
 * Uses the fortran script provided in dc3dm/external to fill in an h-matrix with
 * Okada solutions
 *
*/


#include <iostream>
#include <cmath>
#define _USE_MATH_DEFINES

extern "C"{
void dc3d0_(char* SPACE, double* ALPHA, double* X, double* Y, double* Z,
              double* DEPTH, double* DIP, double* AL1, double* AL2, double* AW1,
              double* AW2, double* DISL1, double* DISL2, double* DISL3, double* UX,
              double* UY, double* UZ, double* UXX, double* UYX, double* UZX,
              double* UXY, double* UYY, double* UZY, double* UXZ, double* UYZ,
              double* UZZ);
}


class GreensFnOkada : public ImplGreensFn {
public:
  virtual void Init(const KeyValueFile* kvf) throw (Exception);
  virtual Hd* ComputeHd (double eta) {return NewHd(_x, _x, NULL, eta); }
  virtual bool Call(const CompressBlockInfo& cbi, const vector<UInt>& rs,
                    const vector<UInt>& cs, double* B) const;

private:
  //geometry
  Matd _x;

  // halfspace
  char _h;

  // elastic properties
  double _mu;
  double _nu;
  double _alpha;

  // size of blocks
  double _dz;

  // angle of fault
  double _dip;

  double _depth;

  // length/width of fault
  double _L;
  double _W;

  // dislocations
  double _d1;
  double _d2;
  double _d3;

  double Eval(UInt i, UInt j) const;
};

inline double GreensFnOkada::Eval (UInt i, UInt j) const {
  //find individual values to the hmatrix here
  // i is the receiver, j is the source

  // these correspond to AL1 and AW1 respectively
  double zL = 0;
  double zW = 0;

  // for source depth; source measured from top left
  double srcdepth = (double)_x(2,j);

  // getting parameters
  // for observer position; obs/rec is measured from center
  double obsx = (double)_x(1,i) - 0.5*_dz - (double)_x(1,j);
  double obsy = (double)_x(2,i) - 0.5*_dz - (double)_x(2,j); // rel to src depth
  double obsz = (double)_x(3,i) - 0.5*_dz - (double)_x(3,j);
  //printf("dz: %f, i: %d, _x: %f, obsy: %f\n", _dz, i, _x(2,i), obsy);

  // pointer business
  double tmp1 = 0.5*_L;
  double tmp2 = 0.5*_W;
  char * ph = const_cast<char*>(&_h);
  double * pAlpha = const_cast<double*>(&_alpha);
  double * pDip = const_cast<double*>(&_dip);
  double * pL = const_cast<double*>(&_L);
  double * pW = const_cast<double*>(&_W);
  double * pD1 = const_cast<double*>(&_d1);
  double * pD2 = const_cast<double*>(&_d2);
  double * pD3 = const_cast<double*>(&_d3);

  // outputs
  double ux;
  double uy;
  double uz;
  double uxx;
  double uyx;
  double uzx;
  double uxy;
  double uyy;
  double uzy;
  double uxz;
  double uyz;
  double uzz;

  //printf("\n===============\n");
  //printf("alpha: %f\n", *pAlpha);
  //printf("obsx: %f, obsy: %f, obsz: %f\n", obsx, obsy, obsz);
  //printf("srcdepth: %f, dip: %f\n", srcdepth, _dip);
  //printf("d1: %f, d2: %f, d3: %f\n", *pD1, *pD2, *pD3);
  //printf("================\n");

  /*
  dc3d0_(ph, pAlpha,
        &obsx, &obsy, &obsz,
        &srcdepth, pDip,
        &zL, pL, &zW, pW,
        pD1, pD2, pD3,
        &ux, &uy, &uz, &uxx, &uyx, &uzx, &uxy, &uyy, &uzy, &uxz, &uyz, &uzz);

  double out = _mu*(uxy + uyx);
  printf("out: %f\n", out);

  return out;
  */
  return 1;
}

void GreensFnOkada::Init(const KeyValueFile* kvf) throw (Exception) {
  const Matd* m;
  double d;
  double tmp;

  if (!kvf->GetMatd("X", m)) throw Exception("Missing X.");
  _x = *m;
  if (_x.Size(1) != 3) throw Exception("X must be 3xN.");

  _h = 'f';
  kvf->GetDouble("halfspace", tmp);
  if (tmp == 1.0) _h = 'h';
  printf("h: %c\n", _h);

  kvf->GetDouble("mu", _mu);
  if (_mu < 1) throw Exception("mu must be greater than 0.");
  printf("mu: %f\n", _mu);

  kvf->GetDouble("nu", _nu);
  if (_nu <=0) throw Exception("nu must be greater than 0.");
  printf("nu: %f\n", _nu);

  tmp = (2*_mu*_nu)/(1-2*_nu);
  _alpha = (tmp+_mu)/(tmp+2*_mu);
  printf("alpha: %f\n", _alpha);

  kvf->GetDouble("dz", _dz);
  if (_dz <=0) throw Exception("dz must be greater than 0.");
  printf("dz: %f\n", _dz);

  kvf->GetDouble("dip", _dip);
  if (_dip <0 || _dip>360) throw Exception("dip must be 0<dip<360.");
  printf("dip: %f\n", _dip);

  kvf->GetDouble("depth", _depth);
  printf("depth: %f\n", _depth);

  kvf->GetDouble("L", _L);
  if (_L<0) throw Exception("L must be greater than 0.");
  printf("L: %f\n", _L);

  kvf->GetDouble("W", _W);
  if (_W < 0 ) throw Exception("W must be greater than 0.");
  printf("W: %f\n", _W);

  kvf->GetDouble("d1", _d1);
  kvf->GetDouble("d2", _d2);
  kvf->GetDouble("d3", _d3);
  printf("d1: %f, d2: %f, d3: %f", _d1, _d2, _d3);

}

bool GreensFnOkada::
Call (const CompressBlockInfo& cbi, const vector<UInt>& rs,
      const vector<UInt>& cs, double* B) const {
  for (UInt k = 0, ic = 0; ic < cs.size(); ic++)
    for (UInt ir = 0; ir < rs.size(); ir++, k++)
      B[k] = Eval(rs[ir], cs[ic]);
  return true;
}
