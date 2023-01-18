#include "stressVertShear_s12.h"

class GreensFnTiming : public ImplGreensFn {
public:
  virtual void Init(const KeyValueFile* kvf) throw (Exception);
  virtual Hd* ComputeHd (double eta) { return NewHd(_x, _x, NULL, eta); }
  virtual bool Call(const CompressBlockInfo& cbi, const vector<UInt>& rs,
                    const vector<UInt>& cs, double* B) const;

private:

  // src and receiver
  Matd _x;
  Matd _y;

  // mesh sizing
  Matd _L;

  double Eval(UInt i, UInt j) const;
};

inline double GreensFnTiming::Eval (UInt i, UInt j) const {
  // i is the receiver, j is the source; both start at 1
  // i is row, j column

  return s12::stressVertShear_s12( _x(1,i),_x(2,i),_x(3,i),
                                   _y(1,j),_y(2,j),_y(3,j),
                                   _L(1,j), _L(2,j), _L(3,j), 0,
                                   0, 1, 0, 0, 0, 0,
                                   30*10*10*10, 0.25);

}

void GreensFnTiming::Init (const KeyValueFile* kvf) throw (Exception) {
  const Matd* m;
  const Matd* n;

  const Matd* l;

  // mesh
  if (!kvf->GetMatd("X", m)) throw Exception("Missing X.");
  _x = *m;
  if (_x.Size(1) != 3) throw Exception("X must be 3xN.");

  if (!kvf->GetMatd("Y", n)) throw Exception("Missing Y.");
  _y = *n;
  if (_y.Size(1) != 3) throw Exception("Y must be 3xN.");

  // sizing
  if (!kvf->GetMatd("L", l)) throw Exception("Missing L.");
  _L = *n;
  if (_L.Size(1) != 3) throw Exception("L must be 3xN.");
}

bool GreensFnTiming::
Call (const CompressBlockInfo& cbi, const vector<UInt>& rs,
      const vector<UInt>& cs, double* B) const {
  for (UInt k = 0, ic = 0; ic < cs.size(); ic++)
    for (UInt ir = 0; ir < rs.size(); ir++, k++)
      B[k] = Eval(rs[ir], cs[ic]);
  return true;
}
