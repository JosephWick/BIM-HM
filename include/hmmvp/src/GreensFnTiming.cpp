

class GreensFnTest : public ImplGreensFn {
public:
  virtual void Init(const KeyValueFile* kvf) throw (Exception);
  virtual Hd* ComputeHd (double eta) { return NewHd(_x, _x, NULL, eta); }
  virtual bool Call(const CompressBlockInfo& cbi, const vector<UInt>& rs,
                    const vector<UInt>& cs, double* B) const;

private:
  // geometry
  Matd _x;

  double _N;

  double Eval(UInt i, UInt j) const;
};

inline double GreensFnTest::Eval (UInt i, UInt j) const {
  // i is the receiver, j is the source; both start at 1
  // keep in mind that i/j are the cell number not location
  //printf("ij: %d, %d\n", i, j);

  // take distance between points i and j
  // if dist = 0, kernel value is negative
  // if dist > 0, kernel value is 1/d

  // declaration

  double dist;

  dist = sqrt( pow(_x(1,i)-_y(1,j),2) +
                pow(_x(2,i)-_y(2,j),2) +
                pow(_x(3,i)-_y(3,j),2) );

  if (dist == 0){
    return -2;
  }
  return 1/dist;

}

void GreensFnTest::Init (const KeyValueFile* kvf) throw (Exception) {
  const Matd* m;
  const Matd* n;
  const Matd* o;

  const Matd* l;
  const Matd* w;

  if (!kvf->GetMatd("X", m)) throw Exception("Missing X.");
  _x = *m;
  if (_x.Size(1) != 3) throw Exception("X must be 3xN.");

  kvf->GetDouble("N", _N);
  if (_N <= 0) throw Exception("N must be greater than 0.");
}

bool GreensFnTest::
Call (const CompressBlockInfo& cbi, const vector<UInt>& rs,
      const vector<UInt>& cs, double* B) const {
  for (UInt k = 0, ic = 0; ic < cs.size(); ic++)
    for (UInt ir = 0; ir < rs.size(); ir++, k++)
      B[k] = Eval(rs[ir], cs[ic]);
  return true;
}
