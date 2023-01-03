

class GreensFnTiming : public ImplGreensFn {
public:
  virtual void Init(const KeyValueFile* kvf) throw (Exception);
  virtual Hd* ComputeHd (double eta) { return NewHd(_x, _x, NULL, eta); }
  virtual bool Call(const CompressBlockInfo& cbi, const vector<UInt>& rs,
                    const vector<UInt>& cs, double* B) const;

private:

  // receiver
  Matd _X;
  // source
  Matd _Q;

  // sizing
  Matd L;
  Matd T;
  Matd W;

  double Eval(UInt i, UInt j) const;
};

inline double GreensFnTiming::Eval (UInt i, UInt j) const {
  // i is the receiver, j is the source; both start at 1
  // i is row, j column
  // take the kernel passed in as a parameter

  // printf("%f\n", _k(i,j));
  // return _k(i,j);

  return stressVertShear_s12( _X(1,i),_X(2,i),_X(3,i), _Q(1,j),_Q(2,j),_Q(3,j),
          L(j), T(j), W(j), 0,
          0, 1, 0, 0, 0, 0,
          30*10*10*10, 0.25);

}

void GreensFnTiming::Init (const KeyValueFile* kvf) throw (Exception) {
  const Matd* m;
  const Matd* n;

  const Matd* l;
  const Matd* w;

  if (!kvf->GetMatd("X", m)) throw Exception("Missing X.");
  _X = *m;
  if (_X.Size(1) != 3) throw Exception("X must be 3xN.");

  if (!kvf->GetMatd("Q", n)) throw Exception("Missing Q.");
  _Q = *n;

}

bool GreensFnTiming::
Call (const CompressBlockInfo& cbi, const vector<UInt>& rs,
      const vector<UInt>& cs, double* B) const {
  for (UInt k = 0, ic = 0; ic < cs.size(); ic++)
    for (UInt ir = 0; ir < rs.size(); ir++, k++)
      B[k] = Eval(rs[ir], cs[ic]);
  return true;
}
