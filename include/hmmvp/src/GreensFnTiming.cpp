

class GreensFnTiming : public ImplGreensFn {
public:
  virtual void Init(const KeyValueFile* kvf) throw (Exception);
  virtual Hd* ComputeHd (double eta) { return NewHd(_k, _k, NULL, eta); }
  virtual bool Call(const CompressBlockInfo& cbi, const vector<UInt>& rs,
                    const vector<UInt>& cs, double* B) const;

private:

  // kernel
  Matd _k;

  double Eval(UInt i, UInt j) const;
};

inline double GreensFnTiming::Eval (UInt i, UInt j) const {
  // i is the receiver, j is the source; both start at 1
  // i is row, j column
  // take the kernel passed in as a parameter

  return _k(i,j);

}

void GreensFnTiming::Init (const KeyValueFile* kvf) throw (Exception) {
  const Matd* m;
  const Matd* n;
  const Matd* o;

  const Matd* l;
  const Matd* w;

  if (!kvf->GetMatd("K", n)) throw Exception("Missing K.");
  _k = *n;

}

bool GreensFnTiming::
Call (const CompressBlockInfo& cbi, const vector<UInt>& rs,
      const vector<UInt>& cs, double* B) const {
  for (UInt k = 0, ic = 0; ic < cs.size(); ic++)
    for (UInt ir = 0; ir < rs.size(); ir++, k++)
      B[k] = Eval(rs[ir], cs[ic]);
  return true;
}
