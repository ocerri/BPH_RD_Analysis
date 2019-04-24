#ifndef __VTK_FITTER_HH__
#define __VTK_FITTER_HH__

class BtaCandidate ;
class HepSymMatrix ;
class HepVector ;
class BbrDoubleErr ;
class BbrLorentzVectorErr ;
class BtaFitParams ;

#include <vector>
#include "VtxTreeFitter/VtkErrCode.hh"

namespace vtxtreefit
{

  class DecayChain ;
  class FitParams ;
  class Upsilon ;
  class ParticleBase ;
  
  extern int vtxverbose ;

  class Fitter
  {
  public:
    Fitter() : _decaychain(0), _fitparams(0) {} 
    Fitter(const BtaCandidate& bc, double prec=0.01) ;
    ~Fitter() ;
    void fit() ;
    void print() const ;
    void printConstraints(std::ostream& os=std::cout) const ;
    const HepSymMatrix& cov() const ;
    const HepVector& par() const ;
    HepSymMatrix cov(const std::vector<int>& indexVec) const ;
    HepVector par(const std::vector<int>& indexVec) const ;
    //const DecayChain* decayChain() const { return _decaychain; }
    int posIndex(const BtaCandidate* bc) const ;
    int momIndex(const BtaCandidate* bc) const ;
    int tauIndex(const BtaCandidate* bc) const ;
    
    double chiSquare() const { return _chiSquare ; }
    double globalChiSquare() const ;
    int    nDof()      const ; 
    int status() const { return _status ; }
    int nIter() const { return _niter ; }
    const ErrCode& errCode() { return _errCode ; }
    
    // must be moved to derived class or so ...
    double add(const BtaCandidate& cand) ; 
    double remove(const BtaCandidate& cand) ;
    void updateIndex() ;
    void fitOneStep() ;

    // interface to beta
    BbrDoubleErr decayLength(const BtaCandidate& cand) const ;
    BbrDoubleErr lifeTime(const BtaCandidate& cand) const ;
    BbrDoubleErr decayLengthSum(const BtaCandidate&,const BtaCandidate&) const ;

    BtaCandidate getFitted() const ;
    BtaCandidate getFitted(const BtaCandidate& cand) const ;
    BtaCandidate getFittedTree() const ;
    BtaCandidate* fittedCand(const BtaCandidate& cand, BtaCandidate* headoftree) const ;

    BtaFitParams btaFitParams(const BtaCandidate& cand) const ;

    void updateCand(BtaCandidate& cand) const ;
    void updateTree(BtaCandidate& cand) const ;

    static void setVerbose(int i) { vtxverbose = i ; }

  public:
    BtaFitParams btaFitParams(const ParticleBase* pb) const ;
    BbrDoubleErr decayLength(const ParticleBase* pb) const ;
    BbrDoubleErr decayLengthSum(const ParticleBase*,const ParticleBase*) const ;

    DecayChain* decaychain() { return _decaychain ; }
    FitParams* fitparams() { return _fitparams ; }
    const DecayChain* decaychain() const { return _decaychain ; }
    const FitParams* fitparams() const { return _fitparams ; }
    const BtaCandidate* bc() const { return _bc ; }
    static BbrDoubleErr decayLength(const ParticleBase* pb, const FitParams* ) ;
  private:
    const BtaCandidate* _bc ;
    DecayChain* _decaychain ;
    FitParams* _fitparams ;
    int _status ;
    double _chiSquare ;
    int _niter ;
    double _prec ;
    ErrCode _errCode ;
  } ;

}

#endif
