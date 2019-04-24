#ifndef _VTK_RECOPHOTON_HH__
#define _VTK_RECOPHOTON_HH__

#include "VtxTreeFitter/VtkRecoParticle.hh"

namespace vtxtreefit
{

  class RecoPhoton : public RecoParticle
  {
  public:
    
    RecoPhoton(const BtaCandidate* bc, const ParticleBase* mother) ;
    virtual ~RecoPhoton() ;
    
    virtual int dimM() const { return _useEnergy ? 3 : 2 ; }
    virtual ErrCode initPar2(FitParams*) ;
    virtual ErrCode initCov(FitParams*) const ;
    virtual int type()     const { return kRecoPhoton ; }
    virtual ErrCode projectRecoConstraint(const FitParams&,Projection&) const ;
    ErrCode updCache() ;
    
    virtual void addToConstraintList(constraintlist& alist, int depth) const {
      alist.push_back( Constraint(this,Constraint::photon,depth,dimM()) ) ; }
    
    static bool useEnergy(const BtaCandidate& cand) ;

  private:
    bool _init ;
    bool _useEnergy ;
    HepVector _m ;
    HepSymMatrix _matrixV ;
  } ;

}
#endif
