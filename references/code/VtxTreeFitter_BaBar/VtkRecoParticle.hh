#ifndef __VTK_RECOPARTICLE_HH__
#define __VTK_RECOPARTICLE_HH__

#include "VtxTreeFitter/VtkParticleBase.hh"

class HepVector ;
class HepSymMatrix ;
class HepMatrix ;

namespace vtxtreefit
{

  class RecoParticle : public ParticleBase
  {
  public:
    RecoParticle(const BtaCandidate* bc, const ParticleBase* mother) ;
    virtual ~RecoParticle() ;

    virtual int dimM() const = 0; // dimension of the measurement    
    virtual ErrCode initPar1(FitParams*) { return ErrCode::success ; } 
    //virtual ErrCode initCov(FitParams*) const ;
    virtual std::string parname(int index) const ;
    virtual int dim() const { return 3; }   //(px,py,pz)
 
    virtual int momIndex() const { return index() ; }
    virtual bool hasEnergy() const { return false ; }

    virtual ErrCode projectRecoConstraint(const FitParams& fitparams, Projection& p) const = 0 ;
    virtual ErrCode projectConstraint(Constraint::Type, const FitParams&, Projection&) const ;
    virtual double chiSquare(const FitParams* fitparams) const ;

  protected:
    //     // a few variables that are worth caching
    //     HepVector* _m ;    // measurement vector
    //     HepSymMatrix* _V ; // variance in measurement
    //     HepSymMatrix* _W ; // inverse of _V
  private:
  } ;

}
#endif
