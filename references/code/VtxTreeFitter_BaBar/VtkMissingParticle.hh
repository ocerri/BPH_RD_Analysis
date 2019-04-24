#ifndef __VTK_MISSINGPARTICLE_HH__
#define __VTK_MISSINGPARTICLE_HH__

#include "VtxTreeFitter/VtkParticleBase.hh"

class HepVector ;
class HepSymMatrix ;
class HepMatrix ;

namespace vtxtreefit
{

  class MissingParticle : public ParticleBase
  {
  public:
    MissingParticle(const BtaCandidate* bc, const ParticleBase* mother) ;
    virtual ~MissingParticle() ;

    virtual ErrCode initPar1(FitParams*) ;
    virtual ErrCode initPar2(FitParams*) { return ErrCode::success ; }

    virtual std::string parname(int index) const ;
    virtual int dim() const  { return _constrainMass ? 3 : 4; }   //(px,py,pz,E)
    virtual int momIndex() const { return index() ; }
    virtual bool hasEnergy() const { return _constrainMass ? false : true ; }
    virtual int type() const { return kMissingParticle ; }
    virtual void addToConstraintList(constraintlist& alist, int depth) const {
      //if(_constrainMass) alist.push_back(Constraint(this,Constraint::mass,depth)) ; 
    }
  private:
    bool _constrainMass ;
  } ;

}
#endif
