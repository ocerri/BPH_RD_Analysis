#include "BaBar/BaBar.hh"
#include "VtxTreeFitter/VtkFitParams.hh"
#include "VtxTreeFitter/VtkMissingParticle.hh"
#include <Beta/BtaCandidate.hh>


namespace vtxtreefit
{

  extern int vtxverbose ;

  MissingParticle::MissingParticle(const BtaCandidate* bc, const ParticleBase* mother)
    : ParticleBase(bc,mother),_constrainMass(false)
  {
    // this will be one of the very few particles for which we adjust
    // the dimension if there is a constraint
    // copy constraints
    _constrainMass = bc && bc->constraint(BtaConstraint::Mass) ;
  }

  MissingParticle::~MissingParticle() {}

  ErrCode MissingParticle::initPar1(FitParams* fitpar)
  {
    // take them from the bc
    HepLorentzVector p4 = bc() ? bc()->p4() : HepLorentzVector(0,0,1,1) ;
    int momindex = momIndex();
    fitpar->par()(momindex+1) = p4.x() ;
    fitpar->par()(momindex+2) = p4.y() ;
    fitpar->par()(momindex+3) = p4.z() ;
    if(hasEnergy()) fitpar->par()(momindex+4) = p4.t() ;
    return ErrCode() ;
  }

  std::string MissingParticle::parname(int index) const
  {
    return ParticleBase::parname(index+4) ;
  }
}
