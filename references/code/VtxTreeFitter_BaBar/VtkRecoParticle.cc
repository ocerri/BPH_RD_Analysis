#include "BaBar/BaBar.hh"
#include <stdio.h>

#include <Beta/BtaCandidate.hh>
#include <BetaRecoAdapter/BtaAbsRecoObject.hh>
#include <TrkBase/TrkFit.hh>
#include <TrkBase/TrkRecoTrk.hh>

#include "VtxTreeFitter/VtkRecoParticle.hh"
#include "VtxTreeFitter/VtkFitParams.hh"
#include "VtxTreeFitter/VtkHelixUtils.hh"

namespace vtxtreefit
{
  extern int vtxverbose ;
    
  RecoParticle::RecoParticle(const BtaCandidate* bc, const ParticleBase* mother) 
    : ParticleBase(bc,mother)
  {
  }
  
  RecoParticle::~RecoParticle() 
  {
//     delete _m ;
//     delete _V ;
//     delete _W ;
  }

  std::string RecoParticle::parname(int index) const
  {
    return ParticleBase::parname(index+4) ;
  }

  ErrCode
  RecoParticle::projectConstraint(Constraint::Type type, 
				  const FitParams& fitparams, 
				  Projection& p) const 
  {
    ErrCode status ;
    switch(type) {
    case Constraint::track:
    case Constraint::photon:
      status |= projectRecoConstraint(fitparams,p) ;
      break ;
    default:
      status |= ParticleBase::projectConstraint(type,fitparams,p) ;
    }
    return status ;
  }

  double RecoParticle::chiSquare(const FitParams* fitparams) const
  {
    // project
    Projection p(fitparams->dim(),dimM()) ;
    projectRecoConstraint(*fitparams,p) ;
    return p.chiSquare() ;
  }
}
