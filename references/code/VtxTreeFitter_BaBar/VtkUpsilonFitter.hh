#ifndef __VTK_UPSILONFITTER_HH__
#define __VTK_UPSILONFITTER_HH__

#include "VtxTreeFitter/VtkFitter.hh"

class EventInfo ;

namespace vtxtreefit
{

  class ParticleBase ;

  class UpsilonFitter : public Fitter
  {
  public:
    // constructor from a B
    UpsilonFitter(const BtaCandidate& cand) ;

    // constructor from an upsilon
    UpsilonFitter(const BtaCandidate& ups,const BtaCandidate& recB) ;

    // use this one to set all constraint correctly and add missing particle
    static UpsilonFitter* createFromUpsilon(const BtaCandidate& ups,const BtaCandidate& recB, 
					    const EventInfo* eventInfo) ;
    
    Upsilon* getUpsilon() ;
    const Upsilon* getUpsilon() const ;
    BtaCandidate getTaggingB() const ;
    BbrDoubleErr deltaT() const ;
    BbrDoubleErr deltaZ() const ;
    BbrDoubleErr lifeTimeSum(bool difference=false) const ;
    void replaceTagB(Fitter& vertex, bool lifetimesumconstraint) ;
    
  private:
    const ParticleBase* _recB ;
    const ParticleBase* _tagB ;
  } ;
  
}


#endif
