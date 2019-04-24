#ifndef __VTK_BTAINTERFACE_HH__
#define __VTK_BTAINTERFACE_HH__

// this class imlements functionality that needs access to the
// BtaCandBase in a BtaCandidate: it is derived from BtaOperatorBase.

#include <Beta/BtaOperatorBase.hh>
class BbrLorentzVectorErr ;
class BbrDoubleErr ;
class BtaCandidate ;
class BtaFitParams ;

namespace vtxtreefit
{
  class ParticleBase ;

  class BtaInterface : public BtaOperatorBase
  {
  public:
    void addMissingParticle(BtaCandidate& tagB) ;
    BtaCandidate* createTree(const ParticleBase& pb) ;
  private:
  } ;

}
#endif
