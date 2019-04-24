#ifndef __VTK_CLUSTERER_HH__
#define __VTK_CLUSTERER_HH__

#include <vector>
#include <CLHEP/Alist/AList.h>
#include <Beta/BtaCandidate.hh>

class BtaCandidate ;
class EventInfo ;
class BtaMcAssoc ;
class AbsEvent ;
class BbrDoubleErr ;

namespace vtxtreefit
{
  class TrkAbsFitWrapper ;
  class VtkPairPoca ;
  class FastVertex ;
   
  class Clusterer
  {
  public:
    Clusterer(const HepAList<BtaCandidate>& inputlist,
	      const TrkAbsFitWrapper* tagcand=0,
	      const EventInfo* eventinfo=0,
	      const AbsEvent* absevent=0,
	      bool excludeconversions=false,
	      bool excludeghosts=false) ;
    ~Clusterer() ;
    
  private:
    static const unsigned int gMaxNVertices ;        // the maximum number of vertex candidates 
    static const double gMaxChisqTrackAdd ;
    static const double gMaxSeedChisq ;
    static const unsigned int gMaxCompositeCharge ;
  } ;
}

#endif
