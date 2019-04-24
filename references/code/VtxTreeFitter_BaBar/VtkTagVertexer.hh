#ifndef __VTK_TAGVERTEXER_HH__
#define __VTK_TAGVERTEXER_HH__

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
  class TagTrkAbsFitWrapper ;
  class VtkPairPoca ;
  class TagFastVertex ;
  class UpsilonFitter ;
  
  class TagVertexer
  {
  public:
    TagVertexer(const BtaCandidate& bc,
		const HepAList<BtaCandidate>& bTagList,
		const EventInfo* eventinfo,
		const HepAList<BtaCandidate>& mclist,
		const BtaMcAssoc* truthMap,
		const AbsEvent* absevent) ;
    ~TagVertexer() ;

    BbrDoubleErr deltaT() const ;
    BbrDoubleErr deltaTtauB() const ;
    BbrDoubleErr deltaZ() const ;
    int nTaggingTracks() const ;
    const UpsilonFitter* vertex() const { return _vertex ; }
    BtaCandidate* createUpsilonCandidate() const  ;
    

  private:
    static const unsigned int gMaxNVertices ;        // the maximum number of vertex candidates 
    static const double gMaxChisqTrackAdd ;
    static const double gMaxSeedDeltaZ ;
    static const double gMaxSeedChisq ;
    UpsilonFitter* _vertex ;
    TagFastVertex* _fastvertex ;
    BtaCandidate _recoB ;
    int _verbosity ;

    void clustererA(const std::vector<TagTrkAbsFitWrapper*>& absfitlist,
		    const std::vector<VtkPairPoca*>&  pocalist,
		    double zCP,
		    vector<TagFastVertex*>& vertexlist) ;

    void clustererB(const std::vector<TagTrkAbsFitWrapper*>& absfitlist,
		    const std::vector<VtkPairPoca*>&  pocalist,
		    double zCP,
		    vector<TagFastVertex*>& vertexlist) ;
   } ;
}

#endif
