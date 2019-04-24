#ifndef __VTK_TrkABSFITWRAPPER_HH__

#include <Beta/BtaCandidate.hh>
#include <TrkBase/TrkAbsFit.hh>

class BField ;

namespace vtxtreefit
{

  class TrkAbsFitWrapper
  {
    // utility class that creates TrkAbsFit object if it doesn't yet exist
  public:
    TrkAbsFitWrapper(const BtaCandidate* cand) ;
    TrkAbsFitWrapper(const BtaCandidate* cand, const HepVector& par,
		  const HepSymMatrix& cov, double charge) ;
    ~TrkAbsFitWrapper() { delete _fit ; }
    TrkAbsFitWrapper(const TrkAbsFitWrapper& rhs) 
      : _cand(rhs._cand), _fit(0) {
      if( rhs._fit ) _fit = new TrkCompTrk(*rhs._fit) ;
    }
    const BtaCandidate* cand() const { return _cand ; }
    const TrkAbsFit*  trkFit() const { return _fit ? _fit : _cand->trkAbsFit() ; }
    const TrkDifTraj& traj() const { return trkFit()->traj(); }
    double charge() const { return _cand->charge() ; }
  private:
    const BtaCandidate* _cand ;
    TrkCompTrk* _fit ;
    static const BField* _gbField ;
  } ;
 
}

#endif
