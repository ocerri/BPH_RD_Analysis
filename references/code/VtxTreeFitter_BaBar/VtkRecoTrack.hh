#ifndef _VTK_RECOTRACK_HH__
#define _VTK_RECOTRACK_HH__

#include "VtxTreeFitter/VtkRecoParticle.hh"

class HepVector ;
class HepSymMatrix ;
class HepMatrix ;
class BField ;
class TrkFit ;

namespace vtxtreefit
{

  class RecoTrack : public RecoParticle
  {
  public:
    RecoTrack(const BtaCandidate* bc, const ParticleBase* mother) ;
    virtual ~RecoTrack() ;

    virtual ErrCode initPar2(FitParams*) ;
    virtual ErrCode initCov(FitParams*) const ;
    virtual int dimM() const { return 5 ; }
    virtual int type() const { return kRecoTrack ; }

    virtual ErrCode projectRecoConstraint(const FitParams&, Projection&) const ;
    ErrCode updCache(double flt) ;
    static void setApplyCovCorrection(bool b=true) { gApplyCovCorrection = b ; }
    
    static void correctCov(HepSymMatrix& V) ;
    
    virtual int nFinalChargedCandidates() const { return 1 ; }
    
    virtual void addToConstraintList(constraintlist& alist, int depth) const {
      alist.push_back(Constraint(this,Constraint::track,depth,dimM()) ) ; 
    }
    ErrCode updFltToMother(const FitParams& fitparams) ;
    void setFlightLength(double flt) { _flt = flt ; }
    const TrkFit* trkFit() const { return _trkFit ; }
  private:
    static bool   gApplyCovCorrection ;
    static double gCovCorrection[5] ;
    const BField* _bfield ;
    const TrkFit* _trkFit ;
    bool _cached ;
    double _flt ;
    HepVector    _m ;
    HepSymMatrix _matrixV ;
  } ;

}
#endif
