#include "BaBar/BaBar.hh"
#include <stdio.h>

#include <Beta/BtaCandidate.hh>
#include <BetaRecoAdapter/BtaAbsRecoObject.hh>
#include <TrkBase/TrkFit.hh>
#include <TrkBase/TrkRecoTrk.hh>
#include <TrkBase/TrkPoca.hh>
#include <TrkBase/TrkDifTraj.hh>

#include "VtxTreeFitter/VtkRecoTrack.hh"
#include "VtxTreeFitter/VtkFitParams.hh"
#include "VtxTreeFitter/VtkHelixUtils.hh"
using std::cout;
using std::endl;

namespace vtxtreefit
{

  extern int vtxverbose ;

  // this is just the width of the pulls (d0,phi0,omega,z0,tandip)
  double RecoTrack::gCovCorrection[5] = {1.07,1.10,1.32,1.06,1.10} ;
  bool   RecoTrack::gApplyCovCorrection = false ;
    
  RecoTrack::RecoTrack(const BtaCandidate* cand, const ParticleBase* mother) 
    : RecoParticle(cand,mother),_bfield(0),_trkFit(0),_cached(false),
      _flt(0),_m(7),_matrixV(7)
  {
    _bfield = &(cand->recoTrk()->bField()) ;
    PdtPid::PidType pidType = PdtPid::pion ; 
    if( bc()->pdtEntry() ) pidType = bc()->pdtEntry()->pidId() ;
    _trkFit = bc()->recoTrk()->fitResult(pidType) ;
    if( _trkFit==0 ) _trkFit = bc()->recoTrk()->fitResult() ;
  }
  
  RecoTrack::~RecoTrack() {}

  ErrCode
  RecoTrack::initPar2(FitParams* fitparams)
  {
    if(_flt==0) const_cast<RecoTrack*>(this)->updFltToMother(*fitparams) ;
    Hep3Vector recoP = _trkFit->momentum(_flt) ;
    int momindex = momIndex() ;
    fitparams->par(momindex+1) = recoP.x() ;
    fitparams->par(momindex+2) = recoP.y() ;
    fitparams->par(momindex+3) = recoP.z() ;
    return ErrCode() ;
  }

  ErrCode
  RecoTrack::initCov(FitParams* fitparams) const 
  {
    // we only need a rough estimate of the covariance
    BbrError p4Err = bc()->p4Err() ;
    int momindex = momIndex() ;
    for(int row=1; row<=3; ++row)
      fitparams->cov()(momindex+row,momindex+row) = 1000*p4Err(row,row) ;
    return ErrCode() ;
  }

  ErrCode
  RecoTrack::updFltToMother(const FitParams& fitparams)
  {
    int posindexmother = mother()->posIndex() ;
    HepPoint pt(fitparams.par()(posindexmother+1),
		fitparams.par()(posindexmother+2),
		fitparams.par()(posindexmother+3)) ;
    TrkPoca poca(_trkFit->traj(),_flt,pt) ;
    _flt = poca.status().success() ? poca.flt1() : 0 ;
    // FIX ME: use helix poca to get estimate of flightlength first
    double lowrange  = _trkFit->traj().lowRange() ;
    double highrange = _trkFit->traj().hiRange() ;
    if(     _flt < lowrange)  _flt = lowrange ;
    else if(_flt > highrange) _flt = highrange ;
    return ErrCode() ;
  } ;

  ErrCode
  RecoTrack::updCache(double flt)
  {
    if(vtxverbose>=2)
      cout << "RecoTrack::updCache: " << name().c_str() 
	   << " from " << _flt << " to " << flt << endl ; 
    _flt = flt ;
    const TrkExchangePar recotrackpars = _trkFit->helix(_flt) ;
    _m   = recotrackpars.params() ;
    // FIX ME: bring z0 in the correct domain ...
    _matrixV   = recotrackpars.covariance() ;
    _cached = true ;

    //if(gApplyCovCorrection) correctCov(V) ;
    return ErrCode() ;
  }

  HepVector symdiag(const HepSymMatrix& m) {
    HepVector rc(m.num_row()) ;
    for(int i=1; i<=m.num_row(); ++i)
      rc(i) = sqrt(m.fast(i,i)) ;
    return rc ;
  }

  ErrCode
  RecoTrack::projectRecoConstraint(const FitParams& fitparams, Projection& p) const
  {
    ErrCode status ;
    // create HepVector with parameters
    HepVector vertexpars(6) ;
    int posindexmother = mother()->posIndex() ;
    for(int row=1; row<=3; ++row)
      vertexpars(row) = fitparams.par()(posindexmother+row) ;
    int momindex = momIndex() ;
    for(int row=1; row<=3; ++row)
      vertexpars(3+row) = fitparams.par()(momindex+row) ;

    // translate into trackparameters
    HepVector helixpars(6) ;
    HepMatrix jacobian(6,6) ;
    VtkHelixUtils::helixFromVertex(vertexpars,charge(),*_bfield,helixpars,jacobian) ;
    
    // get the measured track parameters at the poca to the mother
    if(!_cached) {
      RecoTrack* nonconst =  const_cast<RecoTrack*>(this) ;
      if(_flt==0) nonconst->updFltToMother(fitparams) ;
      nonconst->updCache(_flt) ;
    }

    if( vtxverbose>=5) {
      cout << "vertexpars = " << vertexpars.T() << endl ;
      cout << "pred = " << helixpars.T() << endl ;
      cout << "m   = " << _m.T() << endl ;
      cout << "sig = " << symdiag(_matrixV).T() << endl ;
    }
    
    // get the measured track parameters at the flightlength of the vertex
    // double flt = helixpars(6) ;
    //const double fltupdcut = 1000 ; //cm
    //if( fabs(flt - _flt) > fltupdcut ) 
    //  status |= const_cast<RecoTrack*>(this)->updCache(flt) ;
    
    // fill the residual and cov matrix
    for(int row=1; row<=5; ++row) {
      p.r(row) = helixpars(row) - _m(row) ;
      for(int col=1; col<=row; ++col) 
	p.Vfast(row,col) = _matrixV.fast(row,col) ;
    }

    // bring phi-residual in the correct domain ([-pi,pi])
    p.r(2) = VtkHelixUtils::phidomain(p.r(2)) ;
    // FIX ME: bring z0 residual in the correct domain --> this needs some thinking

    // calculate the full projection matrix from the jacobian
    // assumes that H is reset !
    for(int row=1; row<=5; ++row) {
      // the position
      for(int col=1; col<=3; ++col) 
	p.H(row,posindexmother+col) = jacobian(row,col) ;
      // the momentum
      for(int col=1; col<=3; ++col) 
	p.H(row,momindex+col) = jacobian(row,col+3) ;
    }

    return status ;
  }
  
  void
  RecoTrack::correctCov(HepSymMatrix& V)
  {
    if(gApplyCovCorrection)
      for(int row=1; row<=5; ++row)
	for(int col=1; col<=row; ++col)
	  V.fast(row,col) = V.fast(row,col)*gCovCorrection[row-1]*gCovCorrection[col-1] ;
  }

}
