#include "BaBar/BaBar.hh"
#include <AbsEnv/AbsEnv.hh>
#include <BtaEnv/BtaEnv.hh>
#include <TrkBase/TrkCompTrk.hh>

#include "VtxTreeFitter/VtkTrkAbsFitWrapper.hh"


namespace vtxtreefit
{

  const BField* TrkAbsFitWrapper::_gbField = 0 ;
  
  TrkAbsFitWrapper::TrkAbsFitWrapper(const BtaCandidate* cand) 
    : _cand(cand),_fit(0)
  {
    if(!_cand->trkAbsFit()) {
      const BtaAbsVertex* dVtx=cand->decayVtx();
      BbrPointErr pos(dVtx->point(),dVtx->xxCov());
      BbrLorentzVectorErr p4  = cand->p4WCov() ;
      BbrVectorErr mom(p4,p4.covMatrix().sub(1,3));
      if(!_gbField) _gbField = gblEnv->getBta()->bField();
      int icharge = cand->charge()>0? 1: -1 ;
      _fit = new TrkCompTrk( pos,mom,dVtx->xpCov(),icharge,
			     dVtx->chiSquared(),dVtx->nDof(), _gbField);
    }
  }

  TrkAbsFitWrapper::TrkAbsFitWrapper(const BtaCandidate* cand,
			       const HepVector& par,
			       const HepSymMatrix& cov, 
			       double charge)
    : _cand(cand),_fit(0)
  {
    HepPoint pt(par(1),par(2),par(3)) ;
    BbrPointErr  pos(pt,cov.sub(1,3)) ;
    Hep3Vector dir(par(4),par(5),par(6)) ;
    BbrVectorErr mom(dir,cov.sub(4,6)) ;
    HepMatrix xpcov(3,3) ;
    for(int row=1; row<=3; ++row)   // position row
      for(int col=4; col<=6; ++col) // momentum column
	xpcov(row,col) = cov.fast(col,row) ;
    double chisq(1);
    int ndof(1) ;
    int icharge = charge>0? 1: -1 ;
    _fit = new TrkCompTrk(pos,mom,xpcov,icharge,chisq,ndof,_gbField) ;
  }
}
