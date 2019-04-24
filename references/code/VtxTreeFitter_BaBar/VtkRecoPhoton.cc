#include "BaBar/BaBar.hh"
#include <stdio.h>

#include "TrkBase/TrkFit.hh"
#include "TrkBase/TrkRecoTrk.hh"
#include "Beta/BtaCandidate.hh"
#include "AbsCalo/AbsRecoCalo.hh"
#include <ErrLogger/ErrLog.hh>

#include "VtxTreeFitter/VtkRecoPhoton.hh"
#include "VtxTreeFitter/VtkFitParams.hh"
#include "VtxTreeFitter/VtkHelixUtils.hh"

namespace vtxtreefit
{
  extern int vtxverbose ;

  bool RecoPhoton::useEnergy(const BtaCandidate& bc)
  {
    bool rc = true ;
    if( bc.pdtEntry() && bc.pdtEntry()->lundId() != PdtLund::gamma &&
	bc.pdtEntry()->lundId() != PdtLund::pi0 ) {
      rc = false ;
    }
    return rc ;
  }

  RecoPhoton::RecoPhoton(const BtaCandidate* bc, const ParticleBase* mother) 
    : RecoParticle(bc,mother),_init(false),_useEnergy(useEnergy(*bc)),
      _m(dimM()),_matrixV(dimM())
  {
    updCache() ;
  }
  
  RecoPhoton::~RecoPhoton() {}

  ErrCode
  RecoPhoton::initPar2(FitParams* fitparams)
  {
    // calculate the direction
    int posindexmother = mother()->posIndex() ;
    HepVector deltaX(3) ;
    double deltaX2(0) ;
    for(int row=1; row<=3; ++row) {
      double dx = _m(row) - fitparams->par(posindexmother+row) ;
      deltaX(row) = dx ;
      deltaX2 += dx*dx ;
    }
    
    // get the energy
    double energy = _useEnergy ? _m(4) : bc()->fitParams().e() ;
    
    // assign the momentum
    int momindex = momIndex() ;
    for(int row=1; row<=3; ++row)
      fitparams->par(momindex+row) = energy*deltaX(row)/sqrt(deltaX2)  ;
    return ErrCode() ;
  }

  ErrCode
  RecoPhoton::initCov(FitParams* fitparams) const 
  {    
    int momindex = momIndex() ;
    double varEnergy =  _useEnergy ? _matrixV.fast(4,4) : 1 ;
    const double factor = 1000;
    for(int row=1; row<=3; ++row) 
      fitparams->cov().fast(momindex+row,momindex+row) = factor * varEnergy ;
    return ErrCode() ;
  }

  ErrCode
  RecoPhoton::updCache()
  {
    const AbsRecoCalo* recoCalo=bc()->recoCalo();
    HepPoint centroid = recoCalo->position();
    double energy = recoCalo->energy() ;
    double correctedenergy = bc()->fitParams(HepPoint(0,0,0)).p() ;
    static int printit=10 ;
    if(fabs(energy-correctedenergy)>0.0001 && --printit>=0) 
      ErrMsg(warning) << "Neutral candidate energy (" << correctedenergy
		      << ") different from AbsRecoCalo energy (" << energy << "). "
		      << "Either you are applying an energy correct/smearing correction, or something is wrong. " 
		      << "Will use corrected energy. " << endmsg ;
    energy = correctedenergy ;
    _init = true ;
    _matrixV = recoCalo->errorMatrixXYZ(HepPoint(0,0,0),bc()->pdtEntry()) ;
    _m(1) = centroid.x() ;
    _m(2) = centroid.y() ;
    _m(3) = centroid.z() ;
    if(_useEnergy) _m(4) = energy ;

    if(!_useEnergy && vtxverbose>=3) {
      std::cout << "K-long candidate: " << recoCalo->dynamic_cast_IfrAbs3D() << std::endl
		<< "point is: (" << centroid.x() << "," << centroid.y() << ","
		<< centroid.z() << ")" << std::endl
		<< "cov matrix is: " << _matrixV << std::flush ;
    }

    return ErrCode() ;
  }

  ErrCode
  RecoPhoton::projectRecoConstraint(const FitParams& fitparams, Projection& p) const
  { 
    // residual of photon:
    // r(1-3) = motherpos + mu * photon momentum - cluster position
    // r(4)   = |momentum| - cluster energy
    // mu is calculated from the 'chi2-closest approach' (see below)

    ErrCode status ;

    // calculate the total momentum and energy:
    int momindex  = momIndex() ;
    HepVector mom = fitparams.par().sub(momindex+1,momindex+3) ;
    double mom2 = mom.normsq() ;
    double mass = pdtMass() ;
    double energy = sqrt(mass*mass + mom2) ;

    // calculate dX = Xc - Xmother
    int posindex  = mother()->posIndex() ;
    HepVector dX(3) ; 
    for(int row=1; row<=3; ++row)
      dX(row) = _m(row) - fitparams.par(posindex+row) ;
    
    // the constraints we will use are (dX = Xc - Xmother)
    //  I) r(1) = py * dX - px * dY + pz * dX - px * dZ
    // II) r(2) = px * dZ - pz * dX + py * dZ - pz * dY
    //III) r(3) = Ec - energy
     //
    // We will need two projection matrices:
    // a) the matrix that projects on the measurement parameters (=P)
    // b) the matrix that projects on the fit parameters (=H)
    //
   
    // create the matrix that projects the measurement in the constraint equations
    // this would all be easier if we had tensors. no such luck.
    HepMatrix P(3,4,0) ;
    P(1,1) = mom(2) + mom(3) ; P(1,2) = -mom(1) ; P(1,3) = -mom(1) ;
    P(2,1) = -mom(3) ;         P(2,2) = -mom(3) ; P(2,3) = mom(1) + mom(2) ;
    P(3,4) = 1 ;

    // now get the residual. start in four dimensions
    HepVector residual4(4);
    for(int row=1; row<=3; ++row) residual4(row) = dX(row) ;
    residual4(4) = _m(4) - energy ;

    // project out the 3D part
    HepVector       r = P*residual4 ;
    HepSymMatrix    V = _matrixV.similarity(P) ;
    
    // calculate the parameter projection matrix
    // first the 'position' part
    HepMatrix H(3,7,0) ;
    for(int irow=1; irow<=3; ++irow)
      for(int icol=1; icol<=3; ++icol) 
	H(irow,icol) = - P(irow,icol) ;

    // now the 'momentum' part
    H(1,4) = -dX(2)-dX(3) ;    H(1,5) = dX(1) ;   H(1,6) = dX(1) ;
    H(2,4) =  dX(3) ;          H(2,5) = dX(3) ;   H(2,6) = -dX(1)-dX(2) ; 
    for(int col=1; col<=3; ++col)
      H(3,3+col) = -mom(col)/energy ;

    // done. copy everything back into the 'projection' object
    int dimm = dimM() ; // is we don't use the energy, this is where it will drop out
    
    for(int row=1; row<=dimm; ++row)
      p.r(row) = r(row) ;
  
    // fill the error matrix
    for(int row=1; row<=dimm; ++row)
      for(int col=1; col<=row; ++col)
	p.Vfast(row,col) = V.fast(row,col) ;
    
    // fill the projection
    for(int row=1; row<=dimm; ++row) {
      for(int col=1; col<=3; ++col)
	p.H(row,posindex+col) = H(row,col) ;
      for(int col=1; col<=3; ++col)
	p.H(row,momindex+col) = H(row,col+3) ;
    }

    return status ;
  }
}
