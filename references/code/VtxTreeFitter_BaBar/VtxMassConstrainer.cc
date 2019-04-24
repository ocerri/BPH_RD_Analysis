#include "BaBar/BaBar.hh"
#include "VtxTreeFitter/VtxMassConstrainer.hh"
#include "VtxTreeFitter/VtkFitParams.hh"
#include "VtxTreeFitter/VtkKalmanCalculator.hh"
#include <Beta/BtaAbsVertex.hh>
#include <BbrGeom/BbrLorentzVectorErr.hh>
using std::cout;
using std::endl;

using namespace vtxtreefit ;
namespace vtxtreefit {
  extern int vtxverbose ;
}

VtxMassConstrainer::VtxMassConstrainer(const BbrLorentzVectorErr& p4in,
				       const double pdtmass)
  : _status(BtaAbsVertex::UnFitted), _fitpar(4)
{
  init( p4in, p4in.covMatrix(), pdtmass) ;
}

VtxMassConstrainer::VtxMassConstrainer(const HepLorentzVector& p4in,
				       const HepSymMatrix& cov,
				       const double pdtmass)
  : _status(BtaAbsVertex::UnFitted), _fitpar(4)
{
  init( p4in, cov, pdtmass) ;
}

void VtxMassConstrainer::init(const HepLorentzVector& p4in,
			      const HepSymMatrix& cov, const double pdtmass)
{
  const int dim(1) ;  // dimension of the constraint
  const int maxniter(5) ; // max number of iterations
  const double dchisqconverged = 0.01 ; // chi2 change for convergence

  // initialize the parameter set
  _fitpar.par(1) = p4in.x() ;
  _fitpar.par(2) = p4in.y() ;
  _fitpar.par(3) = p4in.z() ;
  _fitpar.par(4) = p4in.t() ;
  _fitpar.cov()  = cov ;

  // save the prediction
  const HepVector pred = _fitpar.par() ;

  // filter the mass constraint
  HepVector residual(dim) ;
  HepMatrix projmatrix(dim,_fitpar.dim()) ;
  KalmanCalculator kalman ;
  double chisq(0) ;
  double px,py,pz,E;

  ErrCode status ;
  _status = BtaAbsVertex::NonConverged ;
  for(_niter=0; _niter<maxniter && _status==BtaAbsVertex::NonConverged;
      ++_niter) {

    // project the constraint
    px = _fitpar.par(1) ;
    py = _fitpar.par(2) ;
    pz = _fitpar.par(3) ;
    E  = _fitpar.par(4) ;
    residual(1) = E*E-px*px-py*py-pz*pz-pdtmass*pdtmass ;
    projmatrix(1,1) = -2*px ;
    projmatrix(1,2) = -2*py ;
    projmatrix(1,3) = -2*pz ;
    projmatrix(1,4) =  2*E ;

    status |= kalman.init( residual, projmatrix, &_fitpar ) ;
    if(_niter==0) {
      kalman.updatePar( &_fitpar ) ;
    } else {
      kalman.updatePar( pred, &_fitpar ) ;
    }

    double newchisq = kalman.chisq() ;
    double dchisq = newchisq - chisq ;
    chisq = newchisq ;

    if(vtxtreefit::vtxverbose>=2) 
      cout << "VtxMassConstrainer:: niter, status, chisq: "
	   << _niter << " " << status << " " << chisq << " " << dchisq << endl ; 
    
    if( fabs(dchisq) < dchisqconverged )  // converged 
      _status = BtaAbsVertex::Success ;
    else if( _niter > 0 && dchisq>0  )   // diverging
      _status = BtaAbsVertex::Failed ;
  }
  
  // update the covariance matrix and the chisquare
  kalman.updateCov( &_fitpar ) ;
}
