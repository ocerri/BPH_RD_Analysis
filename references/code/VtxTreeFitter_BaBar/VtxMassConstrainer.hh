#ifndef __VTX_MASSCONSTRAINER_HH__
#define __VTX_MASSCONSTRAINER_HH__

////////////////////////////////////////////////////////////////////
// 
// Utility class that forces a mass constraint on a four vector. It
// uses an extended Kalman filter to do so.
//
///////////////////////////////////////////////////////////////////

class BbrLorentzVectorErr ;
class HepLorentzVector ;
class HepSymMatrix ;

#include "VtxTreeFitter/VtkFitParams.hh"
#include <CLHEP/Vector/LorentzVector.h>

class VtxMassConstrainer
{
public:
  VtxMassConstrainer(const BbrLorentzVectorErr& p4in, const double pdtmass) ;
  VtxMassConstrainer(const HepLorentzVector& p4in, const HepSymMatrix& cov, const double pdtmass) ;
  ~VtxMassConstrainer() {}

  HepLorentzVector p4()    const { 
    return HepLorentzVector(_fitpar.par(1),_fitpar.par(2),_fitpar.par(3),_fitpar.par(4)) ; }
  HepSymMatrix     p4Err() const { return _fitpar.cov() ; }
  double chiSquare()       const { return _fitpar.chiSquare() ; }
  int    nDof()            const { return _fitpar.nDof() ; }
  int    status()          const { return _status; }
  int    niter()           const { return _niter ; }
private:
  void init( const HepLorentzVector& p4in, const HepSymMatrix& cov, const double pdtmass) ;
private:
  int _status ;
  vtxtreefit::FitParams _fitpar ;
  int _niter ;
} ;

#endif
