#ifndef __VTK_KALMANCALCULATOR_HH__
#define __VTK_KALMANCALCULATOR_HH__

#include <CLHEP/Matrix/Matrix.h>
#include <CLHEP/Matrix/SymMatrix.h>
#include "VtxTreeFitter/VtkFitParams.hh"
#include "VtxTreeFitter/VtkErrCode.hh"

namespace vtxtreefit
{
  
  class KalmanCalculator
  {
  public:
    ErrCode init(const HepVector& value, const HepMatrix& G, 
		 const FitParams* fitparams, const HepSymMatrix* V=0, int weight=1) ;
    void updatePar(FitParams* fitparams) ;
    void updatePar(const HepVector& prediction, FitParams* fitparams) ;
    void updateCov(FitParams* fitparams, double chisq=-1) ;
    double chisq() const { return _chisq ; }
  private:
    int _nconstraints ; // dimension of the constraint
    int _nparameters  ; // dimension of the state
    const HepVector* _value ;
    const HepMatrix* _matrixG ;
    HepSymMatrix _matrixR;    // cov of residual
    HepSymMatrix _matrixRinv; // inverse of cov of residual
    HepMatrix _matrixK   ;    // kalman gain matrix
    double _chisq ;
    int _ierr ;
    // some temporary results
    HepMatrix _matrixCGT ;
  } ;
}

#endif
