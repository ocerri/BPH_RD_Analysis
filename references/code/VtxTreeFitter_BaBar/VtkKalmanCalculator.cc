#include "BaBar/BaBar.hh"
#include <ErrLogger/ErrLog.hh>
#include "VtxTreeFitter/VtkKalmanCalculator.hh"

#define SLOWBUTSAFE 1
#undef SLOWBUTSAFE
#undef SKIPHIGHACCURACYCORRECTION

namespace vtxtreefit
{
  
  inline double fastsymmatrixaccess(double* m, int row, int col)
  {
    return *(m+(row*(row-1))/2+(col-1));
  }

  inline double symmatrixaccess(double* m, int row, int col)
  {
    return (row>=col? fastsymmatrixaccess(m,row,col) : fastsymmatrixaccess(m,col,row)) ;
  }

  ErrCode
  KalmanCalculator::init(const HepVector& value, const HepMatrix& G, 
			 const FitParams* fitparams, const HepSymMatrix* V, 
			 int weight)
  {
    ErrCode status ;
    _nconstraints = value.num_row() ;  // dimension of the constraint
    _nparameters  = fitparams->dim() ; // dimension of the state

    int valdim  = value.num_row() ; // dimension of the constraint
    int statdim = fitparams->par().num_row() ; // dimension of the state

#ifdef VTK_BOUNDSCHECKING
    assert( G.num_row() == valdim && G.num_col() == statdim &&
	    (!V || V->num_row()==valdim) ) ;
#endif
    _value = &value ;
    _matrixG     = &G ;
    const HepSymMatrix& C = fitparams->cov() ;
    // calculate C*G.T()
#ifdef SLOWBUTSAFE
    _matrixCGT = C * G.T() ;
#else
    double tmp ;
    _matrixCGT = HepMatrix(statdim,valdim,0) ;
    for(int col=1; col<=_nconstraints; ++col)
      for(int k=1; k<=_nparameters; ++k)
 	if( (tmp=G(col,k)) !=0 ) {
	  for(int row=1; row<k; ++row)
	    _matrixCGT(row,col) += C.fast(k,row) * tmp ;
	  for(int row=k; row<=statdim; ++row)
	    _matrixCGT(row,col) += C.fast(row,k) * tmp ;
	}
#endif

    // calculate the error in the predicted residual R = G*C*GT + V
    // slow:
#ifdef SLOWBUTSAFE
    _matrixRinv = fitparams->cov().similarity(G) ;
    if(V) _matrixRinv += weight*(*V) ;
#else
    if(V) {
      _matrixRinv = *V ;
      if(weight!=1) _matrixRinv *= weight ;
    } else _matrixRinv = HepSymMatrix(valdim,0) ;
    
    for(int row=1; row<=_nconstraints; ++row)
      for(int k=1; k<=_nparameters; ++k)
	if( (tmp=G(row,k)) != 0 )
	  for(int col=1; col<=row; ++col)
	    _matrixRinv.fast(row,col) += tmp*_matrixCGT(k,col) ;
#endif    
    _matrixR = _matrixRinv ;
    _matrixRinv.invert(_ierr) ;
    if(_ierr) {
      status |= ErrCode::inversionerror; 
      ErrMsg(warning) << "Error inverting matrix. Vertex fit fails." << endmsg ;
    }
    
    // calculate the gain matrix
    _matrixK = _matrixCGT * _matrixRinv ;
    _chisq = -1 ;
//     // let's see if we get same results using sparce matrices
//     VtkSparseMatrix Gs(G) ; 
//     VtkSparseMatrix CGT = Gs.transposeAndMultiplyRight(fitparams->cov()) ;
//     HepSymMatrix Rs(value.numrow()) ;
//     Gs.multiplyLeft(CGT,Rs) ;
//     if(V) Rs += (*V) ;
//     Rs.invert(_ierr) ;
//     VtkSparseMatrix Ks = CGT*Rs ;
    return status ;
  }

  void 
  KalmanCalculator::updatePar(FitParams* fitparams)
  {
    //fitparams->par() -= fitparams->cov() * (G.T() * (R * value) ) ;
    fitparams->par() -= _matrixK * (*_value) ;
    _chisq = _matrixRinv.similarity(*_value) ;
  }

  void 
  KalmanCalculator::updatePar(const HepVector& pred, FitParams* fitparams)
  {
    // this is still very, very slow !
    HepVector valueprime = (*_value) + (*_matrixG) * (pred-fitparams->par()) ;
    fitparams->par() = pred - _matrixK*valueprime ;
    _chisq = _matrixRinv.similarity( valueprime ) ;
  }
  
  void 
  KalmanCalculator::updateCov(FitParams* fitparams, double chisq)
  {

#ifdef SLOWBUTSAFE
    HepSymMatrix deltaCov = _matrixRinv.similarityT(*_matrixG).similarity(fitparams->cov()) ;
    fitparams->cov() -= deltaCov ;
#else

    // There are two expessions for updating the covariance
    // matrix.
    // slow: deltaCov = - 2*C*GT*KT +  K*R*KT
    // fast: deltaCov = - C*GT*KT
    // The fast expression is very sensitive to machine accuracy. The
    // following piece of code effectively invokes the slow
    // expression. I couldn't write it faster than this.

    double tmp ;
#ifndef SKIPHIGHACCURACYCORRECTION
    // substitute C*GT --> 2*C*GT - K*R. of course, this invalidates
    // C*GT, but we do not need it after this routine.
    
    // we use the fact that _in principle_ C*GT = K*R, such that
    // they have the same zero elements
    for(int row=1; row<=_nparameters; ++row)
      for(int col=1; col<=_nconstraints; ++col) 
	if( (tmp =2*_matrixCGT(row,col))!=0 ) {
	  for(int k=1; k<=_nconstraints; ++k)
	    tmp -= _matrixK(row,k) * _matrixR(k,col) ;
	  _matrixCGT(row,col) = tmp ;
	}
#endif
    
//     HepMatrix KR = _matrixK*_matrixR ;
//     double tmp ;
//     for(int row=1; row<=_nparameters; ++row)
//       for(int k=1; k<=_nconstraints; ++k) 
// 	if( (tmp= (KR(row,k) - 2*_matrixCGT(row,k))) != 0 )
// 	  for(int col=1; col<=row; ++col) 
// 	    fitparams->cov().fast(row,col) += tmp * _matrixK(col,k) ;

    // deltaCov = - C*GT*KT 
    for(int row=1; row<=_nparameters; ++row)
      for(int k=1; k<=_nconstraints; ++k) 
	if( (tmp = -(_matrixCGT(row,k))) != 0 )  // they have same size, and same 'emptiness'
	  for(int col=1; col<=row; ++col) 
	    fitparams->cov().fast(row,col) += tmp * _matrixK(col,k) ;

#endif
    fitparams->addChiSquare(chisq>0 ? chisq : _chisq, _value->num_row()) ;
    for(int col=1; col<=_nconstraints; ++col)
      for(int k=1; k<=_nparameters; ++k)
 	if( (*_matrixG)(col,k) !=0 ) ++(fitparams->nConstraintsVec(k)) ;
  }
}
