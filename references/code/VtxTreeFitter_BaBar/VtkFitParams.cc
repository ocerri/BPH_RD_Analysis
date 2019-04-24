#include "BaBar/BaBar.hh"
#include <iostream>
#include <iomanip>
#include "VtxTreeFitter/VtkFitParams.hh"
#include "VtxTreeFitter/VtkParticleBase.hh"
using std::cout;
using std::endl;
using std::setprecision;
using std::setw;

namespace vtxtreefit {

  FitParams::FitParams(int dim) 
    : _dim(dim),_par(dim,0),_cov(dim,0),_scale(dim,1),
      _chiSquare(0),_nConstraints(0),_nConstraintsVec(dim,0) {}

  FitParams::~FitParams() {}

  void FitParams::resetPar() {
    for(int row=1; row<=_dim; ++row) 
      _par(row) = 0 ;
  }
  
  void FitParams::resetCov(double scale) {
    for(int row=1; row<=_dim; ++row) {
      for(int col=1; col<row; ++col) 
	_cov.fast(row,col) = 0 ;
      _cov.fast(row,row) *= scale ;
      if(_cov.fast(row,row) < 0 ) _cov.fast(row,row)*=-1 ;
    }
    _chiSquare=0 ;
    _nConstraints=0 ;
    for(int row=1; row<=_dim; ++row)
      nConstraintsVec(row) = 0 ;
  }
  
  bool FitParams::testCov() const {
    bool okay=true ;
    for(int row=1; row<=_dim && okay; ++row) 
      okay = _cov.fast(row,row)>0 ;
    return okay ;
  }

  void FitParams::print() const {
    cout << setw(3) << "index" << setw(15) << "val" << setw(15) << "err" << endl ;
    cout << setprecision(5) ;
    for(int row=1; row<=_dim; ++row) 
      cout << setw(3) << row-1
	   << setw(15) << _par(row) 
	   << setw(15) << sqrt(_cov(row,row)) << endl ;
  } ;
  
  HepSymMatrix FitParams::cov(const std::vector<int>& indexVec) const {
    int nrow = indexVec.size() ;
    HepSymMatrix thecov(nrow,0) ;
    for(int row=1; row<=nrow; ++row)
      for(int col=1; col<=row ; ++col)
	thecov(row,col) = _cov(indexVec[row-1]+1,indexVec[col-1]+1) ;
    return thecov ;
  }

  HepVector FitParams::par(const std::vector<int>& indexVec) const {
    int nrow = indexVec.size() ;
    HepVector thepar(nrow,0) ;
    for(int row=1; row<=nrow; ++row)
      thepar(row) = _par(indexVec[row-1]+1) ;
    return thepar ;
  }

  void FitParams::resize(int newdim)
  {
    if( newdim > _dim ) {
      _dim = newdim ;
      // very expensive, but okay ...
      HepVector newpar(newdim,0) ;
      newpar.sub(1,_par);

      HepSymMatrix newcov(newdim,0) ;
      newcov.sub(1,_cov) ;
      
      //      HepVector newpar(newdim,0) ;
      //       HepSymMatrix newcov(newdim,0) ;
      //       cout << newpar << endl ;
      //       for(int row=1; row<=_dim ; ++row) {
      // 	newpar(row) = _par(row) ;
      // 	for(int col=1; col<=row; ++col)
      // 	// 	  newcov(row,col) = _cov(row,col) ;
      //       }
      //      cout << _par << " " << newpar << endl ;

      _par = newpar ;
      _cov = newcov ;
      _dim = newdim ;
      _nConstraintsVec.resize(newdim,0) ;
    }
  }

  void FitParams::copy(const FitParams& rhs, 
		       const indexmap& anindexmap) 
  {
    for(indexmap::const_iterator it = anindexmap.begin() ; 
	it != anindexmap.end(); ++it) {
      int idim =     it->first->dim() ;
      int indexrhs = it->second ;
      int indexlhs = it->first->index() ;
      
#ifdef VTK_BOUNDSCHECKING
      assert( idim + indexlhs <= dim() ) ;
      assert( idim + indexrhs <= rhs.dim() ) ;
#endif

      for(int i=1; i<=idim; ++i)
	_par(indexlhs+i) = rhs._par(indexrhs+i) ;
      
      for(indexmap::const_iterator it2 = it ; 
	  it2 != anindexmap.end(); ++it2) {
	int jdim     = it2->first->dim() ;
	int jndexrhs = it2->second ;
	int jndexlhs = it2->first->index() ;
	
	for(int i=1; i<=idim; ++i)
	  for(int j=1; j<=jdim ; ++j)
	    _cov( indexlhs+i, jndexlhs+j) = rhs._cov( indexrhs+i, jndexrhs+j) ;
      }

    }
  }
  
}
