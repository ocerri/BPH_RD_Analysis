#ifndef __VTK_FITPARAMS_HH__
#define __VTK_FITPARAMS_HH__

#include <vector> 
#include <CLHEP/Matrix/Vector.h>
#include <CLHEP/Matrix/SymMatrix.h>
#include <CLHEP/Matrix/DiagMatrix.h>

namespace vtxtreefit
{
  class ParticleBase ;

  class FitParams
  {
  public:
    // Class that contains the parameters and covariance for the
    // vertex fit.
    FitParams(int dim) ;
    ~FitParams() ;
    
    HepSymMatrix& cov() { return _cov ; }
    HepVector& par() { return _par ; }
    HepDouble& par(int row) { return _par(row) ; }
    
    HepSymMatrix cov(const std::vector<int>& indexVec) const ;
    HepVector par(const std::vector<int>& indexVec) const ;
    
    const HepSymMatrix& cov() const { return _cov ; }
    const HepVector& par() const { return _par ; }
    const HepDouble& par(int row) const { return _par(row) ; }

    HepDiagMatrix& scale() { return _scale ; }

    int& nConstraintsVec(int row) { return _nConstraintsVec[row-1] ; }

    //int dim() const { return _par.num_row() ; }
    int dim() const { return _dim ; }
    double chiSquare() const { return _chiSquare ; }

    int nConstraints() const { return _nConstraints ; }
    int nDof() const { return nConstraints() - dim() ; }
    HepDouble err(int row) const { return sqrt(_cov(row,row)) ; }

    void resize(int newdim) ;
    void resetPar() ;
    void resetCov(double scale=100) ;
    void print() const ;
    bool testCov() const ;
    void addChiSquare(double chisq, int nconstraints) {
      _chiSquare += chisq ; _nConstraints += nconstraints ; }

    typedef std::vector< std::pair<const ParticleBase*,int> > indexmap ;
    void copy(const FitParams& rhs, const indexmap& anindexmap) ;
  protected:
    FitParams() {}
  private:
    int _dim ;
    HepVector    _par ;
    HepSymMatrix _cov ;
    HepDiagMatrix _scale ;
    double _chiSquare ;
    int _nConstraints ;
    std::vector<int> _nConstraintsVec ; // vector with number of constraints per parameter
  } ;
} 

#endif
