#include "BaBar/BaBar.hh"
#include <ErrLogger/ErrLog.hh>
#include "VtxTreeFitter/VtkFitParams.hh"
#include "VtxTreeFitter/VtkParticleBase.hh"
#include "VtxTreeFitter/VtkMergedConstraint.hh"
using std::endl;
using std::ostream;

namespace vtxtreefit
{

  ErrCode
  MergedConstraint::project(const FitParams& fitpar, 
			    Projection& p) const 
  {
    ErrCode status ;
    for(constraintlist::const_iterator it = _list.begin() ;
	it != _list.end() ; ++it) {
      status |= (*it)->project(fitpar,p) ;
      p.incrementOffset((*it)->dim()) ;
    }
 
    return status ;
  }

  void MergedConstraint::print(ostream& os) const
  {
    os << "Merged constraint: " << endl ;
    for(constraintlist::const_iterator it = _list.begin() ;
	it != _list.end() ; ++it) {
      os << "          " ;
      (*it)->print(os) ;
    }
  }

}

