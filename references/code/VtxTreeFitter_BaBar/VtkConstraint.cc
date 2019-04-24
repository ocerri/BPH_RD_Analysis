#include "BaBar/BaBar.hh"
#include <iomanip>
#include <ErrLogger/ErrLog.hh>
#include "VtxTreeFitter/VtkFitParams.hh"
#include "VtxTreeFitter/VtkParticleBase.hh"
#include "VtxTreeFitter/VtkConstraint.hh"
#include "VtxTreeFitter/VtkKalmanCalculator.hh"
using std::cout;
using std::endl;
using std::ostream;
using std::setprecision;
using std::setw;

namespace vtxtreefit
{

  extern int vtxverbose ;

  bool Constraint::operator<(const Constraint& rhs) const
  { 
    // the simple way
    return _type < rhs._type ||
      (_type == rhs._type && _depth < rhs._depth) ; 

    // this is probably the second most complicated routine: how do we
    // order the constraints. there is one very special case:
    // Ks->pipi0 requires the pi0 mass constraints at the very
    // end. otherwise, it just doesn't work. in all other cases, we
    // prefer to fit 'down' the tree'. the 'external' constraints must
    // be filtered first, but soft pions must be fitted after the
    // geometric constraints of the D. You see, this is horrible.

    // if either of the two is external, or either of the two is a
    // mass constraint, we order by _type_
    if( (_type <= Constraint::btacomposite ||
	 rhs._type <= Constraint::btacomposite ) ||
	(_type >= Constraint::mass ||
	 rhs._type >= Constraint::mass ) ) {
      return _type < rhs._type ||
	(_type == rhs._type && _depth < rhs._depth) ;
    } 
    // if not, we order by depth
    return _depth < rhs._depth  ||
      (_depth == rhs._depth && _type < rhs._type ) ;
    
  }
  
  ErrCode
  Constraint::project(const FitParams& fitpar, Projection& p) const
  {
    // this one will be overruled by the MergedConstraint
    return _node->projectConstraint(_type,fitpar,p) ;
  }

  ErrCode
  Constraint::filter(FitParams* fitpar) const
  {
    ErrCode status ;
    if(_type<=Constraint::unknown || _type>=Constraint::ntypes) {
      ErrMsg(error) << "VtkConstraint: unknown constraint: " << _type << endmsg ;
      status |= ErrCode::badsetup ;
    } else if (_type!=merged && !_node) {
      ErrMsg(error) << "VtkConstraint: filter constraint without a node" << endmsg ;
      status |= ErrCode::badsetup ;
    } else {
      if(vtxverbose>=3) { cout << "filtering "  ; print() ;}
      // save the unfiltered ('predicted') parameters. we need to
      // store them if we want to iterate constraints.
      const HepVector* pred(0) ;
      if(_maxNIter>1) pred = new HepVector(fitpar->par()) ;

      Projection p(fitpar->dim(),_dim) ;
      KalmanCalculator kalman ;
      double chisq(0) ;
      int iter(0) ;
      bool finished(false) ;
      while (!finished && !status.failure()) {
	p.reset() ;
	status |= project(*fitpar,p) ;
	if(!status.failure()) {
	  status |= kalman.init( p.r(), p.H(), fitpar, &p.V() ) ;
	  if( !status.failure()) {
	    if(iter==0 || !pred) {
	      kalman.updatePar( fitpar ) ;
	    } else {
	      kalman.updatePar( *pred, fitpar ) ;
	    }
	    const double dchisqconverged = 0.001 ;
	    double newchisq = kalman.chisq() ;
	    double dchisq = newchisq - chisq ;
	    bool diverging = iter > 0 && dchisq>0  ;
	    bool converged = fabs(dchisq) < dchisqconverged ;
	    finished  = ++iter >= _maxNIter || diverging || converged ;
	    
	    if(vtxverbose>=3) { 
	      cout << "chi2,niter: " 
		   << iter << " "<< setprecision(7)
		   << setw(12) << chisq << " " 
		   << setw(12)<< newchisq << " "
		   << setw(12)<< dchisq << " " 
		   << diverging << " " 
		   << converged << " "
		   << status << endl ; 
	    }
	    chisq = newchisq ;
	  }
	}
      }
      if(!status.failure()) {
	kalman.updateCov( fitpar ) ;
	if(_nHidden>0) fitpar->addChiSquare(0,-_nHidden) ;
      } 
      if(pred) delete pred ;
      if(vtxverbose>=4 &&_node&&
	 _node->mother() ) { _node->mother()->print(fitpar) ; }
    }
    return status ;
  }

  void Constraint::print(ostream& os) const
  {
    os << _node->index() << " "  
       <<_node->name().c_str() << " " 
       << name().c_str() << " "
       << _type << " " << _depth << endl ;
  }

  std::string Constraint::name() const
  {
    std::string rc = "unknown constraint!" ;
    switch(_type) 
      {
      case beamspot:     rc = "beamspot" ; break ;
      case beamenergy:   rc = "beamenergy" ; break ;
      case btacomposite: rc = "btacomposite" ; break ;
      case btaresonance: rc = "btaresonance" ; break ;
      case track:        rc = "track" ; break ;
      case photon:       rc = "photon" ; break ;
      case kinematic:    rc = "kinematic" ; break ;
      case geometric:    rc = "geometric" ; break ;
      case mass:         rc = "mass" ; break ;
      case massEnergy:   rc = "massEnergy" ; break ;
      case lifetime:     rc = "lifetime" ; break ;
      case merged:       rc = "merged" ; break ;
      case conversion:   rc = "conversion" ; break ;
      case ntypes:
      case unknown: 
	break ;
      }
    return rc ;
  }
}
