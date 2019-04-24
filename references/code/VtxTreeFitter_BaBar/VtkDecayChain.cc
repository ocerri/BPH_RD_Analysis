#include "BaBar/BaBar.hh"
#include <algorithm>
#include "VtxTreeFitter/VtkFitParams.hh"
#include "VtxTreeFitter/VtkParticleBase.hh"
#include "VtxTreeFitter/VtkDecayChain.hh"
using std::cout;
using std::endl;
using std::ostream;

namespace vtxtreefit
{
  
  extern int vtxverbose ;
  
  void
  DecayChain::printConstraints(ostream& os) const
  {
    os << "Constraints in decay tree: " << endl ;
    for( ParticleBase::constraintlist::const_iterator 
	   it = _constraintlist.begin() ;
	 it != _constraintlist.end(); ++it) 
      it->print(os) ;
  }

  DecayChain::DecayChain(const BtaCandidate* bc, bool forceFitAll) 
    : _dim(0),_mother(0),_isOwner(true)
  {
    _mother = ParticleBase::createParticle(bc,0,forceFitAll) ;
    _mother->updateIndex(_dim) ;
    _cand   = locate(bc) ;
    initConstraintList() ;
    
    if(vtxverbose>=2) {
      cout << "In DecayChain constructor: _dim = " << _dim << endl ;
      printConstraints() ;
    }
  }
  
  DecayChain::~DecayChain()
  { 
    if(_mother && _isOwner) delete _mother ; 
  }

  void
  DecayChain::initConstraintList()
  {
    _constraintlist.clear() ;
    _mother->addToConstraintList(_constraintlist,0) ;
    // the order of the constraints is a rather delicate thing
    std::sort(_constraintlist.begin(),_constraintlist.end()) ;

    // merge all non-lineair constraints
    _mergedconstraintlist.clear() ;
    mergedconstraint = MergedConstraint() ;
    for( ParticleBase::constraintlist::iterator it =  _constraintlist.begin() ;
	 it != _constraintlist.end(); ++it) {
      if( it->isLineair() ) _mergedconstraintlist.push_back(&(*it)) ;
      else  mergedconstraint.push_back(&(*it)) ;
    }
    
    if( mergedconstraint.dim()>0 )
      _mergedconstraintlist.push_back(&mergedconstraint) ;
  }

  ErrCode
  DecayChain::init(FitParams* par) 
  {
    ErrCode status ;

    // set everything to 0
    par->resetPar() ;
    status |= _mother->initPar1(par) ;

    // let the mother do it
    par->resetCov() ;
    status |= _mother->initCov(par) ;
   
    if(vtxverbose>=2) cout << "status after init: " << status << endl ;
 
    return status ;
  }

  ErrCode
  DecayChain::filter(FitParams* par, bool firstpass) 
  {
    ErrCode status ;
    par->resetCov(1000) ;
    if(firstpass || !par->testCov()) status |= _mother->initCov(par) ;


    if( vtxverbose>=3 || (vtxverbose>=2&&firstpass) ) {
      cout << "DecayChain::filter, after initialization: "  << endl ; 
      _mother->print(par) ;
    }
    
#ifdef THEOLDWAY
    for( ParticleBase::constraintlist::const_iterator it = _constraintlist.begin() ;
     	 it != _constraintlist.end(); ++it) {
      status |= it->filter(par) ;
      if( vtxverbose>=2 && status.failure() ) {
	cout << "status is failure after parsing constraint: " ;
	it->print() ;
      }
    }
#else
    for( std::vector<Constraint*>::const_iterator it = _mergedconstraintlist.begin() ;
	  it != _mergedconstraintlist.end(); ++it) {
      status |= (*it)->filter(par) ;
      if( vtxverbose>=2 && status.failure() ) {
	cout << "status is failure after parsing constraint: " ;
	(*it)->print() ;
      }
    }
#endif
    
    if(vtxverbose>=3) cout << "VtkDecayChain::filter: status = " << status << endl ;
      
    return status ;
  }

  double 
  DecayChain::chiSquare(const FitParams* par) const
  {
    return _mother->chiSquare(par) ;
  }

  const ParticleBase* 
  DecayChain::locate(const BtaCandidate* bc) const {
    const ParticleBase* rc(0) ;
    // return _mother->locate(bc) ;
    ParticleMap::const_iterator it = _particleMap.find(bc) ;
    if( it == _particleMap.end() ) {
      rc = _mother->locate(bc) ;
      // used to add every candidate here, but something goes wrong
      // somewhere. now only cache pointers that we internally reference.
      if(rc && rc->bc())
	const_cast<DecayChain*>(this)->_particleMap[rc->bc()] = rc ;
    } else {
      rc = it->second ;
    }
    return rc ;
  }
  
  int 
  DecayChain::index(const BtaCandidate* bc) const {
    int rc = -1 ;
    const ParticleBase* base = locate(bc) ;
    if(base) rc = base->index() ;
    return rc ;
  }

  int 
  DecayChain::posIndex(const BtaCandidate* bc) const {
    int rc = -1 ;
    const ParticleBase* base = locate(bc) ;
    if(base) rc = base->posIndex() ; 
    return rc ;
  }

  int 
  DecayChain::momIndex(const BtaCandidate* bc) const {
    int rc = -1 ;
    const ParticleBase* base = locate(bc) ;
    if(base) rc = base->momIndex() ; 
    return rc ;
  }

  int 
  DecayChain::tauIndex(const BtaCandidate* bc) const {
    int rc = -1 ;
    const ParticleBase* base = locate(bc) ;
    if(base) rc = base->tauIndex() ;  
    return rc ;
  }
}
