#include "BaBar/BaBar.hh"
#include "VtxTreeFitter/VtkFitParams.hh"
#include "VtxTreeFitter/VtkInteractionPoint.hh"
#include <Beta/BtaCandidate.hh>
#include <ErrLogger/ErrLog.hh>
using std::cout;
using std::endl;

namespace vtxtreefit
{

  extern int vtxverbose ;

  InteractionPoint::InteractionPoint(const BtaCandidate* bc, bool forceFitAll)
    : InternalParticle(bc,0,forceFitAll),
      _constrainXY(false), _constrainXYZ(false), _ipPos(3,0), _ipCov(3,1), _ipCovInv(3,1) 
  {
    //assert(bc->pdtEntry()->lundId()%1000 == PdtLund::Upsilon ) ;
    initBeamSpot(bc) ;
  }

  InteractionPoint::InteractionPoint(const BtaCandidate* bc, bool forceFitAll, bool addupsilon) 
    : InternalParticle(0,0,forceFitAll),
      _constrainXY(false), _constrainXYZ(false), _ipPos(3,0), _ipCov(3,1), _ipCovInv(3,1) 
  {
    addDaughter(bc,forceFitAll) ;
    initBeamSpot(bc) ;
  }

  ErrCode
  InteractionPoint::initBeamSpot(const BtaCandidate* bc)
  {
    ErrCode status ;
    // find the constraint
    const BtaConstraint* btaconstraint = bc->constraint(BtaConstraint::Beam) ;
    bool success = false ;
    if( btaconstraint ) {
      success = 
	btaconstraint->getParmValue("bx",_ipPos(1)) &&
	btaconstraint->getParmValue("by",_ipPos(2)) &&
	btaconstraint->getParmValue("bz",_ipPos(3)) &&
	btaconstraint->getParmValue("sxx",_ipCov.fast(1,1)) &&
	btaconstraint->getParmValue("syy",_ipCov.fast(2,2)) &&
	btaconstraint->getParmValue("szz",_ipCov.fast(3,3)) &&
	btaconstraint->getParmValue("sxy",_ipCov.fast(2,1)) &&
	btaconstraint->getParmValue("sxz",_ipCov.fast(3,1)) &&
	btaconstraint->getParmValue("syz",_ipCov.fast(3,2)) ;
      _constrainXY = _constrainXYZ = success ;
      if(success) {
	int ierr ;
	_ipCovInv = _ipCov.inverse(ierr) ;
      }
    }
    if(!success) {
      ErrMsg(error) << "failed to get beam spot data. constraint will not be applied."
		    << endmsg ;
    } else {
      if(vtxverbose>=2)
	cout << "VtkInteractionPoint: initial beam spot = (" 
	     <<_ipPos(1) << "," << _ipPos(2) << "," << _ipPos(3) << ")" << endl ;
    }
    return status ;
  }

  InteractionPoint::~InteractionPoint() {}

  ErrCode 
  InteractionPoint::initPar1(FitParams* fitparams)
  {
    ErrCode status ;
    int posindex = posIndex() ;
    for(int row=1; row<=3; ++row)
      fitparams->par()(posindex+row) = _ipPos(row) ;

    for(daucontainer::const_iterator it = daughters().begin() ;
	it != daughters().end() ; ++it) {
      status |= (*it)->initPar1(fitparams) ;
      status |= (*it)->initPar2(fitparams) ;
    }

    return status ;
  }

  ErrCode 
  InteractionPoint::initPar2(FitParams* fitparams)
  {
    // nothing left to do: actually, should never be called
    assert(0) ;
    return ErrCode::success ;
  }

  ErrCode
  InteractionPoint::initCov(FitParams* fitpar) const 
  {
    ErrCode status ;
    int posindex = posIndex() ;
    for(int row=1; row<=3; ++row)
      fitpar->cov().fast(posindex+row,posindex+row) 
	= 1000*_ipCov.fast(row,row) ;
    
    for(daucontainer::const_iterator it = daughters().begin() ;
	it != daughters().end() ; ++it)
      status |= (*it)->initCov(fitpar) ;

    return status ;
  }

  ErrCode
  InteractionPoint::projectIPConstraint(const FitParams& fitparams, 
					Projection& p) const
  { 
    int posindex = posIndex() ;
    int maxrow = _constrainXYZ ? 3 : (_constrainXY ? 2 : 0 ) ;
    for(int row=1; row<=maxrow; ++row) {
      p.r(row) = fitparams.par()(posindex+row) - _ipPos(row) ;
      p.H(row,posindex+row) = 1 ;
      for(int col=1; col<=row; ++col)
	p.Vfast(row,col) = _ipCov.fast(row,col) ;
    }
    return ErrCode::success ;
  }
  
  ErrCode 
  InteractionPoint::projectConstraint(Constraint::Type type, 
				      const FitParams& fitparams, 
				      Projection& p) const 
  {
    ErrCode status ;
    switch(type) {
    case Constraint::beamspot:
      status |= projectIPConstraint(fitparams,p) ;
      break ;
    default:
      status |= InternalParticle::projectConstraint(type,fitparams,p) ;
    }
    return status ;
  }
  
  double 
  InteractionPoint::chiSquare(const FitParams* fitparams) const
  {
    // calculate the chi2
    int posindex = posIndex() ;
    HepVector residual = _ipPos - fitparams->par().sub(posindex+1,posindex+3) ;
    double chisq = _ipCovInv.similarity(residual) ;

    // add the daughters
    chisq += InternalParticle::chiSquare(fitparams) ;

    return chisq ;
  }

  void InteractionPoint::addToConstraintList(constraintlist& alist, int depth) const
  {
    // first the beamspot
    int dim = _constrainXYZ ? 3 : (_constrainXY ? 2 : 0 ) ;
    alist.push_back(Constraint(this,Constraint::beamspot,depth,dim)) ;

    // then the base class
    InternalParticle::addToConstraintList(alist,depth) ;
  }

}
