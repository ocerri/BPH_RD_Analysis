#include "BaBar/BaBar.hh"
#include <Beta/BtaCandidate.hh>
#include <PDT/Pdt.hh>

#include "VtxTreeFitter/VtkUpsilon.hh"
#include "VtxTreeFitter/VtkMissingParticle.hh"
#include "VtxTreeFitter/VtkFitParams.hh"
using std::cout;
using std::endl;

namespace vtxtreefit
{

  extern int vtxverbose ;

  Upsilon::Upsilon(const BtaCandidate* bc, bool forceFitAll)
    : InteractionPoint(bc,forceFitAll),
      _constrainBeam(false),_constrainLifetimeSum(false),
      _beamMom(4,0), _beamCov(4,0), _beamCovInv(4,1), _tagBcand(0)
  {
    initBeamEnergy(bc) ;
    _constrainLifetimeSum = bc->constraint(BtaConstraint::Life)!=0 ;
  }

  Upsilon::Upsilon(const BtaCandidate* bc, bool forceFitAll, bool addupsilon) 
    : InteractionPoint(bc,forceFitAll,addupsilon),
      _constrainBeam(false),_constrainLifetimeSum(false),
      _beamMom(4,0), _beamCov(4,0), _beamCovInv(4,1), _tagBcand(0)
  {
    initBeamEnergy(bc) ;
    // create a missing particle for the tag B
    const PdtEntry* pdt = bc->pdtEntry()->conjugate() ;
    HepLorentzVector beammom(_beamMom[0],_beamMom[1],_beamMom[2],_beamMom[3]) ;
    HepLorentzVector missingmom = beammom - bc->p4() ;
    _tagBcand = new BtaCandidate(missingmom,pdt) ;
    _tagBcand->addConstraint(BtaConstraint::Mass) ;
    daughters().push_back(new MissingParticle(_tagBcand,this)) ;
  }
 
  ErrCode Upsilon::initBeamEnergy(const BtaCandidate* bc)
  {
    const BtaConstraint* btaconstraint =  bc->constraint(BtaConstraint::BeamEnergy) ;
    ErrCode status ;
    bool success = false ;
    if( btaconstraint ) {
      HepVector    beamp3[2]    = {HepVector(3),HepVector(3)} ;
      HepSymMatrix beamp3cov[2] = {HepSymMatrix(3),HepSymMatrix(3)} ;
      //std::vector<HepVector>    beamp3(2,HepVector(3)) ;
      //std::vector<HepSymMatrix> beamp3cov(2,HepSymMatrix(3)) ;
      success = 
	btaconstraint->getParmValue("eMinusPx", beamp3[0](1)) &&
	btaconstraint->getParmValue("eMinusPy", beamp3[0](2)) &&
	btaconstraint->getParmValue("eMinusPz", beamp3[0](3)) &&
	btaconstraint->getParmValue("eMinusCovPxPx", beamp3cov[0].fast(1,1)) &&
	btaconstraint->getParmValue("eMinusCovPxPy", beamp3cov[0].fast(2,1)) &&
	btaconstraint->getParmValue("eMinusCovPxPz", beamp3cov[0].fast(3,1)) &&
	btaconstraint->getParmValue("eMinusCovPyPy", beamp3cov[0].fast(2,2)) &&
	btaconstraint->getParmValue("eMinusCovPyPz", beamp3cov[0].fast(3,2)) &&
	btaconstraint->getParmValue("eMinusCovPzPz", beamp3cov[0].fast(3,3)) &&
	btaconstraint->getParmValue("ePlusPx", beamp3[1](1)) &&
	btaconstraint->getParmValue("ePlusPy", beamp3[1](2)) &&
	btaconstraint->getParmValue("ePlusPz", beamp3[1](3)) &&
	btaconstraint->getParmValue("ePlusCovPxPx", beamp3cov[1].fast(1,1)) &&
	btaconstraint->getParmValue("ePlusCovPxPy", beamp3cov[1].fast(2,1)) &&
	btaconstraint->getParmValue("ePlusCovPxPz", beamp3cov[1].fast(3,1)) &&
	btaconstraint->getParmValue("ePlusCovPyPy", beamp3cov[1].fast(2,2)) &&
	btaconstraint->getParmValue("ePlusCovPyPz", beamp3cov[1].fast(3,2)) &&
	btaconstraint->getParmValue("ePlusCovPzPz", beamp3cov[1].fast(3,3)) ;
      
      _constrainBeam = success ;
      if(success) {
	// make 4 vectors of these three vectors and add them up
	const double emass  = 0.0005 ;
	for(int i=0; i<2; ++i) {
	  // add to the momentum
	  double energy2 = emass*emass ;
	  for(int row=1; row<=3; ++row) {
	    double px = beamp3[i](row) ;
	    _beamMom(row) += px ;
	    energy2 += px*px ;
	  }
	  double energy = sqrt(energy2) ;
	  _beamMom(4) += energy ;
	  
	  // calculate the jacobian
	  HepMatrix jacobian(4,3,0) ;
	  jacobian(1,1) = 1 ;
	  jacobian(2,2) = 1 ;
	  jacobian(3,3) = 1  ;
	  for(int col=1; col<=3; ++col)
	    jacobian(4,col) = beamp3[i](col)/energy ;

	  _beamCov += beamp3cov[i].similarity(jacobian) ;
	}
	
	// calculate the weight matrix
	int ierr ;
	_beamCovInv = _beamCov.inverse(ierr) ;
      }
    }
    
    if(!success) {
      cout << "WARNING: failed to get beam energy data. constraint will not be applied."
	   << endl ;
    } else {
      if(vtxverbose>=2)
	cout << "VtkUpsilon: initial beam energy = (" 
	     <<_beamMom(1) << "," << _beamMom(2) << "," << _beamMom(3) << "," << _beamMom(4) << ")" << endl ;
    }
    return status ;
  }

  Upsilon::~Upsilon() 
  {
    delete _tagBcand ;
  }

  ErrCode
  Upsilon::initPar1(FitParams* fitpar)
  {
    int momindex = momIndex() ;
    for(int row=1; row<=4; ++row)
      fitpar->par()(momindex+row) = _beamMom(row) ;

    return InteractionPoint::initPar1(fitpar) ;
  }

  ErrCode 
  Upsilon::initCov(FitParams* fitpar) const 
  {
    int momindex = momIndex() ;
    for(int row=1; row<=4; ++row)
      fitpar->cov().fast(momindex+row,momindex+row) 
	= 1000*_beamCov.fast(row,row) ;
    return InteractionPoint::initCov(fitpar) ;
  } 

  ErrCode
  Upsilon::projectBeamEnergyConstraint(const FitParams& fitparams, 
				       Projection& p) const
  {
    int momindex = momIndex() ;
    for(int row=1; row<=4; ++row) {
      p.r(row)   = fitparams.par()(momindex+row) - _beamMom(row) ;
      p.H(row,momindex+row) = 1 ;
      for(int col=1; col<=row; ++col)
	p.Vfast(row,col) = _beamCov.fast(row,col) ;
    }
    return ErrCode::success ;
  }

  ErrCode
  Upsilon::projectLifetimeSumConstraint(const FitParams& fitparams, 
					Projection& p) const
  {
    double pdttausum(0),pdttau2sum(0) ;
    double tausum(0) ;
    for(daucontainer::const_iterator it = daughters().begin() ;
	it != daughters().end() ; ++it) {
      int tauindex = (*it)->tauIndex() ;
      assert(tauindex>=0) ;
      p.H(1,tauindex+1) = 1 ;
      tausum    += fitparams.par()(tauindex+1) ;
      double pdttaudau = (*it)->pdtTau() ;
      pdttausum  += pdttaudau ;
      pdttau2sum += pdttaudau * pdttaudau ;
    }
    p.r(1)       = tausum - pdttausum ;
    p.Vfast(1,1) = pdttau2sum ;
    return ErrCode::success ;
  }

  ErrCode 
  Upsilon::projectConstraint(Constraint::Type type, 
			     const FitParams& fitparams, 
			     Projection& p) const 
  {
    ErrCode status ;
    switch(type) {
    case Constraint::beamenergy:
      status |= projectBeamEnergyConstraint(fitparams,p) ;
      break ;
    case Constraint::lifetime:
      status |= projectLifetimeSumConstraint(fitparams,p) ;
      break ;
    default:
      status |= InteractionPoint::projectConstraint(type,fitparams,p) ;
    }
    return status ;
  }

  double 
  Upsilon::chiSquare(const FitParams* fitparams) const
  {
    // calculate the chi2 to the beam momentum
    int momindex = momIndex() ;
    HepVector residual = _beamMom - fitparams->par().sub(momindex+1,momindex+4) ;
    double chisq = _beamCovInv.similarity(residual) ;
    
    // add the rest
    chisq += InteractionPoint::chiSquare(fitparams) ;
    return chisq ;
  }
  
  void Upsilon::addToConstraintList(constraintlist& alist, int depth) const
  {
    // first the beamenergy
    alist.push_back(Constraint(this,Constraint::beamenergy,depth,4)) ;

    // now the lifetime
    if(_constrainLifetimeSum) 
      alist.push_back(Constraint(this,Constraint::lifetime,depth,1)) ;
    
    // then the base class
    InteractionPoint::addToConstraintList(alist,depth) ;
  } 

//   void Upsilon::replaceTagB(BtaCandidate* newtagB)
//   {
//     cout << "DO NOT CALL" << endl ;
//     assert(0) ;
//     // remove the tagging B
//     daughters().pop_back() ;
//     delete _tagBcand ;
   
//     // copy the new one
//     _tagBcand = newtagB ;

//     // add a mass constraint to the tagging B
//     _tagBcand->addConstraint(BtaConstraint::Mass) ;

//     // create a new particle
//     InternalParticle* pb = new InternalParticle(_tagBcand,this,true) ;

//     // add a particle for the missing momentum
//     pb->daughters().push_back( new MissingParticle(0,pb)) ;
    
    
//     // add to the daughters
//     daughters().push_back( pb ) ;
//   }

}
