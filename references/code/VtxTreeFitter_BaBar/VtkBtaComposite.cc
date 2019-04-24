#include "BaBar/BaBar.hh"
#include <Beta/BtaCandidate.hh>
#include "VtxTreeFitter/VtkBtaComposite.hh"
#include "VtxTreeFitter/VtkFitParams.hh"
using std::cout;
using std::endl;

namespace vtxtreefit
{

  extern int vtxverbose ;
  BtaComposite::BtaComposite(const BtaCandidate* bc, const ParticleBase* mother)
    : ParticleBase(bc,mother),_m(),_matrixV(),_hasEnergy(true) 
  { 
    bool massconstraint = bc && bc->constraint(BtaConstraint::Mass) ;
    _hasEnergy = !massconstraint ;
    updCache() ;
  }

  void BtaComposite::updCache()
  {
    // cache par7 (x,y,z,px,py,pz,E) cov7

    // the simplest is
    const BtaFitParams btafitparams = bc()->fitParams() ;
    HepPoint pos = btafitparams.pos() ;
    HepLorentzVector mom = btafitparams.p4() ;
    _m = HepVector(dimM()) ;
    _m(1) = pos.x() ;
    _m(2) = pos.y() ;
    _m(3) = pos.z() ;
    _m(4) = mom.x() ;
    _m(5) = mom.y() ;
    _m(6) = mom.z() ;
    _m(7) = mom.t() ;
    _matrixV = btafitparams.cov7().sub(1,dimM()) ; // so either 7 or 6, depending on mass constraint
    
    if(vtxverbose>=4) {
      cout << "cov matrix of external candidate: " << name().c_str() 
	   << " " << dimM() << " " << _matrixV << endl ;
    }
  }

  BtaComposite::~BtaComposite() {}

  ErrCode
  BtaComposite::initPar1(FitParams* fitparams)
  {
    int posindex = posIndex() ;
    int momindex = momIndex() ;

    //quick map for parameters
    int indexmap[7]  ;
    for(int i=0; i<3; ++i) indexmap[i]   = posindex+i ; 
    for(int i=0; i<4; ++i) indexmap[i+3] = momindex+i ;
    // copy the 'measurement'
    for(int row=1; row<=dimM(); ++row) 
      fitparams->par()(indexmap[row-1]+1) = _m(row) ;
    
    return ErrCode::success ;
  }

  ErrCode
  BtaComposite::initPar2(FitParams* fitparams)
  {
    // call default lifetime initialization
    return initTau(fitparams) ;
  }

  ErrCode
  BtaComposite::projectBtaComposite(const FitParams& fitparams, 
				    Projection& p) const
  {
    int posindex = posIndex() ;
    int momindex = momIndex() ;
    
    // quick map for parameters
    int indexmap[7]  ;
    for(int i=0; i<3; ++i) indexmap[i]   = posindex+i ; 
    for(int i=0; i<4; ++i) indexmap[i+3] = momindex+i ;
    for(int row=1; row<=dimM(); ++row) {
      p.r(row)                   = fitparams.par()(indexmap[row-1]+1) - _m(row) ;
      p.H(row,indexmap[row-1]+1) = 1 ;
      for(int col=1; col<=row; ++col)
	p.Vfast(row,col) = _matrixV.fast(row,col) ;
    }
 
    return ErrCode::success ;
  }
  
  ErrCode 
  BtaComposite::projectConstraint(Constraint::Type type, 
				  const FitParams& fitparams, 
				  Projection& p) const 
  {
    ErrCode status ;
    switch(type) {
    case Constraint::btacomposite:
      status |= projectBtaComposite(fitparams,p) ;
      break ;
    case Constraint::geometric:
      status |= projectGeoConstraint(fitparams,p) ;
      break ;
    default:
      status |= ParticleBase::projectConstraint(type,fitparams,p) ;
    }
    return status ;
  }

  double 
  BtaComposite::chiSquare(const FitParams* fitparams) const
  {
    Projection p(fitparams->dim(),dimM()) ;
    projectBtaComposite(*fitparams,p) ;
    return p.chiSquare() ;
  }

}
