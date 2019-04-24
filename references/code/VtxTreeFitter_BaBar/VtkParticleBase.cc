#include "BaBar/BaBar.hh"
#include <iomanip>
#include <float.h>

#include <AbsEnv/AbsEnv.hh>
#include <BaBar/Constants.hh>
#include <BField/BField.hh>
#include <BbrGeom/BbrDoubleErr.hh>
#include <Beta/BtaCandidate.hh>
#include <Beta/BtaMcTruth.hh>
#include <BetaCoreTools/BtaMcAssoc.hh>
#include <ErrLogger/ErrLog.hh>
#include <G3Data/GTrack.hh>
#include <G3Data/GVertex.hh>
#include <TrkEnv/TrkEnv.hh>
#include "VtxTreeFitter/VtkParticleBase.hh"
#include "VtxTreeFitter/VtkInternalParticle.hh"
#include "VtxTreeFitter/VtkBtaComposite.hh"
#include "VtxTreeFitter/VtkBtaResonance.hh"
#include "VtxTreeFitter/VtkRecoTrack.hh"
#include "VtxTreeFitter/VtkRecoPhoton.hh"
#include "VtxTreeFitter/VtkResonance.hh"
#include "VtxTreeFitter/VtkInteractionPoint.hh"
#include "VtxTreeFitter/VtkUpsilon.hh"
#include "VtxTreeFitter/VtkMissingParticle.hh"
#include "VtxTreeFitter/VtkFitParams.hh"
using std::cout;
using std::endl;
using std::flush;
using std::setprecision;
using std::setw;

namespace vtxtreefit
{

  //template<class T>
  //inline T sqr(T x) { return x*x ; }

  int vtxverbose=0 ;

  ParticleBase::ParticleBase(const BtaCandidate* bc, const ParticleBase* mother)
    : _bc(bc),_mother(mother),_index(0),_pdtMass(0),_pdtLifeTime(0),_charge(0) 
  {
    if(bc) {
      const PdtEntry* pdt(0) ;
      if( (pdt = bc->pdtEntry()) ) {
	_pdtMass     = pdt->mass() ;
	_pdtLifeTime = pdtLifeTime(bc) ;
	double fltcharge =  pdt->charge() ;
	// round to nearest integer
	_charge = fltcharge < 0 ? int(fltcharge-0.5) : int(fltcharge+0.5) ;
      } else {
	_charge = bc->charge()>0 ? 1 : (bc->charge()<0 ? -1 : 0) ;
      }
    }
  }

  ParticleBase::~ParticleBase() {} ;

  void ParticleBase::updateIndex(int& offset)
  {
    _index = offset ;
    offset += dim() ;
  }

  ParticleBase*
  ParticleBase::createParticle(const BtaCandidate* bc, const ParticleBase* mother, 
			       bool forceFitAll)
  {
    // This routine interpretes a beta candidate as one of the
    // 'Particles' used by the fitter.

    if(vtxverbose>=2)
      cout << "ParticleBase::createParticle: " << forceFitAll << endl ;
    ParticleBase* rc=0 ;
    bool isupsilon    = bc->pdtEntry() && bc->pdtEntry()->lundId()%1000 == PdtLund::Upsilon ;
    bool bsconstraint = bc->constraint(BtaConstraint::Beam) !=0 ;
    bool beconstraint = bc->constraint(BtaConstraint::BeamEnergy) != 0;
    // We refit invalid fits, kinematic fits and composites with beamspot
    // constraint if not at head of tree.
    bool validfit     = !forceFitAll && mother && bc->decayVtx() && 
      (bc->decayVtx()->status() == BtaAbsVertex::Success ) &&
      (bc->decayVtx()->type() == BtaAbsVertex::Geometric ) &&
      !(bsconstraint||beconstraint) ;
    if(validfit && (bc->fitParams().cov7()(7,1)==0 && 
		    bc->fitParams().cov7()(6,1)!=0) ) {
      // this is an awfull hack to detect whether the covariance
      // matrix contains a x-E row. if not, then the matrix is
      // unusable. GeoKin does not fill the xE row.
      static int printit=10 ;
      if(--printit>0)
	ErrMsg(warning) << "ParticleBase: this btacandidate (" << bc->pdtEntry()->name() 
			<< ") has a bad cov matrix. I'll refit it. " << endmsg ;
      validfit = false ;
    }

    // leave this one for now
    if(bc->pdtEntry() && bc->pdtEntry()->lundId() == PdtLund::pi0 &&
       validfit) {
      static int printit=10 ;
      if(--printit>=0)
	ErrMsg(warning) << "VtkParticleBase::createParticle: found pi0 with valid fit." << endl
			<< "This is likely a configuration error." << endmsg ;
      validfit = false ;
    }
    
    if(!mother) { // 'head of tree' particles
      if( beconstraint )
	rc = isupsilon ? new Upsilon(bc,forceFitAll) : new Upsilon(bc,forceFitAll,true) ;
      else if ( bsconstraint ) 
	rc = isupsilon ? new InteractionPoint(bc,forceFitAll) : new InteractionPoint(bc,forceFitAll,true) ;
      else if( bc->isComposite() ) 
	rc = new InternalParticle(bc,0,forceFitAll) ; // still need proper head-of-tree class
      else {
	ErrMsg(error) << "VtkParticleBase::createParticle: You are fitting a decay tree that exists of "
		      << "a single, non-composite particle and which does not have a beamconstraint."
		      << "I do not understand what you want me to do with this." << endmsg ;
	rc = new InternalParticle(bc,0,forceFitAll) ; // still need proper head-of-tree class
      }
    } else if( !(bc->isComposite()) ) { // external particles
      if( bc->recoTrk() )
	rc = new RecoTrack(bc,mother) ;  // reconstructed track
      else if( bc->recoCalo() ) 
	rc = new RecoPhoton(bc,mother) ; // reconstructed photon 
      else if(!validfit || bc->constraint(BtaConstraint::MissingMass) )
	rc = new MissingParticle(bc,mother) ; // missing particle
      else if( isAResonance(bc) )
	rc = new BtaResonance(bc,mother) ;
      else
	rc = new BtaComposite(bc,mother) ;
    } else { // 'internal' particles
      if( validfit /*|| isconversion*/ ) {  // fitted composites
	if( isAResonance(bc) )
	  rc = new BtaResonance(bc,mother) ;
	else
	  rc = new BtaComposite(bc,mother) ;
      } else {         // unfited composites
	if( isAResonance(bc) ) 
	  rc = new Resonance(bc,mother,forceFitAll) ;
	else
	  rc = new InternalParticle(bc,mother,forceFitAll) ;
      }
    } 
    
    if(vtxverbose>=2)
      cout << "ParticleBase::createParticle returns " << rc->type() 
	   << " " << rc->index() << endl ;
    return rc ;
  }

  std::string 
  ParticleBase::name() const 
  {
    return (bc() && bc()->pdtEntry()) ? bc()->pdtEntry()->name() : "Unknown" ;
  }

  double 
  ParticleBase::pdtLifeTime(const BtaCandidate* bc)
  {
    double lifetime = 0;
    if( bc->pdtEntry() && bc->pdtEntry()->lifetime()!=FLT_MAX)
      lifetime =  bc->pdtEntry()->lifetime() ;
    return lifetime ;
  }

  bool 
  ParticleBase::isAResonance(const BtaCandidate* bc) {
    bool rc = false ;
    const PdtEntry* pdt = bc->pdtEntry()  ;
    if( pdt && bc->isComposite() ) {
      switch(pdt->lundId()) {
      case PdtLund::gamma:  // conversions are not treated as a resonance
	rc = false; 
	break ;
      case PdtLund::e_plus: // bremstrahlung is treated as a resonance
      case PdtLund::e_minus:
	rc = true ;
	break ;
      default: // this should take care of the pi0
	rc = bc->isAResonance() || (bc->pdtEntry() && pdtLifeTime(bc)<1.e-8) ;
      }
    }
    return rc ;
  }

  ErrCode
  ParticleBase::initCov(FitParams* fitparams) const
  {
    ErrCode status ;

    // position
    int posindex = posIndex() ;
    if( posindex>=0 ) {
      //double decaylength = pdtLifeTime() ;
      const double sigpos = 20 ; // cm
      // that's how good the initalization should be
      for(int row=posindex+1; row<=posindex+3 ; ++row)
	fitparams->cov().fast(row,row) = sigpos*sigpos ;
    }

    // momentum
    int momindex = momIndex() ;   
    if(momindex>=0) {
      const double sigmom = 1 ; // GeV
      int maxrow = hasEnergy() ? 4 : 3 ;
      for(int row=momindex+1; row<=momindex+maxrow; ++row)
	fitparams->cov().fast(row,row) = sigmom*sigmom ;
    }
    
    // lifetime
    int tauindex = tauIndex() ;
    if(tauindex>=0) {
      double tau = pdtTau() ;
      double sigtau = tau>0 ? 20*tau : 999 ;
      const double maxdecaylength = 20; // [cm] (okay for Ks->pi0pi0)
      double bcP = bc()->p();
      if ( bcP > 0.0 )
	sigtau = std::min( maxdecaylength/bcP, sigtau ) ;
      fitparams->cov().fast(tauindex+1,tauindex+1) = sigtau*sigtau ;
    }

    return status ;
  }

  std::string ParticleBase::parname(int thisindex) const
  {
    std::string rc = name() ;
    switch(thisindex) {
    case 0: rc += "_x  " ; break ;
    case 1: rc += "_y  " ; break ;
    case 2: rc += "_z  " ; break ;
    case 3: rc += "_tau" ; break ;
    case 4: rc += "_px " ; break ;
    case 5: rc += "_py " ; break ;
    case 6: rc += "_pz " ; break ;
    case 7: rc += "_E  " ; break ;
    default: ;
    }
    return rc ;
  }

  void 
  ParticleBase::print(const FitParams* fitpar) const
  {
    cout << setw(5) << "[" << type() << "]" << setw(15) << flush << name().c_str() 
	 << " val" << setw(15) << "err" << endl ;
    cout << setprecision(5) ;
    for(int i=0; i<dim(); ++i) {
      int theindex = index()+i ;
      cout << setw(2) << theindex << " "
	   << setw(20) << parname(i).c_str() 
	   << setw(15) << fitpar->par()(theindex+1)
	   << setw(15) << sqrt(fitpar->cov()(theindex+1,theindex+1)) 
	   << setw(15) << fitpar->cov()(theindex+1,theindex+1) <<endl ;
    }
    if( hasEnergy() ) {
      int momindex = momIndex() ;
      double E  = fitpar->par()(momindex+4) ;
      double px = fitpar->par()(momindex+1) ;
      double py = fitpar->par()(momindex+2) ;
      double pz = fitpar->par()(momindex+3) ;
      double mass2 = E*E-px*px-py*py-pz*pz ;
      double mass = mass2>0 ? sqrt(mass2) : -sqrt(-mass2) ;
      
      HepSymMatrix cov = fitpar->cov().sub(momindex+1,momindex+4) ;
      HepVector G(4,0) ;
      G(1) = -px/mass ;
      G(2) = -py/mass ;
      G(3) = -pz/mass ;
      G(4) =   E/mass ;
      double massvar = cov.similarity(G) ;
      cout << setw(2) << setw(20) << "mass: "
	   << setw(15) << mass
	   << setw(15) << sqrt(massvar) << endl ;
    }
  }

  const 
  ParticleBase* ParticleBase::locate(const BtaCandidate* abc) const
  {
    const ParticleBase* rc = 0;
    if( bc() && ( bc()==abc || bc()->isCloneOf(*abc,true) ) ) rc = this ;
    return rc ;
  }
   
  void ParticleBase::retrieveIndexMap(indexmap& anindexmap) const 
  {
    anindexmap.push_back(std::pair<const ParticleBase*,int>(this,index())) ;
  }

  ErrCode ParticleBase::projectGeoConstraint(const FitParams& fitparams,
					     Projection& p) const
  {
    int posindexmother = mother()->posIndex() ;
    int posindex = posIndex();
    int tauindex = tauIndex() ;
    int momindex = momIndex() ;
 
    double tau =  fitparams.par()(tauindex+1) ;

    // lineair approximation is fine
    for(int row=1; row<=3; ++row) {
      double posxmother = fitparams.par()(posindexmother+row) ;
      double posx       = fitparams.par()(posindex+row) ;
      double momx       = fitparams.par()(momindex+row) ;
      p.r(row) = posxmother - (posx - tau*momx) ;
      p.H(row,posindexmother+row) = 1 ;
      p.H(row,posindex+row)       = -1 ;
      p.H(row,momindex+row)       = tau ;
      p.H(row,tauindex+1)         = momx ;
    }
    
    if( charge()!=0 ) {
      double lambda = bFieldOverC() * charge() ; 
      double px0 = fitparams.par()(momindex+1) ;
      double py0 = fitparams.par()(momindex+2) ;
      double pt0 = sqrt(px0*px0+py0*py0) ;
      const double posprecision = 1e-4 ; // 1mu
      if( fabs(pt0*lambda*tau*tau) > posprecision ) {
	// use the helix, but as if it were a 'correction'
	double sinlt = sin(lambda*tau) ;
	double coslt = cos(lambda*tau) ;
	double px = px0*coslt - py0*sinlt ;
	double py = py0*coslt + px0*sinlt ;
	
	p.r(1) += -tau*px0 + (py-py0)/lambda ;
	p.r(2) += -tau*py0 - (px-px0)/lambda ;

	p.H(1,tauindex+1) += -px0 + px ;
	p.H(1,momindex+1) += -tau + sinlt/lambda ;
	p.H(1,momindex+2) +=        (coslt-1)/lambda ;
	p.H(2,tauindex+1) += -py0 + py ;
	p.H(2,momindex+1) +=      - (coslt-1)/lambda ;
	p.H(2,momindex+2) += -tau + sinlt/lambda ;

	if(vtxverbose>=2)
	  cout << "Using helix for position of particle: " << name().c_str() << " "
	       << lambda << " " << lambda*tau  
	       << "  delta-x,y: " << -tau*px0 + (py-py0)/lambda << "  "
	       << -tau*py0 - (px-px0)/lambda << endl ;
      }
    }
    return ErrCode::success ;
  }

  ErrCode ParticleBase::projectMassConstraint(const FitParams& fitparams,
					      Projection& p) const
  {
    double mass = pdtMass() ;
    double mass2 = mass*mass ;
    int momindex = momIndex() ;

    // initialize the value
    double px = fitparams.par()(momindex+1) ;
    double py = fitparams.par()(momindex+2) ;
    double pz = fitparams.par()(momindex+3) ;
    double E  = fitparams.par()(momindex+4) ;
    p.r(1) = E*E-px*px-py*py-pz*pz-mass2 ;
      
    // calculate the projection matrix
    p.H(1,momindex+1) = -2.0*px ;
    p.H(1,momindex+2) = -2.0*py ;
    p.H(1,momindex+3) = -2.0*pz ;
    p.H(1,momindex+4) =  2.0*E ;
    
    return ErrCode::success ;
  }

  ErrCode 
  ParticleBase::projectConstraint(Constraint::Type atype, const FitParams&, Projection&) const 
  {
    ErrMsg(error) << "no method to project this constaint: " 
		  << name().c_str() << " " << type() << " " << atype << endmsg ;
    return ErrCode::badsetup ;
  }

  double ParticleBase::bFieldOverC()
  {
    static const BField* bfield =  gblEnv->getTrk()->magneticField();
    static const double Bz = BField::cmTeslaToGeVc*bfield->bFieldNominal() ;
    return Bz ;
  }

  ErrCode
  ParticleBase::initTau(FitParams* fitparams) const
  {
    int tauindex = tauIndex() ;
    if(tauindex>=0 && hasPosition() ) {
      const ParticleBase* amother = mother() ;
      int momposindex = amother ? amother->posIndex() : -1 ;
      int posindex = posIndex() ;
      int momindex = momIndex() ;
      assert(momposindex>=0) ; // check code logic: no mother -> no tau
      //assert(fitparams->par(momposindex+1)!=0 ||fitparams->par(momposindex+2)!=0
      //	     ||fitparams->par(momposindex+3)!=0) ; // mother must be initialized

      HepVector dX(3),mom(3) ;
      double mom2(0) ;
      for(int irow=1; irow<=3; ++irow) {
	dX(irow)  = fitparams->par(posindex+irow) - fitparams->par(momposindex+irow) ;
	double px = fitparams->par(momindex+irow) ;
	mom(irow) = px ;
	mom2 += px*px ;
      }
      double tau = dot(dX,mom)/mom2 ;
      // we don't like 0 and we don't like very negative values
      if( tau==0 ) tau=pdtTau() ;
      //tau = tau==0 ? pdtTau() : std::max(tau,-pdtTau()) ;
      fitparams->par(tauindex+1) = tau ;
    }
    return ErrCode::success ;
  }

  bool
  ParticleBase::initFromTruth(const BtaMcAssoc& truthmap, FitParams& fitparams) const
  {
    // returns true if mc match was found
    bool rc = true ;
    // first the daughters
    for( const_iterator dauit = begin() ; dauit!= end() ; ++dauit)
      rc = rc && (*dauit)->initFromTruth(truthmap,fitparams) ;

    // only if all daughters thruthmatched
    if( rc && bc() ) {
      const BtaCandidate* mcmatch = truthmap.mcFromReco(bc()) ;
      rc = mcmatch && mcmatch->pdtEntry() == bc()->pdtEntry() ;
      
      if( rc ) {
	const GTrack* gtrk = mcmatch->mcTruth()->getGTrack() ;
	int momindex = momIndex() ;
	HepLorentzVector p4 ;
	if(momindex>=0) {
	  // if this is a composite, need p4 at decay point
	  if( bc()->isComposite() ) {
	    for( std::vector<const GTrack*>::const_iterator 
		   it = gtrk->daughterList().begin() ; 
		 it != gtrk->daughterList().end(); ++it) 
	      p4 += (*it)->p4() ;
	  } else {
	    p4 = gtrk->p4() ;
	  }
	  fitparams.par(momindex+1) = p4.x() ;
	  fitparams.par(momindex+2) = p4.y() ;
	  fitparams.par(momindex+3) = p4.z() ;
	  if(hasEnergy()) fitparams.par(momindex+4) = p4.t() ;
	}
	  
	int posindex = posIndex() ;
	if(hasPosition() && posindex>=0) {
	  HepPoint decayvtx = gtrk->terminalVertex()->position() ;
	  fitparams.par(posindex+1) = decayvtx.x() ;
	  fitparams.par(posindex+2) = decayvtx.y() ;
	  fitparams.par(posindex+3) = decayvtx.z() ;

	  int tauindex = tauIndex() ;
	  if( tauindex>=0 ) {
	    // the problem is that the meaning of GVertex::time
	    // depends on whether it originates from the generator
	    // (time=t_CMS in [ns]) or from geant (time=c*t_lab in
	    // [s]). To use the generator time, we'll need the CMS
	    // boost, which we don't have here. Therefore, for decay
	    // vertices from GEANT we neglect the production time. For
	    // others we don't use the time at all.

	    double tau(0) ; // tau = tlab/gamma m = c*c*tlab/E
	    HepPoint prodvtx = gtrk->vertex()->position() ;
	    if( gtrk->terminalVertex()->cause()==GVertex::generator ) {
	      Hep3Vector p3 = p4.vect() ;
	      tau = (decayvtx - prodvtx).dot(p3)/p3.mag2() ;
	    } else {
	      tau = Constants::c * gtrk->terminalVertex()->time()/p4.t() ;
	    }

// 	    double taualt = (decayvtx.z() - prodvtx.z())/p4.z() ;
// 	    cout << "mctau: " 
// 		 << name() << " " << gtrk->terminalVertex()->cause() 
// 		 << " " << gtrk->vertex()->cause() 
// 		 << " " <<  taualt << " " << tau << endl ;

	    fitparams.par(tauindex+1) = tau ; 
	    int posindexmom(-1) ;
	    if(mother() && mother()->bc()==0 && (posindexmom=mother()->posIndex())>=0 ) {
	      fitparams.par(posindexmom+1) = prodvtx.x() ;
	      fitparams.par(posindexmom+2) = prodvtx.y() ;
	      fitparams.par(posindexmom+3) = prodvtx.z() ;
	    }
	  }
	}
      }
    }
    return rc ;
  }

}
