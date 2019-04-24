#include "BaBar/BaBar.hh"
#include <Beta/BtaCandidate.hh>
#include <BbrGeom/BbrDoubleErr.hh>
#include <BaBar/Constants.hh>

#include "VtxTreeFitter/VtkUpsilonFitter.hh"
#include "VtxTreeFitter/VtkUpsilon.hh"
#include "VtxTreeFitter/VtkDecayChain.hh"
#include "VtxTreeFitter/VtkFitParams.hh"
#include "VtxTreeFitter/VtkBtaInterface.hh"
using std::cout;
using std::endl;


namespace vtxtreefit
{

  UpsilonFitter* 
  UpsilonFitter::createFromUpsilon(const BtaCandidate& constups, 
			       const BtaCandidate& constrecB,
			       const EventInfo* eventinfo) 
  {
    // first set the constraints on the upsilon ( the contents can be
    // copied by the setConstraint operation!!)
    BtaCandidate* ups = const_cast<BtaCandidate*>(&constups) ;
    setEnergyConstraint(*ups,eventinfo) ;
    setBeamConstraint(*ups,eventinfo) ;

    // now find the B daughters
    HepAListIterator<BtaCandidate> iter=ups->daughterIterator();
    BtaCandidate* recB = iter() ;
    BtaCandidate* tagB = iter() ;
    if( recB==0 || tagB==0 ) {
      cout << "big trouble: " << recB << " " << tagB << endl ;

    }
    if( !recB->isCloneOf(constrecB) ) std::swap(recB,tagB) ;
    
    // set the B constraints 
    recB->invalidatePresentFit() ;
    tagB->invalidatePresentFit() ;
    // these constraints are not okay if candidate way of in dE/mES
    // sideband, because fits will fail. leave up to user to include
    // them.
    //setMassConstraint(*recB) ;
    //setMassConstraint(*tagB) ;
    
    // remove lifetime constraints from the B
    
    // add a missing particle to the tag B
    BtaInterface op ;
    op.addMissingParticle(*tagB) ;
      
    // now create a vertex
    return new UpsilonFitter(*ups,*recB) ;
  }   
  
  UpsilonFitter::UpsilonFitter(const BtaCandidate& ups,const BtaCandidate& recB) 
    : Fitter(ups),_recB(0),_tagB(0) 
  { 
    // constructor from an upsilon
    _recB = decaychain()->locate(&recB) ;
    Upsilon* vtkups = getUpsilon() ;
    _tagB = vtkups->daughters()[1] ;
    if(_tagB==_recB) {
      // daughters are not in correct order, swap them:
      std::swap(vtkups->daughters()[0],vtkups->daughters()[1]) ;
      _tagB = vtkups->daughters()[1] ;
    }
    assert(_tagB && _recB ) ;
  }


  UpsilonFitter::UpsilonFitter(const BtaCandidate& cand)  
    : Fitter(cand), _recB(0),_tagB(0) 
  {
    // constructor from a reco-B
    const Upsilon* ups = getUpsilon() ;
    _recB = ups->daughters()[0] ;
    _tagB = ups->daughters()[1] ;
  }
 
  Upsilon* UpsilonFitter::getUpsilon()
  {
    Upsilon* rc(0) ;
    if( decaychain()->mother()->type() == ParticleBase::kUpsilon )
      rc = static_cast<Upsilon*>(decaychain()->mother()) ;
    return rc ;
  } 
  
  const Upsilon* UpsilonFitter::getUpsilon() const
  {
    const Upsilon* rc(0) ;
    if( decaychain()->mother()->type() == ParticleBase::kUpsilon )
      rc = static_cast<const Upsilon*>(decaychain()->mother()) ;
    return rc ;
  }  

  BtaCandidate UpsilonFitter::getTaggingB() const
  {
    const Upsilon* ups = getUpsilon() ;
    if(!ups) return BtaCandidate() ; 
    const ParticleBase* tagB = ups->tagB() ;
    BtaCandidate thecand(*(tagB->bc())) ;
    updateCand( thecand ) ;
    return thecand ;
  }

  BbrDoubleErr UpsilonFitter::deltaT() const 
  {
    return lifeTimeSum(true) ;
  }

  BbrDoubleErr UpsilonFitter::deltaZ() const 
  {
    BbrDoubleErr rc(-999,-999) ;
    int indexCP   = _recB ? _recB->posIndex() : -1 ;
    int indexTag  = _tagB ? _tagB->posIndex() : -1 ;
    if( indexCP>=0 && indexTag>=0 ) {
      vector<int> indexvec ;
      indexvec.push_back(indexCP+2) ;
      indexvec.push_back(indexTag+2) ;
      HepVector par    = fitparams()->par(indexvec) ;
      HepSymMatrix cov = fitparams()->cov(indexvec) ;
      HepVector jacobian(2) ;
      jacobian(1) = 1 ;
      jacobian(2) = -1 ;
      double deltaZ    = dot(par,jacobian) ;
      double covDeltaZ = cov.similarity(jacobian) ;
      rc = BbrDoubleErr(deltaZ,covDeltaZ) ;
      //cout << "something is right!" << endl ;
    } else {
      //cout << "something is wrong: " << indexCP << " " << indexTag << endl ;
      //      assert(0) ;
    }
    return rc ;
  }

  BbrDoubleErr UpsilonFitter::lifeTimeSum(bool difference) const 
  {
    BbrDoubleErr rc(-999,-999) ;
    int indexCP   = _recB ? _recB->tauIndex() : -1 ;
    int indexTag  = _tagB ? _tagB->tauIndex() : -1 ;
    if( indexCP>=0 && indexTag>=0 ) {
      vector<int> indexvec ;
      indexvec.push_back(indexCP) ;
      indexvec.push_back(indexTag) ;
      HepVector par    = fitparams()->par(indexvec) ;
      HepSymMatrix cov = fitparams()->cov(indexvec) ;
      HepVector jacobian(2) ;
      jacobian(1) = 1 ;
      jacobian(2) = difference ? -1 : 1 ;
      double tsum    = dot(par,jacobian) ;
      double covtsum = cov.similarity(jacobian) ;
      double mass = _recB->pdtMass() ;
      double convfac = 1e12*mass/Constants::c ; 
      rc = BbrDoubleErr(convfac*tsum,convfac*convfac*covtsum) ;
    }
    return rc ;
  }

  void UpsilonFitter::replaceTagB(Fitter& vertex, bool lifetimesumconstraint)
  {
    // this is very complicated :(

    // 1) identify the tag b in the new vertex (should be an ExternalBtaParticle)
 
    //cout << "BEFORE swapping: " << endl ;
    //vertex.print() ;

    InternalParticle* pb = static_cast<InternalParticle*>(vertex.decaychain()->mother()) ;
    const ParticleBase* theoldb(0) ; 
    for( InternalParticle::daucontainer::const_iterator 
	   it = pb->daughters().begin() ;
	 it != pb->daughters().end() ; ++it ) {
      if( (*it)->type() == ParticleBase::kBtaComposite ||
	  (*it)->type() == ParticleBase::kBtaResonance) 
	theoldb = *it ;
    }
    pb->swapMotherDaughter(vertex.fitparams(), theoldb) ;
    pb->setMassConstraint(true) ;
    pb->setMother( decaychain()->mother() ) ;
    
    ParticleBase::indexmap rhsindexmap ;
    pb->retrieveIndexMap(rhsindexmap) ;
    // remove the mother from the index map (changed dimension ...)
    rhsindexmap.pop_back()  ;

    //cout << "AFTER swapping: " << endl ;
    //vertex.print() ;

    // 2) identify the tag b in the current vertex and delete it (should be a MissingParticle)
    Upsilon* upsilon = getUpsilon() ;
    ParticleBase* lhstagb = upsilon->daughters().back() ;
  
    // 3) create an indexmap for the current vertex and also a parameter set
    ParticleBase::indexmap thisindexmap ;
    upsilon->retrieveIndexMap(thisindexmap) ;
    FitParams thisfitparams = *fitparams() ;

    // 4) add the new tree, transfer ownership
    upsilon->daughters().pop_back() ;
    upsilon->daughters().push_back(pb) ;
    vertex.decaychain()->setOwner(false) ;
    _tagB = pb ;

    // copy both sets of parameters
    updateIndex() ;
    fitparams()->copy(thisfitparams,thisindexmap) ;
    fitparams()->copy(*(vertex.fitparams()),rhsindexmap) ;

    // delete the old tagging B
    delete lhstagb ;
  
    // after all the copying ...
    //     cout << "after copuying: " << endl ;
    //     print() ;

    //vtxverbose=1 ;
    // we have to add the new constraints. this is the easy way:
    upsilon->setLifetimeSumConstraint(lifetimesumconstraint) ;
    decaychain()->initConstraintList() ;
    fit() ;
    //    cout << "after fitting: " << endl ;
    //     print() ;
    
    //vtxverbose=0 ;
  }
}
