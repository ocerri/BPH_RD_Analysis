#include <BaBar/BaBar.hh>
#include <PDT/Pdt.hh>
#include <BbrGeom/BbrLorentzVectorErr.hh>
#include <BbrGeom/BbrDoubleErr.hh>
#include <Beta/BtaCandidate.hh>
#include <BetaCoreTools/BtaOpMakeTree.hh>
#include <VtxBase/VtxVertex.hh>
#include "VtxTreeFitter/VtkBtaInterface.hh"
#include "VtxTreeFitter/VtkFitParams.hh"
#include "VtxTreeFitter/VtkParticleBase.hh"


namespace vtxtreefit
{

  void
  BtaInterface::addMissingParticle(BtaCandidate& tagB)
  {
    // adds a missing particle to a tagB if it isn't there yet
    bool foundmissing(false) ;
    HepAListIterator<BtaCandidate> iter=tagB.daughterIterator();
    BtaCandidate* daughter ;
    while( !foundmissing && (daughter=iter() ) ) {
      foundmissing  = daughter->recoTrk()==0 && daughter->recoCalo()==0 &&
	daughter->isComposite()==false ;
    }
    if(!foundmissing) { 
      HepLorentzVector missingmom(0,0,1,1.1) ;
      //static PdtEntry* pion = Pdt::lookup(PdtLund::pi0) ;
      BtaCandidate* missing = new BtaCandidate(missingmom,0) ;
      addDaughterLink(tagB,missing) ;
      // addDaughterLink makes a copy
      delete missing ; 
      //cout << "this candidate already had a missing daughter!" << endl ;
    }
  }

  BtaCandidate*
  BtaInterface::createTree(const ParticleBase& pb)
  {
    // creates a tree from scratch, using the tree structure in the
    // particle base. only useful if tracks have been added or removed from the vertex.
    // obviously this works recursively. 
    // the main problem is that we want to copy all constraints. in
    // fact, that is sort of a showstopper.

    //we don't actually fill any pars here. too much work. we also
    //don't create vertices. problem?
    const BtaCandidate* cand = pb.bc() ;
    BtaCandidate* rc(0) ;

    //cout << "in BtaInterface::createTree: " << pb.name() << endl ;
    if( pb.begin()==pb.end() ) { // no daughters
      
      if( cand ) {
	// there already exists a BtaCandidate. in this case, just
	// make a shallow copy
	
	//const BtaCandBase* bcb = contentOf(*cand) ;
	//BtaCandBase* newbase = bcb->clone() ;
	rc = new BtaCandidate(*cand) ;
      } else {
	// create a missing particle ..
	rc = new BtaCandidate() ;
	rc->setType(Pdt::lookup(PdtLund::null)) ;
      }

    } else {
      // collect all daughters
      HepAList<BtaCandidate> alist ;
      for(ParticleBase::const_iterator idau = pb.begin() ;
	  idau !=  pb.end() ; ++idau) 
	alist += createTree(**idau) ;
      
      // make the new particle
      static BtaOpMakeTree comb;
      HepAListIterator<BtaCandidate> iter(alist) ;
      rc = comb.createFromList(iter) ;
      
      if(cand) {
	// copy constraints
	for(int i=0; i<cand->nConstraints(); ++i)
	  rc->addConstraint( * (cand->constraint(i)) ) ;
	// copy the type
	contentOf(*rc)->setCharge(cand->pdtEntry()->charge()) ;
	rc->setType(cand->pdtEntry()) ;
      }

      // delete all daughters
      //cout << "big memory leak here!" << endl ;
      HepAListDeleteAll(alist) ;
    }
    return rc ;
  }

}
