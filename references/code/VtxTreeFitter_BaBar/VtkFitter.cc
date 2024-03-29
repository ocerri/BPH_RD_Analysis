#include "BaBar/BaBar.hh"
#include <iomanip>
#include <stdio.h>
#include <sstream>

#include <BaBar/Constants.hh>
#include <BbrGeom/BbrDoubleErr.hh>
#include <BbrGeom/BbrLorentzVectorErr.hh>
#include <Beta/BtaCandidate.hh>
#include <BetaCoreTools/BtaPrintTree.hh>
#include <ErrLogger/ErrLog.hh>
#include <VtxBase/VtxVertex.hh>

#include "VtxTreeFitter/VtkFitter.hh"
#include "VtxTreeFitter/VtkFitParams.hh"
#include "VtxTreeFitter/VtkDecayChain.hh"
#include "VtxTreeFitter/VtkParticleBase.hh"
using std::cout;
using std::endl;
using std::ostream;
using std::setprecision;
using std::setw;

namespace vtxtreefit
{

  extern int vtxverbose ;
  
  Fitter::Fitter(const BtaCandidate& bc, double prec) 
    : _bc(new BtaCandidate(bc)),_decaychain(0),_fitparams(0),_status(BtaAbsVertex::UnFitted),
      _chiSquare(-1),_niter(-1), _prec(prec)
  {
    // build the tree
    if(vtxverbose>=2)
      cout << "VtkFitter::VtkFitter: building the tree" << endl ;
    _decaychain = new DecayChain(&bc,false) ;
    
    // allocate the fit parameters
    if(vtxverbose>=2)
      cout << "allocating fit parameters" << endl ;
    _fitparams  = new FitParams(_decaychain->dim()) ;
  }

  Fitter::~Fitter()
  {
    delete _bc ;
    delete _decaychain ;
    delete _fitparams ;
  }
   
  void
  Fitter::fit()
  {
    const int nitermax=10 ;
    const int maxndiverging=3 ;
    const double dChisqConv = _prec ; 
    //const double dChisqQuit = nDof() ; // if chi2 increases by more than this --> fit failed

    // initialize
    _chiSquare = -1 ;
    _errCode.reset() ;
    if( _status==BtaAbsVertex::UnFitted )
      _errCode = _decaychain->init(_fitparams) ;
    
    if(_errCode.failure()) {
      // the input tracks are too far apart
      _status = BtaAbsVertex::BadInput ;
      
    } else {
      // reset the status flag
      _status = BtaAbsVertex::UnFitted ;

      int ndiverging=0 ;
      bool finished = false ;
   
      for(_niter=0; _niter<nitermax && !finished; ++_niter) {
	HepVector prevpar = _fitparams->par() ;
	bool firstpass = _niter==0 ;
	_errCode = _decaychain->filter(_fitparams,firstpass) ;
	double chisq = _fitparams->chiSquare() ;
	double deltachisq = chisq - _chiSquare ;
	// if chi2 increases by more than this --> fit failed
	const double dChisqQuit = std::max(double(2*nDof()),2*_chiSquare) ;
	
	if(_errCode.failure()) {
	  finished = true ;
	  _status = BtaAbsVertex::Failed ;
	} else {
	  if( _niter>0 ) {
	    if( fabs( deltachisq ) < dChisqConv ) {
	      _chiSquare = chisq ;
	      _status = BtaAbsVertex::Success ;
	      finished = true ; 
	    } else if( _niter>1 && deltachisq > dChisqQuit ) {
	      _fitparams->par() = prevpar ;
	      _status  = BtaAbsVertex::Failed ;
	      _errCode = ErrCode::fastdivergingfit ;
	      finished = true ;
	    } else if( deltachisq > 0 && ++ndiverging>=maxndiverging) {
	      _fitparams->par() = prevpar ;
	      _status = BtaAbsVertex::NonConverged ;
	      _errCode = ErrCode::slowdivergingfit ;
	      finished = true ;
	    } else if( deltachisq > 0 ) {
	      // make a half step and reduce stepsize
	      _fitparams->par()   += prevpar ;
	      _fitparams->par()   *= 0.5 ;
	      //if( _niter>10) _fitparams->scale() *= 0.5 ;
	    } 
	  }
	  if ( deltachisq < 0 ) ndiverging=0 ; // start over with counting
	  if(!finished) _chiSquare = chisq ;
	} 
    	
	if(vtxverbose>=1) {
	  cout << "step, chiSquare: "
	       << setw(3) << _niter 
	       << setw(3) << _status
	       << setw(3) << nDof()
	       << setprecision(6) 
	       << setw(12) << chisq
	       << setw(12) << deltachisq << endl ;
	}
	
	if(vtxverbose>=4) {
	  print() ;
	  cout << "press a key ...." << endl ;
	  getchar() ;
	}
      }
      
      if( _niter == nitermax && _status != BtaAbsVertex::Success )
	_status = BtaAbsVertex::NonConverged ;

      //_decaychain->mother()->forceP4Sum(*_fitparams) ;

      if( !(_fitparams->testCov() ) ) {
	ErrMsg(warning) << "vtxtreefitter::Fitter: Error matrix not positive definite. "
			<< "Changing status to failed." << endmsg ;
	_status = BtaAbsVertex::Failed ;
	//print() ;
      }
    }
  }
  
  void
  Fitter::fitOneStep()
  {   
    bool firstpass = _status==BtaAbsVertex::UnFitted ;
    if( firstpass ) _decaychain->init(_fitparams) ;
    _decaychain->filter(_fitparams,firstpass) ;
    _chiSquare = _fitparams->chiSquare() ;
    if(vtxverbose>=1)
      cout << "In VtkFitter::fitOneStep(): " << _status << " " << firstpass << " " << _chiSquare << endl ;
    _status = BtaAbsVertex::Success ;
 }

  void
  Fitter::print() const
  {
    _decaychain->mother()->print(_fitparams) ;
    cout << "chisq,ndof,ncontr,niter,status: " 
	 << chiSquare() << " "
	 << nDof() << " " << _fitparams->nConstraints() << " "
	 << nIter() << " " << status() << " " << _errCode << endl ;
  } 

  void
  Fitter::printConstraints(ostream& os) const
  {
    _decaychain->printConstraints(os) ;
  }

  const HepSymMatrix& Fitter::cov() const { 
    return _fitparams->cov() ;
  }

  const HepVector& Fitter::par() const { 
    return _fitparams->par() ;
  }

  HepSymMatrix Fitter::cov(const vector<int>& indexVec) const {
    return _fitparams->cov(indexVec) ;
  }

  HepVector Fitter::par(const vector<int>& indexVec) const {
    return _fitparams->par(indexVec) ;
  }

  int
  Fitter::nDof() const {
    return _fitparams->nDof() ;
  }

  int Fitter::posIndex(const BtaCandidate* bc) const {
    return _decaychain->posIndex(bc) ;
  }
  
  int Fitter::momIndex(const BtaCandidate* bc) const {
    return _decaychain->momIndex(bc) ;
  }
  
  int Fitter::tauIndex(const BtaCandidate* bc) const {
    return _decaychain->tauIndex(bc) ;
  }

  double Fitter::add(const BtaCandidate& cand)
  {
    // first obtain a map
    //ParticleBase::indexmap indexmap ;
    //_decaychain->mother()->addToMap( indexmap ) ;
    // add this track

    ParticleBase* bp = _decaychain->mother()->addDaughter(&cand) ;
    int offset = _fitparams->dim() ;
    bp->updateIndex(offset) ;
    double deltachisq(999) ;
    if( bp ) {
      // reassign indices
      //int offset(0) ;
      //_decaychain->mother()->updIndex(offset) ;
      // resize the fitparams
      _fitparams->resize(offset) ;
      // initialize the added track, filter and return the delta chisq
      bp->initPar1(_fitparams) ;
      bp->initPar2(_fitparams) ;
      bp->initCov(_fitparams) ;

      ParticleBase::constraintlist constraints ;
      bp->addToConstraintList(constraints,0) ;
      double chisq = _fitparams->chiSquare() ;
      ErrCode status ;
      for( ParticleBase::constraintlist::const_iterator it = constraints.begin() ;
 	   it != constraints.end(); ++it)
	status |= it->filter(_fitparams) ;

      deltachisq = _fitparams->chiSquare() - chisq ;
      _chiSquare = _fitparams->chiSquare() ;

      // we want this somewhere else, but too much work now
      decaychain()->initConstraintList() ;

      //    print() ;
    } else {
      ErrMsg(warning) << "cannot add track to this vertex ..." 
		      << _decaychain->mother()->type() << endmsg ;
    }
    return deltachisq ;
  }

  double Fitter::remove(const BtaCandidate& cand)
  {
    ParticleBase* pb = const_cast<ParticleBase*>(_decaychain->locate(&cand)) ;
    ErrCode status ;
    double dchisq(0) ;
    if(pb) {
      // filter with negative weight
      ParticleBase::constraintlist constraints ;
      pb->addToConstraintList(constraints,0) ;
      double chisq = _fitparams->chiSquare() ;
      for( ParticleBase::constraintlist::iterator it = constraints.begin() ;
	   it != constraints.end(); ++it) {
	it->setWeight(-1) ;
	status |= it->filter(_fitparams) ;
      }
      dchisq = chisq - _fitparams->chiSquare() ;
      // remove
      ParticleBase* themother = const_cast<ParticleBase*>(pb->mother()) ;
      if(themother) themother->removeDaughter(pb);
    }
    return dchisq ;
  }

  void Fitter::updateIndex()
  {
    int offset=0 ;
    _decaychain->mother()->updateIndex(offset) ;
    _fitparams->resize(offset) ;
  }

  double Fitter::globalChiSquare() const 
  {
    return _decaychain->chiSquare(_fitparams) ;
  }

  BtaFitParams 
  Fitter::btaFitParams(const ParticleBase* pb) const
  {
    int posindex = pb->posIndex() ;
    // hack: for tracks and photons, use the production vertex
    if(posindex<0 && pb->mother()) posindex = pb->mother()->posIndex() ;
    int momindex = pb->momIndex() ;

    HepPoint pos(_fitparams->par()(posindex+1),
		 _fitparams->par()(posindex+2),
		 _fitparams->par()(posindex+3)) ;
    HepLorentzVector p ;
    p.setPx( _fitparams->par()(momindex+1) ) ;
    p.setPy( _fitparams->par()(momindex+2) ) ;
    p.setPz( _fitparams->par()(momindex+3) ) ;
    HepSymMatrix cov7(7,0) ;
    if( pb->hasEnergy() ) {
      // if particle has energy, get full p4 from fitparams
      p.setE( _fitparams->par()(momindex+4) ) ;
      int parmap[7] ;
      for(int i=0; i<3; ++i) parmap[i]   = posindex + i ;
      for(int i=0; i<4; ++i) parmap[i+3] = momindex + i ;    
      for(int row=1; row<=7; ++row)
	for(int col=1; col<=row; ++col)
	  cov7.fast(row,col) = _fitparams->cov()(parmap[row-1]+1,parmap[col-1]+1) ;
    } else {
      // if not, use the pdttable mass

      HepSymMatrix cov6(6,0) ;
      int parmap[6] ;
      for(int i=0; i<3; ++i) parmap[i]   = posindex + i ;
      for(int i=0; i<3; ++i) parmap[i+3] = momindex + i ;
      for(int row=1; row<=6; ++row)
	for(int col=1; col<=row; ++col)
	  cov6.fast(row,col) = _fitparams->cov()(parmap[row-1]+1,parmap[col-1]+1) ;
   
      // now fill the jacobian
      double mass = pb->pdtMass() ;
      double energy2 = mass*mass ;
      for(int row=1; row<=3; ++row) {
	double px = _fitparams->par()(momindex+row) ;
	energy2 += px*px ;
      }
      double energy = sqrt(energy2) ;
      
      HepMatrix jacobian(7,6,0) ;
      for(int col=1; col<=6; ++col)
	jacobian(col,col) = 1;
      for(int col=4; col<=6; ++col)
	jacobian(7,col) = _fitparams->par()(momindex+col)/energy ;
      
      p.setE(energy) ;
      cov7 = cov6.similarity(jacobian) ;
    }
    BtaFitParams btafitparams(pb->charge(),pos,p,cov7) ;
    btafitparams.setDecayLength(decayLength(pb)) ;
    return btafitparams ;
  }
  
  BtaFitParams 
  Fitter::btaFitParams(const BtaCandidate& cand) const 
  {
    const ParticleBase* pb = _decaychain->locate(&cand) ;
    if(pb==0) {
      ErrMsg(error) << "cann't find candidate in tree: " << cand.pdtEntry()->name() 
		    << " head of tree = " << _bc->pdtEntry()->name() 
		    << endmsg;
      return BtaFitParams() ;
    }
    return btaFitParams(pb) ;
  }

  BtaCandidate
  Fitter::getFitted() const
  {
    BtaCandidate thecand = *bc() ;
    updateCand(thecand) ;
    return thecand ;
  }

  BtaCandidate
  Fitter::getFitted(const BtaCandidate& cand) const
  {
    BtaCandidate thecand = cand ;
    updateCand(thecand) ;
    return thecand ;
  }

  BtaCandidate* 
  Fitter::fittedCand(const BtaCandidate& cand, BtaCandidate* headOfTree) const 
  {
    // assigns fitted parameters to candidate in tree
    BtaCandidate* acand = const_cast<BtaCandidate*>(headOfTree->cloneInTree(cand)) ;
    updateCand(*acand) ;
    return acand ;
  }

  BtaCandidate
  Fitter::getFittedTree() const
  {
    BtaCandidate cand = *bc() ;
    updateTree( cand ) ;
    return cand ;
  }

  HepString mybtaprint(const BtaCandidate & cand,
                     const ComIOManip & format) {
    std::ostringstream stream;
    const PdtEntry * pdtEntry = cand.pdtEntry();
    if (0 != pdtEntry){
      stream << pdtEntry->name();
      if(!cand.isComposite()) 
	stream << "(" << cand.uid() << ")" ;
      stream << std::ends;
    }
    HepString result(stream.str());
    //delete [] stream.str();           // must take ownership here
    return result;
  }

  void
  Fitter::updateCand(BtaCandidate& cand) const
  {
    // assigns fitted parameters to a candidate
    const ParticleBase* pb = _decaychain->locate(&cand) ;
    if(pb) {
      assert( pb->bc()->pdtEntry() == cand.pdtEntry() ) ;
      BtaFitParams btapar = btaFitParams(pb) ;
      VtxVertex* vtx(0) ;
      int posindex = pb->posIndex() ;
      if( posindex>=0 ) {
	if( pb ==_decaychain->cand() ) {
	  vtx = new VtxVertex(chiSquare(),nDof(),btapar.pos(),btapar.posCov(),btapar.xp4Cov()) ;
	  vtx->setStatus(BtaAbsVertex::VtxStatus(status())) ;
	  vtx->setType(BtaAbsVertex::Geometric) ;
	} else {
	  // all updated daughters are reset to unfitted, but leave the chisquare
	  double chisq = cand.decayVtx() ? cand.decayVtx()->chiSquared() : 0 ;
	  int ndof     = cand.decayVtx() ? cand.decayVtx()->nDof() : 0 ;
	  vtx = new VtxVertex(chisq,ndof,btapar.pos(),btapar.posCov(),btapar.xp4Cov()) ;
	  vtx->setStatus(BtaAbsVertex::UnFitted) ;
	}
      }
      cand.setTrajectory(btapar,vtx) ;
    } else {
      // this error message does not make sense, because it is also
      // triggered for daughters that were simply not refitted. we
      // have to do something about that.
//       BtaPrintTree printer(&mybtaprint) ;
//       ErrMsg(error) << "cann't find candidate " << endl
// 		    << printer.print(cand)
// 		    << "in tree " << endl
// 		    << printer.print(*_bc)
// 		    << endmsg;
    }
  }

  void
  Fitter::updateTree(BtaCandidate& cand) const
  {
    // assigns fitted parameters to all candidates in a decay tree
    updateCand(cand) ;
    HepAListIterator<BtaCandidate> iter=cand.daughterIterator();
    BtaCandidate* daughter ;
    while( (daughter=iter()) )  updateTree(*daughter) ;
  }

  BbrDoubleErr
  Fitter::lifeTime(const BtaCandidate& cand) const
  {
    // returns the lifetime in the rest frame of the candidate
    BbrDoubleErr rc(0,0) ;
    const ParticleBase* pb = _decaychain->locate(&cand) ;
    if(pb && pb->tauIndex()>=0 && pb->mother()) {
      int tauindex = pb->tauIndex() ;
      double tau    = _fitparams->par()(tauindex+1) ;
      double taucov = _fitparams->cov()(tauindex+1,tauindex+1) ;
      double mass   = pb->pdtMass() ; 
      double convfac = mass/Constants::c ;
      rc = BbrDoubleErr(convfac*tau,convfac*convfac*taucov) ;
    }
    return rc ;
  }

  BbrDoubleErr
  Fitter::decayLength(const ParticleBase* pb) const 
  {
    // returns the decaylength in the lab frame
    return decayLength(pb,_fitparams) ;
  }


  BbrDoubleErr
  Fitter::decayLength(const ParticleBase* pb,
		      const FitParams* fitparams)
  {
    // returns the decaylength in the lab frame
    BbrDoubleErr rc(0,0) ;
    if(pb->tauIndex()>=0 && pb->mother()) {
      // one can calculate the error in many ways. I managed to make
      // them all agree, with a few outliers. this one seems to be
      // most conservative/stable/simple one.
 
      // len = tau |mom|
      int tauindex = pb->tauIndex() ;
      int momindex = pb->momIndex() ;
      double tau    = fitparams->par()(tauindex+1) ;
      double mom2(0) ;
      for(int row=1; row<=3; ++row) {
	double px = fitparams->par()(momindex+row) ;
	mom2 += px*px ;
      }
      double mom = sqrt(mom2) ;
      double len = mom*tau ;

      vector<int> indexvec ;
      indexvec.push_back(tauindex) ;
      indexvec.push_back(momindex+0) ;
      indexvec.push_back(momindex+1) ;
      indexvec.push_back(momindex+2) ;
      
      HepVector jacobian(4) ;
      jacobian(1) = mom ;
      jacobian(2) = tau*fitparams->par()(momindex+1)/mom ;
      jacobian(3) = tau*fitparams->par()(momindex+2)/mom ;
      jacobian(4) = tau*fitparams->par()(momindex+3)/mom ;

      rc = BbrDoubleErr(len,fitparams->cov(indexvec).similarity(jacobian)) ;
    }
    return rc ;
  }

  BbrDoubleErr
  Fitter::decayLength(const BtaCandidate& cand) const
  {
    BbrDoubleErr rc(0,0) ;
    const ParticleBase* pb = _decaychain->locate(&cand) ;
    if(pb && pb->tauIndex()>=0 && pb->mother()) rc = decayLength(pb) ;
    return rc ;
  }  

  BbrDoubleErr
  Fitter::decayLengthSum(const ParticleBase* pbA, const ParticleBase* pbB) const 
  {
    // returns the decaylengthsum in the lab frame
    BbrDoubleErr rc(0,0) ;
    if(pbA->tauIndex()>=0 && pbA->mother() &&
       pbB->tauIndex()>=0 && pbB->mother() ) {

      // len = tau |mom|
      int tauindexA = pbA->tauIndex() ;
      int momindexA = pbA->momIndex() ;
      double tauA    = _fitparams->par()(tauindexA+1) ;
      double mom2A(0) ;
      for(int row=1; row<=3; ++row) {
	double px = _fitparams->par()(momindexA+row) ;
	mom2A += px*px ;
      }
      double momA = sqrt(mom2A) ;
      double lenA = momA*tauA ;

      int tauindexB = pbB->tauIndex() ;
      int momindexB = pbB->momIndex() ;
      double tauB    = _fitparams->par()(tauindexB+1) ;
      double mom2B(0) ;
      for(int row=1; row<=3; ++row) {
	double px = _fitparams->par()(momindexB+row) ;
	mom2B += px*px ;
      }
      double momB = sqrt(mom2B) ;
      double lenB = momB*tauB ;

      vector<int> indexvec ;
      indexvec.push_back(tauindexA) ;
      for(int i=0; i<3; ++i) indexvec.push_back(momindexA+i) ;
      indexvec.push_back(tauindexB) ;
      for(int i=0; i<3; ++i) indexvec.push_back(momindexB+i) ;
      
      HepVector jacobian(8) ;
      jacobian(1) = momA ;
      for(int irow=1; irow<=3; ++irow) 
	jacobian(irow+1) = tauA*_fitparams->par()(momindexA+irow)/momA ;
      jacobian(5) = momB ;
      for(int irow=1; irow<=3; ++irow) 
	jacobian(irow+5) = tauB*_fitparams->par()(momindexB+irow)/momB ;
      
      rc = BbrDoubleErr(lenA+lenB,cov(indexvec).similarity(jacobian)) ;
    }
    return rc ;
  }
  
  BbrDoubleErr
  Fitter::decayLengthSum(const BtaCandidate& candA, const BtaCandidate& candB) const
  {
    BbrDoubleErr rc(0,0) ;
    const ParticleBase* pbA = _decaychain->locate(&candA) ;
    const ParticleBase* pbB = _decaychain->locate(&candB) ;
    if(pbA && pbB)  rc = decayLengthSum(pbA,pbB) ;
    return rc ;
  }  

}
  
