#include "BaBar/BaBar.hh"
#include <Beta/BtaCandidate.hh>
#include <VtxBase/VtxVertex.hh>
#include <BetaCoreTools/BtaExclusiveDecayList.hh>
#include <BetaCoreTools/BtaPrintTree.hh>
#include <ErrLogger/ErrLog.hh>

#include "VtxTreeFitter/VtxTreeFitterAlgorithm.hh"
#include "VtxTreeFitter/VtkFitter.hh"
#include "VtxTreeFitter/VtkFitParams.hh"
using std::cout;
using std::endl;

// constructor
VtxTreeFitterAlgorithm::VtxTreeFitterAlgorithm()
  : VtxAbsAlgorithm(VtxAbsAlgorithm::TreeFitter), _vertex(0)
{
}

VtxTreeFitterAlgorithm::~VtxTreeFitterAlgorithm()
{
  if(_vertex) delete _vertex ;
}

static double determinant(const HepMatrix& m)
{
  int nrow = m.num_row() ;
  if(nrow!=m.num_col()) return 0 ;
  if(nrow==1) return m(1,1) ;
  double sum(0) ;
  int fac=1;
  HepMatrix sub(nrow-1,nrow-1) ;
  //cout << "the matrix: " << m << endl ;
  // calculate determinant by walking along the first column
  for(int row=1; row<=nrow; ++row) {
    // create a matrix without this row
    for(int j=1; j<row; ++j)
      for(int col=1; col<=nrow-1; ++col)
	sub(j,col) = m(j,col+1) ;
    for(int j=row+1; j<=nrow; ++j)
      for(int col=1; col<=nrow-1; ++col)
	sub(j-1,col) = m(j,col+1) ;
    sum += fac * m(row,1) * determinant(sub) ;
    fac *= -1 ;
  }
  return sum ;
}

static int counttrajectories(const BtaCandidate& cand)
{
  int ntraj(0) ;
  
  HepAListIterator<BtaCandidate> iter=cand.daughterIterator();
  BtaCandidate* daughter ;
  while( (daughter=iter()) )
    if( daughter->recoTrk() )
      ++ntraj ;
    else if( daughter->isComposite() ) 
      if( daughter->isAResonance() )
	ntraj += counttrajectories(*daughter) ;
      else
	++ntraj ;
  return ntraj ;
}

BtaAbsVertex* VtxTreeFitterAlgorithm::compute(const BtaCandidate* decayTree) 
{
  //  setVerbose(true) ;
  if( decayTree->constraint( BtaConstraint::Beam )==0 ) {
    // count the number of daughters with a trajectory (tracks or composites)
    // if smaller than 2, add a beam spot constraint and issue an error
    int ntraj = counttrajectories(*decayTree) ;
    if(ntraj<2) {
      std::string name = decayTree->pdtEntry() ? decayTree->pdtEntry()->name() : "unknown particle" ;
      BtaPrintTree treeprinter ;
      ErrMsg(error) << "not enough trajectories for geometric fit of " << name.c_str()
		    << ". please, add beam spot constraint!" 
		    << endl << treeprinter.print( *decayTree ) << endmsg ;
      //setBeamConstraint(*decayTree) ;
      // return 0 ;
    }
  }
  
  //
  if(_vertex) delete _vertex ; 
  _vertex = new vtxtreefit::Fitter(*decayTree) ;

  if( mode()==VtxAbsAlgorithm::Fast ) 
    _vertex->fitOneStep() ;
  else
    _vertex->fit() ;
  
  //bool massconstraint = decayTree->constraint(BtaConstraint::Mass)!=0 ;
  if(verbose()) {
    cout << "verbose = " << verbose() << endl ;
    _vertex->print() ;
  }
#ifdef TREEFITDEBUG
  bool somethingwrong = _vertex->fitparams()->testCov()==false ||
    _vertex->status()==BtaAbsVertex::Failed ;
  if(somethingwrong && lund>500 && lund<600) {
    cout << "fit error code: " <<  _vertex->errCode() << endl ;
    cout << "Something went wrong fitting this vertex. Here comes some debug output." << endl ;
    cout << "fitmode: " << mode() << endl ;
    cout << "current decay point: " 
	 << decayTree->fitParams().pos() << " "
	 << decayTree->fitParams().p4() << " "
	 << decayTree->fitParams().m() << endl ;
    //_vertex->print() ;
    int vtxverbosetmp=1 ;
    std::swap(vtxverbosetmp,vtxtreefit::vtxverbose) ;
    delete _vertex ;
    //cout << "invalidating the fit ..." << endl ;
    //const_cast<BtaCandidate*>(decayTree)->invalidateFit() ;
    _vertex = new vtxtreefit::Fitter(*decayTree) ;
    if( mode()==VtxAbsAlgorithm::Fast ) 
      _vertex->fitOneStep() ;
    else
      _vertex->fit() ;
    //_vertex->print() ;
    //   assert(!somethingwrong) ;
    std::swap(vtxverbosetmp,vtxtreefit::vtxverbose) ;
  }
#endif  

  _btaFitParams = _vertex->btaFitParams(*decayTree) ;
  BtaAbsVertex* vtx = new VtxVertex(_vertex->chiSquare(),
				    _vertex->nDof(),
				    _btaFitParams.pos(),
				    _btaFitParams.posCov(),
				    _btaFitParams.xp4Cov()) ;
  vtx->setStatus(BtaAbsVertex::VtxStatus(_vertex->status())) ;
  vtx->setType(BtaAbsVertex::Geometric) ;

  return vtx ; 
}

VtxAbsAlgorithm* VtxTreeFitterAlgorithm::clone() const 
{
  ErrMsg(fatal) << "Sorry, you cann't clone a VtxTreeFitterAlgorithm yet!" << endmsg ;
  return 0 ;
}

HepLorentzVector VtxTreeFitterAlgorithm::p4() const 
{
  return _btaFitParams.p4() ;
}

BbrError VtxTreeFitterAlgorithm::p4Err()
{
  return _btaFitParams.p4Cov() ;
}

double VtxTreeFitterAlgorithm::chi2Contribution(const BtaCandidate&bc) const 
{ 
  // this operation is too expensive in the kalman filter
  return 0; 
}

BtaFitParams 
VtxTreeFitterAlgorithm::fitParams(const BtaCandidate&bc) const 
{
  return _vertex->btaFitParams(bc) ;
}

BbrDoubleErr
VtxTreeFitterAlgorithm::decayLength(const BtaCandidate&bc) const 
{
  return _vertex->decayLength(bc) ;
}

BbrDoubleErr
VtxTreeFitterAlgorithm::lifeTime(const BtaCandidate&bc) const 
{
  return _vertex->lifeTime(bc) ;
}

const BtaCandidate*
VtxTreeFitterAlgorithm::fittedCand( const BtaCandidate& cand,
				    BtaCandidate* headOfTree ) const
{
  return _vertex->fittedCand(cand,headOfTree) ;
}

BtaCandidate 
VtxTreeFitterAlgorithm::getFitted(const BtaCandidate& c ) const 
{
  return _vertex->getFitted(c) ;
}

BtaCandidate 
VtxTreeFitterAlgorithm::getFittedTree() const 
{
  return _vertex->getFittedTree() ;
}

BtaCandidate 
VtxTreeFitterAlgorithm::getFittedTree( const BtaCandidate* c) const 
{
  if(c==0) return _vertex->getFittedTree() ; 
  BtaCandidate acopy(*c) ;
  _vertex->updateTree(acopy) ;
  return acopy ;
}

const vtxtreefit::Fitter*
VtxTreeFitterAlgorithm::vtxTreeFitter() const
{
  return _vertex ;
}

