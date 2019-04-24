#include "BaBar/BaBar.hh"
#include <algorithm>
#include <iomanip>
#include <list>
//#include <pair>

#include <BbrGeom/BbrDoubleErr.hh>
#include <Beta/BtaCandidate.hh>
#include <TrkBase/TrkCompTrk.hh>
#include <TrkBase/TrkDifTraj.hh>
#include <TrkBase/TrkPoca.hh>
#include <AbsEnv/AbsEnv.hh>
#include <BtaEnv/BtaEnv.hh>
#include <BetaCoreTools/BtaBVariables.hh>
#include <BetaCoreTools/BtaExclusiveDecayTruthComp.hh>
#include <BetaCoreTools/BtaMcAssoc.hh>
#include <BetaCoreTools/BtaOpMakeTree.hh>
#include <ProxyDict/Ifd.hh>
#include <ProxyDict/IfdStrKey.hh>
#include <PDT/Pdt.hh>
#include "AbsEvent/AbsEvent.hh"
#include "AbsEventTag/AbsEventTag.hh"

#include <BbrStdUtils/BbrCollectionUtils.hh>
#include "BetaMicroBase/BtaAttributes.hh"
#include "BetaCoreTools/BtaDeltaTConverter.hh"
#include "Beta/EventInfo.hh"
#include "BetaCoreTools/BtaBooster.hh"
#include "ProbTools/ChisqConsistency.hh"

#include "VtxTreeFitter/VtkTrkAbsFitWrapper.hh"
#include "VtxTreeFitter/VtkTagVertexer.hh"
#include "VtxTreeFitter/VtkUpsilonFitter.hh"
#include "VtxTreeFitter/VtkUpsilon.hh"
#include "VtxTreeFitter/VtkBtaInterface.hh"
using std::cout;
using std::endl;
using std::setprecision;
using std::setw;

namespace vtxtreefit
{
  const unsigned int TagVertexer::gMaxNVertices = 4 ;
  const double TagVertexer::gMaxChisqTrackAdd = 8 ; // in the limit (N+1)/N--> 1, on two dofs
  const double TagVertexer::gMaxSeedDeltaZ = 1; 
  const double TagVertexer::gMaxSeedChisq = 25 ;

  class TagTrkAbsFitWrapper : public TrkAbsFitWrapper
  {
  public:
    TagTrkAbsFitWrapper(const BtaCandidate* cand, int mcmatch) 
      : TrkAbsFitWrapper(cand),_mcmatch(mcmatch) {} ;
    TagTrkAbsFitWrapper(const BtaCandidate* cand, const HepVector& par,
		     const HepSymMatrix& cov, double charge) 
      : TrkAbsFitWrapper(cand,par,cov,charge),_mcmatch(0) {}

    int mcmatch() const { return _mcmatch ; }
  private:
    int _mcmatch ;
  } ;

  class VtkPairPoca
  {
  public:
    VtkPairPoca(const TagTrkAbsFitWrapper& cand1, const TagTrkAbsFitWrapper& cand2) ;
    double chisq() const { return _chisq ; }
    double doca() const { return _doca ; }
    double docaErr() const { return _doca/sqrt(_chisq) ; }
    int success() const { return _poca.status().success() ; }
    const HepVector& vertexPos() const { return  _vertexPos ; }
    const HepSymMatrix& vertexCov() const { return  _vertexCov ; }
    void print() const ;
    template<class T> static HepVector convert(const T& vec) ;
    void calcVertex() ;
    double quality() const { return _quality ; }
    bool isFitted() const { return _isFitted ; }
    const TagTrkAbsFitWrapper* first() const { return _cand[0] ; }
    const TagTrkAbsFitWrapper* second() const { return _cand[1] ; }
    const HepVector& pos1() const { return _pos[0] ; }

  private:
    const TagTrkAbsFitWrapper* _cand[2] ;
    double _flt[2] ;
    HepVector _pos[2] ;
    HepSymMatrix _cov[2] ;
    TrkPoca _poca ;
    double _chisq ;
    double _doca ;
    double _vardoca[2] ;
    HepVector _vertexPos ;
    HepSymMatrix _vertexCov ;
    double _quality ;
    bool _isFitted ;
  } ;

  template<class T>
  HepVector VtkPairPoca::convert(const T& vec) 
  {
    HepVector rc(3) ;
    rc(1) = vec.x() ;
    rc(2) = vec.y() ;
    rc(3) = vec.z() ;
    return rc ;
  }

  VtkPairPoca::VtkPairPoca(const TagTrkAbsFitWrapper& cand1, const TagTrkAbsFitWrapper& cand2)
    : _poca(cand1.traj(),0.0,cand2.traj(),0.0,1e-5), _chisq(-1),
      _vertexPos(3),_vertexCov(3),_quality(-999),_isFitted(false) 
  {
    _cand[0] = &cand1 ;
    _cand[1] = &cand2 ;
    if( _poca.status().success() ) {
      _flt[0] = _poca.flt1() ;
      _flt[1] = _poca.flt2() ;
      for(int i=0; i<2; ++i) {
	BbrPointErr position = _cand[i]->trkFit()->positionErr(_flt[i]) ;
	_cov[i] = position.covMatrix() ;
	_pos[i] = convert<HepPoint>(position) ;
      }
      HepVector residual = _pos[1] - _pos[0] ;
      _doca = sqrt( dot(residual,residual) ) ;
      HepVector& jacobian = residual ;
      jacobian /= _doca ;
      for(int i=0; i<2; ++i)
	_vardoca[i] = _cov[i].similarity(jacobian) ;
      _chisq = _doca*_doca/(_vardoca[0]+_vardoca[1]) ;
      // _vertexPos = 0.5*(_pos[1]+_pos[0]) ;
      // this is VERY wrong:
      //_chisq = residual.determineChisq(Hep3Vector(0,0,0)) ;
    }
  }
  
  void
  VtkPairPoca::calcVertex()
  {
    if( _poca.status().success() ) {
      HepVector dChi2dPos(3,0) ;
      HepSymMatrix d2Chi2dPos2(3,0) ;
      HepVector residual = _pos[1] - _pos[0] ;
      for(int i=0; i<2; ++i) {
	HepVector dir = convert<Hep3Vector>(_cand[i]->trkFit()->direction(_flt[i]) ) ;
	HepSymMatrix P(3,1) ;
	for(int row=1; row<=3; ++row)
	  for(int col=1; col<=row; ++col) P.fast(row,col) -= dir(row)*dir(col) ;
	P /= _vardoca[i] ;
	d2Chi2dPos2 += P ;
	dChi2dPos   += P*_pos[i] ;
      }
      int ierr;
      d2Chi2dPos2.invert(ierr) ;
      _vertexCov = d2Chi2dPos2 ;
      _vertexPos = d2Chi2dPos2*dChi2dPos ;
      _quality = -_chisq - log(_vertexCov.determinant()) ;
      //_quality = -_chisq - log(_vertexCov.fast(3,3)) ;
      _isFitted = true ;
    }
  }


  void VtkPairPoca::print() const {
    cout << setprecision(5) 
	 << "x: " << setw(12) << _vertexPos(1) << " +- " << setw(12) << sqrt(_vertexCov(1,1)) << endl 
	 << "y: " << setw(12) << _vertexPos(2) << " +- " << setw(12) << sqrt(_vertexCov(2,2)) << endl 
	 << "z: " << setw(12) << _vertexPos(3) << " +- " << setw(12) << sqrt(_vertexCov(3,3)) << endl ;
  }

  
  inline bool sortQuality(const VtkPairPoca* lhs, const VtkPairPoca* rhs)
  {
    return lhs->quality() > rhs->quality() ;
  }

  class TrackVertexMatch
  {
  public:
    TrackVertexMatch(const TagTrkAbsFitWrapper* track)
      : _track(track),_chisq(999),_logV(999) {}
    void setChisq(double chisq, double logV) { _chisq = chisq ; _logV = logV ; }
    double chisq() const { return _chisq ; }
    double density() const { return -_chisq -_logV ; }
    bool operator<(const TrackVertexMatch& rhs) const {
      return density() > rhs.density() ; 
    }
    const TagTrkAbsFitWrapper* track() const { return _track ; }
  private:
    const TagTrkAbsFitWrapper* _track ;
    double _chisq ;
    double _logV ;
  } ;

  class TagFastVertex : public BtaOperatorBase
  {
  public:
    TagFastVertex() : _vertex(0) {}
    TagFastVertex(const VtkPairPoca& poca) ;
    ~TagFastVertex() { delete _vertex ; delete _seedcand ; }
  
    int status() const { return  _vertex->status() ; }
    TrackVertexMatch deltaChiSq(const TagTrkAbsFitWrapper& cand) const ;
    double add(const TagTrkAbsFitWrapper& cand) ;
    int size() const { return _tracks.size() ; }
    int mcmatch() const { return _tracks.front()->mcmatch() ; }
    int ngood() const ;
    HepSymMatrix cov() const { return _vertex->cov(_indexvec) ; }
    HepVector    par() const { return _vertex->par(_indexvec) ; }
    double chiSquare() const { return _vertex->chiSquare() ; }
    int nDof() const { return _vertex->nDof() ; }
    double chisqprob() const { 
      ChisqConsistency x(chiSquare(),nDof()) ;
      return x.significanceLevel() ;
      //return TMath::Prob(chiSquare(),nDof()) ; } 
    }
    
    BtaCandidate* createTagCandidate() const ;
    Fitter* vertex() { return _vertex ; }
    const VtkPairPoca* seedpoca() const { return _seedpoca ; }
    double quality() const { 
      //return chisqprob()/sqrt( cov().fast(3,3) ) ; } 
      return chisqprob()/cov().fast(3,3) ; }  //         <-- slightly better (89%)
    //return 1./( cov().fast(3,3) * chiSquare() / nDof() ) ; } // <-- the worst (85%)
    //return chisqprob() ; }
    //return 1/cov().fast(3,3) ; }   <--- worse than the worst
    


  private:
    const VtkPairPoca* _seedpoca ;
    Fitter*      _vertex ;
    BtaCandidate* _seedcand ;
    std::vector<int>     _indexvec ;
    std::vector<const TagTrkAbsFitWrapper*> _tracks ;
  } ;

  TagFastVertex::TagFastVertex(const VtkPairPoca& poca)
    : _seedpoca(&poca),_vertex(0),_seedcand(0)
  {
    const TagTrkAbsFitWrapper* tagBseed = poca.first() ;
    const TagTrkAbsFitWrapper* trkseed  = poca.second() ;
    static BtaOpMakeTree comb;
    _seedcand = comb.create(*(tagBseed->cand()),*(trkseed->cand())) ;
    _seedcand->setType(tagBseed->cand()->pdtEntry());
    _vertex = new Fitter(*_seedcand) ;
    _vertex->fit() ;
    //_vertex->print() ;
    int posindex = _vertex->posIndex(_seedcand) ;
    for(int i=0; i<3; ++i) _indexvec.push_back(posindex+i) ;
    _tracks.push_back(trkseed) ;
  }

 
  TrackVertexMatch TagFastVertex::deltaChiSq(const TagTrkAbsFitWrapper& cand) const
  {
    TrackVertexMatch match(&cand) ;
    HepVector    vertexpos = _vertex->par(_indexvec) ;
    HepSymMatrix vertexcov = _vertex->cov(_indexvec) ;
    HepPoint vertexpt(vertexpos(1),vertexpos(2),vertexpos(3)) ;
    TrkPoca thePoca(cand.traj(),0,vertexpt,1e-5);
    if( thePoca.status().success() ) {
      BbrPointErr trkposerr = cand.trkFit()->positionErr(thePoca.flt1()) ;
      HepVector trkpos = VtkPairPoca::convert<HepPoint>(trkposerr) ;
      const HepSymMatrix& trkcov = trkposerr.covMatrix() ;
      HepVector residual = trkpos-vertexpos ;
      double doca = sqrt(dot(residual,residual)) ;
      HepVector jacobian = residual ;
      jacobian /= doca ;
      double vardoca = (trkcov + vertexcov).similarity(jacobian) ;
      match.setChisq( doca*doca/vardoca , log(vardoca) ) ;
    }
    return match ;
  }

  double TagFastVertex::add(const TagTrkAbsFitWrapper& cand)
  {
    _tracks.push_back(&cand) ;
    return _vertex->add(*(cand.cand())) ;
  }
  
  int TagFastVertex::ngood() const
  {
    int rc(0) ;
    for( std::vector<const TagTrkAbsFitWrapper*>::const_iterator it = _tracks.begin() ;
	 it != _tracks.end() ; ++it)
      if( (*it)->mcmatch()==2) ++rc ;
    return rc ;
  }

  BtaCandidate* TagFastVertex::createTagCandidate() const 
  {
    // collect the tracks
    HepLorentzVector missingmom(_seedpoca->first()->cand()->p4()) ;
    HepAList<BtaCandidate> alist ;
    for( std::vector<const TagTrkAbsFitWrapper*>::const_iterator it = _tracks.begin() ;
	 it != _tracks.end() ; ++it) {
      alist     += const_cast<BtaCandidate*>((*it)->cand()) ;
      missingmom -= (*it)->cand()->p4() ;
    }
    // add a missing particle
    BtaCandidate missingcand(missingmom) ;
    alist += &missingcand ;

    // make a composite for the tagging B
    static BtaOpMakeTree comb;
    HepAListIterator<BtaCandidate> iter(alist) ;
    BtaCandidate* tagB = comb.createFromList(iter) ;
    tagB->setType(_seedpoca->first()->cand()->pdtEntry()) ; 
    return tagB ;
  }

  TagVertexer::~TagVertexer() { 
    if(_vertex) delete _vertex ; 
    if(_fastvertex) delete _fastvertex ; 
  }

  TagVertexer::TagVertexer(const BtaCandidate& bc,
			   const HepAList<BtaCandidate>& bTagList,
			   const EventInfo* eventinfo,
			   const HepAList<BtaCandidate>& mclist,
			   const BtaMcAssoc* truthmap,
			   const AbsEvent* absevent) 
    : _vertex(0), _fastvertex(0), _recoB(bc), _verbosity(0)
  {

    // fit the reco B with the upsilon constraint
    setBeamConstraint(_recoB,eventinfo) ;
    setEnergyConstraint(_recoB,eventinfo) ;
    setMassConstraint(_recoB) ;
    _recoB.invalidateFit() ;

    if(_verbosity>=1)
      cout << "In the tagvertexer: " << endl ;
    _vertex = new UpsilonFitter(_recoB) ;
    _vertex->fit() ;

    //_vertex->print() ;
    double upschisqprob = ChisqConsistency(_vertex->chiSquare(),
					   _vertex->nDof()).significanceLevel() ;
    const double minupschisqprob = 0.001 ;

    BtaCandidate fitted_recoB = _vertex->getFitted(_recoB) ;
    BtaBVariables bVariables(fitted_recoB.p4WCov()) ;
    //double mES    = bVariables.energySubstitutedMass();
    //double deltaE = bVariables.deltaE();
    double zCP    = fitted_recoB.decayVtx()->point().z() ;

    // 
    BtaExclusiveDecayTruthComp truthTool(&_recoB,const_cast<HepAList<BtaCandidate>*>(&mclist)) ;
    const Upsilon* upsilon = _vertex->getUpsilon() ;
    const BtaCandidate* upsbc = truthTool.getUpsilon() ;
    const BtaCandidate* mctagB = truthTool.getTagB() ;
    if(true || (upschisqprob>minupschisqprob &&  upsilon && upsbc && mctagB)) {

      // Get fitted information on Upsilon and the tagging B
      int posindex = upsilon->posIndex() ;
      std::vector<int> indexvec ;
      for(int i=0; i<3; ++i) indexvec.push_back(posindex+i) ;
      int dirindex = upsilon->daughters()[1]->momIndex() ;
      for(int i=0; i<3; ++i) indexvec.push_back(dirindex+i) ;
      
      HepVector    par = _vertex->par(indexvec) ;
      HepSymMatrix cov = _vertex->cov(indexvec)  ;
      
      // get the truth
      HepPoint upspoint   = upsbc->decayVtx()->point();
      HepLorentzVector p4 = mctagB->p4() ;
      
#ifdef VTK_MONITOR
      static TNtuple* nt = 
	new TNtuple("upsnt","","dx:dy:dz:sx:sy:sz:dpx:dpy:dpz:spx:spy:spz:vtkchisq:deltaE") ;
      double dx = par(1) - upspoint.x() ;
      double sx = sqrt(cov.fast(1,1)) ;
      double dy = par(2) - upspoint.y() ;
      double sy = sqrt(cov.fast(2,2)) ;
      double dz = par(3) - upspoint.z() ;
      double sz = sqrt(cov.fast(3,3)) ;
      double dpx = par(4) - p4.x() ;
      double spx = sqrt(cov.fast(4,4)) ;
      double dpy = par(5) - p4.y() ;
      double spy = sqrt(cov.fast(5,5)) ;
      double dpz = par(6) - p4.z() ;
      double spz = sqrt(cov.fast(6,6)) ;
      nt->Fill(dx,dy,dz,sx,sy,sz,dpx,dpy,dpz,spx,spy,spz,
	       _vertex->chiSquare(),deltaE) ;
      
      static TNtuple* bsnt = 
	new TNtuple("bsnt","","dx:dy:dz:sx:sy:sz") ;
      BbrPointErr beamSpot(eventinfo->beamSpot());
      bsnt->Fill( beamSpot.x()-upspoint.x(),
		  beamSpot.y()-upspoint.y(),
		  beamSpot.z()-upspoint.z(),
		  sqrt(beamSpot.covMatrix()(1,1)),
		  sqrt(beamSpot.covMatrix()(2,2)),
		  sqrt(beamSpot.covMatrix()(3,3)) ) ;
#endif


      // let's see what we can learn from the truth
      BtaCandidate tagBbc = _vertex->getTaggingB() ;
      TagTrkAbsFitWrapper tagBabsfit(&tagBbc,par,cov,0) ;
      
      HepPoint tagBvtx = mctagB->decayVtx()->point();
      //cout << "true decaypoint of tagB: " << tagBvtx << endl ;
      HepAListIterator<BtaCandidate> iter(bTagList) ;

      //      cout << "From the tag vertexer: " << endl ;
      //       Clusterer clusterer( bTagList, &tagBabsfit ) ;
  
      BtaCandidate* thiscand ;
      
      int ngood(0) ;
      std::vector<TagTrkAbsFitWrapper*> absfitlist ;
      std::vector<VtkPairPoca*>  pocalist ;
      
      while( (thiscand = iter() )) 
	if(thiscand->charge()!=0) {
	  const BtaCandidate* mccand = truthmap->mcFromReco(thiscand) ;
	  int mcmatch = 0 ;
	  if( mccand ) {
	    mcmatch = 1 ;
	    //cout << "mc particle: " << mccand->pdtEntry()->name() << endl ;
	    if( mccand->theMother() ) {
	      HepPoint prodvtx = mccand->theMother()->decayVtx()->point();
	      // cout << "mother decay vertex: " << mccand->theMother()->pdtEntry()->name() << " " 
// 		   << prodvtx << endl ;
	      if( (prodvtx-tagBvtx).mag()<0.001 ) { mcmatch=2 ; ++ngood ; }
	    } else {
	      //cout << "strange: particle has no mother?" << endl ;
	    }
	  }
	  
	  TagTrkAbsFitWrapper* thisabsfit = new TagTrkAbsFitWrapper(thiscand,mcmatch) ;
	  absfitlist.push_back(thisabsfit) ;
	
	  VtkPairPoca* thispoca = new VtkPairPoca(tagBabsfit,*thisabsfit) ;
	  pocalist.push_back(thispoca) ;
	  if( thispoca->success()  && thispoca->chisq()< gMaxSeedChisq )
	    thispoca->calcVertex() ;
	  
	}
      
      std::sort(pocalist.begin(),pocalist.end(),sortQuality) ;

      // fill my ntuple
      //       BtaBooster abooster(bc) ;
      //       int iseed=0 ;
      //       for( std::vector<VtkPairPoca*>::iterator pocait = pocalist.begin() ;
      // 	   pocait != pocalist.end() ; ++pocait, ++iseed) {
      // 	VtkPairPoca* thispoca = *pocait ;
      // 	const BtaCandidate* thiscand = thispoca->second()->cand() ;
      
      // 	HepPoint thisdecayvtx(thispoca->pos1()(1),thispoca->pos1()(2),thispoca->pos1()(3)) ;
      // 	HepPoint productionvtx(par(1),par(2),par(3)) ;
      // 	Hep3Vector momvec(par(4),par(5),par(6)) ;
      // 	double ctau = momvec.dot(thisdecayvtx-productionvtx)/momvec.mag() ;
      // 	double determinant = thispoca->vertexCov().determinant() ;
      // 	double deltaZ = thispoca->vertexPos()(3)-zCP ;
      // 	double pt    = thiscand->pt() ;
      // 	BtaCandidate boosted(abooster.boostTo(*thiscand)) ;
      // 	double pstar = boosted.p() ;
      // 	int mcmatch = thispoca->second()->mcmatch() ;
      // #ifdef VTK_MONITOR
      // 	static TNtuple* nt = new TNtuple("tagnt","",
      // 					 "mcmatch:succ:chisq:detV:docaerr:quality:ctau:deltaZ:pt:pstar:iseed:nseed") ;
      // 	nt->Fill(mcmatch,thispoca->success(),thispoca->chisq(),determinant,
      // 		 thispoca->docaErr(),thispoca->quality(),ctau,deltaZ,pt,pstar,iseed,pocalist.size()) ;
      // #endif
      //       }

      // now we have a sorted list of seeds. start to create vertices
      // and add the tracks.

      std::vector<TagFastVertex*> tagVertexList ;
      clustererB(absfitlist,pocalist,zCP,tagVertexList) ;
      
      // select the best one
      double bestquality=0 ;
      _fastvertex=0 ;
      for( std::vector<TagFastVertex*>::iterator it = tagVertexList.begin() ;
	   it != tagVertexList.end(); ++it) {
	double quality = (*it)->quality() ;
	if( !_fastvertex || bestquality<quality ) {
	  _fastvertex = *it ;
	  bestquality = quality ;
	}
      }
      
      // find also which vertex is really the best (ie closest to mc truth)
      TagFastVertex* mcbestvertex(0) ;
      double bestmcdz=100 ;
      for( std::vector<TagFastVertex*>::iterator it = tagVertexList.begin() ;
	   it != tagVertexList.end(); ++it) {
	double mcdz = fabs( (*it)->par()(3) = tagBvtx.z()) ;
	if( !mcbestvertex || bestmcdz>mcdz ) {
	  bestmcdz = mcdz ;
	  mcbestvertex = *it ;
	}
      }

      if(_verbosity>=1) {
	static int ntot(0),nbest(0) ;
	++ntot ;
	if( mcbestvertex==_fastvertex) ++nbest ;
 	cout << "Ntot, Nbest: " << ntot << " " << nbest << endl ;
      }

      // dump the result
#ifdef VTK_MONITOR
      if( !tagVertexList.empty() ) {
	int iseed = 0 ;
	for( std::vector<TagFastVertex*>::iterator it = tagVertexList.begin() ;
	     it != tagVertexList.end(); ++it, ++iseed) {
	  HepVector    pos = (*it)->par() ;
	  HepSymMatrix cov = (*it)->cov() ;
	  double dx = pos(1) - tagBvtx.x() ;
	  double sx = sqrt(cov.fast(1,1)) ;
	  double dy = pos(2) - tagBvtx.y() ;
	  double sy = sqrt(cov.fast(2,2)) ;
	  double dz = pos(3) - tagBvtx.z() ;
	  double sz = sqrt(cov.fast(3,3)) ; 
	  double chisqprob = ChisqConsistency((*it)->chiSquare(),
					      (*it)->nDof()).significanceLevel() ;
	  
	  // one more try on ctau
	  HepPoint decayvtx(pos(1),pos(2),pos(3)) ;
	  HepPoint productionvtx(par(1),par(2),par(3)) ;
	  Hep3Vector momvec(par(4),par(5),par(6)) ;
	  double ctau = momvec.dot(decayvtx-productionvtx)/momvec.mag() ;
	  bool isbest   = (*it == _fastvertex) ;
	  bool ismcbest = (*it == mcbestvertex) ;
	  int ntrk = singletrackflag? -1 : (*it)->size() ;
	  static TNtuple* seednt = 
	    new TNtuple("seednt","","iseed:dz:sz:chisq:ntrk:ngood:mcmatch:nmc:isbest:nseed:ismcbest") ;
	  seednt->Fill(iseed, //dx,sx,dy,sy,
		       dz,sz,chisqprob,ntrk,ngood,
		       (*it)->mcmatch(),
		       (*it)->ngood(),float(isbest),tagVertexList.size(),float(ismcbest)) ;
	}
      }
#endif
      
      // now, build a new upsilon, including the tagging B

      if( _fastvertex ) {
	if(_verbosity>=2)
	  cout << "NOW ADDING TAG VERTEX TO UPSILON! " << endl
	       << "chisq/ndof ups  before: " << _vertex->chiSquare()/_vertex->nDof() << endl 
	       << "chisq/ndof tagB before: " << _fastvertex->chiSquare()/_fastvertex->nDof() 
	       << "   ntrk: " << _fastvertex->size() << endl ;
	//_fastvertex->vertex()->print() ;
	//	BtaCandidate* tagCandidate = _fastvertex->createTagCandidate() ;	
	// 	_vertex->getUpsilon()->replaceTagB(tagCandidate) ;
	//  	_vertex->updateIndex() ;
	//  	_vertex->fit() ;
	//  	_vertex->print() ;
	// 	assert(0) ;
	_fastvertex->vertex()->fitOneStep() ; // <--- update kinematic constraint ...
	//_fastvertex->vertex()->print() ;
	_vertex->replaceTagB(*(_fastvertex->vertex()),true) ;
	// remove the  vertex from the list (we don't want to delete it)
	std::vector<TagFastVertex*>::iterator it = 
	  std::find(tagVertexList.begin(),tagVertexList.end(),_fastvertex) ;
	if(it != tagVertexList.end()) tagVertexList.erase(it) ;

	if(_verbosity>=2)
	  cout << "chisq/ndof ups  after: " << _vertex->chiSquare()/_vertex->nDof() << endl ;
	//_vertex->print() ;

 	BbrDoubleErr deltaT = _vertex->deltaT() ;
 	BbrDoubleErr deltaZ = _vertex->deltaZ() ;
#ifdef VTK_MONITOR
 	static TNtuple* dtnt = new TNtuple("dtnt","","dt:sdt:dz:sdz:status:chisq:ndof:ntrk:mcmatch") ;
	dtnt->Fill(deltaT.value()-truthTool.deltaT(),
		   sqrt(deltaT.covariance()),
		   deltaZ.value()-truthTool.deltaZ(),
		   sqrt(deltaZ.covariance()),
		   _vertex->status(),
		   _vertex->chiSquare(),
		   _vertex->nDof(),
		   _fastvertex->size(),_fastvertex->mcmatch()) ;
#endif
	BtaCandidate fittedB = _vertex->getFitted(bc) ;

	BtaDeltaTConverter dzdt(deltaZ);
	BbrDoubleErr deltaTtauB = dzdt.deltaT(&fittedB) ;

	std::vector<int> posindexvec ;
	const Upsilon* thisupsilon = _vertex->getUpsilon() ;
	posindexvec.push_back(thisupsilon->posIndex()) ;
	posindexvec.push_back(thisupsilon->recB()->posIndex()) ;
	posindexvec.push_back(thisupsilon->tagB()->posIndex()) ;
	HepSymMatrix zcov = _vertex->cov(posindexvec) ;
#ifdef VTK_MONITOR
	static TNtuple* zcovnt = new TNtuple("zcovnt","","ntrk:status:szIP:szCP:szTAG") ;
	zcovnt->Fill( _fastvertex->size(),_vertex->status(),
		      sqrt(zcov(1,1)),sqrt(zcov(2,2)),sqrt(zcov(3,3)) ) ;
#endif

	// finally, compare to the default tagger

	// first get the correct upsilon
	
	HepAList<BtaCandidate>* ups4slist = 
	  Ifd< HepAList<BtaCandidate> >::get( const_cast<AbsEvent*>(absevent), IfdStrKey("Y4SPlaintag")) ;
	if( ups4slist ) {
	  HepAListIterator<BtaCandidate> ups4siter(*ups4slist) ;
	  BtaCandidate* y4s;
	  bool done= false ;
	  while( (y4s=ups4siter()) && !done) {
	    // find the Brec cand [note : the first one is assumed to be the rec one]
	    HepAListIterator<BtaCandidate> iterDau(y4s->daughterIterator());
	    BtaCandidate* brecCand = iterDau();
	    BtaCandidate* tagCand = iterDau();
	    // if the tagCand vertex failed we will not do the QA on this event
	    if( (done = (tagCand->decayVtx() && brecCand->isCloneOf(bc)) ) ) {
	      // get the BtaAttributes
	      BtaAttributes* attributes = 
		Ifd<BtaAttributes>::get(const_cast<AbsEvent*>(absevent),IfdStrKey("Default")) ;
	      AbsEventTag * tagInfo = attributes->find(y4s);
	      float deltaz,deltazCov ;
	      tagInfo->getFloat(deltaz,"DeltaZ");
	      tagInfo->getFloat(deltazCov,"CovDeltaZ");
	      
#ifdef VTK_MONITOR
	      static TNtuple* dzcompnt = 
		new TNtuple("dzcompnt","","dztrue:dz1:sdz1:dz2:sdz2:status:chisq:ntrk:mcmatch:singletrack") ;
	      dzcompnt->Fill(truthTool.deltaZ(),
			     deltaZ.value(),sqrt(deltaZ.covariance()),
			     deltaz,sqrt(deltazCov),
			     _vertex->status(),_vertex->chiSquare(),
			     _fastvertex->size(),_fastvertex->mcmatch(),singletrackflag) ;
#endif

	      BtaDeltaTConverter dzdt(deltaz,deltazCov);
	      BbrDoubleErr deltaTdefault = dzdt.deltaT(&bc) ;

#ifdef VTK_MONITOR	      
	      static TNtuple* dtcompnt = 
		new TNtuple("dtcompnt","","dttrue:dt1:sdt1:dt2:sdt2:dt3:sdt3:status:chisq:ntrk1:ntrk2:mcmatch:singletrack") ;
	      dtcompnt->Fill(1.e12*truthTool.deltaT(),
			     1.e12*deltaT.value(),1.e12*sqrt(deltaT.covariance()),
			     deltaTdefault.value(),sqrt(deltaTdefault.covariance()),
			     deltaTtauB.value(),sqrt(deltaTtauB.covariance()),
			     _vertex->status(),_vertex->chiSquare(),
			     _fastvertex->size(),
			     tagCand->nDaughters(),
			     _fastvertex->mcmatch(),singletrackflag) ;
#endif
	      if(_verbosity>=1)
		cout << "GeoKin size: " <<  tagCand->nDaughters() << " " << sqrt(deltazCov) << " "
		     << _fastvertex->size() << " "
		     << sqrt(deltaZ.covariance()) << endl ;
	    }
	  }

	  if(!done) cout << "couldn't find a matching upsilon in default list! " 
			 << ups4slist->length() << endl ;
	} else {
	  static int printit=10 ;
	  if(--printit>0)
	    cout << "couldn't find default upsilon list" << endl ;
	}
      }
          
#ifdef THISPARTISBUGGY
#endif

      // cleanup
      std::for_each(tagVertexList.begin(),tagVertexList.end(),babar::Collection::DeleteObject());
      tagVertexList.clear() ;      
      std::for_each(pocalist.begin(),pocalist.end(),babar::Collection::DeleteObject()) ;
      pocalist.clear() ;
      std::for_each(absfitlist.begin(),absfitlist.end(),babar::Collection::DeleteObject());
      absfitlist.clear() ;

    } 
    
  }



  void TagVertexer::clustererA(const std::vector<TagTrkAbsFitWrapper*>& absfitlist,
			       const std::vector<VtkPairPoca*>&  pocalist,
			       double zCP,
			       std::vector<TagFastVertex*>& vertexlist)
  {
    for( std::vector<VtkPairPoca*>::const_iterator pocait = pocalist.begin() ;
	 pocait != pocalist.end() ; ++pocait) {
      const TagTrkAbsFitWrapper* theTrack = (*pocait)->second() ;
      // try to add to existing vertex
      bool isused=false ;
#define ADD_TO_ALL
#ifdef ADD_TO_BEST
      // first find the seed where this vertex fits best
      TagFastVertex* bestvertex(0) ;
      double bestdchisq = 0 ;
      for( std::vector<TagFastVertex*>::iterator it = vertexlist.begin() ;
	   it != vertexlist.end(); ++it) {
	double dchisq = (*it)->deltaChiSq(*theTrack).chisq() ;
	if( !bestvertex || bestdchisq>dchisq ) {
	  bestvertex = *it ;
	  bestdchisq = dchisq ;
	}
      }
      
      if( bestvertex && bestdchisq<trkaddchisqmax ) {
	bestvertex->add(*theTrack) ;
	isused = true ;
	//cout << "dchisq: " << dchisq << " " << bestdchisq << endl ;
      } else if (bestvertex && bestvertex->mcmatch()==2 && theTrack->mcmatch()==2) {
	cout << "missing this good track: " << bestdchisq << " " 
	     << bestvertex->size() << endl ;
      }
#elseifdef ADD_TO_ALL
      // add this track to all vertices it fits
      for( std::vector<TagFastVertex*>::iterator it = vertexlist.begin() ;
	   it != vertexlist.end(); ++it) {
	double dchisq = (*it)->deltaChiSq(*theTrack).chisq() ;
	if( dchisq < trkaddchisqmax ) {
	  (*it)->add(*theTrack) ;
	  isused = true ;
	}
      }
#else
      // add this track to the first vertex it fits
      for( std::vector<TagFastVertex*>::iterator it = vertexlist.begin() ;
	   it != vertexlist.end() && !isused; ++it) {
	double dchisq = (*it)->deltaChiSq(*theTrack).chisq() ;
	if( dchisq < gMaxChisqTrackAdd ) {
	  (*it)->add(*theTrack) ;
	  isused = true ;
	}
      }
#endif
      
      // create a new seed vertex if we could not use this track
      if( !isused && vertexlist.size()< gMaxNVertices && (*pocait)->isFitted() &&
	  fabs((*pocait)->vertexPos()(3) -zCP) < gMaxSeedDeltaZ ) {
	TagFastVertex* vertex = new TagFastVertex(**pocait) ;
	if( vertex->status()==BtaAbsVertex::Success ) {
	  vertexlist.push_back(vertex) ;
	} else {
	  delete vertex ;
	}
      }
    }
  }
  
  
  void TagVertexer::clustererB(const std::vector<TagTrkAbsFitWrapper*>& tracklist,
			       const std::vector<VtkPairPoca*>&  pocalist,
			       double zCP,
			       std::vector<TagFastVertex*>& vertexlist)
  {
    //    cout << "TagFastVertexer::clustererB" << endl ;
    // copy the vectors into lists
    std::list<VtkPairPoca*>  thepocalist ;
    std::list<TagTrkAbsFitWrapper*> thetracklist ;

    for( std::vector<VtkPairPoca*>::const_iterator it = pocalist.begin() ;
	 it != pocalist.end() ; ++it )
      if( (*it)->isFitted() && 
	  (*it)->chisq() < gMaxSeedChisq &&
	  fabs((*it)->vertexPos()(3) -zCP) < gMaxSeedDeltaZ )
	thepocalist.push_back(*it) ;
    int nseeds = thepocalist.size() ;
    //thepocalist.insert(thepocalist.end(),pocalist.begin(),pocalist.end()) ;
    
    thetracklist.insert(thetracklist.end(),tracklist.begin(),tracklist.end()) ;

    while(!thepocalist.empty() && vertexlist.size() < gMaxNVertices ) {
      VtkPairPoca* thepoca = thepocalist.front() ;
      thepocalist.pop_front() ;
      if(thepoca) {

	// remove this track from the track list
	std::list<TagTrkAbsFitWrapper*>::iterator it=
	  remove_if(thetracklist.begin(),thetracklist.end(),
		    bind2nd(std::equal_to<TagTrkAbsFitWrapper*>(),thepoca->second())) ;
	thetracklist.erase(it,thetracklist.end()) ;

	// create a vertex
	TagFastVertex* vertex = new TagFastVertex(*thepoca) ;
	if( vertex->status()==BtaAbsVertex::Success ) {
	  vertexlist.push_back(vertex) ;


	  // collect candidate tracks
	  std::vector<TrackVertexMatch> trackcandidates ;
	  for(std::list<TagTrkAbsFitWrapper*>::const_iterator it= thetracklist.begin();
	      it!= thetracklist.end() ; ++it) {
	    TrackVertexMatch match = vertex->deltaChiSq(**it) ;
	    if( match.chisq() < gMaxChisqTrackAdd ) trackcandidates.push_back( match ) ;
	  }
	  
	  std::sort( trackcandidates.begin(), trackcandidates.end() ) ;
	  bool first = true ;	  
	  for( std::vector<TrackVertexMatch>::iterator imatch = trackcandidates.begin() ;
	       imatch != trackcandidates.end() ; ++imatch ) {
	    if(first || (vertex->deltaChiSq(*(imatch->track())).chisq() < 
		 gMaxChisqTrackAdd) ) {
	      
	      // add to vertex
	      vertex->add( *(imatch->track()) ) ;
	      first = false ;
	      // remove from track list ?
	      //it = find(thetracklist.begin(),thetracklist.end(),imatch->track()) ;
	      //tracklist.erase(it) ;
	      
	      // remove the seed that contains this track
	      for( std::list<VtkPairPoca*>::iterator it=thepocalist.begin() ;
		   it!=thepocalist.end(); ++it)
		if( (*it) && (*it)->second()==imatch->track() ) {
		  *it = 0;
		  break ;
		}
	    }
	  }
	} else {
	  delete vertex ;
	}
      }
    }

    if( _verbosity>=1 && 
	(vertexlist.size() > 1 || 
	 (vertexlist.size()==1 && vertexlist.front()->size()==1) ) ) {
      
      cout << "----------------------------" << endl ;
      cout << "number of tracks: " << tracklist.size() << " | " ;
      for(std::vector<TagTrkAbsFitWrapper*>::const_iterator it = tracklist.begin(); 
	  it != tracklist.end() ; ++it)
	cout << (*it)->mcmatch() << " " ;
      cout << endl ;
      
      cout << "number of vertex seeds: " << nseeds << endl ;
      
      cout << "number of vertices: " << vertexlist.size() << " | "  ;
      for(std::vector<TagFastVertex*>::const_iterator it = vertexlist.begin() ;
	  it != vertexlist.end(); ++it) {
	cout << (*it)->size() << "(" << (*it)->ngood() <<","
	     << setprecision(3) << (*it)->chisqprob() << ","
	     << setprecision(3) << sqrt((*it)->cov().fast(3,3)) << "," 
	     << setprecision(3) << (*it)->quality() << ") "  ;
      }
      cout << endl ;
    }
  }

  BbrDoubleErr TagVertexer::deltaZ() const {
    return _vertex ? _vertex->deltaZ() : BbrDoubleErr() ; 
  }
  

  BbrDoubleErr TagVertexer::deltaTtauB() const {
    BbrDoubleErr dZ = deltaZ() ;
    BtaDeltaTConverter dzdt(dZ.value(),dZ.covariance());
    return dzdt.deltaT(&_recoB) ;
  }

  BbrDoubleErr TagVertexer::deltaT() const {
    return _vertex ? _vertex->deltaT() : BbrDoubleErr() ; 
  }

  BtaCandidate* TagVertexer::createUpsilonCandidate() const 
  {
    // create the upsilon tree from scratch
    BtaInterface bi ;
    BtaCandidate* upscand = bi.createTree(*_vertex->getUpsilon()) ;
    static const PdtEntry* upspdt = Pdt::lookup(PdtLund::Upsilon_4S) ;
    upscand->setType(upspdt) ;
    // fill the parameters
    _vertex->updateTree(*upscand) ;
    return upscand ;
  }

  int TagVertexer::nTaggingTracks() const 
  {
    return _fastvertex ? _fastvertex->size() : 0 ; 
  }
}

