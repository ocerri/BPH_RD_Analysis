#include "BaBar/BaBar.hh"
#include <algorithm>
#include <iomanip>
#include <list>
#include <set>
#include <vector>
//#include <pair>

#include <BbrGeom/BbrDoubleErr.hh>
#include <Beta/BtaCandidate.hh>
#include <TrkBase/TrkCompTrk.hh>
#include <TrkBase/TrkDifTraj.hh>
#include <TrkBase/TrkPoca.hh>
#include <TrkBase/TrkRecoTrk.hh>
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

#include "VtxTreeFitter/VtkClusterer.hh"
#include "VtxTreeFitter/VtkFitter.hh"
using std::cout;
using std::endl;
using std::flush;
using std::setprecision;
using std::setw;

#if !defined(__SUNPRO_CC) && !defined(__SunOS_5_8)

namespace vtxtreefit
{

  const unsigned int Clusterer::gMaxNVertices = 4 ;
  const double Clusterer::gMaxChisqTrackAdd = 8 ; // 2 DOF
  const double Clusterer::gMaxSeedChisq = 25    ; // 1 DOF
  const unsigned int Clusterer::gMaxCompositeCharge = 99 ;

  class ClusPairPoca ;

  class ClusCandidate : public TrkAbsFitWrapper
  {
  public:
    typedef std::set<const ClusPairPoca*> PocaContainer ;
    ClusCandidate(const BtaCandidate* cand) : TrkAbsFitWrapper(cand) {}
    ClusCandidate(const TrkAbsFitWrapper* abs) : TrkAbsFitWrapper(*abs) {}

    typedef std::set<const ClusPairPoca*>::iterator iterator ;
    typedef std::set<const ClusPairPoca*>::const_iterator const_iterator ;
    
    iterator begin() { return _pairpocas.begin() ; }
    iterator end()   { return _pairpocas.end() ; }
    const_iterator begin() const { return _pairpocas.begin() ; }
    const_iterator end()   const { return _pairpocas.end() ; }

    void attach(const ClusPairPoca* poca) {
      _pairpocas.insert(poca) ; }
    void detach(const ClusPairPoca* poca) {
      _pairpocas.erase(poca) ; }
    
    bool isComposite() const { return cand()->isComposite()  ; }
    
    double quality() const { return 0 ; }
    //return cand()->recoTrk() &&cand()->recoTrk()->hits() ? cand()->recoTrk()->hits()->nHit() : 0 ;}
  private:
    PocaContainer _pairpocas ;
  } ;
  
  class ClusPairPoca
  {
  public:
    ClusPairPoca(ClusCandidate& cand1, ClusCandidate& cand2) ;
    ClusPairPoca(const ClusPairPoca& rhs) ;
    ~ClusPairPoca() ;
    double chisq() const { return _chisq ; }
    double doca() const { return _doca ; }
    double docaErr() const { return _doca/sqrt(_chisq) ; }
    int success() const { return _success ; }
    const HepVector& vertexPos() const { return  _vertexPos ; }
    const HepSymMatrix& vertexCov() const { return  _vertexCov ; }
    void print() const ;
    template<class T> static HepVector convert(const T& vec) ;
    void calcVertex() ;
    double quality() const { return _quality ; }
    bool isFitted() const { return _isFitted ; }
    ClusCandidate* first() const { return _cand[0] ; }
    ClusCandidate* second() const { return _cand[1] ; }
    ClusCandidate* cand(int i) const { return _cand[i] ; }

    const HepVector& pos1() const { return _pos[0] ; }
    double charge() const { return first()->charge() + second()->charge() ; }
    double costheta() const { return _costheta ; }
    
    bool overlaps( const std::vector<const ClusCandidate*>& vec ) const {
      return std::find(vec.begin(),vec.end(),_cand[0])!=vec.end() ||
	std::find(vec.begin(),vec.end(),_cand[1])!=vec.end() ; }

    bool overlaps( const ClusPairPoca& rhs) {
      return _cand[0]==rhs._cand[0] || _cand[1]==rhs._cand[0] ||
	_cand[0]==rhs._cand[1] || _cand[1]==rhs._cand[1] ; }
  private:
    ClusCandidate* _cand[2] ;
    double _flt[2] ;
    HepVector _pos[2] ;
    HepSymMatrix _cov[2] ;
    TrkPoca _poca ;
    double _costheta ;
    double _chisq ;
    double _doca ;
    double _vardoca[2] ;
    HepVector _vertexPos ;
    HepSymMatrix _vertexCov ;
    double _quality ;
    bool _isFitted ;
    bool _success ;
  } ;

  template<class T>
  class CompareQuality
  {
  public:
    CompareQuality() {}
    bool operator()(const T& x, const T& y) {
      return x.quality() > y.quality() ;
    }
    bool operator()(const T* x, const T* y) {
      return x->quality() > y->quality() ;
    }
  } ;

  template<class T, class C> 
  class Overlaps
  {
  public:
    Overlaps(const C& c) : _c(c) {}
    bool operator()(const T& t) { return t.overlaps(_c) ; }
    bool operator()(const T* t) { return t->overlaps(_c) ; }
    const C& _c ;
  } ;

  template<class T>
  HepVector ClusPairPoca::convert(const T& vec) 
  {
    HepVector rc(3) ;
    rc(1) = vec.x() ;
    rc(2) = vec.y() ;
    rc(3) = vec.z() ;
    return rc ;
  }

  ClusPairPoca::ClusPairPoca(ClusCandidate& cand1, 
			     ClusCandidate& cand2)
    : _poca(cand1.traj(),0.0,cand2.traj(),0.0,1e-5), _chisq(-1),
      _vertexPos(3),_vertexCov(3),_quality(-999),_isFitted(false), _success(false) 
  {
    _cand[0] = &cand1 ;
    _cand[1] = &cand2 ;
    for(int i=0; i<2; ++i) _cand[i]->attach(this) ;

    bool verbose = cand1.isComposite() || cand2.isComposite() ;

    if( (_success = _poca.status().success()) ) {
      _flt[0] = _poca.flt1() ;
      _flt[1] = _poca.flt2() ;

      // ensure positive lifetime for composites
      Hep3Vector dir[2] ;
      for(int i=0; i<2; ++i) {
	BbrPointErr position = _cand[i]->trkFit()->positionErr(_flt[i]) ;
	_cov[i] = position.covMatrix() ;
	_pos[i] = convert<HepPoint>(position) ;
	dir[i]  = (_cand[i]->trkFit()->momentum( _flt[i] )).unit() ;
	
	// ensure positive lifetime for composites
	if( _cand[i]->isComposite() ) {
	  BtaFitParams fitpars = _cand[i]->cand()->fitParams() ;
	  Hep3Vector deltaX = HepPoint(fitpars.pos()) - HepPoint(position) ;
	  double costheta = deltaX.dot(fitpars.p3()) ;
	  if( costheta < 0 )  _success = 0 ;
	}
      }
      if( _success ) {
	_costheta = dir[0].dot(dir[1]) ;
	HepVector residual = _pos[1] - _pos[0] ;
	_doca = sqrt( dot(residual,residual) ) ;
	HepVector& jacobian = residual ;
	jacobian /= _doca ;
	for(int i=0; i<2; ++i)
	  _vardoca[i] = _cov[i].similarity(jacobian) ;
	_chisq = _doca*_doca/(_vardoca[0]+_vardoca[1]) ;
	
	
	if(false && verbose) {
	  cout << "     cached position: " << _cand[0]->cand()->fitParams().pos() << " "
	     << _cand[1]->cand()->fitParams().pos() << endl ;
	  
	  cout << "     flt: " << _poca.flt1() << " " << _poca.flt2() << endl ; 
	  cout << "     doca/chisq: " << _doca << " " << sqrt(_vardoca[0]) << " " << sqrt(_vardoca[1]) << " "
	       << _chisq << endl ;
	  cout << _cand[0]->trkFit()->positionErr(_flt[0]) << _cand[1]->trkFit()->positionErr(_flt[1]) ;

	  // collect the tracks
	  HepAList<BtaCandidate> alist ;
	  alist.insert( const_cast<BtaCandidate*>(_cand[0]->cand()) ) ;
	  alist.insert( const_cast<BtaCandidate*>(_cand[1]->cand()) ) ;
	  BtaOpMakeTree comb;
	  HepAListIterator<BtaCandidate> iter(alist) ;
	  BtaCandidate* cand = comb.createFromList(iter) ;
	  Fitter newvertex(*cand) ;
	  newvertex.fit() ;
	  BtaCandidate fittedcand = newvertex.getFitted(*cand) ;
	  cout << "      fitted vertex: " << fittedcand.decayVtx()->chiSquared() << " "
	       << fittedcand.fitParams().pos() << endl ;
	  delete cand ;
	}
      }
      // _vertexPos = 0.5*(_pos[1]+_pos[0]) ;
      // this is VERY wrong:
      //_chisq = residual.determineChisq(Hep3Vector(0,0,0)) ;
    }
  }

  ClusPairPoca::ClusPairPoca(const ClusPairPoca& rhs) 
    : _poca(rhs._poca), _chisq(rhs._chisq),_vertexPos(rhs._vertexPos),
      _vertexCov(rhs._vertexCov),_quality(rhs._quality),_isFitted(rhs._isFitted)
  {
    assert(0) ;
    //*this = rhs ;
    //for(int i=0; i<2; ++i) _cand[i]->attach(this) ;
  }

  void
  ClusPairPoca::calcVertex()
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

  ClusPairPoca::~ClusPairPoca()
  {
    for(int i=0; i<2; ++i) _cand[i]->detach(this) ;
  }

  void ClusPairPoca::print() const {
    for(int i=1; i<=3; ++i)
      cout << setprecision(5) << (char)('x'+i-1) 
	   << setw(12) << _vertexPos(i) << " +- " << setw(12) << sqrt(_vertexCov(i,i)) << endl  ;
  }
  
  inline bool sortQuality(const ClusPairPoca* lhs, const ClusPairPoca* rhs)
  {
    return lhs->quality() > rhs->quality() ;
  }

  class TrackVertexMatch
  {
  public:
    TrackVertexMatch(const ClusCandidate* track)
      : _track(track),_chisq(999),_logV(999) {}
    void setChisq(double chisq, double logV) { _chisq = chisq ; _logV = logV ; }
    double chisq() const { return _chisq ; }
    double density() const { return -_chisq -_logV ; }
    bool operator<(const TrackVertexMatch& rhs) const {
      return density() > rhs.density() ; 
    }
    const ClusCandidate* track() const { return _track ; }
  private:
    const ClusCandidate* _track ;
    double _chisq ;
    double _logV ;
  } ;

  class FastVertex : public BtaOperatorBase
  {
  public:
    typedef std::vector<const ClusCandidate*> trackcontainer ;

    FastVertex() : _vertex(0) {}
    FastVertex(const ClusPairPoca& poca) ;
    ~FastVertex() { delete _vertex ; delete _seedcand ; }
  
    int status() const { return  _vertex->status() ; }
    TrackVertexMatch deltaChiSq(const ClusCandidate& cand) const ;
    double add(const ClusCandidate& cand) ;
    int size() const { return _tracks.size() ; }
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
    
    BtaCandidate* createCandidate() const ;
    Fitter* vertex() { return _vertex ; }
    const ClusPairPoca* seedpoca() const { return _seedpoca ; }
    double quality() const { 
      //return chisqprob()/sqrt( cov().fast(3,3) ) ; } 
      return chisqprob()/cov().fast(3,3) ; }  //         <-- slightly better (89%)
    //return 1./( cov().fast(3,3) * chiSquare() / nDof() ) ; } // <-- the worst (85%)
    //return chisqprob() ; }
    //return 1/cov().fast(3,3) ; }   <--- worse than the worst

    const trackcontainer& tracks() const { return _tracks ; }
  private:
    const ClusPairPoca* _seedpoca ;
    Fitter*      _vertex ;
    BtaCandidate* _seedcand ;
    std::vector<int>     _indexvec ;
    trackcontainer _tracks ;
  } ;

  FastVertex::FastVertex(const ClusPairPoca& poca)
    : _seedpoca(&poca),_vertex(0),_seedcand(0)
  {
    const ClusCandidate* trk1 = poca.first() ;
    const ClusCandidate* trk2  = poca.second() ;
    static BtaOpMakeTree comb;
    _seedcand = comb.create(*(trk1->cand()),*(trk2->cand())) ;
    //    _seedcand->setType(tagBseed->cand()->pdtEntry());
    _vertex = new Fitter(*_seedcand) ;
    _vertex->fit() ;
    //_vertex->print() ;
    int posindex = _vertex->posIndex(_seedcand) ;
    for(int i=0; i<3; ++i) _indexvec.push_back(posindex+i) ;
    _tracks.push_back(trk1) ;
    _tracks.push_back(trk2) ;
  }

 
  TrackVertexMatch FastVertex::deltaChiSq(const ClusCandidate& cand) const
  {
    TrackVertexMatch match(&cand) ;
    HepVector    vertexpos = _vertex->par(_indexvec) ;
    HepSymMatrix vertexcov = _vertex->cov(_indexvec) ;
    HepPoint vertexpt(vertexpos(1),vertexpos(2),vertexpos(3)) ;
    TrkPoca thePoca(cand.traj(),0,vertexpt,1e-5);
    if( thePoca.status().success() ) {
      BbrPointErr trkposerr = cand.trkFit()->positionErr(thePoca.flt1()) ;
      HepVector trkpos = ClusPairPoca::convert<HepPoint>(trkposerr) ;
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

  double FastVertex::add(const ClusCandidate& cand)
  {
    _tracks.push_back(&cand) ;
    return _vertex->add(*(cand.cand())) ;
  }
  
  BtaCandidate* FastVertex::createCandidate() const 
  {
    // collect the tracks
    HepAList<BtaCandidate> alist ;
    double charge(0) ;
    for( std::vector<const ClusCandidate*>::const_iterator it = _tracks.begin() ;
	 it != _tracks.end() ; ++it) {
      alist     += const_cast<BtaCandidate*>((*it)->cand()) ;
      charge += (*it)->charge() ;
    }
    static BtaOpMakeTree comb;
    HepAListIterator<BtaCandidate> iter(alist) ;
    BtaCandidate* cand = comb.createFromList(iter) ;
    Fitter newvertex(*cand) ;
    newvertex.fit() ;
    BtaCandidate fittedcand = newvertex.getFitted(*cand) ;
    iter.rewind() ;
    BtaCandidate* daughter ;
    while ( (daughter= iter() ) )
      if(daughter->isComposite()) 
	cout << "lifetime of daughter: " << newvertex.decayLength(*daughter).value()  
	     << " " << sqrt(newvertex.decayLength(*daughter).covariance()) << endl ;

    delete cand ;
    return new BtaCandidate(fittedcand) ;
  }

  Clusterer::~Clusterer() { 
//     if(_vertex) delete _vertex ; 
//     if(_fastvertex) delete _fastvertex ; 
  }

  void dumppocalist(const std::list<ClusPairPoca*>& pocalist)
  {
    for( std::list<ClusPairPoca*>::const_iterator it = pocalist.begin() ;
	 it != pocalist.end(); ++it) {
      cout << "a poca: " << *it << " " << (*it)->quality() << " "
	   << (*it)->chisq() << " " << (*it)->first()->isComposite() << " "
	   << (*it)->second()->isComposite() << endl ;
    }
  }


  Clusterer::Clusterer(const HepAList<BtaCandidate>& inputlist,
		       const TrkAbsFitWrapper* tagcand,
		       const EventInfo* eventinfo,
		       const AbsEvent* absevent,
		       bool excludeconversions,
		       bool excludeclones) 
  {
    // first make the list of absfitwrappers (our interface to the btacandidate)
    typedef std::list<ClusCandidate> ClusCandContainer ;
    typedef std::list<ClusPairPoca*> ClusPocaContainer ;

    ClusCandContainer absfitlist ;
    
    HepAListIterator<BtaCandidate> inputiter(inputlist) ;
    BtaCandidate* thiscand ;
    while( (thiscand = inputiter() )) 
      if( thiscand->recoTrk() || thiscand->isComposite() ||
	  !thiscand->recoCalo() ) 
	absfitlist.push_back(ClusCandidate(thiscand)) ;

    if( tagcand ) absfitlist.push_back( ClusCandidate(tagcand) ) ;
    
    cout << "number of starting tracks:" << absfitlist.size() << endl ;

    // now make all pairs. this is the most time consuming part. it
    // will define our seeding list.
    ClusPocaContainer  pocalist ;
    std::set<const ClusCandidate*> badabsfitlist ; // will contain conversion tracks and ghosttracks
    bool good1(false),good2(false) ;

    for(ClusCandContainer::iterator trk1 =  absfitlist.begin() ; 
	trk1 != absfitlist.end(); ++trk1) {
      good1 = badabsfitlist.find(&(*trk1))==badabsfitlist.end() ;
      if(good1) {
	std::list<ClusPairPoca*>  tmppocalist ;
	for(ClusCandContainer::iterator trk2 = absfitlist.begin() ; 
	     trk2 != trk1 && good1; ++trk2) {
	  good2 = badabsfitlist.find(&(*trk2))==badabsfitlist.end() ;
	  if( good2 ) {
	    
	    ClusPairPoca* thispoca = new ClusPairPoca(*trk1,*trk2) ;
	    if( thispoca->success()  && thispoca->chisq() < gMaxSeedChisq ) {
	      
	      // is this a possible conversion or a duplicated track?
	      if( 1 - thispoca->costheta() < 0.01 ) {
		cout << "********* parallel tracks: "
		     << thispoca->costheta() << " "  << thispoca->charge() << " "
		     << trk1->quality() << " "
		     << trk2->quality() << endl ;
		if( thispoca->charge()==0 && excludeconversions ) {
 		  // this is a conversion. exclude both.
 		  good1=good2=false ;
 		} else if( abs(int(thispoca->charge()))==2 && excludeclones ) {
 		  // exclude the worst track
		  
 		  if( trk1->quality() > trk2->quality() )
 		    good2 = false ;
		  else
 		    good1 = false ;
 		}
	      }
	      
	      if(!good1) badabsfitlist.insert( &(*trk1) ) ;
	      if(!good2) badabsfitlist.insert( &(*trk2) ) ;
	      
	      if(good1&&good2) {
		thispoca->calcVertex() ;
		tmppocalist.push_back(thispoca) ;
	      }
	      else             delete thispoca ;
	    } else             delete thispoca ;
	  }
	}
	if(good1) {
	  tmppocalist.sort() ;
	  tmppocalist.sort(CompareQuality<ClusPairPoca>()) ;
	  pocalist.merge(tmppocalist,CompareQuality<ClusPairPoca>()) ;
	} else {
	  std::for_each(tmppocalist.begin(),tmppocalist.end(),babar::Collection::DeleteObject());
	}
      }
    }
    
    cout << "number of pocas: " << pocalist.size() << endl ;
    //dumppocalist(pocalist) ;

    // from the pocalist, make a new list with 'vertices' that are disjoint.
    vector< FastVertex* > vertexlist ;
    for( ClusPocaContainer::iterator it = pocalist.begin(); 
	 it != pocalist.end(); ++it) {
      bool found = false ;
      for( vector< FastVertex* >::iterator ivec = vertexlist.begin() ;
	   !found && ivec != vertexlist.end() ; ++ivec ) 
	found =  (*it)->overlaps( *(*ivec)->seedpoca() ) ;
      if(!found) {
	FastVertex* vertex = new FastVertex( **it ) ;
	if( vertex->status()==BtaAbsVertex::Success ) 
	  vertexlist.push_back(vertex) ;
	else
	  delete vertex ;
      }
    }
    cout << "number of vertex seeds: " << vertexlist.size() << endl ;

    // create a new list with tracks
    std::set<ClusCandidate*> tracktoaddlist ;
    for(ClusCandContainer::iterator itrk =  absfitlist.begin() ; 
	itrk != absfitlist.end(); ++itrk)
      if( badabsfitlist.find(&(*itrk))==badabsfitlist.end() ) 
	tracktoaddlist.insert(&(*itrk)) ;
    
    // start creating composites
    vector< BtaCandidate* > compositelist ;
    int tracksincomposites(0) ;

    while( !pocalist.empty() ) {
      ClusPairPoca* thispoca = pocalist.front() ;
      pocalist.pop_front() ;

      // is this a valid poca?
      if( badabsfitlist.find( thispoca->first() )==badabsfitlist.end() &&
 	  badabsfitlist.find( thispoca->second() )==badabsfitlist.end() ) {

	FastVertex* vertex = new FastVertex( *thispoca ) ;
	cout << "vertex: " << thispoca << " "
	     << vertex->chiSquare() << " " 
	     << vertex->tracks().size() << endl ;

	if( vertex->status()==BtaAbsVertex::Success ) {
	  tracktoaddlist.erase( thispoca->first() ) ;
	  tracktoaddlist.erase( thispoca->second() ) ;
	
	  cout << "tracktoaddlist size: " << tracktoaddlist.size() << endl ;
	  
	  // try to add new tracks
	  std::vector<TrackVertexMatch> trackcandidates ;
          for(std::set<ClusCandidate*>::const_iterator it= tracktoaddlist.begin();
              it!= tracktoaddlist.end() ; ++it) {
            TrackVertexMatch match = vertex->deltaChiSq(**it) ;
	    cout << " candidate for adding: " << match.chisq() << endl ;
            if( match.chisq() < 2*gMaxChisqTrackAdd ) trackcandidates.push_back( match ) ;
          }
	  
	  std::sort( trackcandidates.begin(), trackcandidates.end() ) ;
	  for( std::vector<TrackVertexMatch>::iterator imatch = trackcandidates.begin() ;
               imatch != trackcandidates.end() ; ++imatch ) {
	    const ClusCandidate* track = imatch->track() ;
	    double dchisq = vertex->deltaChiSq(*track).chisq() ;
	    if(dchisq < gMaxChisqTrackAdd) {
	      vertex->add( *track ) ;
	      cout << "added track: " << vertex->chiSquare() << " " << dchisq << endl ;
	      tracktoaddlist.erase(const_cast<ClusCandidate*>(track)) ;
	      badabsfitlist.insert(track) ;
	    }
	  }

	  // remove all seeds which have any of these tracks
	  ClusPocaContainer::iterator last = 
	    std::remove_if(pocalist.begin(),pocalist.end(),
			   Overlaps<ClusPairPoca,FastVertex::trackcontainer>(vertex->tracks())) ;
	  pocalist.erase(last,pocalist.end()) ;
	  
	  //cout << "number of tracks in vertex: " << vertex->tracks().size() << endl ;

	  // make a composite
	  BtaCandidate* cand =vertex->createCandidate() ;
	  cout << "new candidate: " 
	       << cand->nDaughters() << " " 
	       << cand->mass() << " "
	       << cand->decayVtx()->chiSquared() << " "
	       << cand->decayVtx()->nDof() << " " << cand->charge() << " "
	       << cand->fitParams().pos() << " " 
	       << cand->fitParams().p3() << endl ;

	  for( FastVertex::trackcontainer::const_iterator itrk = vertex->tracks().begin() ;
	       itrk != vertex->tracks().end() ; ++itrk)
	    if (!(*itrk)->isComposite() ) ++tracksincomposites ;

	  compositelist.push_back(cand) ;
	  // // 	// create a new absfitter
	  
	  double rawcharge = cand->charge() ;
	  if( fabs(rawcharge) <= gMaxCompositeCharge ) {
	    
	    double realcharge = rawcharge > 0.5 ? 1 : rawcharge < -0.5 ? -1 : 0 ;
	    cand->setCharge(realcharge) ;
	    absfitlist.push_back( ClusCandidate( cand ) ) ;
	    ClusCandidate& cand = absfitlist.back() ;
	    
	    // create new pairpocas
	    std::list<ClusPairPoca*>  tmppocalist ;
	    for(std::set<ClusCandidate*>::iterator itrk = tracktoaddlist.begin() ; 
		itrk != tracktoaddlist.end(); ++itrk) {
	      ClusPairPoca* newpoca = new ClusPairPoca(cand,**itrk) ;
	      if( newpoca->success()  && newpoca->chisq() < gMaxSeedChisq ) {
		newpoca->calcVertex() ;
		tmppocalist.push_back( newpoca ) ;
	      } else {
		delete newpoca ;
	      }
	    }
	    tracktoaddlist.insert(&cand) ;
	    cout << "new pocas from composite: " << tmppocalist.size() << endl ;
	    tmppocalist.sort( CompareQuality<ClusPairPoca>() ) ;
	    pocalist.merge(tmppocalist,CompareQuality<ClusPairPoca>()) ;
	    //dumppocalist(pocalist) ;
	  }

	  delete vertex ;
	}
      }
      delete thispoca ;
      cout << "current size: " << pocalist.size() << endl ;
    }
    
    cout << "total number composites: " << compositelist.size() << endl ;
    cout << "number of tracks in composites: " << tracksincomposites << " "
	 << tracktoaddlist.size() << endl ;
    cout << "top level vertices: " << flush ;
    for(std::set<ClusCandidate*>::const_iterator it= tracktoaddlist.begin();
	it!= tracktoaddlist.end() ; ++it) 
      if( (*it)->isComposite() ) cout << (*it)->cand()->fitParams().pos() << endl ;
	 

    cout << "++++++++++++++++++++++++++++++++++++++++" << endl ;

// 	if( vertex->success() ) {

// 	  cout << "chisq: " <<  thispoca->chisq() << " " << vertex->chiSquare() << endl ;
// 	  tracktoaddlist

// 	// good. now collect all neighbours that these candidates have
// 	// in common.
// 	set<ClusCandidate*> candidates ;
// 	set<const ClusPairPoca*>  usedpocas ;
// 	candidates.insert( thispoca->first() ) ;
// 	candidates.insert( thispoca->second() ) ;

// 	for( ClusCandidate::iterator ipoca = thispoca->first()->begin() ;
// 	     ipoca!=thispoca->first()->end(); ++ipoca)
// 	  for( ClusCandidate::iterator jpoca = thispoca->second()->begin() ;
// 	       jpoca!=thispoca->second()->end(); ++jpoca) {
// 	    for(int i=0; i<2; ++i)
// 	      for(int j=0; j<2; ++j)
// 		if( (*ipoca)->cand(i) == (*jpoca)->cand(j) ) {
// 		  candidates.insert( (*ipoca)->cand(i) ) ;
// 		  usedpocas.insert( *ipoca ) ;
// 		  usedpocas.insert( *jpoca ) ;
// 		  for( ClusCandidate::iterator kpoca = (*ipoca)->cand(i)->begin() ;
// 		       kpoca!= (*ipoca)->cand(i)->end(); ++kpoca)
// 		    usedpocas.insert( *kpoca ) ;
// 		}
// 	  }
// 	cout << "Common candidates: " << candidates.size() << endl ;

// 	// delete all these pocas
// 	for(  ClusPocaContainer::iterator it = pocalist.begin() ;
// 	      it != pocalist.end() ; ++it) {
// 	  std::set<const ClusPairPoca*>::iterator jt = usedpocas.find(*it) ;
// 	  if( jt != usedpocas.end() ) {
// 	    pocalist.erase(it) ;
// 	    usedpocas.erase(jt) ;
// 	    delete *it ;
// 	    --it ;
// 	  }
// 	}

// // 	// put all cluscandidates on the 'bad' list
// // 	badabsfitlist.insert( candidates.begin(), candidates.end() ) ;
     
// 	// create a composite
// 	HepAList<BtaCandidate> btacandlist ;
// 	for(set<ClusCandidate*>::iterator it = candidates.begin() ;
// 	    it != candidates.end() ; ++it)
// 	  btacandlist.append( const_cast<BtaCandidate*>((*it)->cand() )) ;
// 	BtaOpMakeTree treemaker ;
// 	HepAListIterator<BtaCandidate> iter( btacandlist) ;
// 	BtaCandidate* composite = treemaker.createFromList( iter ) ;
// 	Vertex fitter(*composite) ;
// 	fitter.fit() ;
// 	BtaCandidate* fittedcomposite = new BtaCandidate(fitter.getFitted(*composite)) ;
// 	compositelist.push_back( fittedcomposite ) ;
// 	cout << "composite: " << fittedcomposite->decayVtx()->chiSquared() << " "
// 	     << fittedcomposite->decayVtx()->nDof() << endl ;
// 	delete composite ;

// // 	// create a new absfitter
// // 	if( fabs(fittedcomposite.charge())<=1 ) {
// // 	  absfitlist.push_back( ClusCandidate( *fittedcomposite ) ) ;
// // 	  ClusCandidate& cand = absfitlist.back() ;
	  
// // 	  // create new pairpocas
// // 	  std::list<ClusPairPoca*>  tmppocalist ;
// // 	  for(std::vector<ClusCandidate>::iterator itrk =  absfitlist.begin() ; 
// // 	      itrk != absfitlist.end(); ++itrk) 
// // 	    if( badabsfitlist.find(&(*itrk))==badabsfitlist.end() ) {
// // 	      ClusPairPoca* thispoca = new ClusPairPoca(cand,*itrk) ;
// // 	      if( thispoca->success()  && thispoca->chisq() < gMaxSeedChisq ) {
// // 		tmppocalist.push_back( thispoca ) ;
// // 	      } else {
// // 		delete thispoca ;
// // 	      }
// // 	    }
// // 	  tmppocalist.sort( CompareQuality<ClusPairPoca*>() ) ;
// // 	  pocalist.merge(tmppocalist,CompareQuality<ClusPairPoca*>()) ;
//       }
      
    std::for_each(vertexlist.begin(),vertexlist.end(),babar::Collection::DeleteObject());
    std::for_each(compositelist.begin(),compositelist.end(),babar::Collection::DeleteObject());
    std::for_each(pocalist.begin(),pocalist.end(),babar::Collection::DeleteObject());
  }




#ifdef MOEILIJK
   
    // the final pocalist is sorted. start building vertices.
    std::vector<FastVertex*> clusterlist ;
    
    while( !pocalist.empty() ) {
      ClusPairPoca& thispoca = pocalist.front() ;
      if( badabsfitlist.find( thispoca->first() )==badabsfitlist.end() &&
	  badabsfitlist.find( thispoca->second() )==badabsfitlist.end() ) {
	


#endif
}


#endif
