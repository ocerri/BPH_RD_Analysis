#ifndef __VTX_PARTICLEBASE_HH__
#define __VTX_PARTICLEBASE_HH__

class BtaCandidate ;
class BtaMcAssoc ;
#include <string>
#include <vector> 
#include "VtxTreeFitter/VtkConstraint.hh"
#include "VtxTreeFitter/VtkProjection.hh"

namespace vtxtreefit
{
  class FitParams ;

  class ParticleBase
  {
  public:
    enum ParticleType {kUpsilon=1,kInteractionPoint,
		       kBtaComposite,kBtaResonance,
		       kInternalParticle,kRecoTrack,
		       kResonance,kRecoPhoton,kMissingParticle} ;

		 
    ParticleBase(const BtaCandidate* bc, const ParticleBase* mother) ;
    virtual ~ParticleBase() ;

    static ParticleBase* createParticle(const BtaCandidate* bc, 
					const ParticleBase* mother,
					bool forceFitAll=false) ;

    virtual int dim() const = 0 ;
    virtual void updateIndex(int& offset) ;
    virtual ErrCode initPar1(FitParams*) = 0 ; // init everything that does not need mother vtx
    virtual ErrCode initPar2(FitParams*) = 0 ; // everything else
    virtual ErrCode initCov(FitParams*) const  ;
    virtual std::string parname(int index) const ;
    virtual void print(const FitParams*) const ; 
    virtual const ParticleBase* locate(const BtaCandidate* bc) const ;
    virtual std::string name() const ;

    const BtaCandidate* bc() const { return _bc ; }
    const int index() const { return _index ; }
    const ParticleBase* mother() const { return _mother ; }
    
    virtual ErrCode projectGeoConstraint(const FitParams&, Projection&) const ;
    virtual ErrCode projectMassConstraint(const FitParams&, Projection&) const ;
    virtual ErrCode projectConstraint(Constraint::Type, const FitParams&, Projection&) const ;
    virtual void forceP4Sum(FitParams&) const {} ; // force p4 conservation all along tree

    // indices to fit parameters
    virtual int type() const = 0 ;
    virtual int posIndex() const { return -1 ; }
    virtual int tauIndex() const { return -1 ; }
    virtual int momIndex() const { return -1 ; }
    
    // does the particle have a 3-momentum or a 4-momentum ?
    virtual bool hasEnergy() const { return false ; }

    // does the particle have is own decay vertex ? (resonances and
    // recoparticles do not)
    virtual bool hasPosition() const { return false ; }

    int eneIndex() const { return hasEnergy() ? momIndex()+3 : -1 ; }
    
    // calculates the global chisquare (pretty useless)
    virtual double chiSquare(const FitParams*) const { return 0 ; }
    
    // access to particle PDT parameters
    double pdtMass() const { return _pdtMass ; }
    double pdtLifeTime() const { return _pdtLifeTime ; }
    double pdtTau() const { return _pdtMass >0 ? _pdtLifeTime/_pdtMass : 0 ; } 
    int charge() const { return _charge ; }

    // access to daughters
    typedef std::vector<ParticleBase*> daucontainer ;
    typedef daucontainer::const_iterator const_iterator ;
    virtual const_iterator begin() const { return const_iterator(); }
    virtual const_iterator end()   const { return const_iterator(); }
    virtual ParticleBase* addDaughter(const BtaCandidate*, bool forceFitAll=false) 
    {  return 0 ; }
    virtual void removeDaughter(const ParticleBase* pb) {}

    typedef std::vector< std::pair<const ParticleBase*,int> > indexmap ;
    virtual void retrieveIndexMap(indexmap& anindexmap) const ;
    void setMother(const ParticleBase* m) { _mother = m ; } 
    
    typedef std::vector<vtxtreefit::Constraint> constraintlist ;
    virtual void addToConstraintList(constraintlist& alist, int depth) const = 0 ;
    virtual void addToDaughterList(daucontainer& list) {} 
    virtual int nFinalChargedCandidates() const { return 0 ; }
    void setCandidate(const BtaCandidate* bc) { _bc = bc ; }
    bool initFromTruth(const BtaMcAssoc& truthmap, FitParams& fitparams) const ;
  protected:
    static double pdtLifeTime(const BtaCandidate* bc)  ;
    static bool isAResonance(const BtaCandidate* bc) ;
    static double bFieldOverC() ; // Bz/c
    ErrCode initTau(FitParams* par) const ; 
  protected:
    void setIndex(int i) { _index = i ; }
  private:
    const BtaCandidate* _bc ;
    const ParticleBase* _mother ;
    int _index ;
    double _pdtMass ;     // cached mass
    double _pdtLifeTime ; // cached lifetime
    int _charge ;      // charge
 } ;

}


#endif
