#ifndef __VTKDECAYCHAIN_HH__
#define __VTKDECAYCHAIN_HH__

#include <map>
#include "VtxTreeFitter/VtkParticleBase.hh"
#include "VtxTreeFitter/VtkMergedConstraint.hh"

class BtaCandidate ;

namespace vtxtreefit {

  class FitParams ;
  class ParticleBase ;

  class DecayChain
  {
  public:
    DecayChain() : _mother(0) {}

    DecayChain(const BtaCandidate* bc, bool forceFitAll=false)  ;
    ~DecayChain() ;
    
    int dim() const { return _dim ; }

    void initConstraintList() ;
    ErrCode init(FitParams* par) ;
    ErrCode filter(FitParams* par, bool firstpass=true) ;
    double chiSquare(const FitParams* par) const ;

    ParticleBase* mother() { return _mother ; }
    const ParticleBase* cand() { return _cand ; }
    const ParticleBase* mother() const { return _mother ; }
    const ParticleBase* locate(const BtaCandidate* bc) const ;

    int index(const BtaCandidate* bc) const ;
    int posIndex(const BtaCandidate* bc) const ;
    int momIndex(const BtaCandidate* bc) const ;
    int tauIndex(const BtaCandidate* bc) const ;
    void setOwner(bool b) { _isOwner=b ;}
    int momIndex() const ;

    void printConstraints(std::ostream& os=std::cout) const ;

  private:
    int _dim ;
    ParticleBase* _mother ;     // head of decay tree
    const ParticleBase* _cand ; // fit candidate (not same to mother in case of bs/be constraint)
    ParticleBase::constraintlist _constraintlist ;
    std::vector<Constraint*> _mergedconstraintlist ;
    MergedConstraint mergedconstraint ;
    typedef std::map<const BtaCandidate*,const ParticleBase*> ParticleMap ;
    ParticleMap _particleMap ;
    bool _isOwner ;
  } ;

}



#endif
