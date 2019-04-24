#ifndef __VTK_INTERNALPARTICLE_HH__
#define __VTK_INTERNALPARTICLE_HH__

#include "VtxTreeFitter/VtkParticleBase.hh"
#include <vector>

namespace vtxtreefit 
{

  class InternalParticle : public ParticleBase
  {
  public:
    InternalParticle(const BtaCandidate* bc, const ParticleBase* mother, 
		     bool forceFitAll) ;
    virtual ~InternalParticle() ;

    virtual int dim() const { return mother() ? 8 : 7 ; }
    virtual void updateIndex(int& offset) ;

    virtual ErrCode initPar1(FitParams*) ; 
    virtual ErrCode initPar2(FitParams*) ; 
    virtual ErrCode initCov(FitParams*) const ; 
    virtual int type() const { return kInternalParticle ; }
    virtual void print(const FitParams*) const ; 
    virtual const ParticleBase* locate(const BtaCandidate* bc) const ;

    // parameter definition
    virtual int posIndex() const { return index()   ; }
    virtual int tauIndex() const { return mother() ? index()+3 : -1 ; }
    virtual int momIndex() const { return mother() ? index()+4 : index() + 3 ; }
    virtual bool hasEnergy() const { return true ; }
    virtual bool hasPosition() const { return true ; }
    virtual std::string parname(int index) const ;

    // constraints
    ErrCode projectKineConstraint(const FitParams&, Projection&) const ;
    ErrCode projectLifeTimeConstraint(const FitParams&, Projection&) const ;
    ErrCode projectConversionConstraint(const FitParams&,Projection& p) const ;
    virtual ErrCode projectConstraint(const Constraint::Type type,
				      const FitParams& fitparams,Projection& p) const ;
    virtual void forceP4Sum(FitParams&) const ; 

    // some of that other stuff
    virtual double chiSquare(const FitParams*) const ;

    virtual ParticleBase* addDaughter(const BtaCandidate*, bool forceFitAll=false) ;
    virtual void removeDaughter(const ParticleBase* pb) ;
    daucontainer& daughters() { return _daughters ; }
    const daucontainer& daughters() const { return _daughters ; }
    virtual const_iterator begin() const  { return _daughters.begin() ; }
    virtual const_iterator end()   const  { return _daughters.end() ; }

    virtual int nFinalChargedCandidates() const ;

    virtual void retrieveIndexMap(indexmap& anindexmap) const ;
    virtual void addToConstraintList(constraintlist& alist, int depth) const ;

    bool swapMotherDaughter(FitParams* fitparams, const ParticleBase* newmother) ;
    void setMassConstraint(bool b) { _massconstraint = b ; }

 protected:
    ErrCode initMom( FitParams* fitparams ) const ;
    virtual void addToDaughterList(daucontainer& list) ;
    typedef daucontainer::iterator iterator ;
    iterator begin() { return _daughters.begin() ; }
    iterator end()   { return _daughters.end() ; }
  private:
    daucontainer _daughters ;
    bool _massconstraint ;
    bool _lifetimeconstraint ;
    bool _isconversion ;
  } ;

}


#endif
