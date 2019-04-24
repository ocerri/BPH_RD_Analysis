#ifndef __VTK_UPSILON_HH__
#define __VTK_UPSILON_HH__

#include <vector>
#include <CLHEP/Matrix/Vector.h>
#include <CLHEP/Matrix/SymMatrix.h>
#include "VtxTreeFitter/VtkInteractionPoint.hh"
class BbrDoubleErr ;

namespace vtxtreefit 
{

  class Upsilon : public InteractionPoint
  {
  public:
    Upsilon(const BtaCandidate* recoB, bool forceFitAll) ;
    Upsilon(const BtaCandidate* recoB, bool forceFitAll, bool addupsilon) ;

    virtual ~Upsilon() ;

    ErrCode initBeamEnergy(const BtaCandidate* recoB) ;
    virtual int dim() const { return 7 ; } // (x,y,z,px,py,pz,E)
    virtual ErrCode initPar1(FitParams*) ;
    virtual ErrCode initCov(FitParams*) const ;
    virtual double chiSquare(const FitParams*) const ;
    virtual int type() const { return kUpsilon ; }
     
    virtual int posIndex() const { return index() ; }
    virtual int tauIndex() const { return -1 ; }
    virtual int momIndex() const { return index()+3; }
    virtual bool hasEnergy() const { return true ; }
    virtual std::string name() const { return "Upsilon" ; }

    ErrCode projectBeamEnergyConstraint(const FitParams&, Projection&) const ;
    ErrCode projectLifetimeSumConstraint(const FitParams &, Projection &) const ;
    virtual ErrCode projectConstraint(Constraint::Type, const FitParams&, Projection&) const ;

    virtual void addToConstraintList(constraintlist& alist, int depth) const ;

    virtual ParticleBase* recB() const { return daughters()[0] ; }
    virtual ParticleBase* tagB() const { return daughters()[1] ; }
    void setLifetimeSumConstraint(bool value) {
      _constrainLifetimeSum = value ; }

    //void replaceTagB(BtaCandidate* newtagB) ;
  private:
    bool _constrainBeam ;
    bool _constrainLifetimeSum ;
    HepVector    _beamMom ;
    HepSymMatrix _beamCov ;
    HepSymMatrix _beamCovInv ;
    BtaCandidate* _tagBcand ;
  } ;
}

#endif
