#ifndef _VTK_EXTERNALBTAPARTICLE_HH_
#define _VTK_EXTERNALBTAPARTICLE_HH_

#include "VtxTreeFitter/VtkParticleBase.hh"
#include <CLHEP/Matrix/Vector.h>
#include <CLHEP/Matrix/SymMatrix.h>

namespace vtxtreefit
{

  class BtaComposite : public ParticleBase
  {
  public:
    BtaComposite(const BtaCandidate* bc, const ParticleBase* mother) ;
    virtual ~BtaComposite() ;

    // the number of parameters
    virtual int dim() const { return _hasEnergy ? 8 : 7 ; }// (x,y,z,t,px,py,pz,(E))

    // the number of 'measurements'
    int dimM() const        { return _hasEnergy ? 7 : 6 ; }
    ErrCode projectBtaComposite(const FitParams&, Projection&) const ;
    virtual ErrCode projectConstraint(Constraint::Type, const FitParams&, Projection&) const ;
 
    virtual ErrCode initPar1(FitParams*) ; 
    virtual ErrCode initPar2(FitParams*) ; 
    virtual int type() const { return kBtaComposite ; }  
    
    virtual int posIndex() const { return index()   ; }
    virtual int tauIndex() const { return index()+3 ; }
    virtual int momIndex() const { return index()+4 ; }

    virtual bool hasEnergy() const { return _hasEnergy ; }
    virtual bool hasPosition() const { return true ; }

    virtual void updCache() ;
    virtual double chiSquare(const FitParams* fitparams) const ; 

    virtual void addToConstraintList(constraintlist& alist, int depth) const {
      alist.push_back( Constraint(this,Constraint::btacomposite,depth,dimM()) ) ;
      alist.push_back( Constraint(this,Constraint::geometric,depth,3) ) ;
    }
    
  protected: // I hate this, so we need to change the design ...
    // cache
    HepVector _m ;    // 'measurement' (x,y,zpx,py,pz,E)
    HepSymMatrix _matrixV ; // covariance in measurement
    bool _hasEnergy ;
  } ;

}

#endif
