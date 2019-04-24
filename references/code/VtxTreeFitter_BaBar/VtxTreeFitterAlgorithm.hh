#ifndef __VTX_TREEFITTERALGORITHM_HH__
#define __VTX_TREEFITTERALGORITHM_HH__

#include <VtxBase/VtxAbsAlgorithm.hh>
#include <Beta/BtaFitParams.hh>

// this class is garbage, but that is because there is no definition.

namespace vtxtreefit 
{
  class Fitter ;
}

class VtxTreeFitterAlgorithm : public VtxAbsAlgorithm
{
public:
   // constructor
  VtxTreeFitterAlgorithm();

  // destructor
  virtual ~VtxTreeFitterAlgorithm();
  
  virtual BtaAbsVertex* compute( const BtaCandidate* decayTree ) ;
  
  // clone
  virtual VtxAbsAlgorithm* clone() const ;

  // momentum of the decaying state
  virtual HepLorentzVector    p4() const ;
  virtual BbrError         p4Err() ;
  virtual double chi2Contribution(const BtaCandidate&bc) const ;
  virtual BtaFitParams fitParams(const BtaCandidate&bc) const ;
  virtual BbrDoubleErr decayLength(const BtaCandidate& bc) const ;
  virtual BbrDoubleErr lifeTime(const BtaCandidate& bc) const ;

  virtual BtaCandidate getFitted(const BtaCandidate& c) const ;
  virtual BtaCandidate getFittedTree() const ;

  virtual BtaCandidate getFitted(const BtaCandidate& c, const BtaCandidate*) const {
    return getFitted(c) ; }
  virtual BtaCandidate getFittedTree( const BtaCandidate*) const ;
  virtual const BtaCandidate* fittedCand( const BtaCandidate& theCand,
					  BtaCandidate* theHeadOfTree ) const ;

  
  // give the client access to the underlying fitter object
  virtual const vtxtreefit::Fitter* vtxTreeFitter() const ;
  
private:
  vtxtreefit::Fitter* _vertex ;
  BtaFitParams _btaFitParams ;
} ;

#endif
