#ifndef __VTK_CONSTRAINT_HH__
#define __VTK_CONSTRAINT_HH__

#include <string>
#include <iostream>
#include "VtxTreeFitter/VtkErrCode.hh"

namespace vtxtreefit
{
  class ParticleBase ;
  class Projection ;
  class FitParams ;

  class Constraint
  {
  public:
    // the order of these constraints is important: it is the order in
    // which they are applied.

    enum Type { unknown=0,
		beamspot,
		beamenergy,
		lifetime,
		btaresonance,
		btacomposite,
		track,
		photon,
		conversion,
		kinematic,
		massEnergy,
		geometric,
		mass,
		merged,
		ntypes} ;
    
    bool operator<(const Constraint& rhs) const ;
    
    bool operator==(const Constraint& rhs) const { 
      return _type == rhs._type ; }
    
    // accessors
    Type type() const { return _type ; }
    unsigned int dim() const { return _dim ; }
    bool isLineair() const { return _maxNIter <=1 ; }
    unsigned int nIter() const { return _maxNIter ; }

    Constraint() : _node(0),_depth(0),_type(unknown) {}

    Constraint( const ParticleBase* node, Type type, int depth, 
		unsigned int dim, unsigned int nhidden=0, 
		int maxniter=1, double precision=1e-5)
      : _node(node), _depth(depth), _type(type), _dim(dim), 
	_nHidden(nhidden), _weight(1), _maxNIter(maxniter) {}
    
    virtual ~Constraint() {}

    virtual ErrCode project(const FitParams& fitpar, Projection& p) const ;
    virtual ErrCode filter(FitParams* fitpar) const ;
    virtual void print(std::ostream& os=std::cout) const ;
    std::string name() const ;

    // set to minus one if constraints needs to be removed on next filter
    void setWeight(int w) { _weight = w<0 ? -1 : 1 ; }
    

  protected:
    Constraint(Constraint::Type type) :
      _node(0),_depth(0),_type(type),_dim(0),_nHidden(0),
      _weight(0),_maxNIter(0) {}
    void setDim(unsigned int d) { _dim = d ; }
    void setNIter(unsigned int d) { _maxNIter = d ; }
  private:
    const ParticleBase* _node ;
    int _depth ;
    Type _type ;
    unsigned int _dim ;
    // the number of hidden 'degrees of freedom'. always zero except for the 'photon' constraint
    unsigned int _nHidden ;
    // the weight: guassian constraint can be 'unfilter'
    int _weight ; 
    int _maxNIter ;     // maximum number of iterations for non-linear constraints
  } ;

}

#endif
