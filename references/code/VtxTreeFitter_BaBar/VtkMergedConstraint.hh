#ifndef __VTK_MERGEDCONSTRAINT_HH__
#define __VTK_MERGEDCONSTRAINT_HH__

#include <vector>
#include "VtxTreeFitter/VtkConstraint.hh"

namespace vtxtreefit
{
  class MergedConstraint : public Constraint
  {
  public:
    typedef std::vector<Constraint*> constraintlist ;

    MergedConstraint() : Constraint(Constraint::merged) {}
    virtual ~MergedConstraint() {}

    MergedConstraint( const constraintlist& list ) :
      Constraint(Constraint::merged),_list(list) {
      int d(0) ;
      for(constraintlist::iterator it = _list.begin() ;
	  it != _list.end(); ++it) d += (*it)->dim() ;
      setDim(d) ;
    }

    virtual ErrCode project(const FitParams& fitpar, Projection& p) const ;
    
    void push_back(Constraint* c) { 
      _list.push_back(c) ; 
      setDim(dim()+c->dim()) ; 
      setNIter(std::max(nIter(),c->nIter())) ;
    }

    virtual void print(std::ostream& os=std::cout) const ;

  private:
    constraintlist _list ;
  } ;

}

#endif
