#ifndef __VTK_HELIXUTILS_HH__
#define __VTK_HELIXUTILS_HH__

#include <TrkBase/TrkExchangePar.hh>
#include <string>

class HepVector ;
class HepMatrix ;
class BField ;
class HepPoint ;

namespace vtxtreefit
{
  
  class VtkHelixUtils
  {
  public:
    enum VertexCoor {in_x=0,in_y,in_z,in_px,in_py,in_pz} ;
    enum HelixCoor  {ex_d0=TrkExchangePar::ex_d0,
		     ex_phi0=TrkExchangePar::ex_phi0, 
		     ex_omega=TrkExchangePar::ex_omega, 
		     ex_z0=TrkExchangePar::ex_z0, 
		     ex_tanDip=TrkExchangePar::ex_tanDip,
		     ex_flt=5} ;
    
    static void vertexFromHelix(const HepVector& helixpar,
				const BField& fieldmap,
				HepVector& vertexpar, int& charge) ;
    
    static void helixFromVertex(const HepVector& vertexpar, int charge,
				const BField& fieldmap,
				HepVector& helixpar, HepMatrix& jacobian) ;

    static void helixFromVertexNumerical(const HepVector& vertexpar, int charge,
					 const BField& fieldmap,
					 HepVector& helixpar, HepMatrix& jacobian) ;
    
    static std::string helixParName(int i) ;
    static std::string vertexParName(int i) ;
    static void printHelixPar(const HepVector& helixpar) ;
    static void printVertexPar(const HepVector& vertexpar, int charge) ;

    static double helixPoca(const HepVector& helixpar1,
			    const HepVector& helixpar2,
			    double& flt1, double& flt2, 
			    HepPoint& v, bool parallel=false) ; 
    static double helixPoca(const HepVector& helixpar,const HepPoint& point,
			    double& flt) ;
    static double phidomain(const double phi) ;
  } ;

}


#endif

