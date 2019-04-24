#include "BaBar/BaBar.hh"
#include "VtxTreeFitter/VtkHelixUtils.hh"
#include <iostream>
#include <CLHEP/Matrix/Vector.h>
#include <CLHEP/Matrix/Matrix.h>
#include <BField/BFieldFixed.hh>
using std::cout;
using std::endl;

using namespace vtxtreefit ;

int main()
{
  const double pi = 3.1415927 ;
  HepVector helixpar(6) ;

  helixpar[VtkHelixUtils::ex_d0]     = 2 ;
  helixpar[VtkHelixUtils::ex_phi0]   = +pi-0.1;
  helixpar[VtkHelixUtils::ex_omega]  = 0.05 ;
  helixpar[VtkHelixUtils::ex_z0]     = 0.5 ;
  helixpar[VtkHelixUtils::ex_tanDip] = -7 ;
  helixpar[VtkHelixUtils::ex_flt]    = 10 ;

  cout << "This goes in: " << endl ;
  VtkHelixUtils::printHelixPar(helixpar) ;

  HepVector vertexpar(6) ;
  int charge ;
  BFieldFixed fieldmap(0,0,1.5) ;
  VtkHelixUtils::vertexFromHelix(helixpar,fieldmap,vertexpar,charge) ;

  cout << "This convertes to: " << endl ;
  VtkHelixUtils::printVertexPar(vertexpar,charge) ;

  HepVector helixparback(6) ;
  HepMatrix jacobian(6,6) ;
  VtkHelixUtils::helixFromVertex(vertexpar,charge,fieldmap,helixparback,jacobian) ;

  


  cout << "We get back: " << endl ;
  VtkHelixUtils::printHelixPar(helixparback) ;
  
  // numeric check of the jacobian
  double delta = 1e-5 ; 
  HepVector vertexpartmp(6) ;
  HepVector helixpartmp(6) ;
  HepMatrix jacobiantmp(6,6) ;
  
  HepMatrix jacobiannum(6,6) ;
  VtkHelixUtils::helixFromVertexNumerical(vertexpar,charge,fieldmap,helixparback,jacobiannum) ;

  for(int jin=0; jin<6; ++jin) {
    cout << VtkHelixUtils::vertexParName(jin) << endl ;
    for(int i=0; i<6; ++i) 
      vertexpartmp[i] = vertexpar[i] ;
    vertexpartmp[jin] += delta ;
    //printvertexpar(vertexpartmp,charge) ;
    VtkHelixUtils::helixFromVertex(vertexpartmp,charge,fieldmap,helixpartmp,jacobiantmp) ;
    for(int iex=0; iex<6; ++iex) {
      double numderiv = (helixpartmp[iex]-helixparback[iex])/delta ;
      double anaderiv = jacobiantmp[iex][jin] ;
      double dumderivB = jacobiannum[iex][jin] ;
      cout << "    " 
	   << VtkHelixUtils::helixParName(iex) 
	   << " " << numderiv << " " << anaderiv << " " << dumderivB << endl ;
    }
  }
}
