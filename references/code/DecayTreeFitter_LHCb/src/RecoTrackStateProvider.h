#ifndef DECAYTREEFITTER_RECOTRACKSTATEPROVIDER_H
#define DECAYTREEFITTER_RECOTRACKSTATEPROVIDER_H

#include "TrackInterfaces/ITrackStateProvider.h"
#include "TrackKernel/TrackTraj.h"

namespace DecayTreeFitter
{
  // wrapper around ITrackStateProviderTool
  class RecoTrackStateProvider
  {
  public:
    RecoTrackStateProvider( const ITrackStateProviderTool& tool,
			    double ztolerance,
			    bool usetraj=false)
      : m_tool(&tool), m_ztolerance(ztolerance) {}

    RecoTrackStateProvider(bool usetraj=false)
      : m_tool(0), m_ztolerance(0) {}

    void state( LHCb::State& state, const LHCb::Track& track, double z) {
      if( m_tool ) {
	if( m_usetraj ) {
	  m_tool->stateFromTrajectory( state, track, z ).ignore() ;
	} else {
	  m_tool->state( state, track, z, m_ztolerance).ignore() ;
	}
      } else {
	if( m_usetraj ) {
	  LHCb::TrackTraj traj(track) ;
	  state = traj.state(z) ;
	} else {
	  state = track.closestState( z ) ;
	}
      }
    }
    
  private:
    const ITrackStateProviderTool* m_tool ;
    double m_ztolerance ;
    bool m_usetraj ;
  } ;
}

#endif
