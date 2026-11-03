#ifndef MdmSimTrackingAction_h
#define MdmSimTrackingAction_h 1

#include "G4UserTrackingAction.hh"

namespace MdmSim
{
    /// Records the track IDs of the tagged reaction/decay particles (see
    /// TrackRole) so EventAction can match detector hits to them.
    class TrackingAction : public G4UserTrackingAction
    {
    public:
        TrackingAction() = default;
        ~TrackingAction() override = default;

        void PreUserTrackingAction(const G4Track *track) override;
    };
}

#endif
