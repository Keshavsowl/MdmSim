#include "TrackingAction.hh"
#include "Analysis.hh"
#include "TrackRole.hh"

#include "G4Track.hh"

namespace MdmSim
{
    void TrackingAction::PreUserTrackingAction(const G4Track *track)
    {
        const TrackRole *role = TrackRole::Of(*track);
        if (!role)
            return;

        Analysis *analysis = Analysis::Instance();
        switch (role->GetRole())
        {
        case TrackRole::kReactionLight:
            analysis->SetAlpha1TrackId(track->GetTrackID());
            break;
        case TrackRole::kPrimaryDecayLight:
            analysis->SetAlpha2TrackId(track->GetTrackID());
            break;
        case TrackRole::kDecayGamma:
            analysis->SetGammaTrackId(track->GetTrackID());
            break;
        default:
            break;
        }
    }
}
