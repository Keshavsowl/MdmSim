#ifndef MdmSimTrackRole_h
#define MdmSimTrackRole_h 1

#include "G4Track.hh"
#include "G4VUserTrackInformation.hh"

namespace MdmSim
{
    // Tags the particles of the reaction/decay chain so that detector hits can
    // be matched to the particle that made them, e.g. for
    // 14N(6Li,alpha1)16O*, 16O* -> alpha2 + 12C*, 12C* -> 12C + gamma:
    //   ReactionLight      = alpha1
    //   ReactionHeavy      = 16O*
    //   PrimaryDecayLight  = alpha2 (product 1 of the decay of ReactionHeavy)
    //   PrimaryDecayHeavy  = 12C*  (product 2 of the decay of ReactionHeavy)
    //   DecayGamma         = any gamma from a configured decay
    class TrackRole : public G4VUserTrackInformation
    {
    public:
        enum Role
        {
            kReactionLight,
            kReactionHeavy,
            kPrimaryDecayLight,
            kPrimaryDecayHeavy,
            kDecayGamma,
            kOtherDecayProduct
        };

        explicit TrackRole(Role role) : fRole(role) {}
        Role GetRole() const { return fRole; }

        static const TrackRole *Of(const G4Track &track)
        {
            return dynamic_cast<const TrackRole *>(track.GetUserInformation());
        }

    private:
        Role fRole;
    };
}

#endif
