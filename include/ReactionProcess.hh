#ifndef MdmSimReactionProcess_h
#define MdmSimReactionProcess_h

#include "G4VDiscreteProcess.hh"

#include <map>
#include <string>

namespace MdmSim
{
    class ReactionProcess : public G4VDiscreteProcess
    {
    public:
        ReactionProcess(const G4String &name = "Reaction");
        ~ReactionProcess();

        G4double GetMeanFreePath(const G4Track &, G4double,
                                 G4ForceCondition *);
        // Returns the exact distance to the sampled reaction point, so the
        // reaction depth is uniform and every selected beam ion reacts in the
        // target (the base class would treat it as an exponential mean free path).
        G4double PostStepGetPhysicalInteractionLength(const G4Track &, G4double,
                                                      G4ForceCondition *) override;
        G4VParticleChange *PostStepDoIt(const G4Track &, const G4Step &);

        void StartTracking(G4Track *);

        void SetParams(std::map<std::string, G4double>);

    private:
        void TwoBody(const G4Track &aTrack, const G4ThreeVector &reactionPosition);
        G4double TargetPathLength(const G4Track &aTrack) const;
        G4double TargetDepthAlongTrack(const G4Track &aTrack) const;

        G4int fBeamZ = 0;
        G4int fBeamA = 0;
        G4int fTargetZ = 0;
        G4int fTargetA = 0;
        G4int fLightProductZ = 0;
        G4int fLightProductA = 0;
        G4int fLightProductCharge = 0;
        G4double fLightProductEx = 0.;
        G4int fHeavyProductZ = 0;
        G4int fHeavyProductA = 0;
        G4int fHeavyProductCharge = 0;
        G4double fHeavyProductEx = 0.;
        G4double fReactionFraction = 0.;
        G4double fReactionProbability = 1.;
        G4bool fWillReact = true;
        // Light-product cone: when fLightConeHalfAngle > 0, the CM direction is
        // resampled until the lab light product lies within this half-angle of
        // the MDM axis. Analysis records the number of tries per reaction so the
        // efficiency can be normalised to isotropic CM emission.
        G4double fLightConeHalfAngle = 0.;
        G4double fMdmAngle = 0.;
    };
}
#endif
