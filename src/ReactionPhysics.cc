#include "ReactionPhysics.hh"
#include "ReactionProcess.hh"

#include "G4GenericIon.hh"
#include "G4Alpha.hh"
#include "G4Deuteron.hh"
#include "G4He3.hh"
#include "G4IonTable.hh"
#include "G4ParticleTable.hh"
#include "G4Proton.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4Triton.hh"

#include "globals.hh"

#include "G4PhysicsConstructorFactory.hh"

#include <cmath>

namespace MdmSim
{
    G4_DECLARE_PHYSCONSTR_FACTORY(ReactionPhysics);

    namespace
    {
        G4ParticleDefinition *FindConfiguredIon(const std::map<std::string, G4double> &params)
        {
            const auto zIt = params.find("BeamZ");
            const auto aIt = params.find("BeamA");
            if (zIt == params.end() || aIt == params.end())
                return nullptr;

            const auto z = static_cast<G4int>(std::lround(zIt->second));
            const auto a = static_cast<G4int>(std::lround(aIt->second));
            auto *particleTable = G4ParticleTable::GetParticleTable();
            if (z == 1 && a == 1)
                return particleTable->FindParticle("proton");
            if (z == 1 && a == 2)
                return particleTable->FindParticle("deuteron");
            if (z == 1 && a == 3)
                return particleTable->FindParticle("triton");
            if (z == 2 && a == 3)
                return particleTable->FindParticle("He3");
            if (z == 2 && a == 4)
                return particleTable->FindParticle("alpha");
            if (z > 0 && a > 0)
                return G4IonTable::GetIonTable()->GetIon(z, a, 0.);
            return nullptr;
        }
    }

    ReactionPhysics::ReactionPhysics(G4int)
        : G4VPhysicsConstructor("ReactionPhysics")
    {
    }

    ReactionPhysics::ReactionPhysics(const G4String &name)
        : G4VPhysicsConstructor(name)
    {
    }

    ReactionPhysics::~ReactionPhysics()
    {
    }

    void ReactionPhysics::ConstructParticle()
    {
        G4GenericIon::GenericIon();
        G4Proton::Proton();
        G4Deuteron::Deuteron();
        G4Triton::Triton();
        G4He3::He3();
        G4Alpha::Alpha();
    }

    void ReactionPhysics::ConstructProcess()
    {
        auto attachReaction = [&](G4ParticleDefinition *particle) {
            if (!particle)
                return;
            G4ProcessManager *pm = particle->GetProcessManager();
            if (!pm)
                return;

            const G4ProcessVector *plist = pm->GetProcessList();
            for (size_t i = 0; i < plist->size(); ++i)
            {
                if ((*plist)[i]->GetProcessName() == "Reaction")
                    return;
            }

            auto *reactionProcess = new ReactionProcess("Reaction");
            reactionProcess->SetParams(fReactionParams);
            pm->AddDiscreteProcess(reactionProcess);
        };

        attachReaction(G4GenericIon::GenericIon());
        attachReaction(G4Proton::Proton());
        attachReaction(G4Deuteron::Deuteron());
        attachReaction(G4Triton::Triton());
        attachReaction(G4He3::He3());
        attachReaction(G4Alpha::Alpha());
        attachReaction(FindConfiguredIon(fReactionParams));
    }

}
