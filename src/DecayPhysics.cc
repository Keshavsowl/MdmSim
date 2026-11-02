#include "DecayPhysics.hh"
#include "DecayProcess.hh"
#include "G4GenericIon.hh"
#include "G4Gamma.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

namespace MdmSim
{

	DecayPhysics::DecayPhysics(G4int verbose): G4VPhysicsConstructor("DecayPhysics"){
		SetVerboseLevel(verbose);
	}

	DecayPhysics::DecayPhysics(const G4String &name): G4VPhysicsConstructor(name){
	}
	
	DecayPhysics::~DecayPhysics() = default;

	void DecayPhysics::ConstructParticle()
	{
		G4GenericIon::GenericIonDefinition();
		G4Gamma::GammaDefinition();
	}

	void DecayPhysics::ConstructProcess()
	{
		auto *decayProcess = new DecayProcess();
		decayProcess->SetDecayChannels(fDecayChannels);
		auto *genericIon = G4GenericIon::GenericIonDefinition();
		auto *processManager = genericIon->GetProcessManager();

		if (processManager)
		{
			processManager->AddDiscreteProcess(decayProcess);
		 }
}
}
