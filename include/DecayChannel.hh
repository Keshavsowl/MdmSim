#ifndef MdmSimDecayChannel_h
#define MdmSimDecayChannel_h 1

#include "globals.hh"

#include<string>

namespace MdmSim
{
	struct DecayProduct
	{
		std::string particleType = "ion";
		G4int Z = 0;
		G4int A = 0;

		G4double charge= 0.;
		G4double excitationEnergy = 0.;
	};

	struct DecayChannel
	{
		G4int parentZ = 0;
		G4int parentA = 0;
		G4double parentExcitationEnergy =0.;

		DecayProduct product1;
		DecayProduct product2;

		G4double probability = 1.;
	};

}
#endif
