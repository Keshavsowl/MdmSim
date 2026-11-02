#ifndef MdmSimDecayPhysics_h
#define MdmSimDecayPhysics_h 1
#include "G4VPhysicsConstructor.hh"
#include "DecayChannel.hh"
#include <vector>

namespace MdmSim
{
	class DecayPhysics : public G4VPhysicsConstructor
	{
		public:
			DecayPhysics(G4int verbose = 1);

			DecayPhysics(const G4String &name);

			virtual ~DecayPhysics();

			virtual void ConstructParticle();

			virtual void ConstructProcess();

			void SetDecayChannels(
					const std::vector<DecayChannel> &channels)

			{
				fDecayChannels = channels;

			}

		private:
			std::vector<DecayChannel> fDecayChannels;
	};
}
#endif


