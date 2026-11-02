#ifndef MdmSimDecayProcess_h

#define MdmSimDecayProcess_h

#include "G4VDiscreteProcess.hh"
#include "DecayChannel.hh"
#include <vector>

namespace MdmSim

{
	class DecayProcess : public G4VDiscreteProcess
	{
		public:
			DecayProcess(const G4String &name = "ConfiguredDecay");

			~DecayProcess();

			G4double GetMeanFreePath( const G4Track &, G4double, G4ForceCondition *);

			G4VParticleChange *PostStepDoIt( const G4Track &, const G4Step &);
			
			void StartTracking(G4Track *);


			void SetDecayChannels( const std::vector<DecayChannel> &channels)
			{
				fDecayChannels = channels;
			}



		private:
			const DecayChannel *FindMatchingChannel(
					const G4Track &aTrack) const;

			void TwoBodyDecay(
					const G4Track &aTrack,
					const DecayChannel &channel);

					std::vector<DecayChannel> fDecayChannels;
					const DecayChannel *fActiveChannel = nullptr;
					G4bool fWillDecay = false;

					};
		}
#endif

