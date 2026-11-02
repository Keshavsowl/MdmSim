#ifndef MdmSimNaIDetectorSD_h
#define MdmSimNaIDetectorSD_h 1
#include "G4VSensitiveDetector.hh"
#include "SiDetectorHit.hh"

class G4Step;
class G4HCofThisEvent;

namespace MdmSim
{
	class NaIDetectorSD : public G4VSensitiveDetector
	{
		public:
			NaIDetectorSD(const G4String &name,
			              const G4String &hitsCollectionName);

			~NaIDetectorSD() override;

			void Initialize(G4HCofThisEvent *hce) override;
			G4bool ProcessHits(G4Step *step,
					G4TouchableHistory *history) override;

			void EndOfEvent(G4HCofThisEvent *hce) override;

			private:

			SiDetectorHitsCollection *fHitsCollection =nullptr;

	};

}
#endif
