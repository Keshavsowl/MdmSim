#include "NaIDetectorSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

namespace MdmSim
{
	NaIDetectorSD::NaIDetectorSD(const G4String &name, const G4String &hitsCollectionName) : G4VSensitiveDetector(name){
		collectionName.insert(hitsCollectionName);
	}

	NaIDetectorSD::~NaIDetectorSD()
	{
	}
	void NaIDetectorSD::Initialize(G4HCofThisEvent *hce)
	{
		fHitsCollection = new SiDetectorHitsCollection( SensitiveDetectorName, collectionName[0]);

		G4int hcID = G4SDManager::GetSDMpointer() -> GetCollectionID(collectionName[0]);
		hce->AddHitsCollection(hcID, fHitsCollection);
	}

	G4bool NaIDetectorSD::ProcessHits(G4Step *aStep, G4TouchableHistory *)
	{
		//Energy deposited in the NaI cystel
		G4double edep = aStep->GetTotalEnergyDeposit();

		if (edep ==0.)
			return false;

		// Create a new hit
		SiDetectorHit *newHit = new SiDetectorHit();

		//Get the parent aluminum housing copy number
		G4int copyNo = aStep->GetPreStepPoint() -> GetTouchableHandle() -> GetCopyNumber(1);

		newHit-> SetCopyNo(copyNo);

		// Track identification
		newHit -> SetTrackID( aStep->GetTrack() -> GetTrackID());

		// Energy Deposition
		newHit -> SetEnergyDeposit(edep);

		//Hit Time
		newHit -> SetTime(aStep->GetTrack()->GetGlobalTime());

		//Store the Hit
		fHitsCollection -> insert(newHit);

		return true;
	}

	void NaIDetectorSD::EndOfEvent(G4HCofThisEvent *)
	{
		return;
	}
}
