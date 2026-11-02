#include "SiDetectorSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ios.hh"

namespace MdmSim
{

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    SiDetectorSD::SiDetectorSD(const G4String &name,
                               const G4String &hitsCollectionName)
        : G4VSensitiveDetector(name)
    {
        collectionName.insert(hitsCollectionName);
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    SiDetectorSD::~SiDetectorSD()
    {
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void SiDetectorSD::Initialize(G4HCofThisEvent *hce)
    {
        // Create hits collection

        fHitsCollection = new SiDetectorHitsCollection(SensitiveDetectorName, collectionName[0]);

        // Add this collection in hce

        G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
        hce->AddHitsCollection(hcID, fHitsCollection);
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    G4bool SiDetectorSD::ProcessHits(G4Step *aStep,
                                     G4TouchableHistory *)
    {
        // energy deposit
        G4double edep = aStep->GetTotalEnergyDeposit();

        if (edep == 0.)
            return false;

        SiDetectorHit *newHit = new SiDetectorHit();
	G4int copyNo = 
		aStep->GetPreStepPoint()
		     ->GetTouchableHandle()
		     ->GetCopyNumber();
	
	newHit->SetCopyNo(copyNo);

        // Track ID
        newHit->SetTrackID(aStep->GetTrack()->GetTrackID());

        // Energy deposition
        newHit->SetEnergyDeposit(edep);

        // Hit time
        newHit->SetTime(aStep->GetTrack()->GetGlobalTime());

        // Geant4 ion identity and mass. Raytrace receives the same mass
        // converted to its PMASS convention by the MdmTrace wrapper.
        const G4ParticleDefinition *particle =
            aStep->GetTrack()->GetParticleDefinition();
        newHit->SetMass(particle->GetPDGMass() / (931.48 * MeV));
        newHit->SetIonMassMeV(particle->GetPDGMass() / MeV);
        newHit->SetMassNumber(particle->GetAtomicMass());
        newHit->SetAtomicNumber(particle->GetAtomicNumber());

        // Particle charge
        newHit->SetCharge(aStep->GetTrack()->GetDynamicParticle()->GetCharge() / eplus);

        // Particle kinetic energy
        newHit->SetKineticEnergy(aStep->GetTrack()->GetDynamicParticle()->GetKineticEnergy());



        // Hit position
        const G4ThreeVector globalPosition =
            aStep->GetPreStepPoint()->GetPosition();

        newHit->SetGlobalPosition(globalPosition);

        // Transform from global to detector-local coordinates
        const auto& transform =
            aStep->GetPreStepPoint()
                 ->GetTouchableHandle()
                 ->GetHistory()
                 ->GetTopTransform();

        const G4ThreeVector localPosition =
            transform.TransformPoint(globalPosition);

        newHit->SetLocalPosition(localPosition);

        // Hit momentum
        const G4ThreeVector globalMomentum =
            aStep->GetTrack()->GetMomentum();

        newHit->SetGlobalMomentum(globalMomentum);

        const G4ThreeVector localMomentum =
            transform.TransformAxis(globalMomentum);

        newHit->SetLocalMomentum(localMomentum);


        // Hit position
	/*
        newHit->SetGlobalPosition(aStep->GetPreStepPoint()->GetPosition());
        G4ThreeVector translatedPosition = aStep->GetPreStepPoint()->GetPosition() - aStep->GetPreStepPoint()->GetPhysicalVolume()->GetObjectTranslation();
        G4ThreeVector rotatedPosition = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetObjectRotation()->inverse().operator()(translatedPosition);
        newHit->SetLocalPosition(rotatedPosition);

        // Hit momentum
        newHit->SetGlobalMomentum(aStep->GetTrack()->GetMomentum());
        G4ThreeVector rotatedMomentum = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetObjectRotation()->inverse().operator()(aStep->GetTrack()->GetMomentum());
        newHit->SetLocalMomentum(rotatedMomentum);
        */
        fHitsCollection->insert(newHit);

        // newHit->Print();

        return true;
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void SiDetectorSD::EndOfEvent(G4HCofThisEvent *)
    {
        return;
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
