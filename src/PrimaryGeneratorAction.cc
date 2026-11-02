#include "PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"
#include "G4ChargedGeantino.hh"

#include <cmath>
#include <stdexcept>

namespace MdmSim
{

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    PrimaryGeneratorAction::PrimaryGeneratorAction()
    {
        G4int n_particle = 1;
        fParticleGun = new G4ParticleGun(n_particle);

        // default particle kinematic
        G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
        G4ParticleDefinition *particle = particleTable->FindParticle("chargedgeantino");
        fParticleGun->SetParticleDefinition(particle);
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    PrimaryGeneratorAction::~PrimaryGeneratorAction()
    {
        delete fParticleGun;
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void PrimaryGeneratorAction::ParseParams(std::map<std::string, G4double> params)
    {
        const bool hasBeamCharge = params.find("BeamCharge") != params.end();
        if (params.find("BeamZ") != params.end())
            fBeamZ = static_cast<G4int>(std::lround(params["BeamZ"]));
        if (params.find("BeamA") != params.end())
            fBeamA = static_cast<G4int>(std::lround(params["BeamA"]));
        if (hasBeamCharge)
            fBeamCharge = static_cast<G4int>(std::lround(params["BeamCharge"]));
        else
            fBeamCharge = fBeamZ;
        if (params.find("BeamEnergyInMeV") != params.end())
            fBeamEnergy = params["BeamEnergyInMeV"] * MeV;
        if (params.find("MdmAngleInDeg") != params.end())
            fMdmAngle = params["MdmAngleInDeg"] * deg;

        if (fBeamZ <= 0 || fBeamA <= 0)
            throw std::runtime_error("BeamZ and BeamA must be positive");
        if (fBeamCharge < 0 || fBeamCharge > fBeamZ)
            throw std::runtime_error("BeamCharge must satisfy 0 <= BeamCharge <= BeamZ");
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void PrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent)
    {
        // this function is called at the begining of ecah event
        //

        // In order to avoid dependence of PrimaryGeneratorAction
        // on DetectorConstruction class we get Envelope volume
        // from G4LogicalVolumeStore.
        G4ParticleDefinition *particle = fParticleGun->GetParticleDefinition();
        if (particle == G4ChargedGeantino::ChargedGeantino())
        {
            G4ParticleDefinition *ion;
            G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
            if (particleTable->GetIonTable()->FindIon(fBeamZ, fBeamA, 0.0))
                ion = particleTable->GetIonTable()->FindIon(fBeamZ, fBeamA, 0.0);
            else
                ion = particleTable->GetIonTable()->GetIon(fBeamZ, fBeamA, 0.0);
            // ion = particleTable->FindParticle("alpha");
            fParticleGun->SetParticleDefinition(ion);
        }
        fParticleGun->SetParticleCharge(fBeamCharge * eplus);

        // G4double xAngle = G4RandFlat::shoot(fMdmAngle - 2. * deg, fMdmAngle + 2. * deg);
        // G4double yAngle = G4RandFlat::shoot(-2. * deg, 2. * deg);
        G4double xAngle = 0.;
        G4double yAngle = 0.;
        G4ThreeVector direction(std::tan(xAngle), std::tan(yAngle) * std::sqrt(std::tan(xAngle) * std::tan(xAngle) + 1.), 1.);
        // G4ThreeVector direction(0., 0., 1.);
        // direction.setTheta(std::acos(G4RandFlat::shoot(std::cos(8 * deg), std::cos(0. * deg))));
        // G4double phi = G4RandFlat::shoot(-180. * deg, 180. * deg);
        // G4double phi = G4RandFlat::shoot() > 0.5 ? 90. * deg : -90. * deg;
        // direction.setPhi(phi);
        fParticleGun->SetParticleMomentumDirection(direction);

        fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., -5. * cm));
        // fParticleGun->SetParticleEnergy(G4RandGauss::shoot(20. * MeV, 0.005 * 20. * MeV / 2.355));
        fParticleGun->SetParticleEnergy(fBeamEnergy);

        fParticleGun->GeneratePrimaryVertex(anEvent);
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
