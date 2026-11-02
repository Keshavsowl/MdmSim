#ifndef MdmSimSiDetectorHit_h
#define MdmSimSiDetectorHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "tls.hh"

namespace MdmSim
{

    /// SiDetector hit class
    ///
    /// It defines data members to store the trackID, energy deposit, kinetic energy, momentum,
    /// and position of charged particles in a selected volume:
    /// - fTrackID, fEnergyDeposit, fKineticEnergy, fMomentum, fPosition

    class SiDetectorHit : public G4VHit
    {
    public:
        SiDetectorHit();
        SiDetectorHit(const SiDetectorHit &) = default;
        ~SiDetectorHit() override;

        // operators
        SiDetectorHit &operator=(const SiDetectorHit &) = default;
        G4bool operator==(const SiDetectorHit &) const;

        inline void *operator new(size_t);
        inline void operator delete(void *);

        // methods from base class
        void Print() override;

        // Set methods
        void SetTrackID(G4int track) { fTrackID = track; };
	void SetCopyNo(G4int copyNo) {fCopyNo = copyNo;};
        void SetEnergyDeposit(G4double de) { fEnergyDeposit = de; };
        void SetTime(G4double t) { fTime = t; };
        void SetMass(G4double mass) { fMass = mass; };
        void SetIonMassMeV(G4double mass) { fIonMassMeV = mass; };
        void SetMassNumber(G4int massNumber) { fMassNumber = massNumber; };
        void SetAtomicNumber(G4int atomicNumber) { fAtomicNumber = atomicNumber; };
        void SetCharge(G4double charge) { fCharge = charge; };
        void SetKineticEnergy(G4double e) { fKineticEnergy = e; };
        void SetGlobalPosition(G4ThreeVector xyz) { fGlobalPosition = xyz; };
        void SetLocalPosition(G4ThreeVector xyz) { fLocalPosition = xyz; };
        void SetGlobalMomentum(G4ThreeVector xyz) { fGlobalMomentum = xyz; };
        void SetLocalMomentum(G4ThreeVector xyz) { fLocalMomentum = xyz; };

        // Get methods
        G4int GetTrackID() const { return fTrackID; };
	G4int GetCopyNo() const {return fCopyNo;};
        G4double GetEnergyDeposit() const { return fEnergyDeposit; };
        G4double GetTime() const { return fTime; };
        G4double GetMass() const { return fMass; };
        G4double GetIonMassMeV() const { return fIonMassMeV; };
        G4int GetMassNumber() const { return fMassNumber; };
        G4int GetAtomicNumber() const { return fAtomicNumber; };
        G4double GetCharge() const { return fCharge; };
        G4double GetKineticEnergy() const { return fKineticEnergy; };
        G4ThreeVector GetLocalMomentum() const { return fLocalMomentum; };
        G4ThreeVector GetGlobalMomentum() const { return fGlobalMomentum; };
        G4ThreeVector GetGlobalPosition() const { return fGlobalPosition; };
        G4ThreeVector GetLocalPosition() const { return fLocalPosition; };

    private:
        G4int fTrackID;
	G4int fCopyNo;
        G4double fEnergyDeposit;
        G4double fTime;
        G4double fMass;
        G4double fIonMassMeV;
        G4int fMassNumber;
        G4int fAtomicNumber;
        G4double fCharge;
        G4double fKineticEnergy;
        G4ThreeVector fGlobalPosition;
        G4ThreeVector fLocalPosition;
        G4ThreeVector fGlobalMomentum;
        G4ThreeVector fLocalMomentum;
    };

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    typedef G4THitsCollection<SiDetectorHit> SiDetectorHitsCollection;

    extern G4ThreadLocal G4Allocator<SiDetectorHit> *SiDetectorHitAllocator;

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    inline void *SiDetectorHit::operator new(size_t)
    {
        if (!SiDetectorHitAllocator)
            SiDetectorHitAllocator = new G4Allocator<SiDetectorHit>;
        return (void *)SiDetectorHitAllocator->MallocSingle();
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    inline void SiDetectorHit::operator delete(void *hit)
    {
        SiDetectorHitAllocator->FreeSingle((SiDetectorHit *)hit);
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}

#endif
