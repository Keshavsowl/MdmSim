#ifndef MdmSimDetectorConstruction_h
#define MdmSimDetectorConstruction_h 1

#include "G4FieldManager.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4RotationMatrix.hh"
#include "MdmFieldMapMagneticField.hh"
#include "MdmChargeStateMagneticEqRhs.hh"

#include "globals.hh"

class G4VPhysicalVolume;
class G4Material;
class G4VSensitiveDetector;
class G4VisAttributes;
class G4GenericMessenger;
/// Detector construction class to define materials and geometry.

namespace MdmSim
{
    class DetectorConstruction : public G4VUserDetectorConstruction
    {
    public:
        DetectorConstruction();
        ~DetectorConstruction() override;

        G4VPhysicalVolume *Construct() override;
        void ConstructSDandField() override;
        void ParseParams(std::map<std::string, G4double> params);
        void SetFieldMapPaths(const MdmFieldMapPaths &paths);
        G4double GetTargetThicknessLength() const;
        G4double GetTargetRotationAngle() const { return fTargetRotationAngle; }
        // void SetMdmAngle(G4double ang);

    private:
        static G4ThreadLocal MdmFieldMapMagneticField *fDipoleField;
        static G4ThreadLocal MdmFieldMapMagneticField *fFirstMultipoleField;

        static G4ThreadLocal G4FieldManager *fDipoleFieldMgr;
        static G4ThreadLocal G4FieldManager *fFirstMultipoleFieldMgr;

        G4double fTargetRotationAngle;
        G4double fTargetThickness;

        G4bool fUseDeltaE;
        G4double fSiDetectorAngle;

        G4double fMdmAngle;
        G4int fBeamZ = 0;
        G4int fBeamA = 0;
        G4double fBeamCharge = 0.;
        G4bool fReactionEnabled = false;
        G4int fReactionLightZ = 0;
        G4int fReactionLightA = 0;
        G4double fReactionLightCharge = 0.;
        G4double fReactionLightEx = 0.;
        G4int fReactionHeavyZ = 0;
        G4int fReactionHeavyA = 0;
        G4double fReactionHeavyCharge = 0.;
        G4double fReactionHeavyEx = 0.;

        G4double fFirstMultipoleProbe;
        G4double fDipoleProbe;
        MdmFieldMapPaths fFieldMapPaths;

        G4LogicalVolume *fLogicFirstMultipoleField = nullptr;
        G4VPhysicalVolume *fPhysicFirstMultipoleField = nullptr;
        G4ThreeVector fMultipoleFieldPos;

        G4LogicalVolume *fLogicDipoleField = nullptr;
        G4VPhysicalVolume *fPhysicDipoleField = nullptr;
        G4ThreeVector fDipoleFieldPos;

        G4LogicalVolume *fLogicDipoleMagnet = nullptr;
        G4VPhysicalVolume *fPhysicDipoleMagnet = nullptr;
        G4LogicalVolume *fLogicRaytraceFocalPlane = nullptr;

        G4double fPpacGasPressure;
        G4double fPpacLength;
        G4RotationMatrix *fPpacChamberRot = nullptr;
        G4ThreeVector fPpacChamberPos;
    };

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
}
#endif
