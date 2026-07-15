#include "DetectorConstruction.hh"
#include "SiDetectorSD.hh"
#include "PpacSD.hh"
#include "TrackingPlaneSD.hh"
#include "Constants.hh"
#include "MdmChargeStateMagneticEqRhs.hh"

#include "G4FieldManager.hh"
#include "G4ChordFinder.hh"
#include "G4ClassicalRK4.hh"
#include "G4TransportationManager.hh"
#include "G4Field.hh"
#include "G4AutoDelete.hh"

#include "G4GenericMessenger.hh"

#include "G4NistManager.hh"

#include "G4RunManager.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4UserLimits.hh"
#include "G4SDManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4NucleiProperties.hh"

#include <cmath>
#include <vector>

namespace
{
    G4double GetMdmSimIonMass(G4int z, G4int a, G4double excitationEnergy)
    {
        auto *particleTable = G4ParticleTable::GetParticleTable();
        if (excitationEnergy == 0.)
        {
            if (z == 0 && a == 1)
                return particleTable->FindParticle("neutron")->GetPDGMass();
            if (z == 1 && a == 1)
                return particleTable->FindParticle("proton")->GetPDGMass();
            if (z == 1 && a == 2)
                return particleTable->FindParticle("deuteron")->GetPDGMass();
            if (z == 1 && a == 3)
                return particleTable->FindParticle("triton")->GetPDGMass();
            if (z == 2 && a == 3)
                return particleTable->FindParticle("He3")->GetPDGMass();
            if (z == 2 && a == 4)
                return particleTable->FindParticle("alpha")->GetPDGMass();
        }

        if (z == 0)
            return particleTable->FindParticle("neutron")->GetPDGMass() + excitationEnergy;
        return G4NucleiProperties::GetNuclearMass(a, z) + excitationEnergy;
    }
}

namespace MdmSim
{
    G4ThreadLocal MdmFieldMapMagneticField *DetectorConstruction::fDipoleField = nullptr;
    G4ThreadLocal G4FieldManager *DetectorConstruction::fDipoleFieldMgr = nullptr;
    G4ThreadLocal MdmFieldMapMagneticField *DetectorConstruction::fFirstMultipoleField = nullptr;
    G4ThreadLocal G4FieldManager *DetectorConstruction::fFirstMultipoleFieldMgr = nullptr;
    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    DetectorConstruction::DetectorConstruction()
    {
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    DetectorConstruction::~DetectorConstruction()
    {
        delete fDipoleField;
        delete fDipoleFieldMgr;
        delete fFirstMultipoleField;
        delete fFirstMultipoleFieldMgr;
        delete fPpacChamberRot;
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    G4VPhysicalVolume *DetectorConstruction::Construct()
    {
        // Get nist material manager
        G4NistManager *nist = G4NistManager::Instance();

        // Option to switch on/off checking of volumes overlaps
        //
        G4bool checkOverlaps = false;

        //
        // World
        //
        G4double world_sizeXY = 10. * m;
        G4double world_sizeZ = 15. * m;
        G4Material *world_mat = nist->FindOrBuildMaterial("G4_Galactic");

        G4Box *solidWorld =
            new G4Box("World",                                                    // its name
                      0.5 * world_sizeXY, 0.5 * world_sizeXY, 0.5 * world_sizeZ); // its size

        G4LogicalVolume *logicWorld =
            new G4LogicalVolume(solidWorld, // its solid
                                world_mat,  // its material
                                "World");   // its name

        G4VPhysicalVolume *physWorld =
            new G4PVPlacement(0,               // no rotation
                              G4ThreeVector(), // at (0,0,0)
                              logicWorld,      // its logical volume
                              "World",         // its name
                              0,               // its mother  volume
                              false,           // no boolean operation
                              0,               // copy number
                              checkOverlaps);  // overlaps checking

        //
        // Target
        //
        G4Material *targetMaterial = nist->FindOrBuildMaterial("G4_Galactic");
        G4RotationMatrix *targetRot = new G4RotationMatrix;
        targetRot->rotateY(fTargetRotationAngle);
        G4double targetRMin = 0.;
        G4double targetRMax = 3.5 * mm;
        G4double targetDz = fTargetThickness / (2.253 * g / cm3);
        G4cout << "\n---> DetectorConstruction::Construct(): targetDz= " << G4BestUnit(targetDz, "Length") << G4endl;
        G4Tubs *solidTarget = new G4Tubs("Target",                                                      // its name
                                         targetRMin, targetRMax, 0.5 * targetDz, 2. * M_PI, 2. * M_PI); // its size

        G4LogicalVolume *logicTarget =
            new G4LogicalVolume(solidTarget,    // its solid
                                targetMaterial, // its material
                                "Target");      // its name

        new G4PVPlacement(targetRot,       // rotation
                          G4ThreeVector(), // at (0,0,0)
                          logicTarget,     // its logical volume
                          "Target",        // its name
                          logicWorld,      // its mother  volume
                          false,           // no boolean operation
                          0,               // copy number
                          checkOverlaps);  // overlaps checking

        logicTarget->SetUserLimits(new G4UserLimits(0.1 * targetDz));

        //
        // Si detector - Delta E
        //
        G4Material *DeltaEMaterial = fUseDeltaE ? nist->FindOrBuildMaterial("G4_Si") : nist->FindOrBuildMaterial("G4_Galactic");
        G4ThreeVector DeltaEPos = G4ThreeVector(0., 0., 13.724 * cm);
        DeltaEPos.setTheta(fSiDetectorAngle);
        DeltaEPos.setPhi(180. * deg);
        // Rotate the dE detector
        G4RotationMatrix *DeltaERot = new G4RotationMatrix;
        DeltaERot->rotateY(fSiDetectorAngle);
        // Box shape
        G4double DeltaEX = 5. * cm;
        G4double DeltaEY = 5. * cm;
        G4double DeltaEZ = 32. * um;

        G4Box *solidDeltaE =
            new G4Box("DetdE",
                      0.5 * DeltaEX, 0.5 * DeltaEY, 0.5 * DeltaEZ);

        G4LogicalVolume *logicDeltaE =
            new G4LogicalVolume(solidDeltaE,         // its solid
                                DeltaEMaterial,      // its material
                                "SiDetectorDeltaE"); // its name

        new G4PVPlacement(DeltaERot,          // rotation
                          DeltaEPos,          // at position
                          logicDeltaE,        // its logical volume
                          "SiDetectorDeltaE", // its name
                          logicWorld,         // its mother  volume
                          false,              // no boolean operation
                          0,                  // copy number
                          checkOverlaps);     // overlaps checking

        logicDeltaE->SetUserLimits(new G4UserLimits(0.5 * DeltaEZ));

        //
        // Si detector - E
        //
        G4ThreeVector EPos = G4ThreeVector(0., 0., 14.224 * cm);
        EPos.setTheta(fSiDetectorAngle);
        EPos.setPhi(180. * deg);
        // Rotate the dE detector
        G4RotationMatrix *ERot = new G4RotationMatrix;
        ERot->rotateY(fSiDetectorAngle);
        // Box shape
        G4double EX = 5. * cm;
        G4double EY = 5. * cm;
        G4double EZ = 500. * um;

        G4Box *solidE =
            new G4Box("DetdE",
                      0.5 * EX, 0.5 * EY, 0.5 * EZ);

        G4LogicalVolume *logicE =
            new G4LogicalVolume(solidE,                             // its solid
                                nist->FindOrBuildMaterial("G4_Si"), // its material
                                "SiDetectorE");                     // its name

        new G4PVPlacement(ERot,           // rotation
                          EPos,           // at position
                          logicE,         // its logical volume
                          "SiDetectorE",  // its name
                          logicWorld,     // its mother  volume
                          false,          // no boolean operation
                          0,              // copy number
                          checkOverlaps); // overlaps checking

        logicE->SetUserLimits(new G4UserLimits(0.5 * EZ));

        //
        // MDM slit and slit box
        //
        G4ThreeVector slitBoxPos = G4ThreeVector(0., 0., kSlitBoxPos); // position in mother frame
        slitBoxPos.rotateY(fMdmAngle);
        G4RotationMatrix *slitBoxRot = new G4RotationMatrix;
        slitBoxRot->rotateY(fMdmAngle);

        G4VSolid *solidSlit = new G4Box("SlitSolid", 2.27965 * cm, 2.27965 * cm, 0.5 * kSlitBoxDz); // measured on 1/26/2022
        G4LogicalVolume *logicSlit = new G4LogicalVolume(solidSlit, nist->FindOrBuildMaterial("G4_Galactic"), "SlitLogical");
        new G4PVPlacement(G4Transform3D(*slitBoxRot, slitBoxPos), logicSlit, "SlitPhysical", logicWorld, false, 0, checkOverlaps);

        G4VSolid *solidSlitBoxShape = new G4Box("SlitBoxShape", (2.27965 + 50.) * cm, (2.27965 + 10.) * cm, 0.5 * kSlitBoxDz);
        G4VSolid *solidSlitBox = new G4SubtractionSolid("SlitBox", solidSlitBoxShape, solidSlit);
        G4LogicalVolume *logicSlitBox = new G4LogicalVolume(solidSlitBox, nist->FindOrBuildMaterial("G4_Cu"), "SlitBoxLogical");
        new G4PVPlacement(G4Transform3D(*slitBoxRot, slitBoxPos), logicSlitBox, "SlitBoxPhysical", logicWorld, false, 0, checkOverlaps);

        //
        // First multipole field
        //
        G4VSolid *solidFirstMultipoleField = new G4Tubs("FirstMultipoleFieldTubs", 0., kFirstMultipoleAperture, 0.5 * kFirstMultipoleLength, 0., 360. * deg);
        fLogicFirstMultipoleField = new G4LogicalVolume(solidFirstMultipoleField, nist->FindOrBuildMaterial("G4_Galactic"), "FirstMultipoleFieldLogical");
        fMultipoleFieldPos = G4ThreeVector(0., 0., kFirstMultipoleEntrancePos + 0.5 * kFirstMultipoleLength);
        fMultipoleFieldPos.rotateY(fMdmAngle);
        G4RotationMatrix *firstMultipoleFieldRot = new G4RotationMatrix;
        firstMultipoleFieldRot->rotateY(fMdmAngle);
        fPhysicFirstMultipoleField = new G4PVPlacement(G4Transform3D(*firstMultipoleFieldRot, fMultipoleFieldPos), fLogicFirstMultipoleField,
                                                       "FirstMultipoleFieldPhysical", logicWorld,
                                                       false, 0, checkOverlaps);
        fLogicFirstMultipoleField->SetUserLimits(new G4UserLimits(1.0 * mm));

        //
        // First multipole magnet
        //
        G4VSolid *solidFirstMultipoleMagnetShape = new G4Tubs("FirstMultipoleMagnetShape", 0., kFirstMultipoleAperture + 0.1 * cm, 0.5 * kFirstMultipoleLength, 0., 360. * deg);
        auto solidFirstMultipoleMagnet = new G4SubtractionSolid("FirstMultipoleMagnetSolid", solidFirstMultipoleMagnetShape, solidFirstMultipoleField);
        auto logicFirstMultipoleMagnet = new G4LogicalVolume(solidFirstMultipoleMagnet, nist->FindOrBuildMaterial("G4_Al"), "FirstMultipoleMagnetLogical");
        new G4PVPlacement(G4Transform3D(*firstMultipoleFieldRot, fMultipoleFieldPos), logicFirstMultipoleMagnet,
                          "FirstMultipoleMagnetPhysical", logicWorld,
                          false, 0, checkOverlaps);

        //
        // Dipole field
        //
        // The entrance/exit fringe maps overlap the sector map. The field
        // envelope only needs the parts outside the sector; inside the sector,
        // the same field object is called by the sector volume and dispatches
        // to the entrance/exit maps where MDMTrace would.
        G4VSolid *solidDipoleFieldEntrace = new G4Box("DipoleFieldEntrance", 0.5 * kDipoleFieldWidth, 0.5 * kDipoleZ11, 0.5 * kDipoleFieldHeight);
        G4ThreeVector solidDipoleFieldEntraceTrans(kDipoleFieldRadius, -0.5 * kDipoleZ11, 0.);
        G4VSolid *solidDipoleFieldExit = new G4Box("DipoleFieldExit", 0.5 * kDipoleFieldWidth, 0.5 * kDipoleZ22, 0.5 * kDipoleFieldHeight);
        G4RotationMatrix *solidDipoleFieldExitRot = new G4RotationMatrix;
        solidDipoleFieldExitRot->rotateZ(kDipoleDeflectionAngle);
        G4ThreeVector solidDipoleFieldExitTrans1(-kDipoleFieldRadius * std::cos(M_PI - kDipoleDeflectionAngle), kDipoleFieldRadius * std::sin(M_PI - kDipoleDeflectionAngle), 0.);
        G4ThreeVector solidDipoleFieldExitTrans2(-0.5 * kDipoleZ22 * std::cos(kDipoleDeflectionAngle - M_PI / 2.), -0.5 * kDipoleZ22 * std::sin(kDipoleDeflectionAngle - M_PI / 2.), 0.);
        G4ThreeVector solidDipoleFieldExitTrans = solidDipoleFieldExitTrans1 + solidDipoleFieldExitTrans2;
        G4VSolid *solidDipoleFieldUniform = new G4Tubs("DipoleFieldUniform", kDipoleFieldRadius - 0.5 * kDipoleFieldWidth, kDipoleFieldRadius + 0.5 * kDipoleFieldWidth, 0.5 * kDipoleFieldHeight, 0., kDipoleDeflectionAngle);
        G4UnionSolid *solidDipoleFieldTemp = new G4UnionSolid("solidDipoleFieldTemp", solidDipoleFieldUniform, solidDipoleFieldEntrace, nullptr, solidDipoleFieldEntraceTrans);
        G4UnionSolid *solidDipoleField = new G4UnionSolid("solidDipoleField", solidDipoleFieldTemp, solidDipoleFieldExit, G4Transform3D(*solidDipoleFieldExitRot, solidDipoleFieldExitTrans));
        fLogicDipoleField = new G4LogicalVolume(solidDipoleField, nist->FindOrBuildMaterial("G4_Galactic"), "DipoleFieldLogical");
        fDipoleFieldPos.setX(-kDipoleFieldRadius);
        fDipoleFieldPos.setY(0.);
        fDipoleFieldPos.setZ(kFirstArmLength);
        fDipoleFieldPos.rotateY(fMdmAngle);
        G4RotationMatrix *dipoleRot = new G4RotationMatrix;
        dipoleRot->rotateX(90. * deg);
        dipoleRot->rotateY(fMdmAngle);
        fPhysicDipoleField = new G4PVPlacement(G4Transform3D(*dipoleRot, fDipoleFieldPos), fLogicDipoleField,
                                               "DipoleFieldPhysical", logicWorld,
                                               false, 0, checkOverlaps);
        fLogicDipoleField->SetUserLimits(new G4UserLimits(1.0 * mm));

        //
        // Dipole magnet
        //
        auto solidDipoleMagnetShape = new G4Tubs("DipoleMagnetShape", 30. * cm, kDipoleMagnetRadius, 0.5 * m, 0., kDipoleDeflectionAngle);
        auto solidDipoleMagnet = new G4SubtractionSolid("DipoleMagnetSolid", solidDipoleMagnetShape, solidDipoleField);
        fLogicDipoleMagnet = new G4LogicalVolume(solidDipoleMagnet, nist->FindOrBuildMaterial("G4_Cu"), "DipoleMagnetLogical");
        fPhysicDipoleMagnet = new G4PVPlacement(G4Transform3D(*dipoleRot, fDipoleFieldPos), fLogicDipoleMagnet,
                                                "DipoleMagnetPhysical", logicWorld,
                                                false, 0, checkOverlaps);

        //
        // PPAC
        //
        // PPAC chamber
        G4Material *ppacChamberMaterial = new G4Material("Pentane", (0.00385 * fPpacGasPressure / (133.32236842 * pascal) + 1.666e-5) * kg / m3, nist->FindOrBuildMaterial("G4_N-PENTANE"), kStateGas, 300. * kelvin, fPpacGasPressure);
        G4VSolid *solidPpacChamber = new G4Box("PpacChamberSolid", kPpacWidth / 2., kPpacHeight / 2., kPpacChamberLength / 2.);
        G4LogicalVolume *logicPpacChamber = new G4LogicalVolume(solidPpacChamber, ppacChamberMaterial, "PpacChamberLogical");
        G4ThreeVector solidPpacChamberTrans1(-kDipoleFieldRadius - kDipoleFieldRadius * std::cos(M_PI - kDipoleDeflectionAngle), 0., kDipoleFieldRadius * std::sin(M_PI - kDipoleDeflectionAngle));
        solidPpacChamberTrans1 += G4ThreeVector(0., 0., kFirstArmLength);
        G4ThreeVector solidPpacChamberTrans2(-(kSecondArmLength + 60. * cm) * std::cos(kDipoleDeflectionAngle - M_PI / 2.), 0., -(kSecondArmLength + 60. * cm) * std::sin(kDipoleDeflectionAngle - M_PI / 2.));
        fPpacChamberPos = solidPpacChamberTrans1 + solidPpacChamberTrans2;
        fPpacChamberPos.rotateY(fMdmAngle);
        fPpacChamberRot = new G4RotationMatrix;
        fPpacChamberRot->rotateY(-kDipoleDeflectionAngle + fMdmAngle);
        new G4PVPlacement(G4Transform3D(*fPpacChamberRot, fPpacChamberPos), logicPpacChamber,
                          "PpacChamberPhysical", logicWorld,
                          false, 0, checkOverlaps);

        G4VSolid *solidRaytraceFocalPlane = new G4Box("RaytraceFocalPlaneSolid",
                                                    kPpacWidth / 2.,
                                                    kPpacHeight / 2.,
                                                    kRaytraceFocalPlaneThickness / 2.);
        fLogicRaytraceFocalPlane = new G4LogicalVolume(solidRaytraceFocalPlane,
                                                     nist->FindOrBuildMaterial("G4_Galactic"),
                                                     "RaytraceFocalPlaneLogical");
        const G4ThreeVector raytraceFocalPlaneLocal(0., 0., kRaytraceFocalPlaneLocalZ);
        const G4ThreeVector raytraceFocalPlaneWorld = (*fPpacChamberRot)(raytraceFocalPlaneLocal) + fPpacChamberPos;
        new G4PVPlacement(G4Transform3D(*fPpacChamberRot, raytraceFocalPlaneWorld),
                          fLogicRaytraceFocalPlane,
                          "RaytraceFocalPlanePhysical",
                          logicWorld,
                          false,
                          0,
                          checkOverlaps);

        // PPAC entrance window
        G4VSolid *solidEntranceWindow = new G4Box("EntranceWindowSolid", kPpacWidth / 2., kPpacHeight / 2., kPpacEntranceWindowThickness / 2.);
        G4LogicalVolume *logicEntranceWindow = new G4LogicalVolume(solidEntranceWindow, nist->FindOrBuildMaterial("G4_MYLAR"), "EntranceWindowLogical");
        G4PVPlacement *physicalEntranceWindow = new G4PVPlacement(nullptr, G4ThreeVector(0., 0., 0.5 * (-kPpacChamberLength + kPpacEntranceWindowThickness)), logicEntranceWindow,
                                                                  "EntranceWindowPhysical", logicPpacChamber,
                                                                  false, 0, checkOverlaps);
        // G4cout << "physicalEntranceWindow->Pos(): " << G4BestUnit(physicalEntranceWindow->GetObjectTranslation(), "Length") << G4endl;

        // PPAC cathode Al layer
        G4Material *cathodeAlMaterial = new G4Material("Al", 2.70 * g / cm3, nist->FindOrBuildMaterial("G4_Al"));
        G4double cathodeAlDz = kPpacCathodeAlThickness / (2.70 * g / cm3);
        G4VSolid *solidCathodeAl = new G4Box("CathodeAlSolid", kPpacWidth / 2., kPpacHeight / 2., cathodeAlDz / 2.);
        G4LogicalVolume *logicCathodeAl = new G4LogicalVolume(solidCathodeAl, cathodeAlMaterial, "CathodeAlLogical");

        // PPAC cathode mylar layer
        G4Material *cathodeMylarMaterial = nist->FindOrBuildMaterial("G4_MYLAR");
        G4double cathodeMylarDz = kPpacCathodeMylarThickness / cathodeMylarMaterial->GetDensity();
        G4cout << "cathodeMylarDz=" << G4BestUnit(cathodeMylarDz, "Length") << G4endl;
        G4VSolid *solidCathode0Mylar = new G4Box("Cathode0MylarSolid", kPpacWidth / 2., kPpacHeight / 2., cathodeMylarDz / 2.);
        G4LogicalVolume *logicCathode0Mylar = new G4LogicalVolume(solidCathode0Mylar, cathodeMylarMaterial, "Cathode0MylarLogical");
        G4VSolid *solidCathode1Mylar = new G4Box("Cathode1MylarSolid", kPpacWidth / 2., kPpacHeight / 2., cathodeMylarDz / 2.);
        G4LogicalVolume *logicCathode1Mylar = new G4LogicalVolume(solidCathode1Mylar, cathodeMylarMaterial, "Cathode1MylarLogical");

        // place PPAC cathode 1 (Al_0 + Mylar_0 + Al_1)
        G4ThreeVector cathodeMylar0Pos(0., 0., physicalEntranceWindow->GetObjectTranslation().z() + kPpacSpacingWindowCathode);
        new G4PVPlacement(nullptr, cathodeMylar0Pos, logicCathode0Mylar,
                          "Cathode0MylarPhysical",
                          logicPpacChamber, false, 0, checkOverlaps);
        new G4PVPlacement(nullptr, G4ThreeVector(0., 0., cathodeMylar0Pos.z() - 0.5 * (cathodeAlDz + cathodeMylarDz)), logicCathodeAl,
                          "CathodeAlPhysical",
                          logicPpacChamber, false, 0, checkOverlaps);
        new G4PVPlacement(nullptr, G4ThreeVector(0., 0., cathodeMylar0Pos.z() + 0.5 * (cathodeAlDz + cathodeMylarDz)), logicCathodeAl,
                          "CathodeAlPhysical",
                          logicPpacChamber, false, 1, checkOverlaps);

        // place PPAC cathode 2 (Al_2 + Mylar_1 + Al_3)
        G4ThreeVector cathodeMylar1Pos(0., 0., cathodeMylar0Pos.z() + fPpacLength);
        new G4PVPlacement(nullptr, cathodeMylar1Pos, logicCathode1Mylar,
                          "Cathode1MylarPhysical",
                          logicPpacChamber, false, 0, checkOverlaps);
        new G4PVPlacement(nullptr, G4ThreeVector(0., 0., cathodeMylar1Pos.z() - 0.5 * (cathodeAlDz + cathodeMylarDz)), logicCathodeAl,
                          "CathodeAlPhysical",
                          logicPpacChamber, false, 2, checkOverlaps);
        new G4PVPlacement(nullptr, G4ThreeVector(0., 0., cathodeMylar1Pos.z() + 0.5 * (cathodeAlDz + cathodeMylarDz)), logicCathodeAl,
                          "CathodeAlPhysical",
                          logicPpacChamber, false, 3, checkOverlaps);
        //
        // always return the physical World
        //
        return physWorld;
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void DetectorConstruction::ConstructSDandField()
    {
        // Sensitive detector
        SiDetectorSD *aSiDetectorESD = new SiDetectorSD("MdmSim/SiDetectorESD", "SiDetectorEHitsCollection");
        G4SDManager::GetSDMpointer()->AddNewDetector(aSiDetectorESD);
        SetSensitiveDetector("SiDetectorE", aSiDetectorESD, true);

        SiDetectorSD *aSiDetectorDeltaESD = new SiDetectorSD("MdmSim/SiDetectorDeltaESD", "SiDetectorDeltaEHitsCollection");
        G4SDManager::GetSDMpointer()->AddNewDetector(aSiDetectorDeltaESD);
        SetSensitiveDetector("SiDetectorDeltaE", aSiDetectorDeltaESD, true);

        SiDetectorSD *aSlitSD = new SiDetectorSD("MdmSim/SlitSD", "SlitHitsCollection");
        G4SDManager::GetSDMpointer()->AddNewDetector(aSlitSD);
        SetSensitiveDetector("SlitLogical", aSlitSD, true);

        PpacSD *aPpac1SD = new PpacSD("MdmSim/Ppac1SD", "Ppac1HitsCollection");
        aPpac1SD->SetPpacChamberPosRot(fPpacChamberRot, fPpacChamberPos);
        G4SDManager::GetSDMpointer()->AddNewDetector(aPpac1SD);
        SetSensitiveDetector("Cathode0MylarLogical", aPpac1SD, true);

        PpacSD *aPpac2SD = new PpacSD("MdmSim/Ppac2SD", "Ppac2HitsCollection");
        aPpac2SD->SetPpacChamberPosRot(fPpacChamberRot, fPpacChamberPos);
        G4SDManager::GetSDMpointer()->AddNewDetector(aPpac2SD);
        SetSensitiveDetector("Cathode1MylarLogical", aPpac2SD, true);

        TrackingPlaneSD *raytraceFocalPlaneSD = new TrackingPlaneSD("MdmSim/RaytraceFocalPlaneSD", "RaytraceFocalPlaneHitsCollection");
        raytraceFocalPlaneSD->SetReferenceFrame(fPpacChamberRot, fPpacChamberPos);
        G4SDManager::GetSDMpointer()->AddNewDetector(raytraceFocalPlaneSD);
        SetSensitiveDetector("RaytraceFocalPlaneLogical", raytraceFocalPlaneSD, true);

        // Magnetic field ----------------------------------------------------------
        std::vector<MdmChargeOverride> chargeOverrides;
        if (fReactionEnabled)
        {
            auto addChargeOverride = [&](G4int z, G4int a, G4double excitationEnergy, G4double chargeState) {
                const G4double mass = GetMdmSimIonMass(z, a, excitationEnergy);
                chargeOverrides.push_back({mass, chargeState});
                G4cout << "Set: MDM magnetic transport override for "
                       << "Z=" << z << " A=" << a << " = "
                       << chargeState << " e" << G4endl;
            };

            addChargeOverride(fBeamZ, fBeamA, 0., fBeamCharge);
            addChargeOverride(fReactionLightZ, fReactionLightA,
                              fReactionLightEx, fReactionLightCharge);
            addChargeOverride(fReactionHeavyZ, fReactionHeavyA,
                              fReactionHeavyEx, fReactionHeavyCharge);
        }

        // First multipole
        fFirstMultipoleField = MdmFieldMapMagneticField::CreateMultipole(fFieldMapPaths.multipole, fMdmAngle, fMultipoleFieldPos, fDipoleProbe, fFirstMultipoleProbe);
        fFirstMultipoleFieldMgr = new G4FieldManager();
        fFirstMultipoleFieldMgr->SetDetectorField(fFirstMultipoleField);
        auto *firstMultipoleEquation = fReactionEnabled
                                           ? new MdmChargeStateMagneticEqRhs(fFirstMultipoleField, fBeamCharge, chargeOverrides)
                                           : new MdmChargeStateMagneticEqRhs(fFirstMultipoleField, fBeamCharge);
        auto *firstMultipoleStepper = new G4ClassicalRK4(firstMultipoleEquation);
        auto *firstMultipoleChordFinder = new G4ChordFinder(fFirstMultipoleField, 0.5 * mm, firstMultipoleStepper);
        fFirstMultipoleFieldMgr->SetChordFinder(firstMultipoleChordFinder);
        fFirstMultipoleFieldMgr->SetAccuraciesWithDeltaOneStep(1e-3 * mm);
        G4bool forceToAllDaughters = true;
        fLogicFirstMultipoleField->SetFieldManager(fFirstMultipoleFieldMgr, forceToAllDaughters);
        // Dipole
        fDipoleField = MdmFieldMapMagneticField::CreateDipole(fFieldMapPaths, fMdmAngle, fDipoleProbe, fFirstMultipoleProbe);
        fDipoleFieldMgr = new G4FieldManager();
        fDipoleFieldMgr->SetDetectorField(fDipoleField);
        auto *dipoleEquation = fReactionEnabled
                                   ? new MdmChargeStateMagneticEqRhs(fDipoleField, fBeamCharge, chargeOverrides)
                                   : new MdmChargeStateMagneticEqRhs(fDipoleField, fBeamCharge);
        auto *dipoleStepper = new G4ClassicalRK4(dipoleEquation);
        auto *dipoleChordFinder = new G4ChordFinder(fDipoleField, 0.5 * mm, dipoleStepper);
        fDipoleFieldMgr->SetChordFinder(dipoleChordFinder);
        fDipoleFieldMgr->SetAccuraciesWithDeltaOneStep(1e-3 * mm);
        fLogicDipoleField->SetFieldManager(fDipoleFieldMgr, forceToAllDaughters);
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void DetectorConstruction::SetFieldMapPaths(const MdmFieldMapPaths &paths)
    {
        fFieldMapPaths = paths;
        G4cout << "Set: Multipole field map = " << fFieldMapPaths.multipole << G4endl;
        G4cout << "Set: Dipole entrance field map = " << fFieldMapPaths.dipoleEntrance << G4endl;
        G4cout << "Set: Dipole sector field map = " << fFieldMapPaths.dipoleSector << G4endl;
        G4cout << "Set: Dipole exit field map = " << fFieldMapPaths.dipoleExit << G4endl;
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    G4double DetectorConstruction::GetTargetThicknessLength() const
    {
        return fTargetThickness / (2.253 * g / cm3);
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void DetectorConstruction::ParseParams(std::map<std::string, G4double> params)
    {
        G4cout << "\n---> DetectorConstruction::ParseParams():" << G4endl;
        for (auto it : params)
        {
            if (it.first == "TargetRotationAngleInDeg")
            {
                fTargetRotationAngle = it.second * deg;
                G4cout << "Set: Target rotation angle = " << G4BestUnit(fTargetRotationAngle, "Angle") << G4endl;
            }
            else if (it.first == "TargetThicknessInMgCm2")
            {
                fTargetThickness = it.second * mg / cm2;
                G4cout << "Set: Target thickness = " << G4BestUnit(fTargetThickness, "Mass/Surface") << G4endl;
            }
            else if (it.first == "UseDeltaE")
            {
                fUseDeltaE = it.second;
                G4cout << std::boolalpha << "Set: Use Delta-E detector = " << it.second << G4endl;
            }
            else if (it.first == "SiDetectorAngleInDeg")
            {
                fSiDetectorAngle = it.second * deg;
                G4cout << "Set: Si detector angle = " << G4BestUnit(fSiDetectorAngle, "Angle") << G4endl;
            }
            else if (it.first == "MdmAngleInDeg")
            {
                fMdmAngle = it.second * deg;
                G4cout << "Set: MDM angle = " << G4BestUnit(fMdmAngle, "Angle") << G4endl;
            }
            else if (it.first == "BeamZ")
            {
                fBeamZ = static_cast<G4int>(std::lround(it.second));
            }
            else if (it.first == "BeamA")
            {
                fBeamA = static_cast<G4int>(std::lround(it.second));
            }
            else if (it.first == "BeamCharge")
            {
                fBeamCharge = it.second;
                G4cout << "Set: MDM magnetic transport charge state = " << fBeamCharge << " e" << G4endl;
            }
            else if (it.first == "ReactionEnabled")
            {
                fReactionEnabled = it.second != 0.;
                G4cout << "Set: Reaction mode = " << (fReactionEnabled ? "enabled" : "disabled") << G4endl;
            }
            else if (it.first == "ReactionLightZ")
            {
                fReactionLightZ = static_cast<G4int>(std::lround(it.second));
            }
            else if (it.first == "ReactionLightA")
            {
                fReactionLightA = static_cast<G4int>(std::lround(it.second));
            }
            else if (it.first == "ReactionLightCharge")
            {
                fReactionLightCharge = it.second;
            }
            else if (it.first == "ReactionLightExMeV")
            {
                fReactionLightEx = it.second * MeV;
            }
            else if (it.first == "ReactionHeavyZ")
            {
                fReactionHeavyZ = static_cast<G4int>(std::lround(it.second));
            }
            else if (it.first == "ReactionHeavyA")
            {
                fReactionHeavyA = static_cast<G4int>(std::lround(it.second));
            }
            else if (it.first == "ReactionHeavyCharge")
            {
                fReactionHeavyCharge = it.second;
            }
            else if (it.first == "ReactionHeavyExMeV")
            {
                fReactionHeavyEx = it.second * MeV;
            }
            else if (it.first == "FirstMultipoleProbe")
            {
                fFirstMultipoleProbe = it.second;
                G4cout << "Set: Multipole probe = " << fFirstMultipoleProbe << " Gauss" << G4endl;
            }
            else if (it.first == "DipoleProbe")
            {
                fDipoleProbe = it.second;
                G4cout << "Set: Dipole probe = " << fDipoleProbe << " Gauss" << G4endl;
            }
            else if (it.first == "PpacGasPressureInTorr")
            {
                fPpacGasPressure = it.second * 133.32236842 * pascal;
                G4cout << "Set: PPAC gas pressure = " << G4BestUnit(fPpacGasPressure, "Pressure") << G4endl;
            }
            else if (it.first == "PpacLengthInCm")
            {
                fPpacLength = it.second * cm;
                G4cout << "Set: PPAC length = " << G4BestUnit(fPpacLength, "Length") << G4endl;
            }
        }
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    /*
        void DetectorConstruction::DefineCommands()
        {
            // Define /MdmSim/detector command directory using generic messenger class
            fMessenger = new G4GenericMessenger(this,
                                                "/MdmSim/detector/",
                                                "Detector control");

            // MdmAngle command
            auto &MdmAngleCmd = fMessenger->DeclareMethodWithUnit("MdmAngle", "deg",
                                                                  &DetectorConstruction::SetMdmAngle,
                                                                  "Set rotation angle of the MDM.");
            MdmAngleCmd.SetParameterName("angle", true);
            MdmAngleCmd.SetRange("angle>=-50. && angle<=150.");
            MdmAngleCmd.SetDefaultValue("0.");
        }
    */

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    /*
        void DetectorConstruction::SetMdmAngle(G4double val)
    {
        if (!fPhysicDipoleField)
        {
            G4cerr << "\n---> DetectorConstruction::SetMdmAngle(): Detector has not yet been constructed." << G4endl;
            return;
        }

        fMdmAngle = val;
        *fMdmRotation = G4RotationMatrix(); // make it unit vector
        fMdmRotation->rotateY(-fMdmAngle);
        fMdmRotation->rotateX(-90. * deg);
        G4ThreeVector dipolePos(-kDipoleFieldRadius, 0., kFirstArmLength);
        dipolePos.rotateY(-fMdmAngle);
        fPhysicDipoleField->SetTranslation(dipolePos);
        fPhysicDipoleMagnet->SetTranslation(dipolePos);
        printf("\n---> DetectorConstruction::SetMdmAngle(): fMdmAngle: %.2f\n", fMdmAngle);

        // tell G4RunManager that we change the geometry
        G4RunManager::GetRunManager()->GeometryHasBeenModified();
    }
    */
}
