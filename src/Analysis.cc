#include "Analysis.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace MdmSim
{

    Analysis *Analysis::fInstance = nullptr;

    Analysis::Analysis()
        : fFactoryOn(false)
    {
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    Analysis::~Analysis()
    {
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    Analysis *Analysis::Instance()
    {
        if (!fInstance)
        {
            fInstance = new Analysis();
        }
        return fInstance;
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void Analysis::Book(std::string processNumber)
    {

        if (!fFactoryOn)
        {
            // Open a root file
            std::string fileName = "SimData_" + processNumber + ".root";
            fFile = new TFile(fileName.c_str(), "RECREATE");
            if (!fFile)
            {
                G4cerr << "\n---> Analysis::Book(): cannot open "
                       << fileName << G4endl;
                return;
            }

            fFile->cd();

            // Define two trees, one for accurate hit info, the other one for what you will get in experiment.
            fTreeAccurate = new TTree("AccurateData", "Accurate data recorded");
            fTreeExperiment = new TTree("ExperimentalData", "Experimental data recorded by DAQ");

            // Add branches to trees
            fTreeAccurate->Branch("SlitHitPresent", &fSlitHitPresent);
            fTreeAccurate->Branch("SlitHitTransmitted", &fSlitHitTransmitted);
            fTreeAccurate->Branch("SlitHitTrackId", &fSlitHitTrackId);
            fTreeAccurate->Branch("SlitHitLocalPosX", &fSlitHitLocalPosX);
            fTreeAccurate->Branch("SlitHitLocalPosY", &fSlitHitLocalPosY);
            fTreeAccurate->Branch("SlitHitLocalPosZ", &fSlitHitLocalPosZ);
            fTreeAccurate->Branch("RaytracePositionX", &fRaytracePositionX);
            fTreeAccurate->Branch("RaytracePositionY", &fRaytracePositionY);
            fTreeAccurate->Branch("RaytraceAngleX", &fRaytraceAngleX);
            fTreeAccurate->Branch("RaytraceAngleY", &fRaytraceAngleY);

            fTreeAccurate->Branch("RaytraceFocalPlaneAccepted", &fRaytraceFocalPlaneAccepted);
            fTreeAccurate->Branch("RaytraceFocalPlaneHitTrackId", &fRaytraceFocalPlaneHitTrackId);
            fTreeAccurate->Branch("RaytraceFocalPlaneHitTime", &fRaytraceFocalPlaneHitTime);
            fTreeAccurate->Branch("RaytraceFocalPlaneHitGlobalPosX", &fRaytraceFocalPlaneHitGlobalPosX);
            fTreeAccurate->Branch("RaytraceFocalPlaneHitGlobalPosY", &fRaytraceFocalPlaneHitGlobalPosY);
            fTreeAccurate->Branch("RaytraceFocalPlaneHitGlobalPosZ", &fRaytraceFocalPlaneHitGlobalPosZ);
            fTreeAccurate->Branch("RaytraceFocalPlaneHitLocalPosX", &fRaytraceFocalPlaneHitLocalPosX);
            fTreeAccurate->Branch("RaytraceFocalPlaneHitLocalPosY", &fRaytraceFocalPlaneHitLocalPosY);
            fTreeAccurate->Branch("RaytraceFocalPlaneHitLocalPosZ", &fRaytraceFocalPlaneHitLocalPosZ);
            fTreeAccurate->Branch("RaytraceFocalPlaneHitGlobalMomentumX", &fRaytraceFocalPlaneHitGlobalMomentumX);
            fTreeAccurate->Branch("RaytraceFocalPlaneHitGlobalMomentumY", &fRaytraceFocalPlaneHitGlobalMomentumY);
            fTreeAccurate->Branch("RaytraceFocalPlaneHitGlobalMomentumZ", &fRaytraceFocalPlaneHitGlobalMomentumZ);
            fTreeAccurate->Branch("RaytraceFocalPlaneHitLocalMomentumX", &fRaytraceFocalPlaneHitLocalMomentumX);
            fTreeAccurate->Branch("RaytraceFocalPlaneHitLocalMomentumY", &fRaytraceFocalPlaneHitLocalMomentumY);
            fTreeAccurate->Branch("RaytraceFocalPlaneHitLocalMomentumZ", &fRaytraceFocalPlaneHitLocalMomentumZ);

            fTreeAccurate->Branch("SiDeltaEHitAccepted", &fSiDeltaEHitAccepted);
            fTreeAccurate->Branch("SiDeltaEHitTrackId", &fSiDeltaEHitTrackId);
            fTreeAccurate->Branch("SiDeltaEHitEDep", &fSiDeltaEHitEDep);
            fTreeAccurate->Branch("SiDeltaEHitTime", &fSiDeltaEHitTime);

            fTreeAccurate->Branch("SiEHitAccepted", &fSiEHitAccepted);
            fTreeAccurate->Branch("SiEHitTrackId", &fSiEHitTrackId);
            fTreeAccurate->Branch("SiEHitEDep", &fSiEHitEDep);
            fTreeAccurate->Branch("SiEHitTime", &fSiEHitTime);
            fTreeAccurate->Branch("SiEHitGlobalPosX", &fSiEHitGlobalPosX);
            fTreeAccurate->Branch("SiEHitGlobalPosY", &fSiEHitGlobalPosY);
            fTreeAccurate->Branch("SiEHitGlobalPosZ", &fSiEHitGlobalPosZ);
            fTreeAccurate->Branch("SiEHitLocalPosX", &fSiEHitLocalPosX);
            fTreeAccurate->Branch("SiEHitLocalPosY", &fSiEHitLocalPosY);
            fTreeAccurate->Branch("SiEHitLocalPosZ", &fSiEHitLocalPosZ);
            fTreeAccurate->Branch("SiEHitGlobalMomentumX", &fSiEHitGlobalMomentumX);
            fTreeAccurate->Branch("SiEHitGlobalMomentumY", &fSiEHitGlobalMomentumY);
            fTreeAccurate->Branch("SiEHitGlobalMomentumZ", &fSiEHitGlobalMomentumZ);
            fTreeAccurate->Branch("SiEHitLocalMomentumX", &fSiEHitLocalMomentumX);
            fTreeAccurate->Branch("SiEHitLocalMomentumY", &fSiEHitLocalMomentumY);
            fTreeAccurate->Branch("SiEHitLocalMomentumZ", &fSiEHitLocalMomentumZ);

            fTreeAccurate->Branch("Ppac1Accepted", &fPpac1Accepted);
            fTreeAccurate->Branch("Ppac1HitTrackId", &fPpac1HitTrackId);
            fTreeAccurate->Branch("Ppac1HitEDep", &fPpac1HitEDep);
            fTreeAccurate->Branch("Ppac1HitTime", &fPpac1HitTime);
            fTreeAccurate->Branch("Ppac1HitGlobalPosX", &fPpac1HitGlobalPosX);
            fTreeAccurate->Branch("Ppac1HitGlobalPosY", &fPpac1HitGlobalPosY);
            fTreeAccurate->Branch("Ppac1HitGlobalPosZ", &fPpac1HitGlobalPosZ);
            fTreeAccurate->Branch("Ppac1HitLocalPosX", &fPpac1HitLocalPosX);
            fTreeAccurate->Branch("Ppac1HitLocalPosY", &fPpac1HitLocalPosY);
            fTreeAccurate->Branch("Ppac1HitLocalPosZ", &fPpac1HitLocalPosZ);
            fTreeAccurate->Branch("Ppac1HitGlobalMomentumX", &fPpac1HitGlobalMomentumX);
            fTreeAccurate->Branch("Ppac1HitGlobalMomentumY", &fPpac1HitGlobalMomentumY);
            fTreeAccurate->Branch("Ppac1HitGlobalMomentumZ", &fPpac1HitGlobalMomentumZ);
            fTreeAccurate->Branch("Ppac1HitLocalMomentumX", &fPpac1HitLocalMomentumX);
            fTreeAccurate->Branch("Ppac1HitLocalMomentumY", &fPpac1HitLocalMomentumY);
            fTreeAccurate->Branch("Ppac1HitLocalMomentumZ", &fPpac1HitLocalMomentumZ);

            fTreeAccurate->Branch("Ppac2Accepted", &fPpac2Accepted);
            fTreeAccurate->Branch("Ppac2HitTrackId", &fPpac2HitTrackId);
            fTreeAccurate->Branch("Ppac2HitEDep", &fPpac2HitEDep);
            fTreeAccurate->Branch("Ppac2HitTime", &fPpac2HitTime);
            fTreeAccurate->Branch("Ppac2HitGlobalPosX", &fPpac2HitGlobalPosX);
            fTreeAccurate->Branch("Ppac2HitGlobalPosY", &fPpac2HitGlobalPosY);
            fTreeAccurate->Branch("Ppac2HitGlobalPosZ", &fPpac2HitGlobalPosZ);
            fTreeAccurate->Branch("Ppac2HitLocalPosX", &fPpac2HitLocalPosX);
            fTreeAccurate->Branch("Ppac2HitLocalPosY", &fPpac2HitLocalPosY);
            fTreeAccurate->Branch("Ppac2HitLocalPosZ", &fPpac2HitLocalPosZ);
            fTreeAccurate->Branch("Ppac2HitGlobalMomentumX", &fPpac2HitGlobalMomentumX);
            fTreeAccurate->Branch("Ppac2HitGlobalMomentumY", &fPpac2HitGlobalMomentumY);
            fTreeAccurate->Branch("Ppac2HitGlobalMomentumZ", &fPpac2HitGlobalMomentumZ);
            fTreeAccurate->Branch("Ppac2HitLocalMomentumX", &fPpac2HitLocalMomentumX);
            fTreeAccurate->Branch("Ppac2HitLocalMomentumY", &fPpac2HitLocalMomentumY);
            fTreeAccurate->Branch("Ppac2HitLocalMomentumZ", &fPpac2HitLocalMomentumZ);

            fTreeExperiment->Branch("DaqTrigger", &fDaqTrigger, "DaqTrigger/O");
            fTreeExperiment->Branch("SiDeltaEHitEDep", &fSiDeltaEHitEDepExp);
            fTreeExperiment->Branch("SiDeltaEHitTime", &fSiDeltaEHitTimeExp);
            fTreeExperiment->Branch("SiEHitEDep", &fSiEHitEDepExp);
            fTreeExperiment->Branch("SiEHitTime", &fSiEHitTimeExp);
            fTreeExperiment->Branch("SiEHitHorizontalNo", &fSiEHitHorizontalNo);
            fTreeExperiment->Branch("SiEHitVerticalNo", &fSiEHitVerticalNo);

            fFactoryOn = true;
        }

        G4cout << "\n----> Output file is open in "
               << fFile->GetName() << "."
               << G4endl;
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void Analysis::Save()
    {
        if (!fFactoryOn)
            return;
        fFile->cd();
        fTreeAccurate->Write();
        fTreeExperiment->Write();
        fFile->Close();
        fFactoryOn = false;

        G4cout << "\n----> TTrees are saved.\n"
               << G4endl;
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void Analysis::FillTreeAccurate()
    {
        fTreeAccurate->Fill();
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void Analysis::FillTreeExperiment()
    {
        fTreeExperiment->Fill();
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//     void Analysis::SetPpacAccepted(G4bool *acc)
//     {
//         fPpacAccepted[0]=acc[0];
//         fPpacAccepted[1]=acc[2];
//     }
//     void Analysis::SetPpacHitTrackId(std::vector<G4int> *v)
//     {
//         fPpacHitTrackId[0]=v[0];
//         fPpacHitTrackId[1]=v[1];
//     }
//     void Analysis::SetPpacHitEDep(std::vector<G4double> *v)
//     {
//         fPpacHitEDep[0]=v[0];
//         fPpacHitEDep[1]=v[1];
//     }
//     void Analysis::SetPpacHitTime(std::vector<G4double> *v)
//     {
//         fPpacHitTime[0]=v[0];
//         fPpacHitTime[1]=v[1];
//     }
//     void Analysis::SetPpacHitGlobalPosX(std::vector<G4double> *v)
//     {
//         fPpacHitGlobalPosX[0]=v[0];
//         fPpacHitGlobalPosX[1]=v[1];
//     }
//     void Analysis::SetPpacHitGlobalPosY(std::vector<G4double> *v)
//     {
//         fPpacHitGlobalPosY[0]=v[0];
//         fPpacHitGlobalPosY[1]=v[1];
//     }void Analysis::SetPpacHitGlobalPosZ(std::vector<G4double> *v)
//     {
//         fPpacHitGlobalPosZ[0]=v[0];
//         fPpacHitGlobalPosZ[1]=v[1];
//     }
//    void Analysis::SetPpacHitLocalPosX(std::vector<G4double> *v)
//     {
//         fPpacHitLocalPosX[0]=v[0];
//         fPpacHitLocalPosX[1]=v[1];
//     }
//     void Analysis::SetPpacHitLocalPosY(std::vector<G4double> *v)
//     {
//         fPpacHitLocalPosY[0]=v[0];
//         fPpacHitLocalPosY[1]=v[1];
//     }void Analysis::SetPpacHitLocalPosZ(std::vector<G4double> *v)
//     {
//         fPpacHitLocalPosZ[0]=v[0];
//         fPpacHitLocalPosZ[1]=v[1];
//     }
//     void Analysis::SetPpacHitGlobalMomentumX(std::vector<G4double> *v)
//     {
//         fPpacHitGlobalMomentumX[0]=v[0];
//         fPpacHitGlobalMomentumX[1]=v[1];
//     }
//     void Analysis::SetPpacHitGlobalMomentumY(std::vector<G4double> *v)
//     {
//         fPpacHitGlobalMomentumY[0]=v[0];
//         fPpacHitGlobalMomentumY[1]=v[1];
//     }void Analysis::SetPpacHitGlobalMomentumZ(std::vector<G4double> *v)
//     {
//         fPpacHitGlobalMomentumZ[0]=v[0];
//         fPpacHitGlobalMomentumZ[1]=v[1];
//     }
//    void Analysis::SetPpacHitLocalMomentumX(std::vector<G4double> *v)
//     {
//         fPpacHitLocalMomentumX[0]=v[0];
//         fPpacHitLocalMomentumX[1]=v[1];
//     }
//     void Analysis::SetPpacHitLocalMomentumY(std::vector<G4double> *v)
//     {
//         fPpacHitLocalMomentumY[0]=v[0];
//         fPpacHitLocalMomentumY[1]=v[1];
//     }void Analysis::SetPpacHitLocalMomentumZ(std::vector<G4double> *v)
//     {
//         fPpacHitLocalMomentumZ[0]=v[0];
//         fPpacHitLocalMomentumZ[1]=v[1];
//     }
    

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
