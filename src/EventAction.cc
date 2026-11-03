#include "EventAction.hh"
#include "RunAction.hh"
#include "Analysis.hh"
#include "SiDetectorHit.hh"
#include "PpacHit.hh"

#include "MdmTrace.h"

#include "G4SDManager.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include <cmath>

namespace
{
    MdmIon BuildMdmIonFromHit(const MdmSim::SiDetectorHit &hit)
    {
        MdmIon ion;
        ion.massNumber = hit.GetMassNumber();
        ion.atomicNumber = hit.GetAtomicNumber();
        ion.chargeState = static_cast<int>(std::lround(hit.GetCharge()));
        ion.ionMassMeV = hit.GetIonMassMeV();
        return ion;
    }
}

namespace MdmSim
{

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    EventAction::EventAction(RunAction *runAction)
        : fRunAction(runAction)
    {
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    EventAction::~EventAction()
    {
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void EventAction::BeginOfEventAction(const G4Event *)
    {   
	Analysis *analysis = Analysis::Instance();
	

	analysis->ResetAlpha1Generated();
	analysis->ResetAlpha2Generated();
	analysis->ResetGammaGenerated();
	analysis->ResetCoincidence();

        G4SDManager *sdManager = G4SDManager::GetSDMpointer();
        //fHCID_SiDetectorDeltaE = sdManager->GetCollectionID("SiDetectorDeltaEHitsCollection");
        fHCID_SiDetectorE = sdManager->GetCollectionID("SiDetectorEHitsCollection");
	fHCID_NaI = sdManager->GetCollectionID("NaIHitsCollection");
        fHCID_Slit = sdManager->GetCollectionID("SlitHitsCollection");
        fHCID_RaytraceFocalPlane = sdManager->GetCollectionID("RaytraceFocalPlaneHitsCollection");
        fHCID_Ppac1 = sdManager->GetCollectionID("Ppac1HitsCollection");
        fHCID_Ppac2 = sdManager->GetCollectionID("Ppac2HitsCollection");

        fMdmTrace.SetMdmAngle(fMdmAngle); // deg
        fMdmTrace.SetMdmProbe(fDipoleProbe, fFirstMultipoleProbe);
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void EventAction::EndOfEventAction(const G4Event *aEvent)
    {
        Analysis *analysis = Analysis::Instance();

        G4HCofThisEvent *hce = aEvent->GetHCofThisEvent();
        // if (!hce)
        // {
        //     G4ExceptionDescription msg;
        //     msg << "No hits collection of this event found.\n";
        //     G4Exception("\n---> EventAction::EndOfEventAction()",
        //                 "MdmSimCode01", JustWarning, msg);
        //     return;
        // }

        //
        // Slit
        //
        SiDetectorHitsCollection *hcSlit = (SiDetectorHitsCollection *)hce->GetHC(fHCID_Slit);
        G4int nofHitsSlit = hcSlit->GetSize();
        G4bool slitHitPresent;
        G4bool slitHitTransmitted;
        std::vector<G4int> slitHitTrackId;
        std::vector<G4double> slitHitLocalPosX;
        std::vector<G4double> slitHitLocalPosY;
        std::vector<G4double> slitHitLocalPosZ;
        std::vector<G4double> raytracePositionX;
        std::vector<G4double> raytracePositionY;
        std::vector<G4double> raytraceAngleX;
        std::vector<G4double> raytraceAngleY;
        if (nofHitsSlit == 0)
        {
            // G4ExceptionDescription msg;
            // msg << "No hits in slit of this event found.\n";
            // G4Exception("\n---> EventAction::EndOfEventAction()",
            //             "MdmSimCode01", JustWarning, msg);
            slitHitPresent = false;
            slitHitTransmitted = false;
            slitHitTrackId.push_back(-999);
            slitHitLocalPosX.push_back(-999.);
            slitHitLocalPosY.push_back(-999.);
            slitHitLocalPosZ.push_back(-999.);
            raytracePositionX.push_back(-999.);
            raytracePositionY.push_back(-999.);
            raytraceAngleX.push_back(-999.);
            raytraceAngleY.push_back(-999.);
        }
        else
        {
            slitHitPresent = true;
            slitHitTransmitted = false;
            slitHitTrackId.push_back((*hcSlit)[0]->GetTrackID());
            G4ThreeVector slitLocalPos(0., 0., 0.);
            G4double mdmPosX = 0.;
            G4double mdmPosY = 0.;
            G4double mdmAngX = 0.;
            G4double mdmAngY = 0.;
            G4int itTimes = 0;
            for (int i = 0; i < nofHitsSlit; i++)
            {
                if ((*hcSlit)[i]->GetTrackID() == slitHitTrackId.back())
                {
                    slitLocalPos += (*hcSlit)[i]->GetLocalPosition();
                    G4ThreeVector momentum = (*hcSlit)[i]->GetLocalMomentum();
                    G4ThreeVector direction = momentum / momentum.mag();
                    G4double xAngle = std::atan(direction.x() / direction.z()) * 180. / M_PI;
                    G4double yAngle = std::atan(direction.y() / std::sqrt(std::pow(direction.x(), 2.) + std::pow(direction.z(), 2.))) * 180. / M_PI;
                    fMdmTrace.SetScatteredAngle(xAngle + fMdmAngle, yAngle);
                    fMdmTrace.SetScatteredIon(BuildMdmIonFromHit(*(*hcSlit)[i]));
                    fMdmTrace.SetScatteredEnergy((*hcSlit)[i]->GetKineticEnergy());
                    fMdmTrace.SendRay();
                    const G4double raytracePosX = fMdmTrace.GetFirstWireX() * 10.;
                    const G4double raytracePosY = fMdmTrace.GetFirstWireY() * 10.;
                    mdmPosX += raytracePosX;
                    mdmPosY += raytracePosY;
                    mdmAngX += fMdmTrace.GetFirstWireXAngle();
                    mdmAngY += fMdmTrace.GetFirstWireYAngle();
                    if (std::abs(raytracePosX) < 200. && std::abs(raytracePosY) < 50.)
                    {
                        slitHitTransmitted = true;
                    }
                    itTimes += 1;
                }
                else
                {
                    slitHitTrackId.push_back((*hcSlit)[i]->GetTrackID());
                    slitHitLocalPosX.push_back(slitLocalPos.x() / (double)itTimes);
                    slitHitLocalPosY.push_back(slitLocalPos.y() / (double)itTimes);
                    slitHitLocalPosZ.push_back(slitLocalPos.z() / (double)itTimes);
                    raytracePositionX.push_back(mdmPosX / (double)itTimes);
                    raytracePositionY.push_back(mdmPosY / (double)itTimes);
                    raytraceAngleX.push_back(mdmAngX / (double)itTimes);
                    raytraceAngleY.push_back(mdmAngY / (double)itTimes);
                    slitLocalPos = (*hcSlit)[i]->GetLocalPosition();
                    G4ThreeVector momentum = (*hcSlit)[i]->GetLocalMomentum();
                    G4ThreeVector direction = momentum / momentum.mag();
                    G4double xAngle = std::atan(direction.x() / direction.z()) * 180. / M_PI;
                    G4double yAngle = std::atan(direction.y() / std::sqrt(std::pow(direction.x(), 2.) + std::pow(direction.z(), 2.))) * 180. / M_PI;
                    fMdmTrace.SetScatteredAngle(xAngle + fMdmAngle, yAngle);
                    fMdmTrace.SetScatteredIon(BuildMdmIonFromHit(*(*hcSlit)[i]));
                    fMdmTrace.SetScatteredEnergy((*hcSlit)[i]->GetKineticEnergy());
                    fMdmTrace.SendRay();
                    mdmPosX = fMdmTrace.GetFirstWireX() * 10.;
                    mdmPosY = fMdmTrace.GetFirstWireY() * 10.;
                    mdmAngX = fMdmTrace.GetFirstWireXAngle();
                    mdmAngY = fMdmTrace.GetFirstWireYAngle();
                    if (std::abs(mdmPosX) < 200. && std::abs(mdmPosY) < 50.)
                    {
                        slitHitTransmitted = true;
                    }
                    itTimes = 1;
                }
            }
            slitHitLocalPosX.push_back(slitLocalPos.x() / (double)itTimes);
            slitHitLocalPosY.push_back(slitLocalPos.y() / (double)itTimes);
            slitHitLocalPosZ.push_back(slitLocalPos.z() / (double)itTimes);
            raytracePositionX.push_back(mdmPosX / (double)itTimes);
            raytracePositionY.push_back(mdmPosY / (double)itTimes);
            raytraceAngleX.push_back(mdmAngX / (double)itTimes);
            raytraceAngleY.push_back(mdmAngY / (double)itTimes);
        }
        analysis->SetSlitHitPresent(slitHitPresent);
        analysis->SetSlitHitTransmitted(slitHitTransmitted);
        analysis->SetSlitHitTrackId(slitHitTrackId);
        analysis->SetSlitHitLocalPosX(slitHitLocalPosX);
        analysis->SetSlitHitLocalPosY(slitHitLocalPosY);
        analysis->SetSlitHitLocalPosZ(slitHitLocalPosZ);
        analysis->SetRaytracePositionX(raytracePositionX);
        analysis->SetRaytracePositionY(raytracePositionY);
        analysis->SetRaytraceAngleX(raytraceAngleX);
        analysis->SetRaytraceAngleY(raytraceAngleY);

        auto collectTrackingPlaneHits = [](PpacHitsCollection *hits,
                                           G4bool &accepted,
                                           std::vector<G4int> &trackId,
                                           std::vector<G4double> &time,
                                           std::vector<G4double> &globalPosX,
                                           std::vector<G4double> &globalPosY,
                                           std::vector<G4double> &globalPosZ,
                                           std::vector<G4double> &localPosX,
                                           std::vector<G4double> &localPosY,
                                           std::vector<G4double> &localPosZ,
                                           std::vector<G4double> &globalMomentumX,
                                           std::vector<G4double> &globalMomentumY,
                                           std::vector<G4double> &globalMomentumZ,
                                           std::vector<G4double> &localMomentumX,
                                           std::vector<G4double> &localMomentumY,
                                           std::vector<G4double> &localMomentumZ)
        {
            const G4int nofHits = hits ? hits->GetSize() : 0;
            if (nofHits == 0)
            {
                accepted = false;
                trackId.push_back(-99);
                time.push_back(-99.);
                globalPosX.push_back(-999.);
                globalPosY.push_back(-999.);
                globalPosZ.push_back(-999.);
                localPosX.push_back(-999.);
                localPosY.push_back(-999.);
                localPosZ.push_back(-999.);
                globalMomentumX.push_back(-999.);
                globalMomentumY.push_back(-999.);
                globalMomentumZ.push_back(-999.);
                localMomentumX.push_back(-999.);
                localMomentumY.push_back(-999.);
                localMomentumZ.push_back(-999.);
                return;
            }

            accepted = true;
            G4int currentTrackId = (*hits)[0]->GetTrackID();
            G4double sumTime = 0.;
            G4ThreeVector sumGlobalPos(0., 0., 0.);
            G4ThreeVector sumLocalPos(0., 0., 0.);
            G4ThreeVector sumGlobalMomentum(0., 0., 0.);
            G4ThreeVector sumLocalMomentum(0., 0., 0.);
            G4int count = 0;

            const auto flush = [&]()
            {
                if (count == 0)
                {
                    return;
                }
                const G4double scale = 1. / static_cast<G4double>(count);
                trackId.push_back(currentTrackId);
                time.push_back(sumTime * scale);
                globalPosX.push_back(sumGlobalPos.x() * scale);
                globalPosY.push_back(sumGlobalPos.y() * scale);
                globalPosZ.push_back(sumGlobalPos.z() * scale);
                localPosX.push_back(sumLocalPos.x() * scale);
                localPosY.push_back(sumLocalPos.y() * scale);
                localPosZ.push_back(sumLocalPos.z() * scale);
                globalMomentumX.push_back(sumGlobalMomentum.x() * scale);
                globalMomentumY.push_back(sumGlobalMomentum.y() * scale);
                globalMomentumZ.push_back(sumGlobalMomentum.z() * scale);
                localMomentumX.push_back(sumLocalMomentum.x() * scale);
                localMomentumY.push_back(sumLocalMomentum.y() * scale);
                localMomentumZ.push_back(sumLocalMomentum.z() * scale);
            };

            for (G4int i = 0; i < nofHits; ++i)
            {
                PpacHit *hit = (*hits)[i];
                if (hit->GetTrackID() != currentTrackId)
                {
                    flush();
                    currentTrackId = hit->GetTrackID();
                    sumTime = 0.;
                    sumGlobalPos = G4ThreeVector(0., 0., 0.);
                    sumLocalPos = G4ThreeVector(0., 0., 0.);
                    sumGlobalMomentum = G4ThreeVector(0., 0., 0.);
                    sumLocalMomentum = G4ThreeVector(0., 0., 0.);
                    count = 0;
                }

                sumTime += hit->GetTime();
                sumGlobalPos += hit->GetGlobalPosition();
                sumLocalPos += hit->GetLocalPosition();
                sumGlobalMomentum += hit->GetGlobalMomentum();
                sumLocalMomentum += hit->GetLocalMomentum();
                ++count;
            }
            flush();
        };

        PpacHitsCollection *hcRaytraceFocalPlane = (PpacHitsCollection *)hce->GetHC(fHCID_RaytraceFocalPlane);
        G4bool raytraceFocalPlaneAccepted;
        std::vector<G4int> raytraceFocalPlaneHitTrackId;
        std::vector<G4double> raytraceFocalPlaneHitTime;
        std::vector<G4double> raytraceFocalPlaneHitGlobalPosX;
        std::vector<G4double> raytraceFocalPlaneHitGlobalPosY;
        std::vector<G4double> raytraceFocalPlaneHitGlobalPosZ;
        std::vector<G4double> raytraceFocalPlaneHitLocalPosX;
        std::vector<G4double> raytraceFocalPlaneHitLocalPosY;
        std::vector<G4double> raytraceFocalPlaneHitLocalPosZ;
        std::vector<G4double> raytraceFocalPlaneHitGlobalMomentumX;
        std::vector<G4double> raytraceFocalPlaneHitGlobalMomentumY;
        std::vector<G4double> raytraceFocalPlaneHitGlobalMomentumZ;
        std::vector<G4double> raytraceFocalPlaneHitLocalMomentumX;
        std::vector<G4double> raytraceFocalPlaneHitLocalMomentumY;
        std::vector<G4double> raytraceFocalPlaneHitLocalMomentumZ;
        collectTrackingPlaneHits(hcRaytraceFocalPlane,
                                 raytraceFocalPlaneAccepted,
                                 raytraceFocalPlaneHitTrackId,
                                 raytraceFocalPlaneHitTime,
                                 raytraceFocalPlaneHitGlobalPosX,
                                 raytraceFocalPlaneHitGlobalPosY,
                                 raytraceFocalPlaneHitGlobalPosZ,
                                 raytraceFocalPlaneHitLocalPosX,
                                 raytraceFocalPlaneHitLocalPosY,
                                 raytraceFocalPlaneHitLocalPosZ,
                                 raytraceFocalPlaneHitGlobalMomentumX,
                                 raytraceFocalPlaneHitGlobalMomentumY,
                                 raytraceFocalPlaneHitGlobalMomentumZ,
                                 raytraceFocalPlaneHitLocalMomentumX,
                                 raytraceFocalPlaneHitLocalMomentumY,
                                 raytraceFocalPlaneHitLocalMomentumZ);
        analysis->SetRaytraceFocalPlaneAccepted(raytraceFocalPlaneAccepted);
        analysis->SetRaytraceFocalPlaneHitTrackId(raytraceFocalPlaneHitTrackId);
        analysis->SetRaytraceFocalPlaneHitTime(raytraceFocalPlaneHitTime);
        analysis->SetRaytraceFocalPlaneHitGlobalPosX(raytraceFocalPlaneHitGlobalPosX);
        analysis->SetRaytraceFocalPlaneHitGlobalPosY(raytraceFocalPlaneHitGlobalPosY);
        analysis->SetRaytraceFocalPlaneHitGlobalPosZ(raytraceFocalPlaneHitGlobalPosZ);
        analysis->SetRaytraceFocalPlaneHitLocalPosX(raytraceFocalPlaneHitLocalPosX);
        analysis->SetRaytraceFocalPlaneHitLocalPosY(raytraceFocalPlaneHitLocalPosY);
        analysis->SetRaytraceFocalPlaneHitLocalPosZ(raytraceFocalPlaneHitLocalPosZ);
        analysis->SetRaytraceFocalPlaneHitGlobalMomentumX(raytraceFocalPlaneHitGlobalMomentumX);
        analysis->SetRaytraceFocalPlaneHitGlobalMomentumY(raytraceFocalPlaneHitGlobalMomentumY);
        analysis->SetRaytraceFocalPlaneHitGlobalMomentumZ(raytraceFocalPlaneHitGlobalMomentumZ);
        analysis->SetRaytraceFocalPlaneHitLocalMomentumX(raytraceFocalPlaneHitLocalMomentumX);
        analysis->SetRaytraceFocalPlaneHitLocalMomentumY(raytraceFocalPlaneHitLocalMomentumY);
        analysis->SetRaytraceFocalPlaneHitLocalMomentumZ(raytraceFocalPlaneHitLocalMomentumZ);

        //
        // Micron S2 Annular Silicon detectors
	// S2-1: front detector (CopyNo =0)
	// S2-2: Rear Detector (CopyNo =1)

        //
        // E
        SiDetectorHitsCollection *hcSiDetectorE = (SiDetectorHitsCollection *)hce->GetHC(fHCID_SiDetectorE);
        G4int nofHitsSiDetectorE = hcSiDetectorE->GetSize();
        G4bool siEHitAccepted;
        std::vector<G4int> siEHitTrackId;
	std::vector<G4int> siEHitCopyNo;
        std::vector<G4double> siEHitEDep;
        std::vector<G4double> siEHitTime;
        std::vector<G4double> siEHitGlobalPosX;
        std::vector<G4double> siEHitGlobalPosY;
        std::vector<G4double> siEHitGlobalPosZ;
        std::vector<G4double> siEHitLocalPosX;
        std::vector<G4double> siEHitLocalPosY;
        std::vector<G4double> siEHitLocalPosZ;
        std::vector<G4double> siEHitGlobalMomentumX;
        std::vector<G4double> siEHitGlobalMomentumY;
        std::vector<G4double> siEHitGlobalMomentumZ;
        std::vector<G4double> siEHitLocalMomentumX;
        std::vector<G4double> siEHitLocalMomentumY;
        std::vector<G4double> siEHitLocalMomentumZ;
        std::vector<G4double> siEHitEDepExp;
        std::vector<G4double> siEHitTimeExp;
        std::vector<G4int> siEHitHorizontalNo;
        std::vector<G4int> siEHitVerticalNo;
        if (nofHitsSiDetectorE == 0)
        {
            // G4ExceptionDescription msg;
            // msg << "No hits in Si detector E of this event found.\n";
            // G4Exception("\n---> EventAction::EndOfEventAction()",
            //             "MdmSimCode01", JustWarning, msg);
            siEHitAccepted = false;
            siEHitTrackId.push_back(-99);
	    siEHitCopyNo.push_back(-99);
            siEHitEDep.push_back(-99.);
            siEHitTime.push_back(-99.);
            siEHitGlobalPosX.push_back(-99.);
            siEHitGlobalPosY.push_back(-99.);
            siEHitGlobalPosZ.push_back(-99.);
            siEHitLocalPosX.push_back(-99.);
            siEHitLocalPosY.push_back(-99.);
            siEHitLocalPosZ.push_back(-99.);
            siEHitGlobalMomentumX.push_back(-99.);
            siEHitGlobalMomentumY.push_back(-99.);
            siEHitGlobalMomentumZ.push_back(-99.);
            siEHitLocalMomentumX.push_back(-99.);
            siEHitLocalMomentumY.push_back(-99.);
            siEHitLocalMomentumZ.push_back(-99.);
            siEHitEDepExp.push_back(-99.);
            siEHitTimeExp.push_back(-99.);
            siEHitHorizontalNo.push_back(-99);
            siEHitVerticalNo.push_back(-99);
        }
        else
        {
            siEHitAccepted = true;
            siEHitTrackId.push_back((*hcSiDetectorE)[0]->GetTrackID());
	    siEHitCopyNo.push_back((*hcSiDetectorE)[0]->GetCopyNo());
            G4double siEeDep = 0.;
            G4double siETime = 0.;
            G4ThreeVector siEGlobalPos(0., 0., 0.);
            G4ThreeVector siELocalPos(0., 0., 0.);
            G4ThreeVector siEGlobalMomentum(0., 0., 0.);
            G4ThreeVector siELocalMomentum(0., 0., 0.);
            G4int itTimes = 0;
            for (int i = 0; i < nofHitsSiDetectorE; i++)
            {
                if ((*hcSiDetectorE)[i]->GetTrackID() == siEHitTrackId.back() &&
				(*hcSiDetectorE)[i]->GetCopyNo() == siEHitCopyNo.back())
                {
                    siEeDep += (*hcSiDetectorE)[i]->GetEnergyDeposit();
                    siETime += (*hcSiDetectorE)[i]->GetTime();
                    siEGlobalPos += (*hcSiDetectorE)[i]->GetGlobalPosition();
                    siELocalPos += (*hcSiDetectorE)[i]->GetLocalPosition();
                    siEGlobalMomentum += (*hcSiDetectorE)[i]->GetGlobalMomentum();
                    siELocalMomentum += (*hcSiDetectorE)[i]->GetLocalMomentum();
                    itTimes += 1;
                }
                else
                {
                    siEHitTrackId.push_back((*hcSiDetectorE)[i]->GetTrackID());
		    siEHitCopyNo.push_back((*hcSiDetectorE)[i]->GetCopyNo());
                    siEHitEDep.push_back(siEeDep);
                    siEHitTime.push_back(siETime / (double)itTimes);
                    siEHitGlobalPosX.push_back(siEGlobalPos.x() / (double)itTimes);
                    siEHitGlobalPosY.push_back(siEGlobalPos.y() / (double)itTimes);
                    siEHitGlobalPosZ.push_back(siEGlobalPos.z() / (double)itTimes);
                    siEHitLocalPosX.push_back(siELocalPos.x() / (double)itTimes);
                    siEHitLocalPosY.push_back(siELocalPos.y() / (double)itTimes);
                    siEHitLocalPosZ.push_back(siELocalPos.z() / (double)itTimes);
                    siEHitGlobalMomentumX.push_back(siEGlobalMomentum.x() / (double)itTimes);
                    siEHitGlobalMomentumY.push_back(siEGlobalMomentum.y() / (double)itTimes);
                    siEHitGlobalMomentumZ.push_back(siEGlobalMomentum.z() / (double)itTimes);
                    siEHitLocalMomentumX.push_back(siELocalMomentum.x() / (double)itTimes);
                    siEHitLocalMomentumY.push_back(siELocalMomentum.y() / (double)itTimes);
                    siEHitLocalMomentumZ.push_back(siELocalMomentum.z() / (double)itTimes);
                    siEeDep = (*hcSiDetectorE)[i]->GetEnergyDeposit();
                    siETime = (*hcSiDetectorE)[i]->GetTime();
                    siEGlobalPos = (*hcSiDetectorE)[i]->GetGlobalPosition();
                    siELocalPos = (*hcSiDetectorE)[i]->GetLocalPosition();
                    siEGlobalMomentum = (*hcSiDetectorE)[i]->GetGlobalMomentum();
                    siELocalMomentum = (*hcSiDetectorE)[i]->GetLocalMomentum();
                    itTimes = 1;
                }
            }
            siEHitEDep.push_back(siEeDep);
            siEHitTime.push_back(siETime / (double)itTimes);
            siEHitGlobalPosX.push_back(siEGlobalPos.x() / (double)itTimes);
            siEHitGlobalPosY.push_back(siEGlobalPos.y() / (double)itTimes);
            siEHitGlobalPosZ.push_back(siEGlobalPos.z() / (double)itTimes);
            siEHitLocalPosX.push_back(siELocalPos.x() / (double)itTimes);
            siEHitLocalPosY.push_back(siELocalPos.y() / (double)itTimes);
            siEHitLocalPosZ.push_back(siELocalPos.z() / (double)itTimes);
            siEHitGlobalMomentumX.push_back(siEGlobalMomentum.x() / (double)itTimes);
            siEHitGlobalMomentumY.push_back(siEGlobalMomentum.y() / (double)itTimes);
            siEHitGlobalMomentumZ.push_back(siEGlobalMomentum.z() / (double)itTimes);
            siEHitLocalMomentumX.push_back(siELocalMomentum.x() / (double)itTimes);
            siEHitLocalMomentumY.push_back(siELocalMomentum.y() / (double)itTimes);
            siEHitLocalMomentumZ.push_back(siELocalMomentum.z() / (double)itTimes);
            for (int i = 0; i < siEHitTrackId.size(); i++)
            {
                siEHitEDepExp.push_back(G4RandGauss::shoot(siEHitEDep[i], fSiDetectorEnergyResolutionFwhmFraction * siEHitEDep[i] / 2.355));
                siEHitTimeExp.push_back(siEHitTime[i] + G4RandFlat::shoot(-fTdcResolutionHalfWidth, fTdcResolutionHalfWidth));
                
		
		//Micron S2 annular segmentation
		//Hit coordinates in detector-local frame
		G4double x = siEHitLocalPosX[i];
		G4double y = siEHitLocalPosY[i];

		//Radial distance from S2 center
		G4double r = std::sqrt(x*x + y*y);

		// Azimuthal angle, normalized to [0,360 degree]
		G4double phi = std::atan2(y,x);
		
		if (phi < 0.0)
			phi += 360.0 *deg;

		// S2 active radii
		G4double rInner = 1.10 *cm;
		G4double rOuter = 3.50 *cm;

		// 48 radial rings and 16 sectors
		G4double ringWidth = (rOuter-rInner)/48.0;
		G4double sectorWidth = 360.0 * deg / 16.0;

		G4int ringNo = -1;
		G4int sectorNo = -1;

		if (r >= rInner && r < rOuter)
		{
			ringNo = static_cast<G4int>(
					std::floor((r-rInner)/ringWidth));

			sectorNo = static_cast<G4int>(
					std::floor(phi/sectorWidth));
		}
		siEHitHorizontalNo.push_back(ringNo);
		siEHitVerticalNo.push_back(sectorNo);
		
		
		
		//G4int frontStripNo = std::floor((siEHitLocalPosX[i] + 2.5) / (5. / 16.));
                ///G4int backStripNo = std::floor((siEHitLocalPosY[i] + 2.5) / (5. / 16.));
                //siEHitHorizontalNo.push_back(backStripNo);
                //siEHitVerticalNo.push_back(frontStripNo);
            }
        }
        analysis->SetSiEHitAccepted(siEHitAccepted);
        analysis->SetSiEHitTrackId(siEHitTrackId);
	analysis->SetSiEHitCopyNo(siEHitCopyNo);
        analysis->SetSiEHitEDep(siEHitEDep);
        analysis->SetSiEHitTime(siEHitTime);
        analysis->SetSiEHitGlobalPosX(siEHitGlobalPosX);
        analysis->SetSiEHitGlobalPosY(siEHitGlobalPosY);
        analysis->SetSiEHitGlobalPosZ(siEHitGlobalPosZ);
        analysis->SetSiEHitLocalPosX(siEHitLocalPosX);
        analysis->SetSiEHitLocalPosY(siEHitLocalPosY);
        analysis->SetSiEHitLocalPosZ(siEHitLocalPosZ);
        analysis->SetSiEHitGlobalMomentumX(siEHitGlobalMomentumX);
        analysis->SetSiEHitGlobalMomentumY(siEHitGlobalMomentumY);
        analysis->SetSiEHitGlobalMomentumZ(siEHitGlobalMomentumZ);
        analysis->SetSiEHitLocalMomentumX(siEHitLocalMomentumX);
        analysis->SetSiEHitLocalMomentumY(siEHitLocalMomentumY);
        analysis->SetSiEHitLocalMomentumZ(siEHitLocalMomentumZ);
        analysis->SetDaqTrigger(siEHitAccepted);
        analysis->SetSiEHitEDepExp(siEHitEDepExp);
        analysis->SetSiEHitTimeExp(siEHitTimeExp);
        analysis->SetSiEHitHorizontalNo(siEHitHorizontalNo);
        analysis->SetSiEHitVerticalNo(siEHitVerticalNo);

	//NaI Gamma Detectors
	std::vector<G4double> naiEDep(16,0.0);
	if (hce && fHCID_NaI >= 0)
	{
		SiDetectorHitsCollection *hcNaI = (SiDetectorHitsCollection *)hce->GetHC(fHCID_NaI);
		if(hcNaI)
		{
			G4int nofHitsNaI = hcNaI->GetSize();
			for (G4int i =0; i < nofHitsNaI; i++)
			{
				G4int copyNo = (*hcNaI)[i]->GetCopyNo();
				
				if (copyNo >=0 && copyNo < 16)
				{
					naiEDep[copyNo] += (*hcNaI)[i] ->GetEnergyDeposit();
				}
			}
		}

	}
        analysis->SetNaIEDep(naiEDep);
        //
        // PPAC 1
        PpacHitsCollection *hcPpac1 = (PpacHitsCollection *)hce->GetHC(fHCID_Ppac1);
        G4int nofHitsPpac1 = hcPpac1->GetSize();
        G4bool ppac1HitAccepted;
        std::vector<G4int> ppac1HitTrackId;
        std::vector<G4double> ppac1HitEDep;
        std::vector<G4double> ppac1HitTime;
        std::vector<G4double> ppac1HitGlobalPosX;
        std::vector<G4double> ppac1HitGlobalPosY;
        std::vector<G4double> ppac1HitGlobalPosZ;
        std::vector<G4double> ppac1HitLocalPosX;
        std::vector<G4double> ppac1HitLocalPosY;
        std::vector<G4double> ppac1HitLocalPosZ;
        std::vector<G4double> ppac1HitGlobalMomentumX;
        std::vector<G4double> ppac1HitGlobalMomentumY;
        std::vector<G4double> ppac1HitGlobalMomentumZ;
        std::vector<G4double> ppac1HitLocalMomentumX;
        std::vector<G4double> ppac1HitLocalMomentumY;
        std::vector<G4double> ppac1HitLocalMomentumZ;
        if (nofHitsPpac1 == 0)
        {
            // G4ExceptionDescription msg;
            // msg << "No hits in PPAC1 of this event found.\n";
            // G4Exception("\n---> EventAction::EndOfEventAction()",
            //             "MdmSimCode01", JustWarning, msg);
            ppac1HitAccepted = false;
            ppac1HitTrackId.push_back(-99);
            ppac1HitEDep.push_back(-99.);
            ppac1HitTime.push_back(-99.);
            ppac1HitGlobalPosX.push_back(-999.);
            ppac1HitGlobalPosY.push_back(-999.);
            ppac1HitGlobalPosZ.push_back(-999.);
            ppac1HitLocalPosX.push_back(-999.);
            ppac1HitLocalPosY.push_back(-999.);
            ppac1HitLocalPosZ.push_back(-999.);
            ppac1HitGlobalMomentumX.push_back(-999.);
            ppac1HitGlobalMomentumY.push_back(-999.);
            ppac1HitGlobalMomentumZ.push_back(-999.);
            ppac1HitLocalMomentumX.push_back(-999.);
            ppac1HitLocalMomentumY.push_back(-999.);
            ppac1HitLocalMomentumZ.push_back(-999.);
        }
        else
        {
            ppac1HitAccepted = true;
            ppac1HitTrackId.push_back((*hcPpac1)[0]->GetTrackID());
            G4double ppac1eDep = 0.;
            G4double ppac1Time = 0.;
            G4ThreeVector ppac1GlobalPos(0., 0., 0.);
            G4ThreeVector ppac1LocalPos(0., 0., 0.);
            G4ThreeVector ppac1GlobalMomentum(0., 0., 0.);
            G4ThreeVector ppac1LocalMomentum(0., 0., 0.);
            G4int itTimes = 0;
            for (int i = 0; i < nofHitsPpac1; i++)
            {
                if ((*hcPpac1)[i]->GetTrackID() == ppac1HitTrackId.back())
                {
                    ppac1eDep += (*hcPpac1)[i]->GetEnergyDeposit();
                    ppac1Time += (*hcPpac1)[i]->GetTime();
                    ppac1GlobalPos += (*hcPpac1)[i]->GetGlobalPosition();
                    ppac1LocalPos += (*hcPpac1)[i]->GetLocalPosition();
                    ppac1GlobalMomentum += (*hcPpac1)[i]->GetGlobalMomentum();
                    ppac1LocalMomentum += (*hcPpac1)[i]->GetLocalMomentum();
                    itTimes += 1;
                }
                else
                {
                    ppac1HitTrackId.push_back((*hcPpac1)[i]->GetTrackID());
                    ppac1HitEDep.push_back(ppac1eDep);
                    ppac1HitTime.push_back(ppac1Time / (double)itTimes);
                    ppac1HitGlobalPosX.push_back(ppac1GlobalPos.x() / (double)itTimes);
                    ppac1HitGlobalPosY.push_back(ppac1GlobalPos.y() / (double)itTimes);
                    ppac1HitGlobalPosZ.push_back(ppac1GlobalPos.z() / (double)itTimes);
                    ppac1HitLocalPosX.push_back(ppac1LocalPos.x() / (double)itTimes);
                    ppac1HitLocalPosY.push_back(ppac1LocalPos.y() / (double)itTimes);
                    ppac1HitLocalPosZ.push_back(ppac1LocalPos.z() / (double)itTimes);
                    ppac1HitGlobalMomentumX.push_back(ppac1GlobalMomentum.x() / (double)itTimes);
                    ppac1HitGlobalMomentumY.push_back(ppac1GlobalMomentum.y() / (double)itTimes);
                    ppac1HitGlobalMomentumZ.push_back(ppac1GlobalMomentum.z() / (double)itTimes);
                    ppac1HitLocalMomentumX.push_back(ppac1LocalMomentum.x() / (double)itTimes);
                    ppac1HitLocalMomentumY.push_back(ppac1LocalMomentum.y() / (double)itTimes);
                    ppac1HitLocalMomentumZ.push_back(ppac1LocalMomentum.z() / (double)itTimes);
                    ppac1eDep = (*hcPpac1)[i]->GetEnergyDeposit();
                    ppac1Time = (*hcPpac1)[i]->GetTime();
                    ppac1GlobalPos = (*hcPpac1)[i]->GetGlobalPosition();
                    ppac1LocalPos = (*hcPpac1)[i]->GetLocalPosition();
                    ppac1GlobalMomentum = (*hcPpac1)[i]->GetGlobalMomentum();
                    ppac1LocalMomentum = (*hcPpac1)[i]->GetLocalMomentum();
                    itTimes = 1;
                }
            }
            ppac1HitEDep.push_back(ppac1eDep);
            ppac1HitTime.push_back(ppac1Time / (double)itTimes);
            ppac1HitGlobalPosX.push_back(ppac1GlobalPos.x() / (double)itTimes);
            ppac1HitGlobalPosY.push_back(ppac1GlobalPos.y() / (double)itTimes);
            ppac1HitGlobalPosZ.push_back(ppac1GlobalPos.z() / (double)itTimes);
            ppac1HitLocalPosX.push_back(ppac1LocalPos.x() / (double)itTimes);
            ppac1HitLocalPosY.push_back(ppac1LocalPos.y() / (double)itTimes);
            ppac1HitLocalPosZ.push_back(ppac1LocalPos.z() / (double)itTimes);
            ppac1HitGlobalMomentumX.push_back(ppac1GlobalMomentum.x() / (double)itTimes);
            ppac1HitGlobalMomentumY.push_back(ppac1GlobalMomentum.y() / (double)itTimes);
            ppac1HitGlobalMomentumZ.push_back(ppac1GlobalMomentum.z() / (double)itTimes);
            ppac1HitLocalMomentumX.push_back(ppac1LocalMomentum.x() / (double)itTimes);
            ppac1HitLocalMomentumY.push_back(ppac1LocalMomentum.y() / (double)itTimes);
            ppac1HitLocalMomentumZ.push_back(ppac1LocalMomentum.z() / (double)itTimes);
            // G4cout << "ppac1HitGlobalPosX=" << G4BestUnit(ppac1HitGlobalPosX[0], "Length") << G4endl;
            // G4cout << "ppac1HitGlobalPosY=" << G4BestUnit(ppac1HitGlobalPosY[0], "Length") << G4endl;
            // G4cout << "ppac1HitGlobalPosZ=" << G4BestUnit(ppac1HitGlobalPosZ[0], "Length") << G4endl;
            // G4cout << "ppac1HitLocalPosX=" << ppac1HitLocalPosX[0] << G4endl;
            // G4cout << "ppac1HitLocalPosY=" << ppac1HitLocalPosY[0] << G4endl;
        }
        analysis->SetPpac1Accepted(ppac1HitAccepted);
        analysis->SetPpac1HitTrackId(ppac1HitTrackId);
        analysis->SetPpac1HitEDep(ppac1HitEDep);
        analysis->SetPpac1HitTime(ppac1HitTime);
        analysis->SetPpac1HitGlobalPosX(ppac1HitGlobalPosX);
        analysis->SetPpac1HitGlobalPosY(ppac1HitGlobalPosY);
        analysis->SetPpac1HitGlobalPosZ(ppac1HitGlobalPosZ);
        analysis->SetPpac1HitLocalPosX(ppac1HitLocalPosX);
        analysis->SetPpac1HitLocalPosY(ppac1HitLocalPosY);
        analysis->SetPpac1HitLocalPosZ(ppac1HitLocalPosZ);
        analysis->SetPpac1HitGlobalMomentumX(ppac1HitGlobalMomentumX);
        analysis->SetPpac1HitGlobalMomentumY(ppac1HitGlobalMomentumY);
        analysis->SetPpac1HitGlobalMomentumZ(ppac1HitGlobalMomentumZ);
        analysis->SetPpac1HitLocalMomentumX(ppac1HitLocalMomentumX);
        analysis->SetPpac1HitLocalMomentumY(ppac1HitLocalMomentumY);
        analysis->SetPpac1HitLocalMomentumZ(ppac1HitLocalMomentumZ);

        //
        // PPAC 2
        PpacHitsCollection *hcPpac2 = (PpacHitsCollection *)hce->GetHC(fHCID_Ppac2);
        G4int nofHitsPpac2 = hcPpac2->GetSize();
        G4bool ppac2HitAccepted;
        std::vector<G4int> ppac2HitTrackId;
        std::vector<G4double> ppac2HitEDep;
        std::vector<G4double> ppac2HitTime;
        std::vector<G4double> ppac2HitGlobalPosX;
        std::vector<G4double> ppac2HitGlobalPosY;
        std::vector<G4double> ppac2HitGlobalPosZ;
        std::vector<G4double> ppac2HitLocalPosX;
        std::vector<G4double> ppac2HitLocalPosY;
        std::vector<G4double> ppac2HitLocalPosZ;
        std::vector<G4double> ppac2HitGlobalMomentumX;
        std::vector<G4double> ppac2HitGlobalMomentumY;
        std::vector<G4double> ppac2HitGlobalMomentumZ;
        std::vector<G4double> ppac2HitLocalMomentumX;
        std::vector<G4double> ppac2HitLocalMomentumY;
        std::vector<G4double> ppac2HitLocalMomentumZ;
        if (nofHitsPpac2 == 0)
        {
            // G4ExceptionDescription msg;
            // msg << "No hits in PPAC2 of this event found.\n";
            // G4Exception("\n---> EventAction::EndOfEventAction()",
            //             "MdmSimCode01", JustWarning, msg);
            ppac2HitAccepted = false;
            ppac2HitTrackId.push_back(-99);
            ppac2HitEDep.push_back(-99.);
            ppac2HitTime.push_back(-99.);
            ppac2HitGlobalPosX.push_back(-999.);
            ppac2HitGlobalPosY.push_back(-999.);
            ppac2HitGlobalPosZ.push_back(-999.);
            ppac2HitLocalPosX.push_back(-999.);
            ppac2HitLocalPosY.push_back(-999.);
            ppac2HitLocalPosZ.push_back(-999.);
            ppac2HitGlobalMomentumX.push_back(-999.);
            ppac2HitGlobalMomentumY.push_back(-999.);
            ppac2HitGlobalMomentumZ.push_back(-999.);
            ppac2HitLocalMomentumX.push_back(-999.);
            ppac2HitLocalMomentumY.push_back(-999.);
            ppac2HitLocalMomentumZ.push_back(-999.);
        }
        else
        {
            ppac2HitAccepted = true;
            ppac2HitTrackId.push_back((*hcPpac2)[0]->GetTrackID());
            G4double ppac2eDep = 0.;
            G4double ppac2Time = 0.;
            G4ThreeVector ppac2GlobalPos(0., 0., 0.);
            G4ThreeVector ppac2LocalPos(0., 0., 0.);
            G4ThreeVector ppac2GlobalMomentum(0., 0., 0.);
            G4ThreeVector ppac2LocalMomentum(0., 0., 0.);
            G4int itTimes = 0;
            for (int i = 0; i < nofHitsPpac2; i++)
            {
                if ((*hcPpac2)[i]->GetTrackID() == ppac2HitTrackId.back())
                {
                    ppac2eDep += (*hcPpac2)[i]->GetEnergyDeposit();
                    ppac2Time += (*hcPpac2)[i]->GetTime();
                    ppac2GlobalPos += (*hcPpac2)[i]->GetGlobalPosition();
                    ppac2LocalPos += (*hcPpac2)[i]->GetLocalPosition();
                    ppac2GlobalMomentum += (*hcPpac2)[i]->GetGlobalMomentum();
                    ppac2LocalMomentum += (*hcPpac2)[i]->GetLocalMomentum();
                    itTimes += 1;
                }
                else
                {
                    ppac2HitTrackId.push_back((*hcPpac2)[i]->GetTrackID());
                    ppac2HitEDep.push_back(ppac2eDep);
                    ppac2HitTime.push_back(ppac2Time / (double)itTimes);
                    ppac2HitGlobalPosX.push_back(ppac2GlobalPos.x() / (double)itTimes);
                    ppac2HitGlobalPosY.push_back(ppac2GlobalPos.y() / (double)itTimes);
                    ppac2HitGlobalPosZ.push_back(ppac2GlobalPos.z() / (double)itTimes);
                    ppac2HitLocalPosX.push_back(ppac2LocalPos.x() / (double)itTimes);
                    ppac2HitLocalPosY.push_back(ppac2LocalPos.y() / (double)itTimes);
                    ppac2HitLocalPosZ.push_back(ppac2LocalPos.z() / (double)itTimes);
                    ppac2HitGlobalMomentumX.push_back(ppac2GlobalMomentum.x() / (double)itTimes);
                    ppac2HitGlobalMomentumY.push_back(ppac2GlobalMomentum.y() / (double)itTimes);
                    ppac2HitGlobalMomentumZ.push_back(ppac2GlobalMomentum.z() / (double)itTimes);
                    ppac2HitLocalMomentumX.push_back(ppac2LocalMomentum.x() / (double)itTimes);
                    ppac2HitLocalMomentumY.push_back(ppac2LocalMomentum.y() / (double)itTimes);
                    ppac2HitLocalMomentumZ.push_back(ppac2LocalMomentum.z() / (double)itTimes);
                    ppac2eDep = (*hcPpac2)[i]->GetEnergyDeposit();
                    ppac2Time = (*hcPpac2)[i]->GetTime();
                    ppac2GlobalPos = (*hcPpac2)[i]->GetGlobalPosition();
                    ppac2LocalPos = (*hcPpac2)[i]->GetLocalPosition();
                    ppac2GlobalMomentum = (*hcPpac2)[i]->GetGlobalMomentum();
                    ppac2LocalMomentum = (*hcPpac2)[i]->GetLocalMomentum();
                    itTimes = 1;
                }
            }
            ppac2HitEDep.push_back(ppac2eDep);
            ppac2HitTime.push_back(ppac2Time / (double)itTimes);
            ppac2HitGlobalPosX.push_back(ppac2GlobalPos.x() / (double)itTimes);
            ppac2HitGlobalPosY.push_back(ppac2GlobalPos.y() / (double)itTimes);
            ppac2HitGlobalPosZ.push_back(ppac2GlobalPos.z() / (double)itTimes);
            ppac2HitLocalPosX.push_back(ppac2LocalPos.x() / (double)itTimes);
            ppac2HitLocalPosY.push_back(ppac2LocalPos.y() / (double)itTimes);
            ppac2HitLocalPosZ.push_back(ppac2LocalPos.z() / (double)itTimes);
            ppac2HitGlobalMomentumX.push_back(ppac2GlobalMomentum.x() / (double)itTimes);
            ppac2HitGlobalMomentumY.push_back(ppac2GlobalMomentum.y() / (double)itTimes);
            ppac2HitGlobalMomentumZ.push_back(ppac2GlobalMomentum.z() / (double)itTimes);
            ppac2HitLocalMomentumX.push_back(ppac2LocalMomentum.x() / (double)itTimes);
            ppac2HitLocalMomentumY.push_back(ppac2LocalMomentum.y() / (double)itTimes);
            ppac2HitLocalMomentumZ.push_back(ppac2LocalMomentum.z() / (double)itTimes);
            // G4cout << "ppac2HitGlobalPosX=" << G4BestUnit(ppac2HitGlobalPosX[0], "Length") << G4endl;
            // G4cout << "ppac2HitGlobalPosY=" << G4BestUnit(ppac2HitGlobalPosY[0], "Length") << G4endl;
            // G4cout << "ppac2HitGlobalPosZ=" << G4BestUnit(ppac2HitGlobalPosZ[0], "Length") << G4endl;
            // G4cout << "ppac2HitLocalPosX=" << G4BestUnit(ppac2HitLocalPosX[0], "Length") << G4endl;
            // G4cout << "ppac2HitLocalPosY=" << G4BestUnit(ppac2HitLocalPosY[0], "Length") << G4endl;
            // G4cout << "ppac2HitLocalPosZ=" << G4BestUnit(ppac2HitLocalPosZ[0], "Length") << G4endl;
        }
        analysis->SetPpac2Accepted(ppac2HitAccepted);
        analysis->SetPpac2HitTrackId(ppac2HitTrackId);
        analysis->SetPpac2HitEDep(ppac2HitEDep);
        analysis->SetPpac2HitTime(ppac2HitTime);
        analysis->SetPpac2HitGlobalPosX(ppac2HitGlobalPosX);
        analysis->SetPpac2HitGlobalPosY(ppac2HitGlobalPosY);
        analysis->SetPpac2HitGlobalPosZ(ppac2HitGlobalPosZ);
        analysis->SetPpac2HitLocalPosX(ppac2HitLocalPosX);
        analysis->SetPpac2HitLocalPosY(ppac2HitLocalPosY);
        analysis->SetPpac2HitLocalPosZ(ppac2HitLocalPosZ);
        analysis->SetPpac2HitGlobalMomentumX(ppac2HitGlobalMomentumX);
        analysis->SetPpac2HitGlobalMomentumY(ppac2HitGlobalMomentumY);
        analysis->SetPpac2HitGlobalMomentumZ(ppac2HitGlobalMomentumZ);
        analysis->SetPpac2HitLocalMomentumX(ppac2HitLocalMomentumX);
        analysis->SetPpac2HitLocalMomentumY(ppac2HitLocalMomentumY);
        analysis->SetPpac2HitLocalMomentumZ(ppac2HitLocalMomentumZ);

        //
        // Coincidence flags, matched to the tagged particles by track ID
        //   alpha1 @ MDM : alpha1 crosses the focal plane (AtPpacs: and both PPACs)
        //   alpha2 @ S2  : alpha2 deposits more than kS2ThresholdMeV in the S2 layers
        //   gamma @ NaI  : total NaI deposit above kNaIThresholdMeV in an event with
        //                  a decay gamma (charged products cannot reach the NaI crystals)
        constexpr G4double kS2ThresholdMeV = 0.1;
        constexpr G4double kNaIThresholdMeV = 0.1;
        const G4int alpha1TrackId = analysis->GetAlpha1TrackId();
        const G4int alpha2TrackId = analysis->GetAlpha2TrackId();
        const auto hasTrack = [](PpacHitsCollection *hc, G4int trackId)
        {
            if (!hc || trackId < 0)
                return false;
            for (size_t i = 0; i < hc->GetSize(); i++)
                if ((*hc)[i]->GetTrackID() == trackId)
                    return true;
            return false;
        };
        const auto s2EDep = [&](G4int trackId)
        {
            G4double eDep = 0.;
            for (G4int i = 0; trackId >= 0 && i < nofHitsSiDetectorE; i++)
                if ((*hcSiDetectorE)[i]->GetTrackID() == trackId)
                    eDep += (*hcSiDetectorE)[i]->GetEnergyDeposit();
            return eDep / MeV;
        };
        G4double naiEDepTotal = 0.;
        for (const G4double eDep : naiEDep)
            naiEDepTotal += eDep / MeV;

        const G4bool alpha1AtMdm = hasTrack(hcRaytraceFocalPlane, alpha1TrackId);
        const G4bool alpha1AtPpacs = hasTrack(hcPpac1, alpha1TrackId) && hasTrack(hcPpac2, alpha1TrackId);
        const G4double alpha2S2EDep = s2EDep(alpha2TrackId);
        const G4bool alpha2AtS2 = alpha2S2EDep > kS2ThresholdMeV;
        const G4bool gammaAtNaI = analysis->GetGammaGenerated() && naiEDepTotal > kNaIThresholdMeV;
        analysis->SetAlpha1AtMdm(alpha1AtMdm);
        analysis->SetAlpha1AtPpacs(alpha1AtPpacs);
        analysis->SetAlpha1S2EDep(s2EDep(alpha1TrackId));
        analysis->SetAlpha2AtS2(alpha2AtS2);
        analysis->SetAlpha2S2EDep(alpha2S2EDep);
        analysis->SetGammaAtNaI(gammaAtNaI);
        analysis->SetNaIEDepTotal(naiEDepTotal);
        analysis->SetTripleCoincidence(alpha1AtMdm && alpha2AtS2 && gammaAtNaI);

        //
        // Fill the trees
        analysis->FillTreeAccurate();
        analysis->FillTreeExperiment();
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    void EventAction::ParseParams(std::map<std::string, G4double> params)
    {
        G4cout << "\n---> EventAction::ParseParams():" << G4endl;
        for (auto it : params)
        {
            if (it.first == "SiDetectorEnergyResolutionFwhmFraction")
            {
                fSiDetectorEnergyResolutionFwhmFraction = it.second;
                printf("Set: Si detector energy resolution FWHM = %.2f%%\n", fSiDetectorEnergyResolutionFwhmFraction * 100.);
            }
            else if (it.first == "TdcResolutionHalfWidthInNs")
            {
                fTdcResolutionHalfWidth = it.second * ns;
                G4cout << "Set: TDC resolution half-width = " << G4BestUnit(fTdcResolutionHalfWidth, "Time") << G4endl;
            }
            else if (it.first == "FirstMultipoleProbe")
            {
                fFirstMultipoleProbe = it.second;
                printf("Set: First multipole probe = %.4f\n", fFirstMultipoleProbe);
            }
            else if (it.first == "DipoleProbe")
            {
                fDipoleProbe = it.second;
                printf("Set: Dipole probe = %.4f\n", fDipoleProbe);
            }
            else if (it.first == "BeamCharge")
            {
                fBeamCharge = static_cast<G4int>(std::lround(it.second));
                printf("Set: Primary charge state = %d e; Raytrace uses each slit hit charge\n", fBeamCharge);
            }
            else if (it.first == "MdmAngleInDeg")
            {
                fMdmAngle = it.second;
                printf("Set: Raytrace angle = %.4f deg\n", fMdmAngle);
            }
        }
    }
}
