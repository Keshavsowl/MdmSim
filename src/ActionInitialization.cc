#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"
#include "TrackingAction.hh"

namespace MdmSim
{

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    ActionInitialization::ActionInitialization()
    {
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    ActionInitialization::~ActionInitialization()
    {
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void ActionInitialization::BuildForMaster() const
    {
        RunAction *runAction = new RunAction;
        runAction->SetProcessNumber(std::get<G4int>(fParameters.find("ProcessNumber")->second));
        SetUserAction(runAction);
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void ActionInitialization::Build() const
    {
        PrimaryGeneratorAction *primaryGeneratorAction = new PrimaryGeneratorAction;
        auto primaryParameters = fParameters.find("PrimaryParameters");
        if (primaryParameters != fParameters.end())
            primaryGeneratorAction->ParseParams(std::get<std::map<std::string, G4double>>(primaryParameters->second));
        SetUserAction(primaryGeneratorAction);

        RunAction *runAction = new RunAction;
        runAction->SetProcessNumber(std::get<G4int>(fParameters.find("ProcessNumber")->second));
        SetUserAction(runAction);

        EventAction *eventAction = new EventAction(runAction);
        eventAction->ParseParams(std::get<std::map<std::string, G4double>>(fParameters.find("EventParameters")->second));
        SetUserAction(eventAction);

        SetUserAction(new SteppingAction(eventAction));
        SetUserAction(new TrackingAction);
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
