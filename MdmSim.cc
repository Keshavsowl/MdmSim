#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"
#include "ReactionPhysics.hh"
#include "nlohmann/json.hpp"

#include "G4RunManagerFactory.hh"
#include "G4SteppingVerbose.hh"
#include "G4UImanager.hh"
#include "G4VModularPhysicsList.hh"
#include "QBBC.hh"

#ifdef MDMSIM_USE_UIVIS
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#endif

#include "Randomize.hh"

#include <array>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <map>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>

using namespace MdmSim;
using json = nlohmann::json;

namespace
{
    class ConfigurationError : public std::runtime_error
    {
    public:
        using std::runtime_error::runtime_error;
    };

    std::filesystem::path ResolvePath(const std::filesystem::path &base, const std::filesystem::path &path)
    {
        if (path.is_absolute())
        {
            return path.lexically_normal();
        }
        return (base / path).lexically_normal();
    }

    std::string ResolveMapPath(const json &config,
                               const char *key,
                               const std::filesystem::path &mapDirectory)
    {
        const std::filesystem::path configured = config.at(key).get<std::string>();
        return ResolvePath(mapDirectory, configured).string();
    }

    void RequireFile(const std::filesystem::path &path, const std::string &description)
    {
        if (!std::filesystem::is_regular_file(path))
        {
            throw ConfigurationError(description + " does not exist: " + path.string());
        }
    }

    void ValidateNuclide(G4int z, G4int a, const std::string &name)
    {
        if (z < 0 || a <= 0 || a < z)
        {
            throw ConfigurationError(name + " must satisfy Z >= 0, A > 0, and A >= Z.");
        }
    }

    void ValidateConfiguration(const json &config)
    {
        if (!config.is_object())
        {
            throw ConfigurationError("the top-level JSON value must be an object.");
        }

        const std::array<std::pair<const char *, const char *>, 4> removedKeys{{
            {"MarcoName", "MacroName"},
            {"TdcResolutionInNs", "TdcResolutionHalfWidthInNs"},
            {"SiDetectorEnergyResolution", "SiDetectorEnergyResolutionFwhmFraction"},
            {"PpacVacuumInTorr", "PpacGasPressureInTorr"},
        }};
        for (const auto &[oldName, newName] : removedKeys)
        {
            if (config.contains(oldName))
            {
                throw ConfigurationError(std::string("\"") + oldName + "\" is no longer supported.\nUse \"" + newName + "\".");
            }
        }

        const std::array<const char *, 12> integerKeys{{
            "ProcessNumber",
            "BeamZ",
            "BeamA",
            "BeamCharge",
            "ReactionTargetZ",
            "ReactionTargetA",
            "ReactionLightZ",
            "ReactionLightA",
            "ReactionLightCharge",
            "ReactionHeavyZ",
            "ReactionHeavyA",
            "ReactionHeavyCharge",
        }};
        for (const char *key : integerKeys)
        {
            if (!config.at(key).is_number_integer())
            {
                throw ConfigurationError(std::string(key) + " must be an integer.");
            }
        }

        const G4int beamZ = config.at("BeamZ").get<G4int>();
        const G4int beamA = config.at("BeamA").get<G4int>();
        const G4int beamCharge = config.at("BeamCharge").get<G4int>();
        if (beamZ <= 0)
            throw ConfigurationError("BeamZ must be positive.");
        if (beamA <= 0)
            throw ConfigurationError("BeamA must be positive.");
        if (beamA < beamZ)
            throw ConfigurationError("BeamA must be greater than or equal to BeamZ.");
        if (beamCharge < 0 || beamCharge > beamZ)
            throw ConfigurationError("BeamCharge must satisfy 0 <= BeamCharge <= BeamZ.");
        if (config.at("BeamEnergyInMeV").get<G4double>() < 0.)
            throw ConfigurationError("BeamEnergyInMeV must be nonnegative.\nUnit: MeV.");

        const G4double reactionProbability = config.at("ReactionProbability").get<G4double>();
        if (reactionProbability < 0. || reactionProbability > 1.)
            throw ConfigurationError("ReactionProbability must satisfy 0 <= ReactionProbability <= 1.");

        const G4int targetZ = config.at("ReactionTargetZ").get<G4int>();
        const G4int targetA = config.at("ReactionTargetA").get<G4int>();
        const G4int lightZ = config.at("ReactionLightZ").get<G4int>();
        const G4int lightA = config.at("ReactionLightA").get<G4int>();
        const G4int lightCharge = config.at("ReactionLightCharge").get<G4int>();
        const G4int heavyZ = config.at("ReactionHeavyZ").get<G4int>();
        const G4int heavyA = config.at("ReactionHeavyA").get<G4int>();
        const G4int heavyCharge = config.at("ReactionHeavyCharge").get<G4int>();
        ValidateNuclide(targetZ, targetA, "Reaction target");
        ValidateNuclide(lightZ, lightA, "Reaction light product");
        ValidateNuclide(heavyZ, heavyA, "Reaction heavy product");
        if (lightCharge < 0 || lightCharge > lightZ)
            throw ConfigurationError("ReactionLightCharge must satisfy 0 <= charge <= ReactionLightZ.");
        if (heavyCharge < 0 || heavyCharge > heavyZ)
            throw ConfigurationError("ReactionHeavyCharge must satisfy 0 <= charge <= ReactionHeavyZ.");
        if (config.at("ReactionLightExMeV").get<G4double>() < 0. ||
            config.at("ReactionHeavyExMeV").get<G4double>() < 0.)
            throw ConfigurationError("Reaction excitation energies must be nonnegative.\nUnit: MeV.");

        if (config.at("TargetThicknessInMgCm2").get<G4double>() < 0.)
            throw ConfigurationError("TargetThicknessInMgCm2 must be nonnegative.\nUnit: mg/cm^2.");
        if (config.at("TdcResolutionHalfWidthInNs").get<G4double>() < 0.)
            throw ConfigurationError("TdcResolutionHalfWidthInNs must be nonnegative.\nUnit: ns.");
        if (config.at("SiDetectorEnergyResolutionFwhmFraction").get<G4double>() < 0.)
            throw ConfigurationError("SiDetectorEnergyResolutionFwhmFraction must be nonnegative.");
        if (config.at("PpacGasPressureInTorr").get<G4double>() < 0.)
            throw ConfigurationError("PpacGasPressureInTorr must be nonnegative.\nUnit: Torr.");
        if (config.at("PpacLengthInCm").get<G4double>() <= 0.)
            throw ConfigurationError("PpacLengthInCm must be positive.\nUnit: cm.");
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char **argv)
{
    if (argc != 2)
    {
        G4cerr << "Usage: ./MdmSim <config.json>" << G4endl;
        return 1;
    }

    try
    {
        const std::filesystem::path configFile = argv[1];
        std::ifstream configStream(configFile);
        if (!configStream)
        {
            throw ConfigurationError("cannot open configuration file: " + configFile.string());
        }
        const json config = json::parse(configStream, nullptr, true, true);
        ValidateConfiguration(config);

        const std::filesystem::path configDirectory = std::filesystem::absolute(std::filesystem::path(configFile)).parent_path();
        const std::filesystem::path fieldMapDirectory = ResolvePath(configDirectory, config.at("FieldMapDirectory").get<std::string>());

        const G4bool interactive = config.at("Interactive").get<G4bool>();
        const G4String macroName = config.at("MacroName").get<std::string>();
        const G4int processNumber = config.at("ProcessNumber").get<G4int>();
        const G4bool reactionEnabled = config.at("ReactionEnabled").get<G4bool>();
        const G4String startupMacro = interactive ? "init_vis.mac" : macroName;
        RequireFile(ResolvePath(std::filesystem::current_path(), std::string(startupMacro)), "Macro file");

        std::map<std::string, G4double> primaryParameters; // parameters for PrimaryGeneratorAction
        primaryParameters["BeamZ"] = config.at("BeamZ").get<G4double>();
        primaryParameters["BeamA"] = config.at("BeamA").get<G4double>();
        primaryParameters["BeamCharge"] = config.at("BeamCharge").get<G4double>();
        primaryParameters["BeamEnergyInMeV"] = config.at("BeamEnergyInMeV").get<G4double>();
        primaryParameters["MdmAngleInDeg"] = config.at("MdmAngleInDeg").get<G4double>();

        std::map<std::string, G4double> detectorParameters; // parameters for DectectorConstruction
        detectorParameters["TargetRotationAngleInDeg"] = config.at("TargetRotationAngleInDeg").get<G4double>();
        detectorParameters["TargetThicknessInMgCm2"] = config.at("TargetThicknessInMgCm2").get<G4double>();
        detectorParameters["UseDeltaE"] = config.at("UseDeltaE").get<G4bool>();
        detectorParameters["SiDetectorAngleInDeg"] = config.at("SiDetectorAngleInDeg").get<G4double>();
        detectorParameters["MdmAngleInDeg"] = config.at("MdmAngleInDeg").get<G4double>();
        detectorParameters["BeamZ"] = primaryParameters["BeamZ"];
        detectorParameters["BeamA"] = primaryParameters["BeamA"];
        detectorParameters["BeamCharge"] = primaryParameters["BeamCharge"];
        detectorParameters["ReactionEnabled"] = reactionEnabled ? 1.0 : 0.0;
        detectorParameters["FirstMultipoleProbe"] = config.at("FirstMultipoleProbe").get<G4double>();
        detectorParameters["DipoleProbe"] = config.at("DipoleProbe").get<G4double>();
        detectorParameters["PpacGasPressureInTorr"] = config.at("PpacGasPressureInTorr").get<G4double>();
        detectorParameters["PpacLengthInCm"] = config.at("PpacLengthInCm").get<G4double>();

        MdmFieldMapPaths fieldMapPaths;
        fieldMapPaths.multipole = ResolveMapPath(config, "MultipoleFieldMap", fieldMapDirectory);
        fieldMapPaths.dipoleEntrance = ResolveMapPath(config, "DipoleEntranceFieldMap", fieldMapDirectory);
        fieldMapPaths.dipoleSector = ResolveMapPath(config, "DipoleSectorFieldMap", fieldMapDirectory);
        fieldMapPaths.dipoleExit = ResolveMapPath(config, "DipoleExitFieldMap", fieldMapDirectory);
        RequireFile(fieldMapPaths.multipole, "Multipole field map");
        RequireFile(fieldMapPaths.dipoleEntrance, "Dipole entrance field map");
        RequireFile(fieldMapPaths.dipoleSector, "Dipole sector field map");
        RequireFile(fieldMapPaths.dipoleExit, "Dipole exit field map");

        std::map<std::string, G4double> eventParameters; // parameters for EventAction
        eventParameters["SiDetectorEnergyResolutionFwhmFraction"] = config.at("SiDetectorEnergyResolutionFwhmFraction").get<G4double>();
        eventParameters["TdcResolutionHalfWidthInNs"] = config.at("TdcResolutionHalfWidthInNs").get<G4double>();
        eventParameters["MdmAngleInDeg"] = config.at("MdmAngleInDeg").get<G4double>();
        eventParameters["FirstMultipoleProbe"] = config.at("FirstMultipoleProbe").get<G4double>();
        eventParameters["DipoleProbe"] = config.at("DipoleProbe").get<G4double>();
        eventParameters["BeamCharge"] = primaryParameters["BeamCharge"];

        std::map<std::string, G4double> reactionParameters;
        if (reactionEnabled)
        {
            reactionParameters["BeamZ"] = primaryParameters["BeamZ"];
            reactionParameters["BeamA"] = primaryParameters["BeamA"];
            reactionParameters["BeamCharge"] = primaryParameters["BeamCharge"];
            reactionParameters["ReactionProbability"] = config.at("ReactionProbability").get<G4double>();
            reactionParameters["TargetZ"] = config.at("ReactionTargetZ").get<G4double>();
            reactionParameters["TargetA"] = config.at("ReactionTargetA").get<G4double>();
            reactionParameters["LightProductZ"] = config.at("ReactionLightZ").get<G4double>();
            reactionParameters["LightProductA"] = config.at("ReactionLightA").get<G4double>();
            reactionParameters["LightProductCharge"] = config.at("ReactionLightCharge").get<G4double>();
            reactionParameters["LightProductExMeV"] = config.at("ReactionLightExMeV").get<G4double>();
            reactionParameters["HeavyProductZ"] = config.at("ReactionHeavyZ").get<G4double>();
            reactionParameters["HeavyProductA"] = config.at("ReactionHeavyA").get<G4double>();
            reactionParameters["HeavyProductCharge"] = config.at("ReactionHeavyCharge").get<G4double>();
            reactionParameters["HeavyProductExMeV"] = config.at("ReactionHeavyExMeV").get<G4double>();

            detectorParameters["ReactionLightZ"] = reactionParameters["LightProductZ"];
            detectorParameters["ReactionLightA"] = reactionParameters["LightProductA"];
            detectorParameters["ReactionLightCharge"] = reactionParameters["LightProductCharge"];
            detectorParameters["ReactionLightExMeV"] = reactionParameters["LightProductExMeV"];
            detectorParameters["ReactionHeavyZ"] = reactionParameters["HeavyProductZ"];
            detectorParameters["ReactionHeavyA"] = reactionParameters["HeavyProductA"];
            detectorParameters["ReactionHeavyCharge"] = reactionParameters["HeavyProductCharge"];
            detectorParameters["ReactionHeavyExMeV"] = reactionParameters["HeavyProductExMeV"];
        }

        std::map<std::string, std::variant<G4int, std::map<std::string, G4double>>> actionInitParameters; // parameters for ActionInitialization
        actionInitParameters["ProcessNumber"] = processNumber;
        actionInitParameters["PrimaryParameters"] = primaryParameters;
        actionInitParameters["EventParameters"] = eventParameters;

        G4Random::setTheEngine(new CLHEP::RanecuEngine);
        const std::uint64_t timeValue = static_cast<std::uint64_t>(
            std::chrono::high_resolution_clock::now().time_since_epoch().count());
        std::random_device randomDevice;
        std::seed_seq seedSequence{
            static_cast<std::uint32_t>(randomDevice()),
            static_cast<std::uint32_t>(randomDevice()),
            static_cast<std::uint32_t>(timeValue),
            static_cast<std::uint32_t>(timeValue >> 32),
            static_cast<std::uint32_t>(processNumber)};
        std::array<std::uint32_t, 2> generatedSeeds{};
        seedSequence.generate(generatedSeeds.begin(), generatedSeeds.end());
        std::array<G4long, 3> seeds{
            1 + static_cast<G4long>(generatedSeeds[0] % 2147483562U),
            1 + static_cast<G4long>(generatedSeeds[1] % 2147483398U),
            0};
        G4Random::setTheSeeds(seeds.data());
        G4cout << "Automatic random seeds: " << seeds[0] << " " << seeds[1] << G4endl;

        // use G4SteppingVerboseWithUnits
        G4int precision = 4;
        G4SteppingVerbose::UseBestUnit(precision);

        // Construct the default run manager
        //
        auto *runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::SerialOnly);

        // Set mandatory initialization classes
        //
        // Detector construction
        DetectorConstruction *detector = new DetectorConstruction();
        detector->SetFieldMapPaths(fieldMapPaths);
        detector->ParseParams(detectorParameters);
        runManager->SetUserInitialization(detector);

        // Physics list
        G4VModularPhysicsList *physicsList = new QBBC;
        physicsList->SetVerboseLevel(0);
        if (reactionEnabled)
        {
            auto *reactionPhysics = new ReactionPhysics;
            reactionPhysics->SetReactionParams(reactionParameters);
            physicsList->RegisterPhysics(reactionPhysics);
        }
        runManager->SetUserInitialization(physicsList);

        // User action initialization
        ActionInitialization *action = new ActionInitialization();
        action->SetParameters(actionInitParameters);
        runManager->SetUserInitialization(action);

#ifdef MDMSIM_USE_UIVIS
        G4VisManager *visManager = new G4VisExecutive;
        visManager->Initialize();
#endif

        // Get the pointer to the User Interface manager
        G4UImanager *UImanager = G4UImanager::GetUIpointer();

        // Process macro or start UI session
        //
        if (!interactive)
        {
            // batch mode
            G4String command = "/control/execute ";
            G4String fileName = macroName;
            if (UImanager->ApplyCommand(command + fileName) != 0)
            {
                throw std::runtime_error("Geant4 macro execution failed: " + std::string(fileName));
            }
        }
        else
        {
#ifdef MDMSIM_USE_UIVIS
            // interactive mode
            G4UIExecutive *ui = new G4UIExecutive(argc, argv);
            UImanager->ApplyCommand("/control/execute init_vis.mac");
            ui->SessionStart();
            delete ui;
#else
            throw std::runtime_error("Interactive mode requires WITH_GEANT4_UIVIS=ON");
#endif
        }

        // Job termination
        // Free the store: user actions, physics_list and detector_description are
        // owned and deleted by the run manager, so they should not be deleted
        // in the main() program !

#ifdef MDMSIM_USE_UIVIS
        delete visManager;
#endif
        delete runManager;
        return 0;
    }
    catch (const ConfigurationError &error)
    {
        G4cerr << "Configuration error:\n" << error.what() << G4endl;
        return 1;
    }
    catch (const json::exception &error)
    {
        G4cerr << "Configuration error:\n" << error.what() << G4endl;
        return 1;
    }
    catch (const std::exception &error)
    {
        G4cerr << "MdmSim error:\n" << error.what() << G4endl;
        return 1;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
