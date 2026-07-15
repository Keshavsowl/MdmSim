# MdmSim

MdmSim is an independent Geant4 program for transporting ions from the target
chamber through the MDM spectrometer to the focal-plane detectors. It has its
own `config/config.json`; it does not read an MDMTrace JSON configuration.

The simulation uses four field-map files with embedded text headers followed by
binary field arrays:

```text
field/Multipole.bin
field/DipoleEntrance.bin
field/DipoleSector.bin
field/DipoleExit.bin
```

These files are required before a run. They may be generated with
[MDMTrace](https://github.com/luozf14/MDMTrace), but MdmSim reads only its own
configuration and the four map files.

## Requirements and build

MdmSim requires CMake 3.16 or newer, a C++17 compiler, Geant4 11, ROOT 6, and a
Fortran compiler. Build from the repository root with:

```sh
cmake -S . -B build
cmake --build build -j4
```

Geant4 UI and visualization support is enabled by default. Use
`-DWITH_GEANT4_UIVIS=OFF` for a batch-only build.

CMake copies `init_vis.mac`, `run.mac`, `vis.mac`, and `rayin.dat` into
`build/`.

## Configuration

The configuration remains flat and supports both `//` line comments and
`/* ... */` block comments. The complete checked-in configuration is:

```jsonc
{
  // Run mode.
  "Interactive": false,

  // Geant4 macro used in batch mode.
  // This path is interpreted from the working directory, normally build/,
  // where CMake copies run.mac.
  "MacroName": "run.mac",

  // Job index used in the output filename and in automatic seed mixing.
  "ProcessNumber": 0,

  // Primary beam.
  "BeamZ": 6,
  "BeamA": 12,
  "BeamCharge": 5,                 // charge state Q, in units of e
  "BeamEnergyInMeV": 15.0,         // MeV kinetic energy

  // Forced relativistic two-body reaction.
  "ReactionEnabled": true,
  "ReactionProbability": 1.0,      // range: 0 to 1

  "ReactionTargetZ": 1,
  "ReactionTargetA": 1,

  "ReactionLightZ": 1,
  "ReactionLightA": 1,
  "ReactionLightCharge": 1,        // charge state, in units of e
  "ReactionLightExMeV": 0.0,       // MeV excitation energy

  "ReactionHeavyZ": 6,
  "ReactionHeavyA": 12,
  "ReactionHeavyCharge": 5,        // charge state, in units of e
  "ReactionHeavyExMeV": 0.0,       // MeV excitation energy

  // Virtual target geometry.
  // The current target material is G4_Galactic.
  // TargetThicknessInMgCm2 is converted to a geometrical thickness
  // using 2.253 g/cm^3; real target energy loss is not modeled.
  "TargetRotationAngleInDeg": 0.0, // deg
  "TargetThicknessInMgCm2": 159.0, // mg/cm^2

  // Uniform time smearing:
  // measured time = true time + Uniform(-halfWidth, +halfWidth).
  "TdcResolutionHalfWidthInNs": 0.25, // ns

  "UseDeltaE": false,
  "SiDetectorAngleInDeg": 90.0,    // deg

  // Fractional FWHM of the Gaussian silicon-energy resolution.
  // 0.10 means 10% FWHM; sigma = 0.10 * E / 2.355.
  "SiDetectorEnergyResolutionFwhmFraction": 0.10,

  // Rotation of the MDM geometry around the target.
  // The primary beam itself is generated along +z.
  "MdmAngleInDeg": 0.0,            // deg

  // Requested hall-probe settings. These must agree with the metadata
  // embedded in all four field maps.
  "DipoleProbe": 2338.00,          // Gauss
  "FirstMultipoleProbe": 1659.98,  // Gauss

  // FieldMapDirectory is resolved relative to config/config.json.
  "FieldMapDirectory": "../field",
  "MultipoleFieldMap": "Multipole.bin",
  "DipoleEntranceFieldMap": "DipoleEntrance.bin",
  "DipoleSectorFieldMap": "DipoleSector.bin",
  "DipoleExitFieldMap": "DipoleExit.bin",

  // PPAC working-gas pressure and cathode spacing.
  "PpacGasPressureInTorr": 4.0,     // Torr
  "PpacLengthInCm": 40.0            // cm
}
```

All listed values are required. `MacroName` is interpreted from the process
working directory, normally `build/`. `FieldMapDirectory` is interpreted from
the directory containing the selected configuration file; each field-map name
is then interpreted from that directory.

The primary is an ion with the configured `Z`, `A`, charge state `Q`, and total
kinetic energy. It is generated as a pencil beam along global `+z`.
`MdmAngleInDeg` rotates the MDM geometry around the target; it does not steer
the primary beam.

Detector response uses the following exact models:

```text
measured time = true time
              + Uniform(-TdcResolutionHalfWidthInNs,
                        +TdcResolutionHalfWidthInNs)

silicon sigma = SiDetectorEnergyResolutionFwhmFraction * E / 2.355
```

The silicon setting is a fractional Gaussian FWHM, so `0.10` means 10% FWHM.
`PpacGasPressureInTorr` is the pentane working-gas pressure; renaming it does not
change the existing gas-density calculation.

At startup MdmSim validates physical ranges, opens the macro and four map files,
and rejects missing required values. The program creates a fresh automatic
Ranecu seed state for each invocation by mixing `std::random_device`, a
high-resolution timestamp, and `ProcessNumber`. The two positive engine seeds
are printed once at startup. There is no user seed setting.

## Virtual target and reaction model

The target geometry is a virtual reaction volume made from `G4_Galactic`.
`TargetThicknessInMgCm2` is converted to geometrical length using
`2.253 g/cm3`, but the volume is not a real carbon target and real target
material energy loss is not modeled there.

When `ReactionEnabled` is true, the reaction process:

- considers only the configured primary beam ion;
- selects that primary with `ReactionProbability`;
- chooses a uniformly random reaction fraction along its geometrical target
  path;
- treats the target nucleus as initially at rest;
- obtains particle masses from Geant4 definitions;
- applies relativistic two-body kinematics;
- samples `cos(theta_CM)` uniformly and samples `phi` uniformly;
- throws an error if the configured channel is energetically closed;
- stops the primary and creates the light and heavy secondaries; and
- assigns the configured charge state to each secondary.

The configured charge state `Q` remains fixed. Material energy loss may change
momentum `p`, so magnetic rigidity `B rho = p / Q` may change.

## Field maps and Geant4 integration

Each map header states its map role, coordinate system, grid, units, and the two
hall-probe settings. MdmSim checks the expected role and coordinate-system name,
checks all four maps against `DipoleProbe` and `FirstMultipoleProbe`, and fails
on an incomplete binary payload. Probe comparisons use:

```text
abs(a - b) <= 2.0e-6 * max(1, abs(a), abs(b))
```

The example probe values `2338.00` and `1659.98` match the current MDMTrace
example maps. Those maps request 2.5 mm sampling, stored as 0.25 cm grid spacing
in their headers. MdmSim reads the grid spacing from each map; it does not read
MDMTrace's configuration. If maps are regenerated at different probe settings,
update the two values in MdmSim's own configuration to match their headers.

Map positions are in cm and map fields are in Tesla. Fields are evaluated with
trilinear interpolation and are zero outside the covered map region. Dipole
dispatch checks the entrance map, then the exit map, then the sector map. The
dipole convention uses a nominal angle of 100 degrees and a reference radius of
160 cm.

RAYTRACE uses exactly `c = 3.0e10 cm/s`, so the wrapper applies

```text
(300 mm/ns) / c_light
```

before passing the field to Geant4. Field integration uses `G4ClassicalRK4`, a
0.5 mm chord-finder minimum step, `DeltaOneStep = 1e-3 mm`, and a 1.0 mm user
maximum step in both field volumes.

## Run

Run from the build directory so the default macro path resolves correctly:

```sh
cd build
./MdmSim ../config/config.json
```

In batch mode, `/run/beamOn` in the file selected by `MacroName` controls the
event count.

## ROOT output

For `ProcessNumber = 0`, MdmSim writes `SimData_0.root`. In general the name is:

```text
SimData_<ProcessNumber>.root
```

The `AccurateData` tree contains detailed simulation values. The
`ExperimentalData` tree contains smeared silicon values and `DaqTrigger`, whose
current exact definition is:

```text
DaqTrigger = SiEHitAccepted
```

Important Boolean branches have deliberately narrow meanings:

- `SlitHitPresent` is true when the Geant4 slit hit collection is nonempty.
- `SlitHitTransmitted` is true when at least one embedded legacy RAYTRACE
  prediction made from a slit hit satisfies `|x| < 200 mm` and `|y| < 50 mm`.
  Despite its retained name, it is not a direct Geant4 statement that a
  particle was transmitted through the slit.
- `RaytraceFocalPlaneAccepted` is true when the Geant4 vacuum focal-plane scorer
  has a hit. It is independent of `SlitHitTransmitted`.

Vector branches use sentinel entries such as `-99` and `-999` when no hit is
present. Treat the corresponding acceptance/presence Boolean as authoritative
before consuming those vectors.

ROOT branch units are:

| Quantity | Unit |
| --- | --- |
| Global and local position | mm |
| Hit time | ns |
| Deposited energy | MeV |
| Global and local momentum | MeV/c |
| Embedded RAYTRACE position | mm |
| Embedded RAYTRACE angle | degrees |

## RAYTRACE comparison macro

The embedded legacy RAYTRACE calculation can be compared with Geant4's vacuum
focal-plane scorer or with a projection through the two PPAC planes:

```sh
root -l -b -q 'macros/CompareRaytraceGeant4.C("build/SimData_0.root")'
```

The default output is `raytrace_geant4_compare.root`, with PNG plots for `X`,
`Y`, `AngX`, and `AngY`. Positions are in mm, angles are in degrees, and the
residual convention is RAYTRACE minus Geant4. The macro requires
`SlitHitTransmitted` and distinguishes that legacy window test from the
Geant4 `RaytraceFocalPlaneAccepted` scorer.
