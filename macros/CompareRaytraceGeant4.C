// Compare Raytrace focal-plane output to the Geant4 track.
//
// Usage from the repository root:
//   root -l -b -q 'macros/CompareRaytraceGeant4.C("build/SimData_99.root")'
//
// Source-selection examples:
//   // Use the vacuum Raytrace focal-plane scorer.
//   root -l -b -q 'macros/CompareRaytraceGeant4.C("SimData_93.root","compare_focal.root",-600,-489.99875,-89.99875,"AccurateData","raytrace")'
//
//   // Use projection from PPAC1/PPAC2 hit positions.
//   root -l -b -q 'macros/CompareRaytraceGeant4.C("SimData_93.root","compare_ppac.root",-600,-489.99875,-89.99875,"AccurateData","ppac")'
//
// The Raytrace plane is expressed in the PPAC chamber local z coordinate, in mm.
// In the current geometry the PPAC chamber center is placed 60 cm downstream of
// the Raytrace focal plane, so the default comparison plane is z = -600 mm.
//
// New MdmSim files store a vacuum RaytraceFocalPlane hit at the Raytrace focal
// plane. By default the macro uses that Geant4 branch, avoiding PPAC
// gas/window/cathode scattering in the comparison. The final macro argument can
// choose either "raytrace" focal-plane hits or "ppac" projection.
//
// The stored PPAC local z values are local to each cathode volume, not to the
// PPAC chamber. The PPAC mode therefore uses explicit chamber-local PPAC plane
// z values. Current defaults:
//   PPAC1 cathode: -500 mm + 2.5 um / 2 + 10 mm = -489.99875 mm
//   PPAC2 cathode: PPAC1 + 400 mm = -89.99875 mm
//
// PPAC local X is chamber-centered. In the current detector placement the
// chamber center lies on the Raytrace D-axis, so no empirical X offset is applied.
// The DIPO "Drawing Exit" coefficients in rayin.dat are EFB curvature terms
// (S12..S18), not coordinate-origin shifts. If a future deck uses nonzero
// DIPO XCR2, handle that as geometry, not as a fitted comparison offset.

#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLine.h>
#include <TMath.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TTree.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace {

template <typename T>
bool HasValue(const std::vector<T>* values) {
  return values && !values->empty();
}

double FirstOrNaN(const std::vector<double>* values) {
  return HasValue(values) ? values->front() : std::numeric_limits<double>::quiet_NaN();
}

double AtOrNaN(const std::vector<double>* values, std::size_t index) {
  return values && index < values->size() ? values->at(index) : std::numeric_limits<double>::quiet_NaN();
}

int FirstTrackOr(const std::vector<int>* values, int defaultTrackId) {
  return values && !values->empty() ? values->front() : defaultTrackId;
}

std::size_t IndexForTrack(const std::vector<int>* trackIds, int trackId) {
  if (!trackIds) {
    return 0;
  }
  for (std::size_t i = 0; i < trackIds->size(); ++i) {
    if (trackIds->at(i) == trackId) {
      return i;
    }
  }
  return 0;
}

struct CompareRow {
  double raytraceX = 0.0;
  double raytraceY = 0.0;
  double raytraceAngleX = 0.0;
  double raytraceAngleY = 0.0;
  double geant4X = 0.0;
  double geant4Y = 0.0;
  double geant4AngleX = 0.0;
  double geant4AngleY = 0.0;
};

struct Quantity {
  const char* key;
  const char* title;
  const char* unit;
  std::function<double(const CompareRow&)> raytraceValue;
  std::function<double(const CompareRow&)> geant4Value;
};

enum class Geant4SourceMode {
  RaytraceFocalPlane,
  PpacProjection,
};

std::string Lowercase(const char* value) {
  std::string text = value ? value : "";
  std::transform(text.begin(), text.end(), text.begin(), [](unsigned char c) {
    return static_cast<char>(std::tolower(c));
  });
  return text;
}

bool ParseGeant4SourceMode(const char* value, Geant4SourceMode& mode) {
  const std::string text = Lowercase(value);
  if (text.empty() || text == "raytrace") {
    mode = Geant4SourceMode::RaytraceFocalPlane;
    return true;
  }
  if (text == "ppac") {
    mode = Geant4SourceMode::PpacProjection;
    return true;
  }
  return false;
}

const char* SourceLabel(bool useRaytraceFocalPlane) {
  return useRaytraceFocalPlane ? "RaytraceFocalPlane hit" : "PPAC1/PPAC2 projection";
}

struct LinearFit {
  double intercept = 0.0;
  double slope = 0.0;
  double r2 = 0.0;
};

std::pair<double, double> Range(const std::vector<double>& values,
                                bool includeZero = false) {
  double lo = *std::min_element(values.begin(), values.end());
  double hi = *std::max_element(values.begin(), values.end());
  if (includeZero) {
    lo = std::min(lo, 0.0);
    hi = std::max(hi, 0.0);
  }
  const double span = hi - lo;
  const double pad = span > 0.0 ? 0.08 * span : 0.08 * std::max(1.0, std::abs(lo));
  return {lo - pad, hi + pad};
}

LinearFit FitLine(const std::vector<double>& x, const std::vector<double>& y) {
  LinearFit fit;
  const double n = static_cast<double>(x.size());
  double meanX = 0.0;
  double meanY = 0.0;
  for (std::size_t i = 0; i < x.size(); ++i) {
    meanX += x[i];
    meanY += y[i];
  }
  meanX /= n;
  meanY /= n;

  double sxx = 0.0;
  double sxy = 0.0;
  for (std::size_t i = 0; i < x.size(); ++i) {
    const double dx = x[i] - meanX;
    sxx += dx * dx;
    sxy += dx * (y[i] - meanY);
  }
  fit.slope = sxx > 0.0 ? sxy / sxx : 0.0;
  fit.intercept = meanY - fit.slope * meanX;

  double ssResidual = 0.0;
  double ssTotal = 0.0;
  for (std::size_t i = 0; i < x.size(); ++i) {
    const double fitted = fit.intercept + fit.slope * x[i];
    ssResidual += (y[i] - fitted) * (y[i] - fitted);
    ssTotal += (y[i] - meanY) * (y[i] - meanY);
  }
  fit.r2 = ssTotal > 0.0 ? 1.0 - ssResidual / ssTotal : 1.0;
  return fit;
}

std::string PlotPngPath(const char* outputPath, const char* canvasName) {
  std::string path = outputPath ? outputPath : "";
  const std::size_t slash = path.find_last_of("/\\");
  const std::string directory =
      slash == std::string::npos ? "" : path.substr(0, slash + 1);
  std::string basename =
      slash == std::string::npos ? path : path.substr(slash + 1);
  const std::size_t dot = basename.find_last_of('.');
  if (dot != std::string::npos) {
    basename = basename.substr(0, dot);
  }
  if (basename.empty()) {
    basename = "raytrace_geant4_compare";
  }
  return directory + basename + "_" + canvasName + ".png";
}

void MakePlot(const std::vector<CompareRow>& rows, const Quantity& q,
              const std::string& pngPath) {
  std::vector<double> raytrace;
  std::vector<double> geant4;
  std::vector<double> residual;
  raytrace.reserve(rows.size());
  geant4.reserve(rows.size());
  residual.reserve(rows.size());
  for (const CompareRow& row : rows) {
    raytrace.push_back(q.raytraceValue(row));
    geant4.push_back(q.geant4Value(row));
    residual.push_back(raytrace.back() - geant4.back());
  }

  std::vector<double> both = raytrace;
  both.insert(both.end(), geant4.begin(), geant4.end());
  const auto xyRange = Range(both);
  const auto residualRange = Range(residual, true);
  const LinearFit fit = FitLine(raytrace, geant4);
  double residualSumSq = 0.0;
  for (const double value : residual) {
    residualSumSq += value * value;
  }
  const double residualRms =
      residual.empty() ? 0.0 : std::sqrt(residualSumSq / residual.size());

  TCanvas canvas(q.key, (std::string(q.title) + " Raytrace vs Geant4").c_str(), 900, 900);
  canvas.SetBorderMode(0);
  canvas.SetBorderSize(0);
  canvas.SetFrameBorderMode(0);
  TPad top("top", "", 0.0, 0.32, 1.0, 1.0);
  TPad bottom("bottom", "", 0.0, 0.0, 1.0, 0.32);
  top.SetBorderMode(0);
  top.SetBorderSize(0);
  top.SetFrameBorderMode(0);
  top.SetBottomMargin(0.03);
  top.SetLeftMargin(0.12);
  top.SetRightMargin(0.04);
  bottom.SetBorderMode(0);
  bottom.SetBorderSize(0);
  bottom.SetFrameBorderMode(0);
  bottom.SetTopMargin(0.04);
  bottom.SetBottomMargin(0.30);
  bottom.SetLeftMargin(0.12);
  bottom.SetRightMargin(0.04);
  top.Draw();
  bottom.Draw();

  top.cd();
  TH2D topFrame("topFrame",
                (std::string(q.title) + ";Raytrace " + q.title + " [" + q.unit +
                 "];Geant4 " + q.title + " [" + q.unit + "]")
                    .c_str(),
                1, xyRange.first, xyRange.second, 1, xyRange.first,
                xyRange.second);
  topFrame.SetStats(false);
  topFrame.GetXaxis()->SetLabelSize(0.0);
  topFrame.GetXaxis()->SetTitleSize(0.0);
  topFrame.Draw();
  TGraph scatter(static_cast<int>(rows.size()), raytrace.data(), geant4.data());
  scatter.SetMarkerStyle(20);
  scatter.Draw("P SAME");
  TLine fitLine(xyRange.first, fit.intercept + fit.slope * xyRange.first,
                xyRange.second, fit.intercept + fit.slope * xyRange.second);
  fitLine.SetLineColor(kRed);
  fitLine.SetLineWidth(2);
  fitLine.Draw("SAME");
  TPaveText fitText(0.15, 0.74, 0.62, 0.90, "NDC");
  fitText.SetFillColor(0);
  fitText.SetFillStyle(0);
  fitText.SetBorderSize(0);
  fitText.SetTextAlign(12);
  fitText.SetTextSize(0.045);
  std::ostringstream equation;
  equation << "Fit: y = " << std::setprecision(5) << fit.intercept << " + "
           << fit.slope << " x";
  std::ostringstream r2Text;
  r2Text << "R^{2} = " << std::setprecision(6) << fit.r2;
  fitText.AddText(equation.str().c_str());
  fitText.AddText(r2Text.str().c_str());
  fitText.Draw("SAME");

  bottom.cd();
  TH2D bottomFrame(
      "bottomFrame",
      (std::string(";Raytrace ") + q.title + " [" + q.unit +
       "];Residual [" + q.unit + "]")
          .c_str(),
      1, xyRange.first, xyRange.second, 1, residualRange.first,
      residualRange.second);
  bottomFrame.SetStats(false);
  bottomFrame.GetXaxis()->SetTitleSize(0.11);
  bottomFrame.GetXaxis()->SetLabelSize(0.10);
  bottomFrame.GetYaxis()->SetTitleSize(0.10);
  bottomFrame.GetYaxis()->SetLabelSize(0.09);
  bottomFrame.GetYaxis()->SetTitleOffset(0.55);
  bottomFrame.Draw();
  TGraph residualGraph(static_cast<int>(rows.size()), raytrace.data(), residual.data());
  residualGraph.SetMarkerStyle(20);
  residualGraph.Draw("P SAME");
  TLine zero(xyRange.first, 0.0, xyRange.second, 0.0);
  zero.SetLineColor(kGray + 2);
  zero.SetLineStyle(2);
  zero.Draw("SAME");
  TPaveText residualText(0.58, 0.72, 0.96, 0.90, "NDC");
  residualText.SetFillColor(0);
  residualText.SetFillStyle(0);
  residualText.SetBorderSize(0);
  residualText.SetTextAlign(12);
  residualText.SetTextSize(0.09);
  std::ostringstream residualLabel;
  residualLabel << "Residual RMS = " << std::fixed << std::setprecision(4)
                << residualRms << " " << q.unit;
  residualText.AddText(residualLabel.str().c_str());
  residualText.Draw("SAME");

  canvas.Modified();
  canvas.Update();
  canvas.Write();
  canvas.SaveAs(pngPath.c_str());
}

double RmsResidual(const std::vector<CompareRow>& rows, const Quantity& q) {
  double sumSq = 0.0;
  for (const CompareRow& row : rows) {
    const double residual = q.raytraceValue(row) - q.geant4Value(row);
    sumSq += residual * residual;
  }
  return rows.empty() ? 0.0 : std::sqrt(sumSq / static_cast<double>(rows.size()));
}

double MaxAbsResidual(const std::vector<CompareRow>& rows, const Quantity& q) {
  double maxAbs = 0.0;
  for (const CompareRow& row : rows) {
    const double residual = q.raytraceValue(row) - q.geant4Value(row);
    maxAbs = std::max(maxAbs, std::abs(residual));
  }
  return maxAbs;
}

}  // namespace

void CompareRaytraceGeant4(const char* inputPath = "SimData_0.root",
                           const char* outputPath = "raytrace_geant4_compare.root",
                           double raytracePlaneZLocalMm = -600.0,
                           double ppac1PlaneZLocalMm = -489.99875,
                           double ppac2PlaneZLocalMm = -89.99875,
                           const char* treeName = "AccurateData",
                           const char* geant4Source = "raytrace") {
  std::unique_ptr<TFile> input(TFile::Open(inputPath, "READ"));
  if (!input || input->IsZombie()) {
    std::cerr << "ERROR: cannot open input file: " << inputPath << '\n';
    return;
  }

  TTree* tree = nullptr;
  input->GetObject(treeName, tree);
  if (!tree) {
    std::cerr << "ERROR: cannot find tree '" << treeName << "' in " << inputPath
              << '\n';
    return;
  }

  bool slitHitTransmitted = false;
  bool raytraceFocalPlaneAccepted = false;
  bool ppac1Accepted = false;
  bool ppac2Accepted = false;

  std::vector<int>* slitTrackId = nullptr;
  std::vector<int>* raytraceFocalPlaneTrackId = nullptr;
  std::vector<double>* raytraceX = nullptr;
  std::vector<double>* raytraceY = nullptr;
  std::vector<double>* raytraceAngleX = nullptr;
  std::vector<double>* raytraceAngleY = nullptr;
  std::vector<double>* raytraceFocalPlaneX = nullptr;
  std::vector<double>* raytraceFocalPlaneY = nullptr;
  std::vector<double>* raytraceFocalPlaneMomentumX = nullptr;
  std::vector<double>* raytraceFocalPlaneMomentumY = nullptr;
  std::vector<double>* raytraceFocalPlaneMomentumZ = nullptr;
  std::vector<double>* ppac1X = nullptr;
  std::vector<double>* ppac1Y = nullptr;
  std::vector<double>* ppac2X = nullptr;
  std::vector<double>* ppac2Y = nullptr;

  const bool hasRaytraceFocalPlane =
      tree->GetBranch("RaytraceFocalPlaneAccepted") &&
      tree->GetBranch("RaytraceFocalPlaneHitLocalPosX") &&
      tree->GetBranch("RaytraceFocalPlaneHitLocalPosY") &&
      tree->GetBranch("RaytraceFocalPlaneHitLocalMomentumX") &&
      tree->GetBranch("RaytraceFocalPlaneHitLocalMomentumY") &&
      tree->GetBranch("RaytraceFocalPlaneHitLocalMomentumZ");

  Geant4SourceMode sourceMode = Geant4SourceMode::RaytraceFocalPlane;
  if (!ParseGeant4SourceMode(geant4Source, sourceMode)) {
    std::cerr << "ERROR: unknown geant4Source '" << geant4Source
              << "'. Use \"raytrace\" or \"ppac\".\n";
    return;
  }

  const bool useRaytraceFocalPlane =
      sourceMode == Geant4SourceMode::RaytraceFocalPlane;

  if (sourceMode == Geant4SourceMode::RaytraceFocalPlane &&
      !hasRaytraceFocalPlane) {
    std::cerr << "ERROR: geant4Source=\"raytrace\" was requested, but "
              << "RaytraceFocalPlaneHit* branches are not present in "
              << inputPath << ".\n";
    return;
  }

  tree->SetBranchAddress("SlitHitTransmitted", &slitHitTransmitted);
  tree->SetBranchAddress("SlitHitTrackId", &slitTrackId);
  tree->SetBranchAddress("Ppac1Accepted", &ppac1Accepted);
  tree->SetBranchAddress("Ppac2Accepted", &ppac2Accepted);
  tree->SetBranchAddress("RaytracePositionX", &raytraceX);
  tree->SetBranchAddress("RaytracePositionY", &raytraceY);
  tree->SetBranchAddress("RaytraceAngleX", &raytraceAngleX);
  tree->SetBranchAddress("RaytraceAngleY", &raytraceAngleY);
  if (useRaytraceFocalPlane) {
    tree->SetBranchAddress("RaytraceFocalPlaneAccepted", &raytraceFocalPlaneAccepted);
    tree->SetBranchAddress("RaytraceFocalPlaneHitTrackId", &raytraceFocalPlaneTrackId);
    tree->SetBranchAddress("RaytraceFocalPlaneHitLocalPosX", &raytraceFocalPlaneX);
    tree->SetBranchAddress("RaytraceFocalPlaneHitLocalPosY", &raytraceFocalPlaneY);
    tree->SetBranchAddress("RaytraceFocalPlaneHitLocalMomentumX", &raytraceFocalPlaneMomentumX);
    tree->SetBranchAddress("RaytraceFocalPlaneHitLocalMomentumY", &raytraceFocalPlaneMomentumY);
    tree->SetBranchAddress("RaytraceFocalPlaneHitLocalMomentumZ", &raytraceFocalPlaneMomentumZ);
  }
  tree->SetBranchAddress("Ppac1HitLocalPosX", &ppac1X);
  tree->SetBranchAddress("Ppac1HitLocalPosY", &ppac1Y);
  tree->SetBranchAddress("Ppac2HitLocalPosX", &ppac2X);
  tree->SetBranchAddress("Ppac2HitLocalPosY", &ppac2Y);

  std::unique_ptr<TFile> output(TFile::Open(outputPath, "RECREATE"));
  if (!output || output->IsZombie()) {
    std::cerr << "ERROR: cannot create output file: " << outputPath << '\n';
    return;
  }

  auto* residualTree = new TTree(
      "RaytraceGeant4Residuals",
      "Geant4 PPAC track projected to the Raytrace focal plane");

  Long64_t event = 0;
  double geant4X = 0.0;
  double rawGeant4X = 0.0;
  double geant4Y = 0.0;
  double geant4AngleX = 0.0;
  double geant4AngleY = 0.0;
  double mdmX = 0.0;
  double mdmY = 0.0;
  double mdmAngleX = 0.0;
  double mdmAngleY = 0.0;
  double dx = 0.0;
  double dy = 0.0;
  double dAngleX = 0.0;
  double dAngleY = 0.0;
  double ppac1LocalZ = 0.0;
  double ppac2LocalZ = 0.0;

  residualTree->Branch("Event", &event, "Event/L");
  residualTree->Branch("RaytracePlaneZLocalMm", &raytracePlaneZLocalMm,
                       "RaytracePlaneZLocalMm/D");
  residualTree->Branch("RawGeant4X", &rawGeant4X, "RawGeant4X/D");
  residualTree->Branch("Geant4X", &geant4X, "Geant4X/D");
  residualTree->Branch("Geant4Y", &geant4Y, "Geant4Y/D");
  residualTree->Branch("Geant4AngleXDeg", &geant4AngleX,
                       "Geant4AngleXDeg/D");
  residualTree->Branch("Geant4AngleYDeg", &geant4AngleY,
                       "Geant4AngleYDeg/D");
  residualTree->Branch("RaytraceX", &mdmX, "RaytraceX/D");
  residualTree->Branch("RaytraceY", &mdmY, "RaytraceY/D");
  residualTree->Branch("RaytraceAngleXDeg", &mdmAngleX, "RaytraceAngleXDeg/D");
  residualTree->Branch("RaytraceAngleYDeg", &mdmAngleY, "RaytraceAngleYDeg/D");
  residualTree->Branch("DeltaX", &dx, "DeltaX/D");
  residualTree->Branch("DeltaY", &dy, "DeltaY/D");
  residualTree->Branch("DeltaAngleXDeg", &dAngleX, "DeltaAngleXDeg/D");
  residualTree->Branch("DeltaAngleYDeg", &dAngleY, "DeltaAngleYDeg/D");
  residualTree->Branch("Ppac1LocalZ", &ppac1LocalZ, "Ppac1LocalZ/D");
  residualTree->Branch("Ppac2LocalZ", &ppac2LocalZ, "Ppac2LocalZ/D");

  auto* hDx = new TH1D("hDeltaX", "Raytrace - Geant4;#Delta x [mm];Events", 200, -100.0, 100.0);
  auto* hDy = new TH1D("hDeltaY", "Raytrace - Geant4;#Delta y [mm];Events", 200, -100.0, 100.0);
  auto* hDAngleX = new TH1D("hDeltaAngleX", "Raytrace - Geant4;#Delta x' [deg];Events", 200, -10.0, 10.0);
  auto* hDAngleY = new TH1D("hDeltaAngleY", "Raytrace - Geant4;#Delta y' [deg];Events", 200, -10.0, 10.0);
  hDx->SetTitle("Raytrace - Geant4;#Delta x [mm];Events");
  hDy->SetTitle("Raytrace - Geant4;#Delta y [mm];Events");
  hDAngleX->SetTitle("Raytrace - Geant4;#Delta x' [deg];Events");
  hDAngleY->SetTitle("Raytrace - Geant4;#Delta y' [deg];Events");

  std::vector<CompareRow> rows;
  Long64_t selected = 0;
  Long64_t skipped = 0;
  Long64_t raytraceStopped = 0;
  const Long64_t entries = tree->GetEntries();
  for (Long64_t i = 0; i < entries; ++i) {
    tree->GetEntry(i);

    if (!ppac1Accepted || !ppac2Accepted || !slitHitTransmitted ||
        !HasValue(raytraceX) || !HasValue(raytraceY) || !HasValue(raytraceAngleX) ||
        !HasValue(raytraceAngleY)) {
      ++skipped;
      continue;
    }

    if (useRaytraceFocalPlane) {
      if (!raytraceFocalPlaneAccepted || !HasValue(raytraceFocalPlaneX) ||
          !HasValue(raytraceFocalPlaneY) ||
          !HasValue(raytraceFocalPlaneMomentumX) ||
          !HasValue(raytraceFocalPlaneMomentumY) ||
          !HasValue(raytraceFocalPlaneMomentumZ)) {
        ++skipped;
        continue;
      }
    } else if (!HasValue(ppac1X) || !HasValue(ppac1Y) ||
               !HasValue(ppac2X) || !HasValue(ppac2Y)) {
      ++skipped;
      continue;
    }

    mdmX = FirstOrNaN(raytraceX);
    mdmY = FirstOrNaN(raytraceY);
    mdmAngleX = FirstOrNaN(raytraceAngleX);
    mdmAngleY = FirstOrNaN(raytraceAngleY);
    if (std::abs(mdmX) > 1.0e10 || std::abs(mdmY) > 1.0e10) {
      ++raytraceStopped;
      continue;
    }

    if (useRaytraceFocalPlane) {
      const int trackId = FirstTrackOr(slitTrackId, FirstTrackOr(raytraceFocalPlaneTrackId, -99));
      const std::size_t focalIndex = IndexForTrack(raytraceFocalPlaneTrackId, trackId);
      geant4X = AtOrNaN(raytraceFocalPlaneX, focalIndex);
      geant4Y = AtOrNaN(raytraceFocalPlaneY, focalIndex);
      rawGeant4X = geant4X;
      const double px = AtOrNaN(raytraceFocalPlaneMomentumX, focalIndex);
      const double py = AtOrNaN(raytraceFocalPlaneMomentumY, focalIndex);
      const double pz = AtOrNaN(raytraceFocalPlaneMomentumZ, focalIndex);
      if (!std::isfinite(geant4X) || !std::isfinite(geant4Y) ||
          !std::isfinite(px) || !std::isfinite(py) || !std::isfinite(pz) ||
          std::abs(pz) < 1.0e-12) {
        ++skipped;
        continue;
      }
      geant4AngleX = std::atan2(px, pz) * TMath::RadToDeg();
      geant4AngleY =
          std::atan2(py, std::sqrt(px * px + pz * pz)) * TMath::RadToDeg();
      ppac1LocalZ = raytracePlaneZLocalMm;
      ppac2LocalZ = raytracePlaneZLocalMm;
    } else {
      const double x1 = ppac1X->front();
      const double y1 = ppac1Y->front();
      const double z1 = ppac1PlaneZLocalMm;
      const double x2 = ppac2X->front();
      const double y2 = ppac2Y->front();
      const double z2 = ppac2PlaneZLocalMm;
      const double dzPpac = z2 - z1;
      if (std::abs(dzPpac) < 1.0e-9) {
        ++skipped;
        continue;
      }

      const double t = (raytracePlaneZLocalMm - z1) / dzPpac;
      rawGeant4X = x1 + t * (x2 - x1);
      geant4X = rawGeant4X;
      geant4Y = y1 + t * (y2 - y1);
      geant4AngleX = std::atan2(x2 - x1, dzPpac) * TMath::RadToDeg();
      geant4AngleY =
          std::atan2(y2 - y1,
                     std::sqrt((x2 - x1) * (x2 - x1) + dzPpac * dzPpac)) *
          TMath::RadToDeg();
      ppac1LocalZ = z1;
      ppac2LocalZ = z2;
    }

    dx = mdmX - geant4X;
    dy = mdmY - geant4Y;
    dAngleX = mdmAngleX - geant4AngleX;
    dAngleY = mdmAngleY - geant4AngleY;
    event = i;

    rows.push_back({mdmX, mdmY, mdmAngleX, mdmAngleY,
                    geant4X, geant4Y, geant4AngleX, geant4AngleY});

    residualTree->Fill();
    hDx->Fill(dx);
    hDy->Fill(dy);
    hDAngleX->Fill(dAngleX);
    hDAngleY->Fill(dAngleY);
    ++selected;
  }

  const std::vector<Quantity> quantities{
      {"X", "X", "mm",
       [](const CompareRow& r) { return r.raytraceX; },
       [](const CompareRow& r) { return r.geant4X; }},
      {"Y", "Y", "mm",
       [](const CompareRow& r) { return r.raytraceY; },
       [](const CompareRow& r) { return r.geant4Y; }},
      {"AngX", "AngX", "deg",
       [](const CompareRow& r) { return r.raytraceAngleX; },
       [](const CompareRow& r) { return r.geant4AngleX; }},
      {"AngY", "AngY", "deg",
       [](const CompareRow& r) { return r.raytraceAngleY; },
       [](const CompareRow& r) { return r.geant4AngleY; }},
  };

  output->cd();
  gStyle->SetOptStat(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetFrameBorderMode(0);
  std::vector<std::string> pngPaths;
  if (!rows.empty()) {
    for (const Quantity& quantity : quantities) {
      pngPaths.push_back(PlotPngPath(outputPath, quantity.key));
      MakePlot(rows, quantity, pngPaths.back());
      output->cd();
    }
  }
  residualTree->Write();
  hDx->Write();
  hDy->Write();
  hDAngleX->Write();
  hDAngleY->Write();
  output->Close();

  std::cout << "Input entries: " << entries << '\n'
            << "Compared entries: " << selected << '\n'
            << "Skipped entries: " << skipped << '\n'
            << "Raytrace stopped entries: " << raytraceStopped << '\n'
            << "Geant4 source: " << SourceLabel(useRaytraceFocalPlane) << '\n'
            << "Raytrace plane local z: " << raytracePlaneZLocalMm << " mm\n"
            << "PPAC1 plane local z: " << ppac1PlaneZLocalMm << " mm\n"
            << "PPAC2 plane local z: " << ppac2PlaneZLocalMm << " mm\n"
            << "PPAC-to-Raytrace X origin offset: 0 mm\n"
            << "Wrote " << outputPath << " with 4 canvases and " << rows.size()
            << " rays\n";
  if (!pngPaths.empty()) {
    std::cout << "Saved PNG plots:\n";
    for (const std::string& pngPath : pngPaths) {
      std::cout << "  " << pngPath << '\n';
    }
  }
  std::cout
            << "Residual convention: Raytrace - Geant4\n";
  if (rows.empty()) {
    std::cout << "No comparable rays found; canvases and PNG plots were not written.\n";
    return;
  }
  std::cout << std::fixed << std::setprecision(6);
  for (const Quantity& quantity : quantities) {
    std::cout << std::setw(6) << quantity.key << " RMS="
              << RmsResidual(rows, quantity)
              << " MaxAbs=" << MaxAbsResidual(rows, quantity) << " "
              << quantity.unit << "\n";
  }
}
