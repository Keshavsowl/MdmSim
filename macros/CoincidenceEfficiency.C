// Coincidence efficiencies for the reaction/decay chain, e.g.
// 14N(6Li,alpha1)16O*, 16O* -> alpha2 + 12C*, 12C* -> 12C + gamma.
//
// Usage from the repository root:
//   root -l -b -q 'macros/CoincidenceEfficiency.C("build/SimData_0.root")'
//
// Detector flags are matched to the tagged particles by track ID (see
// EventAction): Alpha1AtMdm, Alpha2AtS2, GammaAtNaI, TripleCoincidence.
//
// Two ways of stating the result:
//
// 1. Given alpha1 through the MDM: the fraction of events with alpha1 at the
//    MDM focal plane that also have alpha2 in S2 and/or gamma in NaI. This is
//    independent of how alpha1 directions were generated.
//
// 2. Effective solid angle (msr) of each condition for alpha1:
//      yield = N_beam * N_target * (dsigma/dOmega)_lab(alpha1 at MDM angle) * dOmega_eff
//    It needs alpha1 generated uniformly over a known lab solid angle, i.e. a run
//    with ReactionLightConeHalfAngleInDeg > 0 and a cone larger than the MDM
//    entrance. Within the cone, alpha1 is generated isotropically in the CM;
//    over a cone of a few degrees this is uniform in the lab to good accuracy.

#include <cmath>
#include <cstdio>

#include "TFile.h"
#include "TMath.h"
#include "TTree.h"

void CoincidenceEfficiency(const char *fileName = "SimData_0.root",
                           double coneHalfAngleDeg = 3.5)
{
    TFile file(fileName);
    TTree *tree = (TTree *)file.Get("AccurateData");
    if (!tree)
    {
        printf("No AccurateData tree in %s\n", fileName);
        return;
    }

    const Long64_t nEvents = tree->GetEntries("Alpha1Generated");
    const Long64_t nMdm = tree->GetEntries("Alpha1Generated && Alpha1AtMdm");
    if (nEvents == 0)
    {
        printf("No reaction events in %s\n", fileName);
        return;
    }
    const double coneSolidAngleMsr =
        2. * TMath::Pi() * (1. - std::cos(coneHalfAngleDeg * TMath::DegToRad())) * 1000.;

    const auto count = [&](const char *selection)
    {
        return tree->GetEntries(TString::Format("Alpha1Generated && (%s)", selection));
    };

    printf("\nFile: %s\n", fileName);
    printf("Reaction events: %lld (alpha1 generated in a %.2f deg cone = %.3f msr)\n",
           nEvents, coneHalfAngleDeg, coneSolidAngleMsr);
    printf("alpha1 through the MDM (focal plane): %lld\n", nMdm);

    printf("\n1. Given alpha1 through the MDM\n");
    printf("  %-28s %8s   %18s\n", "Also detected", "Events", "Fraction");
    const auto given = [&](const char *label, const char *selection)
    {
        const Long64_t n = count(TString::Format("Alpha1AtMdm && (%s)", selection));
        const double p = nMdm > 0 ? double(n) / double(nMdm) : 0.;
        const double e = nMdm > 0 ? std::sqrt(p * (1. - p) / double(nMdm)) : 0.;
        printf("  %-28s %8lld   %7.3f +- %6.3f\n", label, n, p, e);
    };
    given("alpha1 also at both PPACs", "Alpha1AtPpacs");
    given("alpha2 in S2", "Alpha2AtS2");
    given("gamma in NaI", "GammaAtNaI");
    given("alpha2 in S2 AND gamma in NaI", "Alpha2AtS2 && GammaAtNaI");

    printf("\n2. Effective solid angle for alpha1 (valid if the cone covers the MDM entrance)\n");
    printf("  %-34s %8s   %18s\n", "Condition", "Events", "dOmega_eff (msr)");
    const auto solidAngle = [&](const char *label, const char *selection)
    {
        const Long64_t n = count(selection);
        const double p = double(n) / double(nEvents);
        const double e = std::sqrt(p * (1. - p) / double(nEvents));
        printf("  %-34s %8lld   %7.3f +- %6.3f\n", label, n,
               p * coneSolidAngleMsr, e * coneSolidAngleMsr);
    };
    solidAngle("alpha1 @ MDM", "Alpha1AtMdm");
    solidAngle("alpha1 @ MDM && alpha2 @ S2", "Alpha1AtMdm && Alpha2AtS2");
    solidAngle("alpha1 @ MDM && gamma @ NaI", "Alpha1AtMdm && GammaAtNaI");
    solidAngle("TRIPLE (MDM && S2 && NaI)", "TripleCoincidence");
    printf("\n");
}
