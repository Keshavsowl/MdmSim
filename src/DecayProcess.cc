#include "DecayProcess.hh"
#include "Analysis.hh"
#include "TLorentzVector.h"
#include "G4DynamicParticle.hh"
#include "TVector3.h"

#include "G4Gamma.hh"
#include "G4IonTable.hh"
#include "G4ParticleTable.hh"
#include "G4Ions.hh"
#include "G4ParticleChange.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "Randomize.hh"

#include <cfloat>
#include <cmath>
#include <stdexcept>
#include <algorithm>
namespace
{
	G4ParticleDefinition *GetDecayParticleDefinition(
			const MdmSim::DecayProduct &product)
	{
		if (product.particleType == "gamma")
		{
			return G4Gamma::GammaDefinition();
		}
		if (product.particleType != "ion")
		{
			throw std::runtime_error(
				"Unsupported decay product particle type:" + product.particleType);
		}

		auto *particleTable = G4ParticleTable::GetParticleTable();
		
		if(product.excitationEnergy== 0.)
		{
		if(product.Z == 0 && product.A == 1)
			return particleTable->FindParticle("neutron");
		if(product.Z == 1 && product.A ==1)
			return particleTable->FindParticle("proton");
		if(product.Z == 1 && product.A ==2)
		       return particleTable->FindParticle("deuteron");
	        if(product.Z == 1 && product.A ==3)
			return particleTable->FindParticle("triton");
		if(product.Z == 2 && product.A ==3)
			return particleTable->FindParticle("He3");
		if(product.Z == 2 && product.A ==4)
			return particleTable->FindParticle("alpha");
	}
	return G4IonTable::GetIonTable()->GetIon(product.Z, product.A, product.excitationEnergy);

	}
}	
namespace MdmSim
{
	DecayProcess::DecayProcess(const G4String &processName): G4VDiscreteProcess(processName, fDecay)
	{ 
	}

	DecayProcess::~DecayProcess() = default;

	const DecayChannel *DecayProcess::FindMatchingChannel( const G4Track &aTrack) const
        {
		const auto *definition = aTrack.GetParticleDefinition();
		if(!definition)
			return nullptr;
		const G4int z = definition->GetAtomicNumber();
		const G4int a = definition->GetAtomicMass();

		G4double excitationEnergy = 0.;

		const auto *ion = dynamic_cast<const G4Ions *>(definition);

		if(ion)
			excitationEnergy = ion->GetExcitationEnergy();

		for (const auto &channel :fDecayChannels)
		{
			const G4bool sameZ = z == channel.parentZ;

			const G4bool sameA = a == channel.parentA;

			const G4bool sameExcitation = std::abs(excitationEnergy - channel.parentExcitationEnergy) < 1.0 * eV;

					if(sameZ && sameA && sameExcitation)
					return &channel;
					}
					return nullptr;
	
	}

	void DecayProcess::TwoBodyDecay( const G4Track &aTrack, const DecayChannel &channel)
	{
		G4ParticleDefinition *product1Definition = GetDecayParticleDefinition(channel.product1);
		G4ParticleDefinition *product2Definition = GetDecayParticleDefinition(channel.product2);
		
		if(!product1Definition || !product2Definition)
		{
			throw std::runtime_error(
					"Unable to create decay product particle definitions"
					);
		}
		G4cout << "[DecayProcess] Decaying Z="<< channel.parentZ <<"A=" <<channel.parentA<<"Ex="<<channel.parentExcitationEnergy / MeV <<"MeV"<<G4endl;

		const G4double product1Mass = product1Definition->GetPDGMass();
		const G4double product2Mass = product2Definition->GetPDGMass();
		const G4LorentzVector parent4G4 = aTrack.GetDynamicParticle()->Get4Momentum();
		TLorentzVector parent4( parent4G4.px(), parent4G4.py(), parent4G4.pz(), parent4G4.e());

		const G4double s = parent4.M2();
		const G4double sumMass = product1Mass + product2Mass;

		if(s < sumMass * sumMass)
		{
			throw std::runtime_error("Configured decay is energetically closed");
		}

		const G4double diffMass = product1Mass - product2Mass;
		const G4double sqrtS = std::sqrt(s);
		const G4double pRest = std::sqrt(std::max(0.,(s- sumMass * sumMass)*( s - diffMass * diffMass)))/(2. * sqrtS);
		const G4double e1Rest = (s + product1Mass * product1Mass - product2Mass * product2Mass)/ (2. *sqrtS);
		const G4double e2Rest = (s + product2Mass * product2Mass - product1Mass * product1Mass)/(2. *sqrtS);
		TVector3 product1MomentumRest(pRest, 0., 0.);
		product1MomentumRest.SetTheta(std::acos(-1.0 + 2.0 * G4UniformRand()));
		product1MomentumRest.SetPhi(G4UniformRand() * 2. * CLHEP::pi);
		TLorentzVector product1FourRest(product1MomentumRest, e1Rest);
		TLorentzVector product2FourRest(-product1MomentumRest,e2Rest);
		const TVector3 betaParent = parent4.BoostVector();

		const TVector3 restMomentumSum = product1FourRest.Vect() + product2FourRest.Vect();
		const G4double restDeltaE = product1FourRest.E() + product2FourRest.E() - sqrtS;

		G4cout << "Rest-frame check:"
                       << "  |p1*|=" << product1FourRest.P() / MeV
                       << " MeV/c"
                       << "  |p2*|=" << product2FourRest.P() / MeV
                       << " MeV/c"
                      << "  |p1*+ p2*|=" << restMomentumSum.Mag() / MeV
                      << " MeV/c"
                     << "  dE*=" << restDeltaE / MeV
                    << " MeV"
                    << G4endl;



		product1FourRest.Boost(betaParent);
		product2FourRest.Boost(betaParent);

		const G4double product1KE = product1FourRest.E() - product1Mass;
		if (channel.parentZ == 8 && channel.parentA ==16)
		{
			const TVector3 alpha2MomentumLab = product1FourRest.Vect();
			Analysis *analysis = Analysis::Instance();
			analysis->SetAlpha2Generated(true);
			analysis->SetAlpha2Energy(product1KE / MeV);
			analysis->SetAlpha2Theta(alpha2MomentumLab.Theta() / deg);
			analysis->SetAlpha2Phi(alpha2MomentumLab.Phi() /deg);

			G4cout << "alpha2 lab:" << "E="<< product1KE/MeV
				<<"MeV"
				<<"theta="<<alpha2MomentumLab.Theta() / deg
				<<" deg"
				<<" phi ="<<alpha2MomentumLab.Phi() /deg
				<<" deg"
				<<G4endl;
		}
		G4cout << "product1 : " << product1Definition->GetParticleName() <<"KE="<<product1KE / MeV << "MeV" << " |p|=" <<product1FourRest.P() / MeV <<"MeV/c" << G4endl;

                const G4double product2KE = product2FourRest.E() - product2Mass;
		if (product2Definition == G4Gamma::GammaDefinition())
		{
		const TVector3 gammaMomentumLab = product2FourRest.Vect();
		Analysis *analysis = Analysis::Instance();
		
		analysis->SetGammaGenerated(true);
		analysis->SetGammaEnergy(product2FourRest.E() / MeV);
		analysis->SetGammaTheta(gammaMomentumLab.Theta() / deg);
		analysis->SetGammaPhi(gammaMomentumLab.Phi() / deg);

		G4cout << "gamma lab:"
			<< "E=" << product2FourRest.E() / MeV
		       << "theta = "<< gammaMomentumLab.Theta() / deg
	       		<< "deg"
	 		<< "phi="<<gammaMomentumLab.Phi() / deg
			<< " deg"
			<< G4endl;		
		}




	       G4cout << "Product 2:"<< product2Definition->GetParticleName() << " KE=" << product2KE / MeV << " MeV" << " |p|="<< product2FourRest.P() /MeV <<"MeV/c"<<G4endl;
	       TLorentzVector daughterSum4 = product1FourRest + product2FourRest;
	       const G4double deltaE = daughterSum4.E() - parent4.E();
	       const TVector3 deltaP = daughterSum4.Vect() - parent4.Vect();
	       G4cout << "Conservation check:" << " dE= " << deltaE / MeV << "MeV" << " |dP|= "<<deltaP.Mag() /MeV <<" MeV/c "<<G4endl;





	       auto *product1 = new G4DynamicParticle(product1Definition, G4ThreeVector(product1FourRest.Px(), product1FourRest.Py(), product1FourRest.Pz()));
		auto *product2 = new G4DynamicParticle(product2Definition, G4ThreeVector(product2FourRest.Px(), product2FourRest.Py(), product2FourRest.Pz()));
		if (channel.product1.particleType == "ion")
		{
			product1->SetCharge(channel.product1.charge *eplus);
		}
		if (channel.product2.particleType == "ion")
		{
			product2->SetCharge(channel.product2.charge *eplus);
		}

		auto *product1Track = new G4Track(product1, aTrack.GetGlobalTime(), aTrack.GetPosition());
		auto *product2Track = new G4Track(product2, aTrack.GetGlobalTime(), aTrack.GetPosition());
		
		aParticleChange.SetNumberOfSecondaries(2);
		aParticleChange.AddSecondary(product1Track);
		aParticleChange.AddSecondary(product2Track);







	}



	G4double DecayProcess::GetMeanFreePath( const G4Track &, G4double, G4ForceCondition *condition)
	{
		*condition = NotForced;
		if(!fActiveChannel || !fWillDecay)
				return DBL_MAX;
		return 0.;
	}

	G4VParticleChange *DecayProcess::PostStepDoIt(const G4Track &aTrack, const G4Step &)
	{
		aParticleChange.Initialize(aTrack);

		if(!fActiveChannel || !fWillDecay)
			return &aParticleChange;
		TwoBodyDecay(aTrack, *fActiveChannel);
		aParticleChange.ProposeTrackStatus(fStopAndKill);
		aParticleChange.ProposeEnergy(0.);

		return &aParticleChange;
	}

	void DecayProcess::StartTracking(G4Track *track)
	{
		G4VProcess::StartTracking(track);
		fActiveChannel = nullptr;
		fWillDecay = false;

		if (!track)
			return;
		fActiveChannel = FindMatchingChannel(*track);

		if(fActiveChannel)
		{
			fWillDecay = G4UniformRand() < fActiveChannel->probability;

		}
	}





	}
