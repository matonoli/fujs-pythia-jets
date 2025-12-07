// main.cc
#include "Pythia8/Pythia.h"
#include <cmath>
#include <iostream>

#include "fastjet/ClusterSequence.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/GhostedAreaSpec.hh"

#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TString.h"
#include "TVector2.h"

using namespace Pythia8;

int main(int argc, char *argv[]) {


    // jet parameters
  const double jetRadius = 0.4;
  const double jetPtMin = 3.0;
  // particle parameters
  const double particlePtMin = 0.15;
  const double particleEtaMax = 1.5;
  const double jetEtaMax = particleEtaMax - jetRadius;



  if (argc < 3) {
    std::cerr << "Usage: " << argv[0]
              << " pTHatMin pTHatMax|inf [nEvents=100]";
    return 1;
  }
  // Required: pTHatMin
  const double ptHatMin = std::stod(argv[1]);

  // Required: pTHatMax (can be "inf" or negative for open upper bound)
  double ptHatMax;
  std::string s = argv[2];
  if (s == "inf" || s == "Inf" || s == "INF") {
    ptHatMax = -1.0;
  } else {
    ptHatMax = std::stod(s);
  }

  if (ptHatMax > 0.0 && ptHatMax < ptHatMin) {
    std::cerr << "[error] pTHatMax < pTHatMin\n";
    return 1;
  }

  int nEvents = (argc > 3) ? std::atoi(argv[3]) : 100;

  const TString outFile = Form("jets_ptHat_%.0f_%.0f_jetR%.1f.root", 
                                    ptHatMin, ptHatMax > 0.0 ? ptHatMax : 999.0,
                                    jetRadius);

  // --- Pythia setup ---
  Pythia8::Pythia pythia8;
  pythia8.readString("Beams:idA = 2212"); // proton
  pythia8.readString("Beams:idB = 2212"); // proton
  pythia8.readString("Beams:eCM = 200."); // 200 GeV CM energy
  pythia8.readString("HardQCD:all = on"); // enable all hard QCD processes

  // Phase space cuts
  pythia8.readString(Form("PhaseSpace:pTHatMin = %.1f", ptHatMin));  // set pTHat min
  if (ptHatMax > 0)
    pythia8.readString(Form("PhaseSpace:pTHatMax = %.1f", ptHatMax)); 

  // Init
  if (!pythia8.init()) {
    std::cerr << "[error] PYTHIA init() failed.\n";
    return 2;
  }
  // --- ROOT output ---
  TFile *fout = new TFile(outFile, "RECREATE");
  TTree *tree = new TTree("events", "");

  // Jet branches (store up to 100 jets in single event)
  const int kMaxJets = 100;

// jet branches
  int nJets = 0;
  double pt[kMaxJets];
  double eta[kMaxJets];
  double phi[kMaxJets];
  double area[kMaxJets];
  tree->Branch("nJets", &nJets, "nJets/I");
  tree->Branch("pt", pt, "pt[nJets]/D");
  tree->Branch("eta", eta, "eta[nJets]/D");
  tree->Branch("phi", phi, "phi[nJets]/D");
  tree->Branch("area", area, "area[nJets]/D");

  // Event loop
  for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
    if (!pythia8.next())
      continue;

    // Build input particles for jet finding
    std::vector<fastjet::PseudoJet> particlesContainer;

    for (int iParticle = 0; iParticle < pythia8.event.size(); ++iParticle) {
      const auto &particle = pythia8.event[iParticle];
      // final-state, visible (no neutrinos), basic kinematic filter
      if (!particle.isFinal() || !particle.isVisible() || !particle.isCharged())
        continue;
      if (std::abs(particle.eta()) > particleEtaMax)
        continue; // wide acceptance for clustering
      if (particle.pT() < particlePtMin)
        continue;

      fastjet::PseudoJet input(particle.px(), particle.py(), particle.pz(), particle.e());
      input.set_user_index(iParticle); // <— keep Pythia index to recover charge
      particlesContainer.push_back(input);
    }
    //----------------------------
    // FASTJET jet finding PART:
    //---------------------------- 

    fastjet::JetDefinition jetDefinition(fastjet::antikt_algorithm, jetRadius);
    //fastjet::JetDefinition jetDefinition(fastjet::cambridge_aachen_algorithm, jetRadius);

    fastjet::Selector select_eta = fastjet::SelectorAbsEtaMax(jetEtaMax);
    fastjet::Selector select_pt = fastjet::SelectorPtMin(jetPtMin);
    fastjet::Selector select_pt_eta = select_pt && select_eta;

    fastjet::GhostedAreaSpec areaSpec(2.5);
    fastjet::AreaDefinition areaDefinition(fastjet::active_area, areaSpec);
    fastjet::ClusterSequenceArea clusterSeqArea(particlesContainer, jetDefinition, areaDefinition);

    auto inclusive_jets = clusterSeqArea.inclusive_jets();
    auto sorted_jets = fastjet::sorted_by_pt(inclusive_jets); // here
    auto selected_jets = select_pt_eta(sorted_jets);

    nJets = selected_jets.size();

    for (size_t iJet = 0; iJet < selected_jets.size(); ++iJet) {
      pt[iJet] = selected_jets[iJet].pt();
      eta[iJet] = selected_jets[iJet].eta();
      float phiNormalized = TVector2::Phi_mpi_pi(selected_jets[iJet].phi()); // normalize to [-pi, pi]
      phi[iJet] = phiNormalized;
      area[iJet] = selected_jets[iJet].area();
    }

  // Dump all jets for this event
    cout <<"===================================================="<< endl;
    cout << "Event " << iEvent << ": Number of jets = " << nJets << endl;
    for (int k = 0; k < nJets; ++k) {
      cout << "  Jet " << k << ": pT = " << pt[k] << ", eta = " << eta[k] << ", phi = " << phi[k] << endl;
    }


    tree->Fill();
  }

  // Cross sections (mb)
  const double sigmaGen = pythia8.info.sigmaGen();
  const double sigmaErr = pythia8.info.sigmaErr();

  TH1D *stats = new TH1D("stats", "stats", 4, 0, 4);
  vector<TString> statNames = {"nEvents", "sigmaGen_mb", "ptHatMin",
                               "ptHatMax"};
  for (size_t i = 0; i < statNames.size(); ++i)
    stats->GetXaxis()->SetBinLabel(i + 1, statNames[i]);

  stats->SetBinContent(1, nEvents);
  stats->SetBinContent(2, sigmaGen);
  stats->SetBinError(2, sigmaErr);
  stats->SetBinContent(3, ptHatMin);
  stats->SetBinContent(4, ptHatMax);

  // Print and record
  std::cout << "[done] Wrote " << outFile << "\n"
            << "       sigmaGen   = " << sigmaGen << " mb  (± " << sigmaErr
            << ")\n"
            << "ptHatMin = " << ptHatMin
            << "\n"
               "ptHatMax = "
            << ptHatMax << "\n";

  pythia8.stat();

  fout->Write();
  fout->Close();
  delete fout;

  return 0;
}
