#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"
#include "TMath.h"   // for TMath::Pi(), TMath::TwoPi()
#include <iostream>

void anaTree(TString inputFileName="jets_ptHat_30_50_jetR0.4.root") {

  TFile outFile("anaTrees.root", "RECREATE");
  // -------------------------
  // 1) Define histograms
  // -------------------------

  // Event-level
  TH1D hNJetsPerEvent("hNJetsPerEvent", "Number of jets per event;N_{jets};Counts", 12, 0, 12); // 0 to 12 jets

  // Single-jet level
  TH1D hJetPt("hJetPt", ";jet p_{T} (GeV/c);Counts", 100, 0, 100); // 0 to 100 GeV/c
  TH1D hJetEta("hJetEta", ";jet #eta;Counts", 60, -3, 3); // -3 to 3 in eta
  TH1D hJetPhi("hJetPhi", ";jet #phi;Counts", 64, -TMath::Pi(), TMath::Pi()); // -pi to pi in phi
  TH1D hJetArea("hJetArea", ";jet area;Counts", 50, 0, 1); // 0 to 1 in area

  // Dijet level
  TH1D hAsymmetry("hAsymmetry", ";Dijet asymmetry A_{J};Counts", 50, 0, 1); // from 0 to 1
  TH1D hDeltaPhi("hDeltaPhi", ";#Delta#phi;Counts", 200, 0, TMath::TwoPi()); // 0 to 2pi

  // -------------------------
  // 2) Open input file & tree
  // -------------------------


  TFile inputFile(inputFileName);
    if (inputFile.IsZombie()) {
    std::cout << "Cannot open input file " << inputFileName << std::endl;
    return; // exit the program if file cannot be opened
  }

  // Read number of generated events from the "stats" histogram
  // vector<TString> statNames = {"nEvents", "sigmaGen_mb", "ptHatMin", "ptHatMax"};
  TH1D *stats = (TH1D *)inputFile.Get("stats");
  int nEvents = stats->GetBinContent(1);
  std::cout <<"Total number of analyzed events: "<< nEvents << std::endl;

  TTree *tree = (TTree *)inputFile.Get("events");
  // -------------------------
  // 3) Set up branches
  // -------------------------

  const int kMaxJets = 100;
  // jet branches
  int nJets = 0;
  double pt[kMaxJets];
  double eta[kMaxJets];
  double phi[kMaxJets];
  double area[kMaxJets];

  tree->SetBranchAddress("nJets", &nJets);
  tree->SetBranchAddress("pt", pt);
  tree->SetBranchAddress("eta", eta);
  tree->SetBranchAddress("phi", phi);
  tree->SetBranchAddress("area", area);


  // -------------------------
  // 4) Event loop
  // -------------------------


  Long64_t nEntries = tree->GetEntries();

  for (Long64_t iEvent = 0; iEvent < nEntries; ++iEvent) {
    tree->GetEntry(iEvent);

    // Count jets in this event
    hNJetsPerEvent.Fill(nJets);

    // Fill single-jet histograms
    for (int iJet = 0; iJet < nJets; ++iJet) {
      hJetPt.Fill(pt[iJet]);
      hJetEta.Fill(eta[iJet]);
      hJetPhi.Fill(phi[iJet]);  
      hJetArea.Fill(area[iJet]);
    }

    // -------------------------
    // 5) Dijet analysis (exactly 2 jets)
    // -------------------------

    if (nJets == 2) {
      double pt1 = pt[0];
      double pt2 = pt[1];

      // Dijet asymmetry A_J = (pt1 - pt2) / (pt1 + pt2)
      double asymmetry = (pt1 - pt2) / (pt1 + pt2);
      hAsymmetry.Fill(asymmetry);

      double phi1 = phi[0];
      double phi2 = phi[1];

      // calculate delta phi between the two jets and fill into histogram
      double deltaPhi = std::abs(phi1 - phi2);
      hDeltaPhi.Fill(deltaPhi);

    }

  }

  // -------------------------
  // 6) Save everything
  // -------------------------

  outFile.Write();
  outFile.Close();
}