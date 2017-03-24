// -*- C++ -*-
//
// Package:    L1TriggerTools/L1Analyzer
// Class:      L1Analyzer
// 
/**\class L1Analyzer L1Analyzer.cc L1TriggerTools/L1Analyzer/plugins/L1Analyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Devin Taylor
//         Created:  Wed, 25 Jan 2017 22:51:26 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TTree.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/L1TCalorimeter/interface/CaloTower.h"
#include "CalibFormats/CaloTPG/interface/CaloTPGTranscoder.h"
#include "CalibFormats/CaloTPG/interface/CaloTPGRecord.h"
#include "L1Trigger/L1TCaloLayer1/src/UCTGeometry.hh"
#include "L1Trigger/L1TCalorimeter/interface/CaloParamsHelper.h"
#include "CondFormats/L1TObjects/interface/CaloParams.h"
#include "CondFormats/DataRecord/interface/L1TCaloStage2ParamsRcd.h"
#include "DataFormats/L1Trigger/interface/L1Candidate.h"

#include "DataFormats/L1Trigger/interface/EGamma.h"
#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/L1Trigger/interface/Jet.h"
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/L1Trigger/interface/EtSum.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "DataFormats/Math/interface/deltaR.h"

//
// class declaration
//

class L1Analyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit L1Analyzer(const edm::ParameterSet&);
      ~L1Analyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      int getEcalPtBin(double pt);
      int getHcalPtBin(double pt);
      int getHFPtBin(double pt);

      double getEcalSF(double pt, int ieta);
      double getHcalSF(double pt, int ieta);
      double getHFSF(double pt, int ieta);

      void buildMatchToEcal(std::string name);
      void buildMatchToHcal(std::string name);

      void matchObjectToEcal(std::string name, const reco::Candidate& cand);
      void matchObjectToHcal(std::string name, const reco::Candidate& cand);
      void calculateSums(std::string name, int closestIEta, int closestIPhi);

      // ----------member data ---------------------------

      bool verbose_;
      edm::EDGetTokenT<EcalTrigPrimDigiCollection> ecalDigisToken_; 
      edm::EDGetTokenT<HcalTrigPrimDigiCollection> hcalDigisToken_;
      edm::EDGetTokenT<l1t::CaloTowerBxCollection> stage2Layer1DigisToken_;
      edm::EDGetTokenT<l1t::EGammaBxCollection> stage2EGToken_;
      edm::EDGetTokenT<l1t::TauBxCollection> stage2TauToken_;
      edm::EDGetTokenT<l1t::JetBxCollection> stage2JetToken_;
      //edm::EDGetTokenT<l1t::MuonBxCollection> stage2MuonToken_;
      edm::EDGetTokenT<l1t::EtSumBxCollection> stage2EtSumToken_;
      edm::EDGetTokenT<reco::CandidateView> electronsToken_;
      edm::EDGetTokenT<reco::CandidateView> photonsToken_;
      edm::EDGetTokenT<edm::View<pat::Jet>> jetsToken_;
      edm::EDGetTokenT<reco::CandidateView> metToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > electronIdMapToken_;
      edm::EDGetTokenT<reco::CandidateView> electronPairsToken_;
      edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
      bool storeEcal_;
      bool storeHcal_;
      bool storeStage2Layer1_;
      bool storeStage2_;
      double LSB_;
      bool isMC_;

      edm::Handle<EcalTrigPrimDigiCollection> ecalTPs_;
      edm::Handle<HcalTrigPrimDigiCollection> hcalTPs_;

      edm::ESHandle<CaloTPGTranscoder> decoder_;

      edm::Service<TFileService> fs_;
      TTree *tree_;

      std::vector<int> layer1ECalScaleETBins_;
      std::vector<double> layer1ECalScaleFactors_;
      std::vector<int> layer1HCalScaleETBins_;
      std::vector<double> layer1HCalScaleFactors_;
      std::vector<int> layer1HFScaleETBins_;
      std::vector<double> layer1HFScaleFactors_;

      UCTGeometry uctGeometry_;
      int N_IETA_ = l1tcalo::MaxCaloEta;
      int N_IPHI_ = l1tcalo::MaxCaloPhi;
      int N_ECAL_SF_ = 28;
      int N_ECAL_PT_ = 9;
      int N_HCAL_SF_ = 28;
      int N_HCAL_PT_ = 9;
      int N_HF_SF_ = 12;
      int N_HF_PT_ = 5;

      std::vector<std::vector<int> >    ecalTPs_compressedEt_;
      std::vector<std::vector<double> > ecalTPs_et_;
      std::vector<std::vector<double> > ecalTPs_correctedEt_;

      std::vector<std::vector<int> >    hcalTPs_compressedEt_;
      std::vector<std::vector<double> > hcalTPs_et_;
      std::vector<std::vector<double> > hcalTPs_correctedEt_;

      // branch variables
      Int_t run_;
      Int_t lumi_;
      ULong64_t event_;

      std::vector<int> ecalDigiEta_;
      std::vector<int> ecalDigiPhi_;
      std::vector<int> ecalDigiCompressedEt_;
      std::vector<bool> ecalDigiFineGrain_;
      std::vector<double> ecalDigiEt_;
      std::vector<double> ecalDigiCorrectedEt_;

      std::vector<int> hcalDigiEta_;
      std::vector<int> hcalDigiPhi_;
      std::vector<int> hcalDigiCompressedEt_;
      std::vector<bool> hcalDigiFineGrain0_;
      std::vector<bool> hcalDigiFineGrain1_;
      std::vector<bool> hcalDigiFineGrain2_;
      std::vector<bool> hcalDigiFineGrain3_;
      std::vector<int> hcalDigiVersion_;
      std::vector<double> hcalDigiEt_;
      std::vector<double> hcalDigiCorrectedEt_;

      std::vector<int> stage2Layer1DigiHwEtEm_;
      std::vector<int> stage2Layer1DigiHwEtHad_;
      std::vector<int> stage2Layer1DigiHwEtRatio_;
      std::vector<int> stage2Layer1DigiHwEta_;
      std::vector<int> stage2Layer1DigiHwPhi_;
      std::vector<int> stage2Layer1DigiHwPt_;
      std::vector<int> stage2Layer1DigiHwQual_;

      uint32_t expectedTotalET_;

      // EGamma
      std::vector<float> l1egPt_;
      std::vector<float> l1egEta_;
      std::vector<float> l1egPhi_;
      std::vector<int> l1egHwEta_;
      std::vector<int> l1egHwIso_;
      std::vector<int> l1egHwPhi_;
      std::vector<int> l1egHwPt_;
      std::vector<int> l1egHwQual_;
      std::vector<int> l1egFootprintEt_;
      std::vector<int> l1egIsoEt_;
      std::vector<int> l1egNTT_;
      std::vector<int> l1egRawEt_;
      std::vector<int> l1egShape_;
      std::vector<int> l1egTowerIEta_;
      std::vector<int> l1egTowerIPhi_;

      // Tau
      std::vector<float> l1tauPt_;
      std::vector<float> l1tauEta_;
      std::vector<float> l1tauPhi_;
      std::vector<int> l1tauHwEta_;
      std::vector<int> l1tauHwIso_;
      std::vector<int> l1tauHwPhi_;
      std::vector<int> l1tauHwPt_;
      std::vector<int> l1tauHwQual_;
      std::vector<int> l1tauHasEM_;
      std::vector<int> l1tauIsMerged_;
      std::vector<int> l1tauIsoEt_;
      std::vector<int> l1tauNTT_;
      std::vector<int> l1tauRawEt_;
      std::vector<int> l1tauTowerIEta_;
      std::vector<int> l1tauTowerIPhi_;

      // Jet
      std::vector<float> l1jetPt_;
      std::vector<float> l1jetEta_;
      std::vector<float> l1jetPhi_;
      std::vector<int> l1jetHwEta_;
      std::vector<int> l1jetHwIso_;
      std::vector<int> l1jetHwPhi_;
      std::vector<int> l1jetHwPt_;
      std::vector<int> l1jetHwQual_;
      std::vector<int> l1jetPUEt_;
      std::vector<int> l1jetRawEt_;
      std::vector<int> l1jetSeedEt_;
      std::vector<int> l1jetTowerIEta_;
      std::vector<int> l1jetTowerIPhi_;

      // EtSum
      std::vector<float> l1TotalEt_;
      std::vector<int> l1TotalEtHwPt_;
      std::vector<float> l1TotalHt_;
      std::vector<int> l1TotalHtHwPt_;
      std::vector<float> l1MissingEt_;
      std::vector<float> l1MissingEtPhi_;
      std::vector<int> l1MissingEtHwPhi_;
      std::vector<int> l1MissingEtHwPt_;

      // met
      std::vector<float> metEt_;
      std::vector<float> metPhi_;

      std::map<std::string, std::vector<float> > branches_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
L1Analyzer::L1Analyzer(const edm::ParameterSet& iConfig):
  ecalDigisToken_(consumes<EcalTrigPrimDigiCollection>(iConfig.getParameter<edm::InputTag>("ecalDigis"))),
  hcalDigisToken_(consumes<HcalTrigPrimDigiCollection>(iConfig.getParameter<edm::InputTag>("hcalDigis"))),
  stage2Layer1DigisToken_(consumes<l1t::CaloTowerBxCollection>(iConfig.getParameter<edm::InputTag>("stage2Layer1Digis"))),
  stage2EGToken_(consumes<l1t::EGammaBxCollection>(iConfig.getParameter<edm::InputTag>("stage2EG"))),
  stage2TauToken_(consumes<l1t::TauBxCollection>(iConfig.getParameter<edm::InputTag>("stage2Tau"))),
  stage2JetToken_(consumes<l1t::JetBxCollection>(iConfig.getParameter<edm::InputTag>("stage2Jet"))),
  //stage2MuonToken_(consumes<l1t::MuonBxCollection>(iConfig.getParameter<edm::InputTag>("stage2Muon"))),
  stage2EtSumToken_(consumes<l1t::EtSumBxCollection>(iConfig.getParameter<edm::InputTag>("stage2EtSum"))),
  electronsToken_(consumes<reco::CandidateView>(iConfig.getParameter<edm::InputTag>("electrons"))),
  photonsToken_(consumes<reco::CandidateView>(iConfig.getParameter<edm::InputTag>("photons"))),
  jetsToken_(consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jets"))),
  metToken_(consumes<reco::CandidateView>(iConfig.getParameter<edm::InputTag>("met"))),
  electronIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronIdMap"))),
  electronPairsToken_(consumes<reco::CandidateView>(iConfig.getParameter<edm::InputTag>("electronPairs"))),
  genParticlesToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"))),
  storeEcal_(iConfig.exists("storeEcal") ? iConfig.getParameter<bool>("storeEcal") : true),
  storeHcal_(iConfig.exists("storeHcal") ? iConfig.getParameter<bool>("storeHcal") : true),
  storeStage2Layer1_(iConfig.exists("storeStage2Layer1") ? iConfig.getParameter<bool>("storeStage2Layer1") : true),
  storeStage2_(iConfig.exists("storeStage2") ? iConfig.getParameter<bool>("storeStage2") : true),
  LSB_(iConfig.exists("LSB") ? iConfig.getParameter<double>("LSB") : 0.5),
  isMC_(iConfig.getParameter<bool>("isMC"))
{
  // intialize tfileservice
  usesResource("TFileService");

  // create tree_
  tree_ = fs_->make<TTree>("L1TTree", "L1TTree");

  // initialize branches
  tree_->Branch("run", &run_, "run/I");
  tree_->Branch("lumi", &lumi_, "lumi/I");
  tree_->Branch("event", &event_, "event/l");
  
  // ecalTPs
  if (storeEcal_) {
    tree_->Branch("ecalDigi_ieta", &ecalDigiEta_);
    tree_->Branch("ecalDigi_iphi", &ecalDigiPhi_);
    tree_->Branch("ecalDigi_compressedEt", &ecalDigiCompressedEt_);
    tree_->Branch("ecalDigi_fineGrain", &ecalDigiFineGrain_);
    tree_->Branch("ecalDigi_et", &ecalDigiEt_);
    tree_->Branch("ecalDigi_correctedEt", &ecalDigiCorrectedEt_);
  }

  // hcalTPs
  if (storeHcal_) {
    tree_->Branch("hcalDigi_ieta", &hcalDigiEta_);
    tree_->Branch("hcalDigi_iphi", &hcalDigiPhi_);
    tree_->Branch("hcalDigi_compressedEt", &hcalDigiCompressedEt_);
    tree_->Branch("hcalDigi_fineGrain0", &hcalDigiFineGrain0_);
    tree_->Branch("hcalDigi_fineGrain1", &hcalDigiFineGrain1_);
    tree_->Branch("hcalDigi_fineGrain2", &hcalDigiFineGrain2_);
    tree_->Branch("hcalDigi_fineGrain3", &hcalDigiFineGrain3_);
    tree_->Branch("hcalDigi_version", &hcalDigiVersion_);
    tree_->Branch("hcalDigi_et", &hcalDigiEt_);
    tree_->Branch("hcalDigi_correctedEt", &hcalDigiCorrectedEt_);
  }

  // stage2layer1 digis
  if (storeStage2Layer1_) {
    tree_->Branch("stage2Layer1Digi_hwEtEm", &stage2Layer1DigiHwEtEm_);
    tree_->Branch("stage2Layer1Digi_hwEtHad", &stage2Layer1DigiHwEtHad_);
    tree_->Branch("stage2Layer1Digi_hwEtRatio", &stage2Layer1DigiHwEtRatio_);
    tree_->Branch("stage2Layer1Digi_hwEta", &stage2Layer1DigiHwEta_);
    tree_->Branch("stage2Layer1Digi_hwPhi", &stage2Layer1DigiHwPhi_);
    tree_->Branch("stage2Layer1Digi_hwPt", &stage2Layer1DigiHwPt_);
    tree_->Branch("stage2Layer1Digi_hwQual", &stage2Layer1DigiHwQual_);

    // sums
    tree_->Branch("expectedTotalET", &expectedTotalET_);
  }

  if (storeStage2_) {
    // EGamma
    tree_->Branch("l1eg_pt",          &l1egPt_);
    tree_->Branch("l1eg_eta",         &l1egEta_);
    tree_->Branch("l1eg_phi",         &l1egPhi_);
    tree_->Branch("l1eg_hwEta",       &l1egHwEta_);
    tree_->Branch("l1eg_hwIso",       &l1egHwIso_);
    tree_->Branch("l1eg_hwPhi",       &l1egHwPhi_);
    tree_->Branch("l1eg_hwPt",        &l1egHwPt_);
    tree_->Branch("l1eg_hwQual",      &l1egHwQual_);
    tree_->Branch("l1eg_footprintEt", &l1egFootprintEt_);
    tree_->Branch("l1eg_isoEt",       &l1egIsoEt_);
    tree_->Branch("l1eg_nTT",         &l1egNTT_);
    tree_->Branch("l1eg_rawEt",       &l1egRawEt_);
    tree_->Branch("l1eg_shape",       &l1egShape_);
    tree_->Branch("l1eg_towerIEta",   &l1egTowerIEta_);
    tree_->Branch("l1eg_towerIPhi",   &l1egTowerIPhi_);

    // Tau
    tree_->Branch("l1tau_pt",        &l1tauPt_);
    tree_->Branch("l1tau_eta",       &l1tauEta_);
    tree_->Branch("l1tau_phi",       &l1tauPhi_);
    tree_->Branch("l1tau_hwEta",     &l1tauHwEta_);
    tree_->Branch("l1tau_hwIso",     &l1tauHwIso_);
    tree_->Branch("l1tau_hwPhi",     &l1tauHwPhi_);
    tree_->Branch("l1tau_hwPt",      &l1tauHwPt_);
    tree_->Branch("l1tau_hwQual",    &l1tauHwQual_);
    tree_->Branch("l1tau_hasEM",     &l1tauHasEM_);
    tree_->Branch("l1tau_isMerged",  &l1tauIsMerged_);
    tree_->Branch("l1tau_isoEt",     &l1tauIsoEt_);
    tree_->Branch("l1tau_nTT",       &l1tauNTT_);
    tree_->Branch("l1tau_rawEt",     &l1tauRawEt_);
    tree_->Branch("l1tau_towerIEta", &l1tauTowerIEta_);
    tree_->Branch("l1tau_towerIPhi", &l1tauTowerIPhi_);

    // Jet
    tree_->Branch("l1jet_pt",        &l1jetPt_);
    tree_->Branch("l1jet_eta",       &l1jetEta_);
    tree_->Branch("l1jet_phi",       &l1jetPhi_);
    tree_->Branch("l1jet_hwEta",     &l1jetHwEta_);
    tree_->Branch("l1jet_hwIso",     &l1jetHwIso_);
    tree_->Branch("l1jet_hwPhi",     &l1jetHwPhi_);
    tree_->Branch("l1jet_hwPt",      &l1jetHwPt_);
    tree_->Branch("l1jet_hwQual",    &l1jetHwQual_);
    tree_->Branch("l1jet_puEt",      &l1jetPUEt_);
    tree_->Branch("l1jet_rawEt",     &l1jetRawEt_);
    tree_->Branch("l1jet_seedEt",    &l1jetSeedEt_);
    tree_->Branch("l1jet_towerIEta", &l1jetTowerIEta_);
    tree_->Branch("l1jet_towerIPhi", &l1jetTowerIPhi_);

    // EtSum
    tree_->Branch("l1TotalEt_et",       &l1TotalEt_);
    tree_->Branch("l1TotalEt_hwPt",     &l1TotalEtHwPt_);
    tree_->Branch("l1TotalHt_et",       &l1TotalHt_);
    tree_->Branch("l1TotalHt_hwPt",     &l1TotalHtHwPt_);
    tree_->Branch("l1MissingEt_et",     &l1MissingEt_);
    tree_->Branch("l1MissingEt_phi",    &l1MissingEtPhi_);
    tree_->Branch("l1MissingEt_hwPhi",  &l1MissingEtHwPhi_);
    tree_->Branch("l1MissingEt_hwPt",   &l1MissingEtHwPt_);
  }

  // electrons
  buildMatchToEcal("electron");
  buildMatchToEcal("electronFromZ");

  // jets
  buildMatchToHcal("jet");

  // pions
  if (isMC_) buildMatchToHcal("genPion");

  // met
  tree_->Branch("met_et",  &metEt_);
  tree_->Branch("met_phi", &metPhi_);
}


L1Analyzer::~L1Analyzer()
{
 
}


//
// member functions
//

int
L1Analyzer::getEcalPtBin(double pt) {
  int p=0;
  for (p=0; p<N_ECAL_PT_; ++p) {
    if (pt<layer1ECalScaleETBins_[p]) break;
  }
  return p;
}

int
L1Analyzer::getHcalPtBin(double pt) {
  int p=0;
  for (p=0; p<N_HCAL_PT_; ++p) {
    if (pt<layer1HCalScaleETBins_[p]) break;
  }
  return p;
}

int
L1Analyzer::getHFPtBin(double pt) {
  int p=0;
  for (p=0; p<N_HF_PT_; ++p) {
    if (pt<layer1HFScaleETBins_[p]) break;
  }
  return p;
}

double
L1Analyzer::getEcalSF(double pt, int ieta) {
  if (fabs(ieta)>N_ECAL_SF_) return 1.;     // goes to 1-28 +/-
  int pbin = getEcalPtBin(pt);
  int sfbin = pbin*N_ECAL_SF_+fabs(ieta)-1; // move 1-28 to 0-27
  return layer1ECalScaleFactors_[sfbin]; 
}

double
L1Analyzer::getHcalSF(double pt, int ieta) {
  if (fabs(ieta)==29) return 1.;                       // skip 29
  if (fabs(ieta)>N_HCAL_SF_) return getHFSF(pt, ieta); // goes to 1-28 +/-, if outside try HF
  int pbin = getHcalPtBin(pt);
  int sfbin = pbin*N_HCAL_SF_+fabs(ieta)-1;            // move 1-28 to 0-27
  return layer1HCalScaleFactors_[sfbin]; 
}

double
L1Analyzer::getHFSF(double pt, int ieta) {
  if (fabs(ieta)>N_HCAL_SF_+N_HF_SF_+1) return 1.;   // goes from 30-41 (skip 29) +/-
  int pbin = getHFPtBin(pt);
  int sfbin = pbin*N_HF_SF_+fabs(ieta)-N_HCAL_SF_-2; // move 30-41 to 0-11
  return layer1HFScaleFactors_[sfbin]; 
}

void
L1Analyzer::buildMatchToEcal(std::string name) {
  std::vector<std::string> vars;
  vars.push_back("pt");
  vars.push_back("eta");
  vars.push_back("phi");
  vars.push_back("pdgId");
  vars.push_back("matchedDR");
  vars.push_back("matchedDEt");
  vars.push_back("matchedIEta");
  vars.push_back("matchedIPhi");
  vars.push_back("matchedEt");
  vars.push_back("matchedCorrectedEt");
  vars.push_back("matchedPtBin");
  vars.push_back("matchedEcalSum3x3");
  vars.push_back("matchedEcalSum5x5");
  vars.push_back("matchedEcalSum7x7");
  vars.push_back("matchedEcalCorrectedSum3x3");
  vars.push_back("matchedEcalCorrectedSum5x5");
  vars.push_back("matchedEcalCorrectedSum7x7");
  vars.push_back("matchedHcalSum3x3");
  vars.push_back("matchedHcalSum5x5");
  vars.push_back("matchedHcalSum7x7");
  vars.push_back("matchedHcalCorrectedSum3x3");
  vars.push_back("matchedHcalCorrectedSum5x5");
  vars.push_back("matchedHcalCorrectedSum7x7");
  vars.push_back("matchedSum3x3");
  vars.push_back("matchedSum5x5");
  vars.push_back("matchedSum7x7");
  vars.push_back("matchedCorrectedSum3x3");
  vars.push_back("matchedCorrectedSum5x5");
  vars.push_back("matchedCorrectedSum7x7");

  for (auto var: vars) {
    std::string bname = name+"_"+var;
    branches_[bname] = std::vector<float>();
    tree_->Branch(bname.c_str(), &branches_[bname]);
  }
}

void
L1Analyzer::buildMatchToHcal(std::string name) {
  buildMatchToEcal(name);
}

void
L1Analyzer::matchObjectToEcal(std::string name, const reco::Candidate& cand) {
  float pt = cand.pt();
  if (pt==0) return;
  float eta = cand.eta();
  float phi = cand.phi();
  int pdgId = cand.pdgId();

  // find closest max ecaltp
  double closestDR = 1e6;
  double closestDEt = 1e6;
  double closestEt = 0;
  double closestCorrectedEt = 0;
  int closestIEta = -999;
  int closestIPhi = -999;
  int match = -1;
  for ( const auto& ecalTP : *ecalTPs_ ) {
    int tp_ieta = ecalTP.id().ieta();
    int tp_iphi = ecalTP.id().iphi();
    int tp_compressedEt = ecalTP.compressedEt();
    double tp_eta = uctGeometry_.getUCTTowerEta(tp_ieta);
    double tp_phi = uctGeometry_.getUCTTowerPhi(tp_iphi);
    double tp_et = tp_compressedEt*LSB_;
    double correctedEt = tp_et*getEcalSF(tp_et,tp_ieta);
    double dr = deltaR(tp_eta,tp_phi,eta,phi);
    double det = fabs(pt-tp_et);
    if (dr<0.2 && dr<closestDR && pt<127 && det<closestDEt) {
      closestDR = dr;
      closestDEt = det;
      closestIEta = tp_ieta;
      closestIPhi = tp_iphi;
      closestEt = tp_et;
      closestCorrectedEt = correctedEt;
      match = 1;
    }
  }

  int closestPtBin = getHcalPtBin(closestEt);

  branches_[name+"_pt"].push_back(pt);
  branches_[name+"_eta"].push_back(eta);
  branches_[name+"_phi"].push_back(phi);
  branches_[name+"_pdgId"].push_back(pdgId);

  // if match found and electron is good
  if (match>0) {
    calculateSums(name,closestIEta,closestIPhi);

    branches_[name+"_matchedDR"].push_back(closestDR);
    branches_[name+"_matchedDEt"].push_back(closestDEt);
    branches_[name+"_matchedIEta"].push_back(closestIEta);
    branches_[name+"_matchedIPhi"].push_back(closestIPhi);
    branches_[name+"_matchedEt"].push_back(closestEt);
    branches_[name+"_matchedCorrectedEt"].push_back(closestCorrectedEt);
    branches_[name+"_matchedPtBin"].push_back(closestPtBin);
  }
  else {
    branches_[name+"_matchedEcalSum3x3"].push_back(-1);
    branches_[name+"_matchedEcalSum5x5"].push_back(-1);
    branches_[name+"_matchedEcalSum7x7"].push_back(-1);
    branches_[name+"_matchedEcalCorrectedSum3x3"].push_back(-1);
    branches_[name+"_matchedEcalCorrectedSum5x5"].push_back(-1);
    branches_[name+"_matchedEcalCorrectedSum7x7"].push_back(-1);
    branches_[name+"_matchedHcalSum3x3"].push_back(-1);
    branches_[name+"_matchedHcalSum5x5"].push_back(-1);
    branches_[name+"_matchedHcalSum7x7"].push_back(-1);
    branches_[name+"_matchedHcalCorrectedSum3x3"].push_back(-1);
    branches_[name+"_matchedHcalCorrectedSum5x5"].push_back(-1);
    branches_[name+"_matchedHcalCorrectedSum7x7"].push_back(-1);
    branches_[name+"_matchedSum3x3"].push_back(-1);
    branches_[name+"_matchedSum5x5"].push_back(-1);
    branches_[name+"_matchedSum7x7"].push_back(-1);
    branches_[name+"_matchedCorrectedSum3x3"].push_back(-1);
    branches_[name+"_matchedCorrectedSum5x5"].push_back(-1);
    branches_[name+"_matchedCorrectedSum7x7"].push_back(-1);

    branches_[name+"_matchedDR"].push_back(-1);
    branches_[name+"_matchedDEt"].push_back(-1);
    branches_[name+"_matchedIEta"].push_back(-1);
    branches_[name+"_matchedIPhi"].push_back(-1);
    branches_[name+"_matchedEt"].push_back(-1);
    branches_[name+"_matchedCorrectedEt"].push_back(-1);
    branches_[name+"_matchedPtBin"].push_back(-1);
  }
}

void
L1Analyzer::matchObjectToHcal(std::string name, const reco::Candidate& cand) {
  float pt = cand.pt();
  if (pt==0) return;
  float eta = cand.eta();
  float phi = cand.phi();
  int pdgId = cand.pdgId();

  // find closest max hcaltp
  double closestDR = 1e6;
  double closestDEt = 1e6;
  double closestEt = 0;
  double closestCorrectedEt = 0;
  int closestIEta = -999;
  int closestIPhi = -999;
  int match = -1;
  for ( const auto& hcalTP : *hcalTPs_ ) {
    int tp_ieta = hcalTP.id().ieta();
    int tp_iphi = hcalTP.id().iphi();
    double tp_eta = uctGeometry_.getUCTTowerEta(tp_ieta);
    double tp_phi = uctGeometry_.getUCTTowerPhi(tp_iphi);
    double tp_et = decoder_->hcaletValue(hcalTP.id(), hcalTP.t0());
    double correctedEt = tp_et*getHcalSF(tp_et,tp_ieta);
    double dr = deltaR(tp_eta,tp_phi,eta,phi);
    double det = fabs(pt-tp_et);
    if (abs(tp_ieta)>20) det = fabs(pt-2*tp_et);
    if (dr<0.2 && det<closestDEt && pt<127) {
      closestDR = dr;
      closestDEt = det;
      closestIEta = tp_ieta;
      closestIPhi = tp_iphi;
      closestEt = tp_et;
      closestCorrectedEt = correctedEt;
      match = 1;
    }
  }

  int closestPtBin = getHcalPtBin(closestEt);

  branches_[name+"_pt"].push_back(pt);
  branches_[name+"_eta"].push_back(eta);
  branches_[name+"_phi"].push_back(phi);
  branches_[name+"_pdgId"].push_back(pdgId);

  // if match found
  if (match>0) {
    calculateSums(name,closestIEta,closestIPhi);

    //std::cout << name << " dr " << closestDR << " det " << closestDEt << " ieta " << closestIEta << " iphi " << closestIPhi << " et " << closestEt << " corrEt " << closestCorrectedEt << " ptb " << closestPtBin << std::endl;

    branches_[name+"_matchedDR"].push_back(closestDR);
    branches_[name+"_matchedDEt"].push_back(closestDEt);
    branches_[name+"_matchedIEta"].push_back(closestIEta);
    branches_[name+"_matchedIPhi"].push_back(closestIPhi);
    branches_[name+"_matchedEt"].push_back(closestEt);
    branches_[name+"_matchedCorrectedEt"].push_back(closestCorrectedEt);
    branches_[name+"_matchedPtBin"].push_back(closestPtBin);
  }
  else {
    branches_[name+"_matchedEcalSum3x3"].push_back(-1);
    branches_[name+"_matchedEcalSum5x5"].push_back(-1);
    branches_[name+"_matchedEcalSum7x7"].push_back(-1);
    branches_[name+"_matchedEcalCorrectedSum3x3"].push_back(-1);
    branches_[name+"_matchedEcalCorrectedSum5x5"].push_back(-1);
    branches_[name+"_matchedEcalCorrectedSum7x7"].push_back(-1);
    branches_[name+"_matchedHcalSum3x3"].push_back(-1);
    branches_[name+"_matchedHcalSum5x5"].push_back(-1);
    branches_[name+"_matchedHcalSum7x7"].push_back(-1);
    branches_[name+"_matchedHcalCorrectedSum3x3"].push_back(-1);
    branches_[name+"_matchedHcalCorrectedSum5x5"].push_back(-1);
    branches_[name+"_matchedHcalCorrectedSum7x7"].push_back(-1);
    branches_[name+"_matchedSum3x3"].push_back(-1);
    branches_[name+"_matchedSum5x5"].push_back(-1);
    branches_[name+"_matchedSum7x7"].push_back(-1);
    branches_[name+"_matchedCorrectedSum3x3"].push_back(-1);
    branches_[name+"_matchedCorrectedSum5x5"].push_back(-1);
    branches_[name+"_matchedCorrectedSum7x7"].push_back(-1);

    branches_[name+"_matchedDR"].push_back(-1);
    branches_[name+"_matchedDEt"].push_back(-1);
    branches_[name+"_matchedIEta"].push_back(-1);
    branches_[name+"_matchedIPhi"].push_back(-1);
    branches_[name+"_matchedEt"].push_back(-1);
    branches_[name+"_matchedCorrectedEt"].push_back(-1);
    branches_[name+"_matchedPtBin"].push_back(-1);
  }
}

void
L1Analyzer::calculateSums(std::string name, int closestIEta, int closestIPhi) {
  // calculate sums
  double sum3x3 = 0.;
  double sum5x5 = 0.;
  double sum7x7 = 0.;
  double correctedSum3x3 = 0.;
  double correctedSum5x5 = 0.;
  double correctedSum7x7 = 0.;
  double ecalSum3x3 = 0.;
  double ecalSum5x5 = 0.;
  double ecalSum7x7 = 0.;
  double ecalCorrectedSum3x3 = 0.;
  double ecalCorrectedSum5x5 = 0.;
  double ecalCorrectedSum7x7 = 0.;
  double hcalSum3x3 = 0.;
  double hcalSum5x5 = 0.;
  double hcalSum7x7 = 0.;
  double hcalCorrectedSum3x3 = 0.;
  double hcalCorrectedSum5x5 = 0.;
  double hcalCorrectedSum7x7 = 0.;
  for (int e=-3; e<4; ++e) {
    int ie = N_IETA_+closestIEta;
    if (ie>=N_IETA_) ie-=1;                 // -41 .. -1 1 .. 41 -> 0 .. 81
    ie+=e;
    if (ie<0 || ie>=2*N_IETA_) continue;    // max tower eta 41
    for (int p=-3; p<4; ++p) {
      int ip = closestIPhi-1;               // 1 .. 72 -> 0 .. 71
      ip+=p;
      if (ip<0) ip+=N_IPHI_;                // wrap -1 -> 71, etc
      if (ip>=N_IPHI_) ip-=N_IPHI_;         // wrap 72 -> 0, etc
      if (e>-2 && e<2 && p>-2 && p<2) {
        sum3x3 += ecalTPs_et_[ie][ip];
        correctedSum3x3 += ecalTPs_correctedEt_[ie][ip];
        ecalSum3x3 += ecalTPs_et_[ie][ip];
        ecalCorrectedSum3x3 += ecalTPs_correctedEt_[ie][ip];
        sum3x3 += hcalTPs_et_[ie][ip];
        correctedSum3x3 += hcalTPs_correctedEt_[ie][ip];
        hcalSum3x3 += hcalTPs_et_[ie][ip];
        hcalCorrectedSum3x3 += hcalTPs_correctedEt_[ie][ip];
      }
      if (e>-3 && e<3 && p>-3 && p<3) {
        sum5x5 += ecalTPs_et_[ie][ip];
        correctedSum5x5 += ecalTPs_correctedEt_[ie][ip];
        ecalSum5x5 += ecalTPs_et_[ie][ip];
        ecalCorrectedSum5x5 += ecalTPs_correctedEt_[ie][ip];
        sum5x5 += hcalTPs_et_[ie][ip];
        correctedSum5x5 += hcalTPs_correctedEt_[ie][ip];
        hcalSum5x5 += hcalTPs_et_[ie][ip];
        hcalCorrectedSum5x5 += hcalTPs_correctedEt_[ie][ip];
      }
      sum7x7 += ecalTPs_et_[ie][ip];
      correctedSum7x7 += ecalTPs_correctedEt_[ie][ip];
      ecalSum7x7 += ecalTPs_et_[ie][ip];
      ecalCorrectedSum7x7 += ecalTPs_correctedEt_[ie][ip];
      sum7x7 += hcalTPs_et_[ie][ip];
      correctedSum7x7 += hcalTPs_correctedEt_[ie][ip];
      hcalSum7x7 += hcalTPs_et_[ie][ip];
      hcalCorrectedSum7x7 += hcalTPs_correctedEt_[ie][ip];
    }
  }

  branches_[name+"_matchedEcalSum3x3"].push_back(ecalSum3x3);
  branches_[name+"_matchedEcalSum5x5"].push_back(ecalSum5x5);
  branches_[name+"_matchedEcalSum7x7"].push_back(ecalSum7x7);
  branches_[name+"_matchedEcalCorrectedSum3x3"].push_back(ecalCorrectedSum3x3);
  branches_[name+"_matchedEcalCorrectedSum5x5"].push_back(ecalCorrectedSum5x5);
  branches_[name+"_matchedEcalCorrectedSum7x7"].push_back(ecalCorrectedSum7x7);
  branches_[name+"_matchedHcalSum3x3"].push_back(hcalSum3x3);
  branches_[name+"_matchedHcalSum5x5"].push_back(hcalSum5x5);
  branches_[name+"_matchedHcalSum7x7"].push_back(hcalSum7x7);
  branches_[name+"_matchedHcalCorrectedSum3x3"].push_back(hcalCorrectedSum3x3);
  branches_[name+"_matchedHcalCorrectedSum5x5"].push_back(hcalCorrectedSum5x5);
  branches_[name+"_matchedHcalCorrectedSum7x7"].push_back(hcalCorrectedSum7x7);
  branches_[name+"_matchedSum3x3"].push_back(sum3x3);
  branches_[name+"_matchedSum5x5"].push_back(sum5x5);
  branches_[name+"_matchedSum7x7"].push_back(sum7x7);
  branches_[name+"_matchedCorrectedSum3x3"].push_back(correctedSum3x3);
  branches_[name+"_matchedCorrectedSum5x5"].push_back(correctedSum5x5);
  branches_[name+"_matchedCorrectedSum7x7"].push_back(correctedSum7x7);
}

// ------------ method called for each event  ------------
void
L1Analyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  event_ = iEvent.id().event();
  run_   = iEvent.run();
  lumi_  = iEvent.luminosityBlock();

  // grab collections
  iEvent.getByToken(ecalDigisToken_, ecalTPs_);

  iEvent.getByToken(hcalDigisToken_, hcalTPs_);

  edm::Handle<l1t::CaloTowerBxCollection> caloTowerBXs;
  iEvent.getByToken(stage2Layer1DigisToken_, caloTowerBXs);

  edm::Handle<l1t::EGammaBxCollection> l1egs;
  iEvent.getByToken(stage2EGToken_, l1egs);

  edm::Handle<l1t::TauBxCollection> l1taus;
  iEvent.getByToken(stage2TauToken_, l1taus);

  edm::Handle<l1t::JetBxCollection> l1jets;
  iEvent.getByToken(stage2JetToken_, l1jets);

  edm::Handle<l1t::EtSumBxCollection> l1sums;
  iEvent.getByToken(stage2EtSumToken_, l1sums);

  edm::Handle<reco::CandidateView> electrons;
  iEvent.getByToken(electronsToken_, electrons);

  edm::Handle<reco::CandidateView> photons;
  iEvent.getByToken(photonsToken_, photons);

  edm::Handle<edm::View<pat::Jet>> jets;
  iEvent.getByToken(jetsToken_, jets);

  edm::Handle<reco::CandidateView> met;
  iEvent.getByToken(metToken_, met);

  edm::Handle<edm::ValueMap<bool> > electronIdMap;
  iEvent.getByToken(electronIdMapToken_, electronIdMap);

  edm::Handle<reco::CandidateView> electronPairs;
  iEvent.getByToken(electronPairsToken_, electronPairs);

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticlesToken_, genParticles);

  iSetup.get<CaloTPGRecord>().get(decoder_);

  edm::ESHandle<l1t::CaloParams> paramsHandle;
  iSetup.get<L1TCaloStage2ParamsRcd>().get(paramsHandle);
  l1t::CaloParamsHelper caloParams(*paramsHandle.product());

  layer1ECalScaleETBins_ = caloParams.layer1ECalScaleETBins();
  layer1ECalScaleFactors_ = caloParams.layer1ECalScaleFactors();
  layer1HCalScaleETBins_ = caloParams.layer1HCalScaleETBins();
  layer1HCalScaleFactors_ = caloParams.layer1HCalScaleFactors();
  layer1HFScaleETBins_ = caloParams.layer1HFScaleETBins();
  layer1HFScaleFactors_ = caloParams.layer1HFScaleFactors();

  N_ECAL_PT_ = layer1ECalScaleETBins_.size();
  N_HCAL_PT_ = layer1HCalScaleETBins_.size();
  N_HF_PT_ = layer1HFScaleETBins_.size();

  // clear branches
  ecalDigiEta_.clear();
  ecalDigiPhi_.clear();
  ecalDigiCompressedEt_.clear();
  ecalDigiFineGrain_.clear();
  ecalDigiEt_.clear();
  ecalDigiCorrectedEt_.clear();

  hcalDigiEta_.clear();
  hcalDigiPhi_.clear();
  hcalDigiCompressedEt_.clear();
  hcalDigiFineGrain0_.clear();
  hcalDigiFineGrain1_.clear();
  hcalDigiFineGrain2_.clear();
  hcalDigiFineGrain3_.clear();
  hcalDigiVersion_.clear();
  hcalDigiEt_.clear();
  hcalDigiCorrectedEt_.clear();

  stage2Layer1DigiHwEtEm_.clear();
  stage2Layer1DigiHwEtHad_.clear();
  stage2Layer1DigiHwEtRatio_.clear();
  stage2Layer1DigiHwEta_.clear();
  stage2Layer1DigiHwPhi_.clear();
  stage2Layer1DigiHwPt_.clear();
  stage2Layer1DigiHwQual_.clear();

  expectedTotalET_ = 0;

  l1egPt_.clear();
  l1egEta_.clear();
  l1egPhi_.clear();
  l1egHwEta_.clear();
  l1egHwIso_.clear();
  l1egHwPhi_.clear();
  l1egHwPt_.clear();
  l1egHwQual_.clear();
  l1egFootprintEt_.clear();
  l1egIsoEt_.clear();
  l1egNTT_.clear();
  l1egRawEt_.clear();
  l1egShape_.clear();
  l1egTowerIEta_.clear();
  l1egTowerIPhi_.clear();

  l1tauPt_.clear();
  l1tauEta_.clear();
  l1tauPhi_.clear();
  l1tauHwEta_.clear();
  l1tauHwIso_.clear();
  l1tauHwPhi_.clear();
  l1tauHwPt_.clear();
  l1tauHwQual_.clear();
  l1tauHasEM_.clear();
  l1tauIsMerged_.clear();
  l1tauIsoEt_.clear();
  l1tauNTT_.clear();
  l1tauRawEt_.clear();
  l1tauTowerIEta_.clear();
  l1tauTowerIPhi_.clear();

  l1jetPt_.clear();
  l1jetEta_.clear();
  l1jetPhi_.clear();
  l1jetHwEta_.clear();
  l1jetHwIso_.clear();
  l1jetHwPhi_.clear();
  l1jetHwPt_.clear();
  l1jetHwQual_.clear();
  l1jetPUEt_.clear();
  l1jetRawEt_.clear();
  l1jetSeedEt_.clear();
  l1jetTowerIEta_.clear();
  l1jetTowerIPhi_.clear();

  l1TotalEt_.clear();
  l1TotalEtHwPt_.clear();
  l1TotalHt_.clear();
  l1TotalHtHwPt_.clear();
  l1MissingEt_.clear();
  l1MissingEtPhi_.clear();
  l1MissingEtHwPhi_.clear();
  l1MissingEtHwPt_.clear();

  metEt_.clear();
  metPhi_.clear();

  for (auto& entry: branches_) {
    entry.second.clear();
  }

  // build objects to store TP energy
  ecalTPs_compressedEt_.clear();
  ecalTPs_et_.clear();
  ecalTPs_correctedEt_.clear();
  hcalTPs_compressedEt_.clear();
  hcalTPs_et_.clear();
  hcalTPs_correctedEt_.clear();

  for (int ie=0; ie<N_IETA_*2; ++ie) {
    ecalTPs_compressedEt_.push_back(std::vector<int>(N_IPHI_));
    ecalTPs_et_.push_back(std::vector<double>(N_IPHI_));
    ecalTPs_correctedEt_.push_back(std::vector<double>(N_IPHI_));
    hcalTPs_compressedEt_.push_back(std::vector<int>(N_IPHI_));
    hcalTPs_et_.push_back(std::vector<double>(N_IPHI_));
    hcalTPs_correctedEt_.push_back(std::vector<double>(N_IPHI_));
  }

  // set branch values

  // ecalTPs
  for ( const auto& ecalTp : *ecalTPs_ ) {
    int caloEta = ecalTp.id().ieta();
    int caloPhi = ecalTp.id().iphi();
    int compressedEt = ecalTp.compressedEt();
    double et = compressedEt*LSB_;
    bool fgVeto = ecalTp.fineGrain();
    int ie = caloEta>0 ? N_IETA_+caloEta-1 : N_IETA_+caloEta;
    int ip = caloPhi-1;
    double correctedEt = et*getEcalSF(et,caloEta);
    ecalTPs_compressedEt_[ie][ip] = compressedEt;
    ecalTPs_et_[ie][ip] = et;
    ecalTPs_correctedEt_[ie][ip] = correctedEt;
    if (compressedEt>0) {
      ecalDigiEta_.push_back(caloEta);
      ecalDigiPhi_.push_back(caloPhi);
      ecalDigiCompressedEt_.push_back(compressedEt);
      ecalDigiFineGrain_.push_back(fgVeto);
      ecalDigiEt_.push_back(et);
      ecalDigiCorrectedEt_.push_back(correctedEt);
    }
    expectedTotalET_ += compressedEt;
  }

  // hcalTPs
  for ( const auto& hcalTp : *hcalTPs_ ) {
    int caloEta = hcalTp.id().ieta();
    uint32_t absCaloEta = abs(caloEta);
    int hcalVersion = hcalTp.id().version();
    // Tower 29 is not used by Layer-1
    if(absCaloEta == 29) {
      continue;
    }
    // Prevent usage of HF TPs with Layer-1 emulator if HCAL TPs are old style
    //else if(hcalTp.id().version() == 0 && absCaloEta > 29) {
    //  continue;
    //}
    else if(absCaloEta <= 41) {
      int caloPhi = hcalTp.id().iphi();
      int compressedEt = hcalTp.SOI_compressedEt();
      double et = decoder_->hcaletValue(hcalTp.id(), hcalTp.t0());
      bool fg  = hcalTp.t0().fineGrain(0);
      bool fg2 = hcalTp.t0().fineGrain(1);
      bool fg3 = hcalTp.t0().fineGrain(2);
      bool fg4 = hcalTp.t0().fineGrain(3);
      int ie = caloEta>0 ? N_IETA_+caloEta-1 : N_IETA_+caloEta;
      int ip = caloPhi-1;
      double correctedEt = et*getHcalSF(et,caloEta);
      hcalTPs_compressedEt_[ie][ip] = compressedEt;
      hcalTPs_et_[ie][ip] = et;
      hcalTPs_correctedEt_[ie][ip] = correctedEt;
      if(caloPhi <= 72) {
        uint32_t featureBits = 0;
        if(fg)  featureBits |= 0b01;
        // fg2 should only be set for HF
        if(absCaloEta > 29 && fg2) featureBits |= 0b10;
        if (compressedEt>0) {
          if (compressedEt==10 && caloPhi>61 && caloPhi<64 && abs(caloEta)>16 && abs(caloEta)<29) {
            //std::cout << "iphi: " << caloPhi << " ieta: " << caloEta << " compressedEt: " << compressedEt << " et: " << et << " correctedEt: " << correctedEt << std::endl;
          }
          hcalDigiEta_.push_back(caloEta);
          hcalDigiPhi_.push_back(caloPhi);
          hcalDigiCompressedEt_.push_back(compressedEt);
          hcalDigiFineGrain0_.push_back(fg);
          hcalDigiFineGrain1_.push_back(fg2);
          hcalDigiFineGrain2_.push_back(fg3);
          hcalDigiFineGrain3_.push_back(fg4);
          hcalDigiVersion_.push_back(hcalVersion);
          hcalDigiEt_.push_back(et);
          hcalDigiCorrectedEt_.push_back(correctedEt);
        }
        expectedTotalET_ += compressedEt;
      }
    }
  }

  // stage 2 layer 1 digis
  for ( const l1t::CaloTower& caloTower : *caloTowerBXs ) { // iterate through towers
    int hwEtEm = caloTower.hwEtEm();
    int hwEtHad = caloTower.hwEtHad();
    int hwEtRatio = caloTower.hwEtRatio();
    int hwEta = caloTower.hwEta();
    int hwPhi = caloTower.hwPhi();
    int hwPt = caloTower.hwPt();
    int hwQual = caloTower.hwQual();

    if (hwPt>0) {
      stage2Layer1DigiHwEtEm_.push_back(hwEtEm);
      stage2Layer1DigiHwEtHad_.push_back(hwEtHad);
      stage2Layer1DigiHwEtRatio_.push_back(hwEtRatio);
      stage2Layer1DigiHwEta_.push_back(hwEta);
      stage2Layer1DigiHwPhi_.push_back(hwPhi);
      stage2Layer1DigiHwPt_.push_back(hwPt);
      stage2Layer1DigiHwQual_.push_back(hwQual);
    }
  }

  // L1 EGamma
  for ( const auto& l1eg : *l1egs ) {
    l1egPt_.push_back(         l1eg.pt());
    l1egEta_.push_back(        l1eg.eta());
    l1egPhi_.push_back(        l1eg.phi());
    l1egHwEta_.push_back(      l1eg.hwEta());
    l1egHwIso_.push_back(      l1eg.hwIso());
    l1egHwPhi_.push_back(      l1eg.hwPhi());
    l1egHwPt_.push_back(       l1eg.hwPt());
    l1egHwQual_.push_back(     l1eg.hwQual());
    l1egFootprintEt_.push_back(l1eg.footprintEt());
    l1egIsoEt_.push_back(      l1eg.isoEt());
    l1egNTT_.push_back(        l1eg.nTT());
    l1egRawEt_.push_back(      l1eg.rawEt());
    l1egShape_.push_back(      l1eg.shape());
    l1egTowerIEta_.push_back(  l1eg.towerIEta());
    l1egTowerIPhi_.push_back(  l1eg.towerIPhi());
  }

  // L1 Tau
  for ( const auto& l1tau : *l1taus ) {
    l1tauPt_.push_back(       l1tau.pt());
    l1tauEta_.push_back(      l1tau.eta());
    l1tauPhi_.push_back(      l1tau.phi());
    l1tauHwEta_.push_back(    l1tau.hwEta());
    l1tauHwIso_.push_back(    l1tau.hwIso());
    l1tauHwPhi_.push_back(    l1tau.hwPhi());
    l1tauHwPt_.push_back(     l1tau.hwPt());
    l1tauHwQual_.push_back(   l1tau.hwQual());
    l1tauHasEM_.push_back(    l1tau.hasEM());
    l1tauIsMerged_.push_back( l1tau.isMerged());
    l1tauIsoEt_.push_back(    l1tau.isoEt());
    l1tauNTT_.push_back(      l1tau.nTT());
    l1tauRawEt_.push_back(    l1tau.rawEt());
    l1tauTowerIEta_.push_back(l1tau.towerIEta());
    l1tauTowerIPhi_.push_back(l1tau.towerIPhi());
  }

  // L1 Jet
  for ( const auto& l1jet : *l1jets ) {
    l1jetPt_.push_back(       l1jet.pt());
    l1jetEta_.push_back(      l1jet.eta());
    l1jetPhi_.push_back(      l1jet.phi());
    l1jetHwEta_.push_back(    l1jet.hwEta());
    l1jetHwIso_.push_back(    l1jet.hwIso());
    l1jetHwPhi_.push_back(    l1jet.hwPhi());
    l1jetHwPt_.push_back(     l1jet.hwPt());
    l1jetHwQual_.push_back(   l1jet.hwQual());
    l1jetPUEt_.push_back(     l1jet.puEt());
    l1jetRawEt_.push_back(    l1jet.rawEt());
    l1jetSeedEt_.push_back(   l1jet.seedEt());
    l1jetTowerIEta_.push_back(l1jet.towerIEta());
    l1jetTowerIPhi_.push_back(l1jet.towerIPhi());
  }

  // L1 Jet
  for ( const auto& l1sum : *l1sums ) {
    if ( l1sum.getType() == l1t::EtSum::kTotalEt ) {
      l1TotalEt_.push_back(     l1sum.et());
      l1TotalEtHwPt_.push_back( l1sum.hwPt());
    }
    if ( l1sum.getType() == l1t::EtSum::kTotalHt ) {
      l1TotalHt_.push_back(     l1sum.et());
      l1TotalHtHwPt_.push_back( l1sum.hwPt());
    }
    if ( l1sum.getType() == l1t::EtSum::kMissingEt ) {
      l1MissingEt_.push_back(     l1sum.et());
      l1MissingEtPhi_.push_back(  l1sum.phi());
      l1MissingEtHwPhi_.push_back(l1sum.hwPhi());
      l1MissingEtHwPt_.push_back( l1sum.hwPt());
    }
  }

  // electrons
  for ( size_t i = 0; i < electrons->size(); ++i ) {
    const auto electron = electrons->ptrAt(i);
    if (!(*electronIdMap)[electron]) continue;
    matchObjectToEcal("electron",*electron);
  }

  // electrons from Z
  for (reco::CandidateView::const_iterator it = electronPairs->begin(), ed = electronPairs->end(); it != ed; ++it) {
    if (it->numberOfDaughters()<2) continue;
    const auto electron = it->daughter(1);
    matchObjectToEcal("electronFromZ",*electron);
  }

  // jets
  for ( const auto& jet : *jets ){
    double NHF = jet.neutralHadronEnergyFraction();
    double NEMF = jet.neutralEmEnergyFraction();
    double CHF = jet.chargedHadronEnergyFraction();
    //double MUF = jet.muonEnergyFraction();
    double CEMF = jet.chargedEmEnergyFraction();
    int NumConst = jet.chargedMultiplicity()+jet.neutralMultiplicity();
    int NumNeutralParticles = jet.neutralMultiplicity();
    int CHM = jet.chargedMultiplicity();
    double eta = jet.eta();
    bool looseJetID = false;
    if (abs(eta)<=2.7) {
      looseJetID = (NHF<0.99 && NEMF<0.99 && NumConst>1) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(eta)>2.4) && abs(eta)<=2.7;
    }
    else if (abs(eta)<=3.0) {
      looseJetID = (NHF<0.98 && NEMF>0.01 && NumNeutralParticles>2 && abs(eta)>2.7 && abs(eta)<=3.0 );
    }
    else {
      looseJetID = (NEMF<0.90 && NumNeutralParticles>10 && abs(eta)>3.0 );
    }
    if (!looseJetID) continue;
    matchObjectToHcal("jet",jet);
  }

  // gen partices
  if (genParticles.isValid() && isMC_) {
    // pions
    for ( const auto& pion : *genParticles ){
      matchObjectToHcal("genPion",pion);
    }
  }

  // met
  for (const auto& m: *met) {
    metEt_.push_back(m.et());
    metPhi_.push_back(m.phi());
  }

  // fill tree
  tree_->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void 
L1Analyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
L1Analyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
L1Analyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1Analyzer);
