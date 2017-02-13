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

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

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

      // ----------member data ---------------------------

      bool verbose_;
      edm::EDGetTokenT<EcalTrigPrimDigiCollection> ecalDigisToken_; 
      edm::EDGetTokenT<HcalTrigPrimDigiCollection> hcalDigisToken_;
      edm::EDGetTokenT<l1t::CaloTowerBxCollection> stage2Layer1DigisToken_;
      edm::EDGetTokenT<reco::CandidateView> electronsToken_;
      edm::EDGetTokenT<reco::CandidateView> photonsToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > electronIdMapToken_;
      edm::EDGetTokenT<reco::CandidateView> electronPairsToken_;
      edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
      bool storeEcal_;
      bool storeHcal_;
      bool storeStage2Layer1_;
      double LSB_;

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
  electronsToken_(consumes<reco::CandidateView>(iConfig.getParameter<edm::InputTag>("electrons"))),
  photonsToken_(consumes<reco::CandidateView>(iConfig.getParameter<edm::InputTag>("photons"))),
  electronIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronIdMap"))),
  electronPairsToken_(consumes<reco::CandidateView>(iConfig.getParameter<edm::InputTag>("electronPairs"))),
  genParticlesToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"))),
  storeEcal_(iConfig.exists("storeEcal") ? iConfig.getParameter<bool>("storeEcal") : true),
  storeHcal_(iConfig.exists("storeHcal") ? iConfig.getParameter<bool>("storeHcal") : true),
  storeStage2Layer1_(iConfig.exists("storeStage2Layer1") ? iConfig.getParameter<bool>("storeStage2Layer1") : true),
  LSB_(iConfig.exists("LSB") ? iConfig.getParameter<double>("LSB") : 0.5)
{
  // intialize tfileservice
  usesResource("TFileService");

  // create tree_
  tree_ = fs_->make<TTree>("L1TTree", "L1TTree");

  // initialize branches
  
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


  // electrons
  buildMatchToEcal("electron");
  buildMatchToEcal("electronFromZ");

  // pions
  buildMatchToHcal("genPion");
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

  for (auto var: vars) {
    std::string bname = name+"_"+var;
    branches_[bname] = std::vector<float>();
    tree_->Branch(bname.c_str(), &branches_[bname]);
  }
}

void
L1Analyzer::buildMatchToHcal(std::string name) {
  buildMatchToEcal(name);

  std::vector<std::string> vars;
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
L1Analyzer::matchObjectToEcal(std::string name, const reco::Candidate& cand) {
  float pt = cand.pt();
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
    double tp_phi = uctGeometry_.getUCTTowerEta(tp_iphi);
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

  // if match found and electron is good
  if (pt>0 && match>0) {
    // calculate sums
    double ecalSum3x3 = 0.;
    double ecalSum5x5 = 0.;
    double ecalSum7x7 = 0.;
    double ecalCorrectedSum3x3 = 0.;
    double ecalCorrectedSum5x5 = 0.;
    double ecalCorrectedSum7x7 = 0.;
    for (int e=-3; e<4; ++e) {
      int ie = N_IETA_+closestIEta;
      if (ie>=N_IETA_) ie-=1;                 // -41 .. -1 1 .. 41 -> 0 .. 81
      ie+=e;
      if (ie<0 || ie>=2*N_IETA_) continue;    // max tower eta 41
      for (int p=-3; p<4; ++p) {
        int ip = closestIPhi-1;               // i .. 72 -> 0 .. 71
        ip+=p;
        if (ip<0) ip+=N_IPHI_;                // wrap -1 -> 71, etc
        if (ip>=N_IPHI_) ip-=N_IPHI_;         // wrap 72 -> 0, etc
        if (e>-2 && e<2 && p>-2 && p<2) {
          ecalSum3x3 += ecalTPs_et_[ie][ip];
          ecalCorrectedSum3x3 += ecalTPs_correctedEt_[ie][ip];
        }
        if (e>-3 && e<3 && p>-3 && p<3) {
          ecalSum5x5 += ecalTPs_et_[ie][ip];
          ecalCorrectedSum5x5 += ecalTPs_correctedEt_[ie][ip];
        }
        ecalSum7x7 += ecalTPs_et_[ie][ip];
        ecalCorrectedSum7x7 += ecalTPs_correctedEt_[ie][ip];
      }
    }

    branches_[name+"_pt"].push_back(pt);
    branches_[name+"_eta"].push_back(eta);
    branches_[name+"_phi"].push_back(phi);
    branches_[name+"_pdgId"].push_back(pdgId);
    branches_[name+"_matchedDR"].push_back(closestDR);
    branches_[name+"_matchedDEt"].push_back(closestDEt);
    branches_[name+"_matchedIEta"].push_back(closestIEta);
    branches_[name+"_matchedIPhi"].push_back(closestIPhi);
    branches_[name+"_matchedEt"].push_back(closestEt);
    branches_[name+"_matchedCorrectedEt"].push_back(closestCorrectedEt);
    branches_[name+"_matchedPtBin"].push_back(closestPtBin);
    branches_[name+"_matchedEcalSum3x3"].push_back(ecalSum3x3);
    branches_[name+"_matchedEcalSum5x5"].push_back(ecalSum5x5);
    branches_[name+"_matchedEcalSum7x7"].push_back(ecalSum7x7);
    branches_[name+"_matchedEcalCorrectedSum3x3"].push_back(ecalCorrectedSum3x3);
    branches_[name+"_matchedEcalCorrectedSum5x5"].push_back(ecalCorrectedSum5x5);
    branches_[name+"_matchedEcalCorrectedSum7x7"].push_back(ecalCorrectedSum7x7);
  }
}

void
L1Analyzer::matchObjectToHcal(std::string name, const reco::Candidate& cand) {
  float pt = cand.pt();
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
    double tp_phi = uctGeometry_.getUCTTowerEta(tp_iphi);
    double tp_et = decoder_->hcaletValue(hcalTP.id(), hcalTP.t0());
    double correctedEt = tp_et*getHcalSF(tp_et,tp_ieta);
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

  // if match found and electron is good
  if (pt>0 && match>0) {
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
        int ip = closestIPhi-1;               // i .. 72 -> 0 .. 71
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
          hcalCorrectedSum3x3 += hcalTPs_correctedEt_[ie][ip];
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

    branches_[name+"_pt"].push_back(pt);
    branches_[name+"_eta"].push_back(eta);
    branches_[name+"_phi"].push_back(phi);
    branches_[name+"_pdgId"].push_back(pdgId);
    branches_[name+"_matchedDR"].push_back(closestDR);
    branches_[name+"_matchedDEt"].push_back(closestDEt);
    branches_[name+"_matchedIEta"].push_back(closestIEta);
    branches_[name+"_matchedIPhi"].push_back(closestIPhi);
    branches_[name+"_matchedEt"].push_back(closestEt);
    branches_[name+"_matchedCorrectedEt"].push_back(closestCorrectedEt);
    branches_[name+"_matchedPtBin"].push_back(closestPtBin);
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
}

// ------------ method called for each event  ------------
void
L1Analyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // grab collections
  iEvent.getByToken(ecalDigisToken_, ecalTPs_);

  iEvent.getByToken(hcalDigisToken_, hcalTPs_);

  edm::Handle<l1t::CaloTowerBxCollection> caloTowerBXs;
  iEvent.getByToken(stage2Layer1DigisToken_, caloTowerBXs);

  edm::Handle<reco::CandidateView> electrons;
  iEvent.getByToken(electronsToken_, electrons);

  edm::Handle<reco::CandidateView> photons;
  iEvent.getByToken(photonsToken_, photons);

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

  // gen partices
  if (genParticles.isValid()) {
    // pions
    for ( const auto& pion : *genParticles ){
      matchObjectToHcal("genPion",pion);
    }
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
