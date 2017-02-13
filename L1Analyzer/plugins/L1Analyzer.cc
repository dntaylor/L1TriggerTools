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

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/ValueMap.h"

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

      // ----------member data ---------------------------

      bool verbose_;
      edm::EDGetTokenT<EcalTrigPrimDigiCollection> ecalDigisToken_; 
      edm::EDGetTokenT<HcalTrigPrimDigiCollection> hcalDigisToken_;
      edm::EDGetTokenT<l1t::CaloTowerBxCollection> stage2Layer1DigisToken_;
      edm::EDGetTokenT<reco::CandidateView> electronsToken_;
      edm::EDGetTokenT<reco::CandidateView> photonsToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > electronIdMapToken_;
      edm::EDGetTokenT<reco::CandidateView> electronPairsToken_;
      bool storeEcal_;
      bool storeHcal_;
      bool storeStage2Layer1_;
      double LSB_;

      edm::Service<TFileService> fs_;
      TTree *tree_;

      UCTGeometry uctGeometry_;
      int N_IETA_ = l1tcalo::MaxCaloEta;
      int N_IPHI_ = l1tcalo::MaxCaloPhi;

      // branch variables
      std::vector<int> ecalDigiEta_;
      std::vector<int> ecalDigiPhi_;
      std::vector<int> ecalDigiCompressedEt_;
      std::vector<bool> ecalDigiFineGrain_;
      std::vector<double> ecalDigiEt_;

      std::vector<int> hcalDigiEta_;
      std::vector<int> hcalDigiPhi_;
      std::vector<int> hcalDigiCompressedEt_;
      std::vector<bool> hcalDigiFineGrain0_;
      std::vector<bool> hcalDigiFineGrain1_;
      std::vector<bool> hcalDigiFineGrain2_;
      std::vector<bool> hcalDigiFineGrain3_;
      std::vector<int> hcalDigiVersion_;
      std::vector<double> hcalDigiEt_;

      std::vector<int> stage2Layer1DigiHwEtEm_;
      std::vector<int> stage2Layer1DigiHwEtHad_;
      std::vector<int> stage2Layer1DigiHwEtRatio_;
      std::vector<int> stage2Layer1DigiHwEta_;
      std::vector<int> stage2Layer1DigiHwPhi_;
      std::vector<int> stage2Layer1DigiHwPt_;
      std::vector<int> stage2Layer1DigiHwQual_;

      uint32_t expectedTotalET_;

      std::vector<float> electronPt_;
      std::vector<float> electronEta_;
      std::vector<float> electronPhi_;
      std::vector<float> electronMatchedDR_;
      std::vector<float> electronMatchedDEt_;
      std::vector<float> electronMatchedIEta_;
      std::vector<float> electronMatchedIPhi_;
      std::vector<float> electronMatchedEt_;
      std::vector<float> electronMatchedSum3x3_;
      std::vector<float> electronMatchedSum5x5_;
      std::vector<float> electronMatchedSum7x7_;

      std::vector<float> photonPt_;
      std::vector<float> photonEta_;
      std::vector<float> photonPhi_;

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
  tree_->Branch("electron_pt", &electronPt_);
  tree_->Branch("electron_eta", &electronEta_);
  tree_->Branch("electron_phi", &electronPhi_);
  tree_->Branch("electron_matchedDR", &electronMatchedDR_);
  tree_->Branch("electron_matchedDEt", &electronMatchedDEt_);
  tree_->Branch("electron_matchedIEta", &electronMatchedIEta_);
  tree_->Branch("electron_matchedIPhi", &electronMatchedIPhi_);
  tree_->Branch("electron_matchedEt", &electronMatchedEt_);
  tree_->Branch("electron_matchedSum3x3", &electronMatchedSum3x3_);
  tree_->Branch("electron_matchedSum5x5", &electronMatchedSum5x5_);
  tree_->Branch("electron_matchedSum7x7", &electronMatchedSum7x7_);

  // photons
  tree_->Branch("photon_pt", &photonPt_);
  tree_->Branch("photon_eta", &photonEta_);
  tree_->Branch("photon_phi", &photonPhi_);
}


L1Analyzer::~L1Analyzer()
{
 
}


//
// member functions
//

// ------------ method called for each event  ------------
void
L1Analyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // grab collections
  edm::Handle<EcalTrigPrimDigiCollection> ecalTPs;
  iEvent.getByToken(ecalDigisToken_, ecalTPs);

  edm::Handle<HcalTrigPrimDigiCollection> hcalTPs;
  iEvent.getByToken(hcalDigisToken_, hcalTPs);

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

  edm::ESHandle<CaloTPGTranscoder> decoder;
  iSetup.get<CaloTPGRecord>().get(decoder);

  // clear branches
  ecalDigiEta_.clear();
  ecalDigiPhi_.clear();
  ecalDigiCompressedEt_.clear();
  ecalDigiFineGrain_.clear();
  ecalDigiEt_.clear();

  hcalDigiEta_.clear();
  hcalDigiPhi_.clear();
  hcalDigiCompressedEt_.clear();
  hcalDigiFineGrain0_.clear();
  hcalDigiFineGrain1_.clear();
  hcalDigiFineGrain2_.clear();
  hcalDigiFineGrain3_.clear();
  hcalDigiVersion_.clear();
  hcalDigiEt_.clear();

  stage2Layer1DigiHwEtEm_.clear();
  stage2Layer1DigiHwEtHad_.clear();
  stage2Layer1DigiHwEtRatio_.clear();
  stage2Layer1DigiHwEta_.clear();
  stage2Layer1DigiHwPhi_.clear();
  stage2Layer1DigiHwPt_.clear();
  stage2Layer1DigiHwQual_.clear();

  expectedTotalET_ = 0;

  electronPt_.clear();
  electronEta_.clear();
  electronPhi_.clear();
  electronMatchedDR_.clear();
  electronMatchedDEt_.clear();
  electronMatchedIEta_.clear();
  electronMatchedIPhi_.clear();
  electronMatchedEt_.clear();
  electronMatchedSum3x3_.clear();
  electronMatchedSum5x5_.clear();
  electronMatchedSum7x7_.clear();

  photonPt_.clear();
  photonEta_.clear();
  photonPhi_.clear();

  // set branch values

  std::vector<std::vector<int> > ecalTPs_compressedEt(N_IETA_*2, std::vector<int>(N_IPHI_));
  std::vector<std::vector<double> > ecalTPs_et(N_IETA_*2, std::vector<double>(N_IPHI_));

  std::vector<std::vector<int> > hcalTPs_compressedEt(N_IETA_*2, std::vector<int>(N_IPHI_));
  std::vector<std::vector<double> > hcalTPs_et(N_IETA_*2, std::vector<double>(N_IPHI_));

  // ecalTPs
  //std::cout << "Number of ecal tps: " << ecalTPs->size() << std::endl;
  for ( const auto& ecalTp : *ecalTPs ) {
    int caloEta = ecalTp.id().ieta();
    int caloPhi = ecalTp.id().iphi();
    int compressedEt = ecalTp.compressedEt();
    double et = compressedEt*LSB_;
    bool fgVeto = ecalTp.fineGrain();
    int ie = caloEta>0 ? N_IETA_+caloEta-1 : N_IETA_+caloEta;
    int ip = caloPhi-1;
    ecalTPs_compressedEt[ie][ip] = compressedEt;
    ecalTPs_et[ie][ip] = et;
    if (compressedEt>0) {
      ecalDigiEta_.push_back(caloEta);
      ecalDigiPhi_.push_back(caloPhi);
      ecalDigiCompressedEt_.push_back(compressedEt);
      ecalDigiFineGrain_.push_back(fgVeto);
      ecalDigiEt_.push_back(et);
    }
    expectedTotalET_ += compressedEt;
  }

  // hcalTPs
  //std::cout << "Number of hcal tps: " << hcalTPs->size() << std::endl;
  for ( const auto& hcalTp : *hcalTPs ) {
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
      double et = decoder->hcaletValue(hcalTp.id(), hcalTp.t0());
      bool fg  = hcalTp.t0().fineGrain(0);
      bool fg2 = hcalTp.t0().fineGrain(1);
      bool fg3 = hcalTp.t0().fineGrain(2);
      bool fg4 = hcalTp.t0().fineGrain(3);
      int ie = caloEta>0 ? N_IETA_+caloEta-1 : N_IETA_+caloEta;
      int ip = caloPhi-1;
      hcalTPs_compressedEt[ie][ip] = compressedEt;
      hcalTPs_et[ie][ip] = et;
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
  //for (reco::CandidateView::const_iterator it = electronPairs->begin(), ed = electronPairs->end(); it != ed; ++it) {
  //  const auto electron = it->daughter(1);
    if (!(*electronIdMap)[electron]) continue;

    float pt = electron->pt();
    float eta = electron->eta();
    float phi = electron->phi();

    // find closest max ecaltp
    double closestDR = 1e6;
    double closestDEt = 1e6;
    double closestEt = 0;
    int closestIEta = -999;
    int closestIPhi = -999;
    int match = -1;
    for ( const auto& ecalTP : *ecalTPs ) {
      int tp_ieta = ecalTP.id().ieta();
      int tp_iphi = ecalTP.id().iphi();
      int tp_compressedEt = ecalTP.compressedEt();
      double tp_eta = uctGeometry_.getUCTTowerEta(tp_ieta);
      double tp_phi = uctGeometry_.getUCTTowerEta(tp_iphi);
      double tp_et = tp_compressedEt*LSB_;
      double dr = deltaR(tp_eta,tp_phi,eta,phi);
      double det = fabs(pt-tp_et);
      if (dr<0.2 && dr<closestDR && pt<127 && det<closestDEt) {
        closestDR = dr;
        closestDEt = det;
        closestIEta = tp_ieta;
        closestIPhi = tp_iphi;
        closestEt = tp_et;
        match = 1;
      }
    }

    // if match found and electron is good
    if (pt>0 && match>0) {
      // calculate sums
      double sum3x3 = 0.;
      double sum5x5 = 0.;
      double sum7x7 = 0.;
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
          if (e>-2 && e<2 && p>-2 && p<2) sum3x3 += ecalTPs_et[ie][ip];
          if (e>-3 && e<3 && p>-3 && p<3) sum5x5 += ecalTPs_et[ie][ip];
          sum7x7 += ecalTPs_et[ie][ip];
        }
      }

      // fill branches
      electronPt_.push_back(pt);
      electronEta_.push_back(eta);
      electronPhi_.push_back(phi);
      electronMatchedDR_.push_back(closestDR);
      electronMatchedDEt_.push_back(closestDEt);
      electronMatchedIEta_.push_back(closestIEta);
      electronMatchedIPhi_.push_back(closestIPhi);
      electronMatchedEt_.push_back(closestEt);
      electronMatchedSum3x3_.push_back(sum3x3);
      electronMatchedSum5x5_.push_back(sum5x5);
      electronMatchedSum7x7_.push_back(sum7x7);
    }
  }

  // photons
  for ( const auto& photon : *photons ) {
    float pt = photon.pt();
    float eta = photon.eta();
    float phi = photon.phi();

    if (pt>0) {
      photonPt_.push_back(pt);
      photonEta_.push_back(eta);
      photonPhi_.push_back(phi);
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
