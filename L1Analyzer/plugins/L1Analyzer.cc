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

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TTree.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/L1TCalorimeter/interface/CaloTower.h"

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

      edm::Service<TFileService> fs_;
      TTree *tree_;

      // branch variables
      std::vector<int> ecalDigiEta_;
      std::vector<int> ecalDigiPhi_;
      std::vector<int> ecalDigiCompressedEt_;
      std::vector<bool> ecalDigiFineGrain_;

      std::vector<int> hcalDigiEta_;
      std::vector<int> hcalDigiPhi_;
      std::vector<int> hcalDigiCompressedEt_;
      std::vector<bool> hcalDigiFineGrain0_;
      std::vector<bool> hcalDigiFineGrain1_;
      std::vector<bool> hcalDigiFineGrain2_;
      std::vector<bool> hcalDigiFineGrain3_;
      std::vector<int> hcalDigiVersion_;

      std::vector<int> stage2Layer1DigiHwEtEm_;
      std::vector<int> stage2Layer1DigiHwEtHad_;
      std::vector<int> stage2Layer1DigiHwEtRatio_;
      std::vector<int> stage2Layer1DigiHwEta_;
      std::vector<int> stage2Layer1DigiHwPhi_;
      std::vector<int> stage2Layer1DigiHwPt_;
      std::vector<int> stage2Layer1DigiHwQual_;

      uint32_t expectedTotalET_;

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
  stage2Layer1DigisToken_(consumes<l1t::CaloTowerBxCollection>(iConfig.getParameter<edm::InputTag>("stage2Layer1Digis")))
{
   // intialize tfileservice
   usesResource("TFileService");

   // create tree_
   tree_ = fs_->make<TTree>("L1TTree", "L1TTree");

   // initialize branches
   
   // ecalTPs
   tree_->Branch("ecalDigi_ieta", &ecalDigiEta_);
   tree_->Branch("ecalDigi_iphi", &ecalDigiPhi_);
   tree_->Branch("ecalDigi_compressedEt", &ecalDigiCompressedEt_);
   tree_->Branch("ecalDigi_fineGrain", &ecalDigiFineGrain_);

   // hcalTPs
   tree_->Branch("hcalDigi_ieta", &hcalDigiEta_);
   tree_->Branch("hcalDigi_iphi", &hcalDigiPhi_);
   tree_->Branch("hcalDigi_compressedEt", &hcalDigiCompressedEt_);
   tree_->Branch("hcalDigi_fineGrain0", &hcalDigiFineGrain0_);
   tree_->Branch("hcalDigi_fineGrain1", &hcalDigiFineGrain1_);
   tree_->Branch("hcalDigi_fineGrain2", &hcalDigiFineGrain2_);
   tree_->Branch("hcalDigi_fineGrain3", &hcalDigiFineGrain3_);
   tree_->Branch("hcalDigi_version", &hcalDigiVersion_);

   // stage2layer1 digis
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

  // clear branches
  ecalDigiEta_.clear();
  ecalDigiPhi_.clear();
  ecalDigiCompressedEt_.clear();
  ecalDigiFineGrain_.clear();

  hcalDigiEta_.clear();
  hcalDigiPhi_.clear();
  hcalDigiCompressedEt_.clear();
  hcalDigiFineGrain0_.clear();
  hcalDigiFineGrain1_.clear();
  hcalDigiFineGrain2_.clear();
  hcalDigiFineGrain3_.clear();
  hcalDigiVersion_.clear();

  stage2Layer1DigiHwEtEm_.clear();
  stage2Layer1DigiHwEtHad_.clear();
  stage2Layer1DigiHwEtRatio_.clear();
  stage2Layer1DigiHwEta_.clear();
  stage2Layer1DigiHwPhi_.clear();
  stage2Layer1DigiHwPt_.clear();
  stage2Layer1DigiHwQual_.clear();

  expectedTotalET_ = 0;

  // set branch values

  // ecalTPs
  for ( const auto& ecalTp : *ecalTPs ) {
    int caloEta = ecalTp.id().ieta();
    int caloPhi = ecalTp.id().iphi();
    int et = ecalTp.compressedEt();
    bool fgVeto = ecalTp.fineGrain();
    if (et>0) {
      ecalDigiEta_.push_back(caloEta);
      ecalDigiPhi_.push_back(caloPhi);
      ecalDigiCompressedEt_.push_back(et);
      ecalDigiFineGrain_.push_back(fgVeto);
    }
    expectedTotalET_ += et;
  }

  // hcalTPs
  for ( const auto& hcalTp : *hcalTPs ) {
    int caloEta = hcalTp.id().ieta();
    uint32_t absCaloEta = abs(caloEta);
    int hcalVersion = hcalTp.id().version();
    // Tower 29 is not used by Layer-1
    if(absCaloEta == 29) {
      continue;
    }
    // Prevent usage of HF TPs with Layer-1 emulator if HCAL TPs are old style
    else if(hcalTp.id().version() == 0 && absCaloEta > 29) {
      continue;
    }
    else if(absCaloEta <= 41) {
      int caloPhi = hcalTp.id().iphi();
      int et = hcalTp.SOI_compressedEt();
      bool fg  = hcalTp.t0().fineGrain(0);
      bool fg2 = hcalTp.t0().fineGrain(1);
      bool fg3 = hcalTp.t0().fineGrain(2);
      bool fg4 = hcalTp.t0().fineGrain(3);
      if(caloPhi <= 72) {
        uint32_t featureBits = 0;
        if(fg)  featureBits |= 0b01;
        // fg2 should only be set for HF
        if(absCaloEta > 29 && fg2) featureBits |= 0b10;
        if (et>0) {
          hcalDigiEta_.push_back(caloEta);
          hcalDigiPhi_.push_back(caloPhi);
          hcalDigiCompressedEt_.push_back(et);
          hcalDigiFineGrain0_.push_back(fg);
          hcalDigiFineGrain1_.push_back(fg2);
          hcalDigiFineGrain2_.push_back(fg3);
          hcalDigiFineGrain3_.push_back(fg4);
          hcalDigiVersion_.push_back(hcalVersion);
        }
        expectedTotalET_ += et;
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
