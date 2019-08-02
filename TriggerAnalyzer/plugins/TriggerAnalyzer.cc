// Package:    HLTAnalysis/TriggerAnalyzer
// Class:      TriggerAnalyzer
// 
/**\class TriggerAnalyzer TriggerAnalyzer.cc HLTAnalysis/TriggerAnalyzer/plugins/TriggerAnalyzer.cc

 Author: Vasileios (Vasilis) Belis

 Description & Notes: 
 //Analyzer made for pTrel bPurity analysis, mainly for MC events. Includes some code remnants, "#include's" and variables of an old Analyzer (e.g.HLTAnalyze function)
 //made by G. Karathanasis. These remnants are not needed for the function of the current Analyzer. They are kept for potential future functionalities of the Analyzer.
          The original Analyzer code of G. Karathanasis  was critical for my work and the creation of this Analyzer.

 Implementation:
     [Notes on implementation]
*/

// system include files
#include <memory>
#include <iostream>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "HLTrigger/HLTcore/interface/defaultModuleLabel.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/L1Trigger/interface/EtSumHelper.h"
#include "DataFormats/L1Trigger/interface/EtSum.h"
#include "DataFormats/Common/interface/ValueMap.h"


#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateIsolation.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include "HLTrigger/HLTcore/interface/HLTFilter.h"
#include "HLTrigger/Egamma/plugins/HLTGenericFilter.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "CondFormats/DataRecord/interface/L1TUtmTriggerMenuRcd.h"
#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/ParametrizedEngine/src/OAEParametrizedMagneticField.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "L1Trigger/L1TNtuples/interface/MuonID.h"
#include <vector>
#include "TTree.h"
#include <string>
#include <iostream>
#include "TMath.h"
#include <cmath>
#include "DataFormats/Common/interface/Ref.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "TLorentzVector.h"
#include "TVector2.h"
#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/CombinedKinematicConstraint.h"
#include <RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h>
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/KinematicFit/interface/MultiTrackPointingKinematicConstraint.h"
#include "CondFormats/DataRecord/interface/L1TUtmTriggerMenuRcd.h"
#include "CondFormats/L1TObjects/interface/L1TUtmAlgorithm.h"
#include "CondFormats/L1TObjects/interface/L1TUtmTriggerMenu.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "DataFormats/L1TGlobal/interface/GlobalExtBlk.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

//VB:
#include "SimDataFormats/JetMatching/interface/JetFlavourInfo.h"

//#include "NtupleContent.h"

/*namespace edm {
  class ConfigurationDescriptions;
  }*/



using namespace std;
//using namespace edm;


//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.
template<typename T1>
class TriggerAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {

  typedef std::vector<T1> T1Collection;
  typedef edm::Ref<T1Collection> T1Ref;
  typedef edm::AssociationMap<edm::OneToValue<std::vector<T1>, float > > T1IsolationMap;

public:
  explicit TriggerAnalyzer(const edm::ParameterSet&);
  ~TriggerAnalyzer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  bool isAncestor(const reco::Candidate* ancestor, const reco::Candidate* particle); 
  
private:
  //----------member class methods: ------------------------
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
 
  std::pair<const reco::GenJet*, const reco::Candidate*> genJetMuAnalyze(const edm::Event& , const edm::EventSetup& );

  std::pair<std::vector<float>,std::vector<std::vector<float>>> L1Analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  std::vector<std::pair<string,int>> L1SeedAnalyze(const edm::Event& iEvent,TString * algoBitToName, std::vector<string> Seed);
  std::pair<std::vector<float>,std::vector<std::vector<std::vector<float>>>> HLTAnalyze(const edm::Event& iEvent, const edm::EventSetup& iSetup,std::vector<string> HLTPath);
  std::vector<std::vector<float>> track_DCA(std::vector<reco::TransientTrack> ttks);
  std::vector<GlobalVector>refit_tracks(TransientVertex myVertex,std::vector<reco::TransientTrack> tracks);

  float Dphi(const float& phi1, const float& phi2);
  float DR(const float& eta1, const float& phi1, const float& eta2, const float& phi2);

  // ----------member data types: ---------------------------
  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> vtxToken_;
  edm::EDGetToken electronsToken_;
  edm::EDGetToken muonsToken_;
  edm::EDGetToken jetsToken_;
  edm::EDGetToken metToken_;
  edm::EDGetToken photonToken_;
  edm::EDGetToken PFCands_;
  edm::EDGetToken LostTracks_;
//  edm::EDGetToken Tracks_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleIdMapVetoToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleIdMapSoftToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleIdMapMediumToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleIdMapTightToken_;
  edm::EDGetTokenT<edm::ValueMap<int> > elIdMapValueToken_;
  edm::EDGetTokenT<GlobalAlgBlkBxCollection> l1resultToken_;
  edm::EDGetToken l1MuonsToken_;
  edm::EDGetToken l1JetsToken_;
  edm::EDGetToken l1MetToken_;
  vector<string> Seed_;
  edm::EDGetTokenT<edm::TriggerResults> trgresultsToken_; 
  edm::EDGetTokenT<vector<pat::TriggerObjectStandAlone>> trigobjectsToken_;
  vector<string> HLTPath_;

  edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
  edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedGenToken_;
  edm::EDGetTokenT<edm::View<reco::GenJet>> GenJetToken_;
 
//  edm::ParameterSet const& conf;
 

  edm::Service<TFileService> fs;
  TTree * t1;
//  NtupleContent nt;
 int JetFlavourTag=0;
 int totb1=0,totb2=0,totb3=0;
  int trigger1=0,trigger2=0,trigger3=0,trigger4=0,trigger5=0,trigger6=0,trigger7=0,trigger8=0;
  int l1_seed1=0,l1_seed2=0,l1_seed3=0,l1_seed4=0,l1_seed5=0,l1_seed6=0;
  int event=0,nmuons=0,njets=0,total_triggers=0,good_vertex=1,nel=0,njpsi_mumu=0,njpsi_ee=0;
  float beam_x=0,beam_y=0,beam_z=0,beam_ex=0,beam_ey=0,beam_ez=0;
  double ptmet=0,phimet=-999,l1met=0,l1met_phi=-999,l1ht=0,l1hf_met=0,l1hf_met_phi=-999;
  int trg_counter=0,fire_counter=0,l3_counter=0,ntracks=0,nkaons=0;
  std::vector<float> muon_pt,muon_eta,muon_phi,muon_qual,muon_charge,muon_global_flag,muon_tracker_flag,muon_standalone_flag,muon_dxy,muon_dz,muon_edxy,muon_edz,muon_d0,muon_ed0,muon_vx,muon_vz,muon_vy,muon_iso,l1muon_eta,l1muon_iso,l1muon_pt,l1muon_phi,l1muon_qual,l1jet_pt,l1jet_phi,l1jet_eta,l1jet_iso,l1jet_qual,muon_trkpt,muon_trketa,muon_trkphi;
  std::vector<std::vector<float>> tr1_obj_pt_eta_phi,tr2_obj_pt_eta_phi,tr3_obj_pt_eta_phi,tr4_obj_pt_eta_phi,tr5_obj_pt_eta_phi,tr6_obj_pt_eta_phi,tr7_obj_pt_eta_phi,tr8_obj_pt_eta_phi;
  std::vector<float> el_pt,el_eta,el_phi,el_charge,el_vx,el_vy,el_vz,el_dxy,el_dz,el_mva_out,el_mva_iso,el_iso,el_trkpt,el_trketa,el_trkphi,el_edxy,el_edz;
  float pvertex_x,pvertex_y,pvertex_z,pvertex_ex,pvertex_ey,pvertex_ez;
  std::vector<float> vertex_x,vertex_y,vertex_z,vertex_ex,vertex_ey,vertex_ez,vertex_chi,vertex_ndof;
  std::vector<bool>el_veto,el_soft,el_medium,el_tight;
  std::vector<float>el_mva_map_value;
  std::vector<float> jet_pt,jet_eta,jet_phi,jet_cEmEF,jet_cHEF,jet_cHMult,jet_cMuEF,jet_cMult,jet_MuEF,jet_nEmEF,jet_nHEF,jet_nMult,jet_pEF,jet_eEF;
  std::vector<bool> muon_medium,muon_loose,muon_tight,muon_soft;
  std::vector<float> charged_pf_pt,charged_pf_eta,charged_pf_phi;
  std::vector<float> track_pt,track_eta,track_phi,track_norm_chi2,track_charge,track_dxy,track_dz,track_validhits,track_losthits,track_fromPV,track_highPurity;
 
  std::vector<std::vector<int>>muon_vtx_IsValid,el_vtx_IsValid;
  int run_number=0,ls=0;
  int npv=0;

;



  //Gen event data:
  std::vector<float> genJet_mass,bgenJet_genQuark_pt_ratio,cgenJet_genQuark_pt_ratio,lightgenJet_genQuark_pt_ratio,genJet_pt,genJet_eta,genJet_y,genJet_phi,genJet_flavour;
  std::vector<float> genJet_genMu_pt_ratio, bJet_genMu_pt_ratio, cJet_genMu_pt_ratio,lightJet_genMu_pt_ratio;
  std::vector<float> genMuB_pt,genMuB_eta,genMuB_phi,genMuC_pt,genMuC_eta,genMuC_phi,genMuUDS_pt,genMuUDS_eta,genMuUDS_phi,dR_bJet_matching,dR_cJet_matching,dR_lightJet_matching;
  std::vector<float> genMu_pt,genMu_eta,genMu_phi,genMu_quarkMO_pdgId,genMu_quarkMO_pt,genMu_quarkMO_eta,genMu_quarkMO_phi;
  std::vector<bool> map_quark_jet_muon;//true when quark is parent to both Jet and Muon.
  std::vector<float> bFragmentation_function_pT,cFragmentation_function_pT,lightFragmentation_function_pT;
  std::vector<float> primHadron_pdgId,primHadron_pt,primHadron_y,primHadron_phi,primHadron_m,ancQuark_pdgId,ancQuark_pt,ancQuark_y,ancQuark_phi,ancQuark_m;
  //pat data (25.06.19):
  std::vector<float> bjet_pt,bjet_eta,bjet_phi,cjet_pt,cjet_eta,cjet_phi,qjet_pt,qjet_eta,qjet_phi;
  std::vector<float> pTb_rel,pTc_rel,pTq_rel,pTb_rel_dir,pTb_rel_indir;
  std::vector<float> mub_dir_pt,mub_dir_eta,mub_dir_phi,mub_seq_pt,mub_seq_eta,mub_seq_phi;

  ///Run Parameters options:
  bool data=true; bool saveTracks=true;
  bool saveHLT=true; bool saveL1=true; bool saveOnlyHLTFires=false;
 //internal
 ParticleMass part_mass = 0.1056583; float part_sigma = 0.0000001;
 ParticleMass kaon_mass = 0.493677; float kaon_sigma = 0.000016;
 float chi = 0.; float ndf = 0.;
 std::pair<std::vector<float>,std::vector<std::vector<float>>> l1objects;
 std::vector<std::pair<string,int>> l1seeds;
 std::vector<std::shared_ptr<reco::TransientTrack>> KTrack;
 std::vector<unsigned int> KTrack_index;
 unsigned int nmupairs=0;
 std::vector<std::shared_ptr<reco::TransientTrack>> muTrack1,muTrack2;//,eTrack1,eTrack2;
 std::vector<std::pair<unsigned int,unsigned int>> used_muTrack_index,used_eTrack_index;
 std::vector<reco::CandidatePtr> footprint;
 std::vector<pat::PackedCandidate> tracks;


 TString * algoBitToName = new TString[512];
 int count=0;
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
template<typename T1>
TriggerAnalyzer<T1>::TriggerAnalyzer(const edm::ParameterSet& iConfig): 

   beamSpotToken_(consumes<reco::BeamSpot>(iConfig.getParameter <edm::InputTag>("beamSpot"))),
   vtxToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices"))),
   electronsToken_(consumes<std::vector<pat::Electron>>(iConfig.getParameter<edm::InputTag>  ("electrons"))),
   muonsToken_(consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
   jetsToken_(consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>  ("jets"))),
   metToken_(consumes<std::vector<pat::MET>>(iConfig.getParameter<edm::InputTag>("met"))),
   photonToken_(consumes<std::vector<pat::Photon>>(iConfig.getParameter<edm::InputTag>("photons"))),
   PFCands_(consumes<std::vector<pat::PackedCandidate> >(iConfig.getParameter<edm::InputTag>("PFCands"))),
   LostTracks_(consumes<std::vector<pat::PackedCandidate> >(iConfig.getParameter<edm::InputTag>("losttracks"))),
  // Tracks_(consumes<std::vector<reco::Track> >(iConfig.getParameter<edm::InputTag>("tracks"))),
   eleIdMapVetoToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleIdMapVeto"))),
   eleIdMapSoftToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleIdMapSoft"))),
   eleIdMapMediumToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleIdMapMedium"))),
   eleIdMapTightToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleIdMapTight"))),
   elIdMapValueToken_(consumes<edm::ValueMap<int> >(iConfig.getParameter<edm::InputTag>("eleIdMapValue"))),
   l1resultToken_(consumes<GlobalAlgBlkBxCollection>(iConfig.getParameter<edm::InputTag>("l1seed"))),
   l1MuonsToken_(consumes<l1t::MuonBxCollection>(iConfig.getParameter<edm::InputTag>("l1muons"))),
   l1JetsToken_(consumes<l1t::JetBxCollection>(iConfig.getParameter<edm::InputTag>("l1jets"))),
   l1MetToken_(consumes<BXVector<l1t::EtSum> >(iConfig.getParameter<edm::InputTag>("l1met"))),
   Seed_(iConfig.getParameter<vector<string> >("Seed")),
   trgresultsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag> ("triggerresults"))),
   trigobjectsToken_(consumes<vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag> ("triggerobjects"))),
   HLTPath_(iConfig.getParameter<vector<string> >("HLTPath")),
   prunedGenToken_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned"))),
   packedGenToken_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packed"))),
   GenJetToken_(consumes<edm::View<reco::GenJet>>(iConfig.getParameter<edm::InputTag>  ("genjets"))) 
  
  {  
   edm::ParameterSet runParameters=iConfig.getParameter<edm::ParameterSet>("RunParameters");

   data=runParameters.getParameter<bool>("Data");
   saveTracks=runParameters.getParameter<bool>("SaveTracks");
   saveHLT=runParameters.getParameter<bool>("SaveHLT");
   saveL1=runParameters.getParameter<bool>("SaveL1");
   saveOnlyHLTFires=runParameters.getParameter<bool>("SaveResultsOnlyIfAPathFired"); 
  }

template<typename T1>
TriggerAnalyzer<T1>::~TriggerAnalyzer()
{
  // cout<<"total "<<trg_counter<<" fires "<<fire_counter<<" l3 "<<l3_counter<<endl;
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
 // cout<<"total trg2(Mu17) "<<totb1<<" trg3(Mu20) "<<totb2<<" trg5(Mu27) "<<totb3<<endl;  
}


//
// member functions
//
  template<typename T1>
  bool TriggerAnalyzer<T1>::isAncestor(const reco::Candidate* ancestor, const reco::Candidate* particle){
  //particle is already the ancestor
      if(ancestor == particle ) return true;
      //otherwise loop on mothers, if any and return true if the ancestor is found
     for(size_t i=0;i< particle->numberOfMothers();i++){
         if(isAncestor(ancestor,particle->mother(i))) return true;
                                                       }
      //if we did not return true yet, then particle and ancestor are not relatives
       return false;
                                                                                                        }
  
// ------------ method called for each event  ------------

    template<typename T1>
    std::pair<const reco::GenJet*, const reco::Candidate*> TriggerAnalyzer<T1>::genJetMuAnalyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

         using namespace std;
         using namespace edm;
         using namespace reco;
         using namespace trigger;
	 using namespace pat;

        JetFlavourTag=0;
        genJet_pt.clear(); genJet_eta.clear();genJet_y.clear();genJet_phi.clear();genJet_flavour.clear();genMuB_pt.clear();genMuB_eta.clear(); genMuB_phi.clear(); 
        genMuC_pt.clear();genMuC_eta.clear(); genMuC_phi.clear();genMuUDS_pt.clear();genMuUDS_eta.clear(); genMuUDS_phi.clear();
        genJet_mass.clear();bgenJet_genQuark_pt_ratio.clear();cgenJet_genQuark_pt_ratio.clear();lightgenJet_genQuark_pt_ratio.clear();
        dR_bJet_matching.clear(); dR_cJet_matching.clear(); dR_lightJet_matching.clear(); 
        genJet_genMu_pt_ratio.clear(), bJet_genMu_pt_ratio.clear(), cJet_genMu_pt_ratio.clear(),lightJet_genMu_pt_ratio.clear();
        genMu_pt.clear(),genMu_eta.clear(),genMu_phi.clear(),genMu_quarkMO_pdgId.clear(),genMu_quarkMO_pt.clear(),genMu_quarkMO_eta.clear(),genMu_quarkMO_phi.clear();
        bFragmentation_function_pT.clear(),cFragmentation_function_pT.clear(),lightFragmentation_function_pT.clear();
        map_quark_jet_muon.clear();

        primHadron_pdgId.clear();
        primHadron_pt.clear();
        primHadron_y.clear();
        primHadron_phi.clear();
        primHadron_m.clear();
        ancQuark_pdgId.clear();
        ancQuark_pt.clear();
        ancQuark_y.clear();
        ancQuark_phi.clear();
        ancQuark_m.clear();

        edm::Handle<edm::View<reco::GenParticle>> genPruned;
        iEvent.getByToken(prunedGenToken_,genPruned);
        edm::Handle<edm::View<pat::PackedGenParticle>> genPacked;
        iEvent.getByToken(packedGenToken_,genPacked);
        edm::Handle<edm::View<reco::GenJet>> genJet;
        iEvent.getByToken(GenJetToken_,genJet);

   
        cout<<"==================START OF EVENT ("<<event<<") GENJET/GENMU COLLECTION ================"<<endl;
/***********************************************************************         
        cout<<"!!!!!PrunedGenParticles Collection:\n"<<endl;
        const reco::Candidate *bquark=nullptr,*primHadron=nullptr;
         for(typename edm::View<reco::GenParticle>::const_iterator prunedPart = genPruned->begin(); prunedPart != genPruned->end();++prunedPart){
            if(abs(prunedPart->pdgId())>=500 && abs(prunedPart->pdgId())<600) primHadron=&*prunedPart;
            if(abs(prunedPart->pdgId())==5){
            bquark=&*prunedPart;
            const reco::Candidate* mo = prunedPart->mother();
            printf("prunedPart: status/pdgId/pt/eta/MOpdgId= %i/%i/%f/%f/%i\n", prunedPart->status(),prunedPart->pdgId(),prunedPart->pt(),prunedPart->eta(),mo->pdgId());              
                for(unsigned int imom=0;imom<mo->numberOfMothers();++imom){
                 const reco::Candidate* recursive_mo = mo->mother(imom); 
                 printf("recursive mother of PRUNED muon: pdgId/status/pdgId/pt= %i/%i/%f\n",recursive_mo->pdgId(),recursive_mo->status(),recursive_mo->pt());
                                                                          }
                                            }
                                                                                                                                                  }
************************************************************************/

//         cout<<"@@@@@PackedGenParticles Collection:\n"<<endl;
         std::pair<const reco::GenJet*, const reco::Candidate*> genJet_genMu_pair;
         genJet_genMu_pair = std::make_pair(nullptr,nullptr);
         float mu_pt=-1.;
         const reco::Candidate *max_mu=nullptr;  
         for(typename edm::View<pat::PackedGenParticle>::const_iterator packedPart = genPacked->begin(); packedPart != genPacked->end();++packedPart){
            if(abs(packedPart->pdgId())==13 && packedPart->pt()>5. && abs(packedPart->eta())<2.1){
//            const reco::Candidate* mo = packedPart->mother(0);
         //   printf("%li packedPart: status/pdgId/pt/eta/MOpdgId= %i/%i/%f/%f/%i\n", packedPart-genPacked->begin(),packedPart->status(),packedPart->pdgId(),packedPart->pt(),packedPart->eta(),mo->pdgId());              
//                for(unsigned int imom=0;imom<mo->numberOfMothers();++imom){
//                 const reco::Candidate* recursive_mo = mo->mother(imom); 
//                 printf("recursive mother of Packed muon: pdgId/status/pt= %i/%i/%f\n",recursive_mo->pdgId(),recursive_mo->status(),recursive_mo->pt());
//                                                            }
              if(mu_pt<packedPart->pt()){
                max_mu = &(*packedPart);
                mu_pt=packedPart->pt();
                                        }
                                                                                                     }
                                                                                                                                                     }
      if(max_mu != nullptr){
       printf("muon pT_max= %f\n",max_mu->pt());
       float min_deltaR=9999;
         for(typename edm::View<reco::GenJet>::const_iterator genjet=genJet->begin(); genjet !=genJet->end();++genjet){
          if(fabs(genjet->eta())>2.1 || genjet->pt()<10.) continue;//first lets have no genJet cut. in MINI-AOD only pT>8 Jets are saved/merged.
          if(min_deltaR<= deltaR(max_mu->eta(),max_mu->phi(),genjet->eta(),genjet->phi())) continue;
             min_deltaR= deltaR(max_mu->eta(),max_mu->phi(),genjet->eta(),genjet->phi());  
             genJet_genMu_pair = std::make_pair(&*genjet,max_mu);
                                                                                                                       }
///           
        cout<<"=================END OF EVENT ("<<event<<") GENJET/GENMU COLLECTION ==============="<<endl;
      
      //Checking on how to access mu's which are paired with collected GenJets...
//            std::string ancestor_name = "mother_";
//           const reco::Candidate* mo= genJet_genMu_pair.second->mother(0);
//           float genjet_eta=genJet_genMu_pair.first->eta(),genjet_phi=genJet_genMu_pair.first->phi();
//           printf("---CHECK: GenMu: /status/pt/DR(mu-jet)= /%i/%f/%f \n",genJet_genMu_pair.second->status(),genJet_genMu_pair.second->pt(),deltaR(genJet_genMu_pair.second->eta(),genJet_genMu_pair.second->phi(),genjet_eta,genjet_phi));
//           for(unsigned int imom=0;imom<mo->numberOfMothers();++imom){
//                  const reco::Candidate* recursive_mo = mo->mother(imom);
//                  cout<<"---^"+ancestor_name+to_string(imom+1)+"-> "<<recursive_mo->pdgId()<<"/"<<recursive_mo->status()<<"/"<<recursive_mo->pt()<<"/"<<deltaR(recursive_mo->eta(),recursive_mo->phi(),genjet_eta,genjet_phi)<<"\n";
//                                                                      }
      if(genJet_genMu_pair.first == nullptr) cout<<"No genJet passed the cuts"<<endl;
      else{
         cout<<"=================START OF EVENT ("<<event<<") GEN: QUARK-MU-JET MATCHING ==============="<<endl;
          vector<const reco::GenParticle *> genQuark_collection;
          vector<unsigned int> genQuark_collection_pdgId;
          bool btag=false,ctag=false,udstag=false;
      
           for(typename edm::View<reco::GenParticle>::const_iterator genq=genPruned->begin(); genq != genPruned->end() ; ++genq){
      
              if(fabs(genq->pdgId())>5 || genq->status()!=71 || genq->isLastCopy()==false || deltaR(genq->eta(),genq->phi(),genJet_genMu_pair.first->eta(),genJet_genMu_pair.first->phi())>0.4) continue;
                     genQuark_collection.push_back(&*genq);
                     genQuark_collection_pdgId.push_back(fabs(genq->pdgId()));
        //             cout<<"GenQuark ('"<<genQuark_collection.size()<<"') merged into GenJet (collected index'"<<q<<"')-cone:|pdgId|/pT/[pT(Quark)/pT(Jet)]/= "<<fabs(genq->pdgId())<<"/"<<genq->pt()<<"/["<<genq->pt()/genJet_genMu_pair.first->pt()<<"]/" << endl;
                                                                                                                            }
                 if(genQuark_collection.empty())  printf("Mu-paired GenJet  in event ('%i') did not have a genQuark in cone\n",event);   
                 else{
                   unsigned int max_pdgId = *std::max_element(genQuark_collection_pdgId.begin(),genQuark_collection_pdgId.end());
                   unsigned int max_index = std::distance(genQuark_collection_pdgId.begin(),std::max_element(genQuark_collection_pdgId.begin(),genQuark_collection_pdgId.end()));     
                   //27/1/19: iterate on genQuark_collection, and compare the pdgId's of quarks not placed on max_index position, then take the one with the minimum 1-pT(quark)/pT(jet)
                   for(unsigned int i_quark=0;i_quark<genQuark_collection.size();++i_quark){
                     if(i_quark==max_index) continue;
                     if(genQuark_collection_pdgId.at(i_quark)== max_pdgId){
        //               cout<<"!-!-! Found two GenQuark (c.i.:('"<<max_index<<"' & '"<<i_quark<<"')"<<" in GenJet (c.i.:'"<<q<<"') with the max_pdgId. Comparing their DR & 1-pT(quark)/pT(jet)."<<endl;
                       float pt_ratio_max_index= 1.-genQuark_collection.at(max_index)->pt()/genJet_genMu_pair.first->pt(),pt_ratio_i_quark= 1.-genQuark_collection.at(i_quark)->pt()/genJet_genMu_pair.first->pt();
      /*
                       float DR_jet_quark_max_index = deltaR(genQuark_collection.at(max_index)->eta(),genQuark_collection.at(max_index)->phi(),genJet_genMu_pair.first->eta(),genJet_genMu_pair.first->phi());
                       float DR_jet_quark_i_quark= deltaR(genQuark_collection.at(i_quark)->eta(),genQuark_collection.at(i_quark)->phi(),genJet_genMu_pair.first->eta(),genJet_genMu_pair.first->phi());
      */
          //             printf("1-pT(quark'%u')/pT(jet)= %f, 1-pT(quark '%u')/pT(jet)= %f, DR(quark'%u'-jet)= %f, DR(quark'%u'-jet)= %f \n",max_index,pt_ratio_max_index,i_quark,pt_ratio_i_quark,max_index,DR_jet_quark_max_index,i_quark,DR_jet_quark_i_quark);
                       if(pt_ratio_max_index<pt_ratio_i_quark) continue;
                        max_index=i_quark;
                                                                           }
                                                                                           }
             
      float genQuark_genJet_pt_ratio = genQuark_collection.at(max_index)->pt()/genJet_genMu_pair.first->pt(); 
                   if(max_pdgId==5){cout<<"b-tag: pT(Quark)/pT(Jet) = "<<genQuark_genJet_pt_ratio<<endl;
                                    JetFlavourTag = genQuark_collection.at(max_index)->pdgId();
                                     btag=true;
//                                    printf("GenJet: pT/eta/phi= %f/%f/%f\n",genJet_genMu_pair.first->pt(),genJet_genMu_pair.first->eta(),genJet_genMu_pair.first->phi());
                                   }
                   else if(max_pdgId==4){cout<<"c-tag: pT(Quark)/pT(Jet)= "<<genQuark_genJet_pt_ratio<<endl;
                                         JetFlavourTag = genQuark_collection.at(max_index)->pdgId();
                                          ctag=true;
//                                         printf("GenJet: pT/eta/phi= %f/%f/%f\n",genJet_genMu_pair.first->pt(),genJet_genMu_pair.first->eta(),genJet_genMu_pair.first->phi());
                                        }
                   else{cout<<"uds-tag: pT(Quark)/pT(Jet)= "<<genQuark_genJet_pt_ratio<<endl;
                        JetFlavourTag = genQuark_collection.at(max_index)->pdgId();
//                        printf("GenJet: pT/eta/phi= %f/%f/%f\n",genJet_genMu_pair.first->pt(),genJet_genMu_pair.first->eta(),genJet_genMu_pair.first->phi());
                        udstag=true;
                       }

	 //(17.07.19) Collection of outgoin quarks from hard process. Before any gluon radiation. PYTHIA status code 23:
            for(typename edm::View<reco::GenParticle>::const_iterator genq=genPruned->begin(); genq != genPruned->end() ; ++genq){
               if(fabs(genq->pdgId())>5 ||  deltaR(genq->eta(),genq->phi(),genJet_genMu_pair.first->eta(),genJet_genMu_pair.first->phi())>0.4) continue;
	       printf("####genQuark: pdgId/status/pT/p/E/m=%i/%i/%f/%f/%f/%f\n",genq->pdgId(),genq->status(),genq->pt(),genq->p(),genq->energy(),genq->mass());

	    }
	 //(28.07.19) Collection of outgoin quarks from hard process. Before any gluon radiation. PYTHIA status code 23:
            for(typename edm::View<reco::GenParticle>::const_iterator genh=genPruned->begin(); genh != genPruned->end() ; ++genh){
               if(std::to_string(abs(genh->pdgId())).length()<3 ||  deltaR(genh->eta(),genh->phi(),genJet_genMu_pair.first->eta(),genJet_genMu_pair.first->phi())>0.4) continue;
	       printf("####genHadron: pdgId/status/pT/p/E/mass=%i/%i/%f/%f/%f/%f\n",genh->pdgId(),genh->status(),genh->pt(),genh->p(),genh->energy(),genh->mass());

	    }

        const reco::Candidate* recursive_mo = genQuark_collection.at(max_index)->mother();
        const reco::Candidate* hardQuark = genQuark_collection.at(max_index);
	//while(recursive_mo != nullptr && recursive_mo->pdgId() != 2212 && recursive_mo->status() != 23){
	while(recursive_mo != nullptr && recursive_mo->pdgId() != 2212 && recursive_mo->pdgId() != 21){
		printf("loop, recursive_quark: pdgId/status/pT=%i/%i/%f\n",recursive_mo->pdgId(),recursive_mo->status(),recursive_mo->pt());
		recursive_mo = recursive_mo->mother();
		hardQuark = hardQuark->mother();
	}
        printf("last recursive_quark: pdgId/status/pT=%i/%i/%f\n",recursive_mo->pdgId(),recursive_mo->status(),recursive_mo->pt());
        printf("hardQuark: pdgId/status/pT = %i/%i/%f\n",hardQuark->pdgId(),hardQuark->status(),hardQuark->pt());
        recursive_mo = genJet_genMu_pair.second->mother();
        const reco::Candidate* first_fragmentation_meson = genJet_genMu_pair.second;
        while(recursive_mo != nullptr && recursive_mo->pdgId() != 2212 && std::to_string(abs(recursive_mo->pdgId())).length()>=3){//length of pdgId>=3 Hadrons, not working for negative pdgId....
//              cout<<"recursive_mo->pdgId()= "<<recursive_mo->pdgId()<<endl; 
//              cout<<"recursive_mo->status()= "<<recursive_mo->status()<<endl; 
              recursive_mo = recursive_mo->mother();                               
              first_fragmentation_meson = first_fragmentation_meson->mother();
                                                                                   }        
printf("first_fragmentation_meson: pdgId/status/pT = %i/%i/%f\n",first_fragmentation_meson->pdgId(),first_fragmentation_meson->status(),first_fragmentation_meson->pt());
           bool successful_match = isAncestor(genQuark_collection.at(max_index),genJet_genMu_pair.second) && isAncestor(genQuark_collection.at(max_index),first_fragmentation_meson);
           map_quark_jet_muon.push_back(successful_match);

         if(successful_match){
//	     if((int)(abs(first_fragmentation_meson->pdgId())/100) ==0){//Let's check if there is a flaw in first_fragmentation_meson selection
//             cout<<"$$$$$$$$$$$$first_fragmentation_meson->pdgId()= "<<first_fragmentation_meson->pdgId()<<endl;
//             cout<<"$$$$$$$$$$$$pT(meson)/pT(quark)= "<<first_fragmentation_meson->pt()/genQuark_collection.at(max_index)->pt()<<endl;
//	                                                               }
             float dR = deltaR(genJet_genMu_pair.second->eta(),genJet_genMu_pair.second->phi(),genJet_genMu_pair.first->eta(),genJet_genMu_pair.first->phi());           
              genMu_pt.push_back(genJet_genMu_pair.second->pt());
              genMu_eta.push_back(genJet_genMu_pair.second->eta());
              genMu_phi.push_back(genJet_genMu_pair.second->phi());

              primHadron_pdgId.push_back(first_fragmentation_meson->pdgId());
	      primHadron_pt.push_back(first_fragmentation_meson->pt());
	      primHadron_y.push_back(first_fragmentation_meson->y());
	      primHadron_phi.push_back(first_fragmentation_meson->phi());
	      primHadron_m.push_back(first_fragmentation_meson->mass());
	      ancQuark_pdgId.push_back(genQuark_collection.at(max_index)->pdgId());
	      ancQuark_pt.push_back(genQuark_collection.at(max_index)->pt());
	      ancQuark_y.push_back(genQuark_collection.at(max_index)->y());
	      ancQuark_phi.push_back(genQuark_collection.at(max_index)->phi());
	      ancQuark_m.push_back(genQuark_collection.at(max_index)->mass());

              genMu_quarkMO_pdgId.push_back(genQuark_collection.at(max_index)->pdgId());
              genMu_quarkMO_pt.push_back(genQuark_collection.at(max_index)->pt());
              genMu_quarkMO_eta.push_back(genQuark_collection.at(max_index)->eta());
              genMu_quarkMO_phi.push_back(genQuark_collection.at(max_index)->phi());

              genJet_genMu_pt_ratio.push_back(genJet_genMu_pair.second->pt()/genJet_genMu_pair.first->pt());
              genJet_pt.push_back(genJet_genMu_pair.first->pt());
              genJet_eta.push_back(genJet_genMu_pair.first->eta());
              genJet_y.push_back(genJet_genMu_pair.first->y());
              genJet_phi.push_back(genJet_genMu_pair.first->phi());
              genJet_mass.push_back(genJet_genMu_pair.first->mass());
              genJet_flavour.push_back(genQuark_collection.at(max_index)->pdgId());

             if(btag){
                      bFragmentation_function_pT.push_back(first_fragmentation_meson->pt()/genQuark_collection.at(max_index)->pt());
//		      if(bFragmentation_function_pT.at(0)<first_fragmentation_meson->pt()/genQuark_collection.at(max_index)->pt())printf("bFragmentation function= %f, pT(meson)/pT(quark)= %f\n",bFragmentation_function_pT.at(0),first_fragmentation_meson->pt()/genQuark_collection.at(max_index)->pt());
                      genMuB_pt.push_back(genJet_genMu_pair.second->pt());   
                      genMuB_eta.push_back(genJet_genMu_pair.second->eta());   
                      genMuB_phi.push_back(genJet_genMu_pair.second->phi());   
                      dR_bJet_matching.push_back(dR);    
                      bgenJet_genQuark_pt_ratio.push_back(genQuark_genJet_pt_ratio);
                      bJet_genMu_pt_ratio.push_back(genJet_genMu_pair.second->pt()/genJet_genMu_pair.first->pt());
                     }
             if(ctag){
                      cFragmentation_function_pT.push_back(first_fragmentation_meson->pt()/genQuark_collection.at(max_index)->pt());
//		      if(cFragmentation_function_pT.at(0)<first_fragmentation_meson->pt()/genQuark_collection.at(max_index)->pt())printf("cFragmentation function= %f, pT(meson)/pT(quark)= %f\n",cFragmentation_function_pT.at(0),first_fragmentation_meson->pt()/genQuark_collection.at(max_index)->pt());
                      genMuC_pt.push_back(genJet_genMu_pair.second->pt());   
                      genMuC_eta.push_back(genJet_genMu_pair.second->eta());   
                      genMuC_phi.push_back(genJet_genMu_pair.second->phi());   
                      dR_cJet_matching.push_back(dR);    
                      cgenJet_genQuark_pt_ratio.push_back(genQuark_genJet_pt_ratio);
                      cJet_genMu_pt_ratio.push_back(genJet_genMu_pair.second->pt()/genJet_genMu_pair.first->pt());
                     }
             if(udstag){
                      lightFragmentation_function_pT.push_back(first_fragmentation_meson->pt()/genQuark_collection.at(max_index)->pt());
//		      if(lightFragmentation_function_pT.at(0)<first_fragmentation_meson->pt()/genQuark_collection.at(max_index)->pt())printf("lightFragmentation function= %f, pT(meson)/pT(quark)= %f\n",lightFragmentation_function_pT.at(0),first_fragmentation_meson->pt()/genQuark_collection.at(max_index)->pt());
                      genMuUDS_pt.push_back(genJet_genMu_pair.second->pt());   
                      genMuUDS_eta.push_back(genJet_genMu_pair.second->eta());   
                      genMuUDS_phi.push_back(genJet_genMu_pair.second->phi());   
                      dR_lightJet_matching.push_back(dR);    
                      lightgenJet_genQuark_pt_ratio.push_back(genQuark_genJet_pt_ratio);
                      lightJet_genMu_pt_ratio.push_back(genJet_genMu_pair.second->pt()/genJet_genMu_pair.first->pt());
                       }
                              }
         else cout<<std::boolalpha<<"successful_match= "<<successful_match<<endl; 
                }//if_empty_genQuark_collection       
                }//if_genJet_genMu_pair.first == nullptr
                }//if_nullptr_mu'si
        else{ cout<<"No muons passed the cuts"<<endl; }
        cout<<"=================END OF EVENT ("<<event<<") GEN: QUARK-MU-JET MATCHING ==============="<<endl;
             return genJet_genMu_pair;                                                                       }


template<typename T1>
   std::pair<std::vector<float>,std::vector<std::vector<float>>>  TriggerAnalyzer<T1>::L1Analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
     using namespace std;
     using namespace edm;
     using namespace reco;
     using namespace trigger;
   
    
     edm::Handle<GlobalAlgBlkBxCollection> l1result;
     iEvent.getByToken(l1resultToken_,l1result);
     edm::Handle<l1t::MuonBxCollection> l1Muons;
     iEvent.getByToken(l1MuonsToken_,l1Muons);
     edm::Handle<l1t::JetBxCollection> l1Jets;
     iEvent.getByToken(l1JetsToken_,l1Jets);
     edm::Handle<BXVector<l1t::EtSum> > l1Met;
     iEvent.getByToken(l1MetToken_,l1Met);
    
   
     std::vector<float> l1met;
     std::vector< std::vector<float>> l1event;
     std::vector<float> l1muon_pt,l1muon_eta,l1muon_phi,l1muon_iso,l1muon_qual,l1jet_pt,l1jet_eta,l1jet_phi,l1jet_iso,l1jet_qual;
     for(typename std::vector< l1t::Muon >::const_iterator mu=l1Muons->begin(0); mu !=l1Muons->end(0); mu++){
       l1muon_pt.push_back(mu->et());
       l1muon_eta.push_back(mu->eta());
       l1muon_phi.push_back(mu->phi());
       l1muon_iso.push_back(mu->hwIso());
       l1muon_qual.push_back(mu->hwQual());
   }
     l1event.push_back(l1muon_pt); l1event.push_back(l1muon_eta); l1event.push_back(l1muon_phi);  l1event.push_back(l1muon_iso);  l1event.push_back(l1muon_qual); 
      for(typename std::vector< l1t::Jet >::const_iterator jet=l1Jets->begin(0); jet!=l1Jets->end(0); jet++){
         l1jet_pt.push_back(jet->et());
         l1jet_eta.push_back(jet->eta());
         l1jet_phi.push_back(jet->phi());
         l1jet_iso.push_back(jet->hwIso());
         l1jet_qual.push_back(jet->hwQual()); 
   }
     l1event.push_back(l1jet_pt); l1event.push_back(l1jet_eta); l1event.push_back(l1jet_phi);  l1event.push_back(l1jet_iso);  l1event.push_back(l1jet_qual); 
     for(typename std::vector< l1t::EtSum >::const_iterator et=l1Met->begin(0); et !=l1Met->end(0); et++){
     if (et->getType()==l1t::EtSum::kMissingEt) {l1met.push_back(et->et()); l1met.push_back(et->phi()); }
     if (et->getType()==l1t::EtSum::kTotalHt) l1met.push_back(et->et()); 
     if (et->getType()==l1t::EtSum::kMissingEtHF){ l1met.push_back(et->et()); l1met.push_back(et->phi()); }
       }
     
         return std::make_pair(l1met,l1event);
   }

template<typename T1>
    std::vector<std::pair<string,int>>  TriggerAnalyzer<T1>::L1SeedAnalyze(const edm::Event& iEvent, TString * algoBitToName,std::vector< string> Seed){
      using namespace std;
      using namespace edm;
      using namespace reco;
      using namespace trigger;
       edm::Handle<GlobalAlgBlkBxCollection> l1result;
      iEvent.getByToken(l1resultToken_,l1result); 
      std::vector<std::pair<string,int>> res;
    if (l1result.isValid()) {
    GlobalAlgBlk const &result = l1result->at(0, 0);
    for (unsigned int iseed=0; iseed<Seed.size(); ++iseed){
     bool sfire=false;
     for (unsigned int itrg=0; itrg<result.maxPhysicsTriggers; ++itrg){
       if (result.getAlgoDecisionFinal(itrg)!=1) continue;
       std::string l1trigName = static_cast<const char *>(algoBitToName[itrg]);
       if (l1trigName!=Seed[iseed]) continue;
       sfire=true;
       }
     if (sfire) {
       std::pair<string,int> temp(Seed[iseed],1);
       res.push_back(temp);
      }  
     else {
        std::pair<string,int> temp(Seed[iseed],0);
        res.push_back(temp);
      }
     //  std::cout<<"bit "<<itrg<<" algo "<<algoBitToName[itrg]<<" reult "<<result.getAlgoDecisionFinal(itrg)<<endl;   
    }
    }
    return res;
    }



   template<typename T1>
   std::pair<std::vector<float>,std::vector<std::vector<std::vector<float>>>>  TriggerAnalyzer<T1>::HLTAnalyze(const edm::Event& iEvent, const edm::EventSetup& iSetup,std::vector< string> HLTPath ){
      using namespace std;
     using namespace edm;
     using namespace reco;
     using namespace trigger;
    
     edm::Handle<vector<pat::TriggerObjectStandAlone>> triggerObjects;
     iEvent.getByToken(trigobjectsToken_ ,triggerObjects);
     edm::Handle<edm::TriggerResults> trigResults;
     iEvent.getByToken(trgresultsToken_, trigResults);
     std::vector<float> fires;
     std::vector<std::vector<std::vector<float>>> trg_event; 
     const edm::TriggerNames & trigName = iEvent.triggerNames(*trigResults);
     
   for (unsigned int ifilter=0; ifilter<HLTPath.size(); ifilter++){    
       std::vector<std::vector<float>> tot_tr_obj_pt_eta_phi;
        for(pat::TriggerObjectStandAlone itrg :*triggerObjects){
         if(!itrg.id(83)) continue;
          bool save=false;
          itrg.unpackPathNames(trigName);
          std::vector<std::string> const & pathnames = itrg.pathNames(true,true);
          for(std::string name : pathnames){
             if (name.find(HLTPath[ifilter]) != std::string::npos) save=true;
   	   }
           std::vector<float> tr_obj_pt_eta_phi;
   	if(!save) continue;
           tr_obj_pt_eta_phi.push_back(itrg.pt());
           tr_obj_pt_eta_phi.push_back(itrg.eta());
           tr_obj_pt_eta_phi.push_back(itrg.phi());
           tr_obj_pt_eta_phi.push_back(itrg.charge());
           tot_tr_obj_pt_eta_phi.push_back( tr_obj_pt_eta_phi);
   	}
      if(tot_tr_obj_pt_eta_phi.size()>0){
        trg_event.push_back(tot_tr_obj_pt_eta_phi);}            
      else{
        std::vector<float> def; def.push_back(-1); 
        std::vector<std::vector<float>> def2; def2.push_back(def);
        trg_event.push_back(def2);
      }
     }
     //paths
     float fire0=0,fire1=0,fire2=0,fire3=0,fire4=0,fire5=0,fire6=0,fire7=0;
       
       if( !trigResults.failedToGet() ) {
       int N_Triggers = trigResults->size();
       for( int i_Trig = 0; i_Trig < N_Triggers; ++i_Trig ) {
         TString TrigPath =trigName.triggerName(i_Trig);
         if (!trigResults->accept(i_Trig)) continue;
         if (TrigPath.Contains(HLTPath[0])) fire0=1;
         if (TrigPath.Contains(HLTPath[1])) fire1=1;       
         if (TrigPath.Contains(HLTPath[2])) fire2=1; 
         if (TrigPath.Contains(HLTPath[3])) fire3=1;
         if (TrigPath.Contains(HLTPath[4])) fire4=1;
         if (TrigPath.Contains(HLTPath[5])) fire5=1;
         if (TrigPath.Contains(HLTPath[6])) fire6=1;
         if (TrigPath.Contains(HLTPath[7])) fire7=1;
          } 
         }
         							   	
       fires.push_back(fire0); fires.push_back(fire1); fires.push_back(fire2);
       fires.push_back(fire3); fires.push_back(fire4); fires.push_back(fire5);
       fires.push_back(fire6); fires.push_back(fire7);
       //  return std::make_pair(fires,tot_tr_obj_pt_eta_phi);
        return std::make_pair(fires,trg_event);
   }

  template<typename T1> 
  std::vector<std::vector<float>> TriggerAnalyzer<T1>::track_DCA(std::vector<reco::TransientTrack> ttks) {
   std::vector<std::vector<float>> dca;
   std::vector<float> def;
   def.push_back(-9999999);
   if(ttks.size()<2) {dca.push_back(def); return dca; }
   for(unsigned int tk1=0; tk1<ttks.size(); tk1++){
     TrajectoryStateClosestToPoint mu1TS = ttks[tk1].impactPointTSCP();
     std::vector<float> temp1;
     for(unsigned int tk2=tk1+1; tk2<ttks.size(); tk2++){
        TrajectoryStateClosestToPoint mu2TS = ttks[tk2].impactPointTSCP();      
        if (mu1TS.isValid() && mu2TS.isValid()) {
          ClosestApproachInRPhi cdca;
  	cdca.calculate(mu1TS.theState(), mu2TS.theState());
          if (cdca.status()) temp1.push_back(cdca.distance());
          else temp1.push_back(-999999999);
        }
        else temp1.push_back(-99999999);
     }
     dca.push_back(temp1);         
   }  
   return dca;
   }

template<typename T1>
   std::vector<GlobalVector>
   TriggerAnalyzer<T1>::refit_tracks(TransientVertex myVertex,std::vector<reco::TransientTrack> tracks){
       std::auto_ptr<TrajectoryStateClosestToPoint> traj1;
       std::auto_ptr<TrajectoryStateClosestToPoint> traj2;
       GlobalPoint vtxPos(myVertex.position().x(), myVertex.position().y(), myVertex.position().z());
        GlobalVector gvmu1,gvmu2;
        if(myVertex.hasRefittedTracks()){
           std::vector<reco::TransientTrack> refited;
           refited=myVertex.refittedTracks();
           reco::TransientTrack* Track1 = &refited[0];
           reco::TransientTrack* Track2= &refited[1];
           traj1.reset(new TrajectoryStateClosestToPoint(Track1->trajectoryStateClosestToPoint(vtxPos)));
           traj2.reset(new TrajectoryStateClosestToPoint(Track2->trajectoryStateClosestToPoint(vtxPos)));    
           gvmu1=traj1->momentum();  gvmu2=traj2->momentum(); 
        }         
        else {
           traj1.reset(new TrajectoryStateClosestToPoint(tracks[0].trajectoryStateClosestToPoint(vtxPos)));
           traj2.reset(new TrajectoryStateClosestToPoint(tracks[1].trajectoryStateClosestToPoint(vtxPos)));
           gvmu1=traj1->momentum();  gvmu2=traj2->momentum();
        }                       
      std::vector<GlobalVector> gvmu; gvmu.push_back(gvmu1); gvmu.push_back(gvmu2);
      return gvmu;
   }

template<typename T1>
void
TriggerAnalyzer<T1>::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  using namespace std;
  using namespace edm;
  using namespace reco;
  using namespace trigger;

  ++event;

  //Get a few collections to apply basic electron ID
  //Get data
  edm::Handle<reco::BeamSpot> theBeamSpot;
  iEvent.getByToken(beamSpotToken_,theBeamSpot); 
   edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  //continue if there are no vertices
  if (vertices->size()==0) return;
  edm::Handle<std::vector<pat::Electron>> electrons;
  iEvent.getByToken(electronsToken_, electrons);
  edm::Handle<std::vector<pat::Jet>> jets;
  iEvent.getByToken(jetsToken_, jets);
  edm::Handle<std::vector<pat::Muon>> muons;
  iEvent.getByToken(muonsToken_,muons);
  edm::Handle<std::vector<pat::MET>> met;
  iEvent.getByToken(metToken_,met);
  edm::Handle<std::vector<pat::Photon>> photons;
  iEvent.getByToken(photonToken_,photons);
  edm::Handle<vector<pat::PackedCandidate>> tracks1;
  iEvent.getByToken(PFCands_, tracks1);
  edm::Handle<vector<pat::PackedCandidate>> tracks2;
  iEvent.getByToken(LostTracks_, tracks2);
//  edm::Handle<vector<reco::Track>> tracks3;
//  iEvent.getByToken(Tracks_, tracks3);
  edm::Handle<edm::ValueMap<bool> > ele_veto_id;
  iEvent.getByToken(eleIdMapVetoToken_ ,ele_veto_id);
  edm::Handle<edm::ValueMap<bool> > ele_soft_id;
  iEvent.getByToken(eleIdMapSoftToken_ ,ele_soft_id);
  edm::Handle<edm::ValueMap<bool> > ele_medium_id;
  iEvent.getByToken(eleIdMapMediumToken_ ,ele_medium_id);
  edm::Handle<edm::ValueMap<bool> > ele_tight_id;
  iEvent.getByToken(eleIdMapTightToken_ ,ele_tight_id);
  edm::Handle<edm::ValueMap<int> > ele_mva_id_value;
  iEvent.getByToken( elIdMapValueToken_ ,ele_mva_id_value);
  edm::ESHandle<MagneticField> bFieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);
  // GeneratorBTree gen(prunedGenToken_,packedGenToken_,iEvent);
edm::Handle<GlobalAlgBlkBxCollection> l1result;
  iEvent.getByToken(l1resultToken_,l1result);
if (count==0){
edm::ESHandle<L1TUtmTriggerMenu> menu;
  iSetup.get<L1TUtmTriggerMenuRcd>().get(menu);

if (l1result.isValid()) {
  for (auto const & keyval: menu->getAlgorithmMap()) {
   std::string const & trigName  = keyval.second.getName();
   unsigned int index = keyval.second.getIndex();  
   int itrig = index;
   algoBitToName[itrig] = TString( trigName );
 }
  count++;
  }
 }

  

 vertex_x.clear();  vertex_y.clear();  vertex_z.clear();  vertex_ex.clear();  vertex_ey.clear();  vertex_ez.clear();  vertex_chi.clear();  vertex_ndof.clear(); npv=0;
  pvertex_x=-9999; pvertex_y=-9999; pvertex_z=-9999; pvertex_ex=-9999; pvertex_ey=-9999; pvertex_ez=-9999,beam_x=-9999,beam_y=-9999,beam_z=-99999,beam_ex=-9999,beam_ey=-9999,beam_ez=-99999;
  footprint.clear(); 
  trigger1=0; trigger2=0; trigger3=0; trigger4=0; trigger5=0; trigger6=0;
   trigger7=0; trigger8=0;
 // event=iEvent.id().event();
  nmuons=0; njets=0; nel=0; ntracks=0;
  //muons
  muon_pt.clear(); muon_eta.clear(); muon_phi.clear(); muon_qual.clear();
  muon_charge.clear(); muon_dz.clear(); muon_dxy.clear(); muon_global_flag.clear();
  muon_standalone_flag.clear(); muon_tracker_flag.clear();  muon_edz.clear();
  muon_edxy.clear(); muon_d0.clear(); muon_ed0.clear(); muon_vx.clear();
  muon_vy.clear(); muon_vz.clear(); muon_iso.clear(); muon_soft.clear();
  muon_loose.clear(); muon_medium.clear(); muon_tight.clear();
  muon_trkpt.clear(); muon_trketa.clear(); muon_trkphi.clear();
  //e
  el_pt.clear(); el_eta.clear(); el_phi.clear(); el_charge.clear();
  el_vx.clear(); el_vy.clear(); el_vz.clear(); el_dxy.clear(); el_dz.clear();
  el_edxy.clear(); el_edz.clear(); el_mva_out.clear(); el_mva_iso.clear();
  el_iso.clear(); el_veto.clear(); el_soft.clear(); el_medium.clear();
  el_tight.clear();
  //hlt
  tr1_obj_pt_eta_phi.clear(); tr2_obj_pt_eta_phi.clear(); tr3_obj_pt_eta_phi.clear();
  tr4_obj_pt_eta_phi.clear(); tr5_obj_pt_eta_phi.clear(); tr6_obj_pt_eta_phi.clear();
   tr7_obj_pt_eta_phi.clear(); tr8_obj_pt_eta_phi.clear();
  //l1
  l1muon_pt.clear(); l1muon_eta.clear(); l1muon_phi.clear(); l1muon_iso.clear();
  l1muon_qual.clear(); l1jet_pt.clear(); l1jet_eta.clear(); l1jet_phi.clear();
  l1jet_iso.clear(); l1jet_qual.clear();


  ptmet=-1; phimet=-99; jet_pt.clear(); jet_phi.clear(); jet_eta.clear();
  l1met=-1; l1met_phi=-999; l1ht=-1; l1hf_met=-1; l1hf_met_phi=-999;
 
  
  bjet_pt.clear(),bjet_eta.clear(),bjet_phi.clear(),cjet_pt.clear(),cjet_eta.clear(),cjet_phi.clear(),qjet_pt.clear(),qjet_eta.clear(),qjet_phi.clear();
  pTb_rel.clear(),pTc_rel.clear(),pTq_rel.clear(),pTb_rel_dir.clear(),pTb_rel_indir.clear();
  mub_dir_pt.clear(),mub_dir_eta.clear(),mub_dir_phi.clear(),mub_seq_pt.clear(),mub_seq_eta.clear(),mub_seq_phi.clear();

  jet_cEmEF.clear(); jet_cHEF.clear(); jet_cHMult.clear(); jet_cMuEF.clear();
  jet_cMult.clear(); jet_MuEF.clear(); jet_nEmEF.clear(); jet_nHEF.clear();
  jet_nMult.clear(); jet_pEF.clear(); jet_eEF.clear();
 track_pt.clear(); track_eta.clear(); track_phi.clear(); track_norm_chi2.clear(); track_charge.clear();  track_dxy.clear(); track_dz.clear(); track_validhits.clear(); track_losthits.clear(); track_fromPV.clear(); track_highPurity.clear();
 
  el_trkpt.clear(); el_trketa.clear(); el_trkphi.clear();

 el_mva_map_value.clear();
l1_seed1=0; l1_seed2=0; l1_seed3=0; l1_seed4=0; l1_seed5=0; l1_seed6=0;


  run_number=iEvent.id().run();
  ls=iEvent.luminosityBlock();

  beam_x= theBeamSpot->x0();  beam_y= theBeamSpot->y0(); beam_z= theBeamSpot->z0();
  beam_ex= theBeamSpot->x0Error(); beam_ey= theBeamSpot->y0Error(); beam_ez= theBeamSpot->z0Error();

     const  reco::Vertex firstGoodVertex=vertices->front();
     int firstGoodVertexIdx = 0;
     float pvx=-99,pvy=-99,pvz=-99,pevx=-99,pevy=-99,pevz=-99;
     std::vector<reco::TransientTrack> muttks,ettks,jpsimuttks,jpsiettks;
    for (const reco::Vertex &vtx : *vertices) {
         bool isFake = vtx.isFake();
         if ( isFake || !vtx.isValid () ) continue;
         if (firstGoodVertexIdx==0){
              firstGoodVertexIdx=1; 
	      pvx=vtx.x(); pvy=vtx.y(); pvz=vtx.z(); pevx=vtx.xError(); pevy=vtx.yError(); pevz=vtx.zError();    }						                      vertex_x.push_back(vtx.x());  vertex_y.push_back(vtx.y());  vertex_z.push_back(vtx.z()); 
	      vertex_ex.push_back(vtx.xError());  vertex_ey.push_back(vtx.yError());  vertex_ez.push_back(vtx.zError()); 
	      vertex_chi.push_back(vtx.chi2()); vertex_ndof.push_back(vtx.ndof());
    }

      pvertex_x=pvx; pvertex_y=pvy; pvertex_z=pvz; 
      pvertex_ex=pevx; pvertex_ey=pevy; pvertex_ez=pevz;
      reco::TrackBase::Point vertex_point; vertex_point.SetCoordinates(pvx,pvy,pvz);


  if(!data){
    std::pair<const reco::GenJet*, const reco::Candidate*> genJet_genMu_pair = genJetMuAnalyze(iEvent,iSetup);
    if(genJet_genMu_pair.first != nullptr && genJet_genMu_pair.second != nullptr){
    cout<<"=================START OF EVENT ("<<event<<") pat::JET & MUON ANALYSIS ==============="<<endl;
      for (const pat::Jet &jet : *jets){
        if(jet.pt()<10. || abs(jet.eta())>2.1) continue;
        if(jet.neutralHadronEnergyFraction()>=0.90 || jet.neutralEmEnergyFraction()>=0.90 || (jet.neutralMultiplicity()+jet.chargedMultiplicity())<=1 || jet.chargedHadronMultiplicity()<=0 || jet.chargedMultiplicity()<=0) continue;
        const reco::GenJet* genJet_matched = jet.pat::Jet::genJet();
           if(genJet_matched != nullptr){
                 if(genJet_genMu_pair.first == genJet_matched){
                 bool found_muon = false;
		 bool flavourTag_success = false;
                    cout<<"PAT::JET HAS CORRECTLY MATCHED TO RECO::GENJET"<<endl;
//                    printf("matched GenJet: pT/eta/phi= %f/%f/%f\n", genJet_matched->pt(),genJet_matched->eta(),genJet_matched->phi());                    
//                    printf("associated GenMu: pT/eta/phi= %f/%f/%f\n", genJet_genMu_pair.second->pt(),genJet_genMu_pair.second->eta(),genJet_genMu_pair.second->phi());
                    cout<<"jet.pat::Jet::jetFlavourInfo()->getPartonFlavour()= "<< jet.pat::Jet::jetFlavourInfo().getPartonFlavour()<<endl;
                    cout<<"My JetFlavourTag = "<<JetFlavourTag<<endl; 
//                    if(jet.genParton() != nullptr) printf("jet.genParton: pT/pdgId/status/ = %f/%i/%i/\n",jet.genParton()->pt(),jet.genParton()->pdgId(),jet.genParton()->status());
//		    printf("getcHadrons.size(),getbHadrons().size(),getLeptons().size(),getPartons().size() = %lu,%lu,%lu,%lu\n",jet.jetFlavourInfo().getcHadrons().size(),jet.jetFlavourInfo().getbHadrons().size(),jet.jetFlavourInfo().getLeptons().size(),jet.jetFlavourInfo().getPartons().size());
//	            for(unsigned int icHadron=0; icHadron<jet.jetFlavourInfo().getcHadrons().size();++icHadron) printf("Lepton inside Jet: pT/eta/pdgId/ = %f/%f/%i\n",jet.jetFlavourInfo().getcHadrons().at(icHadron).pt(),jet.jetFlavourInfo().getcHadrons().at(icHadron).eta(),jet.jetFlavourInfo().getcHadrons().at(icHadron).pdgId());
//	            for(unsigned int ibHadron=0; ibHadron<jet.jetFlavourInfo().getbHadrons().size();++ibHadron) printf("Lepton inside Jet: pT/eta/pdgId/ = %f/%f/%i\n",jet.jetFlavourInfo().getbHadrons().at(ibHadron).pt(),jet.jetFlavourInfo().getbHadrons().at(ibHadron).eta(),jet.jetFlavourInfo().getbHadrons().at(ibHadron).pdgId());
//	            for(unsigned int iParton=0; iParton<jet.jetFlavourInfo().getPartons().size();++iParton) printf("Lepton inside Jet: pT/eta/pdgId/ = %f/%f/%i\n",jet.jetFlavourInfo().getPartons().at(iParton).pt(),jet.jetFlavourInfo().getPartons().at(iParton).eta(),jet.jetFlavourInfo().getPartons().at(iParton).pdgId());
	            //for(unsigned int iLepton=0; iLepton<jet.jetFlavourInfo().getLeptons().size();++iLepton) printf("Lepton inside Jet: pT/eta/pdgId/ = %f/%f/%i\n",jet.jetFlavourInfo().getLeptons().at(0).at(iLepton).pt(),jet.jetFlavourInfo().getLeptons().at(iLepton).eta(),jet.jetFlavourInfo().getLeptons().at(iLepton).pdgId());
                    cout<<">>>>>>> pat::Mu analysis:"<<endl;
                    for (const pat::Muon &mu : *muons){   
                    if(!mu.pat::Muon::isSoftMuon(firstGoodVertex)) continue;
                    if(mu.pat::Lepton<reco::Muon>::genLepton() == nullptr ) continue; 
                    float DR_genLepton_pairedMu = deltaR(mu.genLepton()->eta(),mu.genLepton()->phi(),genJet_genMu_pair.second->eta(),genJet_genMu_pair.second->phi());
                    float DR_patMu_pairedMu = deltaR(mu.eta(),mu.phi(),genJet_genMu_pair.second->eta(),genJet_genMu_pair.second->phi());

                    //printf("mu.pat::Lepton<reco::Muon>::genLepton(): pT/eta/phi/deltaR(this,paired) = %f/%f/%f/%f\n",mu.genLepton()->pt(),mu.genLepton()->eta(),mu.genLepton()->phi(),DR_genLepton_pairedMu);
                   // printf("pat::Muon: pT/eta/phi/deltaR(this,paired) = %f/%f/%f/%f\n",mu.pt(),mu.eta(),mu.phi(),deltaR(mu.eta(),mu.phi(),genJet_genMu_pair.second->eta(),genJet_genMu_pair.second->phi()));
                    printf("mu.simPdgId()/mu.simFlavour()/mu.simHeaviestMotherFlavour()/mu.simMotherPdgId()/mu.simType() = %i/%i/%i/%i/%i\n",mu.simPdgId(),mu.simFlavour(),mu.simHeaviestMotherFlavour(),mu.simMotherPdgId(),mu.simType());
                       if(DR_genLepton_pairedMu>=0.005 || DR_patMu_pairedMu>=0.005 || abs((mu.genLepton()->pt()-mu.pt())/mu.genLepton()->pt() )> 0.075) continue;
                          cout <<"GOOD PAT::MUON GENMUON MATCH"<<endl;
                          if(abs(JetFlavourTag)==abs(jet.pat::Jet::jetFlavourInfo().getPartonFlavour()) && abs(JetFlavourTag)==mu.simHeaviestMotherFlavour()){
			      flavourTag_success=true;
                              TLorentzVector p4_mu;
                              p4_mu.SetPtEtaPhiM(mu.pt(),mu.eta(),mu.phi(),0.10566);
                              TVector3 p_jet(jet.px(),jet.py(),jet.pz()),p_mu = p4_mu.Vect();
		      	    if(abs(JetFlavourTag)==5){
				 cout<<"bJet tag pTb_rel"<<endl;
		      		 bjet_pt.push_back(jet.pt());bjet_eta.push_back(jet.eta());bjet_phi.push_back(jet.phi());
                                 pTb_rel.push_back((p_mu.Cross(p_jet)).Mag()/p_jet.Mag());		  
                                 if(abs(JetFlavourTag)==mu.simFlavour()){pTb_rel_dir.push_back((p_mu.Cross(p_jet)).Mag()/p_jet.Mag()); cout<<"Direct b->muon"<<endl;}
		      		 else if(mu.simHeaviestMotherFlavour() == 5 && mu.simFlavour()<=4){pTb_rel_indir.push_back((p_mu.Cross(p_jet)).Mag()/p_jet.Mag());cout<<"Indirect b->c->muon"<<endl;}
		      		                     }
		      	    if(abs(JetFlavourTag)==4){
				 cout<<"cJet tag pTc_rel"<<endl;
		      		 cjet_pt.push_back(jet.pt());cjet_eta.push_back(jet.eta());cjet_phi.push_back(jet.phi());
                                 pTc_rel.push_back((p_mu.Cross(p_jet)).Mag()/p_jet.Mag());		  
		      		                }
		      	    if(abs(JetFlavourTag)<4 && abs(JetFlavourTag)>0){
				 cout<<"qJet tag pTq_rel"<<endl;
		      		 qjet_pt.push_back(jet.pt());qjet_eta.push_back(jet.eta());qjet_phi.push_back(jet.phi());
                                 pTq_rel.push_back((p_mu.Cross(p_jet)).Mag()/p_jet.Mag());		  
		      		                }
		      	                                                                        }//if_jetFlavourTagging is consistent between custom and formal
	                  else if(abs(jet.pat::Jet::jetFlavourInfo().getPartonFlavour())==21 && mu.simHeaviestMotherFlavour()<=5 && mu.simHeaviestMotherFlavour()>0){
		      	    //Include also g jets with on flight K/pi->muon decays
			    cout<<"We got one gJet and pTq_rel"<<endl;
			    flavourTag_success=true;
                            TLorentzVector p4_mu;
                            p4_mu.SetPtEtaPhiM(mu.pt(),mu.eta(),mu.phi(),0.10566);
                            TVector3 p_jet(jet.px(),jet.py(),jet.pz()),p_mu = p4_mu.Vect();
		      	    qjet_pt.push_back(jet.pt());qjet_eta.push_back(jet.eta());qjet_phi.push_back(jet.phi());
                            pTq_rel.push_back((p_mu.Cross(p_jet)).Mag()/p_jet.Mag());		  
		        	                       													      }
	                  if(flavourTag_success==false) continue;
			    cout<<"Saving muon data"<<endl;
                          muon_pt.push_back(mu.pt());  muon_phi.push_back(mu.phi());
                          muon_eta.push_back(mu.eta()); muon_charge.push_back(mu.charge());
                          const Track * mutrack= mu.bestTrack();
                          muon_dz.push_back(mutrack->dz(vertex_point));
                          muon_dxy.push_back(mutrack->dxy(vertex_point));
                          muon_global_flag.push_back(mu.isGlobalMuon());
                          muon_standalone_flag.push_back(mu.isStandAloneMuon());
                          muon_tracker_flag.push_back(mu.isTrackerMuon());
                          muon_vx.push_back(mu.vx());  muon_vy.push_back(mu.vy());
                          muon_vz.push_back(mu.vz());  muon_edz.push_back(mu.dzError());
                          muon_edxy.push_back(mu.dxyError());  muon_d0.push_back(mutrack->d0());
                          muon_ed0.push_back(mutrack->d0Error());
                          muon_medium.push_back(mu.isMediumMuon());
                          muon_loose.push_back(mu.isLooseMuon());
                          muon_trkpt.push_back(mutrack->pt()); muon_trketa.push_back(mutrack->eta());
                          muon_trkphi.push_back(mutrack->phi());
                          muon_tight.push_back(mu.isTightMuon(firstGoodVertex));
                          muon_soft.push_back(mu.isSoftMuon(firstGoodVertex));
                         //if (mu.isSoftMuon(firstGoodVertex))snmu++;
                          const MuonPFIsolation&  isol=mu.pfIsolationR04();
                          double mu_iso=(isol.sumChargedHadronPt+max(0.,isol.sumNeutralHadronEt+isol.sumPhotonEt-0.5*isol.sumPUPt))/mu.pt();
                          muon_iso.push_back(mu_iso);
                          muttks.push_back(reco::TransientTrack(*mutrack,&(*bFieldHandle)));   
                          nmuons++;
                          found_muon = true;
                          break;//for_pat::Muon
                                                      }
                    if(!found_muon) continue;
                    cout<<"Saving Jet data"<<endl;
                    jet_pt.push_back(jet.pt()); jet_eta.push_back(jet.eta());
                    jet_phi.push_back(jet.phi());
                    jet_cEmEF.push_back(jet.chargedEmEnergyFraction());
                    jet_cHEF.push_back(jet.chargedHadronEnergyFraction());
                    jet_cHMult.push_back(jet.chargedHadronMultiplicity());
                    jet_cMuEF.push_back(jet.chargedMuEnergyFraction());
                    jet_cMult.push_back(jet.chargedMultiplicity());
                    jet_MuEF.push_back(jet.muonEnergyFraction());
                    jet_eEF.push_back(jet.electronEnergyFraction());
                    jet_nEmEF.push_back(jet.neutralEmEnergyFraction());
                    jet_nHEF.push_back(jet.neutralHadronEnergyFraction());
                    jet_nMult.push_back(jet.neutralMultiplicity());
                    jet_pEF.push_back(jet.photonEnergyFraction());
                    njets++;
                     break;//for_pat::jet
                                                              }//if_genJet_genMu_pair.first == jet.pat::genJet()
                                        }//if_genJet_matched != nullptr
                                       }//for_pat::Jets
                                              }//if_genJet_genMu_pair = nullptr's
    cout<<"=================END OF EVENT ("<<event<<") pat::JET & MUON ANALYSIS ==============="<<endl;
/*****ONLY MUONS*****
           for (const pat::Muon &mu : *muons){
  	   cout<<"iEvent= "<<event<<endl;
	   if(mu.pt()<5 || mu.simHeaviestMotherFlavour() !=5) continue;
           if(mu.pat::Lepton<reco::Muon>::genLepton() == nullptr ) continue; 
           printf("mu.simPdgId()/mu.simFlavour()/mu.simHeaviestMotherFlavour()/mu.simMotherPdgId()/mu.simType() = %i/%i/%i/%i/%i\n",mu.simPdgId(),mu.simFlavour(),mu.simHeaviestMotherFlavour(),mu.simMotherPdgId(),mu.simType());
	   if(mu.simHeaviestMotherFlavour() == mu.simFlavour()){
		   cout<<"b->muon"<<endl;
		   mub_dir_pt.push_back(mu.pt()); mub_dir_eta.push_back(mu.eta());mub_dir_phi.push_back(mu.phi());
		   }  
	   else if(mu.simHeaviestMotherFlavour()>mu.simFlavour() && mu.simFlavour()>0){
		   cout<<"b->...->muon"<<endl;
		   mub_seq_pt.push_back(mu.pt()); mub_seq_eta.push_back(mu.eta());mub_seq_phi.push_back(mu.phi());
		  }
	   }
*****ONLY MUONS*****/

           }//if_!data
  else{
      for (const pat::Jet &jet : *jets){ 
        if(jet.pt()<10. || abs(jet.eta())>2.1) continue;
        jet_pt.push_back(jet.pt()); jet_eta.push_back(jet.eta());
        jet_phi.push_back(jet.phi());
        jet_cEmEF.push_back(jet.chargedEmEnergyFraction());
        jet_cHEF.push_back(jet.chargedHadronEnergyFraction());
        jet_cHMult.push_back(jet.chargedHadronMultiplicity());
        jet_cMuEF.push_back(jet.chargedMuEnergyFraction());
        jet_cMult.push_back(jet.chargedMultiplicity());
        jet_MuEF.push_back(jet.muonEnergyFraction());
        jet_eEF.push_back(jet.electronEnergyFraction());
        jet_nEmEF.push_back(jet.neutralEmEnergyFraction());
        jet_nHEF.push_back(jet.neutralHadronEnergyFraction());
        jet_nMult.push_back(jet.neutralMultiplicity());
        jet_pEF.push_back(jet.photonEnergyFraction());
        njets++;
                                       }

      }

  if(saveL1){
      l1objects=L1Analyze(iEvent, iSetup);
      l1muon_pt=l1objects.second[0];    l1muon_eta=l1objects.second[1];
      l1muon_phi=l1objects.second[2];   l1muon_iso=l1objects.second[3];
      l1muon_qual=l1objects.second[4];  l1jet_pt=l1objects.second[5];
      l1jet_eta=l1objects.second[6];    l1jet_phi=l1objects.second[7];
      l1jet_iso=l1objects.second[8];    l1jet_qual=l1objects.second[9];
      l1met=l1objects.first[0];         l1met_phi=l1objects.first[1]; 
      l1ht=l1objects.first[2];          l1hf_met=l1objects.first[3];
      l1hf_met_phi=l1objects.first[4];
      l1seeds=L1SeedAnalyze(iEvent,algoBitToName,Seed_);
      l1_seed1=l1seeds[0].second; l1_seed2=l1seeds[1].second; 
      l1_seed3=l1seeds[2].second; l1_seed4=l1seeds[3].second; 
      l1_seed5=l1seeds[4].second; l1_seed6=l1seeds[5].second;
   }
  if (saveHLT){  
    std::pair<std::vector<float>,std::vector<std::vector<std::vector<float>>>> trgresult=HLTAnalyze(iEvent,iSetup,HLTPath_);     
     trigger1=trgresult.first[0]; trigger2=trgresult.first[1]; trigger3=trgresult.first[2];  
     trigger4=trgresult.first[3]; trigger5=trgresult.first[4]; trigger6=trgresult.first[5];   
     trigger7=trgresult.first[6]; trigger8=trgresult.first[7];
      if(trigger1==1) tr1_obj_pt_eta_phi=trgresult.second[0];  
      if(trigger2==1) tr2_obj_pt_eta_phi=trgresult.second[1];  
      if(trigger3==1) tr3_obj_pt_eta_phi=trgresult.second[2];   
      if(trigger4==1) tr4_obj_pt_eta_phi=trgresult.second[3]; 
      if(trigger5==1) tr5_obj_pt_eta_phi=trgresult.second[4]; 
      if(trigger6==1) tr6_obj_pt_eta_phi=trgresult.second[5];  
      if(trigger7==1) tr7_obj_pt_eta_phi=trgresult.second[6];
      if(trigger8==1) tr8_obj_pt_eta_phi=trgresult.second[7];
    }
    totb1+=trigger2; totb2+=trigger3; totb3+=trigger5;
    //cout<<"B trg1 "<<trigger1<<" trg2 "<<trigger2<<" trg3 "<<trigger3<<endl;

    t1->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
template<typename T1>
void 
TriggerAnalyzer<T1>::beginJob()
{
  t1=fs->make<TTree>("mytree","mytree");
  //  nt.SetTree(t1);
//  nt.SetNtupleVariables("ALL");
  t1->Branch("event",&event); t1->Branch("run_number",&run_number);
  t1->Branch("ls",&ls);
  
  t1->Branch("mub_dir_pt",&mub_dir_pt);
  t1->Branch("mub_dir_eta",&mub_dir_eta);
  t1->Branch("mub_dir_phi",&mub_dir_phi);
  t1->Branch("mub_seq_pt",&mub_seq_pt);
  t1->Branch("mub_seq_eta",&mub_seq_eta);
  t1->Branch("mub_seq_phi",&mub_seq_phi);
  t1->Branch("bjet_pt",&bjet_pt);
  t1->Branch("bjet_eta",&bjet_eta);
  t1->Branch("bjet_phi",&bjet_phi);
  t1->Branch("cjet_pt",&cjet_pt);
  t1->Branch("cjet_eta",&cjet_eta);
  t1->Branch("cjet_phi",&cjet_phi);
  t1->Branch("qjet_pt",&qjet_pt);
  t1->Branch("qjet_eta",&qjet_eta);
  t1->Branch("qjet_phi",&qjet_phi);
  t1->Branch("pTb_rel",&pTb_rel);
  t1->Branch("pTb_rel_dir",&pTb_rel_dir);
  t1->Branch("pTb_rel_indir",&pTb_rel_indir);
  t1->Branch("pTc_rel",&pTc_rel);
  t1->Branch("pTq_rel",&pTq_rel);
  
  t1->Branch("primHadron_pdgId",&primHadron_pdgId);
  t1->Branch("primHadron_pt",&primHadron_pt);
  t1->Branch("primHadron_y",&primHadron_y);
  t1->Branch("primHadron_phi",&primHadron_phi);
  t1->Branch("primHadron_m",&primHadron_m);
  t1->Branch("ancQuark_pdgId",&ancQuark_pdgId);
  t1->Branch("ancQuark_pt",&ancQuark_pt);
  t1->Branch("ancQuark_y",&ancQuark_y);
  t1->Branch("ancQuark_phi",&ancQuark_phi);
  t1->Branch("ancQuark_m",&ancQuark_m);

  t1->Branch("map_quark_jet_muon",&map_quark_jet_muon);
  t1->Branch("bgenJet_genQuark_pt_ratio",&bgenJet_genQuark_pt_ratio);
  t1->Branch("dR_bJet_matching",&dR_bJet_matching);
  t1->Branch("cgenJet_genQuark_pt_ratio",&cgenJet_genQuark_pt_ratio);
  t1->Branch("dR_cJet_matching",&dR_cJet_matching);
  t1->Branch("lightgenJet_genQuark_pt_ratio",&lightgenJet_genQuark_pt_ratio);
  t1->Branch("dR_lightJet_matching",&dR_lightJet_matching);
  t1->Branch("genJet_genMu_pt_ratio",&genJet_genMu_pt_ratio);
  t1->Branch("bJet_genMu_pt_ratio",&bJet_genMu_pt_ratio);
  t1->Branch("cJet_genMu_pt_ratio",&cJet_genMu_pt_ratio);
  t1->Branch("lightJet_genMu_pt_ratio",&lightJet_genMu_pt_ratio);
  t1->Branch("genJet_pt",&genJet_pt);
  t1->Branch("genJet_eta",&genJet_eta);
  t1->Branch("genJet_y",&genJet_y);
  t1->Branch("genJet_mass",&genJet_mass);
  t1->Branch("genJet_phi",&genJet_phi);
  t1->Branch("genJet_flavour",&genJet_flavour); 
  t1->Branch("bFragmentation_function_pT",&bFragmentation_function_pT);
  t1->Branch("cFragmentation_function_pT",&cFragmentation_function_pT);
  t1->Branch("lightFragmentation_function_pT",&lightFragmentation_function_pT);;
  t1->Branch("genMu_pt",&genMu_pt);
  t1->Branch("genMu_eta",&genMu_eta);
  t1->Branch("genMu_phi",&genMu_phi);
  t1->Branch("genMu_quarkMO_pdgId",&genMu_quarkMO_pdgId);
  t1->Branch("genMu_quarkMO_pt",&genMu_quarkMO_pt);
  t1->Branch("genMu_quarkMO_eta",&genMu_quarkMO_eta);
  t1->Branch("genMu_quarkMO_phi",&genMu_quarkMO_phi);
  t1->Branch("genMuB_pt",&genMuB_pt);
  t1->Branch("genMuB_eta",&genMuB_eta);
  t1->Branch("genMuB_phi",&genMuB_phi);
  t1->Branch("genMuC_pt",&genMuC_pt);
  t1->Branch("genMuC_eta",&genMuC_pt);
  t1->Branch("genMuC_phi",&genMuC_phi);
  t1->Branch("genMuUDS_pt",&genMuUDS_pt);
  t1->Branch("genMuUDS_eta",&genMuUDS_eta);
  t1->Branch("genMuUDS_phi",&genMuUDS_phi);

  t1->Branch("HLT_path1",&trigger1); t1->Branch("HLT_path2",&trigger2);
  t1->Branch("HLT_path3",&trigger3); t1->Branch("HLT_path4",&trigger4);
  t1->Branch("HLT_path5",&trigger5); t1->Branch("HLT_path6",&trigger6);
  t1->Branch("HLT_path7",&trigger5); t1->Branch("HLT_path8",&trigger6);
  t1->Branch("L1_seed1",&l1_seed1); t1->Branch("L1_seed2",&l1_seed2);
  t1->Branch("L1_seed3",&l1_seed3); t1->Branch("L1_seed4",&l1_seed4);
  t1->Branch("L1_seed5",&l1_seed5); t1->Branch("L1_seed6",&l1_seed6);
  t1->Branch("met_pt",&ptmet); t1->Branch("met_phi",&phimet);
  t1->Branch("nmuon",&nmuons); t1->Branch("muon_pt",&muon_pt);
  t1->Branch("muon_eta",&muon_eta); t1->Branch("muon_phi",&muon_phi);
  t1->Branch("muon_charge",&muon_charge); t1->Branch("muon_qual",&muon_qual);
  t1->Branch("muon_dxy",&muon_dxy); t1->Branch("muon_dz",&muon_dz);
  t1->Branch("muon_edxy",&muon_edxy); t1->Branch("muon_edz",&muon_edz);
  t1->Branch("muon_d0",&muon_d0); t1->Branch("muon_ed0",&muon_ed0);
  t1->Branch("muon_vx",&muon_vx); t1->Branch("muon_vy",&muon_vy);
  t1->Branch("muon_vz",&muon_vz); t1->Branch("muon_iso",&muon_iso);
  t1->Branch("muon_global_flag",&muon_global_flag);
  t1->Branch("muon_tracker_flag",&muon_tracker_flag);
  t1->Branch("muon_standalone_flag",&muon_standalone_flag);
  t1->Branch("muon_soft",&muon_soft); t1->Branch("muon_loose",&muon_loose);
  t1->Branch("muon_medium",&muon_medium);  t1->Branch("muon_tight",&muon_tight);
  t1->Branch("jet_pt",&jet_pt); t1->Branch("jet_eta",&jet_eta);
  t1->Branch("jet_phi",&jet_phi); t1->Branch("jet_cEmEF",&jet_cEmEF);
  t1->Branch("jet_cHEF",&jet_cHEF); t1->Branch("jet_cHMult",&jet_cHMult);
  t1->Branch("jet_cMuEF",&jet_cMuEF); t1->Branch("jet_cMult",&jet_cMult);
  t1->Branch("jet_MuEF",&jet_MuEF); t1->Branch("jet_nEmEF",&jet_nEmEF);
  t1->Branch("jet_nHEF",&jet_nHEF); t1->Branch("jet_nMult",&jet_nMult);
  t1->Branch("jet_eEF",&jet_eEF); t1->Branch("jet_pEF",&jet_pEF);
  t1->Branch("njets",&njets); t1->Branch("nelectron",&nel);
  t1->Branch("el_pt",&el_pt); t1->Branch("el_eta",&el_eta);
  t1->Branch("el_phi",&el_phi); t1->Branch("el_charge",&el_charge);
  t1->Branch("el_dxy",&el_dxy); t1->Branch("el_dz",&el_dz);
  t1->Branch("el_edxy",&el_edxy); t1->Branch("el_edz",&el_edz);
  t1->Branch("el_vx",&el_vx); t1->Branch("el_vy",&el_vy);
  t1->Branch("el_vz",&el_vz); t1->Branch("el_mva_out",&el_mva_out);
  t1->Branch("el_mva_iso",&el_mva_iso); t1->Branch("el_iso",&el_iso);
  t1->Branch("el_veto",&el_veto); t1->Branch("el_soft",&el_soft);
  t1->Branch("el_medium",&el_medium); t1->Branch("el_tight",&el_tight);
  t1->Branch("el_mva_map_value",&el_mva_map_value);
  t1->Branch("HLT1Obj_pt_eta_phi_charge",&tr1_obj_pt_eta_phi);
  t1->Branch("HLT2Obj_pt_eta_phi_charge",&tr2_obj_pt_eta_phi);
  t1->Branch("HLT3Obj_pt_eta_phi_charge",&tr3_obj_pt_eta_phi);
  t1->Branch("HLT4Obj_pt_eta_phi_charge",&tr4_obj_pt_eta_phi);
  t1->Branch("HLT5Obj_pt_eta_phi_charge",&tr5_obj_pt_eta_phi);
  t1->Branch("HLT6Obj_pt_eta_phi_charge",&tr6_obj_pt_eta_phi);
  t1->Branch("HLT7Obj_pt_eta_phi_charge",&tr7_obj_pt_eta_phi);
  t1->Branch("HLT8Obj_pt_eta_phi_charge",&tr8_obj_pt_eta_phi);
  t1->Branch("l1muon_pt",&l1muon_pt); t1->Branch("l1muon_eta",&l1muon_eta);
  t1->Branch("l1muon_phi",&l1muon_phi); t1->Branch("l1muon_iso",&l1muon_iso);
  t1->Branch("l1muon_qual",&l1muon_qual);
  t1->Branch("l1jet_pt",&l1jet_pt); t1->Branch("l1jet_eta",&l1jet_eta);
  t1->Branch("l1jet_phi",&l1jet_phi); t1->Branch("l1jet_iso",&l1jet_iso);
  t1->Branch("l1jet_qual",&l1jet_qual); t1->Branch("l1met",&l1met);
  t1->Branch("l1met_phi",&l1met_phi); t1->Branch("l1ht",&l1ht);
  t1->Branch("l1hf_met",&l1hf_met); t1->Branch("l1hf_met_phi",&l1hf_met_phi);
  //  t1->Branch("muon_DCA",&muon_DCA);
  t1->Branch("muon_trkpt",&muon_trkpt);
  t1->Branch("muon_trketa",&muon_trketa);
  t1->Branch("muon_trkphi",&muon_trkphi);
  t1->Branch("el_trkpt",&el_trkpt);
  t1->Branch("el_trketa",&el_trketa);
  t1->Branch("el_trkphi",&el_trkphi);
  
  t1->Branch("pvertex_x",&pvertex_x);
  t1->Branch("pvertex_y",&pvertex_y);
  t1->Branch("pvertex_z",&pvertex_z);
  t1->Branch("pvertex_ex",&pvertex_ex);
  t1->Branch("pvertex_ey",&pvertex_ey);
  t1->Branch("pvertex_ez",&pvertex_ez);
  t1->Branch("ntracks",&ntracks);
  t1->Branch("track_pt",&track_pt);
  t1->Branch("track_eta",&track_eta);
  t1->Branch("track_phi",&track_phi);
  t1->Branch("track_norm_chi2",&track_norm_chi2);
  t1->Branch("track_charge",&track_charge);
  t1->Branch("track_dxy",&track_dxy);
  t1->Branch("track_dz",&track_dz);
  t1->Branch("track_validhits",&track_validhits);
  t1->Branch("track_losthits",&track_losthits);
  t1->Branch("track_fromPV",&track_fromPV);
  t1->Branch("track_highPurity",&track_highPurity);

}

// ------------ method called once each job just after ending the event loop  ------------
template<typename T1>
void 
TriggerAnalyzer<T1>::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
template<typename T1>
void
TriggerAnalyzer<T1>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
typedef TriggerAnalyzer<reco::RecoEcalCandidate> TriggerAnalyzerb;
DEFINE_FWK_MODULE(TriggerAnalyzerb);
