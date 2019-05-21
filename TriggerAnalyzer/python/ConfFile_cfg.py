
import FWCore.ParameterSet.Config as cms
IsData=False
Run="A"
output_path="/afs/cern.ch/work/v/vbelis/private/test_for_no04DR_210519.root"
RecoBtoKLepLep=False
RecoBtoKstarLepLep=True
SkipEventWithNoRecoB=False
MuonsOnly=True
ElectronsOnly=False
addlostTrk=True
saveTrk=False
#Nentries=629771
Nentries=10
File=['/store/data/Run2018B/ParkingBPH5/MINIAOD/PromptReco-v1/000/317/650/00000/321646CB-F76E-E811-91FF-FA163EE936A8.root']
############
if RecoBtoKLepLep : 
   print "reconstructing B->Kll channel"
if RecoBtoKstarLepLep :
   print "reconstructing B->K*ll->Kpill channel"

Addel=True
Onlyel=False
if MuonsOnly and not ElectronsOnly:
  Addel=False 
elif ElectronsOnly and not MuonsOnly:
  Onlyel=True
elif ElectronsOnly and MuonsOnly:
   print "warning both chanels will be kept"


SkipNoKLL=False
if RecoBtoKLepLep and SkipEventWithNoRecoB :
   SkipNoKLL=True
#print SkipNoKLL
SkipNoKsLL=False
if RecoBtoKstarLepLep and SkipEventWithNoRecoB :
   SkipNoKsLL=True

if Run=="A":
   n1="HLT_Mu9_IP6_part" ; n2="HLT_Mu8p5_IP3p5" ; n3 ="HLT_Mu10p5_IP3p5" ; 
   n4="HLT_Mu8_IP3"; n5="empty" ; n6="empty" ; n7="empty" ; n8="empty" 
elif Run=="B":
   n1="HLT_Mu9_IP6_part" ; n2="HLT_Mu9_IP5" ; n3 ="HLT_Mu7_IP4" ; 
   n4="HLT_Mu8_IP3"; n5="HLT_Mu12_IP6" ; n6="empty" ; n7="empty" ; n8="empty"
elif Run=="D":
   n1="HLT_Mu9_IP6_part" ; n2="HLT_Mu9_IP5" ; n3="HLT_Mu7_IP4"; n4="HLT_Mu8_IP3"; n5="HLT_Mu12_IP6" ; n6="HLT_Mu9_IP4" ; n7="HLT_Mu8_IP6"; n8="HLT_Mu8_IP5"
else:
   n1="empty" ; n2=n1 ;n3=n1 ; n4=n1 ; n5=n1 ; n6=n1 ; n7=n1 ; n8=n1
   

globaltag='102X_upgrade2018_realistic_v15' 
L1save=False ; HLTsave=False ; HLTfired=False
if IsData:
   print "We have established we Run on data"
   globaltag='101X_dataRun2_Prompt_v11'
   L1save=True ; HLTsave=True ; HLTfired=True
else:
   print "We have established we Run on MC"
print "Run parameters ",globaltag," save L1 ",L1save," HLT ",HLTsave," save ev. if path fired ",HLTfired

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
#tracks from pf
process.load('PhysicsTools.PatAlgos.slimming.unpackedTracksAndVertices_cfi')

#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
#    ignoreTotal = cms.untracked.int32(1)
#)

process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag,globaltag, '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(Nentries) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#    File
'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/120000/FF72630A-46DB-324E-86D7-1503611633F5.root',
'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/120000/FF41E446-BAC2-0F4B-BA11-455B299A1D30.root',
'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/120000/FED7DFCE-B7C5-D847-9FE0-06017DA28314.root',
'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/120000/FEA92678-B29A-C84B-A127-DF98462ED421.root',
'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/120000/FDB6508C-26DE-A444-8CED-1F8193AF803A.root',
'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/120000/FD73A93B-3584-2844-B108-8AC2B309785D.root',
'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/120000/FBE5B7D7-1BAE-2046-8018-56A201153B01.root',
'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/120000/FB50A629-BDC6-E841-A6A6-5B9FB162AF21.root',
'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/120000/FB466A4D-DC5F-D44E-BE47-15AD57072741.root',
'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/120000/FAFC97C6-7716-DA4C-BC25-0528B0E57EEB.root',
'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/120000/FAC93801-1BF4-3C49-9A7A-BBCFF396DA80.root',
'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/120000/F9F46D55-107D-AE43-948A-CD7875489C82.root',
'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/120000/F7B51987-22BC-9740-B5B4-6F8A047A9A89.root',
'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/120000/F7060754-93D8-3A48-8274-38A8CBB1CE26.root',
'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/120000/F6D29DDE-6A17-1649-B9FD-0C5568AF7EB3.root',
'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/120000/F67AF596-2279-D84E-A7CA-79D25560824D.root',
'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/120000/F6134C44-7F43-5C4E-AAC0-C06B29EE3EA6.root',
'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/120000/F5815F2A-85BC-0D4A-82E1-852941270B3E.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/70000/D79A9D9B-E255-9B4F-8CAC-C9937829D398.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/70000/BA314D29-FE76-5C42-8563-F04C4AA80771.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/60000/FE478DDA-CDA2-9C41-9A8C-18DAAA5193CF.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/60000/DEE2E14E-8A89-2A40-9CC2-DBAFC51449C3.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/60000/DC918904-9517-EF4B-AFEE-B0D2F5F5F5AB.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/60000/D0624E72-FBB1-5D4F-A7E8-30D9CB25AC11.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/60000/BC4D607E-911B-4A4F-B19A-04A919D73FD2.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/60000/B91947AB-34AE-E641-9029-C79F59AC3572.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/60000/AD4C6F5F-9BDE-EF4F-BDEF-1602190C05E7.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/60000/ABF90A6C-1BCD-D843-9839-FB7C38D98048.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/60000/AB5A9476-3D6F-9A47-8B36-384E51CB0A79.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/60000/A382EDB7-5E0D-084C-9EF6-59B1AA5C833F.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/60000/8E06C481-EC16-8740-854A-F96D07A36418.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/60000/89D8A3DA-FBCA-BC4B-97E8-EFE087DE20BE.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/60000/8280654F-16EB-2944-B475-C524AC666837.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/60000/7F8CF172-C038-044A-BAD6-CF0511F455D5.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/60000/7604A01F-4365-F04C-9035-922BDE22D722.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/60000/6AD174AC-FEA7-5445-A7D2-4B797D4DA3AD.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/60000/60CEB336-A170-934D-A9CF-56B43460024F.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/60000/5B5E3617-66F3-F647-BDB0-50DA256D637C.root',
#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-20to30_MuEnrichedPt5_TuneCP5_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v4/120000/07E8D448-ABDE-B447-AC2C-B828ED0E4443.root'
#for aod comparison
# '/store/data/Run2018A/ParkingBPH1/MINIAOD/14May2018-v1/30000/F637E740-F259-E811-B4BB-782BCB3B1A58.root'
#'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_MINIAODSIM_18_03_21/180319_172427/0001/BToKee_MINIAODSIM_1727.root'
#aod comp 2
#'/store/data/Run2018B/ParkingBPH5/MINIAOD/PromptReco-v1/000/317/696/00000/F216061E-AF70-E811-B98A-FA163EC375AA.root'
#'/store/user/tstreble/BToKmm_Pythia/BToKmm_Pythia_MINIAODSIM_18_03_26/180326_093119/0000/BToKmm_MINIAODSIM_99.root'
#fastsim
#"file:/afs/cern.ch/work/g/gkaratha/private/SUSYCMG/BtoKppmunu_production/CMSSW_9_3_6/src/step3_PAT.root"
#'/store/group/cmst3/user/gkaratha/PAT_FastSim_GenBToKEE_NonD_muFilter_ForProbe_try1/CRAB_UserFiles/crab_PAT_FastSim_GenBToKEE_NonD_muFilter_ForProbe_try1/190129_202905/0000/step3_PAT_2.root'
#'/store/mc/RunIIAutumn18MiniAOD/BuToK_Toee_MuFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PUPoissonAve20_102X_upgrade2018_realistic_v15-v3/90000/FE81DFA9-EA47-764C-842D-1F7A31327500.root'

  ),
   secondaryFileNames=cms.untracked.vstring(
),
 # eventsToProcess=cms.untracked.VEventRange('316187:796:MIN-316187:796:MAX'),
 #  eventsToProcess=cms.untracked.VEventRange('317696:399:MIN-317696:399:MAX')
  #eventsToProcess=cms.untracked.VEventRange('317696:252:329504213-317696:252:329504215'),
   inputCommands=cms.untracked.vstring(
                  'keep *',
                  'drop *_ctppsPixelClusters_*_*',
                  
          )

)
'''process.selectedPFCandidatesHP = cms.EDFilter("PATPackedCandidateSelector",
    src = cms.InputTag("packedPFCandidates"),
    cut = cms.string("pt() > 0.8 && abs(eta()) < 2.5 && trackHighPurity() > 0")
)'''
#taskB0.add(process.selectedPFCandidatesHP)


from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)
my_id_modules = [
        'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff', 
 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V1_cff', 
]
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)


process.demo = cms.EDAnalyzer('TriggerAnalyzerb',
                              beamSpot = cms.InputTag('offlineBeamSpot'),
                              electrons    = cms.InputTag("slimmedElectrons"),
                              vertices     = cms.InputTag("offlineSlimmedPrimaryVertices"),
                              jets = cms.InputTag("slimmedJets"),
                              photons = cms.InputTag("slimmedPhotons"),
                              eleIdMapVeto = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-veto"),
                              eleIdMapSoft = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-noIso-V1-wpLoose"),
                              eleIdMapMedium = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-noIso-V1-wp90"),
                              eleIdMapTight = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-noIso-V1-wp80"),
                              eleIdMapValue = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV1Categories"),
                              #If you want no L1_Seed, write "default" in the first element and the tree will write the value -100
                               Seed=cms.vstring("L1_SingleMu7er1p5","L1_SingleMu8er1p5","L1_SingleMu9er1p5","L1_SingleMu10er1p5","L1_SingleMu12er1p5","L1_SingleMu22"),
                              HLTPath=cms.vstring(n1,n2,n3,n4,n5,n6,n7,n8),
################################NORMALLY USE THIS####################### 
                              triggerresults = cms.InputTag("TriggerResults::HLT"),
                              triggerobjects = cms.InputTag('slimmedPatTrigger','','RECO'),                   
                              muons=cms.InputTag("slimmedMuons"),
                              met=cms.InputTag("slimmedMETs"),
                              l1seed=cms.InputTag("gtStage2Digis::RECO"),
                              l1met=cms.InputTag('caloStage2Digis','EtSum','RECO'), 
                              l1muons=cms.InputTag("gmtStage2Digis","Muon","RECO"),
                              l1jets=cms.InputTag('caloStage2Digis','Jet','RECO'),                            
                              packed = cms.InputTag("packedGenParticles"),
                              pruned = cms.InputTag("prunedGenParticles"),
                              PFCands=cms.InputTag("packedPFCandidates"),
                              losttracks=cms.InputTag("lostTracks"),
                              genjets=cms.InputTag("slimmedGenJets"),

 #                              tracks=cms.InputTag("unpackedTracksAndVertices"),

                               RunParameters = cms.PSet(
      Data= cms.bool(IsData),SaveTracks=cms.bool(saveTrk),
      SaveHLT=cms.bool(HLTsave),SaveL1=cms.bool(L1save),
      SaveResultsOnlyIfAPathFired=cms.bool(HLTsave),
      ReconstructBMuMuK=cms.bool(RecoBtoKLepLep),
      ReconstructBMuMuKstar=cms.bool(RecoBtoKstarLepLep),
      MuonPtCutForB=cms.double(1.5),TrackPtCutForB=cms.double(1.3),
      PointingConstraint=cms.bool(False),CosThetaCutPointCons=cms.bool(False),
      CosThetaCut=cms.double(-1),UseClosestVertex=cms.bool(False),
      ProbBMuMuKcut=cms.double(-0.001),
      SkipEventWithNoBToMuMuK=cms.bool(SkipNoKLL),UseBeamspot=cms.bool(False),
      AddeeK=cms.bool(Addel),MLLmax_Cut=cms.double(5),MLLmin_Cut=cms.double(0),
      MBmin_Cut=cms.double(4.5),
      MBmax_Cut=cms.double(6),EtaTrk_Cut=cms.double(2.5),
      MKstarMin_Cut=cms.double(0.742),MKstarMax_Cut=cms.double(1.042),
      LepTrkExclusionCone=cms.double(0.005),AddLostTracks=cms.bool(addlostTrk),
      RefitTracks=cms.bool(False),RefitMuTracksOnly=cms.bool(True),
      OnlyKee=cms.bool(Onlyel),UsePFeForCos=cms.bool(True),
      SkipEventWithNoBToMuMuKstar=cms.bool(SkipNoKsLL)
  ),
)

process.load( "HLTrigger.HLTanalyzers.hlTrigReport_cfi" )
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    #SkipEvent = cms.untracked.vstring('ProductNotFound')
)
process.hlTrigReport.HLTriggerResults   = cms.InputTag("TriggerResults", "", "HLT")


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(output_path)
                                   )
process.fevt = cms.OutputModule("PoolOutputModule",
   # SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("path")),
    outputCommands = cms.untracked.vstring(#"drop *",
    ),
    fileName = cms.untracked.string("edm_output.root"))

#process.p = cms.Path(process.egmGsfElectronIDSequence)#* process.demo)
process.p = cms.Path(
   process.egmGsfElectronIDSequence   
  # +process.unpackedTracksAndVertices
#+process.SecondaryVerticesFromHighPurityTracks
 #  +process.selectedPFCandidatesHP
   +process.demo
   #+process.hlTrigReport
  
   )
   
#process.endjob=cms.EndPath(process.fevt)
    #samples
    #'/store/mc/RunIIAutumn18MiniAOD/BuToK_ToMuMu_MuFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PUPoissonAve20_102X_upgrade2018_realistic_v15-v2/80000/FBEB6F2C-3302-9C4E-9E9D-F253120EC027.root'

