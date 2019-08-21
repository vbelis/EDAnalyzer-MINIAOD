# EDAnalyzer-MINIAOD
## Introduction
EDAnalyzer package for MINIAOD file format, runs on MC and data. Used in CMS experiment software framework (CMSSW_10_2_6), compiled/built for architecture slc6_amd64_gcc700 in SLC 6 OS, ran on CERN's lxplus cluster.

## Basic Functionalities
Basic Functions used for the pTrel b-Purity study:
### _genJetMuAnalyze_:
Creates a pair of a Jet and a Muon at Generation level. For the pair to be created, the following must be true:
- deltaR(Jet,Muon)<0.4.
- The quark whose fragmentation/hadronization is responsible for the Jet, is an ancestor of the muon.
### _analyze_:
The main function of every CMSSW analyzer, which loops over the events and calls all other defined functions. For the b-Purity study, this function did only two things. The first is to call genJetMuAnalyze in every event and the second is to save all Muon and Jet information of each event at Reconstuction level. Then all info is saved in an output .root file -_MC Sample_- in a TTree form. From there, the analysis continues with the [MC Template Producer & Fitter](https://github.com/vbelis/MC_templateProducer_fitter).

There is also a section of code that is under development. Its purpose is to create a preliminary _MC Template_ with all the necessary information needed for the analysis (e.g. b-Purity study). Then the final MC Template can be created by another program (e.g. the, similarly, under development [Direct MC Template Producer/Analyzer](https://github.com/vbelis/DataAnalyzer_directMCAnalyzer/tree/master/newMC_direct_TriggerAnalyzer_analysis)) , which, just, applies the kinematic cuts and othe criteria.  Thus the purpose is to "eliminate" the need to use a seperate MC Template producer which takes as an input the MC Sample.  However, up until now,for the b-Purity study, this section of code was not used. The refered section of code is:
```c++
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
```

