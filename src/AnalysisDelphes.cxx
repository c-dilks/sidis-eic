#include "AnalysisDelphes.h"

ClassImp(AnalysisDelphes)

using std::map;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;

// constructor
AnalysisDelphes::AnalysisDelphes(
  TString infileName_,
  Double_t eleBeamEn_,
  Double_t ionBeamEn_,
  Double_t crossingAngle_,
  TString outfilePrefix_
) : Analysis(
  infileName_,
  eleBeamEn_,
  ionBeamEn_,
  crossingAngle_,
  outfilePrefix_
) {
  // delphes-specific settings defaults
  /* ... none defined yet ... */
};


//=============================================
// perform the analysis
//=============================================
void AnalysisDelphes::Execute() {

  // setup
  Prepare();

  // read delphes tree
  TChain *chain = new TChain("Delphes");
  for(Int_t idx=0; idx<infiles.size(); ++idx) {
    for(std::size_t idxF=0; idxF<infiles[idx].size(); ++idxF) {
      std::cout << "Adding " << infiles[idx][idxF] << " with " << inEntries[idx][idxF] << std::endl;
      chain->Add(infiles[idx][idxF].c_str(), inEntries[idx][idxF]);
    }
  }
  ExRootTreeReader *tr = new ExRootTreeReader(chain);
  ENT = tr->GetEntries();

  // calculate cross section
  if(maxEvents>0) ENT = maxEvents; // limiter

  // branch iterators
  TObjArrayIter itTrack(tr->UseBranch("Track"));
  TObjArrayIter itElectron(tr->UseBranch("Electron"));
  TObjArrayIter itParticle(tr->UseBranch("Particle"));
  TObjArrayIter itEFlowTrack(tr->UseBranch("EFlowTrack"));
  TObjArrayIter itEFlowPhoton(tr->UseBranch("EFlowPhoton"));
  TObjArrayIter itEFlowNeutralHadron(tr->UseBranch("EFlowNeutralHadron"));
  TObjArrayIter itpfRICHTrack(tr->UseBranch("pfRICHTrack"));
  TObjArrayIter itDIRCepidTrack(tr->UseBranch("barrelDIRC_epidTrack"));
  TObjArrayIter itDIRChpidTrack(tr->UseBranch("barrelDIRC_hpidTrack"));
  TObjArrayIter itBTOFepidTrack(tr->UseBranch("BTOF_eTrack"));
  TObjArrayIter itBTOFhpidTrack(tr->UseBranch("BTOF_hTrack"));
  TObjArrayIter itdualRICHagTrack(tr->UseBranch("dualRICHagTrack"));
  TObjArrayIter itdualRICHcfTrack(tr->UseBranch("dualRICHcfTrack"));

  CalculateEventQ2Weights();

  // event loop =========================================================
  cout << "begin event loop..." << endl;
  int errorCount=0;
  for(Long64_t e=0; e<ENT; e++) {
    if(e>0&&e%10000==0) cout << (Double_t)e/ENT*100 << "%" << endl;
    tr->ReadEntry(e);

    // electron loop
    // - finds max-momentum electron
    itElectron.Reset();
    maxEleP = 0;
    while(Electron *ele = (Electron*) itElectron()) {
      eleP = ele->PT * TMath::CosH(ele->Eta);
      if(eleP>maxEleP) {
        maxEleP = eleP;
	kin->vecElectron.SetPtEtaPhiM(
            ele->PT,
            ele->Eta,
            ele->Phi,
            Kinematics::ElectronMass()
            );
      };
    };
    if(maxEleP<0.001) continue; // no scattered electron found

    // - repeat for truth electron
    itParticle.Reset();
    maxElePtrue = 0;
    bool found_elec = false;
    bool found_ion = false;
    while(GenParticle *part = (GenParticle*) itParticle()){
      if(part->PID == 11 && part->Status == 1){
        elePtrue = part->PT * TMath::CosH(part->Eta);
        if(elePtrue > maxElePtrue){
          maxElePtrue = elePtrue;
          kinTrue->vecElectron.SetPtEtaPhiM(
              part->PT,
              part->Eta,
              part->Phi,
              part->Mass
              );
	  kin->vecElectron.SetPtEtaPhiM(
              part->PT,
              part->Eta,
              part->Phi,
              part->Mass
              );
        };
      };
      if(part->PID == 11 && part->Status == 4){
        if(!found_elec){
          found_elec = true;
          kinTrue->vecEleBeam.SetPtEtaPhiM(
              part->PT,
              part->Eta,
              part->Phi,
              part->Mass
              );
	  kin->vecEleBeam.SetPtEtaPhiM(
              part->PT,
              part->Eta,
              part->Phi,
              part->Mass
              );
        }else{
          if(++errorCount<100) cerr << "ERROR: Found two beam electrons in one event" << endl;
        };
      };
      if(part->PID != 11 && part->Status == 4){
        if(!found_ion){
          found_ion = true;
          kinTrue->vecIonBeam.SetPtEtaPhiM(
              part->PT,
              part->Eta,
              part->Phi,
              part->Mass
              );
	  kin->vecIonBeam.SetPtEtaPhiM(
              part->PT,
              part->Eta,
              part->Phi,
              part->Mass
              );
        }else{
          if(++errorCount<100) cerr << "ERROR: Found two beam ions in one event" << endl;
        };
      };
    };
    if(!found_elec){
      if(++errorCount<100) cerr << "ERROR: Didn't find beam electron in event" << endl;
    };
    if(!found_ion){
      if(++errorCount<100) cerr << "ERROR: Didn't find beam ion in event" << endl;
    }
    if(errorCount>=100 && errorCount<1000) { cerr << "ERROR: .... suppressing beam finder errors ...." << endl; errorCount=1000; };

    // get hadronic final state variables

    kin->GetTrueHFS(itParticle);
    kinTrue->GetTrueHFS(itParticle);

    // calculate DIS kinematics
    if(!(kin->CalculateDIS(reconMethod))) continue; // reconstructed
    if(!(kinTrue->CalculateDIS("ele"))) continue; // generated (truth)

    //    cout << "rec Q2: " << kin->Q2 << " rec x: " << kin->x << " rec y: " << kin->y << " rec s: " << kin->s << endl;
    // cout << "true Q2: " << kinTrue->Q2 << " true x: "	<< kinTrue->x << " true y: " << kinTrue->y << " true s: " << kinTrue->s << endl;
    //cout << "rec Q: " << kin->vecQ.Px() << " " << kin->vecQ.Py() << " " << kin->vecQ.Pz() << " " <<kin->vecQ.E() << endl;
    //cout << "true Q: " << kinTrue->vecQ.Px() << " " << kinTrue->vecQ.Py() << " " << kinTrue->vecQ.Pz() << " " << kinTrue->vecQ.E() << endl;


    //    cout << endl;
    // get vector of jets
    // TODO: should this have an option for clustering method?
    kin->GetJets(itEFlowTrack, itEFlowPhoton, itEFlowNeutralHadron, itParticle);
   
    // track loop - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    itTrack.Reset();
    while(Track *trk = (Track*) itTrack()) {
      //cout << e << " " << trk->PID << endl;

      // final state cut
      // - check PID, to see if it's a final state we're interested in for
      //   histograms; if not, proceed to next track
      // pid = trk->PID; //NOTE: trk->PID is currently not smeared so it just returns the truth-level PID
      pid = kin->getTrackPID( // get smeared PID
          trk,
          itpfRICHTrack,
          itDIRCepidTrack, itDIRChpidTrack,
          itBTOFepidTrack, itBTOFhpidTrack,
          itdualRICHagTrack, itdualRICHcfTrack
          );
      auto kv = PIDtoFinalState.find(pid);
      if(kv!=PIDtoFinalState.end()) finalStateID = kv->second; else continue;
      if(activeFinalStates.find(finalStateID)==activeFinalStates.end()) continue;

      // get parent particle, to check if pion is from vector meson
      GenParticle *trkParticle = (GenParticle*)trk->Particle.GetObject();
      TObjArray *brParticle = (TObjArray*)itParticle.GetCollection();
      GenParticle *parentParticle = (GenParticle*)brParticle->At(trkParticle->M1);
      int parentPID = (parentParticle->PID); // TODO: this is not used yet...

      // calculate hadron kinematics
      kin->hadPID = pid;
      kin->vecHadron.SetPtEtaPhiM(
          trk->PT,
          trk->Eta,
          trk->Phi,
          trk->Mass /* TODO: do we use track mass here ?? */
          );
      GenParticle* trkPart = (GenParticle*)trk->Particle.GetObject();
      kinTrue->hadPID = pid;
      kinTrue->vecHadron.SetPtEtaPhiM(
          trkPart->PT,
          trkPart->Eta,
          trkPart->Phi,
          trkPart->Mass /* TODO: do we use track mass here ?? */
          );
      
      kin->CalculateHadronKinematics();
      kinTrue->CalculateHadronKinematics();

      // asymmetry injection
      //kin->InjectFakeAsymmetry(); // sets tSpin, based on reconstructed kinematics
      //kinTrue->InjectFakeAsymmetry(); // sets tSpin, based on generated kinematics
      //kin->tSpin = kinTrue->tSpin; // copy to "reconstructed" tSpin
  
      // Get index of file that the event comes from.
      Double_t Q2weightFactor = GetEventQ2Weight(kinTrue->Q2, inLookup[chain->GetTreeNumber()]);
      wTrack = Q2weightFactor * weight->GetWeight(*kinTrue);
      wTrackTotal += wTrack;

      // fill track histograms in activated bins
      FillHistosTracks();

      // fill simple tree
      // - not binned
      // - `activeEvent` is only true if at least one bin gets filled for this track
      if( writeSimpleTree && activeEvent ) ST->FillTree(wTrack);

      // tests
      //kin->ValidateHeadOnFrame();

    }; // end track loop


    // jet loop - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    finalStateID = "jet";
    if(activeFinalStates.find(finalStateID)!=activeFinalStates.end()) {

      #if INCCENTAURO == 1
      if(useBreitJets) kin->GetBreitFrameJets(itEFlowTrack, itEFlowPhoton, itEFlowNeutralHadron, itParticle);
      #endif

      Double_t Q2weightFactor = GetEventQ2Weight(kinTrue->Q2, chain->GetTreeNumber());
      wJet = Q2weightFactor * weightJet->GetWeight(*kinTrue); // TODO: should we separate weights for breit and non-breit jets?
      wJetTotal += wJet;

      Int_t nJets;
      if(useBreitJets) nJets = kin->breitJetsRec.size();
      else      nJets = kin->jetsRec.size();

      for(int i = 0; i < kin->jetsRec.size(); i++){

        if(useBreitJets) {
          #if INCCENTAURO == 1
          jet = kin->breitJetsRec[i];
          kin->CalculateBreitJetKinematics(jet);
          #endif
        } else {
          jet = kin->jetsRec[i];
          kin->CalculateJetKinematics(jet);
        };

        // fill jet histograms in activated bins
        FillHistosJets();

      };
    }; // end jet loop

  };
  cout << "end event loop" << endl;
  // event loop end =========================================================


  // finish execution
  Finish();
  //cout << "DEBUG PID in HFS: nSmeared=" << kin->countPIDsmeared << "  nNotSmeared=" << kin->countPIDtrue << endl;
};


// destructor
AnalysisDelphes::~AnalysisDelphes() {
};

