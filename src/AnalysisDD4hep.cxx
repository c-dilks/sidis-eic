#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TLorentzVector.h>
#include <TMath.h>

#include "AnalysisDD4hep.h"

using std::map;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;

AnalysisDD4hep::AnalysisDD4hep(
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
      // dd4hep-specific settings defaults
      /*// initialize scatt. electron cuts // DEPRECATED; TODO: remove?
        fEThreshold = eleBeamEn_*0.1; // min energy cut
        fIsoR = 1.0;                  // Isolation cone R
        fIsoCut = 0.1;                // 10%
        */
    };

// destructor
AnalysisDD4hep::~AnalysisDD4hep() {
};


void AnalysisDD4hep::process_event()
{
  // setup
  Prepare();

  // read dd4hep tree
  TChain *chain = new TChain("events");
  for(Int_t idx=0; idx<infiles.size(); ++idx) {
    for(std::size_t idxF=0; idxF<infiles[idx].size(); ++idxF) {
      std::cout << "Adding " << infiles[idx][idxF] << " with " << inEntries[idx][idxF] << std::endl;
      chain->Add(infiles[idx][idxF].c_str(), inEntries[idx][idxF]);
    }
  }

  TTreeReader tr(chain);

  // Truth
  TTreeReaderArray<Int_t>    mcparticles_ID(tr,        "mcparticles.ID");
  TTreeReaderArray<Int_t>    mcparticles_pdgID(tr,     "mcparticles.pdgID");
  TTreeReaderArray<Double_t> mcparticles_psx(tr,       "mcparticles.ps.x");
  TTreeReaderArray<Double_t> mcparticles_psy(tr,       "mcparticles.ps.y");
  TTreeReaderArray<Double_t> mcparticles_psz(tr,       "mcparticles.ps.z");
  TTreeReaderArray<Int_t>    mcparticles_status(tr,    "mcparticles.status");
  TTreeReaderArray<Int_t>    mcparticles_genStatus(tr, "mcparticles.genStatus");
  TTreeReaderArray<Double_t> mcparticles_mass(tr,      "mcparticles.mass");

  // Reco
  TTreeReaderArray<Int_t> ReconstructedParticles_pid(tr,   "ReconstructedParticles.pid");
  TTreeReaderArray<float> ReconstructedParticles_energy(tr,  "ReconstructedParticles.energy");
  TTreeReaderArray<float> ReconstructedParticles_p_x(tr,     "ReconstructedParticles.p.x");
  TTreeReaderArray<float> ReconstructedParticles_p_y(tr,     "ReconstructedParticles.p.y");
  TTreeReaderArray<float> ReconstructedParticles_p_z(tr,     "ReconstructedParticles.p.z");
  TTreeReaderArray<float> ReconstructedParticles_p(tr,       "ReconstructedParticles.momentum");
  TTreeReaderArray<float> ReconstructedParticles_th(tr,      "ReconstructedParticles.direction.theta");
  TTreeReaderArray<float> ReconstructedParticles_phi(tr,     "ReconstructedParticles.direction.phi");
  TTreeReaderArray<float> ReconstructedParticles_mass(tr,    "ReconstructedParticles.mass");
  TTreeReaderArray<short> ReconstructedParticles_charge(tr,  "ReconstructedParticles.charge");
  TTreeReaderArray<int>   ReconstructedParticles_mcID(tr,    "ReconstructedParticles.mcID.value");

  TTreeReader::EEntryStatus entrystats = tr.SetEntry(0);

  CalculateEventQ2Weights();

  // counters
  Long64_t nevt, numNoBeam, numEle, numProxMatched;
  nevt = numNoBeam = numEle = numProxMatched = 0;

  // event loop =========================================================
  cout << "begin event loop..." << endl;
  int errorCount=0;
  while(tr.Next())
  {
    if(nevt%10000==0) cout << nevt << " events..." << endl;
    nevt++;
    if(nevt>maxEvents && maxEvents>0) break;

    // test asan ###########################################
    //Node *leak = new Node();
    //###########################################################

    // mcparticles loop (generated truth)
    std::vector<Particles> mcpart;
    double maxP = 0;
    int electronID = 0;
    bool foundBeamElectron = false;
    bool foundBeamIon = false;
    for(int imc=0; imc<mcparticles_pdgID.GetSize(); imc++)
    {
      int pid_ = mcparticles_pdgID[imc];
      int genStatus_ = mcparticles_genStatus[imc]; // genStatus 4: beam particle 1: final state
      double px_ = mcparticles_psx[imc];
      double py_ = mcparticles_psy[imc];
      double pz_ = mcparticles_psz[imc];
      double mass_ = mcparticles_mass[imc]; // in GeV
      double p_ = sqrt(pow(mcparticles_psx[imc],2) + pow(mcparticles_psy[imc],2) + pow(mcparticles_psz[imc],2));

      if(genStatus_ == 1) // final state
      {
        Particles part;
        part.pid = pid_;
        part.vecPart.SetPxPyPzE(px_, py_, pz_, sqrt(p_*p_ + mass_*mass_));
        part.mcID = mcparticles_ID[imc];
        mcpart.push_back(part);

        if(pid_ == 11)
        {
          if(p_ > maxP)
          {
            maxP = p_;
            kinTrue->vecElectron.SetPxPyPzE(
                px_,
                py_,
                pz_,
                sqrt(p_*p_ + mass_*mass_)
                );
            electronID = mcparticles_ID[imc];
          }
        }// if electron
      }//

      else if(genStatus_ == 4) { // beam particles
        if(pid_ == 11) { // electron beam
          if(!foundBeamElectron) {
            foundBeamElectron = true;
            kinTrue->vecEleBeam.SetPxPyPzE(
                px_,
                py_,
                pz_,
                sqrt(p_*p_ + mass_*mass_)
                );
          }
          else { if(++errorCount<100) cerr << "ERROR: Found two beam electrons in one event" << endl; }
        }
        else { // ion beam
          if(!foundBeamIon) {
            foundBeamIon = true;
            kinTrue->vecIonBeam.SetPxPyPzE(
                px_,
                py_,
                pz_,
                sqrt(p_*p_ + mass_*mass_)
                );
          }
          else { if(++errorCount<100) cerr << "ERROR: Found two beam ions in one event" << endl; }
        }
      }
    }//mcparticles loop

    if(!foundBeamElectron || !foundBeamIon) { numNoBeam++; continue; };
    if(errorCount>=100 && errorCount<1000) { cerr << "ERROR: .... suppressing beam finder errors ...." << endl; errorCount=1000; };


    // calculate true DIS kinematics
    kinTrue->CalculateDIS(reconMethod); // generated (truth)

    // collect reconstructed particles
    std::vector<Particles> recopart;
    double hpx=0;
    double hpy=0;
    double hpz=0;
    double hE=0;
    int foundElectron = 0;
    for(int ireco=0; ireco<ReconstructedParticles_pid.GetSize(); ireco++)
    {
      int pid_ = ReconstructedParticles_pid[ireco];

      // pid==0: reconstructed tracks with no matching truth pid
      if(pid_ == 0) continue;


      Particles part;
      part.pid = pid_;
      part.mcID = ReconstructedParticles_mcID[ireco];
      part.charge = ReconstructedParticles_charge[ireco];

      double reco_E = ReconstructedParticles_energy[ireco];
      double reco_px = ReconstructedParticles_p_x[ireco];
      double reco_py = ReconstructedParticles_p_y[ireco];
      double reco_pz = ReconstructedParticles_p_z[ireco];
      double reco_mass = ReconstructedParticles_mass[ireco];
      double reco_p = sqrt(reco_px*reco_px + reco_py*reco_py + reco_pz*reco_pz);

      part.vecPart.SetPxPyPzE(
          reco_px,
          reco_py,
          reco_pz,
          sqrt(reco_p*reco_p + reco_mass*reco_mass));

      recopart.push_back(part);

      hpx += reco_px;
      hpy += reco_py;
      hpz += reco_pz;
      hE += reco_E;

      // find scattered electron
      if(pid_ == 11 && part.mcID == electronID)
      {
        foundElectron = 1;
        kin->vecElectron.SetPxPyPzE(reco_px,
            reco_py,
            reco_pz,
            sqrt(reco_p*reco_p + reco_mass*reco_mass));
      }
    }//reco loop

    // skip the event if the scattered electron is not found
    // and we need it to calculate the DIS kinematics
    if(foundElectron < 1){
      numEle++;
      if(reconMethod != "JB")
        continue;
    }

    kin->vecHadron.SetPxPyPzE(hpx, hpy, hpz, hE);
    kin->vecHadron -= kin->vecElectron;

    //Hadronic reconstruction
    TLorentzVector head_vecElectron;
    TLorentzVector head_vecHadron;
    kin->TransformToHeadOnFrame(kin->vecElectron,head_vecElectron);
    kin->TransformToHeadOnFrame(kin->vecHadron,head_vecHadron);
    kin->sigmah = (head_vecHadron.E() - head_vecHadron.Pz());
    kin->Pxh = head_vecHadron.Px();
    kin->Pyh = head_vecHadron.Py();

    // calculate DIS kinematics
    kin->CalculateDIS(reconMethod); // reconstructed

    // calculate hadron kinematics
    for(auto part : recopart)
    {
      int pid_ = part.pid;
      int mcid_ = part.mcID;

      // final state cut
      // - check PID, to see if it's a final state we're interested in for
      //   histograms; if not, proceed to next track
      auto kv = PIDtoFinalState.find(pid_);
      if(kv!=PIDtoFinalState.end()) finalStateID = kv->second; else continue;
      if(activeFinalStates.find(finalStateID)==activeFinalStates.end()) continue;

      kin->vecHadron = part.vecPart;
      kin->CalculateHadronKinematics();

      // find the matching truth information
      // using mcID
      if(mcid_ > 0){
        for(auto imc : mcpart)
        {
          if(mcid_ == imc.mcID)
          {
            kinTrue->vecHadron = imc.vecPart;
            break;
          }
        }
      }
      else{
        // give it another shot: proximity matching
        double mineta = 4.0;
        numProxMatched++;
        for(int imc=0; imc<(int)mcpart.size(); imc++)
        {
          if(pid_ == mcpart[imc].pid)
          {
            double deta = abs(kin->vecHadron.Eta() - mcpart[imc].vecPart.Eta());
            if( deta < mineta )
            {
              mineta = deta;
              kinTrue->vecHadron = mcpart[imc].vecPart;
            }
          }
        }
      }

      kinTrue->CalculateHadronKinematics();

      // asymmetry injection
      //kin->InjectFakeAsymmetry(); // sets tSpin, based on reconstructed kinematics
      //kinTrue->InjectFakeAsymmetry(); // sets tSpin, based on generated kinematics
      //kin->tSpin = kinTrue->tSpin; // copy to "reconstructed" tSpin

      Double_t Q2weightFactor = GetEventQ2Weight(kinTrue->Q2, inLookup[chain->GetTreeNumber()]);
      wTrack = Q2weightFactor * weight->GetWeight(*kinTrue);
      wTrackTotal += wTrack;

      // fill track histograms in activated bins
      FillHistosTracks();

      // fill simple tree
      // - not binned
      // - `activeEvent` is only true if at least one bin gets filled for this track
      // - TODO [critical]: add a `finalState` cut (also needed in AnalysisDelphes)
      if( writeSimpleTree && activeEvent ) ST->FillTree(wTrack);

    }//hadron loop

  }// tree reader loop

  cout << "end event loop" << endl;
  // event loop end =========================================================

  // finish execution
  Finish();

  // final printout
  cout << "Total number of scattered electrons found: " << numEle << endl;
  if(numNoBeam>0)
    cerr << "WARNING: skipped " << numNoBeam << " events which had no beam particles" << endl;
  if(numProxMatched>0)
    cerr << "WARNING: " << numProxMatched << " recon. electrons were proximity matched to truth (when mcID match failed)" << endl;

}//execute




// DEPRECATED electron finder methods; TODO: remove?
/*
double AnalysisDD4hep::isolation(double cone_theta, double cone_phi, std::vector<Clusters*> cluster_container, double E_threshold)
{
  double cone_eta = -log(tan(cone_theta/2.0));
  double cone_iso = 0.0;

  // loop over clusters
  for(auto icl : cluster_container)
    {
      auto clusters = *icl;
      double clus_e     = clusters.E;
      double clus_phi   = clusters.phi;
      double clus_theta = clusters.theta;
      double clus_eta   = -log(tan(clus_theta/2.0));
      double clus_pt    = clus_e*sin(clus_theta);
      
      double dphi = clus_phi - cone_phi;
      if( dphi > TMath::Pi() ) dphi = dphi - 2*TMath::Pi();
      double dr = sqrt( pow(dphi,2) + pow((clus_eta - cone_eta), 2) );

      // get E_cone
      if(dr < fIsoR)
	{
	  cone_iso += clus_e;
	}
    }// for loop

  return cone_iso;
}

int AnalysisDD4hep::find_electron(std::vector<Clusters*> ecal_cluster, std::vector<Clusters*> hcal_cluster, double e_threshold)
{

  // init 
  double ptmax = 0;
  double electron_iso = 999;
  int index_max = -999;
  
  int index = 0;
  for(auto icl : ecal_cluster)
    {
      index++;
      auto icluster = *icl;
      double clus_E = icluster.E;
      double clus_theta = icluster.theta;
      double clus_phi   = icluster.phi;
      double clus_eta   = -log(tan(clus_theta/2.0));
      double clus_pt    = clus_E * sin(clus_theta);

      if(clus_E < e_threshold)
	continue;
      if(icluster.theta * 180./TMath::Pi() < 2)
	continue;

      // FIXME hard-coded threshold (currently not used)
      double clus_ecal_iso = isolation(clus_theta, clus_phi, ecal_cluster, 0.1);
      double clus_hcal_iso = isolation(clus_theta, clus_phi, hcal_cluster, 0.1);
      double clus_iso = clus_ecal_iso + clus_hcal_iso - clus_E; // subtract the electron energy

      if(clus_iso > clus_E*fIsoCut)
	continue;

      if(clus_pt > ptmax)
	{
	  ptmax = clus_pt;
	  electron_iso = clus_iso;
	  index_max = index;
	}
      }// loop over ecal clusters
  
  return index_max-1;
}
*/
