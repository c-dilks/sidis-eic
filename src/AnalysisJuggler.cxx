#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TLorentzVector.h>
#include <TMath.h>

#include "AnalysisJuggler.h"

using std::cout;
using std::cerr;
using std::endl;

AnalysisJuggler::AnalysisJuggler(
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
					      ) {};

// destructor
AnalysisJuggler::~AnalysisJuggler() {
};


void AnalysisJuggler::process_event()
{
  // setup
  Prepare();

  // read output tree
  TChain *chain = new TChain("events");
  for(Int_t idx=0; idx<infiles.size(); ++idx) {
    chain->Add(infiles[idx], inEntries[idx]);
  }

  TTreeReader tr(chain);

  // Truth
  TTreeReaderArray<Int_t>    mcparticles_pdgID(tr,     "mcparticles.pdgID");
  TTreeReaderArray<Double_t> mcparticles_psx(tr,       "mcparticles.ps.x");
  TTreeReaderArray<Double_t> mcparticles_psy(tr,       "mcparticles.ps.y");
  TTreeReaderArray<Double_t> mcparticles_psz(tr,       "mcparticles.ps.z");
  TTreeReaderArray<Int_t>    mcparticles_status(tr,    "mcparticles.status");
  TTreeReaderArray<Int_t>    mcparticles_genStatus(tr, "mcparticles.genStatus");
  TTreeReaderArray<Double_t> mcparticles_mass(tr,      "mcparticles.mass");
  TTreeReaderArray<Int_t> mcparticles_charge(tr,      "mcparticles.charge");

  // Reco
  TTreeReaderArray<int> ReconstructedParticles_pid(tr, "ReconstructedParticles.pid");
  TTreeReaderArray<float> ReconstructedParticles_energy(tr, "ReconstructedParticles.energy");
  TTreeReaderArray<float> ReconstructedParticles_p_x(tr, "ReconstructedParticles.p.x");
  TTreeReaderArray<float> ReconstructedParticles_p_y(tr, "ReconstructedParticles.p.y");
  TTreeReaderArray<float> ReconstructedParticles_p_z(tr, "ReconstructedParticles.p.z");
  TTreeReaderArray<float> ReconstructedParticles_mass(tr, "ReconstructedParticles.mass");

  // Inclusive Kinematics Truth
  TTreeReaderArray<Float_t> InclusiveKinematicsTruth_x(tr, "InclusiveKinematicsTruth.x");
  TTreeReaderArray<Float_t> InclusiveKinematicsTruth_y(tr, "InclusiveKinematicsTruth.y");
  TTreeReaderArray<Float_t> InclusiveKinematicsTruth_Q2(tr, "InclusiveKinematicsTruth.Q2");
  TTreeReaderArray<Float_t> InclusiveKinematicsTruth_W(tr, "InclusiveKinematicsTruth.W");
  TTreeReaderArray<Float_t> InclusiveKinematicsTruth_nu(tr, "InclusiveKinematicsTruth.nu");
  TTreeReaderArray<Int_t> InclusiveKinematicsTruth_scatID_value(tr, "InclusiveKinematicsTruth.scatID.value");
  TTreeReaderArray<Int_t> InclusiveKinematicsTruth_scatID_source(tr, "InclusiveKinematicsTruth.scatID.source");

  // Inclusive Kinematics Electron
  TTreeReaderArray<Float_t> InclusiveKinematicsElectron_x(tr, "InclusiveKinematicsElectron.x");
  TTreeReaderArray<Float_t> InclusiveKinematicsElectron_y(tr, "InclusiveKinematicsElectron.y");
  TTreeReaderArray<Float_t> InclusiveKinematicsElectron_Q2(tr, "InclusiveKinematicsElectron.Q2");
  TTreeReaderArray<Float_t> InclusiveKinematicsElectron_W(tr, "InclusiveKinematicsElectron.W");
  TTreeReaderArray<Float_t> InclusiveKinematicsElectron_nu(tr, "InclusiveKinematicsElectron.nu");
  TTreeReaderArray<Int_t> InclusiveKinematicsElectron_scatID_value(tr, "InclusiveKinematicsElectron.scatID.value");
  TTreeReaderArray<Int_t> InclusiveKinematicsElectron_scatID_source(tr, "InclusiveKinematicsElectron.scatID.source");

  // Inclusive Kinematics JB
  TTreeReaderArray<Float_t> InclusiveKinematicsJB_x(tr, "InclusiveKinematicsJB.x");
  TTreeReaderArray<Float_t> InclusiveKinematicsJB_y(tr, "InclusiveKinematicsJB.y");
  TTreeReaderArray<Float_t> InclusiveKinematicsJB_Q2(tr, "InclusiveKinematicsJB.Q2");
  TTreeReaderArray<Float_t> InclusiveKinematicsJB_W(tr, "InclusiveKinematicsJB.W");
  TTreeReaderArray<Float_t> InclusiveKinematicsJB_nu(tr, "InclusiveKinematicsJB.nu");
  TTreeReaderArray<Int_t> InclusiveKinematicsJB_scatID_value(tr, "InclusiveKinematicsJB.scatID.value");
  TTreeReaderArray<Int_t> InclusiveKinematicsJB_scatID_source(tr, "InclusiveKinematicsJB.scatID.source");
  
  // Inclusive Kinematics DA 
  TTreeReaderArray<Float_t> InclusiveKinematicsDA_x(tr, "InclusiveKinematicsDA.x");
  TTreeReaderArray<Float_t> InclusiveKinematicsDA_y(tr, "InclusiveKinematicsDA.y");
  TTreeReaderArray<Float_t> InclusiveKinematicsDA_Q2(tr, "InclusiveKinematicsDA.Q2");
  TTreeReaderArray<Float_t> InclusiveKinematicsDA_W(tr, "InclusiveKinematicsDA.W");
  TTreeReaderArray<Float_t> InclusiveKinematicsDA_nu(tr, "InclusiveKinematicsDA.nu");
  TTreeReaderArray<Int_t> InclusiveKinematicsDA_scatID_value(tr, "InclusiveKinematicsDA.scatID.value");
  TTreeReaderArray<Int_t> InclusiveKinematicsDA_scatID_source(tr, "InclusiveKinematicsDA.scatID.source");

  TTreeReader::EEntryStatus entrystats = tr.SetEntry(0);

  // event loop =========================================================
  cout << "begin event loop..." << endl;
  Long64_t nevt = 0;
  while(tr.Next())
    {
      if(nevt%10000==0) cout << nevt << " events..." << endl;
      nevt++;      
      if(nevt>maxEvents) break;

      if (InclusiveKinematicsTruth_x.GetSize() != 1) continue;
      for (int imc=0; imc<mcparticles_pdgID.GetSize(); imc++)
	{
	  int pid_ = mcparticles_pdgID[imc];
	  double px_ = mcparticles_psx[imc];
	  double py_ = mcparticles_psy[imc];
	  double pz_ = mcparticles_psz[imc];
	  double mass_ = mcparticles_mass[imc];
	  double p_ = sqrt(px_*px_ + py_*py_ + pz_*pz_);
	  double energy_ = sqrt(p_*p_ + mass_*mass_);

	  if (mcparticles_genStatus[imc] == 4 && pid_ == 11){
	    kinTrue->vecEleBeam.SetPxPyPzE(px_, py_, pz_, energy_);
	  }
	  if (mcparticles_genStatus[imc] == 4 && pid_ == 2212){
	    kinTrue->vecIonBeam.SetPxPyPzE(px_, py_, pz_, energy_);
	  }
	  if (imc == InclusiveKinematicsTruth_scatID_value[0] && pid_ == 11){
	    kinTrue->vecElectron.SetPxPyPzE(px_, py_, pz_, energy_);
	  }
	}
      kinTrue->x = InclusiveKinematicsTruth_x[0];
      kinTrue->y = InclusiveKinematicsTruth_y[0];
      kinTrue->Q2 = InclusiveKinematicsTruth_Q2[0];
      kinTrue->W = InclusiveKinematicsTruth_W[0];
      kinTrue->Nu = InclusiveKinematicsTruth_nu[0];

      if (reconMethod.CompareTo("Ele", TString::kIgnoreCase) == 0){
	if (InclusiveKinematicsElectron_x.GetSize() != 1) continue;
	kin->x = InclusiveKinematicsElectron_x[0];
	kin->y = InclusiveKinematicsElectron_y[0];
	kin->Q2 = InclusiveKinematicsElectron_Q2[0];
	kin->W = InclusiveKinematicsElectron_W[0];
	kin->Nu = InclusiveKinematicsElectron_nu[0];
      }
      else if (reconMethod.CompareTo("DA", TString::kIgnoreCase) == 0){
	if (InclusiveKinematicsDA_x.GetSize() != 1) continue;
	kin->x = InclusiveKinematicsDA_x[0];
	kin->y = InclusiveKinematicsDA_y[0];
	kin->Q2 = InclusiveKinematicsDA_Q2[0];
	kin->W = InclusiveKinematicsDA_W[0];
	kin->Nu = InclusiveKinematicsDA_nu[0];
      }
      else if (reconMethod.CompareTo("JB", TString::kIgnoreCase) == 0){
	if (InclusiveKinematicsJB_x.GetSize() != 1) continue;
	kin->x = InclusiveKinematicsJB_x[0];
	kin->y = InclusiveKinematicsJB_y[0];
	kin->Q2 = InclusiveKinematicsJB_Q2[0];
	kin->W = InclusiveKinematicsJB_W[0];
	kin->Nu = InclusiveKinematicsJB_nu[0];
      }
      else {
	cerr << "ERROR: unknown reconstruction method" << endl;
	return;
      }





    }// tree reader loop

  // finish execution
  Finish();
}//execute
