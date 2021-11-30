#ifndef __AnalysisJuggler_H__
#define __AnalysisJuggler_H__

#include "Analysis.h"

class Histos;
class SimpleTree;
class Kinematics;
class BinSet;

class AnalysisJuggler : public Analysis
{
 public:
  AnalysisJuggler(
		  TString infileName_="",
		  Double_t eleBeamEn_=5,
		  Double_t ionBeamEn_=41,
		  Double_t crossingAngle_=0,
		  TString outfilePrefix_=""
		  );
  ~AnalysisJuggler();

  // perform the analysis
  void process_event();
  void Execute() override { process_event(); };

  ClassDefOverride(AnalysisJuggler,1); 
};

#endif /* __AnalysisJuggler_H__ */
