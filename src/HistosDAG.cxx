#include "HistosDAG.h"

ClassImp(HistosDAG)

// default constructor
HistosDAG::HistosDAG()
  : debug(false)
{
  InitializeDAG();
};


// build the DAG from specified bin scheme
void HistosDAG::Build(std::map<TString,BinSet*> binSchemes) {

  // build the DAG, given the bin scheme
  BuildDAG(binSchemes);

  // leaf operator, to create Histos objects
  LeafOp([this](NodePath *P){
    TString histosN = CreatePayloadName("histos");
    TString histosT = CreatePayloadTitle("");
    if(debug) {
      std::cout << "At path " << P->PathString() << ": ";
      std::cout << "Create " << histosN << std::endl;
    };
    // instantiate Histos object
    Histos *H = new Histos(histosN,histosT);
    // add CutDefs to Histos object
    for(Node *N : P->GetBinNodes()) { H->AddCutDef(N->GetCut()); };
    // append to `histosMap`
    InsertPayloadData(P,H);
  });

  // execution
  if(debug) std::cout << "Begin Histos instantiation..." << std::endl;
  ExecuteAndClearOps();
};


HistosDAG::~HistosDAG() {
};
