// Implement Histos payload for Adage
#ifndef HistosDAG_
#define HistosDAG_

#include "Adage.h"

class HistosDAG : public Adage
{
  public:
    HistosDAG();
    ~HistosDAG();

    // build the DAG from specified bin scheme
    void Build(std::map<TString,BinSet*> binSchemes) override;

  ClassDefOverride(HistosDAG,1);
};

#endif
