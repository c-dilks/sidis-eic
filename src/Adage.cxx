#include "Adage.h"

ClassImp(Adage)

// default constructor
Adage::Adage()
  : debug(false)
{
  InitializeDAG();
};


// build the DAG from specified bin scheme
void Adage::BuildDAG(std::map<TString,BinSet*> binSchemes) {
  // initialize DAG and payloadHash
  InitializeDAG();
  payloadHash.clear();
  // add the finalState layer first, if it exists
  try { 
    BinSet *finalLayer = binSchemes.at("finalState");
    if(finalLayer->GetNumBins()>0) AddLayer(finalLayer);
  } catch(const std::out_of_range &ex) {
    std::cerr << "WARNING: no finalState bins defined" << std::endl;
  };
  // add one layer for each BinSet with nonzero bins defined
  for(auto kv : binSchemes) {
    if(debug) std::cout << "add BinSet " << kv.first << " to Adage" << std::endl;
    if(kv.first=="finalState") continue;
    BinSet *binScheme = kv.second;
    if(binScheme->GetNumBins()>0) AddLayer(binScheme);
  };
};


// build the DAG from ROOT file
void Adage::Build(TFile *rootFile, TString pattern) {
  // initialize DAG and payloadHash, and read rootFile keys
  InitializeDAG();
  payloadHash.clear();
  TListIter nextKey(rootFile->GetListOfKeys());
  TString keyname;
  // add each BinSet as a new layer
  while(TKey *key = (TKey*)nextKey()) {
    keyname = TString(key->GetName());
    if(keyname.Contains(TRegexp("^binset__"))) {
      if(debug) std::cout << "READ LAYER " << keyname << std::endl;
      BinSet *B = (BinSet*)key->ReadObj();
      if(B->GetNumBins()>0) AddLayer(B);
    };
  };
  nextKey.Reset();
  // add each found to histMap
  while(TKey *key = (TKey*)nextKey()) {
    keyname = TString(key->GetName());
    if(keyname.Contains(TRegexp(TString("^")+pattern)) {
      // get NodePath from matched name
      if(debug) std::cout << "READ HISTOS " << keyname << std::endl;
      NodePath P;
      P.nodes.insert(GetRootNode());
      P.nodes.insert(GetLeafNode());
      TString tokID;
      Ssiz_t tf=0;
      while(keyname.Tokenize(tokID,tf,"__")) {
        if(tokID==pattern) continue;
        Node *N = GetNode(tokID);
        if(N) P.nodes.insert(N);
        else {
          std::cerr << "ERROR: mismatch of Node \"" << tokID << "\" between Payload Object and BinSets" << std::endl;
          return;
        };
      };
      // append to `payloadHash`
      if(debug) std::cout << "-> PATH: " << P.PathString() << std::endl;
      payloadHash.insert(std::pair<std::set<Node*>,PL*>(P.GetBinNodes(),(PL*)key->ReadObj()));
    };
  };
};


// create a unique name or title for a payload object
TString CreatePayloadName(TString keyName, NodePath *P) {
  TString ret = keyName;
  for(Node *N : P->GetSortedBins()) {
    ret += "__" + N->GetID();
  };
  return ret;
};
TString CreatePayloadTitle(TString keyName, NodePath *P) {
  TString ret = keyName;
  for(Node *N : P->GetSortedBins()) {
    ret += N->GetCut()->GetCutTitle() + ", ";
  };
  ret(TRegexp(", $")) = "";
  return ret;
};

// return payload object associated with the given NodePath
PL *Adage::GetPayloadData(NodePath *P) {
  PL *ret;
  try { ret = payloadHash.at(P->GetBinNodes()); }
  catch(const std::out_of_range &ex) {
    std::cerr << "ERROR: no Payload Object associated with NodePath "
              << P->PathString() << std::endl;
    return nullptr;
  };
  return ret;
};

// return payload object associated with the given external NodePath, by ID-matching its Nodes to the local DAG's Nodes
PL *Adage::GetPayloadDataViaID(NodePath *extP) {
  NodePath *intP = new NodePath();
  for(auto extN : extP->nodes) {
    auto intN = this->GetNode(extN->GetID(),true); // find internal node by ID-matching external node
    if(intN!=nullptr) intP->nodes.insert(intN);
  }
  return this->GetPayloadData(intP);
};

Adage::~Adage() {
};

