/* ADAGE backend data structure: a Directed Acyclic Graph (DAG)
 * A - Analysis in a
 * D - Directed
 * A - Acyclic
 * G - Graph
 * E - Environment
 */
#ifndef DAG_
#define DAG_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <functional>

// ROOT
#include "TSystem.h"
#include "TObject.h"
#include "TObjArray.h"
#include "TNamed.h"
#include "TString.h"

// sidis-eic
#include "CutDef.h"
#include "BinSet.h"
#include "Node.h"
#include "NodePath.h"

class DAG : public TObject
{
  public:

    DAG();
    ~DAG();

    // ---------------------------------------------------------------------
    /* operator staging functions:
     * - higher order functions to stage operators (lambdas) on the DAG
     * - some functions will rearrange the DAG layers; use the Print methods
     *   to help verify it is behaving how you want
     * - lambdas arguments can be `Node*` or `NodePath*`, or both, or neither
     * - execute the lambdas by calling `ExecuteOps`
     */

    // `Initial`: initial operator, executed at the beginning of DAG::ExecuteOps
    template<class O> void Initial(O op) { GetRootNode()->StageInboundOp(op); };

    // `Final`: final operator, executed at the end of DAG::ExecuteOps
    template<class O> void Final(O op) { GetRootNode()->StageOutboundOp(op); };

    // `Subloop`: add inbound and outbound operators to a control node, given
    // by the first layer in `layers`
    // - use `BeforeSubloop` or `AfterSubloop` if you only need `opBefore` or `opAfter`
    // - the specified layers will be moved to the leaf node
    // - a control node will be added before the first listed layer, and the
    //   lambdas `opBefore` and `opAfter` will be staged to it
    // - control nodes for all other listed layers will be removed (and 
    //   replaced by full patches)
    template<class O1, class O2>
    void Subloop(std::vector<TString> layers, O1 opBefore, O2 opAfter) {
      if(layers.size()==0) {
        std::cerr << "ERROR: empty layers list in DAG::Subloop" << std::endl;
        return;
      };
      bool first = true;
      TString controlID;
      // loop over layers
      for(TString layer : layers) {
        RepatchToLeaf(layer);
        if(first) controlID = layer+"_control";
        else RepatchToFull(layer+"_control");
        first = false;
      };
      // stage lambdas
      GetNode(controlID)->StageInboundOp(opBefore);
      GetNode(controlID)->StageOutboundOp(opAfter);
    };
    template<class O>
    void Subloop(std::vector<TString> layers, O opBefore) { Subloop(layers,opBefore,[](){}); };
    void Subloop(std::vector<TString> layers) { Subloop(layers,[](){},[](){}); };

    // `BeforeSubloop`: stage inbound lambda on a control node, given by the first layer in `layers`
    // - if the control node does not exist, move specified layers to the leaf node,
    //   then create the control node and stage the lambda
    // - if the control node exists, only the lambda is staged and layers are not changed
    template<class O>
    void BeforeSubloop(std::vector<TString> layers, O op) {
      Node *controlNode = GetNode(layers.at(0)+"_control",true);
      if(controlNode) controlNode->StageInboundOp(op);
      else Subloop(layers,op,[](){});
    };

    // `AfterSubloop`: stage inbound lambda on a control node, given by the first
    // layer in `layers` (see `BeforeSubloop`)
    template<class O>
    void AfterSubloop(std::vector<TString> layers, O op) {
      Node *controlNode = GetNode(layers.at(0)+"_control",true);
      if(controlNode) controlNode->StageOutboundOp(op);
      else Subloop(layers,[](){},op);
    };


    // `LeafOp`: add a lambda to the leaf node
    // - this is the main operator, acting on all full root-to-leaf paths
    // - there is no difference between inbound and outbound at the leaf
    template<class O> void LeafOp(O op) { GetLeafNode()->StageInboundOp(op); };


    // `MultiLeafOp`: make a multiControl node, which will alter the LeafOp to
    // be `opLeaf`; you can make multiple `MultiLeafOp`s controlling the same layer,
    // and they will be executed sequentially
    // - the inbound operator will change the current LeafOp, thus it is not recommended
    //   to use `MultiLeafOp` if you also use `LeafOp`; moreover, do not define `MultiLeafOp`s
    //   to control different layers, or their LeafOps will overwrite each other
    template<class O1, class O2, class O3>
    void MultiLeafOp(std::vector<TString> layers, O1 opLeaf, O2 opBefore, O3 opAfter) {
      if(layers.size()==0) {
        std::cerr << "ERROR: empty layers list in DAG::MultiLeafOp" << std::endl;
        return;
      };
      // count how many multi-control nodes exist for the first layer in `layers`
      TString controlVar = layers.at(0);
      Int_t nMulti = 0;
      for(auto kv : nodeMap) {
        auto N = kv.second;
        if( N->GetNodeType()==NT::multi && N->GetVarName()==controlVar) nMulti++;
      };
      TString multiID;
      // if this is the first multi-control node, call Subloop to create an empty control node,
      // then convert it to a multi-control node
      if(nMulti==0) {
        Subloop(layers);
        Node *C = GetNode(controlVar+"_control",true);
        multiID = controlVar+"_multi_0";
        ModifyNode(C,multiID,NT::multi);
      }
      // if this is not the first multi-control node, add a new multi-control node by
      // copying a previous multi-control node's connections
      else {
        multiID = Form("%s_multi_%d",controlVar.Data(),nMulti);
        TString oldMultiID = Form("%s_multi_%d",controlVar.Data(),nMulti-1);
        AddNode(NT::multi,multiID);
        Node *newMultiNode = GetNode(multiID);
        Node *oldMultiNode = GetNode(oldMultiID);
        for(auto N : oldMultiNode->GetInputs())  AddEdge(N,newMultiNode,true);
        for(auto N : oldMultiNode->GetOutputs()) AddEdge(newMultiNode,N,true);
      };
      // stage lambdas: setting the leaf lambda is combined with the inbound lambda
      Node *multiNode = GetNode(multiID);
      auto opBeforeFormatted = Node::FormatOp(opBefore);
      auto opBeforeModded = [this,opBeforeFormatted,opLeaf](Node* N,NodePath *P){
        opBeforeFormatted(N,P);
        LeafOp(opLeaf);
      };
      multiNode->StageInboundOp(opBeforeModded);
      multiNode->StageOutboundOp(opAfter);
    };
    template<class O1, class O2>
    void MultiLeafOp(std::vector<TString> layers, O1 opLeaf, O2 opBefore) { MultiLeafOp(layers,opLeaf,opBefore,[](){}); };
    template<class O1>
    void MultiLeafOp(std::vector<TString> layers, O1 opLeaf) { MultiLeafOp(layers,opLeaf,[](){},[](){}); };


    // end operator templates
    // ---------------------------------------------------------------------


    // initialize DAG, with only the root node connected to the leaf node
    void InitializeDAG();

    // node accessors
    // - search for a node by ID; return nullptr if not found;
    //   set silence=true for no error print out when node not found
    Node *GetNode(TString id_, Bool_t silence=false);
    // get root, leaf, or any other unique node; non-uniqueness 
    // will cause error printout
    Node *GetRootNode();
    Node *GetLeafNode();
    Node *GetUniqueNode(Int_t type_,TString typeStr);
    // layer accessors
    BinSet *GetBinSet(TString varName); // return the BinSet (layer) associated with specified variable name

    // add nodes and edges, which are directional connectors between nodes
    // - nodes are only added if the ID is unique; if not unique,
    //   do nothing (silence=true suppresses error print out)
    void AddNode(Node *N, Bool_t silence=false);
    void AddNode(Int_t nodeType, TString id, Bool_t silence=false);
    void AddEdge(Node *inN, Node *outN, Bool_t silence=false);
    void ModifyNode(Node *N, TString newName, Int_t newType=-1); // rename or repurpose a node

    // add a layer of nodes, fully connected ("patched") to the last layer of the DAG
    // - primary usage is to add a layer of bins from a BinSet
    // - "patch" refers to any set of edges between two layers
    void AddLayer(BinSet *BS);
    void AddLayer(std::vector<Node*> nodes);

    // DAG traversals: iterate through nodes, executing the lambda on each iteration:
    // -- breadth-first: loop through nodes of each layer; lambda operates on the node
    void TraverseBreadth(Node *N, std::function<void(Node*)> lambda);
    void TraverseBreadth(std::function<void(Node*)> lambda) { TraverseBreadth(GetRootNode(),lambda); };
    // -- depth-first: traverse toward the leaf node, iterating over the possible unique paths;
    //    the lambda operates on the unique path to the current node, in addition to the current node
    void TraverseDepth(Node *N, std::function<void(Node*,NodePath*)> lambda, NodePath P=NodePath());
    void TraverseDepth(std::function<void(Node*,NodePath*)> lambda, NodePath P=NodePath()) { 
      TraverseDepth(GetRootNode(),lambda,P); };
    // -- run each node's staged lambdas, while traversing depth first; if `activeNodesOnly`, all nodes
    //    must have `active==true` (by default all nodes are active)
    void ExecuteOps(Bool_t activeNodesOnly=false, Node *N=nullptr, NodePath P=NodePath());
    // -- clear all staged lambdas
    void ClearOps();
    // -- execute then clear
    void ExecuteAndClearOps(Bool_t activeNodesOnly=false) { ExecuteOps(activeNodesOnly); ClearOps(); };

    // patch operations: manipulate connections ("patches") between layers
    /*   - "patch" refers to a set of edges between two layers, which
     *     can be two types:
     *     - "full patch": each node of one layer is connected to each node of the other
     *     - "control patch": each node of one layer is connected to a control node, which
     *       is then connected to each node of the next layer
     */
    // convert control patch(es) to full patch(es); control node(s) will be removed
    void RepatchToFull(TString id);
    void RepatchToFull(Node *N);
    void RepatchAllToFull();
    // re-patch a layer to the current leaf node's output, and re-patch the
    // adjacent layers together; the current leaf node will become a control
    // node, and a new leaf node is created
    void RepatchToLeaf(TString varName);
    // mark a node as visited
    void SetVisited(TString id_) { visitList.push_back(id_); };

    // remove nodes and edges
    void RemoveNode(Node *N);
    void RemoveEdge(Node *inN, Node *outN);

    // set all nodes to active (the default) or inactive (if active_==false)
    void ActivateAllNodes(Bool_t active_=true);

    // print the whole DAG (breadth first or depth first), or
    // specific parts of it
    void PrintBreadth(TString header="DAG");
    void PrintDepth(TString header="DAG");
    void PrintLeafPaths(TString header="DAG leaf paths");


  protected:
    Bool_t Visited(TString id_);

  private:
    Bool_t debug;
    std::map<TString,Node*> nodeMap;
    std::map<TString,BinSet*> layerMap;
    std::vector<TString> visitList;

  ClassDef(DAG,1);
};

#endif
