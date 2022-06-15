/* ADAGE frontend class
 * A - Analysis in a
 * D - Directed
 * A - Acyclic
 * G - Graph
 * E - Environment
 *
 * - Associates DAG paths to unique objects
 * - User-friendly higher order functions (`Payload`, `MultiPayload`, ...) for operating
 *   on the stored objects
 */
#ifndef Adage_
#define Adage_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>

// ROOT
#include "TSystem.h"
#include "TObject.h"
#include "TNamed.h"
#include "TString.h"
#include "TRegexp.h"
#include "TFile.h"
#include "TKey.h"

// sidis-eic
#include "DAG.h"
#include "BinSet.h"
#include "Node.h"
#include "NodePath.h"


class Adage : public DAG
{
  public:
    Adage();
    ~Adage();

    // build the DAG from specified bin scheme
    /* - payload classes must override this, following this general format:
     *   - first, call BuildDAG(binSchemes)
     *   - next, stage a LeafOp lambda
     *     - instantiate your payload object, give it a unique name (`CreatePayloadName`, `CreatePayloadTitle`), etc.
     *     - last line of LeafOp lambda should be InsertPayloadData
     *   - last, call ExecuteAndClearOps();
     */
    virtual void Build(std::map<TString,BinSet*> binSchemes) = 0;

    // build the DAG from specified bin scheme
    void BuildDAG(std::map<TString,BinSet*> binSchemes);

    // associate a path to payload data
    template<class PL>
      InsertPayloadData(NodePath *P, PL *PLptr) {
        pathHash.insert(std::pair<std::set<Node*>,PL*>(P->GetBinNodes(),PLptr));
      };

    // build the DAG from ROOT file
    /* - all BinSets will become layers
     * - all objects that have names starting with `pattern` will be associated to DAG paths
     */
    void Build(TFile *rootFile, TString pattern);


    // create a unique name or title for a payload object
    TString CreatePayloadName(TString keyName, NodePath *P);
    TString CreatePayloadTitle(TString keyName, NodePath *P);

    // payload operator, executed on the specified object; see `FormatPayload`
    // for allowed arguments of operator `op`
    template<class O>
      void Payload(O op) {
        auto opFormatted = FormatPayload(op);
        LeafOp( [opFormatted,this](NodePath *P){ opFormatted(this->GetPayloadData(P),P); } );
      };

    // multi-payload operators, for defining multiple payloads for a subloop
    template<class O1, class O2, class O3>
      void MultiPayload(std::vector<TString> layers, O1 opPayload, O2 opBefore, O3 opAfter) {
        auto opFormatted = FormatPayload(opPayload);
        MultiLeafOp(
            layers,
            [opFormatted,this](NodePath *P){ opFormatted(this->GetPayloadData(P),P); },
            opBefore,
            opAfter
            );
      };
    template<class O1, class O2>
      void MultiPayload(std::vector<TString> layers, O1 opPayload, O2 opBefore) { MultiPayload(layers,opPayload,opBefore,[](){}); };
    template<class O1>
      void MultiPayload(std::vector<TString> layers, O1 opPayload) { MultiPayload(layers,opPayload,[](){},[](){}); };

    // format payload operators with proper arguments `(PayLoad object pointer, NodePath pointer)`, to allow easy overloading
    template<class PL>
      static std::function<void(PL*,NodePath*)> FormatPayload(std::function<void(PL*,NodePath*)> op) {
        return op;
      };
    template<class PL>
      static std::function<void(PL*,NodePath*)> FormatPayload(std::function<void(NodePath*,PL*)> op) {
        return [op](PL *H, NodePath *P){ op(P,H); };
      };
    template<class PL>
      static std::function<void(PL*,NodePath*)> FormatPayload(std::function<void(PL*)> op) {
        return [op](PL *H, NodePath *P){ op(H); };
      };
    template<class PL>
      static std::function<void(PL*,NodePath*)> FormatPayload(std::function<void(NodePath*)> op) {
        return [op](PL *H, NodePath *P){ op(P); };
      };
    template<class PL>
      static std::function<void(PL*,NodePath*)> FormatPayload(std::function<void()> op) {
        return [op](PL *H, NodePath *P){ op(); };
      };

    // return payload data `PL` associated with the given NodePath
    template<class PL> PL *GetPayloadData(NodePath *P);
    // if you have a NodePath from another DAG that has the same binning scheme, use GetPayloadDataViaID instead
    template<class PL> PL *GetPayloadDataViaID(NodePath *extP);

  private:
    Bool_t debug;
    std::map<std::set<Node*>,PL*> payloadHash; // map DAG path of bin nodes -> PL*

  ClassDefOverride(Adage,1);
};

#endif
