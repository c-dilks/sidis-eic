#include "SimpleTree.h"

ClassImp(SimpleTree)

// constructor
SimpleTree::SimpleTree(TString treeName_, Kinematics *K_, Kinematics *Ktrue_) 
  : treeName(treeName_)
  , K(K_)
  , Ktrue(Ktrue_)
{
  // branch names are set to match `brufit` implementation
  // (see `https://github.com/c-dilks/dispin/tree/master/src`)
  T = new TTree(treeName,treeName);
  T->Branch("TrueX",     &(Ktrue->x),     "TrueX/D");
  T->Branch("TrueZ",     &(Ktrue->z),     "TrueZ/D");
  T->Branch("TrueQSq",   &(Ktrue->Q2),    "TrueQSq/D");
  T->Branch("TruePhPerp",&(Ktrue->pT),    "TruePhPerp/D");
  T->Branch("TruePhiH",  &(Ktrue->phiH),  "TruePhiH/D");
  T->Branch("TruePhiS",  &(Ktrue->phiS),  "TruePhiS/D");
  T->Branch("QSq",       &(K->Q2),        "QSq/D");
  T->Branch("X",         &(K->x),         "X/D");
  T->Branch("Y",         &(K->y),         "Y/D");
  T->Branch("Z",         &(K->z),         "Z/D");
  T->Branch("W",         &(K->W),         "W/D");
  T->Branch("MX",        &(K->mX),        "MX/D");
  T->Branch("PhPerp",    &(K->pT),        "PhPerp/D");
  T->Branch("PhiH",      &(K->phiH),      "PhiH/D");
  T->Branch("PhiS",      &(K->phiS),      "PhiS/D");
  T->Branch("PolT",      &(K->polT),      "PolT/D");
  T->Branch("PolL",      &(K->polL),      "PolL/D");
  T->Branch("PolB",      &(K->polBeam),   "PolB/D");
  T->Branch("Depol1",    &(K->depolP1),   "Depol1/D");
  T->Branch("Depol2",    &(K->depolP2),   "Depol2/D");
  T->Branch("Depol3",    &(K->depolP3),   "Depol3/D");
  T->Branch("Depol4",    &(K->depolP4),   "Depol4/D");
  T->Branch("HadPID",    &(K->hadPID),    "HadPID/I");
  T->Branch("Spin_idx",  &(K->tSpin),     "Spin_idx/I");
  T->Branch("SpinL_idx", &(K->lSpin),     "SpinL_idx/I");
  T->Branch("Weight",    &(weight),       "Weight/D");
};

// destructor
SimpleTree::~SimpleTree() {
};

