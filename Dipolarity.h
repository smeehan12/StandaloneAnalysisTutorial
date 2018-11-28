////////////////////////
/// 
/// Dipolarity.h
///  Originally by David Miller, modified by Sam Meehan
///  Nov 19, 2011
///
////////////////////////

#ifndef DIPOLARITY_H
#define DIPOLARITY_H

const Int_t kMaxEvent = 1;
const Int_t kMaxGenParticle = 994;

TTree          *fChain;

// Declaration of leaf types
Int_t           Event_;
UInt_t          Event_fUniqueID[kMaxEvent];   //[Event_]
UInt_t          Event_fBits[kMaxEvent];   //[Event_]
Long64_t        Event_Number[kMaxEvent];   //[Event_]
Int_t           Event_size;
Int_t           GenParticle_;
UInt_t          GenParticle_fUniqueID[kMaxGenParticle];   //[GenParticle_]
UInt_t          GenParticle_fBits[kMaxGenParticle];   //[GenParticle_]
Int_t           GenParticle_PID[kMaxGenParticle];   //[GenParticle_]
Int_t           GenParticle_Status[kMaxGenParticle];   //[GenParticle_]
Int_t           GenParticle_M1[kMaxGenParticle];   //[GenParticle_]
Int_t           GenParticle_M2[kMaxGenParticle];   //[GenParticle_]
Int_t           GenParticle_D1[kMaxGenParticle];   //[GenParticle_]
Int_t           GenParticle_D2[kMaxGenParticle];   //[GenParticle_]
Double_t        GenParticle_E[kMaxGenParticle];   //[GenParticle_]
Double_t        GenParticle_Px[kMaxGenParticle];   //[GenParticle_]
Double_t        GenParticle_Py[kMaxGenParticle];   //[GenParticle_]
Double_t        GenParticle_Pz[kMaxGenParticle];   //[GenParticle_]
Double_t        GenParticle_PT[kMaxGenParticle];   //[GenParticle_]
Double_t        GenParticle_Eta[kMaxGenParticle];   //[GenParticle_]
Double_t        GenParticle_Phi[kMaxGenParticle];   //[GenParticle_]
Double_t        GenParticle_Rapidity[kMaxGenParticle];   //[GenParticle_]
Double_t        GenParticle_T[kMaxGenParticle];   //[GenParticle_]
Double_t        GenParticle_X[kMaxGenParticle];   //[GenParticle_]
Double_t        GenParticle_Y[kMaxGenParticle];   //[GenParticle_]
Double_t        GenParticle_Z[kMaxGenParticle];   //[GenParticle_]
Int_t           GenParticle_size;

// List of branches
TBranch        *b_Event_;   //!
TBranch        *b_Event_fUniqueID;   //!
TBranch        *b_Event_fBits;   //!
TBranch        *b_Event_Number;   //!
TBranch        *b_Event_size;   //!
TBranch        *b_GenParticle_;   //!
TBranch        *b_GenParticle_fUniqueID;   //!
TBranch        *b_GenParticle_fBits;   //!
TBranch        *b_GenParticle_PID;   //!
TBranch        *b_GenParticle_Status;   //!
TBranch        *b_GenParticle_M1;   //!
TBranch        *b_GenParticle_M2;   //!
TBranch        *b_GenParticle_D1;   //!
TBranch        *b_GenParticle_D2;   //!
TBranch        *b_GenParticle_E;   //!
TBranch        *b_GenParticle_Px;   //!
TBranch        *b_GenParticle_Py;   //!
TBranch        *b_GenParticle_Pz;   //!
TBranch        *b_GenParticle_PT;   //!
TBranch        *b_GenParticle_Eta;   //!
TBranch        *b_GenParticle_Phi;   //!
TBranch        *b_GenParticle_Rapidity;   //!
TBranch        *b_GenParticle_T;   //!
TBranch        *b_GenParticle_X;   //!
TBranch        *b_GenParticle_Y;   //!
TBranch        *b_GenParticle_Z;   //!
TBranch        *b_GenParticle_size;   //!


///////////////////////////
/// Output Variables
///////////////////////////

TTree *outTree;

Int_t                   EventNumber;
//std::vector<Double_t>* temp_E;
//std::vector<Double_t>* temp_px;
//std::vector<Double_t>* temp_py;
//std::vector<Double_t>* temp_pz;



Int_t                   el_n;
std::vector<Double_t>*  el_E;
std::vector<Double_t>*  el_Et;
std::vector<Double_t>*  el_pt;
std::vector<Double_t>*  el_m;
std::vector<Double_t>*  el_eta;
std::vector<Double_t>*  el_phi;

Int_t                   mu_n;
std::vector<Double_t>*  mu_E;
std::vector<Double_t>*  mu_Et;
std::vector<Double_t>*  mu_pt;
std::vector<Double_t>*  mu_m;
std::vector<Double_t>*  mu_eta;
std::vector<Double_t>*  mu_phi;

Int_t                   ph_n;
std::vector<Double_t>*  ph_E;
std::vector<Double_t>*  ph_Et;
std::vector<Double_t>*  ph_pt;
std::vector<Double_t>*  ph_m;
std::vector<Double_t>*  ph_eta;
std::vector<Double_t>*  ph_phi;

Int_t                   jet_AntiKt7_n;
std::vector<Int_t>*                   jet_AntiKt7_content_n;
std::vector<std::vector<Double_t> >*  jet_AntiKt7_content_E;
std::vector<std::vector<Double_t> >*  jet_AntiKt7_content_px;
std::vector<std::vector<Double_t> >*  jet_AntiKt7_content_py;
std::vector<std::vector<Double_t> >*  jet_AntiKt7_content_pz;
std::vector<Double_t>*  jet_AntiKt7_E;
std::vector<Double_t>*  jet_AntiKt7_Et;
std::vector<Double_t>*  jet_AntiKt7_pt;
std::vector<Double_t>*  jet_AntiKt7_m;
std::vector<Double_t>*  jet_AntiKt7_eta;
std::vector<Double_t>*  jet_AntiKt7_phi;

Int_t                   jet_AntiKt10_n;
std::vector<Int_t>*                   jet_AntiKt10_content_n;
std::vector<std::vector<Double_t> >*  jet_AntiKt10_content_E;
std::vector<std::vector<Double_t> >*  jet_AntiKt10_content_px;
std::vector<std::vector<Double_t> >*  jet_AntiKt10_content_py;
std::vector<std::vector<Double_t> >*  jet_AntiKt10_content_pz;
std::vector<Double_t>*  jet_AntiKt10_E;
std::vector<Double_t>*  jet_AntiKt10_Et;
std::vector<Double_t>*  jet_AntiKt10_pt;
std::vector<Double_t>*  jet_AntiKt10_m;
std::vector<Double_t>*  jet_AntiKt10_eta;
std::vector<Double_t>*  jet_AntiKt10_phi;

Double_t  MET_etx;
Double_t  MET_ety;
Double_t  MET_phi;
Double_t  MET_et;

Double_t  MET_nu_etx;
Double_t  MET_nu_ety;
Double_t  MET_nu_phi;
Double_t  MET_nu_et;

// Counting for cuts
Int_t el_tot = 0;
Int_t el_acc = 0;
Int_t mu_tot = 0;
Int_t mu_acc = 0;
Int_t ph_tot = 0;
Int_t ph_acc = 0;
Int_t ak10jet_tot = 0;
Int_t ak10jet_acc = 0;
Int_t ak7jet_tot = 0;
Int_t ak7jet_acc = 0;

#endif
