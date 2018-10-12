#ifndef KnuclHist_h
#define KnuclHist_h 1

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TString.h"
#include "TVector3.h"
#include "TBranch.h"

#include "globals.hh"
#include "KnuclRootData.h"

class TFile;
class KnuclReactionData;
class KnuclHist
{
private:
    TFile *rootfile;

    G4int ReactionID;

    char *output_file;

  public:
    KnuclHist();
    ~KnuclHist();
    void BeginOfRunAction(const char *name);
    void EndOfRunAction();
    void BeginOfEventAction();
    void EndOfEventAction();

    G4int GetReactionID(){ return ReactionID;};
    void SetReactionID(G4int ID){ ReactionID=ID;};

  //----------//
  // for tree
  //----------//
private:
  TTree *tree;  // for eventHeaderMC/detectorData/mcData/reactionData
  TTree *tree2; // for runHeaderMC

  RunHeaderMC* runHeaderMC;
  EventHeaderMC* eventHeaderMC;
  DetectorData* detectorData;
  MCData* mcData;
  ReactionData* reactionData;

public:
  void setRunHeaderMC(RunHeaderMC* val){runHeaderMC=val;}
  void setEventHeaderMC(EventHeaderMC* val){eventHeaderMC=val;}
  void setDetectorData(DetectorData* val){detectorData=val;}
  void setMCData(MCData* val){mcData=val;}
  void setReactionData(ReactionData* val){reactionData=val;}
  void fillTree(){ tree->Fill(); }
  void showTree(){ tree->Show(); }
  void fillTree2(){ tree2->Fill(); }
  void showTree2(){ tree2->Show(); }
};

#endif
