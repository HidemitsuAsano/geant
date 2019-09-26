#include "KnuclHist.hh"
#include "globals.hh"

KnuclHist::KnuclHist()
{
}

KnuclHist::~KnuclHist()
{
}

void KnuclHist::BeginOfRunAction(const char *filename)
{
  rootfile=new TFile(filename,"RECREATE");
  //rootfile->SetCompressionLevel(9);

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
  // tree
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
  tree  = new TTree("tree", "knucl: event data");
  tree2 = new TTree("tree2","knucl: run info.");
  
  runHeaderMC    = new RunHeaderMC();
  eventHeaderMC  = new EventHeaderMC();
  detectorData = new DetectorData();
  mcData       = new MCData();
  reactionData = new ReactionData();
  
  tree2->Branch("RunHeaderMC","RunHeaderMC",&runHeaderMC);
  tree->Branch("EventHeaderMC","EventHeaderMC",&eventHeaderMC);
  tree->Branch("DetectorData","DetectorData",&detectorData);
  tree->Branch("MCData","MCData",&mcData);
  tree->Branch("ReactionData", "ReactionData", &reactionData);

}


void KnuclHist::EndOfRunAction()
{
  tree->Write();
  tree2->Write();
  rootfile->Write();
  rootfile->Close();
}

void KnuclHist::BeginOfEventAction()
{
}

void KnuclHist::EndOfEventAction()
{
}

