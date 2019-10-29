#include <fstream.h>
#include <iostream.h>

#define DALITZ 0

class RunHeader;
class EventHeader;
class DetectorData;
class MCData;
class ReactionData;

class DetectorHit;
class Track;

class CrossSecion;
class CrossSecionTable;

const double Const=29.97;// [cm/ns]

// Counter ID (from knucl4/include/KnuclCommon.h)
enum gCounterID { CID_CDC     = 0,
                  CID_CDH     = 1,
                  CID_BHD     = 2,
                  CID_PA      = 3,
                  CID_T0      = 4,
                  CID_E0      = 5,
                  CID_DEF     = 5,
                  CID_B1      = 6,
                  CID_LC1     = 7,
                  CID_LC2     = 8,
                  CID_AC      = 9,
                  CID_WC      = 10,
                  CID_GC      = 11,
                  CID_Range   = 12,
                  CID_B2      = 13,
                  CID_TOFstop = 14,
                  CID_CVC     = 14,
                  CID_BLC1a   = 15,
                  CID_BLC1b   = 16,
                  CID_BLC2a   = 17,
                  CID_BLC2b   = 18,
                  CID_SDD     = 19,
                  CID_BLC1    = 21,
                  CID_BLC2    = 22,
                  CID_FDC1    = 23,
                  CID_FDC2    = 24,
                  CID_ZVC     = 30,
                  CID_KDV     = 31,
                  CID_NC      = 32,
                  CID_BVC     = 33,
                  CID_PC      = 35,
                  CID_Longbar = 36,
                  CID_LB      = 36,
                  CID_WVC     = 37,
                  CID_BPC     = 40,
                  CID_BPD     = 41,
                  CID_IH      = 42,
                  CID_T0pre   = 51,
                  CID_T0post  = 52,
                  CID_BHDpost = 56,
                  CID_HVC1    = 61,
                  CID_HVC2    = 62,
                  CID_BHDmul  = 81,
                  CID_T0mul   = 82,
                  CID_BVCmul  = 82,
                  CID_HVC1mul = 83,
                  CID_HVC2mul = 84,
                  CID_REFmul  = 85,
                  CID_BD      = 90,
                  CID_BeamDump= 90,
                  CID_TEMP1   = 91,
                  CID_TEMP2   = 92,
                  CID_GPIO    = 97,
                  CID_MISC    = 98,
                  CID_TEMP    = 99,
                  CID_LS      = 150
};

void piLpn()
{
  gSystem->Load("./build/KnuclRootData_cc.so"); //<-- OR read from ./.rootlogon.C
  
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(111111);
  
  gROOT->cd();

  TDatabasePDG *pdg = new TDatabasePDG();
  pdg->ReadPDGTable("pdg_table.txt");
  
  TFile *f = new TFile("simdata015/sim_0000.root");
  TTree *tree = (TTree*)f->Get("tree");
  TTree *tree2 = (TTree*)f->Get("tree2");

  RunHeaderMC*    runHeaderMC=0;
  EventHeaderMC*  eventHeaderMC=0;
  DetectorData*   detectorData=0;
  MCData*         mcData=0;
  ReactionData*   reactionData=0;

  tree2->SetBranchAddress( "RunHeaderMC", &runHeaderMC );
  tree->SetBranchAddress( "EventHeaderMC", &eventHeaderMC );
  tree->SetBranchAddress( "DetectorData", &detectorData );
  tree->SetBranchAddress( "MCData", &mcData );
  tree->SetBranchAddress( "ReactionData", &reactionData );

  tree2->GetEvent(0);
  double totalCS  = runHeaderMC->CStable().TotalCS();
  int totalEv     = 0;
  int generatedEv = 0;
  for( int i=0; i<tree2->GetEntries(); i++ ){
    tree2->GetEvent(i);
    totalEv     += runHeaderMC->numEvent();
    generatedEv += runHeaderMC->numGenerated();
  }
  std::cerr<<"********************************"<<std::endl;
  std::cerr<<" seed                  = "<<runHeaderMC->seed()<<std::endl;
  std::cerr<<" # of total events     = "<<totalEv<<std::endl;
  std::cerr<<" # of events in tree   = "<<tree->GetEntries()<<std::endl;
  std::cerr<<"--------------------------------"<<std::endl;
  std::cerr<<" # of generated events = "<<generatedEv<<std::endl;
  std::cerr<<" total CrossSection    = "<<totalCS<<" mb"<<std::endl;
  std::cerr<<"********************************"<<std::endl;
  //runHeaderMC->CStable().PrintAllCS();


  Int_t nevent = tree->GetEntries();
  //nevent = 10000;

  int nevent_piLpn = 0;
  
  int FermiFlag = 0;
  string initname[2];
  string name[6];
  //------------------------//
  //--- event roop start ---//
  //------------------------//
  cout<<"Start to fill histgrams. Entries = "<<nevent<<endl;
  for (Int_t i=0; i<nevent; i++) {
    tree->GetEvent(i);
    // print information
    int ndecay    = reactionData->NParticle(0);
    int nspec     = reactionData->NParticle(1);
    int nparticle = ndecay+nspec;
    if( i==0 ){
      for (Int_t j=0; j<2; j++) {
	initname[j] = pdg->GetParticle(reactionData->InitPDG(j))->GetName(); 
      }
      for (Int_t j=0; j<nparticle; j++) {
	name[j] = pdg->GetParticle(reactionData->PDG(j))->GetName(); 
      }
      FermiFlag = ( reactionData->NParticle(1) ) ? 1 : 0;
      cerr<<" ##############################################"<<endl;
      cerr<<" ### Reaction ID = "<<reactionData->ReactionID()<<endl;
      cerr<<" ### Initial  Particles = ";
      for (Int_t j=0; j<2; j++) {
	cerr<<initname[j]<<"("<<reactionData->InitPDG(j)<<") ";
      } cerr<<endl;
      cerr<<" ### num of particles   = ("<<reactionData->NParticle(0)<<","
	  <<reactionData->NParticle(1)<<") = (#decays, #spectators)"<<endl;
      cerr<<" ### Final Particles    = ";
      for (Int_t j=0; j<nparticle; j++) {
	cerr<<name[j]<<"("<<reactionData->PDG(j)<<") ";
      } cerr<<endl;
      cerr<<" ### FermiMotion = "<<FermiFlag<<endl;;
      cerr<<" ##############################################"<<endl;
    }

    //### acceptance ###//
    int L_parent = 0;
    //     [4]      = {pi- from L, p from L, p ,n}
    int PDG[4]      = {-211, 2212, 2212, 2112};
    int parentID[4] = {2, 2, 0, 0};
    int ID[4]       = {-1, -1, -1, -1};
    int trackID[4]  = {-1, -1, -1, -1};
    int nparticle   = 0;
    for( int j=0; j<mcData->trackSize(); j++ ){
      int pdgcode = mcData->track(j)->pdgID();
      int parent  = mcData->track(j)->parentTrackID();
      int track   = mcData->track(j)->trackID();
      for( int k=0; k<4; k++ ){
        if( pdgcode==PDG[k] && parent==parentID[k] && ID[k]==-1 ){
          ID[k] = j;
          trackID[k] = track;
          nparticle++;
	}
      }
    }
    int nCDHhit[4]  = {0, 0, 0, 0};    
    for (Int_t j=0; j<detectorData->detectorHitSize(); j++) {
      int cid    = detectorData->detectorHit(j)->detectorID();
      int track  = detectorData->detectorHit(j)->trackID();
      int parent = mcData->track(j)->parentTrackID();
      for( int k=0; k<4; k++ ){
	if( cid==CID_CDH && track==trackID[k] ) nCDHhit[k]++;
      }
      //*** neutron hit search ***//
      if( cid==CID_CDH && parent==trackID[3] ) nCDHhit[3]++;
      //*** neutron hit search ***//
    }
    int piLpn_event = 1;
    for( int k=0; k<4; k++ ){
      piLpn_event *= nCDHhit[k];
    }

    if( piLpn_event ) nevent_piLpn++;

#if 0
    std::cerr<<piLpn_event<<" | ";
    for( int k=0; k<4; k++ ){
      std::cerr<<nCDHhit[k]<<" ";
    }
    std::cerr<<endl;
#endif
    
  }// for (Int_t i=0;i<nevent;i++) {
  cout<<"end of filling"<<endl;
  //----------------------//
  //--- event roop end ---//
  //----------------------//

  std::cerr<<" nevent_piLpn = "<<nevent_piLpn<<std::endl;
}
