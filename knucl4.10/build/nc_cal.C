// 2015, 5/31, F.Sakuma
// for neutron counter with CDH study
//
// Knucl:ProcessID = neutron_efficiency_mode
//
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

void nc_cal()
{
  gSystem->Load("KnuclRootData_cc.so"); //<-- OR read from ./.rootlogon.C
  
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(111111);
  
  gROOT->cd();

  TDatabasePDG *pdg = new TDatabasePDG();
  pdg->ReadPDGTable("pdg_table.txt");
  
  //TFile *f = new TFile("test_mc.root");
  //TFile *f = new TFile("nc_cal_scinti.root");
  //TFile *f = new TFile("test_singlen.root");
  TFile *f = new TFile("test_mod.root");
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

  // scale facor for energy deposit //
  double dE_par1 = 0.2;
  double dE_par2 = 40;
  bool dE_scale_flag = false;//true;

  double threshold = 1.0; // MeVee

  tree2->GetEvent(0);

  Int_t nevent = tree->GetEntries();
  //nevent = 1000;

  const int    BIN = 1000.0;
  const double MIN = 0.0;
  const double MAX = 100.0;
  char com[128];
  sprintf(com, "CDH_dE");
  TH1F* his0 = new TH1F(com, com, BIN, MIN, MAX);
  sprintf(com, "CDH_dE_cut");
  TH1F* his1 = new TH1F(com, com, BIN, MIN, MAX);
  TH1F* his2[36];
  TH1F* his3[36];
  for( int i=0; i<36; i++ ){
    sprintf(com, "his2-%d", i);
    his2[i] = new TH1F(com, "CDH_dE", BIN, MIN, MAX);
    sprintf(com, "his3-%d", i);
    his3[i] = new TH1F(com, "CDH_dE_cut", BIN, MIN, MAX);
  }
  sprintf(com, "dE_max_seg");
  TH1F* his4 = new TH1F(com, com, 36, -0.5, 10.5);

  double mom[15];
  double mom_err[15];
  int num_mom[15];
  int acc_mom[15];
  double eff_mom[15];
  double eff_mom_err[15];
  for( int j=0; j<15; j++ ){
    mom[j] = j*0.1+0.05;
    mom_err[j] = 0.05;
    num_mom[j] = 0;
    acc_mom[j] = 0;
    eff_mom[j] = 0;
    eff_mom_err[j] = 0;
  }
  TH1F* his20[15];
  TH1F* his30[15];
  for( int i=0; i<15; i++ ){
    sprintf(com, "his20-%d", i);
    his20[i] = new TH1F(com, "CDH_dE", BIN, MIN, MAX);
    sprintf(com, "his30-%d", i);
    his30[i] = new TH1F(com, "CDH_dE_cut", BIN, MIN, MAX);
  }

  //------------------------//
  //--- event roop start ---//
  //------------------------//
  int count0 = 0;
  int count1 = 0;
  cout<<"Start to fill histgrams. Entries = "<<nevent<<endl;
  for (Int_t i=0; i<nevent; i++) {
    tree->GetEvent(i);
    //int ndecay    = reactionData->NParticle(0);
    //int nspec     = reactionData->NParticle(1);
    //int nparticle = ndecay+nspec;
    //std::cout << "ndecay " << ndecay << std::endl;
    //std::cout << "nspec " << nspec << std::endl;
    //if(nspec==0) continue;
    double momentum = reactionData->GetParticle(0).P();
    int m = int(momentum/100);
    num_mom[m]++;

    double ene_dep = 0;
    double de_seg[36] = {0};
    for (Int_t j=0; j<detectorData->detectorHitSize(); j++) {
      int cid = detectorData->detectorHit(j)->detectorID();
      if( cid == CID_CDH ){
	int seg = detectorData->detectorHit(j)->channelID();
	double de = detectorData->detectorHit(j)->adc();
	ene_dep += de;
	de_seg[seg] = de;
	//cerr<<i<<" seg = "<<seg<<" dE = "<<de<<endl;
      }
    }
    
    //if(ene_dep) cerr<<ene_dep<<" -> ";
    //if( dE_scale_flag ) ene_dep = ene_dep*dE_par1*(exp(ene_dep/dE_par2));
    //if(ene_dep) cerr<<ene_dep<<endl;

    if( 0 < ene_dep ){
      his0->Fill(ene_dep);
      for( int j=0; j<16; j++ ) his2[j]->Fill(de_seg[j]);
      his20[m]->Fill(ene_dep);
      count0++;
      if( threshold < ene_dep ){
	his1->Fill(ene_dep);
	for( int j=0; j<16; j++ ) his3[j]->Fill(de_seg[j]);
	his30[m]->Fill(ene_dep);
	count1++;
	acc_mom[m]++;
      
	double ene_max = 0;
	double max_seg = -1;
	for (Int_t j=0; j<16; j++) {
	  if( ene_max < de_seg[j] ){
	    max_seg = j;
	    ene_max = de_seg[j];
	  }
	}
	his4->Fill(max_seg);
      } // if( threshold < ene_dep ){
    } // if( 0 < ene_dep ){

  }// for (Int_t i=0;i<nevent;i++) {
  cout<<"end of filling"<<endl;
  //----------------------//
  //--- event roop end ---//
  //----------------------//
  cerr<<"# of hit = "<<count0<<" / " <<count1<<endl;

  for( int j=1; j<15; j++ ){
    eff_mom[j] = (double)acc_mom[j]/num_mom[j];
    eff_mom_err[j] = sqrt(eff_mom[j]*(1-eff_mom[j])/num_mom[j]);
    eff_mom[j] *= 100; // percent
    eff_mom_err[j] *= 100; // percent
  }
  //cerr<<mom[j]<<" "<<eff_mom[j]<<" +/- "<<eff_mom_err[j]<<endl;
  for( int j=1; j<15; j++ ) cerr<<eff_mom[j]<<", ";
  cerr<<endl;
  for( int j=1; j<15; j++ ) cerr<<eff_mom_err[j]<<", ";
  cerr<<endl;

  //--- plot ---//
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(0);

  TCanvas *c1 = new TCanvas("c1", "", 1800, 600);
  c1->Divide(3,1);
  c1->cd(1); his0->Draw("");
  his0->SetXTitle("dE (MeV)");
  c1->cd(2); his1->Draw("");
  his1->SetXTitle("dE (MeV)");
  c1->cd(3); his4->Draw("");
  his4->SetXTitle("dE maximum segment");
  c1->Print("tmp.pdf(");

#if 1
  TCanvas *c3 = new TCanvas("c3", "", 1200, 1200);
  c3->Divide(4,4);
  for( int i=0; i<16; i++ ){
    c3->cd(i+1); his2[i]->Draw("");
    his2[i]->SetXTitle("dE (MeV)");
  }
  c3->Print("tmp.pdf");

  TCanvas *c4 = new TCanvas("c4", "", 1200, 1200);
  c4->Divide(4,4);
  for( int i=0; i<16; i++ ){
    c4->cd(i+1); his3[i]->Draw("");
    his3[i]->SetXTitle("dE (MeV)");
  }
  c4->Print("tmp.pdf");
#endif

#if 1
  TCanvas *c5 = new TCanvas("c5", "", 1200, 1200);
  c5->Divide(4,4);
  for( int i=0; i<15; i++ ){
    c5->cd(i+1); his20[i]->Draw("");
    his20[i]->SetXTitle("dE (MeV)");
  }
  c5->Print("tmp.pdf");

  TCanvas *c6 = new TCanvas("c6", "", 1200, 1200);
  c6->Divide(4,4);
  for( int i=0; i<15; i++ ){
    c6->cd(i+1); his30[i]->Draw("");
    his30[i]->SetXTitle("dE (MeV)");
  }
  c6->Print("tmp.pdf");
#endif

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TCanvas *c2 = new TCanvas("c2", "", 600, 600);
  c2->SetGrid();
  TH2F *frame = new TH2F("frame", "", 150, 0, 1.5, 100, 0, 30);
  frame->Draw();
  frame->SetXTitle("momentum [GeV/c]");
  frame->SetYTitle("efficiency [%]");
  TGraph *gr;
  gr = new TGraphErrors(15, mom, eff_mom, mom_err, eff_mom_err);
  gr->SetMarkerStyle(20); gr->SetMarkerColor(4);
  gr->SetMarkerSize(1); gr->SetLineColor(4);
  gr->Draw("P");
  c2->Print("tmp.pdf)");

}
