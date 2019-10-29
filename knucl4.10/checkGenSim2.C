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

void piS(const char *filename="")
{
  gSystem->Load("./build/KnuclRootData_cc.so"); //<-- OR read from ./.rootlogon.C
  
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(111111);
  
  gROOT->cd();

  TDatabasePDG *pdg = new TDatabasePDG();
  pdg->ReadPDGTable("pdg_table.txt");
  
  //TFile *f = new TFile("ReacID_735.root");
  //TFile *f = new TFile("kpp_50_50_Lp.root");
  TFile *f = new TFile(filename);
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
  
  const int    BIN = 1500.0;
  const double MIN = 0.0;
  const double MAX = 1500.0;
  const double MAXF = 300.0;
  const double MINT = -1.0;
  const double MAXT = 1.0;
  char com[128];
  // for study
  TH1F* his0[2]; TH1F* his00[2]; TH1F* his01[2];
  TH1F* his1[6]; TH1F* his10[6]; TH1F* his11[6];
  TH1F* his2[6]; TH1F* his20[6]; TH1F* his21[6];
  TH1F* his3[6]; TH1F* his30[6]; TH1F* his31[6];
  TH1F* his4[6]; TH1F* his40[6]; TH1F* his41[6];
  for( int i=0; i<2; i++ ){
    sprintf(com, "his0-%d", i);
    his0[i] = new TH1F(com, "Fermi mom (Lab)", BIN, MIN, MAXF);
    sprintf(com, "his00-%d", i);
    his00[i] = new TH1F(com, "Fermi mom (Lab)", BIN, MIN, MAXF);
    sprintf(com, "his01-%d", i);
    his01[i] = new TH1F(com, "Fermi mom (Lab)", BIN, MIN, MAXF);
  }
  for( int i=0; i<6; i++ ){
    sprintf(com, "his1-%d", i);
    his1[i] = new TH1F(com, "paritlce mom (CM)", BIN, MIN, MAX);
    sprintf(com, "his10-%d", i);
    his10[i] = new TH1F(com, "paritlce mom (CM)", BIN, MIN, MAX);
    sprintf(com, "his11-%d", i);
    his11[i] = new TH1F(com, "paritlce mom (CM)", BIN, MIN, MAX);

    sprintf(com, "his2-%d", i);
    his2[i] = new TH1F(com, "paritlce mom (Lab)", BIN, MIN, MAX);
    sprintf(com, "his20-%d", i);
    his20[i] = new TH1F(com, "paritlce mom (Lab)", BIN, MIN, MAX);
    sprintf(com, "his21-%d", i);
    his21[i] = new TH1F(com, "paritlce mom (Lab)", BIN, MIN, MAX);

    sprintf(com, "his3-%d", i);
    his3[i] = new TH1F(com, "angler distribution (CM)", 50, MINT, MAXT);
    sprintf(com, "his30-%d", i);
    his30[i] = new TH1F(com, "angler distribution (CM)", 50, MINT, MAXT);
    sprintf(com, "his31-%d", i);
    his31[i] = new TH1F(com, "angler distribution (CM)", 50, MINT, MAXT);

    sprintf(com, "his4-%d", i);
    his4[i] = new TH1F(com, "angler distribution (Lab)", 50, MINT, MAXT);
    sprintf(com, "his40-%d", i);
    his40[i] = new TH1F(com, "angler distribution (Lab)", 50, MINT, MAXT);
    sprintf(com, "his41-%d", i);
    his41[i] = new TH1F(com, "angler distribution (Lab)", 50, MINT, MAXT);
  }
  TH1F* his99  = new TH1F("his99", "reactionID", 1120, 0, 1120);
  TH1F* his100 = new TH1F("his100", "reactionID", 1120, 0, 1120);
  TH1F* his101 = new TH1F("his101", "reactionID", 1120, 0, 1120);
  TH2F* hisd   = new TH2F("hisd", "Dalitz's plot", 44, -0.65, 0.65, 44, -0.05, 1.05);
  TH2F* his_ipi = new TH2F("his_ipi", "initial pi momentum vs. cos(#theta_{#pi}^{CM})", 80, MINT, MAXT, 150, 0, 1.5);
  TH2F* his_S   = new TH2F("his_S", "#Sigma momentum vs. cos(#theta_{#pi}^{CM})", 80, MINT, MAXT, 150, 0, 1.5);
  TH2F* his_pi  = new TH2F("his_pi", "pi momentum vs. cos(#theta_{#pi}^{CM})", 80, MINT, MAXT, 150, 0, 1.5);
  TH2F* his_pi_acc = new TH2F("his_pi_acc", "pi momentum vs. cos(#theta_{#pi}^{CM})", 80, MINT, MAXT, 150, 0, 1.5);
  TH2F* his_n   = new TH2F("his_n", "neutron momentum vs. cos(#theta_{#pi}^{CM})", 80, MINT, MAXT, 150, 0, 1.5);
  TH2F* his_n_acc   = new TH2F("his_n_acc", "neutron momentum vs. cos(#theta_{#pi}^{CM})", 80, MINT, MAXT, 150, 0, 1.5);
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

    //### neutron from Sigma ###//
    double cospi = reactionData->GetCMParticle(0).CosTheta();
    double mom_S  = -1;
    double mom_pi = -1;
    double mom_n  = -1;
    double theta_n  = -1;
    double theta_pi = -1;
    bool acc_flag = false;
    for (Int_t j=0; j<mcData->trackSize(); j++) {
      int track   = mcData->track(j)->trackID();
      int pdgcode = mcData->track(j)->pdgID();
      int parent  = mcData->track(j)->parentTrackID();
      //cerr<<track<<" "<<pdgcode<<" "<<parent<<endl;
      if( pdgcode==2112 && parent==3 ){ // neutron from Sigma
	mom_n = mcData->track(j)->momentum().Mag()/1000;
	theta_n = mcData->track(j)->momentum().Theta()*180/3.1415; // degrees
      }
      if( abs(pdgcode)==211 && parent==3 ){ // pi from Sigma
	mom_pi = mcData->track(j)->momentum().Mag()/1000;
	theta_pi = mcData->track(j)->momentum().Theta()*180/3.1415; // degrees
      }
    }
    acc_flag = ( 54<theta_n && theta_n<126 && 54<theta_pi && theta_pi<126 );
    his_ipi->Fill(cospi,reactionData->GetParticle(0).P()/1000);
    his_S->Fill(cospi,reactionData->GetParticle(1).P()/1000);
    his_pi->Fill(cospi,mom_pi);
    his_n->Fill(cospi,mom_n);
    if( acc_flag ){
      his_pi_acc->Fill(cospi,mom_pi);
      his_n_acc->Fill(cospi,mom_n);
    }

    //### all events ###//
    for (Int_t j=0; j<nspec; j++) {
       his0[j]->Fill(reactionData->FermiMom(j));
    }
    for (Int_t j=0; j<nparticle; j++) {
      his1[j]->Fill(reactionData->GetCMParticle(j).P());
      his2[j]->Fill(reactionData->GetParticle(j).P());
      his3[j]->Fill(reactionData->GetCMParticle(j).CosTheta());
      his4[j]->Fill(reactionData->GetParticle(j).CosTheta());
    }
    his99->Fill(reactionData->ReactionID());

    //### in CDH acceptance ###//    
    double weight = 1;
    if( acc_flag  ){
      //std::cerr<<CDHhit[0]<<":"<<part[0]<<" "<<CDHhit[1]<<":"<<part[1]<<std::endl;
      for (Int_t j=0; j<nparticle; j++) {
	his10[j]->Fill(reactionData->GetCMParticle(j).P(), weight);
	his20[j]->Fill(reactionData->GetParticle(j).P(), weight);
	his30[j]->Fill(reactionData->GetCMParticle(j).CosTheta(), weight);
	his40[j]->Fill(reactionData->GetParticle(j).CosTheta(), weight);
      }
      his100->Fill(reactionData->ReactionID());
    }

  }// for (Int_t i=0;i<nevent;i++) {
  cout<<"end of filling"<<endl;
  //----------------------//
  //--- event roop end ---//
  //----------------------//



  //--- plot ---//
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(0);

  TLatex *tex;

  //### plot reaction-ID
  TCanvas *c1 = new TCanvas("c1", "", 800, 400);
  c1->Divide(2,1);
  c1->cd(1); his99->Draw("");
  his99->SetXTitle("reactionID [CERN-HERA-83-02]");
  int max = his99->GetMaximum();
  tex = new TLatex(MINT+0.05*(MAXT-MINT), max*0.9, "All events");
  tex->SetTextColor(1);
  tex->SetTextSize(0.05);
  tex->Draw();

  c1->cd(2); his100->Draw("");
  his100->SetXTitle("reactionID [CERN-HERA-83-02]");
  his100->SetLineColor(2);
  int max = his100->GetMaximum();
  tex = new TLatex(MINT+0.05*(MAXT-MINT), max*0.9, "CDH-2charged");
  tex->SetTextColor(2);
  tex->SetTextSize(0.05);
  tex->Draw();
  c1->Print("tmp.pdf(");


  //### plot particle momentum @ CM
  TCanvas *c3 = new TCanvas("c3", "", 800, 400);
  c3->Divide(2,1);
  for( int i=0; i<nparticle; i++ ){
    c3->cd(i+1); his1[i]->Draw("");
    his1[i]->SetXTitle("momentum [MeV/c]");
    int max = his1[i]->GetMaximum();
    tex = new TLatex(MIN+0.05*(MAX-MIN), max*0.9, name[i].c_str());
    tex->SetTextColor(4);
    tex->SetTextSize(0.1);
    tex->Draw();
    his10[i]->Draw("same"); his10[i]->SetLineColor(2);
    his11[i]->Draw("same"); his11[i]->SetLineColor(4);
  }
  c3->Print("tmp.pdf");

  //### plot particle momentum @ LAB
  TCanvas *c4 = new TCanvas("c4", "", 800, 400);
  c4->Divide(2,1);
  for( int i=0; i<nparticle; i++ ){
    c4->cd(i+1); his2[i]->Draw(""); //gPad->SetLogy();
    his2[i]->SetXTitle("momentum [MeV/c]");
    int max = his2[i]->GetMaximum();
    tex = new TLatex(MIN+0.05*(MAX-MIN), max*0.9, name[i].c_str());
    tex->SetTextColor(4);
    tex->SetTextSize(0.1);
    tex->Draw();
    his20[i]->Draw("same"); his20[i]->SetLineColor(2);
    his21[i]->Draw("same"); his21[i]->SetLineColor(4);
  }
  c4->Print("tmp.pdf");

  //### plot angular distribution @ CM
  TCanvas *c6 = new TCanvas("c6", "", 800, 400);
  c6->Divide(2,1);
  for( int i=0; i<nparticle; i++ ){
    c6->cd(i+1); his3[i]->Draw("");
    his3[i]->SetStats(0);
    his3[i]->SetMinimum(0);
    his3[i]->SetXTitle("cos#theta");
    int max = his3[i]->GetMaximum();
    tex = new TLatex(MINT+0.05*(MAXT-MINT), max*0.9, name[i].c_str());
    tex->SetTextColor(4);
    tex->SetTextSize(0.1);
    tex->Draw();
    his30[i]->Draw("same"); his30[i]->SetLineColor(2);
    his31[i]->Draw("same"); his31[i]->SetLineColor(4);
  }
  c6->Print("tmp.pdf");
  //c6->Print("tmp6.png");

  //### plot angular distribution @ LAB
  TCanvas *c7 = new TCanvas("c7", "", 800, 400);
  c7->Divide(2,1);
  for( int i=0; i<nparticle; i++ ){
    c7->cd(i+1); his4[i]->Draw("");
    his4[i]->SetStats(0);
    his4[i]->SetMinimum(0);
    his4[i]->SetXTitle("cos#theta");
    int max = his4[i]->GetMaximum();
    tex = new TLatex(MINT+0.05*(MAXT-MINT), max*0.9, name[i].c_str());
    tex->SetTextColor(4);
    tex->SetTextSize(0.1);
    tex->Draw();
    his40[i]->Draw("same"); his40[i]->SetLineColor(2);
    his41[i]->Draw("same"); his41[i]->SetLineColor(4);
  }
  c7->Print("tmp.pdf");



  TCanvas *c8 = new TCanvas("c8", "", 800, 400);
  c8->Divide(2,1);
  c8->cd(1); his_pi->Draw("colz"); gPad->SetGrid();
  his_pi->SetXTitle("cos#theta_{#pi}^{CM}");
  his_pi->SetYTitle("#pi from #Sigma momentum [GeV/c]");
  c8->cd(2); his_n->Draw("colz"); gPad->SetGrid();
  his_n->SetXTitle("cos#theta_{#pi}^{CM}");
  his_n->SetYTitle("n from #Sigma momentum [GeV/c]");
  c8->Print("tmp.pdf");

  TCanvas *c9 = new TCanvas("c9", "", 800, 400);
  c9->Divide(2,1);
  c9->cd(1); his_pi_acc->Draw("colz"); gPad->SetGrid();
  his_pi_acc->SetXTitle("cos#theta_{#pi}^{CM}");
  his_pi_acc->SetYTitle("#pi from #Sigma momentum [GeV/c]");
  c9->cd(2); his_n_acc->Draw("colz"); gPad->SetGrid();
  his_n_acc->SetXTitle("cos#theta_{#pi}^{CM}");
  his_n_acc->SetYTitle("n from #Sigma momentum [GeV/c]");
  c9->Print("tmp.pdf");


  TCanvas *c10 = new TCanvas("c10", "", 1200, 400);
  c10->SetGrid();
  c10->Divide(3,1);
  int sta = 69;
  int bin = 4;
  for( int i=0; i<3; i++ ){
    c10->cd(i+1);
    double a = sta+i*bin;
    double b = sta+(i+1)*bin-1;
    double aa = a*0.025-1;
    double bb = b*0.025-1;
    sprintf(com, "%0.1f-%0.1f", aa, bb);
    his_n_acc->ProjectionY(com, a, b)->Draw();
  }
  c10->Print("tmp.pdf");


  TCanvas *c11 = new TCanvas("c11", "", 800, 400);
  c11->Divide(2,1);
  c11->cd(1); his_ipi->Draw("colz"); gPad->SetGrid();
  his_ipi->SetXTitle("cos#theta_{#pi}^{CM}");
  his_ipi->SetYTitle("initial #pi momentum [GeV/c]");
  c11->cd(2); his_S->Draw("colz"); gPad->SetGrid();
  his_S->SetXTitle("cos#theta_{#pi}^{CM}");
  his_S->SetYTitle("#Sigma momentum [GeV/c]");
  c11->Print("tmp.pdf)");

}
