#include <fstream>
#include <iostream>
#include <vector>
#include <string>

#include <TApplication.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include <TSystem.h>
#include <TBranch.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TRint.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TGraphErrors.h> 
#include <TDatabasePDG.h>
#include <TRandom.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TColor.h>

#include "/home/sakuma/work/ana/geant/knucl4.10/include/KnuclRootData.h"

class RunHeader;
class EventHeader;
class RunHeaderMC;
class EventHeaderMC;
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

int main()
{
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(111111);
  
  gROOT->cd();

  TDatabasePDG *pdg = new TDatabasePDG();
  pdg->ReadPDGTable("pdg_table.txt");
  
  double Kpp_mass   = 2*pdg->GetParticle("proton")->Mass()+pdg->GetParticle("K-")->Mass();
  double piSp_mass1 = pdg->GetParticle("pi+")->Mass()+pdg->GetParticle("Sigma-")->Mass()+pdg->GetParticle("proton")->Mass();
  double piSp_mass2 = pdg->GetParticle("pi-")->Mass()+pdg->GetParticle("Sigma+")->Mass()+pdg->GetParticle("proton")->Mass();
  double piLn_mass  = pdg->GetParticle("pi+")->Mass()+pdg->GetParticle("Lambda0")->Mass()+pdg->GetParticle("neutron")->Mass();
  double Kp_mass    = pdg->GetParticle("proton")->Mass()+pdg->GetParticle("K-")->Mass();
  double piS_mass1 = pdg->GetParticle("pi+")->Mass()+pdg->GetParticle("Sigma-")->Mass();
  double piS_mass2 = pdg->GetParticle("pi-")->Mass()+pdg->GetParticle("Sigma+")->Mass();
  double piL_mass  = pdg->GetParticle("pi+")->Mass()+pdg->GetParticle("Lambda0")->Mass();

#if 0
  const int    PISIGMA_BIN = 50;
  const double PISIGMA_MIN = 1.0;
  const double PISIGMA_MAX = 2.0;
#else
  const int    PISIGMA_BIN = 66;
  const double PISIGMA_MIN = 1.0;
  const double PISIGMA_MAX = 1.99;
#endif

  const int    PISIGMAP_BIN = 50;
  const double PISIGMAP_MIN = 2.0;
  const double PISIGMAP_MAX = 3.0;

  const int    COSN_BIN = 40;
  const double COSN_MIN = -1.0;
  const double COSN_MAX = 1.0;

  TH2F *Cosn_IMnpipi = new TH2F( "Cosn_IMnpipi", "Cosn_IMnpipi", PISIGMA_BIN, PISIGMA_MIN, PISIGMA_MAX, COSN_BIN, COSN_MIN, COSN_MAX );
  TH2F *Cosp_IMnpipi = new TH2F( "Cosp_IMnpipi", "Cosp_IMnpipi", PISIGMA_BIN, PISIGMA_MIN, PISIGMA_MAX, COSN_BIN, COSN_MIN, COSN_MAX );
  TH2F *CosX_IMnpipi = new TH2F( "CosX_IMnpipi", "CosX_IMnpipi", PISIGMA_BIN, PISIGMA_MIN, PISIGMA_MAX, COSN_BIN, COSN_MIN, COSN_MAX );

  TH2F *Cosn_IMnppipi = new TH2F( "Cosn_IMnppipi", "Cosn_IMnppipi", PISIGMAP_BIN, PISIGMAP_MIN, PISIGMAP_MAX, COSN_BIN, COSN_MIN, COSN_MAX );

  TH2F *Momp_IMnpipi = new TH2F( "Momp_IMnpipi", "Momp_IMnpipi", PISIGMA_BIN, PISIGMA_MIN, PISIGMA_MAX, 50, 0, 1.5 );
  TH2F *Momn_IMnpipi = new TH2F( "Momn_IMnpipi", "Momn_IMnpipi", PISIGMA_BIN, PISIGMA_MIN, PISIGMA_MAX, 50, 0, 1.5 );
  TH2F *MomX_IMnpipi = new TH2F( "MomX_IMnpipi", "MomX_IMnpipi", PISIGMA_BIN, PISIGMA_MIN, PISIGMA_MAX, 50, 0, 1.5 );

  TH2F *Cosn_IMnpipi_acc = new TH2F( "Cosn_IMnpipi_acc", "Cosn_IMnpipi_acc", PISIGMA_BIN, PISIGMA_MIN, PISIGMA_MAX, COSN_BIN, COSN_MIN, COSN_MAX );
  TH2F *Cosn_IMnppipi_acc = new TH2F( "Cosn_IMnppipi_acc", "Cosn_IMnppipi_acc", PISIGMAP_BIN, PISIGMAP_MIN, PISIGMAP_MAX, COSN_BIN, COSN_MIN, COSN_MAX );

  TH2F *IMnpipi_IMnppipi = new TH2F( "IMnpipi_IMnppipi", "IMnpipi_IMnppipi", PISIGMAP_BIN, PISIGMAP_MIN, PISIGMAP_MAX, PISIGMA_BIN, PISIGMA_MIN, PISIGMA_MAX );

  TH2F *Cosn_IMnpipi_Sp = new TH2F( "Cosn_IMnpipi_Sp", "Cosn_IMnpipi_Sp", PISIGMA_BIN, PISIGMA_MIN, PISIGMA_MAX, COSN_BIN, COSN_MIN, COSN_MAX );
  TH2F *Cosn_IMnpipi_Sm = new TH2F( "Cosn_IMnpipi_Sm", "Cosn_IMnpipi_Sm", PISIGMA_BIN, PISIGMA_MIN, PISIGMA_MAX, COSN_BIN, COSN_MIN, COSN_MAX );

  TH2F *Cosn_IMnppipi_Sp = new TH2F( "Cosn_IMnppipi_Sp", "Cosn_IMnppipi_Sp", PISIGMAP_BIN, PISIGMAP_MIN, PISIGMAP_MAX, COSN_BIN, COSN_MIN, COSN_MAX );
  TH2F *Cosn_IMnppipi_Sm = new TH2F( "Cosn_IMnppipi_Sm", "Cosn_IMnppipi_Sm", PISIGMAP_BIN, PISIGMAP_MIN, PISIGMAP_MAX, COSN_BIN, COSN_MIN, COSN_MAX );

  int FermiFlag = 0;
  std::string initname[2];
  std::string name[6];
  char file[64];


  for(int x=0; x<40; x++){

    //sprintf(file, "data17/sim_000%d.root", x);
    if( x<10)       sprintf(file, "simdata001/sim_000%d.root", x);
    else if( x<100) sprintf(file, "simdata001/sim_00%d.root", x);
  //TFile *f = new TFile("data9/sim_0001.root");
  TFile *f = new TFile(file);
  TTree *tree = (TTree*)f->Get("tree");
  TTree *tree2 = (TTree*)f->Get("tree2");

  RunHeaderMC*    runHeaderMC = 0;
  EventHeaderMC*  eventHeaderMC = 0;
  DetectorData*   detectorData = 0;
  MCData*         mcData = 0;
  ReactionData*   reactionData = 0;

  tree2->SetBranchAddress( "RunHeaderMC", &runHeaderMC );
  tree->SetBranchAddress( "EventHeaderMC", &eventHeaderMC );
  tree->SetBranchAddress( "DetectorData", &detectorData );
  tree->SetBranchAddress( "MCData", &mcData );
  tree->SetBranchAddress( "ReactionData", &reactionData );

  Int_t nevent = tree->GetEntries();
  //nevent=100;

  //------------------------//
  //--- event roop start ---//
  //------------------------//
  std::cout<<"Start to fill histgrams. Entries = "<<nevent<<std::endl;
  for (Int_t i=0; i<nevent; i++) {
    tree->GetEvent(i);

    if( i%10000==0 ) std::cout<<" !!! event number "<<i<<std::endl;

    //** print event information **//
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
      std::cerr<<" ##############################################"<<std::endl;
      std::cerr<<" ### Reaction ID = "<<reactionData->ReactionID()<<std::endl;
      std::cerr<<" ### Initial  Particles = ";
      for (Int_t j=0; j<2; j++) {
	std::cerr<<initname[j]<<"("<<reactionData->InitPDG(j)<<") ";
      } std::cerr<<std::endl;
      std::cerr<<" ### num of particles   = ("<<reactionData->NParticle(0)<<","
	  <<reactionData->NParticle(1)<<") = (#decays, #spectators)"<<std::endl;
      std::cerr<<" ### Final Particles    = ";
      for (Int_t j=0; j<nparticle; j++) {
	std::cerr<<name[j]<<"("<<reactionData->PDG(j)<<") ";
      } std::cerr<<std::endl;
      std::cerr<<" ### FermiMotion = "<<FermiFlag<<std::endl;;
      std::cerr<<" ##############################################"<<std::endl;
    }

    TLorentzVector mom_beam   = reactionData->GetInitParticle(0);
    TLorentzVector mom_target = reactionData->GetInitParticle(1);
    TLorentzVector mom_target_beam = mom_target+mom_beam;
    TVector3 boost = mom_target_beam.BoostVector();
    TLorentzVector mom_beam_CM = mom_beam;
    mom_beam_CM.Boost(-boost);

    TLorentzVector mom_pi_CM    = reactionData->GetCMParticle(3);
    TLorentzVector mom_Sigma_CM = reactionData->GetCMParticle(0);
    TLorentzVector mom_p_CM     = reactionData->GetCMParticle(1);
    TLorentzVector mom_n_CM     = reactionData->GetCMParticle(2);
    TLorentzVector mom_X_CM     = mom_pi_CM+mom_Sigma_CM;

    double cos_n = mom_n_CM.Vect().Dot(mom_beam_CM.Vect())/(mom_n_CM.Vect().Mag()*mom_beam_CM.Vect().Mag());
    double cos_p = mom_p_CM.Vect().Dot(mom_beam_CM.Vect())/(mom_p_CM.Vect().Mag()*mom_beam_CM.Vect().Mag());
    double cos_X = mom_X_CM.Vect().Dot(mom_beam_CM.Vect())/(mom_X_CM.Vect().Mag()*mom_beam_CM.Vect().Mag());

    TLorentzVector mom_piS  = mom_pi_CM+mom_Sigma_CM;
    TLorentzVector mom_piSp = mom_pi_CM+mom_Sigma_CM+mom_p_CM;

    Cosn_IMnpipi->Fill(mom_piS.M()/1000, cos_n);
    Cosp_IMnpipi->Fill(mom_piS.M()/1000, cos_p);
    CosX_IMnpipi->Fill(mom_piS.M()/1000, cos_X);

    Momn_IMnpipi->Fill(mom_piS.M()/1000, mom_n_CM.P()/1000);
    Momp_IMnpipi->Fill(mom_piS.M()/1000, mom_p_CM.P()/1000);
    MomX_IMnpipi->Fill(mom_piS.M()/1000, mom_X_CM.P()/1000);

    Cosn_IMnppipi->Fill(mom_piSp.M()/1000, cos_n);

    if( reactionData->ReactionID()==2120 ){ // 2120: K- 3He -> S+ p n p-
      Cosn_IMnpipi_Sp->Fill(mom_piS.M()/1000, cos_n);
      Cosn_IMnppipi_Sp->Fill(mom_piSp.M()/1000, cos_n);
    }
    else if( reactionData->ReactionID()==2130 ){ // 130: K- 3He -> S- p n p+
      Cosn_IMnpipi_Sm->Fill(mom_piS.M()/1000, cos_n);
      Cosn_IMnppipi_Sm->Fill(mom_piSp.M()/1000, cos_n);
    }


    //std::cerr<<"==============="<<std::endl;
    int init_pi = 5;
    int init_n  = 4;
    int init_p  = 3;
    int init_S  = 2;
    int pi_S    = -1;
    int n_S    = -1;
    int trackID[5] = {-1, -1, -1, -1, -1};
    for (Int_t j=0; j<mcData->trackSize(); j++) {
      int pdgcode = mcData->track(j)->pdgID();
      int parent  = mcData->track(j)->parentTrackID();
      int track   = mcData->track(j)->trackID();
      //if( parent==0 ) std::cerr<<" *** "<<pdgcode<<" "<<parent<<" "<<track<<std::endl;
      if( parent==init_S ){
	if( abs(pdgcode)==211 )	pi_S = track;
	if( pdgcode==2112 )	n_S  = track;
      }
    }
    trackID[0] = init_n;
    trackID[1] = init_p;
    trackID[2] = init_pi;
    trackID[3] = pi_S;
    trackID[4] = n_S;
    
    int nCDH[5] = {0, 0, 0, 0, 0};
    for (Int_t j=0; j<detectorData->detectorHitSize(); j++) {
      int cid     = detectorData->detectorHit(j)->detectorID();
      int pdgcode = detectorData->detectorHit(j)->pdg();
      int parent  = detectorData->detectorHit(j)->parentID();
      int track   = detectorData->detectorHit(j)->trackID();
      if( cid==CID_CDH ){
	for(int k=0; k<5; k++ ){
	  if( track==trackID[k] ) nCDH[k]++;
	}
      }
    }
    IMnpipi_IMnppipi->Fill(mom_piSp.M()/1000, mom_piS.M()/1000);
    //for(int k=0; k<5; k++ ){ std::cerr<<nCDH[k]<<" "; }std::cerr<<std::endl;
    if( nCDH[1]*nCDH[2]*nCDH[3]*nCDH[4] ){
      //std::cerr<<" !!! piSp detected !!! "<<std::endl;
      Cosn_IMnpipi_acc->Fill(mom_piS.M()/1000, cos_n);
      Cosn_IMnppipi_acc->Fill(mom_piSp.M()/1000, cos_n);
      //IMnpipi_IMnppipi->Fill(mom_piSp.M()/1000, mom_piS.M()/1000);
    }


  }// for (Int_t i=0;i<nevent;i++) {
  std::cout<<"end of filling"<<std::endl;
  //----------------------//
  //--- event roop end ---//
  //----------------------//

  }

  //--- plot ---//
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(0);

  TLegend *leg;
  TLatex *tex;
  TLine *line;
  TH1F *his;
  double ymax, width;

  TCanvas *c1 = new TCanvas("c1", "", 600, 600);
  gPad->SetLogz();
  Cosn_IMnpipi->Draw("colz");
  Cosn_IMnpipi->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  Cosn_IMnpipi->SetYTitle("cos#theta_{n}^{CM}");
  Cosn_IMnpipi->SetStats(0);
  line = new TLine(Kp_mass, -1, Kp_mass, 1);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->Draw();
  line = new TLine(piS_mass1, -1, piS_mass1, 1);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->Draw();
  line = new TLine(piS_mass2, -1, piS_mass2, 1);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->Draw();
  c1->Print("tmp.pdf(");
  
  TCanvas *c2;
  c2 = new TCanvas("c2", "", 600, 600);
  gPad->SetLogz();
  Cosp_IMnpipi->Draw("colz");
  Cosp_IMnpipi->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  Cosp_IMnpipi->SetYTitle("cos#theta_{p}^{CM}");
  Cosp_IMnpipi->SetStats(0);
  line = new TLine(Kp_mass, -1, Kp_mass, 1);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->Draw();
  line = new TLine(piS_mass1, -1, piS_mass1, 1);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->Draw();
  line = new TLine(piS_mass2, -1, piS_mass2, 1);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->Draw();
  c2->Print("tmp.pdf");

  TCanvas *c3;
  c3 = new TCanvas("c3", "", 600, 600);
  gPad->SetLogz();
  CosX_IMnpipi->Draw("colz");
  CosX_IMnpipi->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  CosX_IMnpipi->SetYTitle("cos#theta_{n#pi^{+}#pi^{-}}^{CM}");
  CosX_IMnpipi->SetStats(0);
  line = new TLine(Kp_mass, -1, Kp_mass, 1);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->Draw();
  line = new TLine(piS_mass1, -1, piS_mass1, 1);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->Draw();
  line = new TLine(piS_mass2, -1, piS_mass2, 1);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->Draw();
  c3->Print("tmp.pdf");

  TCanvas *c4;
  c4 = new TCanvas("c4", "", 600, 600);
  gPad->SetLogz();
  Momp_IMnpipi->Draw("colz");
  Momp_IMnpipi->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  Momp_IMnpipi->SetYTitle("p_{p}^{CM} [GeV/c]");
  Momp_IMnpipi->SetStats(0);
  line = new TLine(Kp_mass, 0, Kp_mass, 1.5);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->Draw();
  line = new TLine(piS_mass1, 0, piS_mass1, 1.5);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->Draw();
  line = new TLine(piS_mass2, 0, piS_mass2, 1.5);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->Draw();
  c4->Print("tmp.pdf");

  TCanvas *c5;
  c5 = new TCanvas("c5", "", 600, 600);
  gPad->SetLogz();
  Momn_IMnpipi->Draw("colz");
  Momn_IMnpipi->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  Momn_IMnpipi->SetYTitle("p_{missing-n}^{CM} [GeV/c]");
  Momn_IMnpipi->SetStats(0);
  line = new TLine(Kp_mass, 0, Kp_mass, 1.5);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->Draw();
  line = new TLine(piS_mass1, 0, piS_mass1, 1.5);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->Draw();
  line = new TLine(piS_mass2, 0, piS_mass2, 1.5);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->Draw();
  c5->Print("tmp.pdf");

  TCanvas *c6;
  c6 = new TCanvas("c6", "", 600, 600);
  gPad->SetLogz();
  MomX_IMnpipi->Draw("colz");
  MomX_IMnpipi->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  MomX_IMnpipi->SetYTitle("p_{n#pi^{+}#pi^{-}}^{CM} [GeV/c]");
  MomX_IMnpipi->SetStats(0);
  line = new TLine(Kp_mass, 0, Kp_mass, 1.5);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->Draw();
  line = new TLine(piS_mass1, 0, piS_mass1, 1.5);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->Draw();
  line = new TLine(piS_mass2, 0, piS_mass2, 1.5);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->Draw();
  c6->Print("tmp.pdf)");

  //==========//

  TCanvas *t1 = new TCanvas("t1", "", 600, 600);
  gPad->SetLogz();
  Cosn_IMnpipi_acc->Draw("colz");
  Cosn_IMnpipi_acc->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  Cosn_IMnpipi_acc->SetYTitle("cos#theta_{n}^{CM}");
  Cosn_IMnpipi_acc->SetStats(0);
  line = new TLine(Kp_mass, -1, Kp_mass, 1);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->Draw();
  line = new TLine(piS_mass1, -1, piS_mass1, 1);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->Draw();
  line = new TLine(piS_mass2, -1, piS_mass2, 1);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->Draw();
  t1->Print("tmp2.pdf(");
  
  his = (TH1F*)Cosn_IMnpipi_acc->ProjectionX();

  TCanvas *t2 = new TCanvas("t2", "", 600, 600);
  his->Draw();
  his->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  width = his->GetBinWidth(0)*1000;
  his->SetYTitle(Form("counts per %.0f MeV/c^{2}", width));
  ymax = 1.05*his->GetMaximum();
  line = new TLine(Kp_mass, 0, Kp_mass, ymax);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->Draw();
  line = new TLine(piS_mass1, 0, piS_mass1, ymax);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->Draw();
  line = new TLine(piS_mass2, 0, piS_mass2, ymax);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->Draw();
  t2->Print("tmp2.pdf");

  TCanvas *t10 = new TCanvas("t10", "", 600, 600);
  gPad->SetLogz();
  Cosn_IMnppipi_acc->Draw("colz");
  //Cosn_IMnppipi_acc->Draw("cont2");
  Cosn_IMnppipi_acc->SetXTitle("IM(np#pi^{+}#pi^{-}) [GeV/c^{2}]");
  Cosn_IMnppipi_acc->SetYTitle("cos#theta_{n}^{CM}");
  Cosn_IMnppipi_acc->SetStats(0);
  line = new TLine(Kpp_mass, -1, Kpp_mass, 1);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->Draw();
  line = new TLine(piSp_mass1, -1, piSp_mass1, 1);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->Draw();
  line = new TLine(piSp_mass2, -1, piSp_mass2, 1);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->Draw();
  t10->Print("tmp2.pdf");

  his = (TH1F*)Cosn_IMnppipi_acc->ProjectionX();

  TCanvas *t20 = new TCanvas("t20", "", 600, 600);
  his->Draw();
  his->SetXTitle("IM(np#pi^{+}#pi^{-}) [GeV/c^{2}]");
  width = his->GetBinWidth(0)*1000;
  his->SetYTitle(Form("counts per %.0f MeV/c^{2}", width));
  ymax = his->GetMaximum();
  line = new TLine(Kpp_mass, 0, Kpp_mass, ymax);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->Draw();
  line = new TLine(piSp_mass1, 0, piSp_mass1, ymax);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->Draw();
  line = new TLine(piSp_mass2, 0, piSp_mass2, ymax);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->Draw();
  t20->Print("tmp2.pdf");
  
  TCanvas *u3 = new TCanvas("u3", "", 600, 600);
  gPad->SetLogz();
  IMnpipi_IMnppipi->Draw("colz");
  IMnpipi_IMnppipi->SetXTitle("IM(np#pi^{+}#pi^{-}) [GeV/c^{2}]");
  IMnpipi_IMnppipi->SetYTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  IMnpipi_IMnppipi->SetStats(0);
  line = new TLine(Kpp_mass, 1, Kpp_mass, 2);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->Draw();
  line = new TLine(piSp_mass1, 1, piSp_mass1, 2);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->Draw();
  line = new TLine(piSp_mass2, 1, piSp_mass2, 2);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->Draw();
  line = new TLine(2, Kp_mass, 3, Kp_mass);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->Draw();
  line = new TLine(2, piS_mass1, 3, piS_mass1);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->Draw();
  line = new TLine(2, piS_mass2, 3, piS_mass2);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->Draw();
  u3->Print("tmp2.pdf)");

  TFile *out = new TFile("out.root", "recreate");
  his = (TH1F*)Cosn_IMnpipi_acc->ProjectionX();
  his->Write();
  Cosn_IMnpipi_Sp->Write();
  Cosn_IMnpipi_Sm->Write();
  Cosn_IMnppipi_Sp->Write();
  Cosn_IMnppipi_Sm->Write();
  out->Close();

  return 1;
}
