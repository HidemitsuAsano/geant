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

#include "include/KnuclRootData.h"

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

int main( int argc, char** argv )
{

  if( argc!=3 ){
    std::cerr<<"usage: "<<argv[0]<<" [input.root] [output.root]"<<std::endl;
    return -1;
  }
  
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

  TLine *line;
  double ymax;
  char com[32];

  int    offset[3] = {0, -1, 1};
  int    id = 0;
  int    PISIGMA_BIN = 46;
  double PISIGMA_MIN = 1.21+offset[id]*0.005;
  double PISIGMA_MAX = 1.90+offset[id]*0.005;

  const int    PISIGMAP_BIN = 50;
  const double PISIGMAP_MIN = 2.0;
  const double PISIGMAP_MAX = 3.0;

  const int    COSN_BIN = 40;
  const double COSN_MIN = -1.0;
  const double COSN_MAX = 1.0;

  const int    QKN_BIN = 30;
  const double QKN_MIN = 0.0;
  const double QKN_MAX = 1.5;

  // [3] = {S+, S-, sum}
  TH2F *Cosn_IMnpipi[3];
  TH2F *Cosn_IMnppipi[3];

  TH2F *Qkn_IMnpipi[3];
  TH2F *Qkn_IMnppipi[3];

  TH2F *IMnpipi_IMnppipi[3];
    
  //=== acceptance study ===//
  TH2F *Cosn_IMnpipi_acc[3];
  TH2F *Cosn_IMnppipi_acc[3];

  TH2F *Qkn_IMnpipi_acc[3];
  TH2F *Qkn_IMnppipi_acc[3];
  //=== acceptance study ===//
  
  for( int x=0; x<3; x++ ){
    sprintf( com, "Cosn_IMnpipi_%d", x );
    Cosn_IMnpipi[x] = new TH2F( com, com, PISIGMA_BIN, PISIGMA_MIN, PISIGMA_MAX, COSN_BIN, COSN_MIN, COSN_MAX );
    Cosn_IMnpipi[x]->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
    Cosn_IMnpipi[x]->SetYTitle("cos#theta_{#Lambda}^{CM}");
    sprintf( com, "Cosn_IMnppipi_%d", x );
    Cosn_IMnppipi[x] = new TH2F( com, com, PISIGMAP_BIN, PISIGMAP_MIN, PISIGMAP_MAX, COSN_BIN, COSN_MIN, COSN_MAX );
    Cosn_IMnppipi[x]->SetXTitle("IM(n#pi^{+}#pi^{-}p) [GeV/c^{2}]");
    Cosn_IMnppipi[x]->SetYTitle("cos#theta_{#Lambda}^{CM}");
    
    sprintf( com, "Qkn_IMnpipi_%d", x );
    Qkn_IMnpipi[x] = new TH2F( com, com, PISIGMA_BIN, PISIGMA_MIN, PISIGMA_MAX, QKN_BIN, QKN_MIN, QKN_MAX );
    Qkn_IMnpipi[x]->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
    Qkn_IMnpipi[x]->SetYTitle("q_{Kn} [GeV/c]");
    sprintf( com, "Qkn_IMnppipi_%d", x );
    Qkn_IMnppipi[x] = new TH2F( com, com, PISIGMAP_BIN, PISIGMAP_MIN, PISIGMAP_MAX, QKN_BIN, QKN_MIN, QKN_MAX );
    Qkn_IMnppipi[x]->SetXTitle("IM(n#pi^{+}#pi^{-}p) [GeV/c^{2}]");
    Qkn_IMnppipi[x]->SetYTitle("q_{Kn} [GeV/c]");

    sprintf( com, "IMnpipi_IMnppipi_%d", x );
    IMnpipi_IMnppipi[x] = new TH2F( com, com, PISIGMAP_BIN, PISIGMAP_MIN, PISIGMAP_MAX, PISIGMA_BIN, PISIGMA_MIN, PISIGMA_MAX );
    IMnpipi_IMnppipi[x]->SetXTitle("IM(n#pi^{+}#pi^{-}p) [GeV/c^{2}]");
    IMnpipi_IMnppipi[x]->SetYTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");

    
    //=== acceptance study ===//
    sprintf( com, "Cosn_IMnpipi_acc_%d", x );
    Cosn_IMnpipi_acc[x] = new TH2F( com, com, PISIGMA_BIN, PISIGMA_MIN, PISIGMA_MAX, COSN_BIN, COSN_MIN, COSN_MAX );
    Cosn_IMnpipi_acc[x]->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
    Cosn_IMnpipi_acc[x]->SetYTitle("cos#theta_{#Lambda}^{CM}");
    sprintf( com, "Cosn_IMnppipi_acc_%d", x );
    Cosn_IMnppipi_acc[x] = new TH2F( com, com, PISIGMAP_BIN, PISIGMAP_MIN, PISIGMAP_MAX, COSN_BIN, COSN_MIN, COSN_MAX );
    Cosn_IMnppipi_acc[x]->SetXTitle("IM(n#pi^{+}#pi^{-}p) [GeV/c^{2}]");
    Cosn_IMnppipi_acc[x]->SetYTitle("cos#theta_{#Lambda}^{CM}");
    
    sprintf( com, "Qkn_IMnpipi_acc_%d", x );
    Qkn_IMnpipi_acc[x] = new TH2F( com, com, PISIGMA_BIN, PISIGMA_MIN, PISIGMA_MAX, QKN_BIN, QKN_MIN, QKN_MAX );
    Qkn_IMnpipi_acc[x]->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
    Qkn_IMnpipi_acc[x]->SetYTitle("q_{Kn} [GeV/c]");
    sprintf( com, "Qkn_IMnppipi_acc_%d", x );
    Qkn_IMnppipi_acc[x] = new TH2F( com, com, PISIGMAP_BIN, PISIGMAP_MIN, PISIGMAP_MAX, QKN_BIN, QKN_MIN, QKN_MAX );
    Qkn_IMnppipi_acc[x]->SetXTitle("IM(n#pi^{+}#pi^{-}p) [GeV/c^{2}]");
    Qkn_IMnppipi_acc[x]->SetYTitle("q_{Kn} [GeV/c]");
    //=== acceptance study ===//
  }


  int FermiFlag = 0;
  std::string initname[2];
  std::string name[6];

  TFile *f = new TFile(argv[1]);
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

    //=== acceptance study ===//
    // beam_K(K+), pi+, pi-, p, n, Lambda(missing)
    int reactionID = reactionData->ReactionID();
    int PDG[6] = {321, 211, -211, 2212, 2112, 3122};
    int parentID[6] = {0,0,0,0,0,0};
    int ID[6]       = {-1,-1,-1,-1,-1,-1};
    int trackID[6]  = {-1,-1,-1,-1,-1,-1};
    nparticle = 0;
    for (Int_t j=0; j<mcData->trackSize(); j++) {
      int pdgcode = mcData->track(j)->pdgID();
      int parent  = mcData->track(j)->parentTrackID();
      int track   = mcData->track(j)->trackID();
      for( int k=0; k<6; k++ ){
        if( pdgcode==PDG[k] && parent==parentID[k] && ID[k]==-1 ){
          ID[k] = j;
          trackID[k] = track;
          nparticle++;
	  break;
	}
      }
    } // for (Int_t j=0; j<mcData->trackSize(); j++) {

    bool flagG4Decay = (nparticle==6) ? true : false;
    int nCDHhit[6] = {0, 0, 0, 0, 0, 0};
    if( flagG4Decay ){
      for (Int_t j=0; j<detectorData->detectorHitSize(); j++) {
	int cid     = detectorData->detectorHit(j)->detectorID();
	int track   = detectorData->detectorHit(j)->trackID();
	for( int k=1; k<5; k++ ){
	  if( cid==CID_CDH && track==trackID[k] ) nCDHhit[k]++;
	}
      }
    }
    bool piSpn_detect = (nCDHhit[1] && nCDHhit[2] && nCDHhit[3] && nCDHhit[4]) ?  true : false;
    //=== acceptance study ===//
    
    TLorentzVector mom_beam   = reactionData->GetInitParticle(0);
    TLorentzVector mom_target = reactionData->GetInitParticle(1);
    TLorentzVector mom_target_beam = mom_target+mom_beam;
    TVector3 boost = mom_target_beam.BoostVector();
    TLorentzVector mom_beam_CM = mom_beam;
    mom_beam_CM.Boost(-boost);

    TLorentzVector mom_pi_CM     = reactionData->GetCMParticle(3); // pi+
    TLorentzVector mom_Lambda_CM = reactionData->GetCMParticle(0);
    TLorentzVector mom_p_CM      = reactionData->GetCMParticle(1);
    TLorentzVector mom_n_CM      = reactionData->GetCMParticle(2);
    TLorentzVector mom_pi2_CM    = reactionData->GetCMParticle(5); // pi-
    TLorentzVector mom_X_CM      = mom_pi_CM+mom_Lambda_CM;

    double cos_n = mom_n_CM.Vect().Dot(mom_beam_CM.Vect())/(mom_n_CM.Vect().Mag()*mom_beam_CM.Vect().Mag());
    double cos_p = mom_p_CM.Vect().Dot(mom_beam_CM.Vect())/(mom_p_CM.Vect().Mag()*mom_beam_CM.Vect().Mag());
    double cos_X = mom_X_CM.Vect().Dot(mom_beam_CM.Vect())/(mom_X_CM.Vect().Mag()*mom_beam_CM.Vect().Mag());

    TLorentzVector mom_n = mom_n_CM;
    mom_n.Boost(boost);
    TLorentzVector qkn = mom_beam-mom_n;
    
    //TLorentzVector mom_piS  = mom_pi_CM+mom_Lambda_CM;
    //TLorentzVector mom_piSp = mom_pi_CM+mom_Lambda_CM+mom_p_CM;
    // change to BG study of piSpn
    TLorentzVector mom_piS  = mom_pi_CM+mom_pi2_CM+mom_n_CM;
    TLorentzVector mom_piSp = mom_pi_CM+mom_pi2_CM+mom_n_CM+mom_p_CM;
    
    Cosn_IMnpipi[2]->Fill(mom_piS.M()/1000, cos_n);
    Cosn_IMnppipi[2]->Fill(mom_piSp.M()/1000, cos_n);
    
    Qkn_IMnpipi[2]->Fill(mom_piS.M()/1000, qkn.P()/1000.0);
    Qkn_IMnppipi[2]->Fill(mom_piSp.M()/1000, qkn.P()/1000.0);

    IMnpipi_IMnppipi[2]->Fill(mom_piSp.M()/1000, mom_piS.M()/1000.0);
    
    if( reactionData->ReactionID()==2120 ){ // 2120: K- 3He -> S+ p n p-
      Cosn_IMnpipi[0]->Fill(mom_piS.M()/1000, cos_n);
      Cosn_IMnppipi[0]->Fill(mom_piSp.M()/1000, cos_n);

      Qkn_IMnpipi[0]->Fill(mom_piS.M()/1000, qkn.P()/1000.0);
      Qkn_IMnppipi[0]->Fill(mom_piSp.M()/1000, qkn.P()/1000.0);

      IMnpipi_IMnppipi[0]->Fill(mom_piSp.M()/1000, mom_piS.M()/1000.0);
    }
    else if( reactionData->ReactionID()==2130 ){ // 130: K- 3He -> S- p n p+
      Cosn_IMnpipi[1]->Fill(mom_piS.M()/1000, cos_n);
      Cosn_IMnppipi[1]->Fill(mom_piSp.M()/1000, cos_n);

      Qkn_IMnpipi[1]->Fill(mom_piS.M()/1000, qkn.P()/1000.0);
      Qkn_IMnppipi[1]->Fill(mom_piSp.M()/1000, qkn.P()/1000.0);
      
      IMnpipi_IMnppipi[1]->Fill(mom_piSp.M()/1000, mom_piS.M()/1000.0);
    }

    //=== acceptance study ===//
    if (piSpn_detect) {
      Cosn_IMnpipi_acc[2]->Fill(mom_piS.M()/1000, cos_n);
      Cosn_IMnppipi_acc[2]->Fill(mom_piSp.M()/1000, cos_n);
      
      Qkn_IMnpipi_acc[2]->Fill(mom_piS.M()/1000, qkn.P()/1000.0);
      Qkn_IMnppipi_acc[2]->Fill(mom_piSp.M()/1000, qkn.P()/1000.0);

      if( reactionData->ReactionID()==2120 ){ // 2120: K- 3He -> S+ p n p-
	Cosn_IMnpipi_acc[0]->Fill(mom_piS.M()/1000, cos_n);
	Cosn_IMnppipi_acc[0]->Fill(mom_piSp.M()/1000, cos_n);

	Qkn_IMnpipi_acc[0]->Fill(mom_piS.M()/1000, qkn.P()/1000.0);
	Qkn_IMnppipi_acc[0]->Fill(mom_piSp.M()/1000, qkn.P()/1000.0);
      }
      else if( reactionData->ReactionID()==2130 ){ // 130: K- 3He -> S- p n p+
	Cosn_IMnpipi_acc[1]->Fill(mom_piS.M()/1000, cos_n);
	Cosn_IMnppipi_acc[1]->Fill(mom_piSp.M()/1000, cos_n);

	Qkn_IMnpipi_acc[1]->Fill(mom_piS.M()/1000, qkn.P()/1000.0);
	Qkn_IMnppipi_acc[1]->Fill(mom_piSp.M()/1000, qkn.P()/1000.0);
      }
    } // if (piSpn_detect) {
    //=== acceptance study ===//

    
  }// for (Int_t i=0;i<nevent;i++) {
  std::cout<<"end of filling"<<std::endl;
  //----------------------//
  //--- event roop end ---//
  //----------------------//

  //** ++ ** ++ ** ++ ** ++ ** ++ ** ++ ** ++ ** ++ **//
  int n;
  double width;
  
  TCanvas *c0;
  c0 = new TCanvas("c0", "", 900, 300);
  c0->Divide(3,1);
  n = 1;
  for( int x=0; x<3; x++ ){
    c0->cd(n); n++;
    gPad->SetLogz();
    Cosn_IMnpipi[x]->Draw("colz");
    Cosn_IMnpipi[x]->SetMinimum(1);
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
  }
  c0->Print("tmp.pdf(");
  TCanvas *c1;
  c1 = new TCanvas("c1", "", 900, 300);
  c1->Divide(3,1);
  n = 1;
  for( int x=0; x<3; x++ ){
    c1->cd(n); n++;
    gPad->SetLogz();
    Qkn_IMnpipi[x]->Draw("colz");
    Qkn_IMnpipi[x]->SetMinimum(1);
    line = new TLine(Kp_mass, 0, Kp_mass, 1.5);
    line->SetLineColor(1);
    line->SetLineStyle(2);
    line->Draw();
  }
  c1->Print("tmp.pdf");
  TCanvas *c2;
  c2 = new TCanvas("c2", "", 900, 300);
  c2->Divide(3,1);
  n = 1;
  for( int x=0; x<3; x++ ){
    c2->cd(n); n++;
    Cosn_IMnpipi[x]->ProjectionX()->Draw("e");
    c2->Update();
    ymax = gPad->GetUymax();
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
    width = Cosn_IMnpipi[x]->ProjectionX()->GetBinWidth(0)*1000;
    Cosn_IMnpipi[x]->ProjectionX()->SetYTitle(Form("counts per %.0f MeV/c^{2}", width));
  }
  c2->Print("tmp.pdf");
  TCanvas *c3;
  c3 = new TCanvas("c3", "", 900, 300);
  c3->Divide(3,1);
  n = 1;
  for( int x=0; x<3; x++ ){
    c3->cd(n); n++;
    sprintf(com,"tmp1_%d",x);
    Cosn_IMnpipi[x]->ProjectionX(com,35,40)->Draw("e");
    c3->Update();
    ymax = gPad->GetUymax();
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
    width = Cosn_IMnpipi[x]->ProjectionX()->GetBinWidth(0)*1000;
    Cosn_IMnpipi[x]->ProjectionX(com,35,40)->SetYTitle(Form("counts per %.0f MeV/c^{2}", width));
  }
  c3->Print("tmp.pdf");
  TCanvas *c4;
  c4 = new TCanvas("c4", "", 900, 300);
  c4->Divide(3,1);
  n = 1;
  for( int x=0; x<3; x++ ){
    c4->cd(n); n++;
    sprintf(com,"tmp2_%d",x);
    Qkn_IMnpipi[x]->ProjectionX(com,1,13)->Draw("e");
    c4->Update();
    ymax = gPad->GetUymax();
    line = new TLine(Kp_mass, 0, Kp_mass, ymax);
    line->SetLineColor(1);
    line->SetLineStyle(2);
    line->Draw();
    width = Qkn_IMnpipi[x]->ProjectionX()->GetBinWidth(0)*1000;
    Qkn_IMnpipi[x]->ProjectionX(com,1,13)->SetYTitle(Form("counts per %.0f MeV/c^{2}", width));
  }
  c4->Print("tmp.pdf)");


  TCanvas *d0;
  d0 = new TCanvas("d0", "", 900, 300);
  d0->Divide(3,1);
  n = 1;
  for( int x=0; x<3; x++ ){
    d0->cd(n); n++;
    gPad->SetLogz();
    Cosn_IMnppipi[x]->Draw("colz");
    Cosn_IMnppipi[x]->SetMinimum(1);
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
  }
  d0->Print("tmp2.pdf(");
  TCanvas *d1;
  d1 = new TCanvas("d1", "", 900, 300);
  d1->Divide(3,1);
  n = 1;
  for( int x=0; x<3; x++ ){
    d1->cd(n); n++;
    gPad->SetLogz();
    Qkn_IMnppipi[x]->Draw("colz");
    Qkn_IMnppipi[x]->SetMinimum(1);
    line = new TLine(Kpp_mass, 0, Kpp_mass, 1.5);
    line->SetLineColor(1);
    line->SetLineStyle(2);
    line->Draw();
  }
  d1->Print("tmp2.pdf");
  TCanvas *d2;
  d2 = new TCanvas("d2", "", 900, 300);
  d2->Divide(3,1);
  n = 1;
  for( int x=0; x<3; x++ ){
    d2->cd(n); n++;
    Cosn_IMnppipi[x]->ProjectionX()->Draw("e");
    d2->Update();
    ymax = gPad->GetUymax();
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
    width = Cosn_IMnppipi[x]->ProjectionX()->GetBinWidth(0)*1000;
    Cosn_IMnppipi[x]->ProjectionX()->SetYTitle(Form("counts per %.0f MeV/c^{2}", width));
  }
  d2->Print("tmp2.pdf");
  TCanvas *d3;
  d3 = new TCanvas("d3", "", 900, 300);
  d3->Divide(3,1);
  n = 1;
  for( int x=0; x<3; x++ ){
    d3->cd(n); n++;
    sprintf(com,"tmp3_%d",x);
    Cosn_IMnppipi[x]->ProjectionX(com,35,40)->Draw("e");
    d3->Update();
    ymax = gPad->GetUymax();
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
    width = Cosn_IMnppipi[x]->ProjectionX()->GetBinWidth(0)*1000;
    Cosn_IMnppipi[x]->ProjectionX(com,35,40)->SetYTitle(Form("counts per %.0f MeV/c^{2}", width));
  }
  d3->Print("tmp2.pdf");
  TCanvas *d4;
  d4 = new TCanvas("d4", "", 900, 300);
  d4->Divide(3,1);
  n = 1;
  for( int x=0; x<3; x++ ){
    d4->cd(n); n++;
    sprintf(com,"tmp4_%d",x);
    Qkn_IMnppipi[x]->ProjectionX(com,1,13)->Draw("e");
    d4->Update();
    ymax = gPad->GetUymax();
    line = new TLine(Kpp_mass, 0, Kpp_mass, ymax);
    line->SetLineColor(1);
    line->SetLineStyle(2);
    line->Draw();
    width = Qkn_IMnppipi[x]->ProjectionX()->GetBinWidth(0)*1000;
    Qkn_IMnppipi[x]->ProjectionX(com,1,13)->SetYTitle(Form("counts per %.0f MeV/c^{2}", width));
  }
  d4->Print("tmp2.pdf)");


  TCanvas *e0;
  e0 = new TCanvas("e0", "", 800, 800);
  e0->Divide(2,2);
  for( int x=2; x<3; x++ ){
    e0->cd(3);
    //gPad->SetLogz();
    Cosn_IMnpipi[x]->Draw("colz");
    Cosn_IMnpipi[x]->SetMinimum(1);
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
  }
  for( int x=2; x<3; x++ ){
    e0->cd(1);
    Cosn_IMnpipi[x]->ProjectionX()->Draw("e");
    e0->Update();
    ymax = gPad->GetUymax();
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
    width = Cosn_IMnpipi[x]->ProjectionX()->GetBinWidth(0)*1000;
    Cosn_IMnpipi[x]->ProjectionX()->SetYTitle(Form("counts per %.0f MeV/c^{2}", width));
  }
  for( int x=2; x<3; x++ ){
    e0->cd(4);
    //gPad->SetLogz();
    Cosn_IMnppipi[x]->Draw("colz");
    Cosn_IMnppipi[x]->SetMinimum(1);
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
  }
  for( int x=2; x<3; x++ ){
    e0->cd(2);
    Cosn_IMnppipi[x]->ProjectionX()->Draw("e");
    e0->Update();
    ymax = gPad->GetUymax();
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
    width = Cosn_IMnppipi[x]->ProjectionX()->GetBinWidth(0)*1000;
    Cosn_IMnppipi[x]->ProjectionX()->SetYTitle(Form("counts per %.0f MeV/c^{2}", width));
  }
  e0->Print("tmp3.pdf");

  
  TFile *out = new TFile(argv[2], "recreate");
  for( int x=0; x<3; x++ ){
    Cosn_IMnpipi[x]->Write();
    Cosn_IMnppipi[x]->Write();
    Qkn_IMnpipi[x]->Write();
    Qkn_IMnppipi[x]->Write();
    IMnpipi_IMnppipi[x]->Write();
    //=== acceptance study ===//
    Cosn_IMnpipi_acc[x]->Write();
    Cosn_IMnppipi_acc[x]->Write();
    Qkn_IMnpipi_acc[x]->Write();
    Qkn_IMnppipi_acc[x]->Write();
    //=== acceptance study ===//
  }
  out->Close();

}
