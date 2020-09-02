#include <fstream>
#include <iostream>
#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TMultiGraph.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TLatex.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>

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

void pi_single2(const char* filename="")
{
  gSystem->Load("KnuclRootData_cc.so"); //<-- OR read from ./.rootlogon.C
  gSystem->Load("~/ana/k18ana/src/lib/libAll.so"); //<-- OR read from ./.rootlogon.C
  
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(111111);
  
  gROOT->cd();

  TDatabasePDG *pdg = new TDatabasePDG();
  pdg->ReadPDGTable("pdg_table.txt");
  
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

  const double threshold = 2.0; // MeVee

  tree2->GetEvent(0);

  Int_t nevent = tree->GetEntries();
  //nevent = 1000;
  TH1::SetDefaultSumw2();

  TH1F* his0 = new TH1F("CDH_dE", "CDH_dE", 1000, 0, 100);
  TH1F* his1 = new TH1F("CDH_dE_cut","CDH_dE_cut",  1000, 0, 100);
  TH1F* CDHmul = new TH1F("CDH_mul","CDH_mul",10,0,10);
  TH1F* diff_CDH = new TH1F("diff_CDH","diff_CDH",73,-36.5,73.5);
  TH1F* his2[36];
  TH1F* his3[36];
  for( int i=0; i<36; i++ ){
    his2[i] = new TH1F(Form("CDH%d_dE",i),Form("CDH%d_dE",i), 1000, 0, 100);
    his3[i] = new TH1F(Form("CDH%d_dE_cut",i), Form("CDH%d_dE_cut",i), 1000, 0, 100);
  }
  TH1F* his4 = new TH1F("dE_max_seg", "dE_max_seg", 36, -0.5, 35.5);

  double mom[10];
  double mom_err[10];
  int num_mom[10];
  int acc_mom[10];
  double eff_mom[10];
  double eff_mom_err[10];
  for( int j=0; j<10; j++ ){
    mom[j] = j*0.1+0.05;
    mom_err[j] = 0.05;
    num_mom[j] = 0;
    acc_mom[j] = 0;
    eff_mom[j] = 0;
    eff_mom_err[j] = 0;
  }
  TH1F* his20[36];
  TH1F* his30[36];
  for( int i=0; i<36; i++ ){
    his20[i] = new TH1F(Form("CDH%d_dE_end",i),Form("CDH%d_dE",i), 1000, 0, 100);
    his30[i] = new TH1F(Form("CDH%d_dE_end_cut",i),Form("CDH%d_dE_cut",i), 1000, 0, 100);
  }
  
  TH1F* hgen = new TH1F("hgen","hgen",50,0,1);
  hgen->SetXTitle("momentum [GeV/c]");
  TH1F* hdet = new TH1F("hdet","hdet",50,0,1);
  hdet->SetXTitle("momentum [GeV/c]");
  TH1F* heff = new TH1F("heff","heff",50,0,1);
  heff->SetTitle("neutron efficiency");
  heff->SetXTitle("momentum [GeV/c]");
  //parent particle comes from vertex.
  TH1F* hdet_pri = new TH1F("hdet_pri","hdet_pri",50,0,1);
  hdet_pri->SetXTitle("momentum [GeV/c]");
  TH1F* heff_pri = new TH1F("heff_pri","heff_pri",50,0,1);
  heff_pri->SetXTitle("momentum [GeV/c]");

  //------------------------//
  //--- event roop start ---//
  //------------------------//
  int count0 = 0;
  int count1 = 0;
  cout<<"Start to fill histgrams. Entries = "<<nevent<<endl;
  for (Int_t i=0; i<nevent; i++) {
  //for (Int_t i=0; i<500; i++) {
    if(i%5000==0) std::cout << i  << std::endl;
    if(i==30000) break;
    tree->GetEvent(i);
    double momentum = reactionData->GetParticle(0).P();
    //std::cout << reactionData->GetParticle(0).Phi() << std::endl; 
    hgen->Fill(momentum/1000.);
    int m = int(momentum/100);
    num_mom[m]++;
 
    double ene_dep = 0;
    double de_seg[36] = {0};
    int cdhmul=0;
    for (Int_t j=0; j<detectorData->detectorHitSize(); j++) {
      int cid = detectorData->detectorHit(j)->detectorID();
      if( cid == CID_CDH ){
        int seg = detectorData->detectorHit(j)->channelID();
        double dE = detectorData->detectorHit(j)->adc();
        DetectorHit *dhit = detectorData->detectorHit(j);
        ene_dep += dE;
        de_seg[seg] = dE;
        if(dE>threshold){
          cdhmul++;
          int pdgid = dhit->pdg();
          int trackid = dhit->trackID();
          Track *track_p  = FindTrackFromMcIndex(mcData,trackid);
          int parentID = dhit->parentID();
          Track *parenttrack_p = FindTrackFromMcIndex(mcData,parentID);

          int parentpdg = 0;
          if(parenttrack_p !=0){
            parentpdg= parenttrack_p->pdgID();
          }
          Tools::H1(Form("hit_pdg"),pdgid,16000,-8000,8000);
          Tools::H1(Form("hit_parentpdg"),parentpdg,16000,-8000,8000);
          TVector3 ParentVtx = FillAncestryVertex(mcData,detectorData->detectorHit(j),  dE);
          hdet->Fill(momentum/1000.0);
          if(fabs(ParentVtx.Mag()<1.0))hdet_pri->Fill(momentum/1000.0);
        }
      }
    }
    CDHmul->Fill(cdhmul);
    //if(ene_dep) cerr<<ene_dep<<" -> ";
    //if( dE_scale_flag ) ene_dep = ene_dep*dE_par1*(exp(ene_dep/dE_par2));
    //if(ene_dep) cerr<<ene_dep<<endl;

    if( 0 < ene_dep ){
      his0->Fill(ene_dep);
      for( int j=0; j<36; j++ ) his2[j]->Fill(de_seg[j]);
      his20[m]->Fill(ene_dep);
      count0++;
      if( threshold < ene_dep ){
        his1->Fill(ene_dep);
        for( int j=0; j<36; j++ ) his3[j]->Fill(de_seg[j]);
        his30[m]->Fill(ene_dep);
        count1++;
        //if(ParentVtxR<10.0)hdet_pri 
        acc_mom[m]++;
      
	double ene_max = 0;
	double max_seg = -1;
	for (Int_t j=0; j<36; j++) {
	  if( ene_max < de_seg[j] ){
	    max_seg = j;
	    ene_max = de_seg[j];
	  }
	}
	his4->Fill(max_seg);
      } // if( threshold < ene_dep )
    } // if( 0 < ene_dep )

  }// for (Int_t i=0;i<nevent;i++) 
  cout<<"end of filling"<<endl;
  //----------------------//
  //--- event roop end ---//
  //----------------------//
  cerr<<"# of hit = "<<count0<<" / " <<count1<<endl;

  for( int j=1; j<10; j++ ){
    eff_mom[j] = (double)acc_mom[j]/num_mom[j];
    eff_mom_err[j] = sqrt(eff_mom[j]*(1-eff_mom[j])/num_mom[j]);
    eff_mom[j] *= 100; // percent
    eff_mom_err[j] *= 100; // percent
  }
  //cerr<<mom[j]<<" "<<eff_mom[j]<<" +/- "<<eff_mom_err[j]<<endl;
  for( int j=1; j<10; j++ ) cerr<<eff_mom[j]<<", ";
  cerr<<endl;
  for( int j=1; j<10; j++ ) cerr<<eff_mom_err[j]<<", ";
  cerr<<endl;

  //--- plot ---//
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(0);

  TCanvas *c1 = new TCanvas("c1", "", 1200, 600);
  c1->Divide(3,1);
  c1->cd(1); his0->Draw("");
  his0->SetXTitle("dE (MeV)");
  c1->cd(2); his1->Draw("");
  his1->SetXTitle("dE (MeV)");
  c1->cd(3); his4->Draw("");
  his4->SetXTitle("dE maximum segment");
  //c1->Print("tmp.pdf(");

#if 1
  TCanvas *c3 = new TCanvas("c3", "", 800, 800);
  c3->Divide(4,4);
  for( int i=0; i<16; i++ ){
    c3->cd(i+1); his2[i]->Draw("");
    his2[i]->SetXTitle("dE (MeV)");
  }
  //c3->Print("tmp.pdf");

  TCanvas *c4 = new TCanvas("c4", "", 800, 800);
  c4->Divide(4,4);
  for( int i=0; i<16; i++ ){
    c4->cd(i+1); his3[i]->Draw("");
    his3[i]->SetXTitle("dE (MeV)");
  }
 // c4->Print("tmp.pdf");
#endif

#if 1
  TCanvas *c5 = new TCanvas("c5", "", 800, 800);
  c5->Divide(4,4);
  for( int i=0; i<10; i++ ){
    c5->cd(i+1); his20[i]->Draw("");
    his20[i]->SetXTitle("dE (MeV)");
  }
 // c5->Print("tmp.pdf");

  TCanvas *c6 = new TCanvas("c6", "", 1200, 1200);
  c6->Divide(4,4);
  for( int i=0; i<10; i++ ){
    c6->cd(i+1); his30[i]->Draw("");
    his30[i]->SetXTitle("dE (MeV)");
  }
//  c6->Print("tmp.pdf");
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
//  c2->Print("tmp.pdf)");
  
  TCanvas *cgen = new TCanvas("cgen","cgen");
  cgen->cd();
  hgen->Draw("HE");

  TCanvas *cdet = new TCanvas("cdet","cdet");
  cdet->cd();
  hdet->Draw("HE");
  hdet_pri->SetLineColor(2);
  hdet_pri->Draw("HEsame");
  
  TCanvas *ceff = new TCanvas("ceff","ceff");
  ceff->cd();
  heff = (TH1F)hdet->Clone();
  heff_pri = (TH1F)hdet_pri->Clone();
  heff->Divide(hdet,hgen,1.0,1.0,"b");
  heff_pri->Divide(hdet_pri,hgen,1.0,1.0,"b");
  heff->Draw("HE");
  heff_pri->Draw("HEsame");
  
  TIter nexthist2(gDirectory->GetList());
  TFile *fout = new TFile("pimout2.root","RECREATE");
  fout->cd();
  while( (obj = (TObject*)nexthist2())!=0  ){
    obj->Write();
  }
  fout->Close();
}

Track* FindTrackFromMcIndex(MCData *mcdata, int trackid)
{
  for( int itr=0;itr<mcdata->trackSize();itr++){
    Track *track=mcdata->track(itr);
    if( track->trackID()==trackid) return track;
  }
  return 0;
}


TVector3 FillAncestryVertex(MCData *mcdata,DetectorHit *dhit, double dE)
{
  if(dE<0.01) return 0;
  Track *parentTr= (Track*)FindTrackFromMcIndex(mcdata, dhit->parentID());
  
  TVector3 vtxp;
  if(parentTr) vtxp = parentTr->vertex();
  else return vtxp;
  Tools::H2(Form("dE_track_vtxr_ncan_1stg"),vtxp.Perp()/10.,dE,1200,0,120,100,0,10);
  int gen=0;
  while( parentTr!=0){
    TVector3 vtx = parentTr->vertex();
    Tools::H2(Form("dE_track_vtxr_ncan"),vtx.Perp()/10.,dE,1200,0,120,100,0,10);
    Tools::H2(Form("time_track_vtxr_ncan"),vtx.Perp()/10.,dhit->time(),1200,0,120,1000,0,100);
    Tools::H2(Form("gen_track_vtxr_ncan"),vtx.Perp()/10.,gen+1,1200,0,120,10,0,10);
    parentTr=FindTrackFromMcIndex(mcdata,parentTr->parentTrackID());
    gen++;
  }
  //return vtxp.Perp()/10.;
  return vtxp;
}

