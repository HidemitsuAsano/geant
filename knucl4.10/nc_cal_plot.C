#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>

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
#include <TProfile.h>
#include <TFractionFitter.h>

int nc_cal_plot()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  double mom[15];
  double mom_err[15];
  for( int j=0; j<15; j++ ){
    mom[j] = j*0.1+0.05;
    mom_err[j] = 0.05;
  }

#if 0
  // 0.0 < dE
  double eff_mom0[15] = {0, 19.7573, 12.3479, 6.72916, 4.77536, 3.8526, 4.04453, 4.23693, 4.14967, 3.8138, 4.55368, 4.0276, 4.47699, 4.39066, 4.65804};
  double eff_mom_err0[15] = {0, 0.467528, 0.386855, 0.29632, 0.253468, 0.227388, 0.233863, 0.239781, 0.236537, 0.226793, 0.246019, 0.233312, 0.24309, 0.243054, 0.249995};
#else
  // 1.0 < dE
  double eff_mom0[15] = {0, 17.9512, 10.6748, 5.94572, 4.22436, 3.64322, 3.67813, 3.81182, 3.8824, 3.4212, 4.26124, 3.68962, 4.00719, 3.92626, 4.30622};
  double eff_mom_err0[15] = {0, 0.450634, 0.363109, 0.279704, 0.239085, 0.221363, 0.223444, 0.227938, 0.229112, 0.215241, 0.238352, 0.223701, 0.230547, 0.230398, 0.240812};
#endif

  // 2.0 < dE
  double eff_mom1[15] = {0, 14.3665, 9.72069, 5.40011, 3.65923, 3.25237, 3.42446, 3.37254, 3.55887, 3.09871, 3.84348, 3.32348, 3.68937, 3.60259, 3.8559};
  double eff_mom_err1[15] = {0, 0.411849, 0.348348, 0.267334, 0.223174, 0.209576, 0.215885, 0.214891, 0.219727, 0.205187, 0.226861, 0.212714, 0.221582, 0.221069, 0.228408};
  
  // 5.0 < dE
  double eff_mom2[15] = {0, 7.39005, 8.25498, 4.63067, 3.17886, 2.81965, 3.02988, 2.9616, 2.9962, 2.80426, 3.30038, 2.85875, 3.21957, 3.13819, 3.36335};
  double eff_mom_err2[15] = {0, 0.30718, 0.323609, 0.248561, 0.208529, 0.195573, 0.203481, 0.201802, 0.202197, 0.195491, 0.210815, 0.197756, 0.207498, 0.206825, 0.213867};

  // 5.0 < dE, scaled
  double eff_mom3[15] = {0, 0.344685, 2.79314, 2.57415, 1.80842, 1.54941, 1.76156, 1.71461, 1.63173, 1.69658, 2.1167, 1.76032, 2.12795, 2.12496, 2.33605};
  double eff_mom_err3[15] = {0, 0.0688181, 0.193761, 0.18731, 0.158391, 0.14592, 0.156164, 0.154532, 0.150262, 0.15292, 0.16986, 0.156055, 0.169641, 0.17108, 0.179183};

  // from Michael's MC (dE> 2 MeV)
  double MKin[9] = {5.5, 8.5, 12.5, 16.5, 21, 41, 61, 82.5, 105}; // MeV
  double Meff[9] = {0.158, 0.134, 0.114, 0.109, 0.108, 0.053, 0.043, 0.041, 0.036};
  double Mmom[9];
  const double nmass = 939.565346;
  for( int i=0; i<9; i++ ){
    Mmom[i] = sqrt((nmass+MKin[i])*(nmass+MKin[i])-nmass*nmass)/1000;
    Meff[i] = Meff[i]*100;
  }



  TCanvas *c1 = new TCanvas("c1", "", 900, 900);
  c1->SetGrid();
  TH2F *frame = new TH2F("frame", "", 100, 0, 1.5, 100, 0, 20);
  frame->Draw();
  frame->SetXTitle("momentum [GeV/c]");
  frame->SetYTitle("efficiency [%]");
  TGraph *gr0 = new TGraphErrors(15, mom, eff_mom0, mom_err, eff_mom_err0);
  gr0->SetMarkerStyle(20); gr0->SetMarkerColor(1);
  gr0->SetMarkerSize(2); gr0->SetLineColor(1);
  gr0->Draw("same,P");
  TGraph *gr1 = new TGraphErrors(15, mom, eff_mom1, mom_err, eff_mom_err1);
  gr1->SetMarkerStyle(20); gr1->SetMarkerColor(3);
  gr1->SetMarkerSize(2); gr1->SetLineColor(3);
  gr1->Draw("same,P");
  TGraph *gr2 = new TGraphErrors(15, mom, eff_mom2, mom_err, eff_mom_err2);
  gr2->SetMarkerStyle(20); gr2->SetMarkerColor(2);
  gr2->SetMarkerSize(2); gr2->SetLineColor(2);
  gr2->Draw("same,P");
#if 0
  TGraph *gr3 = new TGraphErrors(15, mom, eff_mom3, mom_err, eff_mom_err3);
  gr3->SetMarkerStyle(20); gr3->SetMarkerColor(4);
  gr3->SetMarkerSize(2); gr3->SetLineColor(4);
  gr3->Draw("same,P");
#endif
  TGraph *gr4 = new TGraph(9, Mmom, Meff);
  gr4->SetMarkerStyle(24); gr4->SetMarkerColor(3);
  gr4->SetMarkerSize(2); gr4->SetLineColor(3);
  gr4->Draw("same,P");
  gr4->Print();

  TLegend *leg = new TLegend(0.4, 0.6, 0.9, 0.9);
#if 0
  leg->AddEntry(gr0, "0.0 MeVee < dE", "LP");
#else
  leg->AddEntry(gr0, "1.0 MeVee < dE", "LP");
#endif
  leg->AddEntry(gr1, "2.0 MeVee < dE", "LP");
  leg->AddEntry(gr2, "5.0 MeVee < dE", "LP");
  //leg->AddEntry(gr3, "5.0 MeVee < dE, scaled", "LP");
  leg->AddEntry(gr4, "2.0 MeVee < dE, from Michael's (KLOE) MC", "P");
  leg->Draw();


  c1->Print("tmp.pdf)");

  
  return true;
}
