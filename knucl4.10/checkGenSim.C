#define DALITZ 1

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

void checkGenSim(const char *filename="sim2/sim_nSmpip_0000.root")
{
  std::string outname = std::string(filename);
  outname.insert(std::string(filename).size()-5,"_hist");
  std::cout << outname << std::endl;

  if(gSystem->Load("./build/KnuclRootData_cc.so")!=0) return; //<-- OR read from ./.rootlogon.C
  
  gROOT->SetStyle("Plain");
  //gStyle->SetOptStat(111111);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(111111);
  
  gROOT->cd();

  TDatabasePDG *pdg = new TDatabasePDG();
  pdg->ReadPDGTable("pdg_table.txt");
  
  TFile *f = new TFile(filename);
  if(!f) return;
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
  runHeaderMC->CStable().PrintAllCS();


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
  //TH2F* hisd   = new TH2F("hisd", "Dalitz's plot", 44, -0.65, 0.65, 44, -0.05, 1.05);
  TH2F* hiss   = new TH2F("hiss", "Phase Space IM{Sigma}", 500,1,2,300,0,1.5);
  TH2F* hisd   = new TH2F("hisd", "Phase Space IM{Sigma pi}", 500,1,2,300,0,1.5);
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
  bool isstate=false;
  for (Int_t i=0; i<nevent; i++) {
    tree->GetEvent(i);
    // print information
    int ndecay    = reactionData->NParticle(0);
    int nspec     = reactionData->NParticle(1);
    int nparticle = ndecay+nspec;
    if(!isstate){
      std::cout << "nparticle:" << nparticle << std::endl;
      isstate=true;
    }
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

    int flagCDH = 0;
    int flagNC = 0;
    int flagCVCPC = 0;
    for (Int_t j=0; j<detectorData->detectorHitSize(); j++) {
      int cid = detectorData->detectorHit(j)->detectorID();
      int pdgcode = detectorData->detectorHit(j)->pdg();
      if( pdgcode>10000000 ) continue;
      //if( pdgcode == 2112 ) cerr<<cid<<endl;
      //cerr<<cid<<" "<<pdgcode<<endl;
      if( cid == CID_CDH ){
	if( pdg->GetParticle(pdgcode)->Charge() ) flagCDH++;
      }
      else if( cid == CID_NC ){
	if( pdgcode==2112 ) flagNC++;
	//cerr<<pdgcode<<endl;
      }
      else if( cid == CID_CVC || cid == CID_PC ){
	if( pdg->GetParticle(pdgcode)->Charge() ) flagCVCPC++;
      }
    }

    //### CDH-1charged & NC-neutron ###//    
    double weight = 1;
    if( flagCDH && flagNC ){
      for (Int_t j=0; j<nspec; j++) {
	his00[j]->Fill(reactionData->TmpVal(j), weight);
      }
     for (Int_t j=0; j<nparticle; j++) {
       his10[j]->Fill(reactionData->GetCMParticle(j).P(), weight);
       his20[j]->Fill(reactionData->GetParticle(j).P(), weight);
       his30[j]->Fill(reactionData->GetCMParticle(j).CosTheta(), weight);
       his40[j]->Fill(reactionData->GetParticle(j).CosTheta(), weight);
     }
     his100->Fill(reactionData->ReactionID());
    }

    //### CDH-1charged & CVC/PC-charged ###//    
    if( flagCDH && flagCVCPC ){
      for (Int_t j=0; j<nspec; j++) {
	his01[j]->Fill(reactionData->TmpVal(j), weight);
      }
     for (Int_t j=0; j<nparticle; j++) {
       his11[j]->Fill(reactionData->GetCMParticle(j).P(), weight);
       his21[j]->Fill(reactionData->GetParticle(j).P(), weight);
       his31[j]->Fill(reactionData->GetCMParticle(j).CosTheta(), weight);
       his41[j]->Fill(reactionData->GetParticle(j).CosTheta(), weight);
     }
     his101->Fill(reactionData->ReactionID());
    }

#if DALITZ
    /* original 
    // Dalitz's plot
    double T[3];
    for(int j=0; j<nparticle; j++){
      double ene = reactionData->GetCMParticle(j).E();
      double mass = pdg->GetParticle(reactionData->PDG(j))->Mass();
      mass *= 1000;
      T[j] = ene-mass;
    }
    double Q = T[0]+T[1]+T[2];
    double x = (T[1]-T[0])/(sqrt(3)*Q);
    double y = T[2]/Q;
    hisd->Fill(x,y);
    */
    // Dalitz's plot
    double T[3];
    TLorentzVector TL_piSigma;
    TLorentzVector TL_nmiss;
    /*
    for(int j=0; j<nparticle; j++){
      double ene = reactionData->GetCMParticle(j).E();
      double mass = pdg->GetParticle(reactionData->PDG(j))->Mass();
      //TVector3 mom = reactionData->GetCMParticle(j).P();
      mass *= 1000;
    }*/
    //1: sigma ,2 :pion 3 neutron
    TLorentzVector TL_Sigma;
    TL_Sigma = reactionData->GetParticle(1);
    TL_piSigma = TL_Sigma+reactionData->GetParticle(2);
    TL_nmiss = reactionData->GetParticle(0);
    //std::cout << TL_nmiss.M() << std::endl;
    //std::cout << reactionData->GetCMParticle(1).M()<<std::endl;
    TLorentzVector TL_beam;
    TVector3 beammom(0,0,1000);
    TL_beam.SetVectM(beammom, 493);
    double q = (TL_beam.Vect()-TL_nmiss.Vect()).Mag();
    //TLorentzVector TL_sum = TL_gene[0]+TL_gene[1]+TL_gene[2];
    //double Q = T[0]+T[1]+T[2];
    //double x = (T[1]-T[0])/(sqrt(3)*Q);
    //double y = T[2]/Q;
    double mass = TL_piSigma.M();
    //std::cout << mass << " " << q << std::endl;
    hiss->Fill(TL_Sigma.M()/1000.,q/1000.);
    hisd->Fill(mass/1000.,q/1000.);
#endif
  }// for (Int_t i=0;i<nevent;i++) {
  cout<<"end of filling"<<endl;
  //----------------------//
  //--- event roop end ---//
  //----------------------//

  int Nall = his99->GetEntries();
  int Nnc  = his100->GetEntries();
  int Npc  = his101->GetEntries();
  double acc_nc = double(Nnc)/Nall;
  double acc_pc = double(Npc)/Nall;
  cout<<"acceptance:"<<endl;
  cout<<" CDH-1charged & NC-neutron     = "<<acc_nc<<endl;
  cout<<" CDH-1charged & CVC/PC-charged = "<<acc_pc<<endl;





  //--- plot ---//
  //gStyle->SetOptStat(111111);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TLatex *tex;

  //### plot reaction-ID
  TCanvas *c1 = new TCanvas("c1", "", 1200, 400);
  c1->Divide(3,1);
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
  tex = new TLatex(MINT+0.05*(MAXT-MINT), max*0.9, "CDH-1charged & NC-neutron");
  tex->SetTextColor(2);
  tex->SetTextSize(0.05);
  tex->Draw();

  c1->cd(3); his101->Draw("");
  his101->SetXTitle("reactionID [CERN-HERA-83-02]");
  his101->SetLineColor(4);
  int max = his101->GetMaximum();
  tex = new TLatex(MINT+0.05*(MAXT-MINT), max*0.9, "CDH-1charged & CVC/PC-charged");
  tex->SetTextColor(4);
  tex->SetTextSize(0.05);
  tex->Draw();
  //c1->Print("tmp1.pdf");
  //c1->Print("tmp1.png");

  //### plot Fermi Momentum
  FermiFlag = 1;
  if( FermiFlag ){
    TCanvas *c2 = new TCanvas("c2", "", 800, 400);
    c2->Divide(2,1);
    for( int i=0; i<nspec; i++ ){
      c2->cd(i+1); his0[i]->Draw("");
      his0[i]->SetXTitle("momentum [MeV/c]");
      int max = his0[i]->GetMaximum();
      tex = new TLatex(MIN+0.05*(MAXF-MIN), max*0.9, name[i+ndecay].c_str());
      tex->SetTextColor(4);
      tex->SetTextSize(0.1);
      tex->Draw();
      his00[i]->Draw("same"); his00[i]->SetLineColor(2);
      his01[i]->Draw("same"); his01[i]->SetLineColor(4);
    }
    //c2->Print("tmp2.pdf");
    //c2->Print("tmp2.png");
  }

  //### plot particle momentum @ CM
  TCanvas *c3 = new TCanvas("c3", "", 900, 600);
  c3->Divide(3,2);
  for( int i=0; i<nparticle; i++ ){
    c3->cd(i+1); his1[i]->Draw("");
    his1[i]->SetXTitle("momentum [MeV/c]");
    int max = his1[i]->GetMaximum();
    tex = new TLatex(MIN+0.05*(MAX-MIN), max*0.9, name[i].c_str());
    tex->SetTextColor(4);
    tex->SetTextSize(0.1);
    tex->Draw();
    //his10[i]->Draw("same"); his10[i]->SetLineColor(2);
    //his11[i]->Draw("same"); his11[i]->SetLineColor(4);
  }
  //c3->Print("tmp3.pdf");
  //c3->Print("tmp3.png");

  //### plot particle momentum @ LAB
  TCanvas *c4 = new TCanvas("c4", "", 900, 600);
  c4->Divide(3,2);
  for( int i=0; i<nparticle; i++ ){
    c4->cd(i+1); his2[i]->Draw(""); //gPad->SetLogy();
    his2[i]->SetXTitle("momentum [MeV/c]");
    int max = his2[i]->GetMaximum();
    tex = new TLatex(MIN+0.05*(MAX-MIN), max*0.9, name[i].c_str());
    tex->SetTextColor(4);
    tex->SetTextSize(0.1);
    tex->Draw();
    //his20[i]->Draw("same"); his20[i]->SetLineColor(2);
    //his21[i]->Draw("same"); his21[i]->SetLineColor(4);
  }
  //c4->Print("tmp4.pdf");
  //c4->Print("tmp4.png");

#if DALITZ
  // Dalitz's plot
  TLine* line1 = new TLine(0, 1, -1/sqrt(3), 0);
  line1->SetLineColor(4);
  TLine* line2 = new TLine(-1/sqrt(3), 0, 1/sqrt(3), 0);
  line2->SetLineColor(4);
  TLine* line3 = new TLine(1/sqrt(3), 0, 0, 1);
  line3->SetLineColor(4);

  TCanvas *c5 = new TCanvas("c5", "");
  hisd->SetXTitle("IM(#pi#Sigma) [GeV/c^{2}]");
  hisd->SetYTitle("mom. transfer [GeV/c]");
  hisd->GetXaxis()->CenterTitle();
  hisd->GetYaxis()->CenterTitle();
  hisd->Draw("colz");
  //hisd->SetXTitle("(T_{decay2}-T_{decay1})/#sqrt{3}Q");
  //hisd->SetYTitle("T_{deuteron}/Q");
  //line1->Draw(); line2->Draw(); line3->Draw();
  //c5->Print("tmp5.pdf");
  //c5->Print("tmp5.png");
#endif

  //### plot angular distribution @ CM
  TCanvas *c6 = new TCanvas("c6", "", 900, 600);
  c6->Divide(3,2);
  for( int i=0; i<nparticle; i++ ){
    c6->cd(i+1); his3[i]->Draw("");
    his3[i]->SetMinimum(0);
    his3[i]->SetXTitle("cos#theta");
    int max = his3[i]->GetMaximum();
    tex = new TLatex(MINT+0.05*(MAXT-MINT), max*0.9, name[i].c_str());
    tex->SetTextColor(4);
    tex->SetTextSize(0.1);
    tex->Draw();
    //his30[i]->Draw("same"); his30[i]->SetLineColor(2);
    //his31[i]->Draw("same"); his31[i]->SetLineColor(4);
  }
  //c6->Print("tmp6.pdf");
  //c6->Print("tmp6.png");

  //### plot angular distribution @ CM
  TCanvas *c7 = new TCanvas("c7", "", 900, 600);
  c7->Divide(3,2);
  for( int i=0; i<nparticle; i++ ){
    c7->cd(i+1); his4[i]->Draw("");
    his4[i]->SetMinimum(0);
    his4[i]->SetXTitle("cos#theta");
    int max = his4[i]->GetMaximum();
    tex = new TLatex(MINT+0.05*(MAXT-MINT), max*0.9, name[i].c_str());
    tex->SetTextColor(4);
    tex->SetTextSize(0.1);
    tex->Draw();
    //his40[i]->Draw("same"); his40[i]->SetLineColor(2);
    //his41[i]->Draw("same"); his41[i]->SetLineColor(4);
  }
  TFile *fout = new TFile(outname.c_str(),"RECREATE");
  fout->cd();
  hiss->Write();
  hisd->Write();
  fout->Close();
  return;
  //c7->Print("tmp7.pdf");
  //c7->Print("tmp7.png");
}
