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

Track* FindTrackFromMcIndex(MCData *mcdata, int trackid)
{
  for( int itr=0;itr<mcdata->trackSize();itr++){
    Track *track=mcdata->track(itr);
    if( track->trackID()==trackid) return track;
  }
  return 0;
}

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
   
  TH2F *IMnpim_IMnpip1 = new TH2F("IMnpim_IMnpip1","IMnpim_IMnpip1",2000,0,2,2000,0,2); 
  TH2F *IMnpim_IMnpip2 = new TH2F("IMnpim_IMnpip2","IMnpim_IMnpip2",2000,0,2,2000,0,2); 
  TH1F *pimom = new TH1F("pimom","pimom",200,-1,1); 
  TH2F *MMnmiss_nmom = new TH2F("MMnmiss_nmom","MMnmiss_nmom",100,0,1.0, 100,0.4,1.9);
  TH2F *MMnmiss_nmom_det = new TH2F("MMnmiss_nmom_det","MMnmiss_nmom_det",100,0,1.0, 100,0.4,1.9);
  TH2F *MMnmiss_nmom_comb = new TH2F("MMnmiss_nmom_comb","MMnmiss_nmom_comb",100,0,1.0, 100,0.4,1.9);


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
    int EventType=0;
    struct pimInfo{
      int gen;
      int processID;
      int parentProcessID;
      int seg;
      double dE;
      double mom;
      TVector3 mom_vec;
      double time;
      pimInfo(){
        gen=-1;
        processID=-1;
        parentProcessID=-1;
        dE=0.0;
        seg=-1;
        mom=0.0;
        time=0.0;
        mom_vec.SetXYZ(-999,-999,-999);
      }
    };

    struct pipInfo{
      int gen;
      int processID;
      int parentProcessID;
      int seg;
      double dE;
      double mom;
      TVector3 mom_vec;
      double time;
      pipInfo(){
        gen=-1;
        processID=-1;
        parentProcessID=-1;
        dE=0.0;
        seg=-1;
        mom=0.0;
        time=0.0;
        mom_vec.SetXYZ(-999,-999,-999);
      }
    };

    struct NcanInfo{
      int pdg;
      int parentpdg;
      int gen;
      int processID;
      int parentProcessID;
      TVector3 mom_vec;
      int seg;
      double dE;
      double mom;
      double vtx_r;
      double vtx_r_g1parent;
      double vtx_r_g2parent;
      DetectorHit *dhitncan;
      double time;
      int ntype;
      NcanInfo(){
        pdg=-9999;
        parentpdg=-9999;
        gen=-1;
        processID=-1;
        parentProcessID=-1;
        mom_vec.SetXYZ(-999.,-999.,-999.);
        seg=-1;
        dE=0.0;
        mom=0.0;
        vtx_r=0.0;
        vtx_r_g1parent=0.0;
        vtx_r_g2parent=0.0;
        dhitncan=NULL;
        time=0.0;
        ntype=-1;
      }
    };
  
    pimInfo piminfo;
    pipInfo pipinfo;
    NcanInfo ncaninfo;

    //** print event information **//
    int ndecay    = reactionData->NParticle(0);
    int nspec     = reactionData->NParticle(1);
    int nparticle = ndecay+nspec;
    //std::cout << "nparticle " << nparticle << std::endl;
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
    // beam_K(K+), prompt pi+, prompt pi-,  prompt n, Lambda(missing), p from L decay, pi- from L decay
    // beam_K(K+), prompt pi+, prompt pi-,  prompt n, Lambda(missing), n from L decay, pi0 from L decay
    int reactionID = reactionData->ReactionID();
    const int npart=7;
    //int PDG[npart] = {321, 211, -211, 2112, 3122,-211,2212};//br0, L->pi-
    int PDG[npart] = {321, 211, -211, 2112, 3122,2112,111};//br1, L->npi0
    int parentID[npart] = {0,0,0,0,0,0,0};
    int ID[npart]       = {-1,-1,-1,-1,-1,-1};
    int trackID[npart]  = {-1,-1,-1,-1,-1,-1};
    nparticle = 0;
    for (Int_t itrack=0; itrack<mcData->trackSize(); itrack++) {
      int pdgcode = mcData->track(itrack)->pdgID();
      int parent  = mcData->track(itrack)->parentTrackID();
      int track   = mcData->track(itrack)->trackID();
      for( int ip=0; ip<npart; ip++ ){
        if( pdgcode==PDG[ip] && parent==parentID[ip] && ID[ip]==-1 ){
       // if( pdgcode==PDG[ip] && parent==parentID[ip] ){
          ID[ip] = itrack;
          trackID[ip] = track;
          nparticle++;
          //std::cout << "gen ip " << ip << std::endl; 
          if(ip==4){ //if parent is Lambda
            parentID[5]=track;
            parentID[6]=track;
          }
          break;
        }
      }
    } // for itrack

    if(nparticle!=6){
      //std::cout << "nparticle " << nparticle << std::endl; 
      //continue;
    }
    //bool flagG4Decay = (nparticle==(npart-1)) ? true : false;
    int nCDHhit[npart] = {0, 0, 0, 0, 0,0,0};
    //Since CDH hit is made by proton, heavy ions, generated by the neutron
    int CDHfired = 0;
    int npim=0;
    int npip=0;
    for (Int_t j=0; j<detectorData->detectorHitSize(); j++) {
      DetectorHit *dhit = detectorData->detectorHit(j);
      int cid     = dhit->detectorID();
      if(cid!=CID_CDH){
        continue;
      }
      int pdg    = dhit->pdg();
      int track   = dhit->trackID();
      int seg  = dhit->channelID();//0 origin
      double dE   = dhit->adc();
      Track *track_p  = FindTrackFromMcIndex(mcData,track);
      for( int ip=1; ip<npart; ip++ ){
        if( cid==CID_CDH && track==trackID[ip]){
          nCDHhit[ip]++;
        }
      }//ip
      double truemom = (track_p->momentum()).Mag()/1000.;
      TVector3 momVec = track_p->momentum();
      if(pdg==-211){
         npim++;
        // int parentpdg = 0;
        // int parentId = dhit->parentID();
        // Track *parenttrack_p = FindTrackFromMcIndex(mcData,parentId);

         piminfo.dE = dE;
         piminfo.seg = seg;
         piminfo.mom = truemom;
         piminfo.mom_vec = momVec;
      }
      else if(pdg==211){
         npip++;
         pipinfo.dE = dE;
         pipinfo.seg = seg;
         pipinfo.mom = truemom;
         pipinfo.mom_vec = momVec;
      }else{
         int parentpdg = 0;
         int parentId = dhit->parentID();
         Track *parenttrack_p = FindTrackFromMcIndex(mcData,parentId);
         if(parenttrack_p !=0){
           parentpdg= parenttrack_p->pdgID();
         }
        
         while(parenttrack_p!=0){
           if(parenttrack_p->trackID()==trackID[3]) ncaninfo.ntype=1;
           if(parenttrack_p->trackID()==trackID[5]) ncaninfo.ntype=2;
           ncaninfo.mom_vec = parenttrack_p->momentum();
           parenttrack_p=FindTrackFromMcIndex(mcData,parenttrack_p->parentTrackID());
         }

         ncaninfo.pdg = pdg;
         ncaninfo.seg = seg;
         ncaninfo.parentpdg = parentpdg;
         ncaninfo.mom = truemom;
      }

      if(dE>2.0)CDHfired++;
    }
                              //prompt  pi+, prompt pi-,   prompt n        
    bool pippimn_detect1 =  (nCDHhit[1] && nCDHhit[2] && nCDHhit[3]) ? true : false;
    //if(pippimn_detect1) std::cout << "yes 1" << std::endl;
                              //prompt pi+, prompt pi-,   decay n
    bool pippimn_detect2 =  (nCDHhit[1] && nCDHhit[2] && nCDHhit[5]) ? true : false;
    //if(pippimn_detect2) std::cout << "yes 2" << std::endl;

    bool pippim_detect1 = (nCDHhit[1] && nCDHhit[2]) ? true : false;
    bool pippim_detect2 = (nCDHhit[1] && nCDHhit[2]) ? true : false;
    if(!(npim==1 && npip==1 && abs(ncaninfo.seg - pipinfo.seg)!=1 && abs(ncaninfo.seg - pipinfo.seg)!=35 
                            && abs(ncaninfo.seg - piminfo.seg)!=1 && abs(ncaninfo.seg - piminfo.seg)!=35 
                            && CDHfired==3 )) continue;


    //if(CDHfired==3 &&  pippim_detect1) std::cout << "yes" << std::endl;
    //=== acceptance study ===//
    TLorentzVector TL_mcData[npart];
    for(int ip=0;ip<npart;ip++){
      if(ip) TL_mcData[ip].SetVectM(mcData->track(ID[ip])->momentum(),pdg->GetParticle(PDG[ip])->Mass()*1000.0);
      else  TL_mcData[ip].SetVectM(mcData->track(ID[ip])->momentum()*-1.0,pdg->GetParticle(PDG[ip])->Mass()*1000.0);
      //std::cout << "IP " << ip << " Mass " << TL_mcData[ip].M() << "  Mom " << TL_mcData[ip].P() <<  std::endl;
    }
    
    TLorentzVector TL_detData[4];
    TL_detData[1].SetVectM(piminfo.mom_vec,pdg->GetParticle(PDG[1])->Mass()*1000.0);
    TL_detData[2].SetVectM(pipinfo.mom_vec,pdg->GetParticle(PDG[2])->Mass()*1000.0);
    TL_detData[3].SetVectM(ncaninfo.mom_vec,pdg->GetParticle(PDG[3])->Mass()*1000.0);

    TLorentzVector TL_mcData_npim;
    //if(pippimn_detect1) TL_mcData_npim = TL_mcData[2] + TL_mcData[3];
    //else if(pippimn_detect2) TL_mcData_npim = TL_mcData[6] + TL_mcData[3];
    if(pippim_detect1) TL_mcData_npim = TL_mcData[2] + TL_mcData[5];
    //TL_mcData_npim = TL_mcData[2] + TL_mcData[5];
    //else if(pippim_detect2) TL_mcData_npim = TL_mcData[6] + TL_mcData[3];
   
    TLorentzVector TL_mcData_npip;
    //if(pippimn_detect1 || pippimn_detect2) TL_mcData_npip = TL_mcData[1] + TL_mcData[2];
    //if(pippim_detect1 || pippim_detect2) TL_mcData_npip = TL_mcData[1] + TL_mcData[5];
    TL_mcData_npip = TL_mcData[1] + TL_mcData[5];
     
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
    pimom->Fill(TL_mcData[1].P()/1000.0);
    pimom->Fill(-1.*TL_mcData[2].P()/1000.0);
    if(cos_n<0.7)IMnpim_IMnpip1->Fill(TL_mcData_npim.M()/1000.0,TL_mcData_npip.M()/1000.);
    //TLorentzVector mom_piS  = mom_pi_CM+mom_Lambda_CM;
    //TLorentzVector mom_piSp = mom_pi_CM+mom_Lambda_CM+mom_p_CM;
    // change to BG study of piSpn
    TLorentzVector mom_piS  = mom_pi_CM+mom_pi2_CM+mom_n_CM;
    TLorentzVector mom_piSp = mom_pi_CM+mom_pi2_CM+mom_n_CM+mom_p_CM;
    
    TLorentzVector mom_miss = mom_target+mom_beam-TL_mcData[1]-TL_mcData[2]-TL_mcData[3];
  //  TLorentzVector mom_miss_det = mom_target+mom_beam-TL_detData[1]-TL_detData[2]-TL_detData[3];
    TLorentzVector mom_miss_det = mom_target+mom_beam-TL_mcData[1]-TL_mcData[2]-TL_detData[3];
    TLorentzVector mom_miss_comb = mom_target+mom_beam-TL_mcData[1]-TL_mcData[2]-TL_mcData[5];

    Cosn_IMnpipi[2]->Fill(mom_piS.M()/1000, cos_n);
    Cosn_IMnppipi[2]->Fill(mom_piSp.M()/1000, cos_n);
    
    Qkn_IMnpipi[2]->Fill(mom_piS.M()/1000, qkn.P()/1000.0);
    Qkn_IMnppipi[2]->Fill(mom_piSp.M()/1000, qkn.P()/1000.0);

    IMnpipi_IMnppipi[2]->Fill(mom_piSp.M()/1000, mom_piS.M()/1000.0);
    
    if(CDHfired==3 && pippim_detect1 && ncaninfo.ntype==1){
      MMnmiss_nmom->Fill(TL_mcData[3].M()/1000.0,mom_miss.M()/1000.);
      MMnmiss_nmom_det->Fill(TL_mcData[3].M()/1000.0,mom_miss_det.M()/1000.);
    }
    else if(CDHfired==3 && pippim_detect1 && ncaninfo.ntype==2){
      MMnmiss_nmom_comb->Fill(TL_mcData[3].M()/1000.0,mom_miss_det.M()/1000.);
    }
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
    if (pippimn_detect1) {
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
    } // if (pippimn_detect) {
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
   
  TCanvas *ctest = new TCanvas("ctest","ctest");
  IMnpim_IMnpip1->Draw("colz");

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
    Cosn_IMnpipi_acc[x]->Write();
    Cosn_IMnppipi_acc[x]->Write();
    Qkn_IMnpipi_acc[x]->Write();
    Qkn_IMnppipi_acc[x]->Write();
    MMnmiss_nmom->Write();
    MMnmiss_nmom_det->Write();
    MMnmiss_nmom_comb->Write();
  }
  pimom->Write();
  IMnpim_IMnpip1->Write();
  out->Close();

}

