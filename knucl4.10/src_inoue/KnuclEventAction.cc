// $Id: KnuclEventAction.cc,v 1.1.1.1 2013/12/25 01:25:10 sakuma Exp $
// GEANT4 tag $Name:  $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
#include "KnuclEventAction.hh"
#include "KnuclRunAction.hh"
#include "KnuclHit.hh"
#include "KnuclAnaManager.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4TrajectoryPoint.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4UnitsTable.hh"
#include "G4UImanager.hh"
#include "G4LorentzVector.hh"

#include "Randomize.hh"

//#include "KnuclEventData.h"
#include "KnuclCommon.h"

#include "stdlib.h"
#include <string.h>
#include <ctype.h>
#include <list>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
KnuclEventAction::KnuclEventAction(KnuclAnaManager* ana)
  : AnaManager(ana), counterCollID(-1), chamberCollID(-1)
{
  t0=time(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
KnuclEventAction::~KnuclEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void KnuclEventAction::BeginOfEventAction(const G4Event* evt)
{
  
  G4int evtNb = evt->GetEventID();
  //G4cout << "\n---> Begin of event: " << evtNb << G4endl;

#if 0
  if( evtNb<21000 ) {
    G4UImanager * UImanager = G4UImanager::GetUIpointer();
    UImanager->ApplyCommand("/event/abort");
  }
#endif

  if (evtNb<10) 
    G4cout << "\n---> Begin of event: " << evtNb << G4endl;
  else if (evtNb>=10 && evtNb<100) {
    int t1=time(0);
    if (evtNb%10 == 0) G4cout << "\n---> Begin of event: " << evtNb << "   time : "<<t1-t0<<G4endl;
  }
  else if (evtNb>=100 && evtNb<1000) {
    int t1=time(0);
    if (evtNb%100 == 0) G4cout << "\n---> Begin of event: " << evtNb << "   time : "<<t1-t0<<G4endl;
  }
  else if (evtNb>=1000) {
    int t1=time(0);
    if (evtNb%1000 == 0) G4cout << "\n---> Begin of event: " << evtNb << "   time : "<<t1-t0<<G4endl;
  }

  if (counterCollID == -1) {
    G4SDManager * SDman = G4SDManager::GetSDMpointer();
    counterCollID = SDman->GetCollectionID("counterSD");
    chamberCollID = SDman->GetCollectionID("chamberSD");
    G4cout << "counterCollID:" << counterCollID << G4endl;
    G4cout << "chamberCollID:" << chamberCollID << G4endl;
  }

  AnaManager->BeginOfEventAction();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void KnuclEventAction::EndOfEventAction(const G4Event* evt)
{
  G4int evtNb = evt->GetEventID();

  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  KnuclHitsCollection* CounterHC = 0;
  KnuclHitsCollection* ChamberHC = 0;
  if (HCE) CounterHC = (KnuclHitsCollection*)(HCE->GetHC(counterCollID));
  if (HCE) ChamberHC = (KnuclHitsCollection*)(HCE->GetHC(chamberCollID));

  G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();

  G4int counter_n_hit = 0;
  G4int chamber_n_hit = 0;
 
  // std::cout<<"debug flag = " <<AnaManager->GetDebugFlag()<<std::endl;
  // std::cout<<"debug flag = " <<AnaManager->GetComment()<<std::endl;
  // getchar();

  G4ThreeVector vtx = AnaManager->GetVertexPosition();

  //--------------------------//
  //--- class DetectorData ---//
  //--------------------------//
  DetectorData* detectorData = new DetectorData();
  std::vector <DetectorHit> detectorHit;
  Int_t nd = 0;
  // CounterHit
  std::list   <int> tracklist_det;
  if (CounterHC) {
    counter_n_hit = CounterHC->entries();
    for (G4int i=0;i<counter_n_hit;i++){
      DetectorHit tmp;
      if( (*CounterHC)[i]->GetEdepAbs()<AnaManager->GetADCThreshold((*CounterHC)[i]->GetDetectorID()) ) continue;
      tmp.setHitID(nd);
      tmp.setDetectorID((*CounterHC)[i]->GetDetectorID());
      //if( (*CounterHC)[i]->GetDetectorID() == CID_NC )  std::cerr<<"--- NC hit --- "<<(*CounterHC)[i]->GetPDGCode()<<" "<<(*CounterHC)[i]->GetEdepAbs()<<std::endl;
      tmp.setLayerID((*CounterHC)[i]->GetLayerID());
      tmp.setChannelID((*CounterHC)[i]->GetChannelID());
      tmp.setAdc((*CounterHC)[i]->GetEdepAbs());
      tmp.setTdc((*CounterHC)[i]->GetTime());
      tmp.setPos(ConvVecGT((*CounterHC)[i]->GetPos()));
      tmp.setMomentum(ConvVecGT((*CounterHC)[i]->GetMomentum()));
      tmp.setTrackID((*CounterHC)[i]->GetTrackID());
      tmp.setTime((*CounterHC)[i]->GetTime());
      tmp.setDe((*CounterHC)[i]->GetEdepAbs());
      tmp.setDt(0);
      tmp.setDx(0);
      tmp.setPDG((*CounterHC)[i]->GetPDGCode());
      tmp.setParentID((*CounterHC)[i]->GetParentID());
      detectorHit.push_back(tmp);
      nd++;

      if( (*CounterHC)[i]->GetDetectorID()==CID_NC || (*CounterHC)[i]->GetDetectorID()==CID_CVC ||  
	  (*CounterHC)[i]->GetDetectorID()==CID_PC || (*CounterHC)[i]->GetDetectorID()==CID_CDH ){
	tracklist_det.push_back((*CounterHC)[i]->GetTrackID());
	int id=(*CounterHC)[i]->GetParentID();
	while( id>0 ){
	  tracklist_det.push_back(id);
	  id=GetTrajectory(evt, id)->GetParentID();
	}
      }
    }
  }

  // ChamberHit
  if (ChamberHC) {
    chamber_n_hit = ChamberHC->entries();
    for (G4int i=0;i<chamber_n_hit;i++){
      DetectorHit tmp;
      tmp.setHitID(nd);
      tmp.setDetectorID((*ChamberHC)[i]->GetDetectorID());
      tmp.setLayerID((*ChamberHC)[i]->GetLayerID());
      tmp.setChannelID((*ChamberHC)[i]->GetChannelID());
      tmp.setAdc((*ChamberHC)[i]->GetEdepAbs());
      tmp.setTdc((*ChamberHC)[i]->GetTime());
      tmp.setPos(ConvVecGT((*ChamberHC)[i]->GetPos()));
      tmp.setMomentum(ConvVecGT((*ChamberHC)[i]->GetMomentum()));
      tmp.setTrackID((*ChamberHC)[i]->GetTrackID());
      tmp.setTime((*ChamberHC)[i]->GetTime());
      tmp.setDe((*ChamberHC)[i]->GetEdepAbs());
      tmp.setDt(0);
      tmp.setDx((*ChamberHC)[i]->GetDx());
      tmp.setPDG((*ChamberHC)[i]->GetPDGCode());
      tmp.setParentID((*ChamberHC)[i]->GetParentID());
      detectorHit.push_back(tmp);
      nd++;
    }
  }

  bool trigger_flag=trigger(detectorHit, AnaManager->GetTrigger());
  //  std::cout<<std::boolalpha<<trigger_flag<<std::endl;
  if( trigger_flag ) detectorData->setDetectorHit(detectorHit);

  AnaManager->setHistDetectorData(detectorData);
  //std::cout<<counter_n_hit<<" "<<chamber_n_hit<<std::endl;

  //-- select NC/PC fired events ---//
  bool NCflag = false;
  bool CVCflag = false;
  bool PCflag = false;
  for (int j=0; j<detectorData->detectorHitSize(); j++) {
    if      (detectorData->detectorHit(j)->detectorID() == CID_NC ) NCflag = true;
    else if (detectorData->detectorHit(j)->detectorID() == CID_PC ) PCflag = true;
    else if (detectorData->detectorHit(j)->detectorID() == CID_CVC) CVCflag = true;
  }
  //-- select NC/PC fired events ---//

  std::list   <int> tracklist;
  std::list   <int> tracklist_trj;
 
  std::list   <int> parentID;
  std::list   <int> daughter;
  std::list   <int> granddaug;
  std::list   <int> ggranddaug;

  std::list   <int>::iterator it_find;

  //find primarys ( i.e. pID==0 trajectries )
  for (int i=0; i<n_trajectories; i++) {
     KnuclTrajectory* trj = (KnuclTrajectory*)((*(evt->GetTrajectoryContainer()))[i]);
     G4int tID    = trj->GetTrackID(); G4int PartID = trj->GetParentID(); 
     if ( PartID==0 )                                 {
       parentID .push_back(tID); 
     }
  }

  //find secondary( i.e. pID == parentID trajectries )
  for (int i=0; i<n_trajectories; i++) {
     KnuclTrajectory* trj = (KnuclTrajectory*)((*(evt->GetTrajectoryContainer()))[i]);
     G4int tID    = trj->GetTrackID(); G4int PartID = trj->GetParentID(); G4int PID    = trj->GetPDGEncoding();
     for ( it_find = parentID.begin(); it_find != parentID.end(); ++it_find ) {
       if ( PartID == *it_find ) {
         if ( PID!=11 && abs(PID)!=22 )
           daughter.push_back(tID); 
       }
     }
  }

  //find secondary( i.e. pID == daugtherID trajectries )
  for (int i=0; i<n_trajectories; i++) {
     KnuclTrajectory* trj = (KnuclTrajectory*)((*(evt->GetTrajectoryContainer()))[i]);
     G4int tID    = trj->GetTrackID(); G4int PartID = trj->GetParentID(); G4int PID    = trj->GetPDGEncoding();
     for ( it_find = daughter.begin(); it_find != daughter.end(); ++it_find ) {
       if ( PartID == *it_find ) {
         if ( PID!=11 && abs(PID)!=22 )
           granddaug.push_back(tID); 
       }
     }
  }

  //find secondary( i.e. pID == granddaugID trajectries )
  for (int i=0; i<n_trajectories; i++) {
     KnuclTrajectory* trj = (KnuclTrajectory*)((*(evt->GetTrajectoryContainer()))[i]);
     G4int tID    = trj->GetTrackID(); G4int PartID = trj->GetParentID(); G4int PID    = trj->GetPDGEncoding();
     for ( it_find = granddaug.begin(); it_find != granddaug.end(); ++it_find ) {
       if ( PartID == *it_find ) {
         if ( PID!=11 && abs(PID)!=22 )
           ggranddaug.push_back(tID); 
       }
     }
  }

  std::vector <Track> track;

  int FillLevel = AnaManager->GetFillLevel();
  std::list  <int>::iterator it_v;
  //std::cerr<<tracklist.size()<<" "<<parentID.size()<<" "<<daughter.size()<<" "<<granddaug.size()<<" "<<ggranddaug.size()<<std::endl;

  if( FillLevel<5 ){
    // add parent  
    if( FillLevel >= 1 ){
      for ( it_v = parentID.begin(); it_v != parentID.end(); ++it_v ) {
	tracklist.push_back(*it_v);
      }
    }
    // add daugther
    if( FillLevel >= 2 ){
      for ( it_v = daughter.begin(); it_v != daughter.end(); ++it_v ) {
	tracklist.push_back(*it_v);
      }
    }
    // add grand-daugther
    if( FillLevel >= 3 ){
      for ( it_v = granddaug.begin(); it_v != granddaug.end(); ++it_v ) {
	tracklist.push_back(*it_v);
      }
    }
    // add grand-grand-daugther
    if( FillLevel >= 4 ){
      for ( it_v = ggranddaug.begin(); it_v != ggranddaug.end(); ++it_v ) {
	tracklist.push_back(*it_v);
      }
    }
  }
  // user selection
  else if( FillLevel==5 ){
    std::vector<KnuclTrajectory*> parent;
    for (int i=0; i<n_trajectories; i++) {
      KnuclTrajectory* trj = (KnuclTrajectory*)((*(evt->GetTrajectoryContainer()))[i]);
      G4int tID    = trj->GetTrackID(); G4int PartID = trj->GetParentID();
      if ( PartID==0 ){
        //      std::cout<<"> Primary Particle  "<<std::setw(10)<<std::left<<trj->GetParticleName()<<" tID : "<<trj->GetTrackID()<<std::endl; 
        tracklist.push_back(tID);
        if( IsFill(trj) ) parent.push_back(trj);
      }
    }
    while( true ){
      std::vector<KnuclTrajectory*> daug;
      for( int i=0; i<n_trajectories; i++ ){
        KnuclTrajectory* trj = (KnuclTrajectory*)((*(evt->GetTrajectoryContainer()))[i]);
        for( int j=0; j<parent.size(); j++ ){
          int level = IsFill(parent[j], trj);
          if( level==1 ) tracklist.push_back(trj->GetTrackID());
          if( level==2 ){
            tracklist.push_back(trj->GetTrackID());
            daug.push_back(trj);
          }
        }
      }
      if( daug.empty() ) break;
      parent = daug;
    }
  }
  for ( std::list<int>::iterator it_l = tracklist_det.begin(); it_l != tracklist_det.end(); ++it_l ) tracklist.push_back(*it_l);

  tracklist.sort(); tracklist.unique();
  //std::cerr<<tracklist.size()<<std::endl;

  std::list<int>::iterator it_l;
  for ( std::list<int>::iterator it_l = tracklist.begin(); it_l != tracklist.end(); ++it_l ) {
    for (G4int j=0; j<n_trajectories; j++) {
      KnuclTrajectory* trj = (KnuclTrajectory*)((*(evt->GetTrajectoryContainer()))[j]);
      G4int tmp_tID = trj->GetTrackID();
      if (tmp_tID == *it_l ) tracklist_trj.push_back(j);  
    }
  }
    
  tracklist_trj.sort();  tracklist_trj.unique();

  for ( it_l = tracklist_trj.begin(); it_l != tracklist_trj.end(); ++it_l ) {
    KnuclTrajectory* trj = (KnuclTrajectory*)((*(evt->GetTrajectoryContainer()))[*it_l]);
    G4int tID = trj->GetTrackID();
    G4int pID = trj->GetParentID();
    int reacID = GetProcessID(trj->GetProcessName());
    
    Track tmp;
    tmp.setTrackID(tID);
    tmp.setProcessID(reacID);
    tmp.setParentTrackID(pID);
    tmp.setPdgID(trj->GetPDGEncoding());
    tmp.setVertex(ConvVecGT(trj->GetVertexPosition()));
    tmp.setMomentum(ConvVecGT(trj->GetInitialMomentum()));
    for (int x=0; x<(G4int)detectorHit.size(); x++){
      if( detectorHit[x].trackID() == tID )
        tmp.setDetectorHitLink(detectorHit[x].hitID());
    }
    //
    // search CDH hits
    //
    // Flight Length calc.
    //
    bool fl_flag=false;
    double   FlightLength = 0.0;
    TVector3 prev_pos;
    TVector3 curr_pos;
    std::vector <int> parent_list;
    
    for (int ii=0; ii<trj->GetPointEntries()-1; ii++) {
      for (int ih=0; ih<detectorHit.size(); ih++ ){
	double hit_x = detectorHit[ih].pos().X();
	double hit_y = detectorHit[ih].pos().Y();
	double delx = fabs (hit_x-(trj->GetPoint(ii)->GetPosition()).x());
	double dely = fabs (hit_y-(trj->GetPoint(ii)->GetPosition()).y());
	if (delx==0 && dely==0 ) {
	  if( detectorHit[ih].detectorID()==CID_CDH || detectorHit[ih].detectorID()==CID_CVC || 
	      detectorHit[ih].detectorID()==CID_PC || detectorHit[ih].detectorID()==CID_NC ){
	    fl_flag=true;
	    break;
	  }
	}
      }
      if( fl_flag ){
	break;
      }
      FlightLength+=(trj->GetPoint(ii)->GetPosition()-trj->GetPoint(ii+1)->GetPosition()).mag();
    }
    
    if( fl_flag ) tmp.setFlightLength( FlightLength );
    track.push_back(tmp);
    //printf("%5d %5d %6d %8.3f\n",trj->GetTrackID(),trj->GetParentID(),trj->GetPDGEncoding(),FlightLength);
  }
  
  if(!strcmp(AnaManager->GetDebugFlag(),"true")){
    //strcmp: if two char * are the same, return 0!

    // std::cout<<"debug flag = " <<AnaManager->GetDebugFlag()<<std::endl;
    // std::cout<<"debug flag = " <<AnaManager->GetComment()<<std::endl;
    // getchar();
  }


  //-----------------//
  //--- fill tree ---//
  //-----------------//
  EventHeaderMC* eventHeaderMC = new EventHeaderMC();
  eventHeaderMC->setEventID(evtNb);
  AnaManager->setHistEventHeaderMC(eventHeaderMC);

  MCData* mcData = new MCData();
  if( trigger_flag ){
    mcData->setTrack(track);
    AnaManager->setHistMCData(mcData);
  }
  else AnaManager->setHistMCData(mcData);

  if( !AnaManager->GetFWDTrigger() )
    AnaManager->fillHistTree();
  else if( AnaManager->GetFWDTrigger() && (NCflag || PCflag) )
    AnaManager->fillHistTree();

  AnaManager->EndOfEventAction();
  //AnaManager->showHistTree();
  
  delete detectorData;
  delete eventHeaderMC;
  delete mcData;
}

KnuclTrajectory* KnuclEventAction::GetTrajectory(const G4Event *evt, int id){
  for (int i=0; i<evt->GetTrajectoryContainer()->entries(); i++) {
    KnuclTrajectory* trj = (KnuclTrajectory*)((*(evt->GetTrajectoryContainer()))[i]);
    if( trj->GetTrackID()==id ) return trj;
  }

  throw std::runtime_error("KnuclTrajectory not found");
}

G4int KnuclEventAction::GetProcessID(const G4String &name)
{
  if( name=="" ) return 0;
  else if( name=="Decay"                 ) return 1;
  else if( name=="conv"                  ) return 2;
  else if( name=="Transportation"        ) return 3;
  else if( name=="phot"                  ) return 4;
  else if( name=="annihil"               ) return 5;
  else if( name=="compt"                 ) return 6;
  else if( name=="eBrem"                 ) return 7;
  else if( name=="hadElastic"            ) return 8;
  else if( name=="CoulombScat"           ) return 9;
  else if( name=="nKiller"               ) return 10;
  else if( name=="photonNuclear"         ) return 11;
  else if( name=="msc"                   ) return 12;
  else if( name=="pi-Inelastic"          ) return 100;
  else if( name=="pi+Inelastic"          ) return 101;
  else if( name=="kaon-Inelastic"        ) return 102;
  else if( name=="kaon+Inelastic"        ) return 103;
  else if( name=="kaon0LInelastic"       ) return 104;
  else if( name=="kaon0SInelastic"       ) return 105;
  else if( name=="lambdaInelastic"       ) return 106;
  else if( name=="sigma+Inelastic"       ) return 107;
  else if( name=="sigma-Inelastic"       ) return 108;
  else if( name=="sigma0Inelastic"       ) return 109;
  else if( name=="protonInelastic"       ) return 110;
  else if( name=="neutronInelastic"      ) return 111;
  else if( name=="dInelastic"            ) return 112;
  else if( name=="tInelastic"            ) return 113;
  else if( name=="He3Inelastic"          ) return 114;
  else if( name=="alphaInelastic"        ) return 115;
  else if( name.find("Inelastic")!=std::string::npos ) return 199;
  else if( name=="eIoni"                 ) return 200;
  else if( name=="hIoni"                 ) return 201;
  else if( name=="ionIoni"               ) return 202;
  else if( name=="muIoni"                ) return 203;
  else if( name=="hBertiniCaptureAtRest" ) return 204;
  else if( name=="nCapture"              ) return 205;
  else if( name=="muMinusCaptureAtRest"  ) return 206;
  else if( name=="unknown"          ) return -999;
  else return -1;
}

bool KnuclEventAction::IsFill(KnuclTrajectory *trj)
{
  G4ParticleDefinition *def = trj-> GetParticleDefinition();
  std::string pname = def-> GetParticleName();
  if( pname.find("sigma")!=std::string::npos || pname.find("lambda")!=std::string::npos || pname.find("kaon")!=std::string::npos ||
      pname.find("rho")!=std::string::npos || pname.find("eta")!=std::string::npos || pname.find("star")!=std::string::npos ||
      pname.find("delta")!=std::string::npos || pname=="D+" || pname=="D0" || pname=="D-" ){
    return true;
  }
  return false;
}

int KnuclEventAction::IsFill(KnuclTrajectory *parent, KnuclTrajectory *daughter)
{
  if( !parent ) std::cout<<" !!! parent error !!!"<<std::endl;
  if( !daughter ) std::cout<<" !!! daughter error !!!"<<std::endl;

  if( parent->GetTrackID()!=daughter->GetParentID() ) return 0;
  G4ParticleDefinition *def_p = parent-> GetParticleDefinition();
  G4ParticleDefinition *def_d = daughter-> GetParticleDefinition();
  std::string pname = def_p-> GetParticleName();
  std::string dname = def_d-> GetParticleName();

  if( daughter-> GetPDGEncoding()==11 ) return 0;

  std::string process = daughter->GetProcessName();
  if( daughter-> GetPDGEncoding()==22 && process.find("Decay")==std::string::npos ) return 0;
  //  std::cout<<"> Secondary  : "<<dname<<"  tID: "<<daughter->GetTrackID()<<"  Process : "<<process<<"  parent  : "<<pname<<std::endl;
  if( dname.find("sigma")!=std::string::npos || dname.find("lambda")!=std::string::npos || dname.find("kaon")!=std::string::npos ||
      dname.find("rho")!=std::string::npos || dname.find("eta")!=std::string::npos || dname.find("star")!=std::string::npos ||
      dname.find("delta")!=std::string::npos || pname=="D+" || pname=="D0" || pname=="D-" ) return 2;
  return 1;
}
