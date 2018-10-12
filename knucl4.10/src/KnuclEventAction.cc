// $Id: KnuclEventAction.cc,v 1.6 2016/11/09 08:41:09 inoue Exp $
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
    if (evtNb%10 == 0) G4cout << "\n---> Begin of event: " << evtNb << G4endl;
  }
  else if (evtNb>=100 && evtNb<1000) {
    if (evtNb%100 == 0) G4cout << "\n---> Begin of event: " << evtNb << G4endl;
  }
  else if (evtNb>=1000) {
    if (evtNb%1000 == 0) G4cout << "\n---> Begin of event: " << evtNb << G4endl;
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
  if( trigger_flag ) detectorData->setDetectorHit(detectorHit);

  AnaManager->setHistDetectorData(detectorData);
  //std::cout<<counter_n_hit<<" "<<chamber_n_hit<<std::endl;

  //-- select NC/PC fired events ---//
  bool NCflag = false;
  bool PCflag = false;
  for (int j=0; j<detectorData->detectorHitSize(); j++) {
    if      (detectorData->detectorHit(j)->detectorID() == CID_NC) NCflag = true;
    else if (detectorData->detectorHit(j)->detectorID() == CID_PC) PCflag = true;
  }
  //-- select NC/PC fired events ---//



  //-------------------------//
  //--- track information ---//
  //-------------------------//
  //
  // *** d_cdh[] is a hit container for CDH/CVC/PC ***
  //
  CounterHit d_cdh[N_CDH+N_CVC+N_PC];
  for (int j=0; j<detectorData->detectorHitSize(); j++) {
    int hit_id = -999;
    switch ( detectorData->detectorHit(j)->detectorID() ){
    case CID_CDH:
      hit_id = detectorData->detectorHit(j)->channelID();
      break;
    case CID_CVC:
      hit_id = N_CDH + detectorData->detectorHit(j)->channelID();
      break;
    case CID_PC:
      hit_id = N_CDH + N_CVC + detectorData->detectorHit(j)->channelID();
      break;
    default:
      hit_id = -999;
      break;
    }
    if ( hit_id>0 ){
      d_cdh[hit_id].Append(detectorData->detectorHit(j)->adc(),
                           detectorData->detectorHit(j)->tdc(),
                           detectorData->detectorHit(j)->pdg(),
                           detectorData->detectorHit(j)->trackID(),
                           detectorData->detectorHit(j)->pos() );
    }
  }

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

  // Creating tracklist which contained tracks
  // which creates hits on CDH, CDC PC and NC.
  // 
  //
  if (CounterHC) {
    counter_n_hit = CounterHC->entries();
    for (G4int i=0;i<counter_n_hit;i++){
      if ( 
           (*CounterHC)[i]->GetDetectorID() == CID_CDH || 
           (*CounterHC)[i]->GetDetectorID() == CID_CVC || 
           (*CounterHC)[i]->GetDetectorID() == CID_PC  ||
           (*CounterHC)[i]->GetDetectorID() == CID_NC   
         ) {
        G4int tID = (*CounterHC)[i]->GetTrackID();
        G4int pID = (*CounterHC)[i]->GetParentID();
        tracklist.push_back(tID);
        while (pID!=0){
           for (G4int j=0; j<n_trajectories; j++) {
              KnuclTrajectory* trj = (KnuclTrajectory*)((*(evt->GetTrajectoryContainer()))[j]);
              G4int tmp_tID = trj->GetTrackID();
              G4int tmp_pID = trj->GetParentID();
              if ( tmp_tID == pID ) {
                 tracklist.push_back(tmp_tID);
                 pID = tmp_pID; j=0;
              }
           }
        }
      }
    }
  }

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

  tracklist.sort(); tracklist.unique();
  //std::cerr<<tracklist.size()<<std::endl;

  std::list<int>::iterator it_l;
  for ( it_l = tracklist.begin(); it_l != tracklist.end(); ++it_l ) {
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

    Track tmp;
    tmp.setTrackID(tID);
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
     int hit_CDH_id = -999;
     for (int ii=0; ii<trj->GetPointEntries(); ii++) {
       for (int ih=0; ih<N_CDH; ih++){
         double hit_x = d_cdh[ih].getX();
         double hit_y = d_cdh[ih].getY();
         double delx = fabs (hit_x-(trj->GetPoint(ii)->GetPosition()).x());
         double dely = fabs (hit_y-(trj->GetPoint(ii)->GetPosition()).y());
         if (delx==0 && dely==0 ) {
           hit_CDH_id = ih;
         }
       }
     }
     // Flight Length calc.
     //
     double   FlightLength = -999.0;
     TVector3 prev_pos;
     TVector3 curr_pos;
     std::vector <int> parent_list;

     if (hit_CDH_id!=-999) {
        FlightLength = 0.0;
        for (int iii=0; iii<trj->GetPointEntries(); iii++) {
           double xx = (trj->GetPoint(iii)->GetPosition()).x();
           double yy = (trj->GetPoint(iii)->GetPosition()).y();
           double zz = (trj->GetPoint(iii)->GetPosition()).z();
           if (iii==0) {
              prev_pos.SetXYZ(xx,yy,zz);
              if (xx!=vtx.x() && yy!=vtx.y() && zz!=vtx.z() ) {
                parent_list.push_back(pID);
              }
           } else {
              curr_pos.SetXYZ(xx,yy,zz);
              TVector3 delta = curr_pos - prev_pos;
              double delx = fabs(d_cdh[hit_CDH_id].getX() -  curr_pos.x());
              double dely = fabs(d_cdh[hit_CDH_id].getY() -  curr_pos.y());
              FlightLength += delta.Mag();
              if (delx==0.0 && dely==0.0) break;
              prev_pos.SetXYZ(xx,yy,zz);
           }
        }
     }
     if ( parent_list.size()!=0 ) {
       for (int ii=0; ii<n_trajectories; ii++) {
         KnuclTrajectory* trj_2 = (KnuclTrajectory*)((*(evt->GetTrajectoryContainer()))[ii]);
         G4int tID_2 = trj_2->GetTrackID();
         std::vector<int>::iterator it_parent;
         for (it_parent = parent_list.begin(); it_parent != parent_list.end(); ++it_parent )
         {
            if ( tID_2 == *it_parent ) {
              for (int iii=0; iii<trj_2->GetPointEntries(); iii++) {
                double xx_2 = (trj_2->GetPoint(iii)->GetPosition()).x();
                double yy_2 = (trj_2->GetPoint(iii)->GetPosition()).y();
                double zz_2 = (trj_2->GetPoint(iii)->GetPosition()).z();
                if (iii==0) {
                  prev_pos.SetXYZ(xx_2,yy_2,zz_2);
                } else {
                  curr_pos.SetXYZ(xx_2,yy_2,zz_2);
                  TVector3 delta = curr_pos - prev_pos;
                  FlightLength += delta.Mag();
                  prev_pos.SetXYZ(xx_2,yy_2,zz_2);
                }
              }
            }
         }
       }
     }
     tmp.setFlightLength( FlightLength );
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

