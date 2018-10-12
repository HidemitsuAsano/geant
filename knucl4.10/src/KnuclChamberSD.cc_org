//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: KnuclChamberSD.cc,v 1.1.1.1 2013/12/25 01:25:10 sakuma Exp $
// GEANT4 tag $Name:  $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "KnuclChamberSD.hh"
#include "KnuclHit.hh"
#include "KnuclDetectorConstruction.hh"
#include "KnuclCommon.h"

#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"

#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include "G4ios.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// from K18ana, class CDSTrack
void LineToLine( const G4ThreeVector &x1, const G4ThreeVector &a1,
		 const G4ThreeVector &x2, const G4ThreeVector &a2,
		 G4double &dist,
		 G4ThreeVector &n1, G4ThreeVector &n2 )
{
  G4ThreeVector x = x2-x1;   // x = x1 + t*a1
  G4double a =  a1.dot(a1);  // x = x2 + s*a2
  G4double b = -a1.dot(a2);  //    ||
  G4double c =  a2.dot(a1);  //    \/
  G4double d = -a2.dot(a2);  // a*t + b*s = A1
  G4double A1 = a1.dot(x);   // c*t + d*s = A2
  G4double A2 = a2.dot(x);

  G4double D = a*d-b*c;

  G4ThreeVector x2p;
  if( fabs(D)<0.00000000000001 ){
    dist = sqrt(x.mag2()-A1*A1);
    n1 = G4ThreeVector(0,0,0);
    n2 = G4ThreeVector(0,0,0);
  }
  else{
    G4double ss = (a*A2-c*A1)/D;
    G4double tt = (d*A1-b*A2)/D;
    n1 = x1 + tt*a1;
    n2 = x2 + ss*a2;
    dist = (n1-n2).mag();
  }

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

KnuclChamberSD::KnuclChamberSD(G4String name,
			       KnuclDetectorConstruction* det, double param[200][20])
  :G4VSensitiveDetector(name),Detector(det)
{
  detectorname = name;
  collectionName.insert(name);
  CDCLength=param[0][9]*cm;
#if 0
  std::cout<<"CDCLenght: "<<CDCLength<<std::endl;
#endif
  for(int i=0;i<15;i++){
    CDC_CELL_OFFSET[i]=-param[i+1][4]*degree;
    CDC_RADIUS[i]=param[i+1][1]*cm;
    G4ThreeVector posp=CDCLocalWirePos(i,-CDCLength/2.);
    G4ThreeVector pos=CDCLocalWirePos(i,CDCLength/2.);
#if 0
    std::cout<<"Layer "<<i+1<<"  "<<CDC_RADIUS[i]<<"  "<<CDC_CELL_OFFSET[i]/degree
	     <<"  "<<pos <<"  "<<pos.phi()/degree
	     <<"  "<<posp<<"  "<<posp.phi()/degree<<std::endl;
#endif
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

KnuclChamberSD::~KnuclChamberSD()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void KnuclChamberSD::Initialize(G4HCofThisEvent*)
{
  CalCollection = new KnuclHitsCollection
                      (SensitiveDetectorName,collectionName[0]); 
  for (G4int i=0; i<15; i++){
    for (G4int j=0; j<200; j++){
      CDC_HitID[i][j] = -1;
    }
  }
  for (G4int i=0; i<8; i++){
    for (G4int j=0; j<30; j++){
      BLC1_HitID[i][j] = -1;
      BLC2_HitID[i][j] = -1;
    }
  }
  for (G4int i=0; i<6; i++){
    for (G4int j=0; j<64; j++){
      FDC1_HitID[i][j] = -1;
    }
    for (G4int j=0; j<128; j++){
      FDC2_HitID[i][j] = -1;
    }
  }
  for (G4int i=0; i<8; i++){
    for (G4int j=0; j<15; j++){
      BPC_HitID[i][j] = -1;
    }
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool KnuclChamberSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  const G4StepPoint* preStepPoint= aStep-> GetPreStepPoint();

  G4double edep = aStep->GetTotalEnergyDeposit();
  G4double time = aStep->GetTrack()->GetGlobalTime();
  G4int    pdg  = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();  
  
  G4double stepl = 0.;
  if ((aStep->GetTrack()->GetDefinition()->GetPDGCharge()  != 0.) || 
      (pdg==22 || pdg==2112)   )
      stepl = aStep->GetStepLength();

  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
    
  G4VPhysicalVolume* physVol = theTouchable->GetVolume(); 
  G4ThreeVector position = physVol->GetObjectTranslation();

  G4int trackID = aStep->GetTrack()->GetTrackID();
  G4int detectorID = -1;
  G4int layerID = -1;
  G4int channelID = -1;
  G4String name = physVol->GetName();

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
  // chamber-hit reduction is very simple, so it must be imoproved.
  //    TODO :: store all hits => sort hits by layer & channel 
  //               => timing selection => final hit
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
  if (name.find(G4String("BLC2a"))!=G4String::npos){
    //--- only charged particles ---//
    if (aStep->GetTrack()->GetDefinition()->GetPDGCharge()==0) return false;
    detectorID=CID_BLC2a;
    G4String tmp=name;
    layerID=atoi(tmp.erase(0,11).erase(1).c_str())-1;
    channelID=physVol->GetCopyNo();
    if(BLC1_HitID[layerID][channelID]==-1) BLC1_HitID[layerID][channelID]=1;
    else return false;
  } 
  else if (name.find(G4String("BLC2b"))!=G4String::npos){
    //--- only charged particles ---//
    if (aStep->GetTrack()->GetDefinition()->GetPDGCharge()==0) return false;
    detectorID=CID_BLC2b;
    G4String tmp=name;
    layerID=atoi(tmp.erase(0,11).erase(1).c_str())-1;
    channelID=physVol->GetCopyNo();
    if(BLC2_HitID[layerID][channelID]==-1) BLC2_HitID[layerID][channelID]=1;
    else return false;
  }
  else if (name.find(G4String("FDC1"))!=G4String::npos){
    //--- only charged particles ---//
    if (aStep->GetTrack()->GetDefinition()->GetPDGCharge()==0) return false;
    detectorID=CID_FDC1;
    G4String tmp=name;
    layerID   = atoi(tmp.remove(0, 10).remove(1, 2))-1;
    channelID = physVol->GetCopyNo();
    if(FDC1_HitID[layerID][channelID]==-1) FDC1_HitID[layerID][channelID]=1;
    else return false;
  }
  else if (name.find(G4String("FDC2"))!=G4String::npos){
    //--- only charged particles ---//
    if (aStep->GetTrack()->GetDefinition()->GetPDGCharge()==0) return false;
    detectorID=CID_FDC2;
    G4String tmp=name;
    layerID   = atoi(tmp.remove(0, 10).remove(1, 3))-1;
    channelID = physVol->GetCopyNo();
    if(FDC2_HitID[layerID][channelID]==-1) FDC2_HitID[layerID][channelID]=1;
    else return false;
  }
  else if (name.find(G4String("CDC"))!=G4String::npos){
    //--- only charged particles ---//
    if (aStep->GetTrack()->GetDefinition()->GetPDGCharge()==0) return false;
    detectorID=CID_CDC;
    G4String tmp=name;
    layerID=atoi(tmp.erase(0,8).erase(2,4).c_str())-1;
    channelID=physVol->GetCopyNo();
    if(CDC_HitID[layerID][channelID]==-1) CDC_HitID[layerID][channelID]=1;
    else return false;
  }
  else if (name.find(G4String("BPC"))!=G4String::npos){
    //--- only charged particles ---//
    if (aStep->GetTrack()->GetDefinition()->GetPDGCharge()==0) return false;
    detectorID=CID_BPC;
    G4String tmp=name;
    layerID=atoi(tmp.erase(0,9).erase(1,3).c_str())-1;
    channelID=physVol->GetCopyNo();
    //if(   trackID ==1 ) return false; // comment out, 20120709 sakuma
    //std::cout<<"BPC "<<layerID<<" "<<channelID<<" "<<BPC_HitID[layerID][channelID]<<std::endl;
    if(BPC_HitID[layerID][channelID]==-1) BPC_HitID[layerID][channelID]=-1;
    else return false;
  }
  //G4cerr<<detectorID<<" "<<layerID<<" "<<channelID<<G4endl;

  //G4ThreeVector trackpos = aStep->GetTrack()->GetPosition();
  G4ThreeVector trackpos = preStepPoint->GetPosition(); // --> need for CalcLocalVec()
  G4double x = trackpos.x();
  G4double y = trackpos.y();
  G4double z = trackpos.z();

  KnuclHit* calHit = new KnuclHit();
  calHit->SetDetectorID(detectorID);
  calHit->SetLayerID(layerID);
  calHit->SetChannelID(channelID);
  calHit->AddAbs(edep, stepl);
  calHit->SetTime(time);
  calHit->SetPos(x, y, z);
  calHit->SetLocalPos(calHit->CalcLocalPos(aStep));
  calHit->SetMomentum(aStep->GetTrack()->GetMomentum());
  calHit->SetDx(0);
  calHit->SetPDGCode(aStep->GetTrack()->GetDefinition()->GetPDGEncoding());
  calHit->SetParentID(aStep->GetTrack()->GetParentID());
  calHit->SetTrackID(trackID);
  CalCollection->insert(calHit);

  //----- CDC -----//
  if ( detectorID==CID_CDC ){
    //------------------------------------------------------------//
    // find minimum-distance positon from anode-wire
    //------------------------------------------------------------//
    // calc new_lpos
    G4ThreeVector lpos = calHit->GetLocalPos();
    G4ThreeVector mom = calHit->CalcLocalVec(aStep, preStepPoint->GetMomentumDirection());
    G4ThreeVector w  = CDCLocalWirePos(layerID, lpos.z());
    G4ThreeVector w0 = CDCLocalWirePos(layerID, 0);
    G4ThreeVector wdir = w-w0;
    G4double dx;
    G4ThreeVector w1, new_lpos;
    LineToLine(lpos, mom, w0, wdir, dx, new_lpos, w1);
    G4ThreeVector new_pos = calHit->CalcPos(aStep, new_lpos);
    //G4cerr<<layerID<<" "<<channelID<<" "<<lpos<<" "<<new_lpos<<" "<<dx<<G4endl;
    //G4cerr<<layerID<<" "<<channelID<<" "<<calHit->GetPos()<<" "<<new_pos<<" "<<dx<<G4endl;
    calHit->SetDx(dx);
    calHit->SetPos(new_pos);
    calHit->SetLocalPos(new_lpos);
    //G4cerr<<layerID<<" "<<channelID<<" "<<calHit->GetLocalPos()<<" "<<calHit->GetDx()<<G4endl;
  }else{
    G4ThreeVector lpos = calHit->GetLocalPos();
    G4ThreeVector pos  = calHit->CalcPos(aStep, lpos);
    calHit->SetPos(pos);
    calHit->SetLocalPos(lpos);
    calHit->SetDx(lpos.x());
    //G4cerr<<layerID<<" "<<channelID<<" "<<calHit->GetLocalPos()<<" "<<calHit->GetDx()<<G4endl;
  }

#if 1
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
    {
      G4Circle circle(trackpos);
      circle.SetScreenSize(10.);
      //circle.SetFillStyle(G4Circle::filled);
      G4VisAttributes attribs(Yellow);
      circle.SetVisAttributes(attribs);
      pVVisManager->Draw(circle);
    }
#endif
  
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void KnuclChamberSD::EndOfEvent(G4HCofThisEvent* HCE)
{
  static G4int HCID = -1;
  if(HCID<0)
  { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  HCE->AddHitsCollection(HCID,CalCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void KnuclChamberSD::clear()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void KnuclChamberSD::DrawAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void KnuclChamberSD::PrintAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

