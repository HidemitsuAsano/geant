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
// $Id: KnuclCounterSD.cc,v 1.1.1.1 2013/12/25 01:25:10 sakuma Exp $
// GEANT4 tag $Name:  $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "KnuclCounterSD.hh"

#include "KnuclHit.hh"
#include "KnuclDetectorConstruction.hh"

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
#include "DetectorList.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

KnuclCounterSD::KnuclCounterSD(G4String name,
			       KnuclDetectorConstruction* det,
			       G4double th)
  : G4VSensitiveDetector(name), Detector(det)
{
  detectorname = name;
  collectionName.insert(name);
  ADCthreshold = th*MeV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

KnuclCounterSD::~KnuclCounterSD()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void KnuclCounterSD::Initialize(G4HCofThisEvent*)
{
  CalCollection = new KnuclHitsCollection
    (SensitiveDetectorName,collectionName[0]); 
  HitIDContainer.clear();
  ADCthContainer.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool KnuclCounterSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  const G4StepPoint* preStepPoint= aStep-> GetPreStepPoint();

  G4double edep = aStep->GetTotalEnergyDeposit();
  G4double time = aStep->GetTrack()->GetGlobalTime();
  //G4int    pdg  = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();
 
  G4double stepl = 0.;
  //if (aStep->GetTrack()->GetDefinition()->GetPDGCharge() != 0. || (pdg==2112) )
      stepl = aStep->GetStepLength();
      
  if ((edep==0.)&&(stepl==0.)) return false;      

  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
    
  G4VPhysicalVolume* physVol = theTouchable->GetVolume(); 
  G4ThreeVector position = physVol->GetObjectTranslation();

  G4int detectorID = -1;
  G4int layerID = -1;
  G4int channelID = -1;
  G4String name = physVol->GetName();

  DetectorList *dlist=DetectorList::GetInstance();

  G4String tmp=name;
  if(name.find_first_of('_')!=std::string::npos)
    tmp.remove(name.find_first_of('_') , name.size());
  if(tmp.find_first_of('s')!=std::string::npos)
    tmp.remove(tmp.find_first_of('s') , 1);
 
  detectorID=dlist->GetCID(tmp);
  if(detectorID==-1) return false;
  tmp=name;

  layerID=0;
  channelID=physVol->GetCopyNo();
  //  G4cout<<detectorID<<" "<<layerID<<" "<<channelID<<G4endl;
  
  G4ThreeVector trackpos = preStepPoint->GetPosition(); //--> wrong NC ???120607
  time =                   preStepPoint->GetGlobalTime();
  
  G4int tmpkey=KEY(detectorID,layerID,channelID);

  std::map<G4int,G4int>::const_iterator ihitID = HitIDContainer.find(tmpkey);
  std::map<G4int,G4bool>::const_iterator iadcID = ADCthContainer.find(tmpkey);

  G4int hitID=-1;
  if( ihitID==HitIDContainer.end() )
    { 
      KnuclHit* calHit = new KnuclHit();
      calHit->SetDetectorID(detectorID);
      calHit->SetLayerID(layerID);
      calHit->SetChannelID(channelID);
      calHit->AddAbs(edep, stepl);
      calHit->SetParentID(aStep->GetTrack()->GetParentID());
      calHit->SetTrackID(aStep->GetTrack()->GetTrackID());
      calHit->SetTime(time);
      calHit->SetPos(trackpos.x(), trackpos.y(), trackpos.z());
      calHit->SetLocalPos(calHit->CalcLocalPos(aStep));
      //      calHit->SetMomentum(aStep->GetTrack()->GetMomentum());
      calHit->SetMomentum(preStepPoint->GetMomentum());
      calHit->SetPDGCode(aStep->GetTrack()->GetDefinition()->GetPDGEncoding());
      hitID= CalCollection->insert(calHit) - 1;
      HitIDContainer[tmpkey] =hitID;
    }
  else
    { 
      hitID= ihitID->second;
      (*CalCollection)[hitID]->AddAbs(edep,stepl);
    }

  if( iadcID==ADCthContainer.end() && ADCthreshold>=0 &&
      (*CalCollection)[hitID]->GetEdepAbs() > ADCthreshold )
    {
      ADCthContainer[tmpkey]=true;
      (*CalCollection)[hitID]->SetTime(time);
      (*CalCollection)[hitID]->SetPos(trackpos.x(), trackpos.y(), trackpos.z());
      (*CalCollection)[hitID]->SetLocalPos((*CalCollection)[hitID]->CalcLocalPos(aStep));
      //      (*CalCollection)[hitID]->SetMomentum(aStep->GetTrack()->GetMomentum());
      (*CalCollection)[hitID]->SetMomentum(preStepPoint->GetMomentum());
      (*CalCollection)[hitID]->SetPDGCode(aStep->GetTrack()->GetDefinition()->GetPDGEncoding());
      (*CalCollection)[hitID]->SetParentID(aStep->GetTrack()->GetParentID());
      (*CalCollection)[hitID]->SetTrackID(aStep->GetTrack()->GetTrackID());
    }
  
  //  G4cerr<<hitID<<" "<<detectorID<<" "<<layerID<<" "<<channelID<<" "<<aStep->GetTrack()->GetPosition()<<G4endl;

#if 0
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
    {
      G4Circle circle(trackpos);
      circle.SetScreenSize(10.);
      circle.SetFillStyle(G4Circle::filled);
      G4VisAttributes attribs(Yellow);
      if (detectorID==CID_CDH){
	attribs.SetColor(White);
      }
      if (detectorID==CID_IH){
	attribs.SetColor(Blue);
      }
      circle.SetVisAttributes(attribs);
      pVVisManager->Draw(circle);
    }
#endif

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void KnuclCounterSD::EndOfEvent(G4HCofThisEvent* HCE)
{
  static G4int HCID = -1;
  if(HCID<0)
  { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  HCE->AddHitsCollection(HCID,CalCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void KnuclCounterSD::clear()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void KnuclCounterSD::DrawAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void KnuclCounterSD::PrintAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

