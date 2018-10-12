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
// $Id: KnuclSteppingAction.cc,v 1.2 2013/12/25 03:40:48 sakuma Exp $
// GEANT4 tag $Name:  $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "KnuclSteppingAction.hh"

#include "KnuclDetectorConstruction.hh"
#include "KnuclEventAction.hh"
#include "KnuclAnaManager.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4EventManager.hh"

#include <string>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

KnuclSteppingAction::KnuclSteppingAction(KnuclAnaManager* ana)
  : anaManager(ana)
{
  Verbosity = 0;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

KnuclSteppingAction::~KnuclSteppingAction() 
{ 
  Verbosity = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void KnuclSteppingAction::UserSteppingAction(const G4Step* aStep)
{
  G4Track* track = aStep->GetTrack();
  G4VPhysicalVolume* volume = track->GetVolume();
  G4String volname = volume->GetName();
  G4String parname = track->GetDefinition()->GetParticleName();
  G4double length = track->GetTrackLength();

#if 0  
  if(track->GetTrackID()==1){
    std::cout<<volname<<"   "<<track->GetMaterial()->GetName()<<"  "<<length<<"   "<<track->GetMomentum().mag()<<std::endl;
  }
#endif

#if 1
  G4double time   = track->GetLocalTime();
 
  if (time>1000000.0) {
    if(Verbosity){ 
      G4cout << " time out of range "<< G4endl;  
      G4cout << " # track maybe trapped in Magnet?? (time) # " << volname << " " <<  parname << " " << track->GetMomentum() << G4endl;  
    } 
    track->SetTrackStatus(fStopAndKill);
 } 
 
 if (length>100000.0) { // 100m, G4double expHall_z =25.0*m;
   if(Verbosity){
     G4cout << " length out of range "<< G4endl;  
     G4cout << " # track maybe trapped in Magnet?? (length) # " << volname << " " <<  parname << " " << track->GetMomentum() << G4endl;  
   }
   track->SetTrackStatus(fStopAndKill);
 } 
#endif


#if 0
 //----------------------------------------------------//
 // cut of low-momentum gamma, ekectron, and neucleon
 // in magnet yoke, shiled, beam-dump
 //----------------------------------------------------//
 //std::cerr<<parname<<" "<<volname<<" "<<track->GetMomentum().mag()<<std::endl;
 static const double threshold = 100*MeV;
 std::string vol = std::string(volume->GetName());
 if( vol.find("sDoraemon", 0)!=std::string::npos || // yoke or end-caps of Solenoid magnet
     vol.find("sUSWK",     0)!=std::string::npos || // yokes of Sweeping magnet
     vol.find("BeamDump",  0)!=std::string::npos ||
     vol.find("Shield",    0)!=std::string::npos ||
     vol.find("Con",       0)!=std::string::npos ||
     vol.find("Floor",     0)!=std::string::npos ){
   double mom = track->GetMomentum().mag();
   int pid = track->GetDefinition()->GetPDGEncoding();
   if( mom<threshold && (pid==22 || abs(pid)==11 || pid==2112 || pid==2212) ){ // gamma,electron,neucleon
     track->SetTrackStatus(fStopAndKill);
     //std::cerr<<"cut | "<<vol<<" "<<pid<<" "<<mom<<std::endl;
   }
 }
#endif

#if 0
 //----------------------------------------------------//
 // dump all steps
 //----------------------------------------------------//
 static int currentEventID = -1;
 int eventID = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();
 if( currentEventID!=eventID ){
   std::cerr<<"===== eventID:"<<eventID<<" ====="<<std::endl;
   currentEventID = eventID;
 }
 const G4VProcess* processDefinedTheStep
   = aStep->GetPostStepPoint()->GetProcessDefinedStep();
 std::cerr<<track->GetTrackID()<<","<<track->GetParentID()<<" : "
	  <<parname<<" ["<<track->GetMomentum().mag()<<" MeV/c, "<<length<<"mm] "
	  <<processDefinedTheStep->GetProcessName()<<" ("<<volname<<") "
	  <<std::endl;
#endif

#if 0
 //----------------------------------------------------//
 // dump reacltions
 //----------------------------------------------------//
 const G4VProcess* processDefinedTheStep
   = aStep->GetPostStepPoint()->GetProcessDefinedStep();
 if( (parname=="sigma+" || parname=="sigma0" || parname=="sigma-" || parname=="lambda") &&
     (processDefinedTheStep->GetProcessName()=="hBertiniCaptureAtRest" ||
      processDefinedTheStep->GetProcessName()=="CHIPSNuclearCaptureAtRest") )
   std::cerr<<track->GetTrackID()<<","<<track->GetParentID()<<" : "
	    <<parname<<" ["<<track->GetMomentum().mag()<<" MeV/c, "<<length<<"mm] "
	    <<processDefinedTheStep->GetProcessName()<<" ("<<volname<<") "
	    <<std::endl;
#endif

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



