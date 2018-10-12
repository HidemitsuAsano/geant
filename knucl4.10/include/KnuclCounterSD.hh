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
// $Id: KnuclCounterSD.hh,v 1.1.1.1 2013/12/25 01:25:10 sakuma Exp $
// GEANT4 tag $Name:  $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef KnuclCounterSD_h
#define KnuclCounterSD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"
#include <map>

class KnuclDetectorConstruction;
class G4HCofThisEvent;
class G4Step;
#include "KnuclHit.hh"
#include "KnuclAnaManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class KnuclCounterSD : public G4VSensitiveDetector
{
  public:
  
  KnuclCounterSD(G4String, KnuclDetectorConstruction*, G4double);
     ~KnuclCounterSD();

      void Initialize(G4HCofThisEvent*);
      G4bool ProcessHits(G4Step*,G4TouchableHistory*);
      void EndOfEvent(G4HCofThisEvent*);
      void clear();
      void DrawAll();
      void PrintAll();

  private:

  G4double ADCthreshold;

  static const unsigned int SMASK    = 0x01FF;      /* S Mask 9 Bits (0-511) */
  static const unsigned int LMASK    = 0x001F;      /* L Mask 5 Bits (0-31) */
  static const unsigned int CMASK    = 0x007F;      /* C Mask 7 Bits (0-127) */
  static const G4int          SSHIFT   =  4;
  static const G4int          LSHIFT   = 14;
  static const G4int          CSHIFT   = 20;
  G4int KEY(G4int cid,G4int layer,G4int seg){
    return  ( (((cid)&CMASK )<<CSHIFT ) | (((seg)&SMASK)<<SSHIFT) | (((layer)&LMASK)<<LSHIFT) );
  }
  KnuclHitsCollection* CalCollection;      
  KnuclDetectorConstruction* Detector;
  std::map<G4int,G4int> HitIDContainer;
  std::map<G4int,G4bool> ADCthContainer;
  //  G4int HitNCID[7][16], HitCVCID[34], HitPCID[27], HitCDHID[36], HitIHID[24], HitKDVID, HitACID, HitT0ID[5], HitBVCID, HitLC1ID, HitLC2ID, HitBPDID[70], HitLSID;
  G4String detectorname;
};

#endif

