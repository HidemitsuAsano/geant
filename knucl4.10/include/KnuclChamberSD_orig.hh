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
// $Id: KnuclChamberSD.hh,v 1.1.1.1 2013/12/25 01:25:10 sakuma Exp $
// GEANT4 tag $Name:  $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef KnuclChamberSD_h
#define KnuclChamberSD_h 1

#include <vector>

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

#include "G4SystemOfUnits.hh"

class KnuclDetectorConstruction;
class G4HCofThisEvent;
class G4Step;
#include "KnuclHit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// from K18ana, class CDSTrack
void LineToLine( const G4ThreeVector &x1, const G4ThreeVector &a1,
		 const G4ThreeVector &x2, const G4ThreeVector &a2,
		 G4double &dist,
		 G4ThreeVector &n1, G4ThreeVector &n2 );

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class KnuclChamberSD : public G4VSensitiveDetector
{
public:
  
  KnuclChamberSD(G4String, KnuclDetectorConstruction*, double param[200][20]);
  ~KnuclChamberSD();
  
  void Initialize(G4HCofThisEvent*);
  G4bool ProcessHits(G4Step*,G4TouchableHistory*);
  void EndOfEvent(G4HCofThisEvent*);
  void clear();
  void DrawAll();
  void PrintAll();
  G4ThreeVector CDCLocalWirePos(G4int layer, G4double z);
  
private:
  
  KnuclHitsCollection*  CalCollection;      
  KnuclDetectorConstruction*   Detector;
  G4int                   CDC_HitID[15][200];
  G4int                   BLC1_HitID[8][30];
  G4int                   BLC2_HitID[8][30];
  G4int                   FDC1_HitID[6][64];
  G4int                   FDC2_HitID[6][128];
  G4int                   BPC_HitID[8][16];
  G4String                detectorname;

  G4double CDC_CELL_OFFSET[15];
  G4double CDC_RADIUS[15];
  G4double CDCLength;
};

inline G4ThreeVector KnuclChamberSD::CDCLocalWirePos(G4int layer, G4double z)
{
  G4double offset = CDC_CELL_OFFSET[layer]; //degree
  G4double x = CDC_RADIUS[layer]*cos(offset/2.0); //mm
  G4double y = x*tan(offset/2.0)*2.0*z*mm/(CDCLength); //mm
  return G4ThreeVector(x,y,z);
}

#endif

