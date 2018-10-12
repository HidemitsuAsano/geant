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
// $Id: KnuclHit.cc,v 1.1.1.1 2013/12/25 01:25:10 sakuma Exp $
// GEANT4 tag $Name:  $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "KnuclHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"

G4Allocator<KnuclHit> KnuclHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

KnuclHit::KnuclHit()
{
  detectorID = -1;
  layerID = -1;
  channelID = -1;
  EdepAbs = 0.;
  TrackLengthAbs = 0.;
  Time = 0.;
  Pos[XCOORD]=-9999.9;
  Pos[YCOORD]=-9999.9;
  Pos[ZCOORD]=-9999.9;
  LocalPos[XCOORD]=-9999.9;
  LocalPos[YCOORD]=-9999.9;
  LocalPos[ZCOORD]=-9999.9;
  Momentum[XCOORD]=-9999.9;
  Momentum[YCOORD]=-9999.9;
  Momentum[ZCOORD]=-9999.9;
  PDGCode  = -1;
  ParentID = -1;
  TrackID  = -1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

KnuclHit::~KnuclHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

KnuclHit::KnuclHit(const KnuclHit& right)
  : G4VHit()
{
  detectorID = right.detectorID;
  layerID = right.layerID;
  channelID = right.channelID;
  EdepAbs = right.EdepAbs;
  TrackLengthAbs = right.TrackLengthAbs;
  Time = right.Time;
  Pos[XCOORD]=right.Pos[XCOORD];
  Pos[YCOORD]=right.Pos[YCOORD];
  Pos[ZCOORD]=right.Pos[ZCOORD];
  LocalPos[XCOORD]=right.LocalPos[XCOORD];
  LocalPos[YCOORD]=right.LocalPos[YCOORD];
  LocalPos[ZCOORD]=right.LocalPos[ZCOORD];
  Momentum[XCOORD]=right.Momentum[XCOORD];
  Momentum[YCOORD]=right.Momentum[YCOORD];
  Momentum[ZCOORD]=right.Momentum[ZCOORD];
  PDGCode = right.PDGCode;
  ParentID= right.ParentID;
  TrackID = right.TrackID;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const KnuclHit& KnuclHit::operator=(const KnuclHit& right)
{
  detectorID = right.detectorID;
  layerID = right.layerID;
  channelID = right.channelID;
  EdepAbs = right.EdepAbs;
  TrackLengthAbs = right.TrackLengthAbs;
  Time = right.Time;
  Pos[XCOORD]=right.Pos[XCOORD];
  Pos[YCOORD]=right.Pos[YCOORD];
  Pos[ZCOORD]=right.Pos[ZCOORD];
  LocalPos[XCOORD]=right.LocalPos[XCOORD];
  LocalPos[YCOORD]=right.LocalPos[YCOORD];
  LocalPos[ZCOORD]=right.LocalPos[ZCOORD];
  Momentum[XCOORD]=right.Momentum[XCOORD];
  Momentum[YCOORD]=right.Momentum[YCOORD];
  Momentum[ZCOORD]=right.Momentum[ZCOORD];
  PDGCode = right.PDGCode;
  ParentID= right.ParentID;
  TrackID = right.TrackID;

  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int KnuclHit::operator==(const KnuclHit& right) const
{
  return (this==&right) ? 1 : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void KnuclHit::Draw()
{
//  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
//  if(pVVisManager&&detectorID==10&&(EdepAbs>2.0*MeV)){
//    G4Tubs SideTOF_tube("dummy", 40.0*cm, 41.8*cm,
//	                         49.0*cm, 0, twopi/24.0);
//    G4ThreeVector xyzChamber(0.0, 0.0, 0.0);
//    G4RotationMatrix* rotChamber = new G4RotationMatrix;
//    rotChamber->rotateZ(twopi/24.0*(channelID-1.5));
//    G4Transform3D posChamber(*rotChamber, xyzChamber);
//  
//    G4VisAttributes attribs;
//    G4Color color(1.0,1.0,0.0);
//    attribs.SetColor(color);
//    attribs.SetForceSolid(true);
//    pVVisManager->Draw(SideTOF_tube, attribs, posChamber);
//  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void KnuclHit::Print()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

