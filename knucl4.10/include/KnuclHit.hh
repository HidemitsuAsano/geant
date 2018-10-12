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
// $Id: KnuclHit.hh,v 1.1.1.1 2013/12/25 01:25:10 sakuma Exp $
// GEANT4 tag $Name:  $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef KnuclHit_h
#define KnuclHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"

#include "KnuclCommon.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class KnuclHit : public G4VHit
{
private:
  G4int detectorID; // 0-, see KnuclCommon.h
  G4int layerID;    // 0-
  G4int channelID;  // 0-
  G4double EdepAbs, TrackLengthAbs;
  G4double Time;
  G4double Pos[XYZ];
  G4double LocalPos[XYZ];
  G4double Momentum[XYZ];
  G4double Dx;
  G4int PDGCode;
  G4int ParentID;
  G4int TrackID;
  
public:
  
  KnuclHit();
  ~KnuclHit();
  KnuclHit(const KnuclHit&);
  const KnuclHit& operator=(const KnuclHit&);
  G4int operator==(const KnuclHit&) const;
  
  inline void* operator new(size_t);
  inline void  operator delete(void*);
  
  void Draw();
  void Print();
  
public:
  
  void SetDetectorID(G4int i) {detectorID = i;};
  void SetLayerID(G4int i) {layerID = i;};
  void SetChannelID(G4int i) {channelID = i;};
  void AddAbs(G4double de, G4double dl) {EdepAbs += de; TrackLengthAbs += dl;};
  void SetTime(G4double time) {Time =time;};
  void SetPos(G4double x, G4double y, G4double z)
  {Pos[XCOORD]=x;Pos[YCOORD]=y;Pos[ZCOORD]=z;};
  void SetPos(const G4ThreeVector& pos)
  {Pos[XCOORD]=pos.x();Pos[YCOORD]=pos.y();Pos[ZCOORD]=pos.z();};
  void SetLocalPos(G4double x, G4double y, G4double z)
  {LocalPos[XCOORD]=x;LocalPos[YCOORD]=y;LocalPos[ZCOORD]=z;};
  void SetLocalPos(const G4ThreeVector& pos)
  {LocalPos[XCOORD]=pos.x();LocalPos[YCOORD]=pos.y();LocalPos[ZCOORD]=pos.z();};
  void SetMomentum(G4double x, G4double y, G4double z)
  {Momentum[XCOORD]=x;Momentum[YCOORD]=y;Momentum[ZCOORD]=z;};
  void SetMomentum(const G4ThreeVector& momentum)
  {Momentum[XCOORD]=momentum.x();Momentum[YCOORD]=momentum.y();Momentum[ZCOORD]=momentum.z();};
  void SetDx(G4double v){Dx=v;}
  void SetPDGCode(G4int i){PDGCode=i;}
  void SetParentID(G4int i){ParentID=i;}
  void SetTrackID(G4int i){TrackID=i;}

  G4int GetDetectorID() { return detectorID;};
  G4int GetLayerID() { return layerID;};
  G4int GetChannelID() { return channelID;};
  G4double GetEdepAbs()     { return EdepAbs; };
  G4double GetTrakAbs()     { return TrackLengthAbs; };
  G4double GetTime()     { return Time; };
  void GetPos(G4ThreeVector &pos) 
  { pos.set(Pos[XCOORD],Pos[YCOORD],Pos[ZCOORD]);};
  G4ThreeVector GetPos()
  {return G4ThreeVector(Pos[XCOORD],Pos[YCOORD],Pos[ZCOORD]);}
  void GetLocalPos(G4ThreeVector &pos) 
  { pos.set(LocalPos[XCOORD],LocalPos[YCOORD],LocalPos[ZCOORD]);};
  G4ThreeVector GetLocalPos()
  {return G4ThreeVector(LocalPos[XCOORD],LocalPos[YCOORD],LocalPos[ZCOORD]);}
  void GetMomentum(G4ThreeVector &momentum) 
  { momentum.set(Momentum[XCOORD],Momentum[YCOORD],Momentum[ZCOORD]);};
  G4ThreeVector GetMomentum()
  {return G4ThreeVector(Momentum[XCOORD],Momentum[YCOORD],Momentum[ZCOORD]);}
  G4double GetDx(){return Dx;}
  G4int GetPDGCode(){return PDGCode;}
  G4int GetParentID(){return ParentID;}
  G4int GetTrackID(){return TrackID;}

  G4ThreeVector CalcLocalPos(G4Step* aStep);
  G4ThreeVector CalcPos(G4Step* aStep, G4ThreeVector localpos);
  G4ThreeVector CalcLocalVec(G4Step* aStep, G4ThreeVector vec);
  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<KnuclHit> KnuclHitsCollection;

extern G4Allocator<KnuclHit> KnuclHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* KnuclHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) KnuclHitAllocator.MallocSingle();
  return aHit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void KnuclHit::operator delete(void *aHit)
{
  KnuclHitAllocator.FreeSingle((KnuclHit*) aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4ThreeVector KnuclHit::CalcLocalPos(G4Step* aStep)
{
  //G4StepPoint* p1 = aStep->GetPostStepPoint();
  G4StepPoint* p1 = aStep->GetPreStepPoint();
  G4ThreeVector coord1 = p1->GetPosition();
    
  const G4AffineTransform transformation =
    p1->GetTouchable()->
    GetHistory()->GetTopTransform();
  G4ThreeVector local = transformation.TransformPoint(coord1);

  return local;
}

inline G4ThreeVector KnuclHit::CalcPos(G4Step* aStep, G4ThreeVector localpos)
{
  //G4StepPoint* p1 = aStep->GetPostStepPoint();
  G4StepPoint* p1 = aStep->GetPreStepPoint();
    
  const G4AffineTransform transformation =
    p1->GetTouchable()->
    GetHistory()->GetTopTransform();
  G4ThreeVector pos = transformation.Inverse().TransformPoint(localpos);

  return pos;
}

inline G4ThreeVector KnuclHit::CalcLocalVec(G4Step* aStep, G4ThreeVector vec)
{
  //G4StepPoint* p1 = aStep->GetPostStepPoint();
  G4StepPoint* p1 = aStep->GetPreStepPoint();
    
  const G4AffineTransform transformation =
    p1->GetTouchable()->
    GetHistory()->GetTopTransform();
  G4ThreeVector local = transformation.TransformAxis(vec);

  return local;
}

#endif
