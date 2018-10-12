#ifndef KnuclTrajectoryPoint_h
#define KnuclTrajectoryPoint_h 1

#include "G4VTrajectory.hh"
#include "G4Allocator.hh"
#include <stdlib.h>
#include "G4ThreeVector.hh"
#include "G4ios.hh"
#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4TrajectoryPoint.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include <vector>

class KnuclTrajectoryPoint : public G4TrajectoryPoint 
{
public:
  KnuclTrajectoryPoint(){fPosition = G4ThreeVector(0.,0.,0.); DetectorID=-1;}
  KnuclTrajectoryPoint(G4ThreeVector pos,G4int id) {fPosition = pos; DetectorID=id;}
  ~KnuclTrajectoryPoint(){}
  void   SetDetectorID(G4int id){DetectorID=id;}
  G4int  GetDetectorID(){return DetectorID;}
private:
  G4int DetectorID;
}

#endif
