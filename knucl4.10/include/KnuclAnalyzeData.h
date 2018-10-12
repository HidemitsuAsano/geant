#ifndef KnuclAnalyzeData_h
#define KnuclAnalyzeData_h 1

#include "G4ThreeVector.hh"
#include "KnuclCommon.h"

const G4int CounterMultiMax0=20;

struct Counter1ev0 {
  G4int hitLay;
  G4int hitSeg;
  G4double time;
  G4ThreeVector hitPos;
};

class KnuclAnalyzeData{
  public:
    Counter1ev0 Tof1ev[CounterMultiMax0];
    G4double mom_in;
    G4double mom_out;
    G4double theta_out;
    G4double missmass;
};

#endif
