#ifndef KNUCLTRIGGER_HH
#define KNUCLTRIGGER_HH 1

#include "KnuclRootData.h"
#include "KnuclCommon.h"
#include "globals.hh"


bool trigger(const std::vector<DetectorHit> &hits, G4String name);

#endif
