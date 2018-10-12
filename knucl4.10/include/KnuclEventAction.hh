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
// $Id: KnuclEventAction.hh,v 1.3 2016/11/09 08:41:03 inoue Exp $
// GEANT4 tag $Name:  $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
#ifndef KnuclEventAction_h
#define KnuclEventAction_h 1

#include "G4UserEventAction.hh"
#include "G4Trajectory.hh"
#include "G4PrimaryParticle.hh"

#include "KnuclTrajectory.hh"

#include "globals.hh"
#include "KnuclCommon.h"

#include "TVector3.h"

#include "KnuclTrigger.hh"

class G4Event;
class KnuclAnaManager;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


class KnuclEventAction : public G4UserEventAction
{
  public:
  KnuclEventAction(KnuclAnaManager* ana);
   ~KnuclEventAction();
  
  public:
    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);

  private:
    KnuclAnaManager* AnaManager;
    G4int counterCollID;
    G4int chamberCollID;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
