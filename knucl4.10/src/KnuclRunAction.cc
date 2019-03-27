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
// $Id: KnuclRunAction.cc,v 1.3 2019/03/20 21:45:40 asano Exp $
// GEANT4 tag $Name:  $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "KnuclRunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4ios.hh"

#include "TRandom3.h"
#include "Randomize.hh"
#include <time.h>
#include <sys/time.h>

#include "CLHEP/Random/MTwistEngine.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

KnuclRunAction::KnuclRunAction(KnuclAnaManager* ana) : anaManager(ana)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

KnuclRunAction::~KnuclRunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void KnuclRunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

#if 1
  struct timeval tv;
  gettimeofday(&tv, NULL);
  srand(tv.tv_usec);
  //std::cerr<<tv.tv_usec<<" "<<time(0)<<std::endl;
#else
  srand(time(0));
#endif

  long seed;
  if(anaManager->GetSeed()=="random"){
    seed=long((double)rand()/((double)RAND_MAX+1)*10000);
    anaManager->SetSeedNum(seed);
  }else if(anaManager->GetSeed()=="external"){
    std::cout << "seed is given from the argument of knucl" << std::endl;
    seed =  anaManager->GetSeedNum();
    if(seed<0) std::cerr <<  "seed must be >= 0 !!!" << std::endl;
    std::cout << seed << std::endl;
  }else{
    seed=atoi(anaManager->GetSeed().Data());
    anaManager->SetSeedNum(seed);
  }
  CLHEP::HepRandom::setTheEngine(new CLHEP::MTwistEngine(seed));
  CLHEP::HepRandom::showEngineStatus();

  gRandom=new TRandom3(seed);

  anaManager->BeginOfRunAction();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void KnuclRunAction::EndOfRunAction(const G4Run*)
{
  anaManager->EndOfRunAction();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



