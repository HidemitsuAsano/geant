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
// $Id: KnuclPrimaryGeneratorAction.hh,v 1.2 2015/06/18 01:46:32 sakuma Exp $
// GEANT4 tag $Name:  $
//

#ifndef KnuclPrimaryGeneratorAction_h
#define KnuclPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "KnuclAnaManager.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

#include "TGenPhaseSpace.h"

#include "ComCrossSectionTable.hh"

#include <TGraph.h>
#include <TFile.h>
#include <TF1.h>


// ############################################################# //
// ### === 3He ===                                           ### //
// ### momentum distributions for two- and three-body decay  ### //
// ###   measured by quasifree (e,e'p)                       ### //
// ###   -- Ref:: Phys.Rev.Lett. 49 (1982) 974 --            ### //
// ### (4He:: Nucl.Phys. A355 (1981) 333 [not used in here]) ### //
// ### === d ===                                             ### //
// ### momentum distributions measured by quasifree (e,e'p)  ### //
// ###   -- Ref:: Nucl. Phys. A365 (1981) 349 --             ### //
// ### $$$ parameters $$$                                    ### //
// ###   triple-gaussian parameters are obtained by data-fit ### //
// ###   (ag:/home/sakuma/work/ana/calc/double_gauss.C)      ### //
// ###       unit in p^{-3}                                  ### //
// ############################################################# //
//old param until 20150617
//const G4double FermiMotion_3HeTwoBody[6]   = { 241.842, 83.9673, 426.396, 44.488, 0, 0 };
//const G4double FermiMotion_3HeThreeBody[6] = { 65.7717, 99.8646, 134.471, 58.0724, 0, 0 };
const G4double FermiMotion_3HeTwoBody[6]   = { 406.113, 44.0577, 244.472, 79.1843, 16.9276, 125.564 };
const G4double FermiMotion_3HeThreeBody[6] = { 96.6054, 42.8724, 117.396, 84.6603, 4.63981, 151.192 };
const G4double FermiMotion_deuteron[6]     = { 1162.54, 37.6533, 124.49, 79.0454, 4.2879, 156.929};

class G4ParticleGun;
class G4Event;

class KnuclPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    KnuclPrimaryGeneratorAction(KnuclAnaManager* ana);
    ~KnuclPrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event* anEvent);

// %%%%%%%%%%%%%%%%%%%% //
// % common functions % //
// %%%%%%%%%%%%%%%%%%%% //
    float Legendre(int n, float x);
    bool FermiMom_judge(G4int tgtID, G4int nbodyFlag, G4double mom); // nbodyFlag 2: two-body, 3: three-body in 3He
    double FermiMom_gen(G4int tgtID, G4int nbodyFlag); // nbodyFlag 2: two-body, 3: three-body in 3He
    int ReadFromFile(G4Event* anEvent);
    void WriteOutputOscarFile(G4Event* anEvent, int reacID=0);

// %%%%%%%%%%%%%%%%%%%% //
// % K-3He study      % //
// %%%%%%%%%%%%%%%%%%%% //
    int Kminus3He(G4Event* anEvent);
    int KminusReac(G4Event* anEvent, const CrossSection& cs);
    bool ManyBody(G4double CMmass, G4int nBody,
		  const G4double* mass, G4ThreeVector*vec,
		  G4int nFin,
		  G4int nparam, const G4double* param, G4double max);
    bool ManyBody(G4double CMmass, G4int nBody,
		  const G4double* mass, G4ThreeVector* vec);
    G4double GetMass(const G4ParticleDefinition& particle);
    G4ThreeVector RandUnitVec();
    const CrossSection& GenerateReac();
    void PrintCS(const CrossSection& cs);
    void PrintAllCS();


  private:
    KnuclAnaManager* anaManager;
    G4ParticleGun* particleGun;
    
    TGenPhaseSpace *gen;
    CrossSectionTable* csTable;

    ReactionData *reactionData;

    TF1* FermiMotion_3HeTwoBody_dist;
    TF1* FermiMotion_3HeThreeBody_dist;
    TF1* FermiMotion_deuteron_dist;
};

#endif


