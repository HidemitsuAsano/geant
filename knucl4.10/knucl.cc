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
// $Id: knucl.cc,v 1.6 2016/12/21 05:30:56 sakuma Exp $
// GEANT4 tag $Name:  $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - exampleN01 
//
// --------------------------------------------------------------
// Comments
//
// 
// --------------------------------------------------------------

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "G4PhysListFactory.hh"

#include "KnuclDetectorConstruction.hh"
#include "KnuclPrimaryGeneratorAction.hh"
#include "KnuclRunAction.hh"
#include "KnuclEventAction.hh"
#include "KnuclTrackingAction.hh"
#include "KnuclSteppingAction.hh"
#include "KnuclAnaManager.hh"
#include "KnuclVisManager.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "G4DecayTable.hh"
#include "G4PhaseSpaceDecayChannel.hh"
#include "G4VDecayChannel.hh"

//#include "G4SystemOfUnits.hh"

#include "stdlib.h"

void ModDecayBranch();
void AddKpp(KnuclAnaManager* anaManager);
void AddKpn(KnuclAnaManager* anaManager);

class G4Kpp;
class G4Kpn;

#define VISUALIZATION 0

int main(int argc,char** argv)
{
  std::string macro_name="";
  std::string knucl_card="KnuclSetting.card";
  std::string cs_list="";
  std::string shape_file="";
  std::string out_file="";
  for( int i=1; i<argc; i++ ){
    if( strstr(argv[i], ".mac") ) macro_name=argv[i];
    else if( strstr(argv[i], ".card") ) knucl_card=argv[i];
    else if( strstr(argv[i], ".list") ) cs_list=argv[i];
    else if( strstr(argv[i], "shape") && strstr(argv[i], ".root") ) shape_file=argv[i];
    else if( strstr(argv[i], ".root") ) out_file=argv[i];
    else{
      std::cout<<" !!! Invailed input !!! "<<argv[i]<<std::endl;
    }
  }
  std::cout<<"KnuclSetting.card="<<knucl_card<<std::endl;
  if( !macro_name.empty() ) std::cout<<"Macro="<<macro_name<<std::endl;
  if( !cs_list.empty() ) std::cout<<"CrossSectionTable="<<cs_list<<std::endl;
  if( !shape_file.empty() ) std::cout<<"CrossSectionTable="<<cs_list<<std::endl;
  if( !out_file.empty() ) std::cout<<"OutFile="<<out_file<<std::endl;

  // Construct the default run manager
  G4RunManager* runManager = new G4RunManager;

  //### for singleton ###//
  KnuclAnaManager* anaManager= KnuclAnaManager::GetPointer(knucl_card);
  //### for singleton ###//
  if( !out_file.empty() ){
    std::cout<<"OutputROOTFile : "<<out_file<<std::endl;
    anaManager->SetOutputRootFile(out_file);
  }
  if( !shape_file.empty() ){
    std::cout<<"KppShapeFile : "<<shape_file<<std::endl;
    anaManager->SetKppShapeFile(shape_file);
  }
  if( !cs_list.empty() ){
    std::cout<<"CSListFile : "<<cs_list<<std::endl;
    anaManager->SetCSFile(cs_list);
  }

  // set mandatory initialization classes
  runManager->SetUserInitialization(new KnuclDetectorConstruction(anaManager));

  G4String phys_name = anaManager->GetPhysicsList();
  G4int verbose = 1;
  G4PhysListFactory factory;
  G4VModularPhysicsList* physlist = factory.GetReferencePhysList(phys_name);
  factory.SetVerbose(verbose);
  runManager->SetUserInitialization(physlist);
  // === print avalable physics-list ===//
  //std::vector <G4String> list = factory.AvailablePhysLists();
  //for (std::vector<G4String>::iterator it = list.begin(); it != list.end(); ++it) std::cout<<*it<<std::endl;

#if 1
  ModDecayBranch();
  AddKpp(anaManager);
  AddKpn(anaManager);
#endif

  // set mandatory user action class
  runManager->SetUserAction(new KnuclPrimaryGeneratorAction(anaManager));
  runManager->SetUserAction(new KnuclRunAction(anaManager));
  runManager->SetUserAction(new KnuclEventAction(anaManager));
  runManager->SetUserAction(new KnuclTrackingAction);
  runManager->SetUserAction(new KnuclSteppingAction(anaManager));

  // Initialize G4 kernel
  runManager->Initialize();

  // Visualization manager construction
#if VISUALIZATION
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif
  for( int i=0; i<argc; i++ ){
    std::cout<<"argv["<<i<<"] : "<<argv[i]<<std::endl;
  }

  if(macro_name.empty())
  {
    G4UIsession* session =  new  G4UIterminal(new G4UItcsh);
    G4cout << "Session Start!!" << G4endl;
    session->SessionStart();

    delete session;  
  }
  else
  {
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    G4String command = "/control/execute ";
    //G4UIExecutive * ui = new G4UIExecutive(argc,argv); // for GUI
    UImanager->ApplyCommand(command+macro_name);
    //ui->SessionStart(); // for GUI
    //delete ui; // for GUI
  }
  // start a run
  //int numberOfEvent = 100;
  //runManager->BeamOn(numberOfEvent);

  // job termination
#if VISUALIZATION
  delete visManager;
#endif
  delete runManager;
  return 0;
}


void ModDecayBranch()
{
  G4ParticleTable *pTable = G4ParticleTable::GetParticleTable();
  G4VDecayChannel *mode;

  //### lambda(1520) decay mode ### //
  // default setting of G4 is wrong:
  //   G4DecayTable:  lambda(1520)
  //    0:  BR:  0.225  [Phase Space]   :   proton kaon-
  //    1:  BR:  0.225  [Phase Space]   :   neutron anti_kaon0
  //    2:  BR:  0.143333  [Phase Space]   :   sigma+ pi-
  //    3:  BR:  0.143333  [Phase Space]   :   sigma0 pi0
  //    4:  BR:  0.143333  [Phase Space]   :   sigma- pi+
  //    5:  BR:  0.0366667  [Phase Space]   :   sigma(1385)+ pi-
  //    6:  BR:  0.0366667  [Phase Space]   :   sigma(1385)0 pi0
  //    7:  BR:  0.0366667  [Phase Space]   :   sigma(1385)- pi+
  //    8:  BR:  0.01  [Phase Space]   :   lambda gamma
  const G4String name_L1520 = "lambda(1520)";
  // search in particle table
  G4ParticleDefinition *anInstance_L1520 = pTable->FindParticle(name_L1520);
  //create Decay Table
  G4DecayTable *table_L1520 = new G4DecayTable();
  // create a decay channel
  mode = new G4PhaseSpaceDecayChannel("lambda(1520)", 0.225, 2, "proton", "kaon-");
  table_L1520->Insert(mode);
  mode = new G4PhaseSpaceDecayChannel("lambda(1520)", 0.225, 2, "neutron", "anti_kaon0");
  table_L1520->Insert(mode);
  mode = new G4PhaseSpaceDecayChannel("lambda(1520)", 0.14, 2, "sigma+", "pi-");
  table_L1520->Insert(mode);
  mode = new G4PhaseSpaceDecayChannel("lambda(1520)", 0.14, 2, "sigma0", "pi0");
  table_L1520->Insert(mode);
  mode = new G4PhaseSpaceDecayChannel("lambda(1520)", 0.14, 2, "sigma-", "pi+");
  table_L1520->Insert(mode);
  mode = new G4PhaseSpaceDecayChannel("lambda(1520)", 0.06666, 3, "lambda", "pi+", "pi-");
  table_L1520->Insert(mode);
  mode = new G4PhaseSpaceDecayChannel("lambda(1520)", 0.03334, 3, "lambda", "pi0", "pi0");
  table_L1520->Insert(mode);
  mode = new G4PhaseSpaceDecayChannel("lambda(1520)", 0.005, 3, "sigma+", "pi+", "pi0");
  table_L1520->Insert(mode);
  mode = new G4PhaseSpaceDecayChannel("lambda(1520)", 0.005, 3, "sigma0", "pi0", "pi0");
  table_L1520->Insert(mode);
  mode = new G4PhaseSpaceDecayChannel("lambda(1520)", 0.005, 3, "sigma0", "pi+", "pi-");
  table_L1520->Insert(mode);
  mode = new G4PhaseSpaceDecayChannel("lambda(1520)", 0.005, 3, "sigma-", "pi-", "pi0");
  table_L1520->Insert(mode);
  mode = new G4PhaseSpaceDecayChannel("lambda(1520)", 0.01, 2, "lambda", "gamma");
  table_L1520->Insert(mode);
  // set the decay table
  anInstance_L1520->SetDecayTable(table_L1520);

  //### K*+ decay mode ### //
  // default setting of G4 is wrong:
  //   G4DecayTable:  k_star+
  //    0:  BR:  0.5  [Phase Space]   :   kaon+ pi0
  //    1:  BR:  0.5  [Phase Space]   :   kaon0 pi+
  const G4String name_K_star_p = "k_star+";
  G4ParticleDefinition *anInstance_K_star_p = pTable->FindParticle(name_K_star_p);
  G4DecayTable *table_K_star_p = new G4DecayTable();
  mode = new G4PhaseSpaceDecayChannel("k_star+", 0.66666, 2, "kaon0", "pi+");
  table_K_star_p->Insert(mode);
  mode = new G4PhaseSpaceDecayChannel("k_star+", 0.33334, 2, "kaon+", "pi0");
  table_K_star_p->Insert(mode);
  anInstance_K_star_p->SetDecayTable(table_K_star_p);

  //### K*0 decay mode ### //
  // default setting of G4 is wrong:
  //   G4DecayTable:  k_star0
  //   0:  BR:  0.5  [Phase Space]   :   kaon+ pi-
  //   1:  BR:  0.5  [Phase Space]   :   kaon0 pi0
  const G4String name_K_star_0 = "k_star0";
  G4ParticleDefinition *anInstance_K_star_0 = pTable->FindParticle(name_K_star_0);
  G4DecayTable *table_K_star_0 = new G4DecayTable();
  mode = new G4PhaseSpaceDecayChannel("k_star0", 0.66666, 2, "kaon+", "pi-");
  table_K_star_0->Insert(mode);
  mode = new G4PhaseSpaceDecayChannel("k_star0", 0.33334, 2, "kaon0", "pi0");
  table_K_star_0->Insert(mode);
  anInstance_K_star_0->SetDecayTable(table_K_star_0);

  //### K*0bar decay mode ### //
  // default setting of G4 is wrong:
  //   G4DecayTable:  anti_k_star0
  //   0:  BR:  0.5  [Phase Space]   :   kaon- pi+
  //   1:  BR:  0.5  [Phase Space]   :   anti_kaon0 pi0
  const G4String name_anti_K_star_0 = "anti_k_star0";
  G4ParticleDefinition *anInstance_anti_K_star_0 = pTable->FindParticle(name_anti_K_star_0);
  G4DecayTable *table_anti_K_star_0 = new G4DecayTable();
  mode = new G4PhaseSpaceDecayChannel("anti_k_star0", 0.66666, 2, "kaon-", "pi+");
  table_anti_K_star_0->Insert(mode);
  mode = new G4PhaseSpaceDecayChannel("anti_k_star0", 0.33334, 2, "anti_kaon0", "pi0");
  table_anti_K_star_0->Insert(mode);
  anInstance_anti_K_star_0->SetDecayTable(table_anti_K_star_0);

  //### K*- decay mode ### //
  // default setting of G4 is wrong:
  //    G4DecayTable:  k_star-
  //   0:  BR:  0.5  [Phase Space]   :   kaon- pi0
  //   1:  BR:  0.5  [Phase Space]   :   anti_kaon0 pi-
  const G4String name_K_star_m = "k_star-";
  G4ParticleDefinition *anInstance_K_star_m = pTable->FindParticle(name_K_star_m);
  G4DecayTable *table_K_star_m = new G4DecayTable();
  mode = new G4PhaseSpaceDecayChannel("k_star-", 0.66666, 2, "anti_kaon0", "pi-");
  table_K_star_m->Insert(mode);
  mode = new G4PhaseSpaceDecayChannel("k_star-", 0.33334, 2, "kaon-", "pi0");
  table_K_star_m->Insert(mode);
  anInstance_K_star_m->SetDecayTable(table_K_star_m);

#if 1
  // K+ decay is forbidden for beam study
  //### K+ decay mode ### //
  // default setting of G4 is unstable particle
  // G4DecayTable:  kaon+
  //  0:  BR:  0.6355  [Phase Space]   :   mu+ nu_mu
  //  1:  BR:  0.2066  [Phase Space]   :   pi+ pi0
  //  2:  BR:  0.0559  [Phase Space]   :   pi+ pi+ pi-
  //  3:  BR:  0.0507  [KL3 Decay]   :   pi0 e+ nu_e
  //  4:  BR:  0.0335  [KL3 Decay]   :   pi0 mu+ nu_mu
  //  5:  BR:  0.01761  [Phase Space]   :   pi+ pi0 pi0
  const G4String name_Kp = "kaon+";
  G4ParticleDefinition *anInstance_Kp = pTable->FindParticle(name_Kp);
  anInstance_Kp->SetPDGStable(true);
#endif

  G4ParticleDefinition *anInstance_L1405 = pTable->FindParticle("lambda(1405)");
  G4DecayTable *table_L1405 = anInstance_L1405-> GetDecayTable();
  mode = new G4PhaseSpaceDecayChannel("lambda(1405)", 0.0, 3, "pi+", "pi-", "neutron");
  table_L1405-> Insert(mode);
  anInstance_L1405-> SetDecayTable(table_L1405);
}

// ######################################################################
// ### K-pp ###
// ######################################################################
#include "G4VShortLivedParticle.hh"

class G4Kpp : public G4ParticleDefinition{
private:
  static G4Kpp* theInstance;
  G4Kpp() {}
  ~G4Kpp() {}
  
public:
  static G4Kpp* Definition(KnuclAnaManager* anaManager);
  static G4Kpp* KppDefinition(KnuclAnaManager* anaManager);
  static G4Kpp* Kpp(KnuclAnaManager* anaManager);
};

G4Kpp* G4Kpp::theInstance = 0;

G4Kpp* G4Kpp::Definition(KnuclAnaManager* anaManager) {

  if (theInstance != 0) return theInstance;
  const G4String name = "Kpp";

  // search in particle table
  G4ParticleTable *pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition *anInstance = pTable->FindParticle(name);

  G4double BindingEnergy = anaManager->GetBindingEnergy();
  G4double width         = anaManager->GetDecayWidth();

  G4double mass = pTable->FindParticle("kaon-")->GetPDGMass() +
    2*pTable->FindParticle("proton")->GetPDGMass() - BindingEnergy;
  G4double charge    = 1;
  G4int iSpin        = 0;
  G4int iParity      = -1;
  G4int iConjugation = 0; // temporal
  G4int iIsospin     = 1;
  G4int iIsospinZ    = 1; // temporal
  G4int gParity      = 0;
  const G4String& pType = "baryon";
  G4int lepton       = 0;
  G4int baryon       = 2;
  G4int encoding     = 9999; // new
  G4bool stable      = false;
  G4double lifetime  = 0;  // temporal
  G4DecayTable* decaytable = new G4DecayTable();
  G4bool shortlived  = true;
  const G4String&  subType = "";
  G4int anti_encoding = 0;

  if (anInstance == 0)
    {
      // create particle
      //
      //    Arguments for constructor are as follows
      //               name             mass          width         charge
      //             2*spin           parity  C-conjugation
      //          2*Isospin       2*Isospin3       G-parity
      //               type    lepton number  baryon number   PDG encoding
      //             stable         lifetime    decay table
      //             shortlived      subType    anti_encoding
      anInstance = new G4ParticleDefinition(
		    name, mass, width, charge,
		    iSpin, iParity, iConjugation,
		    iIsospin, iIsospinZ, gParity,
		    pType, lepton, baryon, encoding,
		    stable, lifetime, decaytable,
		    shortlived, subType, anti_encoding );
    }

  //create Decay Table
  G4DecayTable *table = new G4DecayTable();

  // create a decay channel
  G4VDecayChannel *mode;
  double BR_Lp   = 1.0;
  double BR_Sp   = 0;
  double BR_piSp = 0;
  double KppMassThreshold = 0;
  switch( anaManager->GetDecayMode() ){
  case 0: // Lambda+p
    BR_Lp   = 1.0;
    BR_Sp   = 0;
    BR_piSp = 0;
    KppMassThreshold = pTable->FindParticle("lambda")->GetPDGMass() +
      pTable->FindParticle("proton")->GetPDGMass();
    break;
  case 1: // Sigma0+p
    BR_Lp   = 0;
    BR_Sp   = 1.0;
    BR_piSp = 0;
    KppMassThreshold = pTable->FindParticle("sigma0")->GetPDGMass() +
      pTable->FindParticle("proton")->GetPDGMass();
    break;
  case 2: // (pi+Sigma)0+p
    BR_Lp   = 0;
    BR_Sp   = 0;
    BR_piSp = 0.33333;
    KppMassThreshold = pTable->FindParticle("pi+")->GetPDGMass() +
      pTable->FindParticle("sigma-")->GetPDGMass() +
      pTable->FindParticle("proton")->GetPDGMass();
    break;
  default: // Lambda+p
    break;
  }
  mode = new G4PhaseSpaceDecayChannel("Kpp", BR_Lp, 2, "lambda", "proton");
  table->Insert(mode);
  mode = new G4PhaseSpaceDecayChannel("Kpp", BR_Sp, 2, "sigma0", "proton");
  table->Insert(mode);
  mode = new G4PhaseSpaceDecayChannel("Kpp", BR_piSp, 3, "pi0", "sigma0", "proton");
  table->Insert(mode);
  mode = new G4PhaseSpaceDecayChannel("Kpp", BR_piSp, 3, "pi+", "sigma-", "proton");
  table->Insert(mode);
  mode = new G4PhaseSpaceDecayChannel("Kpp", BR_piSp, 3, "pi-", "sigma+", "proton");
  table->Insert(mode);
  anInstance->SetDecayTable(table);

  anaManager->SetKppMassThreshold(KppMassThreshold);

  // Print K-pp information //
#if 1
  std::cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<std::endl;
  std::cout<<"@@@  === K-pp is generated ==="<<std::endl;
  std::cout<<"@@@  mass   = "<<mass<<" MeV"<<std::endl;
  std::cout<<"@@@  B.E    = "<<BindingEnergy<<" MeV "<<std::endl;
  std::cout<<"@@@  width  = "<<width<<" MeV"<<std::endl;
  std::cout<<"@@@  BR     = "<<BR_Lp<<" : lambda + proton"<<std::endl;
  std::cout<<"@@@           "<<BR_Sp<<" : sigma0 + proton"<<std::endl;
  std::cout<<"@@@           "<<BR_piSp<<" : pi0 + sigma0 + proton"<<std::endl;
  std::cout<<"@@@           "<<BR_piSp<<" : pi+ + sigma- + proton"<<std::endl;
  std::cout<<"@@@           "<<BR_piSp<<" : pi- + sigma+ + proton"<<std::endl;
  std::cout<<"@@@ mass th = "<<KppMassThreshold<<" MeV "<<std::endl;
  std::cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<std::endl;
#endif

  theInstance = reinterpret_cast<G4Kpp*>(anInstance);
  return theInstance;
}

G4Kpp *G4Kpp::KppDefinition(KnuclAnaManager* anaManager)
{
  return Definition(anaManager);
}

G4Kpp *G4Kpp::Kpp(KnuclAnaManager* anaManager)
{
  return Definition(anaManager);
}

void AddKpp(KnuclAnaManager* anaManager)
{
  G4Kpp::KppDefinition(anaManager);
}

// ######################################################################
// ### K-pn ###
// ######################################################################
class G4Kpn : public G4ParticleDefinition{
private:
  static G4Kpn* theInstance;
  G4Kpn() {}
  ~G4Kpn() {}
  
public:
  static G4Kpn* Definition(KnuclAnaManager* anaManager);
  static G4Kpn* KpnDefinition(KnuclAnaManager* anaManager);
  static G4Kpn* Kpn(KnuclAnaManager* anaManager);
};

G4Kpn* G4Kpn::theInstance = 0;

G4Kpn* G4Kpn::Definition(KnuclAnaManager* anaManager) {

  if (theInstance != 0) return theInstance;
  const G4String name = "Kpn";

  // search in particle table
  G4ParticleTable *pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition *anInstance = pTable->FindParticle(name);

  G4double BindingEnergy = anaManager->GetBindingEnergy();
  G4double width         = anaManager->GetDecayWidth();

  G4double mass = pTable->FindParticle("kaon-")->GetPDGMass() +
    pTable->FindParticle("proton")->GetPDGMass() +
    pTable->FindParticle("neutron")->GetPDGMass() - BindingEnergy;
  G4double charge    = 0;
  G4int iSpin        = 0;
  G4int iParity      = -1;
  G4int iConjugation = 0; // temporal
  G4int iIsospin     = 1;
  G4int iIsospinZ    = 1; // temporal
  G4int gParity      = 0;
  const G4String& pType = "baryon";
  G4int lepton       = 0;
  G4int baryon       = 2;
  G4int encoding     = 9998; // new
  G4bool stable      = false;
  G4double lifetime  = 0;  // temporal
  G4DecayTable* decaytable = new G4DecayTable();
  G4bool shortlived  = true;
  const G4String&  subType = "";
  G4int anti_encoding = 0;

  if (anInstance == 0)
    {
      // create particle
      //
      //    Arguments for constructor are as follows
      //               name             mass          width         charge
      //             2*spin           parity  C-conjugation
      //          2*Isospin       2*Isospin3       G-parity
      //               type    lepton number  baryon number   PDG encoding
      //             stable         lifetime    decay table
      //             shortlived      subType    anti_encoding
      anInstance = new G4ParticleDefinition(
		    name, mass, width, charge,
		    iSpin, iParity, iConjugation,
		    iIsospin, iIsospinZ, gParity,
		    pType, lepton, baryon, encoding,
		    stable, lifetime, decaytable,
		    shortlived, subType, anti_encoding );
    }

  //create Decay Table
  G4DecayTable *table = new G4DecayTable();

  // create a decay channel
  // --- Lambda+n decay mode only !!! --- //
  G4VDecayChannel *mode;
  mode = new G4PhaseSpaceDecayChannel("Kpn", 1.0, 2, "lambda", "neutron");
  table->Insert(mode);
  anInstance->SetDecayTable(table);

  // Print K-pp information //
#if 1
  std::cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<std::endl;
  std::cout<<"@@@  === K-pn is generated ==="<<std::endl;
  std::cout<<"@@@  mass  = "<<mass<<" MeV"<<std::endl;
  std::cout<<"@@@  ( B.E = "<<BindingEnergy<<" MeV )"<<std::endl;
  std::cout<<"@@@  width = "<<width<<" MeV"<<std::endl;
  std::cout<<"@@@  BR    = "<<1.0<<" : lambda + neutron"<<std::endl;
  std::cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<std::endl;
#endif

  theInstance = reinterpret_cast<G4Kpn*>(anInstance);
  return theInstance;
}

G4Kpn *G4Kpn::KpnDefinition(KnuclAnaManager* anaManager)
{
  return Definition(anaManager);
}

G4Kpn *G4Kpn::Kpn(KnuclAnaManager* anaManager)
{
  return Definition(anaManager);
}

void AddKpn(KnuclAnaManager* anaManager)
{
  G4Kpn::KpnDefinition(anaManager);
}
