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
// $Id: KnuclPrimaryGeneratorAction.cc,v 1.4 2015/06/18 01:46:36 sakuma Exp $
// GEANT4 tag $Name:  $
//

#include "KnuclPrimaryGeneratorAction.hh"
#include "KnuclHist.hh"
#include "KnuclAnaManager.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4TransportationManager.hh"
#include "G4Navigator.hh"
#include "G4VPhysicalVolume.hh"
#include "G4UImanager.hh"
#include "globals.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4LorentzVector.hh"
#include "G4LorentzConvertor.hh"
#include "Randomize.hh"
#include <math.h>
#include <stdlib.h>

#include "KnuclCommon.h"
#include "ComCrossSectionTable.hh" 

#define FLAT_ANGL_DIST 0

// ---------------------------------- //
// beam profile from Noumi-san K1.8BR
// ---------------------------------- //
// position @ FF
const G4double beam_z =   0.0  *cm;
// const G4double beam_x =  -0.009*cm;
// const G4double beam_y =   0.024*cm;
// const G4double beam_dx = 0.588*cm;
// const G4double beam_dy = 0.292*cm;
// // direction @ FF (mrad)
// const G4double beam_xp = -0.972*mrad;
// const G4double beam_yp =  0.451*mrad;
// const G4double beam_dxp = 22.676*mrad;
// const G4double beam_dyp = 2.572*mrad;

// from run49c BPC track
// position @ FF
const G4double beam_x =  -0.475*cm;
const G4double beam_y =   0.081*cm;
const G4double beam_dx =  1.971*cm;
const G4double beam_dy =  1.907*cm;
// direction @ FF (mrad)
const G4double beam_xp =   -4.23*mrad;
const G4double beam_yp =  -0.771*mrad;
const G4double beam_dxp =  18.53*mrad;
const G4double beam_dyp =  7.808*mrad;


//////////////////////////////////////////////////////
KnuclPrimaryGeneratorAction::KnuclPrimaryGeneratorAction(KnuclAnaManager* ana)
  : anaManager(ana)
//////////////////////////////////////////////////////
{
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);
  gen = new TGenPhaseSpace();
  reactionData = new ReactionData();
  //temporary using "K3He" flag for E31 mc.from Knucl:ProcessID
  if (anaManager->GetProcessID()==KnuclAnaManager::K3He || 
      anaManager->GetProcessID()==KnuclAnaManager::Single
  ){
    csTable  = new CrossSectionTable(ana->GetCSFileName(), ana->GetBeamMomentum()*GeV);
    PrintAllCS();
    ana->SetCSTable(*csTable);
  }
  int csID = csTable->CS(0).Id();
  std::cout << __FILE__ << " L." << __LINE__ << " csID: " << csID << std::endl;
  TFile *genfile = NULL;
  MakeUniformInqmass = anaManager->GetUniformGenFlag();
  if(MakeUniformInqmass){
    std::cout << "Make Uniform in q vs mass " << std::endl;
    if(csID == 1725)   genfile = new TFile("probSp.root","READ");
    if(csID == 1525)   genfile = new TFile("probSm.root","READ");
    if(csID == 1600)   genfile = new TFile("probLpim.root","READ");
    if(csID == 2006)   genfile = new TFile("probnpipiL.root","READ");
    //G4cout << "File name for making uniform distribution : " << genfile->GetName() << G4endl;
    if(genfile)h2genprob = (TH2D*) genfile->Get("h2prob");
    if(genfile)h2genprob->Print();
  }
  //gErrorIgnoreLevel = 5000;
  // fermi motion with 
  std::string str = "([0]*exp(-x*x/([1]*[1]))+[2]*exp(-x*x/([3]*[3]))+[4]*exp(-x*x/([5]*[5])))*x*x";
  const char* fun = str.c_str();
  const double min = 0;
  const double max = 1000;
  FermiMotion_3HeTwoBody_dist   = new TF1("FermiMotion_3HeTwoBody_dist"  , fun, min, max);
  FermiMotion_3HeThreeBody_dist = new TF1("FermiMotion_3HeThreeBody_dist", fun, min, max);
  FermiMotion_deuteron_dist     = new TF1("FermiMotion_deuteron_dist"    , fun, min, max);
  for ( int i=0; i<6; i++ ){
    FermiMotion_3HeTwoBody_dist->SetParameter(i, FermiMotion_3HeTwoBody[i]);
    FermiMotion_3HeThreeBody_dist->SetParameter(i, FermiMotion_3HeThreeBody[i]);
    FermiMotion_deuteron_dist->SetParameter(i, FermiMotion_deuteron[i]);
  }
}

//////////////////////////////////////////////////////
KnuclPrimaryGeneratorAction::~KnuclPrimaryGeneratorAction()
//////////////////////////////////////////////////////
{
  delete particleGun;
  delete gen;
  delete reactionData;
  if (anaManager->GetProcessID()==KnuclAnaManager::K3He ){
    delete csTable;
  }

  // fermi motion
  delete FermiMotion_3HeTwoBody_dist;
  delete FermiMotion_3HeThreeBody_dist;
  delete FermiMotion_deuteron_dist;
}

//////////////////////////////////////////////////////
void KnuclPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
//////////////////////////////////////////////////////
{
  G4double beam_mom = anaManager->GetBeamMomentum()*GeV; // beam momentum

  anaManager->SetReactionID(-1);
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;

  G4double mass;
  if (anaManager->GetProcessID()==KnuclAnaManager::Single ||
      anaManager->GetProcessID()==KnuclAnaManager::Beam ||
      anaManager->GetProcessID()==KnuclAnaManager::neutron_efficiency_mode ){
    TString name = anaManager->GetSingleParticle();
    if( '0'<=name[1] && name[1]<='9' ){
      mass = particleTable->FindParticle(atoi(name.Data()))->GetPDGMass();
      particleGun->SetParticleDefinition(particleTable->FindParticle(atoi(name.Data())));
    }else{
      mass = particleTable->FindParticle(name.Data())->GetPDGMass();
      particleGun->SetParticleDefinition(particleTable->FindParticle(name.Data()));
    }
  } 
  else if (anaManager->GetProcessID()==KnuclAnaManager::K3He ){ // for reversed beam
    mass = particleTable->FindParticle("kaon+")->GetPDGMass();
    particleGun->SetParticleDefinition(particleTable->FindParticle(particleName="kaon+"));
  }
  
  if (anaManager->GetProcessID()==KnuclAnaManager::neutron_efficiency_mode){
    particleGun->SetParticleDefinition(particleTable->FindParticle(particleName="neutron"));
    mass = particleTable->FindParticle("neutron")->GetPDGMass();
  }

  G4double beam_mom_bite = 0;
  G4double BeamMomentum = 0;
  G4double beam_ene = 0;

  G4double vtx0 = 0;
  G4double vty0 = 0;

  G4double xdir = 0;
  G4double ydir = 0;
  G4double zdir = 0;

  G4double ax = 0; G4double bx = 0;
  G4double ay = 0; G4double by = 0;
 
  G4double vtx = 0;
  G4double vty = 0;
  G4double vtz = 0;

  G4double tdir = 0;
  G4double theta = 0;
  G4double phi = 0;
  G4double px = 0;
  G4double py = 0;
  G4double pz = 0;
  
  /*
  if (anaManager->GetProcessID()==KnuclAnaManager::neutron_efficiency_mode){
    //  if (1){
    beam_ene = sqrt(beam_mom*beam_mom + mass*mass);
    vtx = 0.0;
    vty = 0.0;
    vtz = 10.0;
    px = 0.0;
    py = 0.0;
    pz = beam_mom;

    particleGun->SetParticlePosition(G4ThreeVector(vtx, vty, vtz));
    particleGun->SetParticleMomentumDirection(G4ThreeVector(px, py, pz));
    particleGun->SetParticleTime(0.0);
    particleGun->SetParticleEnergy(beam_ene-mass);
    particleGun->GeneratePrimaryVertex(anEvent);
  }*/
  
  if (anaManager->GetProcessID()==KnuclAnaManager::Single ||
      anaManager->GetProcessID()==KnuclAnaManager::neutron_efficiency_mode ){
    BeamMomentum = CLHEP::RandFlat::shoot(anaManager->GetSingleMomentumMIN(),
					  anaManager->GetSingleMomentumMAX())*GeV;
    beam_ene = sqrt(BeamMomentum*BeamMomentum + mass*mass);

    vtx = CLHEP::RandFlat::shoot(anaManager->GetSingleVertexMIN().x(),
				 anaManager->GetSingleVertexMAX().x())*mm;
    vty = CLHEP::RandFlat::shoot(anaManager->GetSingleVertexMIN().y(),
				 anaManager->GetSingleVertexMAX().y())*mm;
    vtz = CLHEP::RandFlat::shoot(anaManager->GetSingleVertexMIN().z(),
				 anaManager->GetSingleVertexMAX().z())*mm;
    
    G4ThreeVector mom;
    theta = CLHEP::RandFlat::shoot(anaManager->GetSingleThetaMIN(),
				   anaManager->GetSingleThetaMAX());
    phi = 2.0*pi*G4UniformRand()*rad;
    mom.setRThetaPhi(BeamMomentum, acos(theta), phi);

    particleGun->SetParticlePosition(G4ThreeVector(vtx, vty, vtz));
    particleGun->SetParticleMomentumDirection(mom);
    particleGun->SetParticleTime(0.0);
    particleGun->SetParticleEnergy(beam_ene-mass);
    particleGun->GeneratePrimaryVertex(anEvent);

    //std::cerr<<anaManager->GetSingleMomentumMIN()<<std::endl;
    //std::cout<<"(vtx, vty, vtz) = ("<<vtx<<" , "<<vty<<" , "<<vtz<<")"<<std::endl;
    //std::cout<<"(px, py, pz) = ("<<mom.x()<<" , "<<mom.y()<<" , "<<" , "<<mom.z()<<")"<<std::endl;
    //std::cout<<"beam_ene-mass = "<<beam_ene-mass<<std::endl;
    // getchar();
  
    reactionData->Init();
    reactionData->SetPDG(particleGun->GetParticleDefinition()->GetPDGEncoding());
    reactionData->SetParticle( mom.x(), mom.y(), mom.z(), mass );
    anaManager->setHistReactionData(reactionData);
  
  
  }
  else if (anaManager->GetProcessID()==KnuclAnaManager::Beam ||
	   anaManager->GetProcessID()==KnuclAnaManager::K3He ){
    beam_mom_bite = beam_mom*BEAM_MOM_BITE/100./2.35; // FWHM->sigma
    BeamMomentum  = G4RandGauss::shoot(beam_mom, beam_mom_bite);
    beam_ene = sqrt(BeamMomentum*BeamMomentum + mass*mass);
    
    if (anaManager->GetProcessID()==KnuclAnaManager::Beam) {      
      vtx0 = G4RandGauss::shoot(beam_x,beam_dx);
      vty0 = G4RandGauss::shoot(beam_y,beam_dy);      
      xdir = sin(G4RandGauss::shoot(beam_xp,beam_dxp)/rad);
      ydir = sin(G4RandGauss::shoot(beam_yp,beam_dyp)/rad);
#if 0 // ideal condition
      BeamMomentum = beam_mom;
      beam_ene = sqrt(BeamMomentum*BeamMomentum + mass*mass);
      vtx0 = 0;
      vty0 = 0;
      xdir = 0;
      ydir = 0;
#endif
      zdir = sqrt(1.0-xdir*xdir-ydir*ydir);
      vtz = -1.1*m; // before just T0
      ax = xdir/zdir; bx = vtx0;  
      ay = ydir/zdir; by = vty0;        
      vtx = ax*vtz + bx; 
      vty = ay*vtz + by; 
    }else if(anaManager->GetProcessID()==KnuclAnaManager::K3He) {           
      G4TransportationManager *transManager 
	= G4TransportationManager::GetTransportationManager();
      G4Navigator* nav=transManager->GetNavigatorForTracking();
      bool FIDUCIAL=false;
      int count=0;
      do{
	vtx0 = G4RandGauss::shoot(beam_x,beam_dx);
	vty0 = G4RandGauss::shoot(beam_y,beam_dy);      
	xdir = sin(G4RandGauss::shoot(beam_xp,beam_dxp)/rad);
	ydir = sin(G4RandGauss::shoot(beam_yp,beam_dyp)/rad);
	zdir = sqrt(1.0-xdir*xdir-ydir*ydir);
	vtz = ((G4UniformRand()-0.5)*2.0)*(anaManager->GetTargetLength())/2.0*mm;
	ax = xdir/zdir; bx = vtx0;  
	ay = ydir/zdir; by = vty0;        
	vtx = ax*vtz + bx; 
	vty = ay*vtz + by; 
	G4VPhysicalVolume* vol=nav->LocateGlobalPointAndSetup(G4ThreeVector(vtx,vty,vtz));
	//	std::cout<<vtx<<"  "<<vty<<"  "<<vtz<<"  "<<vol->GetName()<<std::endl;
	if(vol->GetName()=="Fiducial") FIDUCIAL=true;
	count++;
	if(count>1000){
	  std::cout<<"maximum loops in defining vertex position !!!" <<std::endl;
	  break;
	}
      } while (!FIDUCIAL);
    }
    tdir  = sqrt(xdir*xdir+ydir*ydir); 
    theta = atan2(tdir,zdir);
    phi   = atan2(ydir,xdir);
    px = BeamMomentum*sin(theta)*cos(phi);
    py = BeamMomentum*sin(theta)*sin(phi);
    pz = BeamMomentum*cos(theta);

    if (anaManager->GetBeamMomentum()){
      particleGun->SetParticlePosition(G4ThreeVector(vtx, vty, vtz));
      if (anaManager->GetProcessID()==KnuclAnaManager::Beam) {
	particleGun->SetParticleMomentumDirection(G4ThreeVector(px, py, pz));
      }else {
	particleGun->SetParticleMomentumDirection(G4ThreeVector(-px, -py, -pz));
      }
      particleGun->SetParticleTime(0.0);
      particleGun->SetParticleEnergy(beam_ene-mass);
      particleGun->GeneratePrimaryVertex(anEvent);    
    }

  }

  //==================================================================================
  // reaction vertex position is now stored for analysis
  //==================================================================================
  anaManager->SetBeamDirection (xdir,ydir,zdir);
  anaManager->SetVertexPosition(vtx, vty, vtz );
  anaManager->SetBeamMomVector (px,  py,  pz  );
  anaManager->SetBeamEne (beam_ene);
  anaManager->SetBeamMass (mass);

  anaManager->SetReactionID(-999);
  //==================================================================================
  // reaction
  //==================================================================================

  int reacID = 0;

  if (anaManager->GetProcessID()==KnuclAnaManager::neutron_efficiency_mode) anaManager->SetReactionID(-9);

  if (anaManager->GetProcessID()==KnuclAnaManager::Beam) anaManager->SetReactionID(0);

  if (anaManager->GetProcessID()==KnuclAnaManager::Single) anaManager->SetReactionID(1);

  if (anaManager->GetProcessID()==KnuclAnaManager::ReadFromFile) reacID = ReadFromFile(anEvent);

  if (anaManager->GetProcessID()==KnuclAnaManager::K3He)   reacID = Kminus3He(anEvent); 

  //--- output data-file ---//
  if( anaManager->GetOutputDataFile() ){
    WriteOutputOscarFile(anEvent, reacID);
  }

}


// %%%%%%%%%%%%%%%%%%%% //
// % common functions % //
// %%%%%%%%%%%%%%%%%%%% //


//////////////////////////////////////////////////////
float KnuclPrimaryGeneratorAction::Legendre(int n, float x){
//////////////////////////////////////////////////////
  float f = 0.0;

  if (n==0){
    f = 1.0;
  } else
  if (n==1){
    f = x;
  } else
  if (n==2){
    f = 0.5*(3.0*x*x-1.0);
  } else
  if (n==3){
    f = 0.5*(5.0*x*x*x-3.0*x);
  } else
  if (n==4){
    f = 0.125*(35.0*x*x*x*x-30.0*x*x + 3.0);
  } else
  if (n==5){
    f = 0.125*(63.0*x*x*x*x*x-70.0*x*x*x+15.0*x);
  } else
  if (n==6){
    f = 0.0625*(231.0*x*x*x*x*x*x - 315.0*x*x*x*x + 105.0*x*x - 5.0);
  } else
  if (n==7){
    f = 0.0625*(429.0*x*x*x*x*x*x*x - 693.0*x*x*x*x*x + 315.0*x*x*x - 35.0*x);
  } else
  if (n==8){
    f = 0.0078125*(6435.0*x*x*x*x*x*x*x*x - 12012.0*x*x*x*x*x*x + 6930.0*x*x*x*x - 1260.0*x*x + 35);
  }

 else {
    G4cout << "n > 8 is not supported!!!" << G4endl;
  }
  return f;
}

//////////////////////////////////////////////////////
bool KnuclPrimaryGeneratorAction::FermiMom_judge(G4int tgtID, G4int nbodyFlag, G4double mom){
//////////////////////////////////////////////////////
// nbodyFlag 2: two-body, 3: three-body in 3He

  if ( mom < 0 ) return true; // for num of spectators = 1 case

  G4double param[6];
  G4double max;
  if( tgtID==1000020030 && nbodyFlag==2 ){
    for (G4int i=0; i<6; i++){ param[i] = FermiMotion_3HeTwoBody[i]; }
    max = 700.;
  }else if( tgtID==1000020030 && nbodyFlag==3 ){
    for (G4int i=0; i<6; i++){ param[i] = FermiMotion_3HeThreeBody[i]; }
    max = 210.;
  }else if( tgtID==1000010020 ){
    for (G4int i=0; i<6; i++){ param[i] = FermiMotion_deuteron[i]; }
    max=1400.;
  }
    
  G4double density = param[0]*exp(-mom*mom/(param[1]*param[1]))
    +param[2]*exp(-mom*mom/(param[3]*param[3]))
    +param[4]*exp(-mom*mom/(param[5]*param[5]));
  G4double ran = max*G4UniformRand();
  //std::cout<<nbodyFlag<<" "<<mom<<" "<<density<<" "<<ran<<std::endl;  

  return ( ran<density ) ? true : false;
}

//////////////////////////////////////////////////////
double KnuclPrimaryGeneratorAction::FermiMom_gen(G4int tgtID, G4int nbodyFlag){
//////////////////////////////////////////////////////
// nbodyFlag 2: two-body, 3: three-body in 3He

  G4double mom;
  if( tgtID==1000020030 && nbodyFlag==2 ){
    mom = FermiMotion_3HeTwoBody_dist->GetRandom()*MeV;
  }else if( tgtID==1000020030 && nbodyFlag==3 ){
    mom = FermiMotion_3HeThreeBody_dist->GetRandom()*MeV;
  }else if( tgtID==1000010020 ){
    mom = FermiMotion_deuteron_dist->GetRandom()*MeV;
  }

  return mom;
}

//////////////////////////////////////////////////////
int KnuclPrimaryGeneratorAction::ReadFromFile(G4Event* anEvent){
//////////////////////////////////////////////////////

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  int      Pid            = 50;

  G4int aInt; G4float aFloat; G4int Mult;
  G4float reacID = 0;

  ifstream* InputFile = anaManager->GetInputData();
  if ( ! (InputFile->fail()) && ! InputFile->eof() ){
    ( *InputFile ) >> aInt ;
    ( *InputFile ) >> Mult;
    ( *InputFile ) >> reacID ;
    ( *InputFile ) >> aFloat ;

    for(Int_t iPart=0; iPart < Mult; iPart++ )
    {
      Int_t nr;
      Int_t id;
      Float_t px,py,pz,E,ma,ptot;
      Float_t x,y,z,t;
      ( *InputFile ) >> nr ;
      ( *InputFile ) >> id ;
      ( *InputFile ) >> px >> py >> pz >> E >> ma ;
      ( *InputFile ) >> x >> y >> z >> t;

      ptot = sqrt(px*px+py*py+pz*pz);

      G4ParticleDefinition *partDef = particleTable->FindParticle(id);
      if (!partDef) continue;
      //G4double mass = partDef->GetPDGMass();
      G4double mass = ma;
      particleGun->SetParticleDefinition(partDef);
      particleGun->SetParticleMomentumDirection(G4ThreeVector(px/ptot,py/ptot,pz/ptot));
      particleGun->SetParticleEnergy(E*GeV-mass*GeV);
      particleGun->SetParticlePosition(G4ThreeVector(x, y, z));
      particleGun->GeneratePrimaryVertex( anEvent );
    }
    anaManager->SetReactionID(Pid);
  }
  return int(reacID);
}



// %%%%%%%%%%%%%%%%%%%% //
// % K-3He study      % //
// %%%%%%%%%%%%%%%%%%%% //

//////////////////////////////////////////////////////
int KnuclPrimaryGeneratorAction::Kminus3He(G4Event* anEvent){
//////////////////////////////////////////////////////
  return KminusReac(anEvent, GenerateReac());
}

//////////////////////////////////////////////////////
int KnuclPrimaryGeneratorAction::KminusReac(G4Event* anEvent, const CrossSection& cs)
//////////////////////////////////////////////////////
{
  G4bool FermiFlag = ( !anaManager->GetFermiMotion() || !cs.SpecPdgSize() ) ? false : true;
  G4int  FermiMode = anaManager->GetFermiMotionMode();

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  //--- initial state particles ---//
  G4ParticleDefinition *Init1 = particleTable->FindParticle(cs.InitPdg(0));
  G4ParticleDefinition *Init2 = FermiFlag ? 
    particleTable->FindParticle(cs.InitPdg(1)) : particleTable->FindParticle(cs.InitPdg(2));
  G4double MassInit1 = Init1->GetPDGMass();
  G4double MassInit2 = Init2->GetPDGMass();
  G4int tgtID = Init2->GetPDGEncoding();

  //--- initial conditions ---//
  G4ThreeVector beam_mom = anaManager->GetBeamMom();
  G4ThreeVector tgt_mom  = G4ThreeVector(0.0, 0.0, 0.0);
  G4LorentzVector beam   = G4LorentzVector(beam_mom, sqrt(beam_mom.mag2()+MassInit1*MassInit1));
  G4LorentzVector tgt    = G4LorentzVector(tgt_mom, sqrt(tgt_mom.mag2()+MassInit2*MassInit2));
  G4ThreeVector boost    = (beam+tgt).boostVector(); // K- + p/3He frame
  G4double CMmass        = (beam+tgt).m();

  //--- final state particles ---//
  G4int nFinl = cs.FinlPdgSize();
  G4int nSpec = cs.SpecPdgSize();
  G4int nBody = FermiFlag ? nFinl+nSpec : nFinl;
  std::vector <G4ParticleDefinition*> particle;
  for( int i=0; i<nFinl; i++ ){
    particle.push_back(particleTable->FindParticle(cs.FinlPdg(i)));
  }
  for( int i=0; i<nSpec; i++ ){
    particle.push_back(particleTable->FindParticle(cs.SpecPdg(i)));
  }

  //--- set parameter of angular distribution ---//
  G4double pol[9] = {0,0,0,0,0,0,0,0,0};
  for(int i=0; i<cs.PolSize(); i++){ pol[i] = cs.Pol(i); }
  G4double fermimom1 = 0;
  G4double fermimom2 = 0;

  //--- for infinite-loop break ---//
 START:
  //--- for infinite-loop break ---//

  //--- many body decay ---//
  G4double      mass[10];
  G4ThreeVector vec[10];
  
  /*
  while( true ){
    G4double totMass=0;
    G4double totGamma=0;
    for( int i=0; i<nFinl; i++ ){
      //should not sum ?
      totMass=particle[i]->GetPDGMass();
      totGamma=particle[i]->GetPDGWidth();	
    }
    //CMmass = beam + tgt
    if( CMmass>totMass ) break;
    else{
      if( (totMass-CMmass)/totGamma>2.0 ){
	std::cout<<"total mass : "<<totMass<<"  W : "<<CMmass<<"  total Gamma : "<<totGamma<<std::endl;
	beam_mom = anaManager->GetBeamMom();
	tgt_mom  = G4ThreeVector(0.0, 0.0, 0.0);
	beam   = G4LorentzVector(beam_mom, sqrt(beam_mom.mag2()+MassInit1*MassInit1));
	tgt    = G4LorentzVector(tgt_mom, sqrt(tgt_mom.mag2()+MassInit2*MassInit2));
	boost    = (beam+tgt).boostVector(); // K- + p/3He frame
	CMmass        = (beam+tgt).m();
      }
      else break;
    }
  }*/

  while( true ){
    G4double sum = 0;
    for( int i=0; i<nFinl+nSpec; i++ ){
      mass[i] = GetMass(*particle[i]);
      vec[i]  = G4ThreeVector(0,0,0);
      if( i<nBody ) sum += mass[i];
    }
    if( sum<CMmass ){ break; }
  }

  //-------------------------------------//
  //--- K-pp shape from the root-file ---//
  //-------------------------------------//
  if( anaManager->GetKppShape() && cs.Id()==3000 ){ // NO spectators
    if( !ManyBody(CMmass, nBody, mass, vec) ) goto START;    
    double P0 = vec[0].mag();
    double P1 = vec[1].mag();
    vec[0] = P0*G4ThreeVector(0,0,1);  // neutron
    vec[1] = P1*G4ThreeVector(0,0,-1); // Kpp
  }
  //--------------------//
  //--- Fermi motion ---//
  //--------------------//
  //asano memo
  //if taking into account Fermi motion
  else if( FermiFlag ){
    if( FermiMode ){
      //==============================//
      //===== on-shell treatment =====//
      //==============================//
      G4LorentzVector sp1 = G4LorentzVector(), sp2 = G4LorentzVector();
      G4double P_sp1 = -1, P_sp2 = -1;
      G4int loop = 0;
      while( true ){
#if FLAT_ANGL_DIST // w/o angular distributon for DEBUG
	if( !ManyBody(CMmass, nBody, mass, vec) ) goto START;
#else // w/ angular distributon
	if( !ManyBody(CMmass, nBody, mass, vec, nFinl, cs.PolSize(), pol, cs.PolMax()) ) goto START;
#endif
	// 3He case only
	sp1.setVectM(vec[nFinl], mass[nFinl]);
	sp1.boost(boost);
	P_sp1 = sp1.rho();
	if ( nSpec == 2 ){
	  sp2.setVectM(vec[nFinl+1], mass[nFinl+1]);
	  sp2.boost(boost);
	  P_sp2 = sp2.rho();
	}
	fermimom1 = P_sp1;
	fermimom2 = P_sp2;
	loop++;
	if ( nSpec == 1 ){
	  if( FermiMom_judge(tgtID, nSpec+1, sp1.rho()) ) break;
	} else {
	  G4LorentzVector sum  = sp1+sp2;
	  G4LorentzVector diff = sp1-sp2;
	  static const double scale = 1.7; // artificial factor
	  if( FermiMom_judge(tgtID, nSpec+1, sum.rho()) && FermiMom_judge(tgtID, nSpec+1, diff.rho()/scale) ) break;
	}
	if( loop >10000000 ){
	  std::cerr<<"KnuclPrimaryGeneratorAction::KminusReac, break infinite loop (num of loop = "
		   <<loop<<"), reactionID = "<<cs.Id()<<std::endl;
	  goto START;
	}
      }//while
    } else { // if( FermiMode ){
      //===============================//
      //===== off-shell treatment =====//
      //===============================//
      G4ThreeVector vec_sum = G4ThreeVector(0);
      if( nSpec == 1 ){
	fermimom1 = FermiMom_gen(tgtID, nSpec+1);
	vec[nFinl]  = fermimom1*RandUnitVec(); // lab frame
	vec_sum = vec[nFinl];
      }else if( nSpec == 2 ){
	static const double scale = 98.7/75.4; // artificial factor
	while( true ){
	  fermimom1 = FermiMom_gen(tgtID, nSpec+1);
	  fermimom2 = FermiMom_gen(tgtID, nSpec+1);
	  vec[nFinl]   = fermimom1*RandUnitVec(); // lab frame
	  vec[nFinl+1] = fermimom2*RandUnitVec(); // lab frame
	  vec_sum = vec[nFinl]+vec[nFinl+1];
	  if( FermiMom_judge(tgtID, nSpec+1, vec_sum.mag()) ){
	    fermimom1 *= scale;
	    fermimom2 *= scale;
	    vec[nFinl]   *= scale;
	    vec[nFinl+1] *= scale;
	    break;
	  }
	}
      }//if nSpec ==2
      
      G4double spec_ene[2] = {0,0};
      G4double ene_sum = 0;
      G4LorentzVector spec[2] = {G4LorentzVector(0), G4LorentzVector(0)};
      for( int i=0; i<nSpec; i++ ){
	spec_ene[i] = sqrt(vec[nFinl+i].mag2()+mass[nFinl+i]*mass[nFinl+i]);
	ene_sum    += spec_ene[i];
	spec[i]     = G4LorentzVector(spec_ene[i], vec[nFinl+i]);
      }
      
      tgt_mom = -1*vec_sum;
      Init2 = particleTable->FindParticle(cs.InitPdg(2));
      G4double tgt_ene = particleTable->FindParticle(cs.InitPdg(1))->GetPDGMass() - ene_sum; // = 3He is stopped
      tgt = G4LorentzVector(tgt_ene, tgt_mom);
      
      boost = (beam+tgt).boostVector(); // K- + p frame
      G4double CMmass_off = (beam+tgt).m();
      if( !ManyBody(CMmass_off, nFinl, mass, vec, nFinl, cs.PolSize(), pol, cs.PolMax()) ) goto START;
      
      for( int i=0; i<nSpec; i++ ){
	spec[i].boost(-1*boost);
	vec[nFinl+i] = spec[i].vect();
      }
    } // if( FermiMode ){
  } else{ // w/o Fermi motion
    if( !ManyBody(CMmass, nFinl, mass, vec, nFinl, cs.PolSize(), pol, cs.PolMax()) ) goto START;
  }


  //------------------------------------------------//
  //--- 2-step reactions (NO momentm dependence) ---//
  //------------------------------------------------//
  if( anaManager->GetTwoStep() && cs.SpecPdgSize() ){
    int mode_ts = anaManager->GetTwoStepMode();
    //--- initializtion ---//
    bool flag_ts = false;
    int in1_id  = -1; // initial particle (-1:temporal value)
    //*** CAUTION: fixed to the first filled spectator ***//
    int in2_id  = 0;  // target particle in spectator(s)
    //*** CAUTION: fixed to the first filled spectator ***//
    int in1_pdg = 0;
    int in2_pdg = cs.SpecPdg(in2_id);
    //asano memo fill virtual particle after 1st step
    std::vector <int> in1_pdg_cand;
    if( mode_ts==0 ){ // Sigma+/0/- -> Lambda conversion
      in1_pdg_cand.push_back(3222); // Sigma+
      in1_pdg_cand.push_back(3212); // Sigma0
      in1_pdg_cand.push_back(3112); // Sigma-
    } else if( mode_ts==1 ){ // Lambda, Sigma+/0/- elastic
      in1_pdg_cand.push_back(3122); // Lambda
      in1_pdg_cand.push_back(3222); // Sigma+
      in1_pdg_cand.push_back(3212); // Sigma0
      in1_pdg_cand.push_back(3112); // Sigma-
    } else if( mode_ts==2 ){ // p, n elastic
      in1_pdg_cand.push_back(2212); // p
      in1_pdg_cand.push_back(2112); // n
    } else if( mode_ts==3 ){ // K-p->K0n => K0ds->Lp only
      in1_pdg_cand.push_back(-311); // K0bar
    } else if( mode_ts==4 ){ // K*-n->K0Xi-
      in1_pdg_cand.push_back(-323); // K*-
    } else if( mode_ts==10  || mode_ts==11){ // K-/0 N elastic
      in1_pdg_cand.push_back(-311); // K0bar
      in1_pdg_cand.push_back(-321); // K-
    } // Add by Inoue
    else if( mode_ts==20 || mode_ts==21 || mode_ts==22 || mode_ts==23 ){
      in1_pdg_cand.push_back(-311); // K0bar
      in1_pdg_cand.push_back(-321); // K-
    }
    //Added by Asano
    else if( mode_ts==30){//K-p -> K0n(1st) -> n-n elastic (2nd) in E31
      in1_pdg_cand.push_back(2112);//n
    }
    //Added by Asano
    else if( mode_ts==40){//K-p->K0n (1st),  K0n -> pi+pi-Lambda (2nd)
      in1_pdg_cand.push_back(-311);//K0bar
    }
    else{
      std::cout<<"!!!!! unknown two step mode !!!!! "<<mode_ts<<std::endl;
      exit(0);
    }
    
    //asano memo.
    //in two-step mode, FindlPdg(i) stores the output particle of the 1st step reaction
    //in1_id stores the index of the output particle
    for( int i=0; i<nFinl; i++ ){
      int pdg = cs.FinlPdg(i);
      for( unsigned int j=0; j<in1_pdg_cand.size(); j++ ){
        if( pdg==in1_pdg_cand[j] ){
          flag_ts = true;
          in1_id  = i;
          in1_pdg = pdg;
        }
      }
    }

    //--- particle selection ---//
    //asano memo
    //set final products produced by 2nd step
    //
    std::vector <G4ParticleDefinition*> particle_ts;
    if( flag_ts ){
      if( mode_ts==0 ){
        particle_ts.push_back(particleTable->FindParticle(3122)); // final hyperon is Lambda
        // --- Sigma+ case --- // *** NO Sigma+ p
        if( in1_pdg==3222 ){
          if( in2_pdg==2112 ){ // Sigma+ n -> Lambda p
            particle_ts.push_back(particleTable->FindParticle(2212));
          }
          else if( in2_pdg==1000010020 ){ // Sigma+ d -> Lambda p p
            particle_ts.push_back(particleTable->FindParticle(2212));
            particle_ts.push_back(particleTable->FindParticle(2212));
          }
        }
        // --- Sigma0 case --- //
        else if( in1_pdg==3212 ){
          particle_ts.push_back(particleTable->FindParticle(in2_pdg)); // Sigma0 X -> Lambda X;
        }
        // --- Sigma- case --- // *** NO Sigma- n
        else if( in1_pdg==3112 ){
          if( in2_pdg==2212 ){ // Sigma- p -> Lambda n
            particle_ts.push_back(particleTable->FindParticle(2112));
          }
          else if( in2_pdg==1000010020 ){ // Sigma- d -> Lambda n n
            particle_ts.push_back(particleTable->FindParticle(2112));
            particle_ts.push_back(particleTable->FindParticle(2112));
          }
        }//Sigma case end
      }
      else if( mode_ts==1 || mode_ts==2 || mode_ts==10 ){
        particle_ts.push_back(particleTable->FindParticle(in1_pdg));
        particle_ts.push_back(particleTable->FindParticle(in2_pdg));
      }
      else if( mode_ts==3 ){ // K-p->K0n => K0ds->Lp only
        particle_ts.push_back(particleTable->FindParticle(3122)); // Lambda
        particle_ts.push_back(particleTable->FindParticle(2212)); // p
      }
      else if( mode_ts==4 ){//  K*-n->K0Xi-
        particle_ts.push_back(particleTable->FindParticle(311)); // K0
        particle_ts.push_back(particleTable->FindParticle(3312)); // Xi-
      }
      //Added by Asano
      //K-n -> K-n (elastic) ps (1st), K-p -> K0n (2nd)
      else if( mode_ts==11){
        particle_ts.push_back(particleTable->FindParticle(-311));
        particle_ts.push_back(particleTable->FindParticle(2112));
      }
      else if( mode_ts==30 ){//K-p -> K0n(1st), n-n elastic in E31  
        particle_ts.push_back(particleTable->FindParticle(2112));
        particle_ts.push_back(particleTable->FindParticle(2112));
      }
      else if( mode_ts==40){//K0n -> pi+pi-Lambda (2nd) 
        particle_ts.push_back(particleTable->FindParticle(-211));
        particle_ts.push_back(particleTable->FindParticle(211));
        particle_ts.push_back(particleTable->FindParticle(3122));
      }
    } // if( flag_ts )

    //--- reaction and fill ---//
    if( 1<particle_ts.size() ){
      G4LorentzVector in1, in2;
      in1.setVectM(vec[in1_id], mass[in1_id]);//output of 1st step, incoming to 2nd step reaction in CM frame (beam + target frame)
      in2.setVectM(vec[nFinl],  mass[nFinl]);//nFinl = # of particles in final state of 1st step, spectator neutron in E31 case in CM frame (beam + target frame)
      G4ThreeVector boost_ts = (in1+in2).boostVector();
      G4double CMmass_ts  = (in1+in2).m();//cm mass of 2nd step reaction
      
      G4int nBody_ts = particle_ts.size();//
      G4double      mass_ts[3];//
      G4ThreeVector vec_ts[3];//mom. vec. in lab frame
      while( true ){
        G4double sum_ts = 0;
        //std::cerr<<nBody_ts<<std::endl;
        for( int i=0; i<nBody_ts; i++ ){
          mass_ts[i] = GetMass(*particle_ts[i]);
          vec_ts[i]  = G4ThreeVector(0,0,0);
          sum_ts += mass_ts[i];
          //std::cerr<<particle_ts[i]->GetParticleName()<<" "<<mass_ts[i]<<std::endl;
        }
	//std::cerr<<"ene: "<<sum_ts<<" "<<CMmass_ts<<std::endl;
        if( CMmass_ts<sum_ts ){ goto START; }
        if( ManyBody(CMmass_ts, nBody_ts, mass_ts, vec_ts) ){ break; }
      }//while 

      std::vector <G4LorentzVector> finl_ts;
      for( int i=0; i<nBody_ts; i++ ){
        G4LorentzVector tmp;
        tmp.setVectM(vec_ts[i], mass_ts[i]);//CM frame of 2nd step particles
        tmp.boost(boost_ts);//CM in 2nd step -> CM in beam + target
        finl_ts.push_back(tmp);// 
      }

      vec[in1_id]      = finl_ts[0].vect();//CM in beam + target
      mass[in1_id]     = mass_ts[0];
      particle[in1_id] = particle_ts[0];
      for( int i=0; i<nBody_ts-1; i++ ){
        vec[nFinl+i]      = finl_ts[i+1].vect();//CM in beam + target
        mass[nFinl+i]     = mass_ts[i+1];
        particle[nFinl+i] = particle_ts[i+1];
      }

      nSpec--;
      nFinl = nFinl-1+nBody_ts;
      nBody = FermiFlag ? nFinl+nSpec : nFinl;

    } // if( particle_ts.size() )

  } // if( anaManager->GetTwoStep() && cs.SpecPdgSize() )

  anaManager->AddCounter(); // for number of generated events

  //---------------------------------------------------//
  // forward acceptance selection for inital particles //
  //---------------------------------------------------//
  //  for "initial" forward neutron
  //     0.0 degrees < theta < 8.0 degrees
  //     0.2 GeV/c   <  mom  < 2.0 GeV/c
  if( anaManager->GetFowardAccept() ){
    const double theta_min = 0.0*degree;
    const double theta_max = 8.0*degree;
    const double mom_min   = 0.2*GeV;
    const double mom_max   = 2.0*GeV;
    double theta, mom;
    for( int i=0; i<nBody; i++ ){
      //asano memo
      //neutron selection (pdg==2112)
      if( particle[i]->GetPDGEncoding()==2112 ){
	G4LorentzVector tmp;
	tmp.setVectM(vec[i], mass[i]);
	tmp.boost(boost);
	theta = tmp.theta();
	mom   = tmp.rho();
	//std::cerr<<theta*180/3.14<<" "<<mom/GeV<<std::endl;
	if( !(theta_min<theta && theta<theta_max &&
	      mom_min<mom && mom<mom_max) ) goto START;
	//std::cerr<<" passed: "<<theta*180/3.14<<" "<<mom/GeV<<std::endl;
      }
    }
  } // if( anaManager->GetFowardAccept() ){
  //std::cerr<<"num of loop = "<<counter<<std::endl;
  
  if(MakeUniformInqmass==1){
    //---------------------------------------------------//
    //uniform distribution for 3-body 
    //---------------------------------------------------//
    G4LorentzVector lvec[3];
    for(int i=0;i<3;i++){
      lvec[i].setVectM(vec[i], mass[i]);//vec is 3-mom. vec in CM frame.
      lvec[i].boost(boost);//boost to the lab frame
    }
    //0: missing neutron
    //1: S+/-
    //2: pi-/+
    G4LorentzVector TL_piSigma = lvec[1]+lvec[2];

    //npipiL simulation
    //0: n
    //1: pi+
    //2: pi-
    //3: missing Lambda
    //G4LorentzVector TL_piSigma = lvec[0]+lvec[1]+lvec[2];
    G4ThreeVector beammom(0,0,1000.);
    G4LorentzVector TL_beam;
    TL_beam.setVectM(beammom,493.);
    //double q = (TL_beam.vect()-lvec[3].vect()).mag()/1000.;//npipiL
    double q = (TL_beam.vect()-lvec[0].vect()).mag()/1000.;//piSigma
    double piSmass = TL_piSigma.m()/1000.;

    double prob = h2genprob->Interpolate(piSmass,q);
    if(  prob <  G4UniformRand()) goto START;
  }
  
  //npipiL sim.
  if(MakeUniformInqmass==2){
    G4LorentzVector lvec[4];
    for(int i=0;i<4;i++){
      lvec[i].setVectM(vec[i], mass[i]);//vec is 3-mom. vec in CM frame.
      lvec[i].boost(boost);//boost to the lab frame
    }
    //0: missing neutron
    //1: S+/-
    //2: pi-/+
    G4LorentzVector TL_piSigma = lvec[1]+lvec[2];

    //npipiL simulation
    //0: n
    //1: pi+
    //2: pi-
    //3: missing Lambda
    //G4LorentzVector TL_piSigma = lvec[0]+lvec[1]+lvec[2];
    G4ThreeVector beammom(0,0,1000.);
    G4LorentzVector TL_beam;
    TL_beam.setVectM(beammom,493.);
  //double q = (TL_beam.vect()-lvec[3].vect()).mag()/1000.;//npipiL
    double q = (TL_beam.vect()-lvec[0].vect()).mag()/1000.;//piSigma
    double piSmass = TL_piSigma.m()/1000.;

    double prob = h2genprob->Interpolate(piSmass,q);
    if(  prob <  G4UniformRand()) goto START;
  }

  //---------------------//
  //--- set particles ---//
  //---------------------//
  G4ThreeVector vtx = anaManager->GetVertexPosition();
  G4PrimaryVertex* vertex= new G4PrimaryVertex(vtx, 0.*ns);
  for( int i=0; i<nBody; i++ ){
    G4LorentzVector tmp;
    tmp.setVectM(vec[i], mass[i]);
    tmp.boost(boost);
    //std::cerr<<tmp.theta()<<std::endl;
    G4PrimaryParticle *primary= new G4PrimaryParticle(particle[i], tmp.px(), tmp.py(), tmp.pz(), tmp.e());
    vertex->SetPrimary(primary);
  }
  //### add the primary vertex ###//
  anEvent->AddPrimaryVertex(vertex);
  //vertex->Print();
  //anEvent->GetPrimaryVertex()->Print();


  //--------------------------//
  //--- class ReactionData ---//
  //--------------------------//
  reactionData->Init();
  reactionData->SetReactionID(cs.Id());
  for( int i=0; i<nBody; i++ ){
    reactionData->SetPDG(particle[i]->GetPDGEncoding());
    G4LorentzVector tmp;
    tmp.setVectM(vec[i], mass[i]);
    reactionData->SetCMParticle( tmp.px(), tmp.py(), tmp.pz(), tmp.m() );
    tmp.boost(boost);
    reactionData->SetParticle( tmp.px(), tmp.py(), tmp.pz(), tmp.m() );
  }
  reactionData->SetInitPDG(Init1->GetPDGEncoding());
  reactionData->SetInitPDG(Init2->GetPDGEncoding());
  reactionData->SetInitParticle( beam.px(), beam.py(), beam.pz(), beam.m() );
  reactionData->SetInitParticle( tgt.px(), tgt.py(), tgt.pz(), tgt.m() );
  reactionData->SetNParticle(nFinl, nSpec);
  reactionData->SetFermiMom(0, fermimom1);
  reactionData->SetFermiMom(1, fermimom2);
  anaManager->setHistReactionData(reactionData);

  return cs.Id();
}
    
//asano memo
//phase space generator with angular distribution
//param[] is the Legendre coefficient
//////////////////////////////////////////////////////
bool KnuclPrimaryGeneratorAction::ManyBody(G4double CMmass, G4int nBody,
					   const G4double* mass, G4ThreeVector*vec,
					   G4int nFinl,
					   G4int nparam, const G4double* param, G4double max)
//////////////////////////////////////////////////////
{
  if( !nparam ){
    if( !ManyBody(CMmass, nBody, mass, vec) ) return false;
  }
  else{
    // -----------------------------------------------------------------------
    // vec[0] must be a target particle for calculaton of the production angle
    // -----------------------------------------------------------------------
    G4LorentzVector fin[2], CM; // nFinl = 2 only
    G4int i;
    G4ThreeVector boost;
    G4double x, S;
    while( true ){
      if( !ManyBody(CMmass, nBody, mass, vec) ) return false;
      
      CM = G4LorentzVector();
      for( i=0; i<nFinl; i++ ){
	fin[i].setVectM(vec[i], mass[i]);
	CM += fin[i];
      }
      boost = CM.boostVector();
      for( i=0; i<nFinl; i++ ){
	fin[i].boost(-boost);
	//std::cout<<fin[i].z()<<" "<<fin[i].y()<<" "<<fin[i].z()<<std::endl;
      }
      x = fin[0].cosTheta();
      
      S = 0;
      for( i=0; i<nparam; i++ ) S += param[i]*Legendre(i, x);
      //std::cout<<x<<", "<<S<<std::endl;
      if ( G4UniformRand()*max < S ) break;
    }
  }
  return true;
}

//////////////////////////////////////////////////////
bool KnuclPrimaryGeneratorAction::ManyBody(G4double CMmass, G4int nBody,
					   const G4double* mass, G4ThreeVector*vec)
//////////////////////////////////////////////////////
{
  G4int total_mass = 0;
  //G4double total_mass = 0;
  for( int i=0; i<nBody; i++ ){
    total_mass += mass[i];
  }
  if( CMmass < total_mass ){
    std::cerr<<"KnuclPrimaryGeneratorAction::ManyBody, CMmass < total_mass!!!"<<std::endl;
    return false;
  }

  G4double WeightMAX = 0;
  switch( nBody ){
  case 2: WeightMAX = 1.0; break;
  case 3: WeightMAX = 0.5; break;
  case 4: WeightMAX = 0.14;  break;
  case 5: WeightMAX = 0.037;  break;
  case 6: WeightMAX = 0.01;  break;
  default: std::cerr<<"KnuclPrimaryGeneratorAction::ManyBody, n>6 body-decay is not supported!!!"<<std::endl;
    return false;
  }

  TLorentzVector W(0.0, 0.0, 0.0, CMmass);
  //static G4double weightMAX = 0;
  gen->SetDecay(W, nBody, mass);
  Double_t weight, ran;
  while( true ){
    weight = (G4double)gen->Generate();
    if ( isnan(weight) ){
      std::cerr<<"KnuclPrimaryGeneratorAction::ManyBody, return value of gen->Generate() is nan!!!"<<std::endl;
      return false;
    }
    ran = WeightMAX*G4UniformRand();
    if( ran<weight  ) {
        break;
    }
    //break; // for debug
  }
  //if ( weightMAX<weight ){ weightMAX = weight; std::cout<<weightMAX<<std::endl; }

  for( int i=0; i<nBody; i++ ){
    vec[i] = ConvVecTG( gen->GetDecay(i)->Vect() );
  }

  return true;
}

//////////////////////////////////////////////////////
G4double KnuclPrimaryGeneratorAction::GetMass(const G4ParticleDefinition& particle)
//////////////////////////////////////////////////////
{
  G4double mass  = particle.GetPDGMass();
  G4int pid = particle.GetPDGEncoding();
  if( particle.GetPDGStable() || !particle.IsShortLived() ){ // stable particles
    return mass;
  }else if( pid==9999 && anaManager->GetKppShape() ){ // Kpp
    G4double th = anaManager->GetKppMassThreshold();
    while( true ){
      G4double dmass = anaManager->GetKppMassSpectrum()->GetRandom()*1000;
      if( th<dmass ) return dmass;
    }
  }else{
    G4ParticleTable *pTable = G4ParticleTable::GetParticleTable();
    G4double th = 0;
    if( pid==9999 ) // Kpp
      th = anaManager->GetKppMassThreshold();
    else if( pid==113 ) // rho0
      th = pTable->FindParticle("pi+")->GetPDGMass() + pTable->FindParticle("pi-")->GetPDGMass();
    else if( pid==2214 || pid==2114 || pid==1114 ) // Delta
      th = pTable->FindParticle("neutron")->GetPDGMass() + pTable->FindParticle("pi+")->GetPDGMass();
    else if( pid==3224 || pid==3214 || pid==3114 ) // Sigma(1385)
      th = pTable->FindParticle("sigma-")->GetPDGMass() + pTable->FindParticle("pi+")->GetPDGMass();
    else if( pid==13122 ) // Lambda(1405)
      th = pTable->FindParticle("sigma-")->GetPDGMass() + pTable->FindParticle("pi+")->GetPDGMass();
    else if( pid==3124 ) // Lambda(1520)
      th = pTable->FindParticle("neutron")->GetPDGMass() + pTable->FindParticle("anti_kaon0")->GetPDGMass();
    else if( pid==13124 ) // Lambda(1690)
      th = pTable->FindParticle("neutron")->GetPDGMass() + pTable->FindParticle("anti_kaon0")->GetPDGMass();
    else if( abs(pid)==313 || abs(pid)==323 ) // K*
      th = pTable->FindParticle("anti_kaon0")->GetPDGMass() + pTable->FindParticle("pi-")->GetPDGMass();
    while( true ){
      G4double dmass = CLHEP::RandBreitWigner::shoot(mass, particle.GetPDGWidth());
      if( th<dmass ) return dmass;
    }
  }
}

//////////////////////////////////////////////////////
G4ThreeVector KnuclPrimaryGeneratorAction::RandUnitVec()
//////////////////////////////////////////////////////
{
  G4double costheta = 2.*G4UniformRand()-1.0;
  G4double sintheta = std::sqrt(1.0-sqr(costheta));
  G4double phi = twopi*G4UniformRand()*rad;
  G4double sinphi = std::sin(phi);
  G4double cosphi = std::cos(phi);
  G4ThreeVector direction0(sintheta*cosphi, sintheta*sinphi, costheta);
  return direction0;
}


//////////////////////////////////////////////////////
const CrossSection& KnuclPrimaryGeneratorAction::GenerateReac()
//////////////////////////////////////////////////////
{
  double ran = csTable->TotalCS()*G4UniformRand();
  int i;
  //std::cout<<"Generate(): ran / i = "<<ran<<" / ";
  for (i=0; i<csTable->CSSize(); i++){
    ran -= csTable->CS(i).Cs()*csTable->CS(i).CsFactor();
    if (ran<0) break;
  }
  //std::cout<<i<<std::endl;
  return csTable->CS(i);
}

//////////////////////////////////////////////////////
void KnuclPrimaryGeneratorAction::PrintCS(const CrossSection& cs)
//////////////////////////////////////////////////////
{
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  std::vector <TString> iname;
  for (int i=0; i<cs.InitPdgSize(); i++){
    iname.push_back(particleTable->FindParticle(cs.InitPdg(i))->GetParticleName());
  }

  std::vector <TString> fname;
  for (int i=0; i<cs.FinlPdgSize(); i++){
    fname.push_back(particleTable->FindParticle(cs.FinlPdg(i))->GetParticleName());
  }

  std::vector <TString> sname;
  for (int i=0; i<cs.SpecPdgSize(); i++){
    sname.push_back(particleTable->FindParticle(cs.SpecPdg(i))->GetParticleName());
  }

  // initial reaction
  G4ParticleDefinition *Init[2] = { particleTable->FindParticle(cs.InitPdg(0)),
				    particleTable->FindParticle(cs.InitPdg(1)) };
  G4double MassInit[2]   = { Init[0]->GetPDGMass(), Init[1]->GetPDGMass() };
  G4ThreeVector beam_mom = G4ThreeVector(0.0, 0.0, cs.InitMom());
  G4ThreeVector tgt_mom  = G4ThreeVector(0.0, 0.0, 0.0);
  G4LorentzVector beam   = G4LorentzVector(beam_mom, sqrt(beam_mom.mag2()+MassInit[0]*MassInit[0]));
  G4LorentzVector tgt    = G4LorentzVector(tgt_mom, sqrt(tgt_mom.mag2()+MassInit[1]*MassInit[1]));
  G4double CMmass        = (beam+tgt).m();

  G4double SumFin = 0;
  for (int i=0; i<cs.FinlPdgSize(); i++){
    SumFin += particleTable->FindParticle(cs.FinlPdg(i))->GetPDGMass();
  }
  G4double SumSpec = 0;
  for (int i=0; i<cs.SpecPdgSize(); i++){
    SumSpec += particleTable->FindParticle(cs.SpecPdg(i))->GetPDGMass();
  }
  G4double Sum = SumFin+SumSpec;
  G4double Qvalue = CMmass-Sum;

  // elementary reaction
  G4ParticleDefinition *eInit[2] = { particleTable->FindParticle(cs.InitPdg(0)),
				     particleTable->FindParticle(cs.InitPdg(2)) };
  G4double eMassInit[2]   = { eInit[0]->GetPDGMass(), eInit[1]->GetPDGMass() };
  G4LorentzVector ebeam   = G4LorentzVector(beam_mom, sqrt(beam_mom.mag2()+eMassInit[0]*eMassInit[0]));
  G4LorentzVector etgt    = G4LorentzVector(tgt_mom, sqrt(tgt_mom.mag2()+eMassInit[1]*eMassInit[1]));
  G4double eCMmass        = (ebeam+etgt).m();
  G4double eQvalue = eCMmass-SumFin;

  std::cout<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<std::endl;
  std::cout<<"  reaction ID:               "<<cs.Id()<<std::endl;
  std::cout<<"  CrossSection [mb]:         "<<cs.Cs()<<std::endl;
  std::cout<<"  CS factor:                 "<<cs.CsFactor()<<std::endl;
  std::cout<<"  initial PDG code:          "<<cs.InitPdg(0)<<", "<<cs.InitPdg(1)<<", "<<cs.InitPdg(2)<<std::endl;
  std::cout<<"  initial PDG name:          "<<iname[0]<<", "<<iname[1]<<", "<<iname[2]<<std::endl;
  std::cout<<"  beam momentum [MeV/c]:     "<<cs.InitMom()<<std::endl;

  std::cout<<"  # of final particles:      "<<cs.FinlPdgSize()<<std::endl;
  std::cout<<"    PDG code:                "<<cs.FinlPdg(0);
  for (int i=1; i<cs.FinlPdgSize(); i++){ std::cout<<", "<<cs.FinlPdg(i); } std::cout<<std::endl;
  std::cout<<"    PDG name:                "<<fname[0];
  for (int i=1; i<cs.FinlPdgSize(); i++){ std::cout<<", "<<fname[i]; } std::cout<<std::endl;

  std::cout<<"  # of spectators:           "<<cs.SpecPdgSize()<<std::endl;
  if( cs.SpecPdgSize() > 0 ){
    std::cout<<"    PDG code:                "<<cs.SpecPdg(0);
    for (int i=1; i<cs.SpecPdgSize(); i++){ std::cout<<", "<<cs.SpecPdg(i); } std::cout<<std::endl;
    std::cout<<"    PDG name:                "<<sname[0];
    for (int i=1; i<cs.SpecPdgSize(); i++){ std::cout<<", "<<sname[i]; } std::cout<<std::endl;
  }

  std::cout<<"  degree of Legendre:        "<<cs.PolSize()<<std::endl;
  if( cs.PolSize() > 0 ){
    std::cout<<"  coefficients of func:      "<<cs.Pol(0);
    for (int i=1; i<cs.PolSize(); i++){ std::cout<<", "<<cs.Pol(i); } std::cout<<std::endl;
    std::cout<<"  maximun random number:     "<<cs.PolMax()<<std::endl;
  }
  std::cout<<"  Q-value of K-3He           "<<Qvalue<<"MeV ("<<CMmass<<"-"<<Sum<<")"<<std::endl;
  std::cout<<"   [elementary reaction:     "<<eQvalue<<"MeV ("<<eCMmass<<"-"<<SumFin<<")]"<<std::endl;
}

//////////////////////////////////////////////////////
void KnuclPrimaryGeneratorAction::PrintAllCS()
//////////////////////////////////////////////////////
{
  // #################################
  // # dump all reactions
  // #################################
  std::cout<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<std::endl;
  for( int i=0; i<csTable->CSSize(); i++ ){
    PrintCS(csTable->CS(i));
  }
  std::cout<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<std::endl;
  //--- total-CS dump ---//
  csTable->CalcTotalCS();
  std::cout<<"======================================================"<<std::endl;
  std::cout<<" ***** input total cross-sections *****"<<std::endl;
  std::cout<<"        "<<csTable->TotalCS()<<" [mb]"<<std::endl;
  std::cout<<"======================================================"<<std::endl;
}

//////////////////////////////////////////////////////
void KnuclPrimaryGeneratorAction::WriteOutputOscarFile(G4Event* anEvent, int reacID)
//////////////////////////////////////////////////////
// format: OSC1997A
// reference: http://karman.physics.purdue.edu/OSCAR-old/docs/file/cascade_output_format/index.html
{
  ofstream* OutputFile = anaManager->GetOutputData();
  int ipart = 0;
  for( int i=0; i<anEvent->GetNumberOfPrimaryVertex(); i++ ){
    ipart += anEvent->GetPrimaryVertex(i)->GetNumberOfParticle();
  }
  *OutputFile<<anEvent->GetEventID()<<" " // event number
	     <<ipart<<" " // number of particles in an event
	     <<reacID<<" " // impact parameter --> cross-section ID
	     <<"0.0"<<std::endl; // azimuthal angle --> temporal
  int n = 0;
  for( int i=0; i<anEvent->GetNumberOfPrimaryVertex(); i++ ){
    for( int j=0; j<anEvent->GetPrimaryVertex(i)->GetNumberOfParticle(); j++ ){
      G4PrimaryParticle* tmp = anEvent->GetPrimaryVertex(i)->GetPrimary(j);
      double E = sqrt(tmp->GetMomentum().mag2()+tmp->GetMass()*tmp->GetMass());
      *OutputFile<<n<<" " // ipart
		 <<tmp->GetPDGcode()<<" " // id
		 <<tmp->GetPx()/GeV<<" " // px
		 <<tmp->GetPy()/GeV<<" " // py
		 <<tmp->GetPz()/GeV<<" " // pz
		 <<E/GeV<<" " // p0
		 <<tmp->GetMass()/GeV<<" " // mass
		 <<anEvent->GetPrimaryVertex(i)->GetX0()<<" " // x
		 <<anEvent->GetPrimaryVertex(i)->GetY0()<<" " // y
		 <<anEvent->GetPrimaryVertex(i)->GetZ0()<<" " // z
		 <<anEvent->GetPrimaryVertex(i)->GetT0()<<std::endl; // t
      n++;
    }
  }
}
