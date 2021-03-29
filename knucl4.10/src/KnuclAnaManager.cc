// Asano memo
// This code hundle analysis parameters from KnuclSetting.card

#include "G4Version.hh"
#include "G4ParticleTable.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"

#include "KnuclAnaManager.hh"
#include "globals.hh"
#include "KnuclCommon.h"
#include "Randomize.hh"

#include "TString.h"
#include "TMacro.h"

#include "DetectorList.hh"

//### for singleton ###//
KnuclAnaManager* KnuclAnaManager::ptrAnaManager= 0;
//### for singleton ###//

//////////////////////////////////////////////////////
KnuclAnaManager::KnuclAnaManager(std::string filename)
//////////////////////////////////////////////////////
{
  G4cout << "!!!! Now began !!!! " << G4endl;

  //--- read card-file ---//
  InitParameters();
  CardFileName = filename;
  ComCardsInput *inp = new ComCardsInput(filename.c_str());
  SetParameters(inp);
  if( strcmp(AddCardFile.Data(),"") ){
    ComCardsInput *inp2 = new ComCardsInput(AddCardFile.Data());
    SetParameters(inp2);
  }
  PrintParameters();

  //--- read conf-file ---//
  const Int_t MAXCHARS = 512;
  char str[MAXCHARS],str1[MAXCHARS];
  std::string DetectorListFileName="DEFAULT";

  ifstream fin(ConfFileName.Data());
  if (!fin.good()){
    std::cerr<<"conf file cant open file: "<<ConfFileName<<std::endl;
    exit(0);
  }
  while (!fin.eof()){
    fin.getline(str, MAXCHARS);
    if( str[0]=='#' ) continue;
    if( sscanf(str,"DetectorList: %s", str1)==1 )
      DetectorListFileName= str1;
  }
  DetectorList *dlist=DetectorList::GetInstance();
  dlist->Initialize(DetectorListFileName.data());

  //--- read OSCAR-file ---//
  if(ProcessID == ReadFromFile) {
    ReadOscarFile();
  }

  //--- output data-file ---//
  if( OutputDataFile ){
    InitOutputOscarFile();
  }

  //-------------------------//
  //--- class RunHeaderMC ---//
  //-------------------------//
  // constructor
  runHeaderMC = new RunHeaderMC();
  // counter list
  std::map<std::string, unsigned int> CIDContainer = dlist->GetList();
  std::map<std::string, unsigned int>::iterator it = CIDContainer.begin();
  while( it != CIDContainer.end() ){
    runHeaderMC->setCounterList((*it).second);
    //    std::cout << (*it).first << ":" << (*it).second << std::endl;
    ++it;
  }
}

//////////////////////////////////////////////////////
KnuclAnaManager::~KnuclAnaManager()
//////////////////////////////////////////////////////
{
  cmOutParticleContainer.clear();
  OutParticleContainer.clear();
  
  //--- read OSCAR-file ---//
  if(ProcessID == ReadFromFile) {
    InputFile->close();
    delete InputFile;
  }

  //--- output data-file ---//
  if( !OutputDataFile ){
    OutputFile->close();
    delete OutputFile;
  }

  //--- read KppShapeFile ---//
  if( KppShape ) {
    shapeFile->Close();
    delete shapeFile;
  }
}


//////////////////////////////////////////////////////
void KnuclAnaManager::BeginOfRunAction()
//////////////////////////////////////////////////////
{
  //--- read KppShapeFile ---//
  if( KppShape ) {
    shapeFile = new TFile(KppShapeFile.Data());
    if(!shapeFile->IsOpen()){
      std::cerr<<" cant open KppShapeFile: "<<KppShapeFile<<std::endl;
      exit(0);
    }
    std::cout<<"KppShapeFile : "<<KppShapeFile<<" is opend."<<std::endl;
    std::cout<<" shape : "<<KppShapeHist.Data()<<std::endl;
    kppMassSpectrum = (TH1F*)shapeFile->Get( Form(KppShapeHist.Data()) );
  }

  hist.BeginOfRunAction(OutputRootFile.Data());
  InitCounter();

  //-- save Card-File as text ---//
  //TMacro cardfile("KnuclSetting.card");
  TMacro cardfile(CardFileName.Data());
  cardfile.Write();
  if( strcmp(AddCardFile.Data(),"") ){
    TMacro cardfile2(AddCardFile.Data());
    cardfile2.Write();
  }
  //-- save CS-File as text ---//
  if( ProcessID==K3He ){
    TMacro csfile(CSFileName);
    csfile.Write();
  }
}


//////////////////////////////////////////////////////
void KnuclAnaManager::EndOfRunAction()
//////////////////////////////////////////////////////
{
  //-------------------------//
  //--- class RunHeaderMC ---//
  //-------------------------//
  // seed number
  runHeaderMC->setSeed(SeedNum);
  // number of total events
  const G4Run* aRun 
    = G4RunManager::GetRunManager()->GetCurrentRun();
  numEvent = aRun->GetNumberOfEvent();
  runHeaderMC->setNumEvent(numEvent);
  // number of generated events
  runHeaderMC->setNumGenerated(counter);
  std::cout<<"========================================="<<std::endl;
  std::cout<<"total number of events    : "<<numEvent<<std::endl;
  std::cout<<"total number of generated : "<<counter<<std::endl;
  std::cout<<"========================================="<<std::endl;
  setHistRunHeaderMC(runHeaderMC);
  fillHistTree2();
  delete runHeaderMC;

  hist.EndOfRunAction();
  //  G4cerr<<" AnaManager::counter = "<<GetCounter()<<G4endl;
}


//////////////////////////////////////////////////////
void KnuclAnaManager::BeginOfEventAction()
//////////////////////////////////////////////////////
{
  hist.BeginOfEventAction();
}


//////////////////////////////////////////////////////
void KnuclAnaManager::EndOfEventAction()
//////////////////////////////////////////////////////
{
}


//////////////////////////////////////////////////////
void KnuclAnaManager::InitParameters()
//////////////////////////////////////////////////////
{
  // -------------------------//
  // parameter initialization //
  // -------------------------//
  // for Knucl
  PhysicsList        = "QGSP_BERT";
  ProcessID          = K3He;
  CSFileName         = "CS.list";
  FermiMotion        = 1;
  FermiMotionMode    = 1;
  TwoStep            = 0;
  TwoStepMode        = 0;
  BeamMomentum       = 1.0;
  BindingEnergy      = 50.0;
  DecayWidth         = 50.0;
  DecayMode          = 0;
  KppShape           = 0;
  KppShapeFile       = "KppShape.root";
  KppShapeHist       = "his_shape";
  FowardAccept       = 0;
  // for Single
  Particle           = "pi-";
  Momentum_MIN       = 0.0;
  Momentum_MAX       = 1.0;
  Theta_MIN          = -1.0;
  Theta_MAX          = 1.0;
  Vertex_MIN_X       = -10.0;
  Vertex_MIN_Y       = -10.0;
  Vertex_MIN_Z       = -30.0;
  Vertex_MAX_X       = 10.0;
  Vertex_MAX_Y       = 10.0;
  Vertex_MAX_Z       = 30.0;
  // for Detector
  ConfFileName       = "conf/Run49c/analyzer.conf";
  DORAFieldMapFlag   = 1;
  USWKFieldMapFlag   = 1;
  FieldInCDC         = 0.716;
  FieldInUshiwaka    = 1.2;
  ADCThreshold       = -1;
  FWDTrigger         = 0;
  TriggerName        = "";
  // for TargetInfo
  TargetLength       = 120.0;
  TargetMaterialID   = 0;
  // for FileInfo
  InputDataFile      = "osc.dat";
  OutputRootFile     = "test.root";
  FillLevel          = 4;
  AddCardFile        = "";
  OutputDataFile     = "";
  // for AnalysisInfo
  Seed               = "random";
  // for DataInfo
  Comment            = "test";
  // for DebugInfo 
  DebugFlag          = "false";
  UniformGenFlag     = 0;
}


//////////////////////////////////////////////////////
void KnuclAnaManager::SetParameters(ComCardsInput *inp)
//////////////////////////////////////////////////////
{
  const char key0[]="Knucl";
  if (inp->FindEntry(key0,"PhysicsList",   "s",1,1)>0) PhysicsList   = inp->GetArg(0);
  if (inp->FindEntry(key0,"ProcessID",     "s",1,1)>0) ProcName      = inp->GetArg(0);
  if (inp->FindEntry(key0,"CSFileName" ,   "s",1,1)>0) CSFileName    = inp->GetArg(0);
  if (inp->FindEntry(key0,"FermiMotion"  , "d",1,1)>0) FermiMotion   = inp->GetArgD(0);
  if (inp->FindEntry(key0,"FermiMotionMode", "d",1,1)>0) FermiMotionMode = inp->GetArgD(0);
  if (inp->FindEntry(key0,"TwoStep"  ,     "d",1,1)>0) TwoStep       = inp->GetArgD(0);
  if (inp->FindEntry(key0,"TwoStepMode"  , "d",1,1)>0) TwoStepMode   = inp->GetArgD(0);
  if (inp->FindEntry(key0,"BeamMomentum",  "f",1,1)>0) BeamMomentum  = inp->GetArgF(0);
  if (inp->FindEntry(key0,"BindingEnergy", "f",1,1)>0) BindingEnergy = inp->GetArgF(0);
  if (inp->FindEntry(key0,"DecayWidth"   , "f",1,1)>0) DecayWidth    = inp->GetArgF(0);
  if (inp->FindEntry(key0,"DecayMode"    , "d",1,1)>0) DecayMode     = inp->GetArgD(0);
  if (inp->FindEntry(key0,"KppShape"    ,  "d",1,1)>0) KppShape      = inp->GetArgD(0);
  if (inp->FindEntry(key0,"KppShapeFile" , "s",1,1)>0) KppShapeFile  = inp->GetArg(0);
  if (inp->FindEntry(key0,"KppShapeHist" , "s",1,1)>0) KppShapeHist  = inp->GetArg(0);
  if (inp->FindEntry(key0,"FowardAccept" , "d",1,1)>0) FowardAccept  = inp->GetArgD(0);

  const char key1[]="Single";
  if (inp->FindEntry(key1,"Particle",     "s",1,1)>0) Particle     = inp->GetArg(0);
  if (inp->FindEntry(key1,"Momentum_MIN", "f",1,1)>0) Momentum_MIN = inp->GetArgF(0);
  if (inp->FindEntry(key1,"Momentum_MAX", "f",1,1)>0) Momentum_MAX = inp->GetArgF(0);
  if (inp->FindEntry(key1,"Theta_MIN",    "f",1,1)>0) Theta_MIN    = inp->GetArgF(0);
  if (inp->FindEntry(key1,"Theta_MAX",    "f",1,1)>0) Theta_MAX    = inp->GetArgF(0);
  if (inp->FindEntry(key1,"Vertex_MIN",   "f",1,1)>0){
    Vertex_MIN_X = inp->GetArgF(0); Vertex_MIN_Y = inp->GetArgF(1); Vertex_MIN_Z = inp->GetArgF(2);
  }
  if (inp->FindEntry(key1,"Vertex_MAX",   "f",1,1)>0){
    Vertex_MAX_X = inp->GetArgF(0); Vertex_MAX_Y = inp->GetArgF(1); Vertex_MAX_Z = inp->GetArgF(2);
  }

  const char key2[]="Detector";
  if (inp->FindEntry(key2,"ConfFileName" ,   "s",1,1)>0)  ConfFileName     = inp->GetArg(0);
  if (inp->FindEntry(key2,"DORAFieldMapFlag",   "f",1,1)>0)  DORAFieldMapFlag    = inp->GetArgF(0);
  if (inp->FindEntry(key2,"USWKFieldMapFlag",   "f",1,1)>0)  USWKFieldMapFlag    = inp->GetArgF(0);
  if (inp->FindEntry(key2,"FieldInCDC",   "f",1,1)>0)  FieldInCDC    = inp->GetArgF(0);
  if (inp->FindEntry(key2,"FieldInUshiwaka","f",1,1)>0)  FieldInUshiwaka = inp->GetArgF(0);
  if (inp->FindEntry(key2,"ADCThreshold","f",1,1)>0) ADCThreshold= inp->GetArgF(0);

  for( int i=0; i<500; i++ ){
    char det_key[512];
    sprintf(det_key, "ADCThreshold_%d", i);
    if( inp->FindEntry(key2, det_key, "f", 1, 1)>0 ) ADCThresholdMap[i]=inp->GetArgF(0);
  }

  if (inp->FindEntry(key2,"FWDTrigger",   "f",1,1)>0)  FWDTrigger    = inp->GetArgF(0);
  if (inp->FindEntry(key2,"Trigger",   "s",1,1)>0)  TriggerName    = inp->GetArg(0);
 
  const char key3[]="TargetInfo";
  if (inp->FindEntry(key3,"TargetLength",  "f",1,1)>0) TargetLength = inp->GetArgF(0);
  if (inp->FindEntry(key3,"TargetMaterial","s",1,1)>0) TargetName   = inp->GetArg(0);

  const char key4[]="FileInfo";
  if (inp->FindEntry(key4,"InputDataFile",   "s",1,1)>0) InputDataFile  = inp->GetArg(0);
  if (inp->FindEntry(key4,"OutputRootFile",  "s",1,1)>0) OutputRootFile = inp->GetArg(0);
  if (inp->FindEntry(key4,"FillLevel",       "f",1,1)>0) FillLevel      = inp->GetArgF(0);
  if (inp->FindEntry(key4,"AddCardFile",     "s",1,1)>0) AddCardFile    = inp->GetArg(0);
  if (inp->FindEntry(key4,"OutputDataFile",  "s",1,1)>0) OutputDataFile = inp->GetArg(0);

  const char key5[]="AnalysisInfo";
  if (inp->FindEntry(key5,"Seed",  "s",1,1)>0) Seed = inp->GetArg(0);

  const char key6[]="DataInfo";
  if (inp->FindEntry(key6,"Comment",  "s",1,1)>0) Comment = inp->GetArg(0);

  const char key7[]="DebugInfo";
  if (inp->FindEntry(key7,"DebugFlag",  "s",1,1)>0) DebugFlag = inp->GetArg(0);

  const char key8[]="UniformGen";
  if (inp->FindEntry(key8,"UniformGen", "f",1,1)>0) UniformGenFlag = inp->GetArgF(0);

  if ( !strcmp(TargetName.Data(),"L3He") ) {
    TargetMaterialID = 0; 
    G4cout << " Target Material is specified as L 3_He " << G4endl;  
  } else
  if ( !strcmp(TargetName.Data(),"CH2") ) {
    TargetMaterialID = 1;
    G4cout << " Target Material is specified as CH2    " << G4endl;  
  } else
  if ( !strcmp(TargetName.Data(),"C") ) {
    TargetMaterialID = 2;
    G4cout << " Target Material is specified as C      " << G4endl;  
  } else
  if ( !strcmp(TargetName.Data(),"Vaccum") ) {
    TargetMaterialID = 999;
    G4cout << " Target Material is specified as Vaccum " << G4endl;  
  } else {
    G4cout << " Target Material unkown " <<TargetName<< G4endl;
  } 

  if ( !strcmp(ProcName.Data(),"ReadFromFile") ) ProcessID = ReadFromFile;
  else 
  if ( !strcmp(ProcName.Data(),"Beam")         ) ProcessID = Beam;
  else
  if ( !strcmp(ProcName.Data(),"K3He")         ) ProcessID = K3He;
  else
  if ( !strcmp(ProcName.Data(),"Single")       ) ProcessID = Single;
  else 
  if ( !strcmp(ProcName.Data(),"neutron_efficiency_mode") ) ProcessID = neutron_efficiency_mode;
  else {
    G4cout << "Event type not specified or not supported " << G4endl;
    exit(0);
  }

}


//////////////////////////////////////////////////////  
void KnuclAnaManager::PrintParameters()
//////////////////////////////////////////////////////
{
  G4cout << " ####################################################### " << G4endl;
  G4cout << " PhysicsList                   = " << PhysicsList << G4endl;
  G4cout << " ProcessID                     = " << ProcName.Data() << G4endl;
  G4cout << " CSFileName                    = " << CSFileName  << " " << G4endl;
  G4cout << " FermiMotion                   = " << FermiMotion << G4endl;
  G4cout << " FermiMotionMode               = " << FermiMotionMode << G4endl;
  G4cout << " TwoStep                       = " << TwoStep << G4endl;
  G4cout << " TwoStepMode                   = " << TwoStepMode << G4endl;
  G4cout << " Beam Momentum                 = " << BeamMomentum   << " GeV/c" << G4endl;
  G4cout << " BindingEnergy                 = " << BindingEnergy << " MeV" << G4endl;
  G4cout << " DecayWidth                    = " << DecayWidth << " MeV" << G4endl;
  G4cout << " DecayMode                     = " << DecayMode << G4endl;
  G4cout << " KppShape                      = " << KppShape << G4endl;
  G4cout << " KppShapeFile                  = " << KppShapeFile << G4endl;
  G4cout << " KppShapeHist                  = " << KppShapeHist << G4endl;
  G4cout << " FowardAccept                  = " << FowardAccept << G4endl;
  G4cout << " ------------------------------------------------------- " << G4endl;
  G4cout << " Single Particle               = " << Particle.Data() << G4endl;  
  G4cout << " Min. Momentum                 = " << Momentum_MIN << " GeV/c" << G4endl;
  G4cout << " Max. Momentum                 = " << Momentum_MAX << " GeV/c" << G4endl;
  G4cout << " Min. Theta Angle              = " << Theta_MIN << " degree" << G4endl;
  G4cout << " Max. Theta Angle              = " << Theta_MAX << " degree" << G4endl;
  G4cout << " Min. Vertex Pos.              = (" << Vertex_MIN_X << ", "
	 << Vertex_MIN_Y << ", "<< Vertex_MIN_Z << ") mm" << G4endl;
  G4cout << " Max. Vertex Pos.              = (" << Vertex_MAX_X << ", "
	 << Vertex_MAX_Y << ", "<< Vertex_MAX_Z << ") mm" << G4endl;
  G4cout << " ------------------------------------------------------- " << G4endl;
  G4cout << " ConfFileName                  = " << ConfFileName  << " " << G4endl;
  G4cout << " DORAFieldMapFlag              = " << DORAFieldMapFlag << " " << G4endl;
  G4cout << " USWKFieldMapFlag              = " << USWKFieldMapFlag << " " << G4endl;
  G4cout << " FieldInCDC                    = " << FieldInCDC    << " T" << G4endl;
  G4cout << " FieldInUshiwaka               = " << FieldInUshiwaka << " T" << G4endl;
  G4cout << " ADCThreshold                  = " << ADCThreshold <<" MeV" << G4endl;
  G4cout << " FWDTrigger                    = " << FWDTrigger <<" " << G4endl;
  G4cout << " ------------------------------------------------------- " << G4endl;
  G4cout << " TargetLength                  = " << TargetLength <<" mm" << G4endl;
  G4cout << " TargetMaterial                = " << TargetName.Data() << G4endl;
  G4cout << " ------------------------------------------------------- " << G4endl;
  G4cout << " input data-file               = " << InputDataFile  << G4endl;
  G4cout << " output root-file              = " << OutputRootFile << G4endl;
  G4cout << " fill level                    = " << FillLevel  << G4endl;
  G4cout << " additional card-file          = " << AddCardFile  << G4endl;
  G4cout << " output data-file              = " << OutputDataFile  << G4endl;
  G4cout << " ------------------------------------------------------- " << G4endl;
  G4cout << " random seed                   = " << Seed << G4endl;
  G4cout << " ------------------------------------------------------- " << G4endl;
  G4cout << " comment                       = " << Comment << G4endl;
  G4cout << " debug flag                    = " << DebugFlag << G4endl;
  G4cout << " ####################################################### " << G4endl;
  G4cout << " " << G4endl;
  G4cout << " " << G4endl;
  G4cout << " " << G4endl;
}


//////////////////////////////////////////////////////
void KnuclAnaManager::ReadOscarFile()
//////////////////////////////////////////////////////
{
  G4cout << " ####################################################### " << G4endl;
  G4cout << " Simulation data will be read from data on disk. "         << G4endl;
  G4cout << " Input  data-file              = " << InputDataFile        << G4endl;
  G4cout << " input file format must be followed OSCAR format "         << G4endl;
  G4cout << " please check " << G4endl;
  G4cout << " https://karman.physics.purdue.edu/OSCAR-old/docs/file/cascade_output_format/node8.html" << G4endl;
  G4cout << " ####################################################### " << G4endl;
  
  InputFile=new ifstream(InputDataFile.Data());
  G4cout <<  InputFile << G4endl;
  if ( InputFile->fail() ) {
    std::cerr << InputDataFile << " doesn't exist" << std::endl;
    exit(-1);
  } else {
    G4cout << " Oscar Input File : " << InputDataFile << " is opened." << G4endl;
  }
  G4cout << " " << G4endl;
  G4cout << " " << G4endl;
  G4cout << " " << G4endl;
  
  char inln[160];
  char c;
  
  std::cout << "#####################################################" << std::endl;
  std::cout << "        Header Information of input data file        " << std::endl;
  std::cout << "#####################################################" << std::endl;
  // ReadOscarHeader
  InputFile->getline(inln,160);
  std::cout << inln << std::endl;
  InputFile->getline(inln,160);
  std::cout << inln << std::endl;
  //
  // read over comment lines beginning with #
  InputFile->get(c);
  while ( c == '#' ) {
    InputFile->getline(inln,160);
    std::cout << inln << std::endl;
    InputFile->get(c);
  }
  InputFile->putback(c);
  // here's the third line with the event information
  TString CodeName;
  TString Version;
  char *cn=new char[80];
  char *ver=new char[80];
  char *ref=new char[80];
  Int_t Aproj,Zproj,Atarg,Ztarg;
  TString RefFrame;
  Float_t Ebeam;
  Int_t NTestPart;
  InputFile->getline(inln,160);
  sscanf(inln,"%s %s (%d,%d)+(%d,%d) %s %f %d",cn,ver,&Aproj,&Zproj,
	 &Atarg,&Ztarg,ref,&Ebeam,&NTestPart);
  CodeName = cn;
  Version  = ver;
  RefFrame = ref;
  
  InputFile->get(c);
  while ( c == '#' ) {
    InputFile->getline(inln,160);
    std::cout << inln << std::endl;
    InputFile->get(c);
  }
  std::cout << "Code:"  << CodeName << " " << Version;
  std::cout << "  ("  << Aproj << "," << Zproj << ") + " ;
  std::cout << "("  << Atarg << "," << Ztarg << ")" ;
  std::cout << " at " << Ebeam << " GeV in " << RefFrame << " frame ";
  std::cout << " #tp:" << NTestPart << std::endl;
  std::cout << "#####################################################" << std::endl;
  
  InputFile->putback(c);
}

//////////////////////////////////////////////////////
void KnuclAnaManager::InitOutputOscarFile()
//////////////////////////////////////////////////////
{
  OutputFile = new ofstream(OutputDataFile.Data());
  G4cout << " Oscar Output File : " << OutputDataFile << " is opened." << G4endl;
  *OutputFile << "OSC1997A" << std::endl;
  *OutputFile << "final_id_p_x" << std::endl;
  std::string vs = G4Version;
  vs = vs.substr(1,vs.size()-2).erase(0,6);
  *OutputFile << vs;
  *OutputFile << " 1.0  ( -1,     1)+(  1,     1)  tar   0.0000E+00"; // temporal
  *OutputFile << std::endl;
}

