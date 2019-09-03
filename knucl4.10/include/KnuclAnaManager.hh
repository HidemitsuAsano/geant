#ifndef KnuclAnaManager_h
#define KnuclAnaManager_h 1

#include "globals.hh"
#include "KnuclHist.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4LorentzVector.hh"

#include <TROOT.h>
#include <TFile.h>
#include <TH1F.h> 

#include <fstream>
#include <cstdio>
#include "ComCardsInput.hh"
#include "ComCrossSectionTable.hh"

class KnuclAnaManager
{
public:
  typedef enum { Single, Beam, K3He, ReadFromFile, neutron_efficiency_mode} Process_t;

//### for singleton ###//
private:
  // this object is a singleton.
  static KnuclAnaManager* ptrAnaManager;
protected:
  KnuclAnaManager(std::string filename);  // constructor should NOT be public.
public:
  static KnuclAnaManager* GetPointer(std::string filename)
  {
    if(!ptrAnaManager) ptrAnaManager= new KnuclAnaManager(filename);
    else{
      std::cout<<"!!! KnuclAnaManager initialized at once !!!"<<std::endl;
      exit(0);
    }
    return ptrAnaManager;
  }
//### for singleton ###//
  
private:
  KnuclHist hist;

public:

  //KnuclAnaManager();
  ~KnuclAnaManager();
  void BeginOfRunAction();
  void EndOfRunAction();
  void BeginOfEventAction();
  void EndOfEventAction();

  void InitParameters();
  void SetParameters(ComCardsInput *inp);
  void PrintParameters();

  void ReadOscarFile();
  void InitOutputOscarFile();

  G4int GetReactionID(){ return hist.GetReactionID();};
  void SetReactionID(G4int ID){ hist.SetReactionID(ID);};

  //----------//
  // for tree
  //----------//
  void setHistRunHeaderMC(RunHeaderMC* val){hist.setRunHeaderMC(val);}
  void setHistEventHeaderMC(EventHeaderMC* val){hist.setEventHeaderMC(val);}
  void setHistDetectorData(DetectorData* val){hist.setDetectorData(val);}
  void setHistMCData(MCData* val){hist.setMCData(val);}
  void setHistReactionData(ReactionData* val){hist.setReactionData(val);}
  void fillHistTree(){hist.fillTree();}
  void showHistTree(){hist.showTree();}
  void fillHistTree2(){hist.fillTree2();}
  void showHistTree2(){hist.showTree2();}

//===========//
// for Knucl //
//===========//
private:
  TString       PhysicsList;
  Process_t     ProcessID;
  TString       CSFileName;
  TString       CardFileName;
  G4int         FermiMotion;
  G4int         FermiMotionMode;
  G4int         TwoStep;  
  G4int         TwoStepMode;  
  G4double      BeamMomentum;
  G4double      BindingEnergy;
  G4double      DecayWidth;
  G4int         DecayMode;
  G4int         KppShape;
  TString       KppShapeFile;
  TString       KppShapeHist;
  G4int         FowardAccept;
public:
  const char*   GetPhysicsList()             {return PhysicsList;     }
  Process_t     GetProcessID()               {return ProcessID;       }
  void          SetProcessID(Process_t id=(Process_t)0) {ProcessID = id;         }
  const char*   GetCSFileName()              {return CSFileName; }
  G4int         GetFermiMotion()             {return FermiMotion;    }
  G4int         GetFermiMotionMode()         {return FermiMotionMode;    }
  G4int         GetTwoStep()                 {return TwoStep;    }
  G4int         GetTwoStepMode()             {return TwoStepMode;    }
  void          SetBeamMomentum(G4double mom){BeamMomentum=mom;}
  G4double      GetBeamMomentum()            {return BeamMomentum;   }
  G4double      GetBindingEnergy()           {return BindingEnergy;   }
  G4double      GetDecayWidth()              {return DecayWidth;      }
  G4int         GetDecayMode()               {return DecayMode;      }
  G4int         GetKppShape()                {return KppShape;      }
  const char*   GetKppShapeFile()            {return KppShapeFile;      }
  const char*   GetKppShapeHist()            {return KppShapeHist;      }
  G4int         GetFowardAccept()            {return FowardAccept;      }
private:
  TString       ProcName;


//============// 
// for Single //
//============// 
private:
  TString       Particle;
  G4double      Momentum_MIN;
  G4double      Momentum_MAX;
  G4double      Theta_MIN;
  G4double      Theta_MAX;
  G4double      Vertex_MIN_X;
  G4double      Vertex_MIN_Y;
  G4double      Vertex_MIN_Z;
  G4double      Vertex_MAX_X;
  G4double      Vertex_MAX_Y;
  G4double      Vertex_MAX_Z;
public:
  TString       GetSingleParticle()    {return Particle;}
  G4double      GetSingleMomentumMIN() {return Momentum_MIN;}
  G4double      GetSingleMomentumMAX() {return Momentum_MAX;}
  G4double      GetSingleThetaMIN()    {return Theta_MIN;}
  G4double      GetSingleThetaMAX()    {return Theta_MAX;}
  G4ThreeVector GetSingleVertexMIN()   {return G4ThreeVector(Vertex_MIN_X, Vertex_MIN_Y, Vertex_MIN_Z);}
  G4ThreeVector GetSingleVertexMAX()   {return G4ThreeVector(Vertex_MAX_X, Vertex_MAX_Y, Vertex_MAX_Z);}

//==============// 
// for Detector //
//==============// 
private:
  TString       ConfFileName;
  G4int         DORAFieldMapFlag;
  G4int         USWKFieldMapFlag;
  G4double      FieldInCDC;
  G4double      FieldInUshiwaka;
  G4double      ADCThreshold;
  std::map<int, double> ADCThresholdMap;
  G4int         FWDTrigger;
  G4String      TriggerName;
public:
  const char*   GetConfFileName()    { return ConfFileName; }
  G4int         GetDORAFieldMapFlag(){ return DORAFieldMapFlag;     }
  G4int         GetUSWKFieldMapFlag(){ return USWKFieldMapFlag;     }
  G4double      GetFieldInCDC()      { return FieldInCDC;     }
  G4double      GetFieldInUshiwaka() { return FieldInUshiwaka;  }
  G4double      GetADCThreshold()    { return ADCThreshold;     }
  G4double      GetADCThreshold(const int &cid){ return ADCThresholdMap.find(cid)==ADCThresholdMap.end() ? GetADCThreshold() : ADCThresholdMap.at(cid); }
  G4int         GetFWDTrigger()      { return FWDTrigger;     }
  G4String      GetTrigger()         { return TriggerName; }

//================// 
// for TargetInfo //
//================// 
private:
  G4double      TargetLength;
  G4int         TargetMaterialID;
public:
  G4double      GetTargetLength()        {return TargetLength;       }
  G4int         GetTargetMaterialID()    {return TargetMaterialID;   }
private:
  TString       TargetName;
//==============// 
// for FileInfo //
//==============//
private:
  TString       InputDataFile;
  TString       OutputRootFile;
  G4int         FillLevel;
  TString       AddCardFile;
  TString       OutputDataFile;
public:
  void          SetOutputRootFile(const std::string &name){ OutputRootFile=name; }
  void          SetCSFile(const std::string &name){ CSFileName=name; }
  void          SetKppShapeFile(const std::string &name){ KppShapeFile=name; }
  TString       GetInputDataFile()    {return InputDataFile;        }
  const char*   GetInputDataFileName(){return InputDataFile.Data(); }
  ifstream*     GetInputData()        {return InputFile;}
  TString       GetOutputRootFile()   {return OutputRootFile;     }
  G4int         GetFillLevel()        {return FillLevel;   }
  TString       GetAddCardFile()      {return AddCardFile;        }
  TString       GetOutputDataFile()   {return OutputDataFile;        }
  ofstream*     GetOutputData()       {return OutputFile;}
private:
  ifstream      *InputFile;
  ofstream      *OutputFile;

//==================// 
// for AnalysisInfo //
//==================//
private:
  TString       Seed;
public:
  TString       GetSeed()   {return Seed;     }
 

//==============// 
// for DataInfo //
//==============// 
private:
  TString       Comment;
public:
  const char*   GetComment() { return Comment.Data(); }


//===============// 
// for DebugInfo //
//===============// 
private:
  TString       DebugFlag;
public:
  const char*   GetDebugFlag() { return DebugFlag.Data(); }


//Uniform Gen
private:
  G4int      UniformGenFlag;
public:
  G4int GetUniformGenFlag() { return UniformGenFlag;}

private:
  G4ThreeVector BeamDirection;
  G4ThreeVector VertexPosition;
  G4ThreeVector BeamMom;
  G4double      BeamEne;
  G4double      BeamMass;
  G4double      KppMassThreshold;

public:
  void SetBeamDirection (G4double fx, G4double fy, G4double fz){BeamDirection.set (fx,fy,fz);}
  void SetVertexPosition(G4double fx, G4double fy, G4double fz){VertexPosition.set(fx,fy,fz);}
  void SetBeamMomVector (G4double fx, G4double fy, G4double fz){BeamMom.set       (fx,fy,fz);}
  void SetBeamEne  (G4double fe){BeamEne = fe;}
  void SetBeamMass (G4double fm){BeamMass = fm;}
  void SetKppMassThreshold(G4double mass){ KppMassThreshold = mass; }

  G4ThreeVector GetBeamDirection() {return BeamDirection;}
  G4ThreeVector GetVertexPosition(){return VertexPosition;} 
  G4ThreeVector GetBeamMom()       {return BeamMom;} 
  G4double      GetBeamEne()       {return BeamEne;} 
  G4double      GetBeamMass()      {return BeamMass;} 
  G4double      GetKppMassThreshold(){ return KppMassThreshold; }

private:
  typedef std::vector<G4LorentzVector> LorentzVectorContainer;
  LorentzVectorContainer cmOutParticleContainer;
  LorentzVectorContainer OutParticleContainer;
  
public:
  void SetCMParticle(G4double px, G4double py, G4double pz, G4double m);
  void SetParticle(G4double px, G4double py, G4double pz, G4double m);


  // for seed number
private:
  G4int SeedNum;
public:
  void  SetSeedNum(G4int v){ SeedNum = v; }
  G4int GetSeedNum()       { return SeedNum; }


  // number of events
  //  used for acceptance study, etc.
private:
  G4int numEvent;
  G4int counter;
public:
  void  InitCounter(){counter = 0;}
  void  AddCounter(){counter++;}
  G4int GetCounter(){return counter;}

  // used for check of kinematical dist., etc.
private:
  G4double tmpValue;
public:
  void     SetTmpValue(G4double v){tmpValue = v;}
  G4double GetTmpValue(){return tmpValue;}

  // RunHeaderMC
private:
  RunHeaderMC* runHeaderMC;
public:
  void SetCSTable(const CrossSectionTable& val){runHeaderMC->setCStable(val);}


// KppShapeHist
private:
  TFile *shapeFile;
  TH1F* kppMassSpectrum;
public:
  const TH1F* GetKppMassSpectrum(){return kppMassSpectrum;}

};
#endif
