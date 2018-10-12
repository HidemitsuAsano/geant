

// ====================================================================
//    KnuclDetectorConstruction.hh
//
//
// ====================================================================

#ifndef KnuclDetectorConstruction_H
#define KnuclDetectorConstruction_H 1

class G4Box;
class G4Para;
class G4Tubs;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4UniformMagField;

#include "G4VUserDetectorConstruction.hh"

#include "G4Transform3D.hh"
#include "G4CSGSolid.hh"

#include "KnuclCounterSD.hh"
#include "KnuclChamberSD.hh"
#include "KnuclAnaManager.hh"
#include "KnuclMaterialManager.hh"

#include "TString.h"
#include "G4Color.hh"

class KnuclMaterialManager;
class KnuclCounterSD;
class KnuclFieldSetup;
class KnuclFieldSetup_uswk;
class KnuclFieldSetup_dora;
class KnuclFieldSetup_combi;

class KnuclDetectorConstruction : public G4VUserDetectorConstruction
{
public:

  KnuclDetectorConstruction(KnuclAnaManager* ana);
  ~KnuclDetectorConstruction();

  G4VPhysicalVolume* Construct();

private:

  KnuclMaterialManager* materialMgr;
  G4double ADCthreshold;

public:
  
  G4Material*         FindChamberGas(const TString& chamberType);
  void                SetMagneticFields();

  G4VPhysicalVolume*  SetupForNCEfficencyStudy();
  
//   void CreateBeamDump         ();
//   void CreateUshiwaka         (double ushiwaka_x,double ushiwaka_y, double ushiwaka_z);
//   void CreateNeutronCounter   ();
//   void CreateChargeVetoCounter();
//   void CreateProtonCounter    ();
//   void CreateDoraemon         (double cdsPos_x,  double cdsPos_y,   double cdsPos_z);
//   void CreateCDC              ();
//   void CreateCDH              ();
//   void CreateInnerHodo        ();
//   void CreateTargetSystem     (double Tgt_Pos_x=0.0, double Tgt_Pos_y=0.0, double Tgt_Pos_z=0.0);
//   void CreateBeamVeto         (double BVC_Pos_x=0.0, double BVC_Pos_y=0.0, double BVC_Pos_z=0.0);
//   void CreateFDC              (double FDC_Pos_x=0.0, double FDC_Pos_y=0.0, double FDC_Pos_z=0.0);
//   void CreateFDC2             (double FDC_Pos_x=0.0, double FDC_Pos_y=0.0, double FDC_Pos_z=0.0);
//   void CreateBLCs             ();
//   void CreateT0               ();
//   void CreateAC               ();
//   void CreateBPC              ();
//   void CreateBPD              ();

private:
  std::vector <Bool_t> vnCounterStat; 	// Detector ON/OFF status flags
					// id follows "gCounterID" in KnuclCommon.h 
  Double_t param[200][200][20];
  void SetGeomParameters(std::string filename);
  void ReadFile(std::string filename);
  G4LogicalVolume* ConstructShape(int cid);
  G4Transform3D MakeTrans(double* par, bool ZFIRST=false);
  G4CSGSolid* MakeShape(double* par, const char* name);

  void ConstructBeamDump();
  void ConstructTarget();
  void CreateCDCCells(G4LogicalVolume* log,bool WIRE);
  void CreateBLCCells(int cid,G4LogicalVolume* log,bool WIRE);
  void CreateFDCCells(int cid,G4LogicalVolume* log,bool WIRE);
  void CreateBPCCells(int cid,G4LogicalVolume* log,bool WIRE);
  void CreateChamberCells(int cid,G4LogicalVolume* log,bool WIRE=false);

public:
  Bool_t GetCounterStatus(int cid) const {return vnCounterStat[cid];}
  void   SetCounterStatus(int cid, Bool_t st){vnCounterStat[cid]=st;} 

private:

  KnuclCounterSD* counterSD;
  KnuclChamberSD* chamberSD;

  //G4UniformMagField*      magField;      // pointer to the magnetic field
  KnuclFieldSetup*        fEmFieldSetupUSHIWAKA;
  KnuclFieldSetup*        fEmFieldSetupCDC;

  KnuclFieldSetup_uswk*        fEmFieldSetupMap_uswk;
  KnuclFieldSetup_dora*        fEmFieldSetupMap_dora;
  KnuclFieldSetup_combi*       fEmFieldSetupMap_combi;  

  G4int                   DORAFieldMapFlag;
  G4int                   USWKFieldMapFlag;
  G4double                FieldInCDC;
  G4double                FieldInUshiwaka;
  
  G4bool                  neutron_efficiency_mode;  

  // World volume
  G4Box*             experimentalHall_box; 
  G4LogicalVolume*   experimentalHall_log;
  G4VPhysicalVolume* experimentalHall_phys;

  // CDS
  G4Tubs*            CDS_tube;
  G4LogicalVolume*   CDS_log;
  G4VPhysicalVolume* CDS_phys;
    
};

//----------------------------------//
// Color Definition
//----------------------------------//
const G4Color Black     = G4Colour::Black();
const G4Color Gray      = G4Colour::Gray();
const G4Color Red       = G4Colour::Red();
const G4Color Blue      = G4Colour::Blue();
const G4Color Green     = G4Colour::Green();
const G4Color Yellow    = G4Colour::Yellow();
const G4Color Orange    = G4Color(1,153/255,0);
const G4Color Purple    = G4Color(153/255,0,1);
const G4Color White     = G4Colour::White();
const G4Color Cyan      = G4Colour::Cyan();
const G4Color Magenta   = G4Colour::Magenta();

#endif

