#ifndef KNUCLCOMMON_H
#define KNUCLCOMMON_H 1

#include "TVector3.h"
#include "G4ThreeVector.hh"

#define XYZ 3
#define XCOORD 0
#define YCOORD 1
#define ZCOORD 2

#define DecayToLambda
#define DecayToPiSigma


//### constant ###//
const Double_t Const=0.299792458; // =c/10^9


//### constant for beam ###//
static const Double_t BEAM_MOM_BITE = 2.5;//(percent,FWHM)
const Double_t CDC_OFFSET = 3.0; // unit :: cell
const int N_CDH = 36;
const int N_CVC = 34;
const int N_PC  = 27;
//const Double_t BPC_POSZ =-197.37;// 159.77+37.6 //mm 
// Counter ID (from k18ana/src/GlobalVariables.h)
enum gCounterID { CID_CDC     = 0,
		  CID_CDH     = 1,
		  CID_BHD     = 2,
		  //		  CID_PA      = 3,
		  CID_T0      = 4,
		  CID_E0      = 5,
		  CID_DEF     = 5,
		  //		  CID_B1      = 6, //
		  CID_LC1     = 7,
		  CID_LC2     = 8,
		  CID_AC      = 9,
		  CID_WC      = 10,
		  CID_GC      = 11,
		  //		  CID_Range   = 12,
		  //		  CID_B2      = 13,
		  CID_TOFstop = 14,
		  CID_CVC     = 14,
		  //		  CID_PDC1    = 15,
		  CID_BLC1a   = 15,
		  //		  CID_PDC2    = 16,
		  CID_BLC1b   = 16,
		  CID_BLC2a   = 17,
		  CID_BLC2b   = 18,
		  CID_SDD     = 19,
		  CID_BLC1    = 21,
		  CID_BLC2    = 22,
		  CID_FDC1    = 23,
		  CID_FDC2    = 24,
		  CID_ZVC     = 30,
		  CID_KDV     = 31,
		  CID_NC      = 32,
		  CID_BVC     = 33,
		  CID_PC      = 35,
		  CID_Longbar = 36,
		  CID_LB      = 36,
		  CID_WVC     = 37,
		  CID_BPC     = 40,
		  CID_BPD     = 41,
		  CID_IH      = 42,
		  CID_T0pre   = 51,
		  CID_T0post  = 52,
		  CID_BHDpost = 56,
		  CID_HVC1    = 61,
		  CID_HVC2    = 62,
		  CID_BHDmul  = 81,
		  CID_T0mul   = 82,
		  CID_BVCmul  = 83,
		  CID_HVC1mul = 84,
		  CID_HVC2mul = 85,
		  CID_REFmul  = 86,
		  CID_BD      = 90,
		  //		  CID_BDC     = 90,
		  CID_TEMP1   = 91,
		  CID_TEMP2   = 92,
		  CID_GPIO    = 97,
		  CID_MISC    = 98,
		  CID_TEMP    = 99,
		  CID_Hall    =100,
		  CID_Floor   =101,
		  CID_BeamDump=110,
		  CID_SideDump=111,
		  CID_NShield =112,
		  CID_SideCon =113,
		  CID_DoorCon =114,
		  CID_Doraemon=120,
		  CID_USWK    =121,
		  CID_CDCCFRP =130,
		  CID_CDCMylar=131,
		  CID_TarChm  =140,
		  CID_RadS    =141,
		  CID_TarCFRP =142,
		  CID_TarCap  =143,
		  CID_TarRing =144,
		  CID_TarCell =150,
		  CID_Target  =151,
		  CID_CellBe  =152,
		  CID_CellAlBe=153,
		  CID_BShield  = 154,
		  CID_BFrange  = 155,
		  CID_CellRing=156,
		  CID_Fiducial = 160 
};

//### utilities ###//
TVector3 ConvVecGT(const G4ThreeVector& v);
G4ThreeVector ConvVecTG(const TVector3& v);

inline TVector3 ConvVecGT(const G4ThreeVector& v)
{
  return TVector3(v.x(),v.y(),v.z());
}

inline G4ThreeVector ConvVecTG(const TVector3& v)
{
  return G4ThreeVector(v.x(),v.y(),v.z());
}


#endif
