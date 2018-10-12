#ifndef KnuclFieldSetup_combi_H
#define KnuclFieldSetup_combi_H

#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

class G4FieldManager;
class G4ChordFinder;
class G4Mag_UsualEqRhs;
class G4MagIntegratorStepper;

class sub_field_uswk
{
  // Storage space for the table
  vector< vector< vector< float > > > xField;
  vector< vector< vector< float > > > yField;
  vector< vector< vector< float > > > zField;
  // The dimensions of the table
  int nx,ny,nz; 
  // The physical limits of the defined region
  double minx, maxx, miny, maxy, minz, maxz;
  // The physical extent of the defined region
  double dx, dy, dz;
  double fXoffset, fYoffset, fZoffset;
  bool invertX, invertY, invertZ;

public:
  sub_field_uswk(const char*, double xOffset, double yOffset, double zOffset);
  // ~combi_field();
  void  GetFieldValue( const  double Point[4],
		       double *Bfield          ) const;
};

class sub_field_dora
{
  // Storage space for the table
  vector< vector< vector< float > > > xField;
  vector< vector< vector< float > > > yField;
  vector< vector< vector< float > > > zField;
  // The dimensions of the table
  int nx,ny,nz; 
  // The physical limits of the defined region
  double minx, maxx, miny, maxy, minz, maxz;
  // The physical extent of the defined region
  double dx, dy, dz;
  double fXoffset, fYoffset, fZoffset;
  bool invertX, invertY, invertZ;
  double norm_factor; 
public:
  sub_field_dora(const char*, double xOffset, double yOffset, double zOffset, double norm);
  // ~combi_field();
  void  GetFieldValue( const  double Point[4],
		       double *Bfield          ) const;
};

class combi_field
 : public G4MagneticField
{
private:
  sub_field_uswk *puswk;
  sub_field_dora *pdora;

public:
  combi_field(const char* uswk_file, 
	      double uswk_xOffset, double uswk_yOffset, double uswk_zOffset,
	      const char* dora_file, 
	      double dora_xOffset, double dora_yOffset, double dora_zOffset,
	      double dora_norm);
  // ~combi_field();
  void  GetFieldValue( const  double Point[4],
		       double *Bfield          ) const;
};

class KnuclFieldSetup_combi
{
public:
  KnuclFieldSetup_combi() ;               //  A zero field
  KnuclFieldSetup_combi(const char* uswk_file, 
	      double uswk_xOffset, double uswk_yOffset, double uswk_zOffset,
	      const char* dora_file, 
	      double dora_xOffset, double dora_yOffset, double dora_zOffset,
	      double dora_norm);

  ~KnuclFieldSetup_combi() ;

  void SetMinStep(G4double step) { fMinStep = step ; }

  G4FieldManager*  GetFieldManager() { return fFieldManager;}

protected:

  G4FieldManager*         fFieldManager ;
  G4ChordFinder*          fChordFinder ;
  G4Mag_UsualEqRhs*       fEquation ; 
  G4MagneticField*        fMagneticField ; 
  G4MagIntegratorStepper* fStepper ;

  G4double                fMinStep ;

};


#endif
