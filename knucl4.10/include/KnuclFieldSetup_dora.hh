#ifndef KnuclFieldSetup_dora_H
#define KnuclFieldSetup_dora_H

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

class dora_field
 : public G4MagneticField
{
  
  // Storage space for the table
  vector< vector< vector< float > > > xField;
  vector< vector< vector< float > > > yField;
  vector< vector< vector< float > > > zField;
  // The dimensions of the table
  int nx,ny,nz; 
  double norm_factor; 
  // The physical limits of the defined region
  double minx, maxx, miny, maxy, minz, maxz;
  // The physical extent of the defined region
  double dx, dy, dz;
  double fXoffset, fYoffset, fZoffset;
  bool invertX, invertY, invertZ;

public:
  dora_field(const char*, double xOffset, double yOffset, double zOffset, double norm_factor);
  // ~dora_field();
  void  GetFieldValue( const  double Point[4],
		       double *Bfield          ) const;
};

class KnuclFieldSetup_dora
{
public:
  KnuclFieldSetup_dora() ;               //  A zero field
  KnuclFieldSetup_dora(G4ThreeVector) ;  //  The value of the field
  KnuclFieldSetup_dora(const char*, double, double, double, double) ;  //  The implementation for Dora field map
  ~KnuclFieldSetup_dora() ;

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
