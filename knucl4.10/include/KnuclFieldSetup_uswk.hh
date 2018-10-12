#ifndef KnuclFieldSetup_uswk_H
#define KnuclFieldSetup_uswk_H

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

class uswk_field
 : public G4MagneticField
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
  uswk_field(const char*, double xOffset, double yOffset, double zOffset);
  // ~uswk_field();
  void  GetFieldValue( const  double Point[4],
		       double *Bfield          ) const;
};

class KnuclFieldSetup_uswk
{
public:
  KnuclFieldSetup_uswk() ;               //  A zero field
  KnuclFieldSetup_uswk(G4ThreeVector) ;  //  The value of the field
  KnuclFieldSetup_uswk(const char*, double, double, double) ;  //  The implementation for Uswk field map
  ~KnuclFieldSetup_uswk() ;

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
