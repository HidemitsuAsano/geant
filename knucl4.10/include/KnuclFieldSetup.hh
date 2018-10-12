#ifndef KnuclFieldSetup_H
#define KnuclFieldSetup_H

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

class KnuclFieldSetup
{
public:
  KnuclFieldSetup() ;               //  A zero field
  KnuclFieldSetup( G4ThreeVector ) ;  //  The value of the field
  ~KnuclFieldSetup() ;

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
