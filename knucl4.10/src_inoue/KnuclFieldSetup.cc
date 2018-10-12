#include "KnuclFieldSetup.hh"

#include "G4UniformMagField.hh"
#include "G4MagneticField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ChordFinder.hh"

#include "G4SystemOfUnits.hh"

#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4ClassicalRK4.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4CashKarpRKF45.hh"
#include "G4RKG3_Stepper.hh"

KnuclFieldSetup::KnuclFieldSetup()
{
}

KnuclFieldSetup::KnuclFieldSetup(G4ThreeVector fieldVector)
 : fChordFinder(0), fStepper(0)
{    
  G4cout << "#### " << fieldVector.x() << " " << fieldVector.y() << " " << fieldVector.z() << G4endl;

  fMinStep     = 0.25*mm ;

  fMagneticField = new G4UniformMagField(fieldVector);
  fEquation      = new G4Mag_UsualEqRhs(fMagneticField);
  fStepper       = new G4ClassicalRK4( fEquation );
  fChordFinder   = new G4ChordFinder( fMagneticField, fMinStep, fStepper);

  fFieldManager  = new G4FieldManager(fMagneticField,fChordFinder);
}

KnuclFieldSetup::~KnuclFieldSetup()
{
  if(fMagneticField) delete fMagneticField;
  if(fChordFinder)   delete fChordFinder;
  if(fStepper)       delete fStepper;
  if(fEquation)      delete fEquation;
  G4cout<<"KnuclFieldSetup::~KnuclFieldSetup()"<<G4endl;
  //  getchar();
}

