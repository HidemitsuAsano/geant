#include "KnuclMaterialManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
// ====================================================================
//    KnuclMaterialManager.cc
//
//                                 F.Sakuma May,2007
// ====================================================================

//#include "KnuclMaterialManager.hh"

// ====================================================================
//
// class description
//
// ====================================================================

//KnuclMaterialManager* KnuclMaterialManager::ptrMaterialManager= 0;

////////////////////////////////////////////
KnuclMaterialManager::KnuclMaterialManager()
////////////////////////////////////////////
{
  elementTable=  G4Element::GetElementTable();
  materialTable= G4Material::GetMaterialTable();
  DefineMaterials();

  // list elements/materials
  G4cout << *elementTable << G4endl;
  G4cout << *materialTable << G4endl;
}

/////////////////////////////////////////////
KnuclMaterialManager::~KnuclMaterialManager()
/////////////////////////////////////////////
{
}

///////////////////////////////////////////////////////////////////////
G4Element* KnuclMaterialManager::GetElement(const G4String& name) const
///////////////////////////////////////////////////////////////////////
{
  G4Element* element= 0;

  G4int i;
  for(i=0; i< (G4int)elementTable-> size(); i++) {
    G4Element* aElement= (*elementTable)[i];
    G4String aName= aElement-> GetName();
    if(aName==name) {
      element= aElement;
      break;
    }
  }

  if(!element) {
    G4String errorMessage= " Element named <" + name + "> is not defined.";
    //new G4Exception definition y.ma@riken.jp 2012/11/08
    G4Exception("KnuclMaterialManager::GetElement",
		"elements not defined",
		FatalException,
		errorMessage);

    //    G4Exception(errorMessage);
    //G4Exception("KnuclMaterialManager::GetElement()",
    //"Invalid Element", FatalException,
    //errorMessage);
  }

  return element;
}


/////////////////////////////////////////////////////////////////////////
G4Material* KnuclMaterialManager::GetMaterial(const G4String& name) const
/////////////////////////////////////////////////////////////////////////
{
  G4Material* material= 0;

  G4int i;
  for(i=0; i< (G4int)materialTable-> size(); i++) {
    G4Material* aMaterial= (*materialTable)[i];
    G4String aName= aMaterial-> GetName();
    if(aName==name) {
      material= aMaterial;
      break;
    }
  }

  if(!material) {
    G4String errorMessage= " Material named <" + name + "> is not defined.";
    //new G4Exception definition y.ma@riken.jp 2012/11/08
    G4Exception("KnuclMaterialManager::GetMaterial",
		"materials not defined",
		FatalException,
		errorMessage);

    //    G4Exception(errorMessage);
    //G4Exception("KnuclMaterialManager::GetMaterial()",
    //"Invalid Material", FatalException,
    //errorMessage);
  } 

  return material;
}

// define elements/materials used in a user program
#include "KnuclMaterials.icc"
