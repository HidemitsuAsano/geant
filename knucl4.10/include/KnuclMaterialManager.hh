// ====================================================================
//    KnuclMaterialManager.hh
//    
//                                 F.Sakuma May,2007
// ====================================================================
#ifndef KNUCL_MATERIAL_MANAGER_H
#define KNUCL_MATERIAL_MANAGER_H

#include "G4Element.hh"
#include "G4Material.hh"

// ====================================================================
//
// class definition
//
// ====================================================================
 

class KnuclMaterialManager {
private:
  // this object is a singleton.
  //static KnuclMaterialManager* ptrMaterialManager;

  //KnuclMaterialManager();  // constructor should NOT be public.

  const G4ElementTable*  elementTable;
  const G4MaterialTable* materialTable;

  // define elements/materials used in a user program
  void DefineMaterials();  
  
public:
  KnuclMaterialManager();
  ~KnuclMaterialManager();

  //static KnuclMaterialManager* GetPointer();

  G4Element* GetElement(const G4String& name) const;
  G4Material* GetMaterial(const G4String& name) const;

};

// ====================================================================
// inline definition
// ====================================================================

//inline PhiG4MaterialManager* PhiG4MaterialManager::GetPointer()
//{
//  if(!ptrMaterialManager) ptrMaterialManager= new PhiG4MaterialManager;
//  return ptrMaterialManager;
//}

#endif
