#ifndef KnuclTrajectory_h
#define KnuclTrajectory_h 1

#include "G4VTrajectory.hh"
#include "G4Allocator.hh"
#include <stdlib.h>
#include "G4ThreeVector.hh"
#include "G4ios.hh"     
#include "globals.hh" 
#include "G4ParticleDefinition.hh" 
#include "G4TrajectoryPoint.hh"   
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4VProcess.hh"
#include <vector>

class G4Polyline;
class G4AttDef;
class G4AttValue;

typedef std::vector<G4VTrajectoryPoint*> KnuclTrajectoryPointContainer;

class KnuclTrajectory : public G4VTrajectory
{
 public:
   KnuclTrajectory();
   KnuclTrajectory(const G4Track* aTrack);
   KnuclTrajectory(KnuclTrajectory &);
   virtual ~KnuclTrajectory();

   inline void* operator new(size_t);
   inline void  operator delete(void*);
   inline int operator == (const KnuclTrajectory& right) const
   {return (this==&right);} 

   virtual G4int GetTrackID() const
   { return fTrackID; }
   virtual G4int GetParentID() const
   { return fParentID; }
   virtual G4String GetParticleName() const
   { return ParticleName; }
   virtual G4String GetCreatorProcess() const
   { return CreatorProcess; }
   virtual G4double GetCharge() const
   { return PDGCharge; }
   virtual G4int GetPDGEncoding() const
   { return PDGEncoding; }
   virtual G4ThreeVector GetInitialMomentum() const
   { return momentum; }
   virtual int GetPointEntries() const
   { return positionRecord->size(); }
   virtual G4VTrajectoryPoint* GetPoint(G4int i) const 
   { return (*positionRecord)[i]; }

   virtual void ShowTrajectory(std::ostream& os=G4cout) const;
   virtual void DrawTrajectory() const;
   virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
   virtual std::vector<G4AttValue>* CreateAttValues() const;
   virtual void AppendStep(const G4Step* aStep);
   virtual void MergeTrajectory(G4VTrajectory* secondTrajectory);

   G4ParticleDefinition* GetParticleDefinition();

   inline G4int GetTrackStatus() const
   { return fTrackStatus; }
   inline const G4ThreeVector& GetVertexPosition() const
   { return vertexPosition; }
   inline G4double GetGlobalTime() const
   { return globalTime; }

   void Extrapolate(G4ThreeVector &vector, G4double z);

  void  SetDrawTrajectoryMode(G4int mode){ DrawTrajectoryMode = mode; };
  G4int GetDrawTrajectoryMode(void){ return DrawTrajectoryMode; };

private:
   KnuclTrajectoryPointContainer*  positionRecord;
   G4int                        fTrackID;
   G4int                        fParentID;
   G4int                        fTrackStatus;
   G4ParticleDefinition*        fpParticleDefinition;
   G4String                     ParticleName;
   G4double                     PDGCharge;
   G4int                        PDGEncoding;
   G4ThreeVector                momentum;
   G4ThreeVector                vertexPosition;
   G4double                     globalTime;
   G4int                        DrawTrajectoryMode;
   G4String                     CreatorProcess;
};

extern G4Allocator<KnuclTrajectory> myTrajectoryAllocator;

inline void* KnuclTrajectory::operator new(size_t)
{
  void* aTrajectory;
  aTrajectory = (void*)myTrajectoryAllocator.MallocSingle();
  return aTrajectory;
}

inline void KnuclTrajectory::operator delete(void* aTrajectory)
{
  myTrajectoryAllocator.FreeSingle((KnuclTrajectory*)aTrajectory);
}

#endif

