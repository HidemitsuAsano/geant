#include "KnuclTrackingAction.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "KnuclTrajectory.hh"

void KnuclTrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
  // Create trajectory only for track in tracking region

  fpTrackingManager->SetStoreTrajectory(true);
  fpTrackingManager->SetTrajectory(new KnuclTrajectory(aTrack));

}

//void KnuclTrackingAction::PostUserTrackingAction(const G4Track* aTrack)
//{}


