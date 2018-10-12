#ifndef KnuclTrackingAction_h 
#define KnuclTrackingAction_h 1

#include "G4UserTrackingAction.hh"

class KnuclTrackingAction : public G4UserTrackingAction {
  public:
    KnuclTrackingAction(){}
    virtual ~KnuclTrackingAction(){}

    virtual void PreUserTrackingAction(const G4Track*);
  //virtual void PostUserTrackingAction(const G4Track*);
};

#endif
