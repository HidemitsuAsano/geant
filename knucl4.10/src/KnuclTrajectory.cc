#include "KnuclTrajectory.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleTypes.hh"
#include "G4Polyline.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UIcommand.hh"
#include "G4VisAttributes.hh"
#include "G4VVisManager.hh"
#include "G4UnitsTable.hh"
#include "G4DynamicParticle.hh"
#include "G4PrimaryParticle.hh"
#include "KnuclDetectorConstruction.hh"

G4Allocator<KnuclTrajectory> myTrajectoryAllocator;

KnuclTrajectory::KnuclTrajectory()
:G4VTrajectory()
{
   fpParticleDefinition = 0;
   ParticleName = "";
   CreatorProcess = "";
   PDGCharge = 0;
   PDGEncoding = 0;
   fTrackID = 0;
   fParentID = 0;
   fTrackStatus = 0;
   positionRecord = 0;
   momentum = G4ThreeVector(0.,0.,0.);
   vertexPosition = G4ThreeVector(0.,0.,0.);
   globalTime = 0.;
   DrawTrajectoryMode = 0;
}

KnuclTrajectory::KnuclTrajectory(const G4Track* aTrack)
:G4VTrajectory()
{
   fpParticleDefinition = aTrack->GetDefinition();
   ParticleName         = fpParticleDefinition->GetParticleName();
   PDGCharge            = fpParticleDefinition->GetPDGCharge();
   PDGEncoding          = fpParticleDefinition->GetPDGEncoding();
   const G4VProcess *pro= aTrack->GetCreatorProcess();
   if(pro){
     CreatorProcess = pro->GetProcessName();
   }
   if(ParticleName=="unknown")
   {
     G4PrimaryParticle*pp = aTrack->GetDynamicParticle()->GetPrimaryParticle();
     if(pp)
     {
       if(pp->GetCharge()<DBL_MAX) PDGCharge = pp->GetCharge();
       PDGEncoding = pp->GetPDGcode();
       if(pp->GetG4code()!=0)
       {
         ParticleName += " : ";
         ParticleName += pp->GetG4code()->GetParticleName();
       }
     }
   }

   fTrackID = aTrack->GetTrackID();
   fParentID = aTrack->GetParentID();
   positionRecord = new KnuclTrajectoryPointContainer();
   positionRecord->push_back(new G4TrajectoryPoint(aTrack->GetPosition()));
   momentum       = aTrack->GetMomentum();
   vertexPosition = aTrack->GetPosition();
   globalTime     = aTrack->GetGlobalTime();
   DrawTrajectoryMode = 0;
}

KnuclTrajectory::KnuclTrajectory(KnuclTrajectory & right)
:G4VTrajectory()
{
  ParticleName = right.ParticleName;
  fpParticleDefinition = right.fpParticleDefinition;
  PDGCharge = right.PDGCharge;
  PDGEncoding = right.PDGEncoding;
  fTrackID = right.fTrackID;
  fParentID = right.fParentID;
  fTrackStatus = right.fTrackStatus;
  positionRecord = new KnuclTrajectoryPointContainer();
  for(size_t i=0;i<right.positionRecord->size();i++)
  {
    G4TrajectoryPoint* rightPoint = (G4TrajectoryPoint*)((*(right.positionRecord))[i]);
    positionRecord->push_back(new G4TrajectoryPoint(*rightPoint));
  }
   momentum = right.momentum;
   vertexPosition = right.vertexPosition;
   globalTime = right.globalTime;
   DrawTrajectoryMode = right.DrawTrajectoryMode;
}

KnuclTrajectory::~KnuclTrajectory()
{
  size_t i;
  for(i=0;i<positionRecord->size();i++){
    delete  (*positionRecord)[i];
  }
  positionRecord->clear();

  delete positionRecord;
}

void KnuclTrajectory::ShowTrajectory(std::ostream& os) const
{
//   os << G4endl << "TrackID =" << fTrackID 
//        << " : ParentID=" << fParentID << " : TrackStatus=" << fTrackStatus << G4endl;
//   os << "Particle name : " << ParticleName << "  PDG code : " << PDGEncoding
//        << "  Charge : " << PDGCharge << G4endl;
//   os << "Original momentum : " <<
//        G4BestUnit(momentum,"Energy") << G4endl;
//   os << "Vertex : " << G4BestUnit(vertexPosition,"Length")
//        << "  Global time : " << G4BestUnit(globalTime,"Time") << G4endl;
//   os << "  Current trajectory has " << positionRecord->size() 
//        << " points." << G4endl;

   os << "TrackID =" << fTrackID << " : " << fParentID << " : " << ParticleName << " " << G4BestUnit(vertexPosition,"Length") << G4endl;

   //for( size_t i=0 ; i < positionRecord->size() ; i++){
   //    G4TrajectoryPoint* aTrajectoryPoint = (G4TrajectoryPoint*)((*positionRecord)[i]);
   //    os << "Point[" << i << "]" 
   //         << " Position= " << aTrajectoryPoint->GetPosition() << G4endl;
   //}
}

void KnuclTrajectory::DrawTrajectory() const
{
  //For k1.8 simulation, mode=0;
  //For neutron efficiency simulation, mode=1;
  int mode=0;

  if(mode==0){
    G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
    G4ThreeVector pos;
    
    G4Polyline pPolyline;
    
    for (size_t i = 0; i < positionRecord->size() ; i++) {
      G4TrajectoryPoint* aTrajectoryPoint = (G4TrajectoryPoint*)((*positionRecord)[i]);
      pos = aTrajectoryPoint->GetPosition();
      if (pos.z()<2.0*m) 
	pPolyline.push_back( pos );
    }
    G4Colour colour(White);
    if(fpParticleDefinition==G4Gamma::GammaDefinition())
      colour = G4Colour(Gray);
    else if(fpParticleDefinition==G4Electron::ElectronDefinition()
	    ||fpParticleDefinition==G4Positron::PositronDefinition())
      colour = G4Colour(Yellow);
    else if(fpParticleDefinition==G4MuonMinus::MuonMinusDefinition()
	    ||fpParticleDefinition==G4MuonPlus::MuonPlusDefinition())
      colour = G4Colour(Green);
    else if(fpParticleDefinition->GetParticleType()=="meson")
      {
	if(PDGCharge>0.)
	  colour = G4Colour(Red);
	else if(PDGCharge<0.)
	  colour = G4Colour(Blue);
	else
	  colour = G4Colour(Gray);
      }
    else if(fpParticleDefinition->GetParticleType()=="baryon")
      {
	if(PDGCharge>0.)
	  colour = G4Colour(Magenta);
	else if(PDGCharge<0.)
	  colour = G4Colour(Cyan);
	else
	  if (PDGEncoding==2112) { 
	    colour = G4Colour(Orange);
	  } else {
	    colour = G4Colour(Cyan);
	  }
      } else {
      G4cout <<  ParticleName << G4endl;
    }
    
    G4cout <<  ParticleName << G4endl;
    G4VisAttributes attribs(colour);
    attribs.SetLineWidth(3.0);
    pPolyline.SetVisAttributes(attribs);
    // except for geantino, neiutrino
    if(!(fpParticleDefinition->GetParticleName()=="geantino" ||
	 fpParticleDefinition->GetParticleName()=="chargedgeantino" ||
	 fpParticleDefinition==G4NeutrinoMu::NeutrinoMuDefinition() ||
	 fpParticleDefinition==G4AntiNeutrinoMu::AntiNeutrinoMuDefinition() ||
	 fpParticleDefinition==G4NeutrinoE::NeutrinoEDefinition() ||
	 fpParticleDefinition==G4AntiNeutrinoE::AntiNeutrinoEDefinition()) )
      if(pVVisManager) pVVisManager->Draw(pPolyline);
  }
  else if(mode==1){
    G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
    G4ThreeVector pos;
    
    G4Polyline pPolyline;
    
    for (size_t i = 0; i < positionRecord->size() ; i++) {
      G4TrajectoryPoint* aTrajectoryPoint = (G4TrajectoryPoint*)((*positionRecord)[i]);
      pos = aTrajectoryPoint->GetPosition();
      //      if (pos.z()<2.0*m) 
	pPolyline.push_back( pos );
    }
    G4Colour colour(White);
    if(fpParticleDefinition==G4Gamma::GammaDefinition())
      colour = G4Colour(Gray);
    else if(fpParticleDefinition==G4Electron::ElectronDefinition()
	    ||fpParticleDefinition==G4Positron::PositronDefinition())
      colour = G4Colour(Yellow);
    else if(fpParticleDefinition==G4MuonMinus::MuonMinusDefinition()
	    ||fpParticleDefinition==G4MuonPlus::MuonPlusDefinition())
      colour = G4Colour(Green);
    else if(fpParticleDefinition->GetParticleType()=="meson")
      {
	if(PDGCharge>0.)
	  colour = G4Colour(Red);
	else if(PDGCharge<0.)
	  colour = G4Colour(Blue);
	else
	  colour = G4Colour(Gray);
      }
    else if(fpParticleDefinition->GetParticleType()=="baryon")
      {
	if(PDGCharge>0.)
	  colour = G4Colour(Magenta);
	else if(PDGCharge<0.)
	  colour = G4Colour(Cyan);
	else
	  if (PDGEncoding==2112) { 
	    colour = G4Colour(Orange);
	  } else {
	    colour = G4Colour(Cyan);
	  }
      } else {
      G4cout <<  ParticleName << G4endl;
    }
    
    G4cout <<  ParticleName << G4endl;
    G4VisAttributes attribs(colour);
    attribs.SetLineWidth(1.0);
    pPolyline.SetVisAttributes(attribs);
    // except for geantino, neiutrino
    if(!(fpParticleDefinition->GetParticleName()=="geantino" ||
	 fpParticleDefinition->GetParticleName()=="chargedgeantino" ||
	 fpParticleDefinition==G4NeutrinoMu::NeutrinoMuDefinition() ||
	 fpParticleDefinition==G4AntiNeutrinoMu::AntiNeutrinoMuDefinition() ||
	 fpParticleDefinition==G4NeutrinoE::NeutrinoEDefinition() ||
	 fpParticleDefinition==G4AntiNeutrinoE::AntiNeutrinoEDefinition()) )
      if(pVVisManager) pVVisManager->Draw(pPolyline);
  }
}

const std::map<G4String,G4AttDef>* KnuclTrajectory::GetAttDefs() const
{
  G4bool isNew;
  std::map<G4String,G4AttDef>* store
    = G4AttDefStore::GetInstance("KnuclTrajectory",isNew);
  if (isNew) {

    G4String ID("ID");
    (*store)[ID] = G4AttDef(ID,"Track ID","Bookkeeping","","G4int");

    G4String PID("PID");
    (*store)[PID] = G4AttDef(PID,"Parent ID","Bookkeeping","","G4int");

    G4String Status("Status");
    (*store)[Status] = G4AttDef(Status,"Track Status","Bookkeeping","","G4int");

    G4String PN("PN");
    (*store)[PN] = G4AttDef(PN,"Particle Name","Bookkeeping","","G4String");

    G4String Ch("Ch");
    (*store)[Ch] = G4AttDef(Ch,"Charge","Physics","e+","G4double");

    G4String PDG("PDG");
    (*store)[PDG] = G4AttDef(PDG,"PDG Encoding","Bookkeeping","","G4int");

    G4String IMom("IMom");
    (*store)[IMom] = G4AttDef(IMom, "Momentum of track at start of trajectory",
			      "Physics","G4BestUnit","G4ThreeVector");

    G4String IMag("IMag");
    (*store)[IMag] = 
      G4AttDef(IMag, "Magnitude of momentum of track at start of trajectory",
	       "Physics","G4BestUnit","G4double");

    G4String VtxPos("VtxPos");
    (*store)[VtxPos] = G4AttDef(VtxPos, "Vertex position",
			      "Physics","G4BestUnit","G4ThreeVector");

    G4String NTP("NTP");
    (*store)[NTP] = G4AttDef(NTP,"No. of points","Bookkeeping","","G4int");

  }
  return store;
}

std::vector<G4AttValue>* KnuclTrajectory::CreateAttValues() const
{
  std::vector<G4AttValue>* values = new std::vector<G4AttValue>;

  values->push_back
    (G4AttValue("ID",G4UIcommand::ConvertToString(fTrackID),""));

  values->push_back
    (G4AttValue("PID",G4UIcommand::ConvertToString(fParentID),""));

  values->push_back
    (G4AttValue("Status",G4UIcommand::ConvertToString(fTrackStatus),""));

  values->push_back(G4AttValue("PN",ParticleName,""));

  values->push_back
    (G4AttValue("Ch",G4UIcommand::ConvertToString(PDGCharge),""));

  values->push_back
    (G4AttValue("PDG",G4UIcommand::ConvertToString(PDGEncoding),""));

  values->push_back
    (G4AttValue("IMom",G4BestUnit(momentum,"Energy"),""));

  values->push_back
    (G4AttValue("IMag",G4BestUnit(momentum.mag(),"Energy"),""));

  values->push_back
    (G4AttValue("VtxPos",G4BestUnit(vertexPosition,"Length"),""));

  values->push_back
    (G4AttValue("NTP",G4UIcommand::ConvertToString(GetPointEntries()),""));

  return values;
}

void KnuclTrajectory::AppendStep(const G4Step* aStep)
{
   positionRecord->push_back( new G4TrajectoryPoint(aStep->GetPostStepPoint()->
                                 GetPosition() ));
}
  
G4ParticleDefinition* KnuclTrajectory::GetParticleDefinition()
{
   return (G4ParticleTable::GetParticleTable()->FindParticle(ParticleName));
}

void KnuclTrajectory::MergeTrajectory(G4VTrajectory* secondTrajectory)
{
  if(!secondTrajectory) return;

  KnuclTrajectory* seco = (KnuclTrajectory*)secondTrajectory;
  G4int ent = seco->GetPointEntries();
  for(int i=1;i<ent;i++) // initial point of the second trajectory should not be merged
  {
    positionRecord->push_back((*(seco->positionRecord))[i]);
  }
  delete (*seco->positionRecord)[0];
  seco->positionRecord->clear();

}

void KnuclTrajectory::Extrapolate(G4ThreeVector &vector, G4double z)
{
  vector.set(-999.0,-999.0,-999.0);

  if ( positionRecord->size()>3 ) {

    G4ThreeVector Position[4]; 

    int icount=0;
    for( size_t i=positionRecord->size() ; i > (positionRecord->size()-4) ; i--){
      G4TrajectoryPoint* aTrajectoryPoint = (G4TrajectoryPoint*)((*positionRecord)[i-1]);
      Position[icount] = aTrajectoryPoint->GetPosition(); icount++;
    }
    
    G4int n = 4;
    G4double sxz = 0, syz = 0;
    G4double sz  = 0;
    G4double sx  = 0, sy  = 0;
    G4double sxx = 0, syy = 0;
    G4double szz = 0;
   
    for (G4int i=0; i<4; i++){
      sxz += Position[i].x()*Position[i].z();
      syz += Position[i].y()*Position[i].z();
      sz  += Position[i].z();
      sx  += Position[i].x();
      sy  += Position[i].y();
      sxx += Position[i].x()*Position[i].x();
      syy += Position[i].y()*Position[i].y();
      szz += Position[i].z()*Position[i].z();
    }

    G4double ax = (n*sxz  - sz  * sx)  /(n * szz-sz*sz); 
    G4double ay = (n*syz  - sz  * sy)  /(n * szz-sz*sz); 
    G4double bx = (szz*sx - sxz * sz)/(n * szz-sz*sz);
    G4double by = (szz*sy - syz * sz)/(n * szz-sz*sz);

    //G4double ax = (Position[0].x()-Position[1].x())/(Position[0].z()-Position[1].z());
    //G4double ay = (Position[0].y()-Position[1].y())/(Position[0].z()-Position[1].z());
    //G4double bx = Position[0].x()-ax*Position[0].z();
    //G4double by = Position[0].y()-ay*Position[0].z();

    vector.setX(ax*z+bx);
    vector.setY(ay*z+by);
    vector.setZ(z);
  }
}

