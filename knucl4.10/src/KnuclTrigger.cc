#include "KnuclTrigger.hh"

bool trigger(const std::vector<DetectorHit> &hits, G4String name)
{
  if( name.empty() ) return true;

  //  std::cout<<"Trigger : "<<name<<std::endl;
  bool trigNC=false;
  if( name.find("NC")!=G4String::npos ) trigNC=true;
  bool trigPC=false;
  if( name.find("PC")!=G4String::npos ) trigPC=true;
  bool trigCVC=false;
  if( name.find("CVC")!=G4String::npos ) trigCVC=true;
  if( name.find("FWD")!=G4String::npos ){ trigNC=true; trigPC=true; trigCVC=true; }

  bool trigCDH1=false;
  if( name.find("CDH1")!=G4String::npos ) trigCDH1=true;
  bool trigCDH2=false;
  if( name.find("CDH2")!=G4String::npos ) trigCDH2=true;
  bool trigCDH3=false;
  if( name.find("CDH3")!=G4String::npos ) trigCDH3=true;

  bool trigBPD1=false;
  if( name.find("BPD1")!=G4String::npos ) trigBPD1=true;
  bool trigBPD2=false;
  if( name.find("BPD2")!=G4String::npos ) trigBPD2=true;

  if( !trigNC && !trigPC && !trigCVC && !trigCDH1 && !trigCDH2 && !trigCDH3 && !trigBPD1 && !trigBPD2 ) return true;

  int nCDH=0;
  int nBPD=0;
  for( int i=0; i<hits.size(); i++ ){

    if( hits[i].detectorID()==CID_NC && trigNC ) return true;
    if( hits[i].detectorID()==CID_PC && trigPC ) return true;
    if( hits[i].detectorID()==CID_CVC && trigCVC ) return true;
    if( hits[i].detectorID()==CID_CDH ){
      nCDH++;
      //      std::cout<<"   CDH hit   pdg="<<std::setw(5)<<hits[i].pdg()<<"  dE : "<<hits[i].adc()<<std::endl;
    }
    if( hits[i].detectorID()==CID_BPD ) nBPD++;
  }

  // std::cout<<" nCDH : "<<nCDH<<std::endl;
  // std::cout<<" nBPD : "<<nBPD<<std::endl;
  // std::cout<<" NC   :"<<trigNC<<std::endl;
  // std::cout<<" CVC  :"<<trigCVC<<std::endl;
  // std::cout<<" PC   :"<<trigPC<<std::endl;

  if( nCDH>=3 && trigCDH3 ) return true;
  if( nCDH>=2 && trigCDH2 ) return true;
  if( nCDH>=1 && trigCDH1 ) return true;

  if( nBPD>=2 && trigBPD2 ) return true;
  if( nBPD>=1 && trigBPD1 ) return true;

  //  std::cout<<" Trigger false"<<std::endl;
  return false;
}
