#include "KnuclTools.h"

std::string PDGToName(const int &pdg_id)
{
  if( pdg_id==11 ) return "e";
  else if( pdg_id==-11 ) return "e+";
  else if( pdg_id==2112 ) return "neutron";
  else if( pdg_id==2212 ) return "proton";
  else if( pdg_id== 211 ) return "pion+";
  else if( pdg_id==-211 ) return "pion-";

  return "unknown";
}

int ProcessNameToProcessID(const std::string &name)
{
  if( name=="" ) return 0; //Initial
  else if( name=="Decay"                 ) return 1;
  else if( name=="conv"                  ) return 2;
  else if( name=="Transportation"        ) return 3;
  else if( name=="phot"                  ) return 4;
  else if( name=="annihil"               ) return 5;
  else if( name=="compt"                 ) return 6;
  else if( name=="eBrem"                 ) return 7;
  else if( name=="hadElastic"            ) return 8;
  else if( name=="CoulombScat"           ) return 9;
  else if( name=="nKiller"               ) return 10;
  else if( name=="photonNuclear"         ) return 11;
  else if( name=="msc"                   ) return 12;
  else if( name=="pi-Inelastic"          ) return 100;
  else if( name=="pi+Inelastic"          ) return 101;
  else if( name=="kaon-Inelastic"        ) return 102;
  else if( name=="kaon+Inelastic"        ) return 103;
  else if( name=="kaon0LInelastic"       ) return 104;
  else if( name=="kaon0SInelastic"       ) return 105;
  else if( name=="lambdaInelastic"       ) return 106;
  else if( name=="sigma+Inelastic"       ) return 107;
  else if( name=="sigma-Inelastic"       ) return 108;
  else if( name=="sigma0Inelastic"       ) return 109;
  else if( name=="protonInelastic"       ) return 110;
  else if( name=="neutronInelastic"      ) return 111;
  else if( name=="dInelastic"            ) return 112;
  else if( name=="tInelastic"            ) return 113;
  else if( name=="He3Inelastic"          ) return 114;
  else if( name=="alphaInelastic"        ) return 115;
  else if( name.find("Inelastic")!=std::string::npos ) return 199;
  else if( name=="eIoni"                 ) return 200;
  else if( name=="hIoni"                 ) return 201;
  else if( name=="ionIoni"               ) return 202;
  else if( name=="muIoni"                ) return 203;
  else if( name=="hBertiniCaptureAtRest" ) return 204;
  else if( name=="nCapture"              ) return 205;
  else if( name=="muMinusCaptureAtRest"  ) return 206;
  else if( name=="unknown"               ) return -999;
  else return -1;
}

std::string ProcessIDToProcessName(const int &id)
{
  if( id==0                 ) return "initial";
  else if( id==1            ) return "Decay";
  else if( id==2            ) return "conv";
  else if( id==3            ) return "Transportation";
  else if( id==4            ) return "phot";
  else if( id==5            ) return "annihil";
  else if( id==6            ) return "compt";
  else if( id==7            ) return "eBrem";
  else if( id==8            ) return "hadElastic";
  else if( id==9            ) return "CoulombScat";
  else if( id==10           ) return "nKiller";
  else if( id==11           ) return "photoNuclear";
  else if( id==12           ) return "msc";
  else if( id==100          ) return "pi-Inelastic";
  else if( id==101          ) return "pi+Inelastic";
  else if( id==102          ) return "kaon-Inelastic";
  else if( id==103          ) return "kaon+Inelastic";
  else if( id==104          ) return "kaon0Inelastic";
  else if( id==105          ) return "kaon0SInelastic";
  else if( id==106          ) return "lambdaInelastic";
  else if( id==107          ) return "sigma-Inelastic";
  else if( id==108          ) return "sigma+Inelastic";
  else if( id==109          ) return "sigma0Inelastic";
  else if( id==110          ) return "protonInelastic";
  else if( id==111          ) return "neutronInelastic";
  else if( id==112          ) return "dInelastic";
  else if( id==113          ) return "tInelastic";
  else if( id==114          ) return "He3Inelastic";
  else if( id==115          ) return "alphaInelastic";
  else if( 100<id && id<200 ) return "Inelastic";
  else if( id==200          ) return "eIoni";
  else if( id==201          ) return "hIoni";
  else if( id==202          ) return "ionIoni";
  else if( id==203          ) return "muIoni";
  else if( id==204          ) return "hBertiniCaptureAtRest";
  else if( id==205          ) return "nCapture";
  else if( id==206          ) return "muMinusCaptureAtRest";
  else if( id==-999         ) return "unknown";
  else return "error";
}
