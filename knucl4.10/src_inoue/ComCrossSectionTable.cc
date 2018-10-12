#include "ComCrossSectionTable.hh"

//-----------------------------------------------//
// class CorssSection
//-----------------------------------------------//
//////////////////////////////////////////////////////
void CrossSection::PrintCS()
//////////////////////////////////////////////////////
{
  std::cout<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<std::endl;
  std::cout<<"  reaction ID:               "<<fId<<std::endl;
  std::cout<<"  CrossSection [mb]:         "<<fCs<<std::endl;
  std::cout<<"  CS factor:                 "<<fCsFactor<<std::endl;
  std::cout<<"  initial PDG code:          "<<fInitPdg[0]<<", "<<fInitPdg[1]<<", "<<fInitPdg[2]<<std::endl;
  std::cout<<"  beam momentum [MeV/c]:     "<<fInitMom<<std::endl;

  std::cout<<"  # of final particles:      "<<fFinlPdg.size()<<std::endl;
  std::cout<<"    PDG code:                "<<fFinlPdg[0];
  for (unsigned int i=1; i<fFinlPdg.size(); i++){ std::cout<<", "<<fFinlPdg[i]; } std::cout<<std::endl;

  std::cout<<"  # of spectators:           "<<fSpecPdg.size()<<std::endl;
  if( fSpecPdg.size() > 0 ){
    std::cout<<"    PDG code:                "<<fSpecPdg[0];
    for (unsigned int i=1; i<fSpecPdg.size(); i++){ std::cout<<", "<<fSpecPdg[i]; } std::cout<<std::endl;
  }

  std::cout<<"  degree of Legendre:        "<<fPol.size()<<std::endl;
  if( fPol.size() > 0 ){
    std::cout<<"  coefficients of func:      "<<fPol[0];
    for (unsigned int i=1; i<fPol.size(); i++){ std::cout<<", "<<fPol[i]; } std::cout<<std::endl;
    std::cout<<"  maximun random number:     "<<fPolMax<<std::endl;
  }
}


//-----------------------------------------------//
// class CorssSectionTable
//-----------------------------------------------//
//////////////////////////////////////////////////////
CrossSectionTable::CrossSectionTable(const char* CSFileName, double initMom)
//////////////////////////////////////////////////////
{
  Init();

  //=== read CS from file & fill ===//
  const int MAXCHAR = 512;
  const int MAXTOKEN = 20;
  const char* DELIMITER = " ";
  char buf[MAXCHAR];

  std::ifstream *CSFile = new std::ifstream(CSFileName);
  if ( CSFile->fail() ) {
    std::cerr << CSFileName << " doesn't exist" << std::endl;
    exit(-1);
  } else {
    std::cout << " CS List File : " << CSFileName << " is opened." << std::endl;
  }

  while ( !CSFile->eof() ||  !CSFile->fail() ){
    CSFile->getline(buf, MAXCHAR);
    if( buf[0]=='#' ) continue;
    int n = 0;
    const char* token1[MAXTOKEN] = {};
    token1[0] = strtok(buf, DELIMITER);

    if (token1[0] != NULL) {
      if( !strcmp(token1[0], "EOF") ) break;
      if( !strcmp(token1[0], "#") ) continue;
      for (n = 1; n < 11; n++){
        token1[n] = strtok( NULL, DELIMITER);
        if (token1[n] == NULL ) break;
      }
    }
    else continue;
    int reacID    = atoi(token1[0]);
    double cs     = atof(token1[1]);
    double csf    = atof(token1[2]);
    int beam      = atoi(token1[4]);
    int target    = atoi(token1[5]);
    int react     = atoi(token1[6]);
    int nFinl     = atoi(token1[7]);
    int nSpec     = atoi(token1[8]);
    int nPol      = atoi(token1[9]);
    double polMax = atof(token1[10]);
    std::vector <int> init;
    init.push_back(beam);
    init.push_back(target);
    init.push_back(react);

    CSFile->getline(buf, MAXCHAR);
    const char* token2[MAXTOKEN] = {};
    token2[0] = strtok(buf, DELIMITER);
    if (token2[0] != NULL) {
      if( !strcmp(token2[0], "#") ) continue;
      for (n = 1; n < nFinl+nSpec+nPol; n++){
        token2[n] = strtok( NULL, DELIMITER);
        if (token2[n] == NULL ) break;
      }
    }
    else continue;
    std::vector <int> finl, spec;
    std::vector <double> pol;
    for ( int i=0; i<nFinl; i++ ){
      finl.push_back(atoi(token2[i]));
    } 
    for ( int i=0; i<nSpec; i++ ){
      spec.push_back(atoi(token2[i+nFinl]));
    } 
    for ( int i=0; i<nPol; i++ ){
      pol.push_back(atof(token2[i+nFinl+nSpec]));
    } 

    CrossSection tmp;
    tmp.Init(initMom);
    tmp.SetId(reacID);
    tmp.SetCs(cs);
    tmp.SetCsFactor(csf);
    tmp.SetInitPdg(init);
    tmp.SetFinlPdg(finl);
    tmp.SetSpecPdg(spec);
    tmp.SetPol(pol);
    tmp.SetPolMax(polMax);
    SetCS(tmp);
  } // while()
  //=== read CS from file & fill ===//

  //--- dump each cs ---//
  //PrintAllCS();

  //--- check overlap in CS-list ---//
  if( !CheckCSList() ){
    std::cerr<<"!!! CrossSectionTable::CrossSectionTable(): overlap in CS-list --> exit !!!"<<std::endl;
    exit(0);
  }

  //--- check total-CS ---//
  if( CalcTotalCS() == 0 ){
    std::cerr<<"!!! CrossSectionTable::CrossSectionTable(): total CS = 0 --> exit !!!"<<std::endl;
    exit(0);
  }
}

//////////////////////////////////////////////////////
double CrossSectionTable::CheckCSList()
//////////////////////////////////////////////////////
{
  std::list <int> reacid;
  std::vector<CrossSection>::const_iterator it;
  for(it = fCS.begin(); it != fCS.end(); ++it){
    reacid.push_back(it->Id());
  }
  int size1 = reacid.size();
  reacid.sort();
  reacid.unique();
  int size2 = reacid.size();
  return (size1==size2);
}

//////////////////////////////////////////////////////
double CrossSectionTable::CalcTotalCS()
//////////////////////////////////////////////////////
{
  fTotalCS = 0;
  std::vector<CrossSection>::const_iterator it;
  for(it = fCS.begin(); it != fCS.end(); ++it){
    fTotalCS += it->Cs()*it->CsFactor();
  }
  return fTotalCS;
}

//////////////////////////////////////////////////////
void CrossSectionTable::PrintAllCS()
//////////////////////////////////////////////////////
{
  // #################################
  // # dump all reactions
  // #################################
  std::cout<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<std::endl;
  std::vector <CrossSection>::iterator it;
  for(it=fCS.begin(); it!=fCS.end(); ++it){
    it->PrintCS();
  }
  std::cout<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<std::endl;
  //--- total-CS dump ---//
  CalcTotalCS();
  std::cout<<"======================================================"<<std::endl;
  std::cout<<" ***** input total cross-sections *****"<<std::endl;
  std::cout<<"        "<<fTotalCS<<" [mb]"<<std::endl;
  std::cout<<"======================================================"<<std::endl;

}
