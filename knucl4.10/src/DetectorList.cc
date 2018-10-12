// DetectorList.cpp
#include "DetectorList.hh"
#include <cstdio>
#include <stdlib.h>
#include <cstdlib>
#include <algorithm>

#if K18ANA
ClassImp(DetectorList)
#endif
const int MAXCHAR = 128;
void DetectorList::Initialize(const char* FileName){
  Clear();
  const char *fn;  
  if (sizeof(FileName) == 0) {
    fn = "Detector.list";
  } else {
    fn = FileName;
  }
  FILE* file = fopen(fn,"r");
  if (file == 0) {
    std::cout<<"List Could not open detector list file "<<fn<<std::endl;
    exit(-1);
  }
  
  char str[MAXCHAR];
  char tmpname[20],matname[20],color[20];
  int cid,type,num1,num2,num3,flag1,flag2;  
  int nd;
  while( fgets(str,MAXCHAR,file)!=0 ){
    if( str[0]=='#' ) continue;
    if( (nd=sscanf(str,"%d %s %d %d %d %d %d %d %s %s",
		   &cid, tmpname, &type, &num1, &num2, &num3, &flag1,&flag2, matname, color)) == 10 ) {
      CIDContainer[tmpname]=cid;
      StringContainer[NAME][cid]=tmpname;
      StringContainer[MATERIAL][cid]=matname;
      StringContainer[COLOR][cid]=color;
      FlagContainer[0][cid]=flag1;
      FlagContainer[1][cid]=flag2;
      DataContainer[TYPE][cid]=type;
      DataContainer[NUM1][cid]=num1;
      DataContainer[NUM2][cid]=num2;
      DataContainer[NUM3][cid]=num3;
      if(type>=0&&type<ntype)
	CIDList[type].push_back(cid);
    }
  }
  fclose(file);
  return;
}

void DetectorList::Clear()
{
  CIDContainer.clear();
  for(int i=0;i<3;i++) StringContainer[i].clear();
  for(int i=0;i<3;i++) DataContainer[i].clear();
  for(int i=0;i<3;i++) FlagContainer[i].clear();
  for(int i=0;i<ntype;i++)   CIDList[i].clear();
}

bool DetectorList::IsList(const int &type,const unsigned int &cid){ 
  if(type>=0&&type<ntype){
    if( find( CIDList[type].begin(), 
	      CIDList[type].end(), 
	      cid ) == CIDList[type].end() ) 
      return false;
    else return true;
  }
  else 
    return false;
}

int DetectorList::GetData(const int &type,const unsigned int &cid){
  if(type<0||type>=ntype) return -1;
  std::map<unsigned int,int>::iterator it = DataContainer[type].find(cid);
  if(it != DataContainer[type].end()) return it->second;
  else return -1;
}

bool DetectorList::GetFlag(const unsigned int &cid,int type){
  if(type<0&&type>1) return false;
  std::map<unsigned int,bool>::iterator it = FlagContainer[type].find(cid);
  if(it != FlagContainer[type].end()) return it->second;
  else return false;
}

std::string DetectorList::GetString(const int &type,const unsigned int &cid){
  std::map<unsigned int,std::string>::iterator it = StringContainer[type].find(cid);
  if(it != StringContainer[type].end()) return it->second;
  else return NULL;
}
unsigned DetectorList::GetCID(const std::string  &name){
  std::map<std::string,unsigned int>::iterator it = CIDContainer.find(name);
  if(it != CIDContainer.end()) return it->second;
  else return -1;
}
