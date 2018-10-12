// DetectorList.h
#ifndef DetectorList_h
#define DetectorList_h 1

#define K18ANA 0
#include <iostream>
#include <string>
#include <vector>
#include <map>
#if K18ANA
#include <TObject.h>
class DetectorList : public TObject
#else
class DetectorList 
#endif
{
 public:
  static DetectorList* GetInstance()
  {
    static DetectorList instance;
    return &instance;
  }
  void Initialize(const char* FileName=NULL);
  unsigned int GetCID(const std::string &name);
  std::string GetName    (const unsigned int &cid) { return GetString(NAME,cid); }
  std::string GetMaterial(const unsigned int &cid) { return GetString(MATERIAL,cid); }
  std::string GetColor   (const unsigned int &cid) { return GetString(COLOR,cid); }
  bool GetFlag(const unsigned int &cid, int num);
  int GetType(const unsigned int &cid) { return GetData(TYPE,cid); }
  int GetNum1(const unsigned int &cid) { return GetData(NUM1,cid); }
  int GetNum2(const unsigned int &cid) { return GetData(NUM2,cid); }
  int GetNum3(const unsigned int &cid) { return GetData(NUM3,cid); }
  bool IsChamber (const unsigned int &cid) { return IsList(CType_Chamber,cid); }
  bool IsCounter (const unsigned int &cid) { return IsList(CType_Counter,cid); }
  bool IsMTDC    (const unsigned int &cid) { return IsList(CType_MTDC,cid); }
  bool IsDetector(const unsigned int &cid) { return !IsList(CType_NODET,cid); }
  
  int GetNlayers (const unsigned int &cid) { return IsChamber(cid) ? GetNum1(cid) : 0 ; } 
  int GetNwires  (const unsigned int &cid) { return IsChamber(cid) ? GetNum2(cid) : 0 ; } 
  int GetNsegs   (const unsigned int &cid) { return !IsChamber(cid) ? GetNum1(cid) : 0 ; } 
  int GetNsensors(const unsigned int &cid) { return IsCounter(cid) ? GetNum2(cid) : 0 ; } 


  
  unsigned int nList()    { return CIDContainer.size(); }
  unsigned int nChamber() { return CIDList[CType_Chamber].size(); }
  unsigned int nCounter() { return CIDList[CType_Counter].size(); }
  unsigned int nMTDC()    { return CIDList[CType_MTDC].size(); }
  unsigned int nNODET()   { return CIDList[CType_NODET].size(); }
  
  std::map<std::string,unsigned int> const GetList(){ return CIDContainer; }

  //  void Dump();

 private:
  DetectorList(){}
  DetectorList(const DetectorList& rhs);
  DetectorList& operator=(const DetectorList& rhs);
  
  void Clear();
  enum gCounterType{ CType_Chamber=0,
		     CType_Counter=1,
		     CType_MTDC=2,
		     CType_NODET=3
  };
  static const int ntype=4;
  enum gStringType{ NAME=0,
		    MATERIAL=1,
		    COLOR=2
  };

  enum gNumberType{ TYPE=0,
		    NUM1=1,
		    NUM2=2,
		    NUM3=3
  };
  
  std::map<std::string,unsigned int> CIDContainer;
  std::map<unsigned int,std::string> StringContainer[3];
  std::map<unsigned int,int> DataContainer[4];
  std::map<unsigned int,bool> FlagContainer[3];
  std::vector<unsigned int> CIDList[ntype];
  
  bool IsList(const int &type,const unsigned int &cid);
  int GetData(const int &type,const unsigned int &cid);
  std::string GetString(const int &type,const unsigned int &cid);
#if K18ANA  
  ClassDef(DetectorList, 1);
#endif
};

#endif
