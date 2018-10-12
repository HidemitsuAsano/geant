#ifndef KNUCLCARDSINPUT_H
#define KNUCLCARDSINPUT_H

#include <stdio.h>
#include <stdlib.h>

class ComCardsInput {
 public:
  ComCardsInput();
  ComCardsInput(const char* fname);
  virtual ~ComCardsInput();
  FILE* OpenFile(const char* fname);
  //
  virtual int FindEntry(const char* key, const char* mod="", 
                        const char* form="", int rew=0, int errmsg=1);
  virtual int NextEntry(const char* key, const char* mod="",const char* form="")
;
  virtual int NextEntry();
  virtual int GetNArgs() {return fNArgs;}
  void  Rewind() {if (fCardsFile) rewind(fCardsFile);}
  void  StepBack() {if (fCardsFile) fseek(fCardsFile, fLastPos, SEEK_SET);}
  char* GetKey() {return fKey;}
  char* GetModifier() {return fModifier;}
  char* GetArg(const int iarg, const char* option="", int *err=0);
  char* GetArg(char* &dest, const int iarg, const char* option="", int *err=0);
  float GetArgF(const int iarg, int *err=0);
  int   GetArgD(const int iarg, int *err=0);
  char  **GetArgs() {return fArgs;} 
  int   CompareKey(const char *key);
  int   CompareModifier(const char *mod);
  int   CompareArgList(const char *form);
  virtual void Print();
  //
 protected:
  virtual void ClearArgs(); 
 protected:
  static const char fgkComment='#';   // comment identifier
  static const char fgkDelimiter=':'; // delimiter between keyword and modifier
  static const char fgkContinuation='\\'; // delimiter between keyword and modifier
  static const int  fgkMaxLen = 2048;  // max. length of the entry
  //
  FILE* fCardsFile;        // pointer on the opened file
  long  fLastPos;            // position in the stream of last record read
  char  fbuffer[fgkMaxLen];// string buffer for current line
  char* fPoint;            // pointer on the beginning of data in the line
  int   fNArgs;            // number of the arguments in the entry
  char** fArgs;            // list of the arguments
  char* fKey;              // current Key
  char* fModifier;         // current Modifier
  //
};

#endif

