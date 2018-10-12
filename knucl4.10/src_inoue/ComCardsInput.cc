#include "ComCardsInput.hh"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>


ComCardsInput::ComCardsInput() 
{
  // Initialize only
  *fbuffer = '\0';
  fPoint = 0;
  fNArgs = 0;
  fArgs = 0;
  fKey = 0;
  fModifier = 0;
  fCardsFile = 0;
}
//--------------------------------------------------
//
ComCardsInput::ComCardsInput(const char* fname) 
{
  // Initialize only
  *fbuffer = '\0';
  fPoint = 0;
  fNArgs = 0;
  fArgs = 0;
  fKey = 0;
  fModifier = 0;
  fCardsFile = 0;
  fCardsFile = OpenFile(fname);
}
//--------------------------------------------------
//
ComCardsInput::~ComCardsInput() 
{
  if (fCardsFile) fclose(fCardsFile);
  ClearArgs();
}
//--------------------------------------------------
//
FILE* ComCardsInput::OpenFile(const char *fname) 
{
  if (fCardsFile) fclose(fCardsFile);
  if ( !(fCardsFile = fopen(fname,"r")) ) { 
    printf("Error in <ComCardsInput::OpenFile>: Did not find file %s\n",fname);
    exit(1);
  }  
  ClearArgs();
  return fCardsFile;
}
//--------------------------------------------------
//
void ComCardsInput::ClearArgs() 
{
  if (fArgs) {
    for (int i=0;i<fNArgs;i++) delete[] fArgs[i];
    delete[] fArgs;
    fArgs = 0;
  }
  if (fKey) {
    delete[] fKey;
    fKey = 0;
  }
  if (fModifier) {
    delete[] fModifier;
    fModifier = 0;
  }
  fNArgs = 0;
}
//--------------------------------------------------
//
int ComCardsInput::FindEntry(const char* ckey, const char* cmod,const char* form,
 int rew, int errmsg) 
{
  // Find and fill arguments list from the line starting with KEY:MODIFIER ...
  // KEY(ckey) is obligatory, while MODIFIER(cmod) can be either 0 (any
  // modifier matches), or "" or " .." - empty modifier is requested.
  // Return number of found arguments (not counting key:modifier)
  // If form is not NULL or nonempty, then looks for the record in the following
  // format according to characters in the "form" string: 
  // 's' - argument should be string (any type of argument mathces)
  // 'a' - argument should be string with at least one non-numeric char
  // 'd' - argument should be integer
  // 'f' - argument should be float (no check for . is done so the integer too matches)
  // '?' - following arguments are optional, but their number should be not less than in 'form' 
  // '*' - following arguments are optional, and their number can be less than in 'form'
  // '|' - no more arguments is allowed
  // If rew != 0, the search is done from the beginning of the file
  // If errmsg is not 0, then the warning is printed in case KEY and MODIFIER were matched
  // but format - not
  // Case insensitive
  // If not found, return -1;
  //
  int narg=-1;
  //
  if ( !ckey ) {
    printf("<ComCardsInput::FindEntry>  KEY is empty\n");
    return -1;
  }
  if (rew) rewind(fCardsFile);
  //
  while ( (narg=NextEntry()) >-1 ) {
    if ( CompareKey(ckey) ) continue; // KEY does not match
    if ( cmod && (cmod[0]=='\0' || cmod[0]==' ') 
         && fModifier[0]!='\0' ) continue; // Empty modifier was requested
    if ( cmod && CompareModifier(cmod) ) continue; // MODIFIER does not match
    //
    if ( CompareArgList(form) ) {  // ArgList form does not match
      if (cmod && errmsg) {
        printf("<ComCardsInput::FindEntry> arguments of %s\n %s to \"%s\"\n",
               fPoint,"do not match to format",form);
        //return -10;
      }
      continue;
    }
    return narg;
  }
  return -1;
}
//--------------------------------------------------
//
int ComCardsInput::NextEntry(const char* ckey, const char* cmod,const char* form)
{
  // Reads next entry in required format
  // KEY(ckey) is obligatory, while MODIFIER(cmod) can be either 0 (any
  // modifier matches), or "" or " .." - empty modifier is requested.
  // Return number of found arguments (not counting key:modifier)
  // If form is not NULL or nonempty, then looks for the record in the following
  // format according to characters in the "form" string: 
  // 's' - argument should be string (any type of argument mathces)
  // 'a' - argument should be string with at least one non-numeric char
  // 'd' - argument should be integer
  // 'f' - argument should be float (no check for . is done so the integer too matches)
  // '?' - following arguments are optional, but their number should be not less than in 'form' 
  // '*' - following arguments are optional, and their number can be less than in 'form'
  // '|' - no more arguments is allowed
  // Case insensitive
  // If next entry does not match, return -1;
  //
  int narg=-1;
  //
  if ( !ckey ) {
    printf("<ComCardsInput::NextEntry(...)>  KEY is empty\n");
    return -1;
  }
  //
  narg=NextEntry();
  if (narg<0) return narg;
  if ( CompareKey(ckey) ) return -1; // KEY does not match
  if (!(cmod && (cmod[0]=='\0'||cmod[0]==' ')&&fModifier[0]!='\0')) // NonEmpty modifier requested
    if ( cmod && CompareModifier(cmod) ) return -1; // MODIFIER does not match
    //
  if ( CompareArgList(form) ) {  // ArgList form does not match
    printf("<ComCardsInput::NextEntry> arguments of %s\n %s to \"%s\"\n",
           fPoint,"do not match to format",form);
    return -1;
  }
  return narg;
  //
}

int ComCardsInput::NextEntry() 
{
  // get next entry, skipping commented lines
  //
  int len;
  char *tmp1,*tmp2;
  char buftmp[fgkMaxLen]; 
  //
  if (!fCardsFile) {
    printf("<ComCardsInput::NextEntry> - No file was opened\n");
    return -1;
  }
  ClearArgs(); 
  while ( 1 ) {
    fLastPos = ftell(fCardsFile);
    if ( !(fPoint=fgets(fbuffer,fgkMaxLen,fCardsFile)) ) break; // end of file reached
    while(*fPoint == ' ') fPoint++; // skip spaces in the beginning
    if (*fPoint == fgkComment) continue; // check if the line is commented
    //
    // Check if there is a line continuation flag
    tmp2 = fPoint;
    do {
      tmp1 = &fPoint[strlen(fPoint)-1];
      if (*tmp1=='\n') {*tmp1--='\0';}
      while(*tmp1 == ' ') *tmp1--='\0'; // skip spaces in the end
      if (*tmp1 == fgkContinuation) {
        if ( !(tmp2=fgets(tmp1,fgkMaxLen-(tmp1-fPoint),fCardsFile)) ) break; // end of file reached
      }
      else
        break;
    } while(1);
    if (!tmp2) break; // end of file reached
    //
    // Check if there is KEY:[MODIFIER]
    tmp1 = fPoint;
    if ( sscanf(tmp1,"%s",buftmp)<=0 ) continue; // Empty line
    if ( (tmp2=strchr(buftmp,fgkDelimiter)) ) { // there is a key/delimiter
      len = (tmp2-buftmp)/sizeof(char);
      fKey = new char[len+1];
      for (int i=0;i<len;i++) fKey[i] = tolower(buftmp[i]);
      fKey[len]='\0';
      //
      len = strlen(++tmp2); // skip delimiter, get modifier length
      fModifier = new char[len+1];
      for (int i=0;i<len;i++) fModifier[i] = tolower(tmp2[i]);
      fModifier[len] = '\0';
      // skip key:mod
      tmp1 += strlen(buftmp);
      while(*tmp1 == ' ') tmp1++; // skip spaces
    }
    else {
      continue; // no delimiter - not valid record
    }

    // First, throw away everything after comment (if any)
    if ( (tmp2=strchr(tmp1,fgkComment)) ) *tmp2 = '\0';
    // now, scan the string to get the number of arguments
    tmp2 = tmp1;
    while( sscanf(tmp2,"%s",buftmp)!= -1 ) {
      while( *tmp2 == ' ') tmp2++;
      tmp2 += strlen(buftmp); 
      fNArgs++;
    }
    // Fill the arguments list
    fArgs = new char*[fNArgs];
    for (int i=0;i<fNArgs;i++) {
      sscanf(tmp1,"%s",buftmp);
      while( *tmp1 == ' ') tmp1++;
      int lgt = strlen(buftmp);
      tmp1 += lgt;
      fArgs[i] = new char[lgt+1];
      strcpy(fArgs[i],buftmp);
    }
    return fNArgs;
  }
  return -1;
}

//--------------------------------------------------
//
char* ComCardsInput::GetArg(const int iarg, const char* option, int *err) 
{
  // Return argument iarg in character string format
  // Optionally converting it to upper or lower case
  // if iarg is wrong, return 0 and set err to 1
  //
  err = 0;
  if (iarg>=fNArgs || iarg<0) { *err = 1; return 0; }
  char *carg = fArgs[iarg];
  int lgt = strlen(carg);
  if (!strncasecmp(option,"l",1))  // convert to lower case
    for (int i=0;i<lgt;i++) 
      carg[i] = tolower(carg[i]);
  else if (!strncasecmp(option,"u",1)) // convert to upper case
    for (int i=0;i<lgt;i++) 
      carg[i] = toupper(carg[i]);
  //
  return carg;
}
//--------------------------------------------------
//
char* ComCardsInput::GetArg(char* &dest, const int iarg, const char* option, int *err) 
{
  // Creates new string at dest and fill it with argument iarg 
  // in character string format.
  // Optionally converting it to upper or lower case
  // if iarg is wrong, return 0 and set err to 1
  // No check of dest == 0 is done!
  err = 0;
  if (iarg>=fNArgs || iarg<0) { *err = 1; return 0; }
  char *carg = fArgs[iarg];
  int lgt = strlen(carg);
  dest = new char[lgt+1];
  if (!strncasecmp(option,"l",1))  // convert to lower case
    for (int i=0;i<lgt;i++) 
      dest[i] = tolower(carg[i]);
  else if (!strncasecmp(option,"u",1)) // convert to upper case
    for (int i=0;i<lgt;i++) 
      dest[i] = toupper(carg[i]);
  else
    for (int i=0;i<lgt;i++) dest[i] = carg[i];
  //
  dest[lgt] = '\0';
  return dest;
}
//--------------------------------------------------
//
//--------------------------------------------------
//
float ComCardsInput::GetArgF(const int iarg, int *err) 
{
  // Get requsted argument in float format, it it is not float, set Error to 1
  err = 0;
  if (iarg>=fNArgs || iarg<0) { *err = 1; return 0.; }
  //
  char *errc = 0;
  float vald = float(strtod(fArgs[iarg],&errc));
  if (*errc) { *err = 1; return 0.; }
  return vald;
  //
}
//--------------------------------------------------
//
int ComCardsInput::GetArgD(const int iarg, int *err) 
{
  // Get requsted argument in integer format, it it is not float, set Error to 1
  err = 0;
  if (iarg>=fNArgs || iarg<0) { *err = 1; return 0; }
  //
  char *errc = 0;
  int vald = int(strtol(fArgs[iarg],&errc,10));
  if (*errc) { *err = 1; return 0; }
  return vald;
  //
}
//--------------------------------------------------
//
int ComCardsInput::CompareKey(const char* key) 
{
  if (fKey) return strcasecmp(key,fKey);
  return 0;
}
//--------------------------------------------------
//
int ComCardsInput::CompareModifier(const char* mod) 
{
  if (fModifier) return strcasecmp(mod,fModifier);
  return 1;
}
//--------------------------------------------------
//

//
int ComCardsInput::CompareArgList(const char *form)
{
  // If form is not NULL or nonempty, check if the record matches following 
  // format according to characters in the "form" string: 
  // 's' - argument should be string (any type of argument mathces)
  // 'a' - argument should be string with at least one non-numeric char
  // 'd' - argument should be integer
  // 'f' - argument should be float (no check for . is done so the integer too matches)
  // '?' - following arguments are optional, but their number should be not less than in 'form' 
  // '*' - following arguments are optional, and their number can be less than in 'form'
  // '|' - no more arguments is allowed
  // If matches, return 0, otherwise -1
  //
  char cf=' ';
  const char *pf = form;
  char *err = 0;
  if ( !form ) return 0; // no for is requested, anything matches
  int iarg = 0;
  float valf=0.0;
  int vald = 0;
  int isoption = 0;
  //
  if (fNArgs == -1 ) {
    printf("<ComCardsInput::CompareArgList> There is no record in the buffer\n");
    return -1;
  }
  while ( (cf=tolower(*pf++)) ) {
    //
    if (cf=='*') { // the following arguments are optional
      isoption = 1;
      continue;
    }
    //
    if (cf=='?') { // the following argument block is optional
      if ( fNArgs > iarg+1 ) {
        isoption = 0; // number of following arguments should be respected
        continue;
      }
      else return 0; // no more arguments, but the block was optional
    }
    //
    if (cf=='|') { // no more arguments are alowed
      if ( fNArgs != iarg ) return-1; // number of arguments > than allowed
      else return 0; // no more arguments, OK
    }
    //
    if ( iarg >= fNArgs ) { // number of agruments is less than requested
      if ( isoption || cf=='|' ) return 0; // but this was allowed -> MATCHED
      else return -1; // too few arguments
    }
    //
    // Analize requsted argument type
    //
    //
    if (cf=='s') { // string argument is requested, anything matches
      iarg++;
      continue;
    }
    //
    if (cf=='d') { // Integer argument is requested
      vald = strtol(fArgs[iarg++],&err,10);
      if (*err) return -1; // Not integer, does not match
      else continue; // this argument matches
    }
    //
    if (cf=='f') { // Float argument is requested
      valf = strtod(fArgs[iarg++],&err);
      if (*err) return -1; // Not float, does not match
      else continue; // this argument matches
    }
    //
    if (cf=='a') { // string argument with non-numeric character is requested
      valf = strtod(fArgs[iarg++],&err);
      if (!(*err)) return -1; // looks like number, does not match
      else continue; // this argument matches
    }
    //
    // unknown format character
    printf("<ComCardsInput::CompareArgList> unknown format requsted \'%c\'\n",cf);
  }
  return 0; // Matched
}
//--------------------------------------------------
//
void ComCardsInput::Print() 
{
  if (fKey) printf("Key = %s\n",fKey);
  else printf("No Key\n");
  if (fModifier) printf("Modifier = %s\n",fModifier);
  else printf("No Modifiers\n");
  for (int i=0;i<fNArgs;i++) {
    printf("Arg %d: %s\n",i+1,fArgs[i]);
  }
}



