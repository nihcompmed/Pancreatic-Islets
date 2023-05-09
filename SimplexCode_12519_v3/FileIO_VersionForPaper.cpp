#include "./FileIO_VersionForPaper.h"

bool file_exists (const std::string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0); 
}

//*********************************************************
//*********************************************************

string getDirectory(const string & s){ //includes last slash
  string directory;
  const size_t last_slash_idx = s.rfind('/');
  if (std::string::npos != last_slash_idx)
    {
      directory = s.substr(0, last_slash_idx+1);
    }
  return directory;
}

//*********************************************************
//*********************************************************

string getFileheader(const string& s) {//strips extension
  
  size_t i = s.rfind('.', s.length());
  if (i != string::npos) {
    return(s.substr(0, i));
  }
  
  return("");
}

//*********************************************************
//*********************************************************

string getFileName(const string& s) {//strips directory info
  char sep = '/';
  size_t i = s.rfind(sep, s.length());
  if (i != string::npos) {
    return(s.substr(i+1, s.length() - i));
  }
  return(s);
}

//*********************************************************
//*********************************************************
string getGroup(const string& s){ //gets up to first period; if no period returns s
  size_t i = s.find('.');
  if (i != string::npos) {
    return(s.substr(0,i ));
  }
  
  return(s);
}

//*********************************************************
//*********************************************************

// string getFileExt(const string& s) {
  
//   size_t i = s.rfind('.', s.length());
//   if (i != string::npos) {
//     return(s.substr(i+1, s.length() - i));
//   }
  
//   return("");
// }




