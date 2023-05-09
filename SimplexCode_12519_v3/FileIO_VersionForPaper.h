#include <iostream>                  // for std::cout
#include <fstream>
#include <utility>                   // for std::pair
#include <algorithm>                 // for std::for_each
#include <vector>
#include <time.h>
#include <cstdlib>
#include <complex>
#include <cmath>
 #include <sys/stat.h>
#include <unistd.h>
 
using namespace std;

  bool file_exists (const std::string& ); 
 string getDirectory(const string & ); 
/* string getFileExt(const string& ); */
 string getFileheader(const string&); 
 string getFileName(const string& ); 
 string getGroup(const string& ); 
