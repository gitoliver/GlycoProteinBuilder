#ifndef IO_H
#define IO_H
#include <string>
#include <fstream>      // std::ifstream
#include <stdio.h>  /* defines FILENAME_MAX */
#include <sys/param.h> // for MIN function
#include <stdio.h>  /* defines FILENAME_MAX */

#ifdef WINDOWS
    #include <direct.h>
    #define GetCurrentDir _getcwd
#else
    #include <unistd.h>
    #define GetCurrentDir getcwd
#endif

bool fileExists(std::__cxx11::string fileName);
std::string SplitFilename (const std::string& str);
std::string Find_Program_Installation_Directory();
std::string Find_Program_Working_Directory();

#endif // IO_H
