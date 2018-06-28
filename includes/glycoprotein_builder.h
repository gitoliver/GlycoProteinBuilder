#ifndef GLYCOPROTEIN_BUILDER_H
#define GLYCOPROTEIN_BUILDER_H

#include <string>
#include <dirent.h>
//#include <unistd.h>
#include <sys/stat.h>
//#include <sys/types.h>
//#include <stdlib.h>     /* getenv */
//#include <fstream>      // std::ifstream
#include "io.h"
#include "glycosylationsite.h"

namespace glycoprotein_builder
{

typedef std::vector<GlycosylationSite> GlycosylationSiteVector;

/*******************************************/
/* Function Declarations                   */
/*******************************************/
void Read_Input_File(GlycosylationSiteVector &glycosites, std::string &proteinPDB, std::string &glycanDirectory, const std::string working_Directory);
void AttachGlycansToGlycosites(MolecularModeling::Assembly &glycoprotein, GlycosylationSiteVector &glycosites, const std::string glycanDirectory);
void PrintDihedralAnglesAndOverlapOfGlycosites(GlycosylationSiteVector &glycosites);
void SetReasonableChi1Chi2Values(GlycosylationSiteVector &glycosites);
void CalculateOverlaps(GlycosylationSiteVector &glycosites);
double GetGlobalOverlap(GlycosylationSiteVector &glycosites);
}

#endif // GLYCOPROTEIN_BUILDER_H
