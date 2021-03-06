#ifndef GLYCOPROTEIN_BUILDER_H
#define GLYCOPROTEIN_BUILDER_H

#include <string>
#include <dirent.h>
//#include <unistd.h>
#include <sys/stat.h>
//#include <sys/types.h>
//#include <stdlib.h>     /* getenv */
//#include <fstream>      // std::ifstream
#include "glycosylationsite.h"
#include "bead_residues.h"
#include "io.h"

typedef std::vector<GlycosylationSite> GlycosylationSiteVector;
typedef std::vector<GlycosylationSite*> GlycosylationSitePointerVector;

// This should become it's own class
namespace glycoprotein_builder
{
/*******************************************/
/* Function Declarations                   */
/*******************************************/
void Read_Input_File(GlycosylationSiteVector &glycosites, std::string &proteinPDB, std::string &glycanDirectory, const std::string working_Directory);
void AttachGlycansToGlycosites(MolecularModeling::Assembly &glycoprotein, GlycosylationSiteVector &glycosites, const std::string glycanDirectory);
void SetOtherGlycosites(GlycosylationSiteVector &glycosites);
void PrintDihedralAnglesAndOverlapOfGlycosites(GlycosylationSiteVector &glycosites);
void SetDefaultDihedralAnglesUsingMetadata(GlycosylationSiteVector &glycosites);
void SetRandomDihedralAnglesUsingMetadata(GlycosylationSiteVector &glycosites);
//double CalculateAtomicOverlaps(GlycosylationSiteVector &glycosites);
//double RandomAngle_360range();
//double RandomAngle_range(int min, int max);
//double RandomAngle_PlusMinusX(double start_point, int max_step_size);
//double GetNewAngleScaledToDegreeOfOverlap(double current_angle, double overlap, int number_of_atoms);
void write_pdb_file(MolecularModeling::Assembly *glycoprotein, int cycle, std::string summary_filename, double overlap);
double GetGlobalOverlap(GlycosylationSiteVector &glycosites);

double CalculateOverlaps(GlycosylationSiteVector &glycosites, OverlapType overlapType = BEAD, MoleculeType moleculeType = ALL, bool recordOverlap = true, bool printOverlap = false);

void PrintOverlaps(GlycosylationSiteVector &glycosites);
void PrintOverlaps(GlycosylationSitePointerVector &glycosites);
void CalculateAndPrintOverlaps(GlycosylationSiteVector &glycosites);
GlycosylationSitePointerVector DetermineSitesWithOverlap(GlycosylationSiteVector &glycosites, double tolerance, OverlapType overlapType = BEAD);
GlycosylationSitePointerVector GetSitesWithOverlap(GlycosylationSiteVector &glycosites, double tolerance);
//GlycosylationSitePointerVector DeleteSitesWithOverlaps(GlycosylationSiteVector &glycosites, double tolerance, std::string overlap_type = "total");
//void DeleteSitesWithOverlapRecordsAboveTolerance(GlycosylationSiteVector &glycosites, double tolerance, std::string overlap_type = "total");
void DeleteSitesIterativelyWithAtomicOverlapAboveTolerance(GlycosylationSiteVector &glycosites, double tolerance);
void UpdateAtomsThatMoveInLinkages(GlycosylationSiteVector &glycosites);
ResidueLinkageVector GetAllFirstAnd1_6Linkages(GlycosylationSiteVector &glycosites);
ResidueLinkageVector GetAllFirstAnd2_XLinkages(GlycosylationSiteVector &glycosites);
void StashCoordinates(GlycosylationSiteVector &glycosites);
void SetStashedCoordinatesWithLowestOverlap(GlycosylationSiteVector &glycosites);
ResidueLinkageVector SplitLinkagesIntoPermutants(ResidueLinkageVector inputLinkages);

//void Overlap_Weighted_Adjust_Torsions_For_X_Cycles(GlycosylationSitePointerVector &sites, GlycosylationSiteVector &glycosites, int max_cycles, double tolerance, std::string overlap_type);
//void Overlap_Weighted_Adjust_Torsions(GlycosylationSitePointerVector &sites);
}

#endif // GLYCOPROTEIN_BUILDER_H
