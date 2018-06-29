#ifndef RESOLVE_OVERLAPS_H
#define RESOLVE_OVERLAPS_H

#include <iostream>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm> // For erase remove
#include <iomanip>
#include <cmath>
#include <ctime>
#include <cstring>

#include "bead_residues.h"
#include "glycosylationsite.h"
#include "glycoprotein_builder.h"
#include "metropolis_criterion.h"

typedef std::vector<GlycosylationSite*> GlycosylationSitePointerVector;

double RandomAngle_360range();
double RandomAngle_range(int min, int max);
double RandomAngle_PlusMinusX(double start_point, int max_step_size);
double GetNewAngleScaledToDegreeOfOverlap(double current_angle, double overlap, int number_of_atoms);
void write_pdb_file(Assembly *glycoprotein, int cycle, std::string summary_filename, double score);
void PrintOverlaps(GlycosylationSiteVector &glycosites);
void PrintOverlaps(GlycosylationSitePointerVector &glycosites);
void CalculateAndPrintOverlaps(GlycosylationSiteVector &glycosites);
void SetBestChi1Chi2(GlycosylationSitePointerVector &glycosites, std::string overlap_type = "total");
GlycosylationSitePointerVector DetermineSitesWithOverlap(GlycosylationSiteVector &glycosites, double tolerance, std::string overlap_type = "total");
GlycosylationSitePointerVector GetSitesWithOverlap(GlycosylationSiteVector &glycosites, double tolerance);
GlycosylationSitePointerVector DeleteSitesWithOverlaps(GlycosylationSiteVector &glycosites, double tolerance, std::string overlap_type = "total");
void DeleteSitesWithOverlapRecordsAboveTolerance(GlycosylationSiteVector &glycosites, double tolerance, std::string overlap_type = "total");
void DeleteSitesIterativelyWithOverlapAboveTolerance(GlycosylationSiteVector &glycosites, double tolerance);
void Overlap_Weighted_Adjust_Torsions_For_X_Cycles(GlycosylationSitePointerVector &sites, GlycosylationSiteVector &glycosites, int max_cycles, double tolerance, std::string overlap_type);
void Overlap_Weighted_Adjust_Torsions(GlycosylationSitePointerVector &sites);

namespace resolve_overlaps
{
    void Resolve_Overlaps_Random_Walk_Scaled_To_Overlap(GlycosylationSiteVector &glycosites, std::string type = "protein", int max_cycles = 500, double strict_tolerance = 0.1, double loose_tolerance = 1.0);
    void weighted_protein_global_overlap_monte_carlo(GlycosylationSiteVector &glycosites);
    void protein_first_random_walk_scaled_to_overlap(GlycosylationSiteVector &glycosites);
    bool dumb_random_walk(GlycosylationSiteVector &glycosites);
}

#endif // RESOLVE_OVERLAPS_H
