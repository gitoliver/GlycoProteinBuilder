#ifndef RESOLVE_OVERLAPS_H
#define RESOLVE_OVERLAPS_H

#include <iostream>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <string>
#include <algorithm> // For erase remove
#include <iomanip>
#include <cmath>
#include <ctime>
#include <cstring>

#include "bead_residues.h"
#include "glycosylationsite.h"

typedef std::vector<GlycosylationSite*> GlycosylationSitePointerVector;

double RandomAngle_360range();
double RandomAngle_PlusMinusX(double start_point, int max_step_size);
double GetNewAngleScaledToDegreeOfOverlap(double current_angle, double overlap, int number_of_atoms);
void write_pdb_file(Assembly &glycoprotein, int cycle, std::string summary_filename, double score);
void PrintOverlaps(GlycosylationSiteVector &glycosites);
void PrintOverlaps(GlycosylationSitePointerVector &glycosites);
void SetBestChi1Chi2(GlycosylationSitePointerVector &glycosites, std::string overlap_type = "total");
GlycosylationSitePointerVector DetermineSitesWithOverlap(GlycosylationSiteVector &glycosites, double tolerance, std::string overlap_type = "total");
GlycosylationSitePointerVector DeleteSitesWithOverlaps(GlycosylationSiteVector &glycosites, double tolerance, std::string overlap_type = "total");
void DeleteSitesWithOverlapRecordsAboveTolerance(GlycosylationSiteVector &glycosites, double tolerance, std::string overlap_type = "total");
void Monte_Carlo_Torsions(GlycosylationSitePointerVector &sites, GlycosylationSiteVector &glycosites, int max_cycles, double tolerance, std::string overlap_type);
void RandomizeTorsions(GlycosylationSitePointerVector &sites);

namespace resolve_overlaps
{
    void protein_first_monte_carlo(GlycosylationSiteVector &glycosites);
    void dumb_monte_carlo(GlycosylationSiteVector &glycosites);
}

#endif // RESOLVE_OVERLAPS_H
