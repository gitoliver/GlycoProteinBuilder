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
void write_pdb_file(Assembly *glycoprotein, int cycle, std::string summary_filename, double score);
void PrintOverlaps(GlycosylationSiteVector *glycosites);
void PrintOverlaps(GlycosylationSitePointerVector &glycosites);
void SetBestProteinChi1Chi2(GlycosylationSitePointerVector &glycosites, Assembly *glycoprotein);
GlycosylationSitePointerVector DetermineSitesWithOverlap(GlycosylationSiteVector *glycosites, double tolerance, std::string returning = "with", std::string type = "total");
GlycosylationSitePointerVector DeleteSitesWithOverlaps(GlycosylationSiteVector &glycosites, double tolerance, std::string type);

namespace resolve_overlaps
{
    void monte_carlo(MolecularModeling::Assembly *glycoprotein, GlycosylationSiteVector *glycosites);
    void dumb_monte_carlo(MolecularModeling::Assembly *glycoprotein, GlycosylationSiteVector *glycosites);
}

#endif // RESOLVE_OVERLAPS_H
