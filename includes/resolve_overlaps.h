#ifndef RESOLVE_OVERLAPS_H
#define RESOLVE_OVERLAPS_H

#include <iostream>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <string>
#include "bead_residues.h"


namespace resolve_overlaps
{
    void monte_carlo(MolecularModeling::Assembly *glycoprotein, GlycosylationSiteVector *glycosites);
    void example_for_Gordon(MolecularModeling::Assembly *glycoprotein, GlycosylationSiteVector *glycosites);
    void genetic_algorithm(MolecularModeling::Assembly *glycoprotein, GlycosylationSiteVector *glycosites);
}

#endif // RESOLVE_OVERLAPS_H
