#ifndef RESOLVE_OVERLAPS_H
#define RESOLVE_OVERLAPS_H

#include "/home/oliver/Programs/gems/gmml/includes/gmml.hpp"
#include "glycosylationsite.h"

//*******************************************

typedef std::vector<GeometryTopology::Coordinate*> CoordinateVector;
typedef std::vector<MolecularModeling::Atom*> AtomVector;
typedef std::vector<Residue*> ResidueVector;
typedef std::vector<MolecularModeling::Assembly*> AssemblyVector;
typedef std::vector<GlycosylationSite*> GlycoSiteVector;

//*******************************************

using namespace MolecularModeling;

namespace resolve_overlaps
{
    void monte_carlo(Assembly glycoprotein, GlycoSiteVector glycosites);
}

#endif // RESOLVE_OVERLAPS_H
