#ifndef BEAD_RESIDUES_H
#define BEAD_RESIDUES_H

#include "glycosylationsite.h"

using namespace MolecularModeling;
typedef std::vector<GlycosylationSite> GlycosylationSiteVector;

void Add_Beads(MolecularModeling::Assembly *glycoprotein, GlycosylationSiteVector *glycosites);
void Remove_Beads(MolecularModeling::Assembly glycoprotein);

double GetMaxDistanceBetweenAtoms(AtomVector atoms);
AtomVector SelectAtomsWithinDistanceOf(MolecularModeling::Atom *query_atom, double distance, AtomVector atoms);
#endif // BEAD_RESIDUES_H

