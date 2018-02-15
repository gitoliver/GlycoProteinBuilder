#ifndef BEAD_RESIDUES_H
#define BEAD_RESIDUES_H

#include "glycosylationsite.h"

//*******************************************
typedef std::vector<GlycosylationSite*> GlycoSiteVector;
typedef std::vector<Residue*> ResidueVector;
typedef std::vector<Atom*> AtomVector;
//*******************************************

AtomVector Add_Beads(Assembly glycoprotein, GlycoSiteVector glycosites);
void Remove_Beads(Assembly glycoprotein);
double Atomwise_CalculateAtomicOverlaps(Atom *atomA, Atom *atomB, double radiusA, double radiusB);
double modified_CalculateAtomicOverlaps(AtomVector atomsA, AtomVector atomsB);
void fast_overlap_calculation();

#endif // BEAD_RESIDUES_H

