#ifndef BEAD_RESIDUES_H
#define BEAD_RESIDUES_H

#include "glycosylationsite.h"

//*******************************************
typedef std::vector<GlycosylationSite> GlycosylationSiteVector;
typedef std::vector<Residue*> ResidueVector;
typedef std::vector<Atom*> AtomVector;
//*******************************************

void Add_Beads(Assembly *glycoprotein, GlycosylationSiteVector *glycosites);
void Remove_Beads(Assembly glycoprotein);
double Calculate_Bead_Overlap(AtomVector beads);


double Atomwise_CalculateAtomicOverlaps(Atom *atomA, Atom *atomB, double radiusA, double radiusB);
double modified_CalculateAtomicOverlaps(AtomVector atomsA, AtomVector atomsB);

#endif // BEAD_RESIDUES_H

