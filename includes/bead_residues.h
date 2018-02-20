#ifndef BEAD_RESIDUES_H
#define BEAD_RESIDUES_H

#include "glycosylationsite.h"
typedef std::vector<GlycosylationSite> GlycosylationSiteVector;

void Add_Beads(Assembly *glycoprotein, GlycosylationSiteVector *glycosites);
void Remove_Beads(Assembly glycoprotein);

#endif // BEAD_RESIDUES_H

