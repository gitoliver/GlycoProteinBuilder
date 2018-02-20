#include "bead_residues.h"

void Add_Beads(MolecularModeling::Assembly *glycoprotein, GlycosylationSiteVector *glycosites)
{
	AtomVector protein_beads; 
    ResidueVector all_residues = glycoprotein->GetAllResiduesOfAssembly();
    // Go through, find all protein residues, add bead on CA atom.
    for (ResidueVector::iterator it1 = all_residues.begin(); it1 != all_residues.end(); ++it1)
    {
        Residue *residue = *it1;
        // std::cout << (*resi_iter)->GetName() << "\t" << (*resi_iter)->CheckIfProtein() << endl;
        if (residue->CheckIfProtein()==1) // the current residue is an amino acid
        {
            //std::cout << residue->GetName() << "\tA\t" << residue->CheckIfProtein() << endl;
            Atom *atomCA;
            AtomVector atoms = residue->GetAtoms();
            for (AtomVector::iterator it2 = atoms.begin(); it2 != atoms.end(); ++it2)
            {
                Atom *atom = *it2;
                if (atom->GetName().compare("CA")==0)
                {
                	//std::cout << "Adding bead to protein " << residue->GetId() << std::endl;
                    Atom* bead_atom = new Atom(residue, "3fat", atom->GetCoordinates());
                    residue->AddAtom(bead_atom);
                    protein_beads.push_back(bead_atom);
                }
            }
        }
    }
    // Go through all glycosite glycans, add bead in center of each residue, attach it to one other atom in residue.
    // Then set protein_beads
    for (GlycosylationSiteVector::iterator it1 = glycosites->begin(); it1 != glycosites->end(); ++it1)
    {
    	GlycosylationSite *glycosite = &(*it1);
    	glycosite->SetProteinBeads(&protein_beads);
    	AtomVector these_beads;
    	ResidueVector glycan_residues = glycosite->GetAttachedGlycan()->GetResidues();
        for (ResidueVector::iterator it2 = glycan_residues.begin(); it2 != glycan_residues.end(); ++it2)
        {
        	Residue *residue = *it2;
        	if (residue->GetName().compare("SUP") !=0) // don't add one to the superimposition atoms
            {
                // std::cout << (*resi_iter)->GetName() << "\tG\t" << (*resi_iter)->CheckIfProtein() << endl;
               // std::cout << "Adding bead to self glycan " << residue->GetId() << std::endl;
                Atom* bead_atom = new Atom(residue, "4fat", residue->GetGeometricCenter());
                residue->AddAtom(bead_atom);
                these_beads.push_back(bead_atom);
                //Bond bead_atom to any other atom in residue so when glycan is moved, bead_atom moves too.
                Atom *any_atom = residue->GetAtoms().at(0); // 0 is arbitrary, any atom would do.
                //std::cout << "Blow here?" << any_atom->GetId() << std::endl;
                any_atom->GetNode()->AddNodeNeighbor(bead_atom);
                AtomVector temp = {any_atom};
                AtomNode *node = new AtomNode(); // DELETE IS FOR LOSERS.
                bead_atom->SetNode(node);
                bead_atom->GetNode()->SetNodeNeighbors(temp);
            }
        }
        glycosite->SetSelfGlycanBeads(&these_beads);
    }
    // Now find beads from other glycans and add them to list of other_glycan_beads for each glycosite
    for (GlycosylationSiteVector::iterator it1 = glycosites->begin(); it1 != glycosites->end(); ++it1)
    {
    	GlycosylationSite *glycosite1 = &(*it1);
    	AtomVector other_glycan_beads;
    	for (GlycosylationSiteVector::iterator it2 = glycosites->begin(); it2 != glycosites->end(); ++it2)
    	{
    		GlycosylationSite *glycosite2 = &(*it2);
    		if(glycosite1 != glycosite2) // Check if same site
    		{
    			// append each other glycosite's beads to list of other_glycan_beads: a.insert(std::end(a), std::begin(b), std::end(b));
    			AtomVector temp = glycosite2->GetSelfGlycanBeads();
    			other_glycan_beads.insert(std::end(other_glycan_beads), std::begin(temp), std::end(temp));
    			//std::cout << "Adding beads of glycosite " << glycosite2->GetResidue()->GetId() << " to " << glycosite1->GetResidue()->GetId() << std::endl;

    		}

    	}
    	glycosite1->SetOtherGlycanBeads(&other_glycan_beads);
    }
}

void Remove_Beads(MolecularModeling::Assembly glycoprotein)
{
    ResidueVector all_residues = glycoprotein.GetAllResiduesOfAssembly();
    for (ResidueVector::iterator it1 = all_residues.begin(); it1 != all_residues.end(); ++it1)
    {
        Residue *residue = *it1;
        AtomVector atoms = residue->GetAtoms();
        for (AtomVector::iterator it2 = atoms.begin(); it2 != atoms.end(); ++it2)
        {
            Atom *atom = *it2;
            if (atom->GetName().find("fat")==1)
            {
                residue->RemoveAtom(atom);
            }
        }
    }
}

