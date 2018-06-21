#include "../includes/bead_residues.h"

void Add_Beads(MolecularModeling::Assembly &glycoprotein, GlycosylationSiteVector &glycosites)
{
	AtomVector protein_beads; 
    ResidueVector all_residues = glycoprotein.GetAllResiduesOfAssembly();
    // Go through, find all protein residues, add bead on CA atom and certain sidechain atoms.
    for (ResidueVector::iterator it1 = all_residues.begin(); it1 != all_residues.end(); ++it1)
    {
        Residue *residue = *it1;
        if (residue->CheckIfProtein()==1) // the current residue is an amino acid
        {
            AtomVector atoms = residue->GetAtoms();
            for (AtomVector::iterator it2 = atoms.begin(); it2 != atoms.end(); ++it2)
            {
                Atom *atom = *it2;
                // Okay so all CA atoms and then I've picked out others by hand in VMD that will cover the sidechains
                if (atom->GetName().compare("CA")==0)
                {
                	//std::cout << "Adding bead to protein " << residue->GetId() << std::endl;
                    Atom* bead_atom = new Atom(residue, "mfat", atom->GetCoordinates());
                    residue->AddAtom(bead_atom);
                    protein_beads.push_back(bead_atom);
                }
                else if ( (atom->GetName().compare("NZ")==0) ||
                          (atom->GetName().compare("CZ")==0) ||
                          (atom->GetName().compare("NE2")==0) ||
                          (atom->GetName().compare("OD1")==0) ||
                          (atom->GetName().compare("SD")==0)  ||
                          ( (atom->GetName().compare("CE2")==0) && residue->GetName().compare("TRP")==0 ) ||
                          ( (atom->GetName().compare("CD1")==0) && ( residue->GetName().compare("LEU")==0 || residue->GetName().compare("ILE")==0 ) ) ||
                          ( (atom->GetName().compare("CD")==0) && residue->GetName().compare("GLU")==0 )
                        )
                { // sfats should move when a chi1, chi2 is moved, so make sure they are connected to something for the SetDihedral function to move them.
                    Atom* bead_atom = new Atom(residue, "sfat", atom->GetCoordinates());
                    residue->AddAtom(bead_atom);
                    protein_beads.push_back(bead_atom);
                }
            }
        }
    }
    // Go through all glycosite glycans, add bead in center of each residue, attach it to one other atom in residue.
    // Then set protein_beads
    for (GlycosylationSiteVector::iterator it1 = glycosites.begin(); it1 != glycosites.end(); ++it1)
    {
    	GlycosylationSite *glycosite = &(*it1);
        Atom *cb_atom = glycosite->GetResidue()->GetAtom("CB");
        // This needs to be updated/scrapped if we start moving the glycan
        double distance = ( GetMaxDistanceBetweenAtoms(glycosite->GetAttachedGlycan()->GetAllAtomsOfAssembly()) + 5); // added five to account for CB-C1(glycan) distance
        AtomVector close_protein_beads = selection::AtomsWithinDistanceOf(cb_atom, distance, protein_beads);
        glycosite->SetProteinBeads(&close_protein_beads);
   //     glycosite->SetProteinBeads(&protein_beads);
    	AtomVector these_beads;
    	ResidueVector glycan_residues = glycosite->GetAttachedGlycan()->GetResidues();
        for (ResidueVector::iterator it2 = glycan_residues.begin(); it2 != glycan_residues.end(); ++it2)
        {
        	Residue *residue = *it2;
            if ( residue->GetName().at(1) != 'S' ) // don't add a bead to a sialic acid in this section (see "else" below). Middle character of resname is always S for sialic acid.
           // if ( residue->GetName().compare(1, 2, "SA") != 0) // don't add one to a sialic acid (see below)
            {
                // std::cout << (*resi_iter)->GetName() << "\tG\t" << (*resi_iter)->CheckIfProtein() << endl;
               // std::cout << "Adding bead to self glycan " << residue->GetId() << std::endl;
                Atom* bead_atom = new Atom(residue, "gfat", residue->GetGeometricCenter());
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
                AtomVector atoms = residue->GetAtoms();
                for (AtomVector::iterator it3 = atoms.begin(); it3 != atoms.end(); ++it3)
                {
                    Atom *atom = *it3;
                    if ( (atom->GetName().compare("C2N") == 0) || (atom->GetName().compare("C6") == 0) )
                    {
                        bead_atom = new Atom(residue, "gfat", atom->GetCoordinates().at(0));
                        residue->AddAtom(bead_atom);
                        these_beads.push_back(bead_atom);
                        any_atom->GetNode()->AddNodeNeighbor(bead_atom);
                        AtomNode *node1 = new AtomNode(); // DELETE IS FOR LOSERS.
                        bead_atom->SetNode(node1);
                        bead_atom->GetNode()->SetNodeNeighbors(temp);
                    }
                }
            }
            else // if it is sialic acid
            {
                AtomVector atoms = residue->GetAtoms();
                for (AtomVector::iterator it3 = atoms.begin(); it3 != atoms.end(); ++it3)
                {
                    Atom *atom = *it3;
                    if ( (atom->GetName().compare("C2") == 0) || (atom->GetName().compare("N5") == 0) || (atom->GetName().compare("C8") == 0) )
                    {
                        Atom* bead_atom = new Atom(residue, "gfat", atom->GetCoordinates().at(0));
                        residue->AddAtom(bead_atom);
                        these_beads.push_back(bead_atom);
                        Atom *any_atom = residue->GetAtoms().at(0);
                        any_atom->GetNode()->AddNodeNeighbor(bead_atom);
                        AtomVector temp = {any_atom};
                        AtomNode *node = new AtomNode(); // DELETE IS FOR LOSERS.
                        bead_atom->SetNode(node);
                        bead_atom->GetNode()->SetNodeNeighbors(temp);
                    }
                }
            }
        }
        glycosite->SetSelfGlycanBeads(&these_beads);
    }
    // Now find beads from other glycans and add them to list of other_glycan_beads for each glycosite
    for (GlycosylationSiteVector::iterator it1 = glycosites.begin(); it1 != glycosites.end(); ++it1)
    {
    	GlycosylationSite *glycosite1 = &(*it1);
    	AtomVector other_glycan_beads;
        for (GlycosylationSiteVector::iterator it2 = glycosites.begin(); it2 != glycosites.end(); ++it2)
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

void Remove_Beads(MolecularModeling::Assembly &glycoprotein)
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

double GetMaxDistanceBetweenAtoms(AtomVector atoms)
{
    double max_distance = 0.0;
    for(AtomVector::iterator it1 = atoms.begin(); it1 != atoms.end(); ++it1)
    {
        Atom *atom1 = (*it1);
        for(AtomVector::iterator it2 = it1; it2 != atoms.end(); ++it2)
        {
            Atom *atom2 = (*it2);
            if (atom1->GetDistanceToAtom(atom2) > max_distance)
            {
                max_distance = atom1->GetDistanceToAtom(atom2);
            }
        }
    }
    return max_distance;
}



