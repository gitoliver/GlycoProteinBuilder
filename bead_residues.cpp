#include "bead_residues.h"

void Add_Beads(MolecularModeling::Assembly glycoprotein, GlycoSiteVector glycosites)
{
	AtomVector protein_beads; 
    ResidueVector all_residues = glycoprotein.GetAllResiduesOfAssembly();
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
                    Atom* bead_atom = new Atom(residue, "3fat", atom->GetCoordinates());
                    residue->AddAtom(bead_atom);
                    protein_beads.push_back(bead_atom);
                }
            }
        }
    }
    // Go through all glycosite glycans, add bead in center of each residue, attach it to one other atom in residue.
    // Then set protein_beads
    for (GlycoSiteVector::iterator it1 = glycosites.begin(); it1 != glycosites.end(); ++it1)
    {
    	GlycosylationSite *glycosite = *it1;
    	glycosite->SetProteinBeads(&protein_beads);
    	AtomVector these_beads;
    	ResidueVector glycan_residues = glycosite->GetAttachedGlycan()->GetAllResiduesOfAssembly();
        for (ResidueVector::iterator it2 = all_residues.begin(); it2 != all_residues.end(); ++it2)
        {
        	Residue *residue = *it2;
        	if (residue->GetName().compare("SUP") !=0) // don't add one to the superimposition atoms
            {
                // std::cout << (*resi_iter)->GetName() << "\tG\t" << (*resi_iter)->CheckIfProtein() << endl;
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
    for (GlycoSiteVector::iterator it1 = glycosites.begin(); it1 != glycosites.end(); ++it1)
    {
    	GlycosylationSite *glycosite1 = *it1;
    	AtomVector other_glycan_beads;
    	for (GlycoSiteVector::iterator it2 = glycosites.begin(); it2 != glycosites.end(); ++it2)
    	{
    		GlycosylationSite *glycosite2 = *it2;
    		if(glycosite1 != glycosite2) // Check if same site
    		{
    			// append each other glycosite's beads to list of other_glycan_beads: a.insert(std::end(a), std::begin(b), std::end(b));
    			AtomVector temp = glycosite2->GetSelfGlycanBeads();

    			other_glycan_beads.insert(std::end(other_glycan_beads), std::begin(temp), std::end(temp));

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

double Calculate_Bead_Overlap(AtomVector beads)
{

}

// most of this is copied from gmml; used by ONLY the fat atom mode score; [TWEAK SETTINGS]
double Atomwise_CalculateAtomicOverlaps(Atom *atomA, Atom *atomB, double radiusA, double radiusB)
{
    double distance = atomA->GetDistanceToAtom(atomB);
    if (radiusA == -0.1) // default value is -0.1, but user can provide.
    {
        // element info not usually set, so I look at first letter of atom name. This may be why you're reading this.
        if (atomA->GetName().at(0) == '3') radiusA = 3.00; // for fat atom mode
        if (atomA->GetName().at(0) == '4') radiusA = 4.00; // for fat atom mode
        if (atomA->GetName().at(0) == '5') radiusA = 5.00; // for fat atom mode
        if (atomA->GetName().at(0) == '6') radiusA = 6.00; // for fat atom mode
    }
    if (radiusB == -0.1) // default value is -0.1, but user can provide.
    {
        if (atomB->GetName().at(0) == '3') radiusB = 3.00; // for fat atom mode
        if (atomB->GetName().at(0) == '4') radiusB = 4.00; // for fat atom mode
        if (atomB->GetName().at(0) == '5') radiusB = 5.00; // for fat atom mode
        if (atomB->GetName().at(0) == '6') radiusB = 6.00; // for fat atom mode
    }
    // std::cout << "Distance: " << distance << " radiusA: " << radiusA << " radiusB: " << radiusB << std::endl;
    double overlap = 0.0;
    if (radiusA + radiusB > distance + 0.6)
    { // 0.6 overlap is deemed acceptable. (Copying chimera:)
        // Eqn 1, Rychkov and Petukhov, J. Comput. Chem., 2006, Joint Neighbours. Each atom against each atom, so overlap can be "double" counted. See paper.
        overlap = ( 2 * (gmml::PI_RADIAN) * radiusA* ( radiusA - distance / 2 - ( ( (radiusA*radiusA) - (radiusB*radiusB) ) / (2 * distance) ) ) );
    }
    // std::cout << "Non-normalized Overlap=" << totalOverlap << std::endl;
    return overlap;
}

// taken from gmml/src/MolecularModeling/overlaps.cc; plan to modify the code to work with fat atoms; [TWEAK SETTINGS]
double modified_CalculateAtomicOverlaps(AtomVector atomsA, AtomVector atomsB)
{
    double distance = 0.0, totalOverlap = 0.0;
    for(AtomVector::iterator it1 = atomsA.begin(); it1 != atomsA.end(); ++it1)
    {
        for(AtomVector::iterator it2 = atomsB.begin(); it2 != atomsB.end(); ++it2)
        {
            Atom *atomA = *it1;
            Atom *atomB = *it2;
            if ( (atomA->GetCoordinates().at(0)->GetX() - atomB->GetCoordinates().at(0)->GetX()) < 2.0 ) // This is faster than calulating distance, and rules out tons of atom pairs.
            {
                distance = atomA->GetDistanceToAtom(atomB);
                if ( ( distance < 8.0 ) && ( distance > 0.0 ) ) //Close enough to overlap, but not the same atom
                {
                    totalOverlap += Atomwise_CalculateAtomicOverlaps(atomA, atomB, -0.1, -0.1); // This calls the overloaded version with default values
                }
            }
        }
    }
    return (totalOverlap / gmml::CARBON_SURFACE_AREA); //Normalise to area of a buried carbon
}