#include <iostream>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <string>

// Change this next line to your PATH:
//#include "/home/ubunter/software/gems/gmml/includes/MolecularModeling/assembly.hpp"
//#include "/home/ubunter/software/gems/gmml/includes/MolecularModeling/overlaps.hpp"
#include "../../../gems/gmml/includes/MolecularModeling/assembly.hpp"
#include "../../../gems/gmml/includes/MolecularModeling/overlaps.hpp"
// #include "/home/oliver/Programs/gems/gmml/includes/MolecularModeling/assembly.hpp"
// #include "/home/oliver/Programs/gems/gmml/includes/MolecularModeling/overlaps.hpp"
#include "resolve_overlaps.h"

using namespace std;
using namespace MolecularModeling;
using namespace GeometryTopology;
using namespace gmml;

void print_coords(Atom* atom)
{
    Coordinate coord = atom->GetCoordinates().at(0);
    std::cout << atom->GetName() << "\t" << coord.GetX() << "\t" << coord.GetY() << "\t" << coord.GetZ() << "\n";
}

// adds the fat atoms to the glycoprotein itself (remember to remove them later or something!)
void Implant_FatAtoms(Assembly glycoprotein, GlycoSiteVector glycosites)
{
    ResidueVector all_residues = glycoprotein.GetAllResiduesOfAssembly();
    for (ResidueVector::iterator resi_iter = all_residues.begin(); resi_iter != all_residues.end(); ++resi_iter)
    {
        // std::cout << (*resi_iter)->GetName() << "\t" << (*resi_iter)->CheckIfProtein() << endl;
        if ((*resi_iter)->CheckIfProtein()==1) // the current residue is an amino acid
        {
            // std::cout << (*resi_iter)->GetName() << "\tA\t" << (*resi_iter)->CheckIfProtein() << endl;
            Atom *atomCA;
            AtomVector atoms = (*resi_iter)->GetAtoms();
            for (AtomVector::iterator atom_iter = atoms.begin(); atom_iter != atoms.end(); ++atom_iter)
            {
                if ( (*atom_iter)->GetName().compare("CA")==0 ) atomCA = *atom_iter;
            }
            Atom* fatom = new Atom(*resi_iter, "3fatom", atomCA->GetCoordinates());
            (*resi_iter)->AddAtom(fatom);
        }
        if ((*resi_iter)->CheckIfProtein()==0) // the current residue is an glycan (or something else?)
        {
            // std::cout << (*resi_iter)->GetName() << "\tG\t" << (*resi_iter)->CheckIfProtein() << endl;
            Atom* fatom = new Atom(*resi_iter, "4fatom", (*resi_iter)->GetRingCenter());
            (*resi_iter)->AddAtom(fatom);
        }
    }
}

// removes the fat atoms from the glycoprotein; do this before trying to spit out a pdb file
void Sacrifice_FatAtoms(Assembly glycoprotein)
{
    ResidueVector all_residues = glycoprotein.GetAllResiduesOfAssembly();
    for (ResidueVector::iterator resi_iter = all_residues.begin(); resi_iter != all_residues.end(); ++resi_iter)
    {
        AtomVector atoms = (*resi_iter)->GetAtoms();
        for (AtomVector::iterator atom_iter = atoms.begin(); atom_iter != atoms.end(); ++atom_iter)
        {
            if ((*atom_iter)->GetName().find("fatom")==1)
            {
                // std::cout << "deleting: "<< (*atom_iter)->GetName() << "\n";
                (*resi_iter)->RemoveAtom(*atom_iter);
            }
        }
    }
}

// leaves only fat atoms behind, created for ResiFilter_ScoreFatAtomOverlap()
AtomVector Filter_FatAtoms(AtomVector atoms)
{
    AtomVector fat_atoms_vector;
    for (AtomVector::iterator atom_iter = atoms.begin(); atom_iter != atoms.end(); ++atom_iter)
    {
        if ((*atom_iter)->GetName().find("fatom")==1) fat_atoms_vector.push_back(*atom_iter);
    }
    return fat_atoms_vector;
}

// taken from gmml/src/MolecularModeling/overlaps.cc; plan to modify the code to work with fat atoms

// the fat atom score function (work in progress) much faster and messy
ResidueVector ResiFilter_ScoreFatAtomOverlap(Assembly* glycoprotein, GlycoSiteVector* glycosites, double* overlap_score)
{
    ResidueVector filtered_residues_list;
    AtomVector protein = Filter_FatAtoms(glycoprotein->GetAllAtomsOfAssemblyWithinProteinResidues());
    AtomVector glycans = Filter_FatAtoms(glycoprotein->GetAllAtomsOfAssemblyNotWithinProteinResidues());

    double total_glycoprotein_overlap = 0.0;

    for(GlycoSiteVector::iterator it1 = glycosites->begin(); it1 != glycosites->end(); ++it1)
    {
        GlycosylationSite *current_glycan = *it1;
        double current_glycan_overlap = gmml::CalculateAtomicOverlaps(protein, current_glycan->GetAttachedGlycan()->GetAllAtomsOfAssembly());
        double glycan_overlap_with_protein = current_glycan_overlap; // overlap of the glycan overlaping against protein
        for(GlycoSiteVector::iterator it2 = glycosites->begin(); it2 != glycosites->end(); ++it2)
        {
            GlycosylationSite *comparison_glycan = *it2;
            if (current_glycan != comparison_glycan) // dont overlap against yourself
            { // overlap of the glycan overlaping against other glycans
                current_glycan_overlap += gmml::CalculateAtomicOverlaps(comparison_glycan->GetAttachedGlycan()->GetAllAtomsOfAssembly(), current_glycan->GetAttachedGlycan()->GetAllAtomsOfAssembly());
            }
        }

        std::cout << current_glycan->GetResidue()->GetId() << "\t\t"; // PRINT the identity of the residue!
        total_glycoprotein_overlap += current_glycan_overlap;

        std::cout << "GLYCAN-PROTEIN: " << glycan_overlap_with_protein << "\t\t";
        // std::cout << "TOTAL:   " << current_glycan_overlap << "\n";
        std::cout << "GLYCAN-GLYCAN: " << current_glycan_overlap - glycan_overlap_with_protein << "\n";
        if (current_glycan_overlap > 5) // set a threshold overlap!?
        { // make residue vector of condition: the full overlap overlap is above threshold
            filtered_residues_list.push_back(current_glycan->GetResidue());
        }
        *overlap_score = total_glycoprotein_overlap;
    }
    return filtered_residues_list;
}


// the original score function and filter; single atom resolution and relatively slow
ResidueVector ResiFilter_ScoreTrueOverlap(Assembly* glycoprotein, GlycoSiteVector* glycosites, double* overlap_score)
{
    ResidueVector filtered_residues_list;
    AtomVector protein = glycoprotein->GetAllAtomsOfAssemblyWithinProteinResidues();
    AtomVector glycans = glycoprotein->GetAllAtomsOfAssemblyNotWithinProteinResidues();

    double total_glycoprotein_overlap = 0.0;

    for(GlycoSiteVector::iterator it1 = glycosites->begin(); it1 != glycosites->end(); ++it1)
    {
        GlycosylationSite *current_glycan = *it1;
        double current_glycan_overlap = gmml::CalculateAtomicOverlaps(protein, current_glycan->GetAttachedGlycan()->GetAllAtomsOfAssembly());
        double glycan_overlap_with_protein = current_glycan_overlap; // overlap of the glycan overlaping against protein
        for(GlycoSiteVector::iterator it2 = glycosites->begin(); it2 != glycosites->end(); ++it2)
        {
            GlycosylationSite *comparison_glycan = *it2;
            if (current_glycan != comparison_glycan) // dont overlap against yourself
            { // overlap of the glycan overlaping against other glycans
                current_glycan_overlap += gmml::CalculateAtomicOverlaps(comparison_glycan->GetAttachedGlycan()->GetAllAtomsOfAssembly(), current_glycan->GetAttachedGlycan()->GetAllAtomsOfAssembly());
            }
        }

        std::cout << current_glycan->GetResidue()->GetId() << "\t\t"; // PRINT the identity of the residue!
        total_glycoprotein_overlap += current_glycan_overlap;

        std::cout << "GLYCAN-PROTEIN: " << glycan_overlap_with_protein << "\t\t";
        // std::cout << "TOTAL:   " << current_glycan_overlap << "\n";
        std::cout << "GLYCAN-GLYCAN: " << current_glycan_overlap - glycan_overlap_with_protein << "\n";
        if (current_glycan_overlap > 5) // set a threshold overlap!?
        { // make residue vector of condition: the full overlap overlap is above threshold
            filtered_residues_list.push_back(current_glycan->GetResidue());
        }
        *overlap_score = total_glycoprotein_overlap;
    }
    return filtered_residues_list;
}

// torsion adjuster function, samples 360 deg for chi1 & 2 (1 degree increments)
void ResiRotor_FullRange(Assembly* glycoprotein, ResidueVector* move_these_guys)
{
    for(ResidueVector::iterator it1 = move_these_guys->begin(); it1!=move_these_guys->end(); ++it1)
    {
        AtomVector atoms = (*it1)->GetAtoms();
        Atom *atom1, *atom2, *atom3, *atom4, *atom5;
        if( (*it1)->GetName().compare("ASN")==0 || (*it1)->GetName().compare("NLN")==0 )
        { // (THIS IS A GOOD CODE FOLDING SPOT!) if your residue is an ASN or NLN: move chi1 and chi2!
            for(AtomVector::iterator atom_iter = atoms.begin(); atom_iter != atoms.end(); ++atom_iter)
            {
                if ( (*atom_iter)->GetName().compare("N"  )==0 ) atom1 = *atom_iter;
                if ( (*atom_iter)->GetName().compare("CA" )==0 ) atom2 = *atom_iter;
                if ( (*atom_iter)->GetName().compare("CB" )==0 ) atom3 = *atom_iter;
                if ( (*atom_iter)->GetName().compare("CG" )==0 ) atom4 = *atom_iter;
                if ( (*atom_iter)->GetName().compare("ND2")==0 ) atom5 = *atom_iter;
            }
            double random_dihedral = (rand() % 360) + 1 - 180;
            glycoprotein->SetDihedral(atom1, atom2, atom3, atom4, random_dihedral); // CHI1
            random_dihedral = (rand() % 360) + 1 - 180;
            glycoprotein->SetDihedral(atom2, atom3, atom4, atom5, random_dihedral); // CHI2
        }
        if( (*it1)->GetName().compare("TYR")==0 || (*it1)->GetName().compare("OLY")==0 )
        { // (THIS IS A GOOD CODE FOLDING SPOT!) if your residue is an TYR or OLY: move chi1 and chi2!
            for(AtomVector::iterator atom_iter = atoms.begin(); atom_iter != atoms.end(); ++atom_iter)
            {
                if ( (*atom_iter)->GetName().compare("N"  )==0 ) atom1 = *atom_iter;
                if ( (*atom_iter)->GetName().compare("CA" )==0 ) atom2 = *atom_iter;
                if ( (*atom_iter)->GetName().compare("CB" )==0 ) atom3 = *atom_iter;
                if ( (*atom_iter)->GetName().compare("CG" )==0 ) atom4 = *atom_iter;
                if ( (*atom_iter)->GetName().compare("CD1")==0 ) atom5 = *atom_iter;
            }
            double random_dihedral = (rand() % 360) + 1 - 180;
            glycoprotein->SetDihedral(atom1, atom2, atom3, atom4, random_dihedral); // CHI1
            random_dihedral = (rand() % 360) + 1 - 180;
            glycoprotein->SetDihedral(atom2, atom3, atom4, atom5, random_dihedral); // CHI2
        }
        if( (*it1)->GetName().compare("THR")==0 || (*it1)->GetName().compare("OLT")==0 )
        { // (THIS IS A GOOD CODE FOLDING SPOT!) if your residue is an THR or OLT: move chi1!
            for(AtomVector::iterator atom_iter = atoms.begin(); atom_iter != atoms.end(); ++atom_iter)
            {
                if ( (*atom_iter)->GetName().compare("N"  )==0 ) atom1 = *atom_iter;
                if ( (*atom_iter)->GetName().compare("CA" )==0 ) atom2 = *atom_iter;
                if ( (*atom_iter)->GetName().compare("CB" )==0 ) atom3 = *atom_iter;
                if ( (*atom_iter)->GetName().compare("OG1")==0 ) atom4 = *atom_iter;
            }
            double random_dihedral = (rand() % 360) + 1 - 180;
            glycoprotein->SetDihedral(atom1, atom2, atom3, atom4, random_dihedral); // CHI1
        }
        if( (*it1)->GetName().compare("SER")==0 || (*it1)->GetName().compare("OLS")==0 )
        { // (THIS IS A GOOD CODE FOLDING SPOT!) if your residue is an SER or OLS: move chi1!
            for(AtomVector::iterator atom_iter = atoms.begin(); atom_iter != atoms.end(); ++atom_iter)
            {
                if ( (*atom_iter)->GetName().compare("N"  )==0 ) atom1 = *atom_iter;
                if ( (*atom_iter)->GetName().compare("CA" )==0 ) atom2 = *atom_iter;
                if ( (*atom_iter)->GetName().compare("CB" )==0 ) atom3 = *atom_iter;
                if ( (*atom_iter)->GetName().compare("OG" )==0 ) atom4 = *atom_iter;
            }
            double random_dihedral = (rand() % 360) + 1 - 180;
            glycoprotein->SetDihedral(atom1, atom2, atom3, atom4, random_dihedral); // CHI1
        }
    }
}

// basically the main function that does all the work
void resolve_overlaps::monte_carlo(Assembly glycoprotein, GlycoSiteVector glycosites)
{
    // glycosites contains pointers to the residues in glycoprotein that have a glycan attached to them. GetResidue()
    // Each glycosite can have multiple "rotamers" aka "glycan shapes" that are attached. This is due to an old design plan.
    std::cout << "----------- start ----------\n";
    /////////////////// SEED THE RANDOMNESS N STUFF ////////////////////////////
    int seed = time(NULL);
    srand(seed);
    std::cout << "USING SEED:    " << seed << "\n";

    AtomVector protein = glycoprotein.GetAllAtomsOfAssemblyWithinProteinResidues();     // WITHOUT fat atoms
    AtomVector glycans = glycoprotein.GetAllAtomsOfAssemblyNotWithinProteinResidues();  // WITHOUT fat atoms
    // ResidueVector residuevector = glycoprotein.GetAllResiduesOfAssembly();
    // Implant_FatAtoms(glycoprotein, glycosites);
    // { // this scope be for testing and stuff
    //     ResidueVector all_residues = glycoprotein.GetAllResiduesOfAssembly();
    //     for (ResidueVector::iterator resi_iter = all_residues.begin(); resi_iter != all_residues.end(); ++resi_iter)
    //     {
    //         std::cout << (*resi_iter)->GetName() << "\n";
    //         AtomVector atoms = (*resi_iter)->GetAtoms();
    //         for (AtomVector::iterator atom_iter = atoms.begin(); atom_iter != atoms.end(); ++atom_iter)
    //         {
    //             if ((*atom_iter)->GetName().find("fatom")==1) print_coords(*atom_iter);
    //         }
    //         std::cout << "\n";
    //     }
    // }

    // Sacrifice_FatAtoms(glycoprotein);

    int cycle = 1, max_tries = 2;
    while (cycle <= max_tries)
    {
        double overlap_score = 0.0;

        std::cout << "=========== cycle\n";

        std::cout << "----- norm\n";
        ResidueVector move_these_guys = ResiFilter_ScoreTrueOverlap(&glycoprotein, &glycosites, &overlap_score);

        std::cout << "\n----- fat atom\n";
        Implant_FatAtoms(glycoprotein, glycosites);
        ResidueVector merb_these_guys = ResiFilter_ScoreFatAtomOverlap(&glycoprotein, &glycosites, &overlap_score);
        Sacrifice_FatAtoms(glycoprotein);

        std::cout << "OVERALL: " << overlap_score << "\n\n";

        ResiRotor_FullRange(&glycoprotein, &move_these_guys);
        cycle++;
    }

}
