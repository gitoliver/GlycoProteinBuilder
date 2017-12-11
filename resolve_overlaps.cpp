#include <iostream>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <string>

// Change this next line to your PATH:
#include "/home/ubunter/software/gems/gmml/includes/MolecularModeling/assembly.hpp"
#include "/home/ubunter/software/gems/gmml/includes/MolecularModeling/overlaps.hpp"
// #include "/home/oliver/Programs/gems/gmml/includes/MolecularModeling/assembly.hpp"
// #include "/home/oliver/Programs/gems/gmml/includes/MolecularModeling/overlaps.hpp"
#include "resolve_overlaps.h"

using namespace std;
using namespace MolecularModeling;
using namespace GeometryTopology;
using namespace gmml;

// // makes a fat atom vector based on residues of assembly; if residue is a glycan, use midpoint of C1 and C4; radius of 4
// AtomVector FatAtomMode_Glycan(ResidueVector residuevector)
// {
//     AtomVector fat_atom_vector;
//     for(ResidueVector::iterator resi_it = residuevector.begin(); resi_it != residuevector.end(); ++resi_it)
//     {
//         AtomVector atoms = (*resi_it)->GetAtoms();
//         Atom *C1, *C4;
//         bool create_fatom = false;
//         for(AtomVector::iterator atom_iter = atoms.begin(); atom_iter != atoms.end(); ++atom_iter)
//         {
//             if ( (*atom_iter)->GetName().compare("C1")==0 )
//             {
//                 C1 = *atom_iter;
//                 create_fatom = true;
//             }
//             if ( (*atom_iter)->GetName().compare("C4")==0 ) C4 = *atom_iter;
//         }
//         if (create_fatom)
//         {
//             Atom fat_atom;
//             Coordinate new_point;
//             CoordinateVector new_point_vector;
//             new_point.SetX((C1->GetCoordinates().at(0)->GetX() + C4->GetCoordinates().at(0)->GetX()) / 2);
//             new_point.SetY((C1->GetCoordinates().at(0)->GetY() + C4->GetCoordinates().at(0)->GetY()) / 2);
//             new_point.SetZ((C1->GetCoordinates().at(0)->GetZ() + C4->GetCoordinates().at(0)->GetZ()) / 2);
//             new_point_vector.push_back(&new_point);
//             fat_atom.SetCoordinates(new_point_vector);
//             fat_atom.SetName("4");
//             // std::cout << C1->GetCoordinates().at(0)->GetX() << "\n";
//             // std::cout << fat_atom.GetCoordinates().at(0)->GetX() << "\n";
//             fat_atom_vector.push_back(&fat_atom);
//         }
//     }
//     // std::cout << fat_atom_vector.size() << "\n";
//     // std::cout << residuevector.size() << "\n";
//     return fat_atom_vector;
// }
//
// // makes a fat atom vector based on residues of assembly; if residue is a amino acid, use atom CA; radius of 4
// AtomVector FatAtomMode_Protein(ResidueVector residuevector)
// {
//     AtomVector fat_atom_vector;
//     for(ResidueVector::iterator resi_it = residuevector.begin(); resi_it != residuevector.end(); ++resi_it)
//     {
//         AtomVector atoms = (*resi_it)->GetAtoms();
//         Atom *CA;
//         bool create_fatom = false;
//         for(AtomVector::iterator atom_iter = atoms.begin(); atom_iter != atoms.end(); ++atom_iter)
//         {
//             if ( (*atom_iter)->GetName().compare("CA")==0 )
//             {
//                 CA = *atom_iter;
//                 create_fatom = true;
//             }
//         }
//         if (create_fatom)
//         {
//             Atom fat_atom;
//             Coordinate new_point;
//             CoordinateVector new_point_vector;
//             new_point.SetX(CA->GetCoordinates().at(0)->GetX());
//             new_point.SetY(CA->GetCoordinates().at(0)->GetY());
//             new_point.SetZ(CA->GetCoordinates().at(0)->GetZ());
//             new_point_vector.push_back(&new_point);
//             fat_atom.SetCoordinates(new_point_vector);
//             fat_atom.SetName("4");
//             // std::cout << C1->GetCoordinates().at(0)->GetX() << "\n";
//             // std::cout << fat_atom.GetCoordinates().at(0)->GetX() << "\n";
//             fat_atom_vector.push_back(&fat_atom);
//         }
//     }
//     // std::cout << fat_atom_vector.size() << "\n";
//     // std::cout << residuevector.size() << "\n";
//                 std::cout << "----------- section ---\n";
//     // std::cout << fat_atom_vector.at(0)->GetName() << "\n";
//     std::cout << fat_atom_vector.at(0)->GetCoordinates().at(0)->GetZ() << "\n";
//                 std::cout << "----------- section ---\n";
//     return fat_atom_vector;
// }

// makes a new vector of fat atoms from the protein CA
AtomVector FatAtomMode_Protein(Assembly protein)
{
    std::cout << "----------- block --\n";
    AtomVector fatom_protein;
    AtomVector atoms = protein.GetAllAtomsOfAssembly();
    // Residue* dummy_res = new Residue();

    for (AtomVector::iterator atom_iter = atoms.begin(); atom_iter != atoms.end(); ++atom_iter)
    {
        if ( (*atom_iter)->GetName().compare("CA")==0 )
        {
            Atom *fatomCA = new Atom( (*atom_iter)->GetResidue(), "3fatom", (*atom_iter)->GetCoordinates());
            fatom_protein.push_back(fatomCA);

            // Atom *fatomCA = new Atom(dummy_res, "3fatom", (*atom_iter)->GetCoordinates());
            // fatom_protein.push_back(fatomCA);
        }
    }
    return fatom_protein;
}

// makes a new vector of fat atoms from the monosaccharide ring center
AtomVector FatAtomMode_Glycan(Assembly glycan)
{
    AtomVector fatom_glycan;
    ResidueVector monosacs = glycan.GetAllResiduesOfAssembly();
    for (ResidueVector::iterator resi_iter = monosacs.begin(); resi_iter != monosacs.end(); ++resi_iter)
    {
        Atom *fatomRINGCENTER = new Atom( (*resi_iter), "4fatom", (*resi_iter)->GetRingCenter());
        fatom_glycan.push_back(fatomRINGCENTER);
    }
    return fatom_glycan;
}

// adds fat atoms to the assembly
void Implant_FatAtoms(Assembly glycoprotein)
{
    ResidueVector glycoprotein_residues = glycoprotein.GetAllResiduesOfAssembly();
    for(ResidueVector::iterator it1 = glycoprotein_residues.begin(); it1!=glycoprotein_residues.end(); ++it1)
    {
        AtomVector atoms = (*it1)->GetAtoms();
        if ( (*it1)->CheckIfProtein()) // it is a protein
        {
            for(AtomVector::iterator atom_iter = atoms.begin(); atom_iter != atoms.end(); ++atom_iter)
            {
                std::cout << "----------- break -\n";
                Atom atomCA;
                bool hasCA = false;
                if ( (*atom_iter)->GetName().compare("CA")==0 )
                {
                    std::cout << "----------- break\n";
                    atomCA = *atom_iter;
                    hasCA = true;
                }
                if (hasCA)
                {
                    GeometryTopology::Coordinate new_point;

                    new_point.SetX(atomCA.GetCoordinates().at(0)->GetX());
                    new_point.SetY(atomCA.GetCoordinates().at(0)->GetY());
                    new_point.SetZ(atomCA.GetCoordinates().at(0)->GetZ());

                    Atom* fatom = new Atom(*it1, "3f", new_point);
                    (*it1)->AddAtom(fatom);
                    std::cout << "pro\n";
                }
            }
        }
        // if (!(*it1)->CheckIfProtein()) // it is a glycan
        // {
        //     for(AtomVector::iterator atom_iter = atoms.begin(); atom_iter != atoms.end(); ++atom_iter)
        //     {
        //         if ( (*atom_iter)->GetName().compare("C1")==0 )
        //         {
        //             C1 = *atom_iter;
        //         }
        //
        //
        //
        //         if ( (*atom_iter)->GetName().compare("CA")==0 )
        //         {
        //             GeometryTopology::Coordinate new_point;
        //
        //             new_point.SetX((*atom_iter)->GetCoordinates().at(0)->GetX());
        //             new_point.SetY((*atom_iter)->GetCoordinates().at(0)->GetY());
        //             new_point.SetZ((*atom_iter)->GetCoordinates().at(0)->GetZ());
        //
        //             Atom *atomOH = new Atom(*it1, "3f", new_point);
        //             std::cout << "gly\n";
        //         }
        //     }
        // }
    }
}

// taken from gmml/src/MolecularModeling/overlaps.cc
// plan to modify the code to work with fat atoms
double ModifiedCalculateAtomicOverlaps(AtomVector atomsA, AtomVector atomsB)
{
    double distance = 0.0, totalOverlap = 0.0;
    for(AtomVector::iterator it1 = atomsA.begin(); it1 != atomsA.end(); ++it1)
    {
        for(AtomVector::iterator it2 = atomsB.begin(); it2 != atomsB.end(); ++it2)
        {
            Atom *atomA = *it1;
            Atom *atomB = *it2;
            // std::cout << "----------- break ---\n";
            // std::cout << atomA->GetName() << "\n";
            std::cout << atomA->GetCoordinates().at(0)->GetX() << "\n";
            // std::cout << "----------- break ---\n";
            if ( (atomA->GetCoordinates().at(0)->GetX() - atomB->GetCoordinates().at(0)->GetX()) < 1.2 ) // This is faster than calulating distance, and rules out tons of atom pairs.
            {
                distance = atomA->GetDistanceToAtom(atomB);
                // std::cout << "----------- break\n";
                if ( ( distance < 3.6 ) && ( distance > 0.0 ) ) //Close enough to overlap, but not the same atom
                {
                    // std::cout << "----------- break\n";
                    totalOverlap += gmml::CalculateAtomicOverlaps(atomA, atomB); // This calls the overloaded version with default values
                }
            }
        }
    }
    return (totalOverlap / CARBON_SURFACE_AREA); //Normalise to area of a buried carbon
}

// the fat atom score function (work in progress) much faster and messy
ResidueVector ResiFilter_ScoreFatAtoms(Assembly* glycoprotein, GlycoSiteVector* glycosites, double* overlap_score)
{
    ResidueVector filtered_residues_list;
    AtomVector protein = glycoprotein->GetAllAtomsOfAssemblyWithinProteinResidues();
    AtomVector glycans = glycoprotein->GetAllAtomsOfAssemblyNotWithinProteinResidues();

    std::cout << "----------- looky here\n";
    AtomVector fatom_protein = FatAtomMode_Protein(glycoprotein);
    std::cout << "----------- looky here\n";

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
    std::cout << "USING SEED:    " << seed << "\n\n";
    /////////////////// Get pointers to protein and glycan parts of assembly ///
    AtomVector protein = glycoprotein.GetAllAtomsOfAssemblyWithinProteinResidues();
    AtomVector glycans = glycoprotein.GetAllAtomsOfAssemblyNotWithinProteinResidues();
    ////////////////////////////////////////////////////////////////////////////

    ResidueVector residuevector = glycoprotein.GetAllResiduesOfAssembly();


    // AtomVector the_atoms = glycoprotein.GetAllAtomsOfAssembly();
    // for (AtomVector::iterator it = fatoms_protein.begin(); it != fatoms_protein.end(); ++it)
    // {
    //     std::cout << (*it)->GetName() << endl;
    // }


    // AtomVector fatoms_protein = FatAtomMode_Protein(glycoprotein);

    // std::cout << fatoms_protein.size() << "\n";
    // Implant_FatAtoms(glycoprotein);


    int cycle = 1, max_tries = 2;
    while (cycle <= max_tries)
    {
        double overlap_score = 0.0;

        std::cout << "----------- cycle\n";
        // ResidueVector move_these_guys = ResiFilter_ScoreTrueOverlap(&glycoprotein, &glycosites, &overlap_score);
        ResidueVector move_these_guys = ResiFilter_ScoreFatAtoms(&glycoprotein, &glycosites, &overlap_score);

        std::cout << "OVERALL: " << overlap_score << "\n\n";

        ResiRotor_FullRange(&glycoprotein, &move_these_guys);
        cycle++;
    }
}
