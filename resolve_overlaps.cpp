#include <iostream>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <string>

// Change this next line to your PATH:
#include "../../gems/gmml/includes/MolecularModeling/assembly.hpp"
#include "../../gems/gmml/includes/MolecularModeling/overlaps.hpp"
#include "resolve_overlaps.h"

using namespace std;
using namespace MolecularModeling;
using namespace GeometryTopology;
using namespace gmml;

// print xyz coords, just used for debugging
void print_coords(Atom* atom)
{
    Coordinate coord = atom->GetCoordinates().at(0);
    std::cout << atom->GetName() << "\t" << coord.GetX() << "\t" << coord.GetY() << "\t" << coord.GetZ() << "\n";
}

// adds the fat atoms to the glycoprotein itself (remember to remove them later)
void Add_FatAtoms(Assembly glycoprotein, GlycoSiteVector glycosites)
{
    ResidueVector all_residues = glycoprotein.GetAllResiduesOfAssembly();
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
                    Atom* fatom = new Atom(residue, "3fat", atom->GetCoordinates());
                    residue->AddAtom(fatom);
                }
            }
        }
        else // the current residue is an glycan (or something else?!?:O)
        {
            if (residue->GetName().compare("SUP") !=0) // don't add one to the superimposition atoms
            {
                // std::cout << (*resi_iter)->GetName() << "\tG\t" << (*resi_iter)->CheckIfProtein() << endl;
                Atom* fatom = new Atom(residue, "4fat", residue->GetGeometricCenter());
                residue->AddAtom(fatom);
                //Bond fatom to any other atom in residue so when glycan is moved, fatom moves too.
                Atom *any_atom = residue->GetAtoms().at(0); // 0 is arbitrary, any atom would do.
                std::cout << "Blow here?" << any_atom->GetId() << std::endl;
                any_atom->GetNode()->AddNodeNeighbor(fatom);
                AtomVector temp = {any_atom};
                AtomNode *node = new AtomNode(); // DELETE IS FOR LOSERS.
                fatom->SetNode(node);
                fatom->GetNode()->SetNodeNeighbors(temp);
            }
        }
    }
}

// removes the fat atoms from the glycoprotein; do this before trying to spit out a pdb file
void Remove_FatAtoms(Assembly glycoprotein)
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

// leaves only fat atoms behind, created for ResiFilter_ScoreFatAtomOverlap()
AtomVector Filter_FatAtoms(AtomVector atoms)
{
    AtomVector fat_atoms_vector;
    for (AtomVector::iterator it1 = atoms.begin(); it1 != atoms.end(); ++it1)
    {
        Atom *atom = *it1;
        if (atom->GetName().find("fat")==1) fat_atoms_vector.push_back(atom);
    }
    return fat_atoms_vector;
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
        overlap = ( 2 * (PI_RADIAN) * radiusA* ( radiusA - distance / 2 - ( ( (radiusA*radiusA) - (radiusB*radiusB) ) / (2 * distance) ) ) );
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
    return (totalOverlap / CARBON_SURFACE_AREA); //Normalise to area of a buried carbon
}

// the fat atom score function (work in progress) much faster and messy
ResidueVector ResiFilter_ScoreFatAtomOverlap(Assembly* glycoprotein, GlycoSiteVector* glycosites, double* overlap_score, double threshold)
{
    ResidueVector filtered_residues_list;
    AtomVector fatom_protein = Filter_FatAtoms(glycoprotein->GetAllAtomsOfAssemblyWithinProteinResidues());
    // AtomVector glycans = Filter_FatAtoms(glycoprotein->GetAllAtomsOfAssemblyNotWithinProteinResidues());
    double total_glycoprotein_overlap = 0.0;

    for(GlycoSiteVector::iterator it1 = glycosites->begin(); it1 != glycosites->end(); ++it1)
    {
        GlycosylationSite *current_glycan = *it1;
        double current_glycan_overlap = modified_CalculateAtomicOverlaps(fatom_protein, Filter_FatAtoms(current_glycan->GetAttachedGlycan()->GetAllAtomsOfAssembly()));
        double glycan_overlap_with_protein = current_glycan_overlap; // overlap of the glycan overlaping against protein

        for(GlycoSiteVector::iterator it2 = glycosites->begin(); it2 != glycosites->end(); ++it2)
        {
            GlycosylationSite *comparison_glycan = *it2;
            if (current_glycan != comparison_glycan) // dont overlap against yourself
            { // overlap of the glycan overlaping against other glycans
                current_glycan_overlap += modified_CalculateAtomicOverlaps(Filter_FatAtoms(comparison_glycan->GetAttachedGlycan()->GetAllAtomsOfAssembly()), Filter_FatAtoms(current_glycan->GetAttachedGlycan()->GetAllAtomsOfAssembly()));
            }
        }

        std::cout << current_glycan->GetResidue()->GetId() << "\t\t"; // PRINT the identity of the residue!
        total_glycoprotein_overlap += current_glycan_overlap;

        std::cout << "GLYCAN-PROTEIN: " << glycan_overlap_with_protein << "\t\t";
        // std::cout << "TOTAL:   " << current_glycan_overlap << "\n";
        std::cout << "GLYCAN-GLYCAN: " << current_glycan_overlap - glycan_overlap_with_protein << "\n";
        if (current_glycan_overlap > threshold) // set a threshold overlap!?
        { // make residue vector of condition: the full overlap overlap is above threshold
            filtered_residues_list.push_back(current_glycan->GetResidue());
        }
        *overlap_score = total_glycoprotein_overlap;
    }
    return filtered_residues_list;
}

// the original score function and filter; single atom resolution and relatively slow
ResidueVector ResiFilter_ScoreTrueOverlap(Assembly* glycoprotein, GlycoSiteVector* glycosites, double* overlap_score, double threshold)
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
        if (current_glycan_overlap > threshold) // set a threshold overlap!?
        { // make residue vector of condition: the full overlap overlap is above threshold
            filtered_residues_list.push_back(current_glycan->GetResidue());
        }
        *overlap_score = total_glycoprotein_overlap;
    }
    return filtered_residues_list;
}

// not actually a score function; simply outputs all glycosylated residues to SHUFFLE with the ResiRotor functions
ResidueVector ResiFilter_Aggregate(Assembly* glycoprotein, GlycoSiteVector* glycosites)
{
    ResidueVector filtered_residues_list;
    for(GlycoSiteVector::iterator it1 = glycosites->begin(); it1 != glycosites->end(); ++it1)
    {
        GlycosylationSite *current_glycan = *it1;
        filtered_residues_list.push_back(current_glycan->GetResidue());
    }
    return filtered_residues_list;
}

// random number generator; allows full range rotation
double RandomAngle_360range()
{
    return (rand() % 360) + 1 - 180;
}

// random number generator; specify a maximum step size relative to a start point
double RandomAngle_PlusMinusX(double start_point, int max_step_size)
{
    return start_point + (rand() % (max_step_size * 2)) - max_step_size;
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
            double random_dihedral = RandomAngle_360range();
            glycoprotein->SetDihedral(atom1, atom2, atom3, atom4, random_dihedral); // CHI1
            random_dihedral = RandomAngle_360range();
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
            double random_dihedral = RandomAngle_360range();
            glycoprotein->SetDihedral(atom1, atom2, atom3, atom4, random_dihedral); // CHI1
            random_dihedral = RandomAngle_360range();
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
            double random_dihedral = RandomAngle_360range();
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
            double random_dihedral = RandomAngle_360range();
            glycoprotein->SetDihedral(atom1, atom2, atom3, atom4, random_dihedral); // CHI1
        }
    }
}

// torsion adjuster function, restricts movement to +/- a given range with RandomAngle_PlusMinusX()
void ResiRotor_BabyStep(Assembly* glycoprotein, ResidueVector* move_these_guys)
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
            double random_dihedral = RandomAngle_PlusMinusX(glycoprotein->CalculateTorsionAngleByAtoms(atom1, atom2, atom3, atom4), 6);
            glycoprotein->SetDihedral(atom1, atom2, atom3, atom4, random_dihedral); // CHI1
            random_dihedral = RandomAngle_PlusMinusX(glycoprotein->CalculateTorsionAngleByAtoms(atom2, atom3, atom4, atom5), 6);
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
            double random_dihedral = RandomAngle_PlusMinusX(glycoprotein->CalculateTorsionAngleByAtoms(atom1, atom2, atom3, atom4), 6);
            glycoprotein->SetDihedral(atom1, atom2, atom3, atom4, random_dihedral); // CHI1
            random_dihedral = RandomAngle_PlusMinusX(glycoprotein->CalculateTorsionAngleByAtoms(atom2, atom3, atom4, atom5), 6);
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
            double random_dihedral = RandomAngle_PlusMinusX(glycoprotein->CalculateTorsionAngleByAtoms(atom1, atom2, atom3, atom4), 6);
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
            double random_dihedral = RandomAngle_PlusMinusX(glycoprotein->CalculateTorsionAngleByAtoms(atom1, atom2, atom3, atom4), 6);
            glycoprotein->SetDihedral(atom1, atom2, atom3, atom4, random_dihedral); // CHI1
        }
    }
}

void write_pdb_file(Assembly glycoprotein, int cycle, string summary_filename, double score)
{
    string pdb_filename = "outputs/pose_" + to_string(cycle) + ".pdb";
    PdbFileSpace::PdbFile *outputPdbFile1 = glycoprotein.BuildPdbFileStructureFromAssembly(-1,0);
    outputPdbFile1->Write(pdb_filename);

    ofstream summary;   // write a file that describes the best conformations found
    summary.open(summary_filename, ios::out | ios::app);
    summary << score << "\t" << "pose_" << cycle << ".pdb\n";
    summary.close();
}

// basically the main function that does all the work
void resolve_overlaps::monte_carlo(Assembly glycoprotein, GlycoSiteVector glycosites)
{
    // glycosites contains pointers to the residues in glycoprotein that have a glycan attached to them. GetResidue()
    std::cout << "----------- start ----------\n";
    /////////////////// SEED THE RANDOMNESS N STUFF ////////////////////////////
    int seed = time(NULL);
    srand(seed);
    std::cout << "USING SEED:    " << seed << "\n";
    std::cout << setprecision(6); // make it so my tabs wont be ruined by the number output
    ////////////////////////////////////////////////////////////////////////////
    AtomVector protein = glycoprotein.GetAllAtomsOfAssemblyWithinProteinResidues();     // WITHOUT fat atoms
    AtomVector glycans = glycoprotein.GetAllAtomsOfAssemblyNotWithinProteinResidues();  // WITHOUT fat atoms

    string summary_filename = "outputs/output_summary.txt";
    double best_score_fat  = -0.1, best_score_normal = -0.1;
    int cycle = 0, cycles_since_last_improvement = 0, max_cycles = 16000;
    std::cout << "\n----- fat atoms\n";
    Add_FatAtoms(glycoprotein, glycosites);
    while (cycle <= max_cycles)
    {
        cycle++;
        cycles_since_last_improvement++;
        std::cout << "============================ CYCLE " << cycle <<"\n";
        bool fine_tune = false;
        ResidueVector move_these_guys;
        double score_to_improve = -0.1, overlap_score = -0.1;

        
        
        move_these_guys = ResiFilter_ScoreFatAtomOverlap(&glycoprotein, &glycosites, &overlap_score, 2.0); // SCORE IT USING FAT ATOMS
        //Remove_FatAtoms(glycoprotein);
        std::cout << "OVERALL: " << overlap_score << "\n\n";
        if (overlap_score < best_score_fat || best_score_fat == -0.1)
        {
            best_score_fat = overlap_score;
            ////////////////////////////////////////////////////////////////////
            std::cout << "----- normal atoms\n";
            ResidueVector move_these_guys_real = ResiFilter_ScoreTrueOverlap(&glycoprotein, &glycosites, &overlap_score, 2.0);
            std::cout << "OVERALL: " << overlap_score << "\n\n";
            write_pdb_file(glycoprotein, cycle, summary_filename, overlap_score);
            best_score_normal = overlap_score;
            ////////////////////////////////////////////////////////////////////
            cycles_since_last_improvement = 0;
            // score_to_improve = overlap_score;
            // fine_tune = true;
        }
        // if (fine_tune = true && cycle > 200)
        // {
        //     while (overlap_score < score_to_improve+120)
        //     {
        //         cycle++;
        //         std::cout << "============================ CYCLE " << cycle <<"\n";
        //         ResidueVector nudge_these_guys = ResiFilter_ScoreTrueOverlap(&glycoprotein, &glycosites, &overlap_score, 2.0);
        //         ResiRotor_BabyStep(&glycoprotein, &nudge_these_guys);
        //         std::cout << "OVERALL: " << overlap_score << "\n\n";
        //         if (overlap_score < score_to_improve)
        //         {
        //             score_to_improve = overlap_score;
        //             write_pdb_file(glycoprotein, cycle, summary_filename, overlap_score);
        //         }
        //     }
        // }
        ResiRotor_FullRange(&glycoprotein, &move_these_guys);
        if (cycles_since_last_improvement > 2700) // allow a certain number of tries before shuffling the protein
        {
            ResidueVector shuffle_these_guys = ResiFilter_Aggregate(&glycoprotein, &glycosites);
            ResiRotor_FullRange(&glycoprotein, &shuffle_these_guys);
        }
    }
}
