#include <iostream>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <string>

// Change this next line to your PATH:
#include "../../gems/gmml/includes/MolecularModeling/assembly.hpp"
#include "../../gems/gmml/includes/MolecularModeling/overlaps.hpp"
#include "resolve_overlaps.h"
#include "bead_residues.h"

using namespace std;
using namespace MolecularModeling;
using namespace GeometryTopology;
using namespace gmml;


// leaves only fat atoms behind, created for ResiFilter_ScoreFatAtomOverlap()
AtomVector Filter_Beads(AtomVector atoms)
{
    AtomVector fat_atoms_vector;
    for (AtomVector::iterator it1 = atoms.begin(); it1 != atoms.end(); ++it1)
    {
        Atom *atom = *it1;
        if (atom->GetName().find("fat")==1) fat_atoms_vector.push_back(atom);
    }
    return fat_atoms_vector;
}



// the fat atom score function (work in progress) much faster and messy
ResidueVector ResiFilter_ScoreFatAtomOverlap(Assembly* glycoprotein, GlycoSiteVector* glycosites, double* overlap_score, double threshold)
{
    ResidueVector filtered_residues_list;
    AtomVector fatom_protein = Filter_Beads(glycoprotein->GetAllAtomsOfAssemblyWithinProteinResidues());
    // AtomVector glycans = Filter_Beads(glycoprotein->GetAllAtomsOfAssemblyNotWithinProteinResidues());
    double total_glycoprotein_overlap = 0.0;

    for(GlycoSiteVector::iterator it1 = glycosites->begin(); it1 != glycosites->end(); ++it1)
    {
        GlycosylationSite *current_glycan = *it1;
        double current_glycan_overlap = modified_CalculateAtomicOverlaps(fatom_protein, Filter_Beads(current_glycan->GetAttachedGlycan()->GetAllAtomsOfAssembly()));
        double glycan_overlap_with_protein = current_glycan_overlap; // overlap of the glycan overlaping against protein

        for(GlycoSiteVector::iterator it2 = glycosites->begin(); it2 != glycosites->end(); ++it2)
        {
            GlycosylationSite *comparison_glycan = *it2;
            if (current_glycan != comparison_glycan) // dont overlap against yourself
            { // overlap of the glycan overlaping against other glycans
                current_glycan_overlap += modified_CalculateAtomicOverlaps(Filter_Beads(comparison_glycan->GetAttachedGlycan()->GetAllAtomsOfAssembly()), Filter_Beads(current_glycan->GetAttachedGlycan()->GetAllAtomsOfAssembly()));
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

void resolve_overlaps::example_for_Gordon(Assembly glycoprotein, GlycoSiteVector glycosites)
{
    double site_total_overlap = 0.0, site_glycan_overlap = 0.0, site_protein_overlap = 0.0, new_dihedral_value = 0.0;
    std::cout << "Site | Total | Protein | Glycan " << std::endl;
    for (GlycoSiteVector::iterator it1 = glycosites.begin(); it1 != glycosites.end(); ++it1)
    {
        GlycosylationSite *current_glycosite = *it1; //I always do this for readability.
        site_total_overlap = current_glycosite->calculate_bead_overlaps(); // Must repeat after rotating chi1, chi2.
        site_glycan_overlap = current_glycosite->GetGlycanOverlap(); // If you wish to have this level of detail
        site_protein_overlap = current_glycosite->GetProteinOverlap(); // If you wish to have this level of detail
        std::cout << current_glycosite->GetResidue()->GetName() << " | " << site_total_overlap << " | " << site_protein_overlap  << " | " << site_glycan_overlap << std::endl;
        new_dihedral_value = RandomAngle_PlusMinusX(current_glycosite->GetChi1Value(), 6);
        current_glycosite->SetChi1Value(new_dihedral_value, &glycoprotein);



    }
    std::cout << "Site | Total | Protein | Glycan " << std::endl;
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
    Add_Beads(glycoprotein, glycosites);
    while (cycle <= max_cycles)
    {
        cycle++;
        cycles_since_last_improvement++;
        std::cout << "============================ CYCLE " << cycle <<"\n";
        bool fine_tune = false;
        ResidueVector move_these_guys;
        double score_to_improve = -0.1, overlap_score = -0.1;

        
        
        move_these_guys = ResiFilter_ScoreFatAtomOverlap(&glycoprotein, &glycosites, &overlap_score, 2.0); // SCORE IT USING FAT ATOMS
        //Remove_Beads(glycoprotein);
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
