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
using namespace gmml;




void resolve_overlaps::monte_carlo(Assembly glycoprotein, GlycoSiteVector glycosites)
{
    std::cout << "----------- start ----------\n";
    // glycosites contains pointers to the residues in glycoprotein that have a glycan attached to them. GetResidue()
    // Each glycosite can have multiple "rotamers" aka "glycan shapes" that are attached. This is due to an old design plan.

    /////////////////// SEED THE RANDOMNESS N STUFF ////////////////////////////
    int seed = time(NULL);
    srand(seed);
    std::cout << "USING SEED:    " << seed << "\n";
    ResidueVector residues = glycoprotein.GetAllResiduesOfAssembly();

    /////////////////// Get pointers to protein and glycan parts of assembly ///
    AtomVector protein = glycoprotein.GetAllAtomsOfAssemblyWithinProteinResidues();
    AtomVector glycans = glycoprotein.GetAllAtomsOfAssemblyNotWithinProteinResidues();

    /////////////////// SCORE IT FIRST /////////////////////////////////////////
    double best_score = gmml::CalculateAtomicOverlaps(protein, glycans); // I've made a function you can pass atomvectors to. USE THIS TO CHECK CLASH
    std::cout << "TOTAL GLYCANS: " << glycosites.size() << "\n";
    std::cout << "INITIAL SCORE: " << best_score << "\n\n\n";

    int cycle = 1, max_tries = 5000;
    while (cycle <= max_tries)
    {
        std::cout << "---------------------\nCYCLE: " << cycle++ << "\n\n";
        /////////////////// SCORING N STUFF ////////////////////////////////////////
        double total_glycoprotein_score =0;
        Residue *worst_residue;
        ResidueVector overall_clash_above_threshold; // the overall clash score is above threshold
        ResidueVector glypro_clash_above_threshold; // the glycan-protein clash score is above threshold
        // ResidueVector glygly_clash_above_threshold; // the glycan-glycan clash score is above threshold
        for(GlycoSiteVector::iterator it1 = glycosites.begin(); it1 != glycosites.end(); ++it1)
        {
            GlycosylationSite *current_glyan = *it1;
            double current_glycan_score = gmml::CalculateAtomicOverlaps(protein, current_glyan->GetAttachedGlycan()->GetAllAtomsOfAssembly());
            double glycan_score_on_protein = current_glycan_score; // score of the glycan clashing against protein
            for(GlycoSiteVector::iterator it2 = glycosites.begin(); it2 != glycosites.end(); ++it2)
            {
                GlycosylationSite *comparison_glyan = *it2;
                if (current_glyan != comparison_glyan) // dont score against yourself
                { // score of the glycan clashing against other glycans
                    current_glycan_score += gmml::CalculateAtomicOverlaps(comparison_glyan->GetAttachedGlycan()->GetAllAtomsOfAssembly(), current_glyan->GetAttachedGlycan()->GetAllAtomsOfAssembly());
                }
            }
            std::cout << current_glyan->GetResidue()->GetId() << "\t\t"; // PRINT the identity of the residue!
            total_glycoprotein_score += current_glycan_score;
            std::cout << "CLASH SURFACE: " << glycan_score_on_protein << "\t\t";
            if (glycan_score_on_protein > 100) // set a threshold score!? (THIS CRITERIA IS CURRENTLY NOT BEING USED)
            { // make residue vector of condition: the glycan-protein clash score is above threshold
                glypro_clash_above_threshold.push_back(current_glyan->GetResidue());
            }
            std::cout << "CLASH OVERALL:   " << current_glycan_score << "\n";
            if (current_glycan_score > 5) // set a threshold score!?
            { // make residue vector of condition: the full clash score is above threshold
                overall_clash_above_threshold.push_back(current_glyan->GetResidue());
            }
        }
        std::cout << "CURRENT POSE:  " << total_glycoprotein_score << "\n";
        std::cout << "VECTOR SIZE-full:  " << overall_clash_above_threshold.size() << "\t";
        std::cout << "VECTOR SIZE-glyp:  " << glypro_clash_above_threshold.size() << "\n\n";

        if (overall_clash_above_threshold.size()==0) // if the score is good, output a pdb file, else do it again!
        {
            std::cout << "BRING IT IN BOIZ, WE GOT OURSELVES A BIG ONE!!!!!!!!!!!!!!!!!\n\n";
            PdbFileSpace::PdbFile *outputPdbFile1 = glycoprotein.BuildPdbFileStructureFromAssembly(-1,0);
            outputPdbFile1->Write("outputs/Glycoprotein_even_better_one.pdb");
            break;
        }
        // PdbFileSpace::PdbFile *outputPdbFile1 = glycoprotein.BuildPdbFileStructureFromAssembly(-1,0);
        // outputPdbFile1->Write("outputs/wayland_stuff/Glycoprotein"+to_string(cycle-1)+".pdb");

        /////////////////// MOVE GLYCAN CHAINS N STUFF /////////////////////////////
        ResidueVector move_these_guys =  overall_clash_above_threshold;
        for(ResidueVector::iterator it1 = move_these_guys.begin(); it1!=move_these_guys.end(); ++it1)
        {
            AtomVector atoms = (*it1)->GetAtoms();
            Atom *atom1, *atom2, *atom3, *atom4, *atom5;
            if( (*it1)->GetName().compare("ASN")==0 || (*it1)->GetName().compare("NLN")==0 )
            { // (THIS IS A GOOD CODE FOLDING SPOT!) if your residue is an ASN or NLN: move chi1 and chi2!
                for(AtomVector::iterator atom_iter = atoms.begin(); atom_iter != atoms.end(); ++atom_iter)
                {
                    if ( (*atom_iter)->GetName().compare("N")==0 )
                    {
                        atom1 = *atom_iter;
                    }
                    if ( (*atom_iter)->GetName().compare("CA")==0 )
                    {
                        atom2 = *atom_iter;
                    }
                    if ( (*atom_iter)->GetName().compare("CB")==0 )
                    {
                        atom3 = *atom_iter;
                    }
                    if ( (*atom_iter)->GetName().compare("CG")==0 )
                    {
                        atom4 = *atom_iter;
                    }
                    if ( (*atom_iter)->GetName().compare("ND2")==0 )
                    {
                        atom5 = *atom_iter;
                    }
                }
                double random_dihedral = (rand() % 360) + 1 - 180;
                glycoprotein.SetDihedral(atom1, atom2, atom3, atom4, random_dihedral); // CHI1
                random_dihedral = (rand() % 360) + 1 - 180;
                glycoprotein.SetDihedral(atom2, atom3, atom4, atom5, random_dihedral); // CHI2
            }
            if( (*it1)->GetName().compare("TYR")==0 || (*it1)->GetName().compare("OLY")==0 )
            { // (THIS IS A GOOD CODE FOLDING SPOT!) if your residue is an TYR or OLY: move chi1 and chi2!
                for(AtomVector::iterator atom_iter = atoms.begin(); atom_iter != atoms.end(); ++atom_iter)
                {
                    if ( (*atom_iter)->GetName().compare("N")==0 )
                    {
                        atom1 = *atom_iter;
                    }
                    if ( (*atom_iter)->GetName().compare("CA")==0 )
                    {
                        atom2 = *atom_iter;
                    }
                    if ( (*atom_iter)->GetName().compare("CB")==0 )
                    {
                        atom3 = *atom_iter;
                    }
                    if ( (*atom_iter)->GetName().compare("CG")==0 )
                    {
                        atom4 = *atom_iter;
                    }
                    if ( (*atom_iter)->GetName().compare("CD1")==0 )
                    {
                        atom5 = *atom_iter;
                    }
                }
                double random_dihedral = (rand() % 360) + 1 - 180;
                glycoprotein.SetDihedral(atom1, atom2, atom3, atom4, random_dihedral); // CHI1
                random_dihedral = (rand() % 360) + 1 - 180;
                glycoprotein.SetDihedral(atom2, atom3, atom4, atom5, random_dihedral); // CHI2
            }
            if( (*it1)->GetName().compare("THR")==0 || (*it1)->GetName().compare("OLT")==0 )
            { // (THIS IS A GOOD CODE FOLDING SPOT!) if your residue is an THR or OLT: move chi1!
                for(AtomVector::iterator atom_iter = atoms.begin(); atom_iter != atoms.end(); ++atom_iter)
                {
                    if ( (*atom_iter)->GetName().compare("N")==0 )
                    {
                        atom1 = *atom_iter;
                    }
                    if ( (*atom_iter)->GetName().compare("CA")==0 )
                    {
                        atom2 = *atom_iter;
                    }
                    if ( (*atom_iter)->GetName().compare("CB")==0 )
                    {
                        atom3 = *atom_iter;
                    }
                    if ( (*atom_iter)->GetName().compare("OG1")==0 )
                    {
                        atom4 = *atom_iter;
                    }
                }
                double random_dihedral = (rand() % 360) + 1 - 180;
                glycoprotein.SetDihedral(atom1, atom2, atom3, atom4, random_dihedral); // CHI1
            }
            if( (*it1)->GetName().compare("SER")==0 || (*it1)->GetName().compare("OLS")==0 )
            { // (THIS IS A GOOD CODE FOLDING SPOT!) if your residue is an SER or OLS: move chi1!
                for(AtomVector::iterator atom_iter = atoms.begin(); atom_iter != atoms.end(); ++atom_iter)
                {
                    if ( (*atom_iter)->GetName().compare("N")==0 )
                    {
                        atom1 = *atom_iter;
                    }
                    if ( (*atom_iter)->GetName().compare("CA")==0 )
                    {
                        atom2 = *atom_iter;
                    }
                    if ( (*atom_iter)->GetName().compare("CB")==0 )
                    {
                        atom3 = *atom_iter;
                    }
                    if ( (*atom_iter)->GetName().compare("OG")==0 )
                    {
                        atom4 = *atom_iter;
                    }
                }
                double random_dihedral = (rand() % 360) + 1 - 180;
                glycoprotein.SetDihedral(atom1, atom2, atom3, atom4, random_dihedral); // CHI1
            }
        }
    }
}
