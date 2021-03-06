#include <iostream>
#include <cstdlib> // for exit()
// Includes for directory reading
#include <string>
#include <dirent.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <stdlib.h>     /* getenv */
#include <fstream>      // std::ifstream
#include <cstdlib>
#include <iomanip>
#include <vector>
#include <stdio.h>
#include <cstring>
#include <algorithm> //erase

#include "../includes/io.h"
#include "../includes/resolve_overlaps.h"
#include "../includes/bead_residues.h"
#include "../includes/glycoprotein_builder.h"

constexpr auto PI = 3.14159265358979323846;

using namespace MolecularModeling;

int main(int argc, char* argv[])
{
    std::string working_Directory = Find_Program_Working_Directory(); // Default behaviour.
    std::string installation_Directory = Find_Program_Installation_Directory();
    if (argc == 2)
    {
        working_Directory = argv[1];
    }
    std::cout << "Working directory is " << working_Directory << "\n";
    std::cout << "Install directory is " << installation_Directory << "\n"; // Assumes folder with glycans inside is present.

    //************************************************//
    // Read input file                                //
    //************************************************//

    GlycosylationSiteVector glycosites;
    std::string proteinPDB, glycanDirectory = installation_Directory + "/../glycans/"; // Default behaviour
   // std::cout << "glycanDirectory: " << glycanDirectory << std::endl;
    std::cout << "Read_Input_File\n";
    glycoprotein_builder::Read_Input_File(glycosites, proteinPDB, glycanDirectory, working_Directory);

    //************************************************//
    // Load Protein PDB file                          //
    //************************************************//
    std::cout << "Build protein structure by distance" << std::endl;
    Assembly glycoprotein( (working_Directory + "/" + proteinPDB), gmml::InputFileType::PDB);
    glycoprotein.BuildStructureByDistance(4, 1.6); // 4 threads, 1.91 cutoff to allow C-S in Cys and Met to be bonded. Nope that did bad things

    //************************************************//
    // Load Glycans and Attach to glycosites          //
    //************************************************//

    std::cout << "AttachGlycansToGlycosites"  << std::endl;
    glycoprotein_builder::AttachGlycansToGlycosites(glycoprotein, glycosites, glycanDirectory);
    PdbFileSpace::PdbFile *outputPdbFileGlycoProteinAll = glycoprotein.BuildPdbFileStructureFromAssembly(-1,0,-1,false);
    outputPdbFileGlycoProteinAll->Write(working_Directory + "/GlycoProtein_Initial.pdb");
    //glycoprotein_builder::SetResonableChi1Chi2DihedralAngles(glycosites);



    //************************************************//
    // Prep for Resolve Overlaps                      //
    //************************************************//

    std::cout << "BuildGlycoproteinStructure"  << std::endl;
   // PdbFileSpace::PdbFile *outputPdbFileGlycoProteinAll = glycoprotein.BuildPdbFileStructureFromAssembly(-1,0,-1,false);
    outputPdbFileGlycoProteinAll->Write(working_Directory + "/GlycoProtein_Initial.pdb");
    // Add beads. They make the overlap calculation faster.
    std::cout << "Add_Beads"  << std::endl;
    beads::Add_Beads(glycoprotein, glycosites);
    glycoprotein_builder::UpdateAtomsThatMoveInLinkages(glycosites); // Must update to include the beads

    //************************************************//
    // Resolve Overlaps                               //
    //************************************************//

    // Fast and stupid:
//    if (!resolve_overlaps::dumb_random_walk(glycosites))
//    {
//        int max_cycles = 500;
//        std::cout << "Could not resolve quickly" << std::endl;
//        glycoprotein_builder::SetDefaultDihedralAnglesUsingMetadata(glycosites); // Reset to reasonable starting points
//        outputPdbFileGlycoProteinAll->Write(working_Directory + "/GlycoProtein_Initial1.pdb");
//        resolve_overlaps::weighted_protein_global_overlap_random_descent(glycosites, max_cycles);
//    }
//    resolve_overlaps::rotamer_permutator(glycosites);
////    glycoprotein_builder::SetDefaultDihedralAnglesUsingMetadata(glycosites); // Reset to reasonable starting points
    bool use_monte_carlo = true;
    int cycles = 100;
    resolve_overlaps::wiggleFirstLinkages(glycosites, BEAD, cycles);
    std::cout << "1. Post WiggleFirst Overlaps Bead: " << glycoprotein_builder::CalculateOverlaps(glycosites, BEAD) << ". Atomic: " << glycoprotein_builder::CalculateOverlaps(glycosites, ATOMIC) << std::endl;
    resolve_overlaps::weighted_protein_global_overlap_random_descent(glycosites, BEAD, cycles, use_monte_carlo);
    std::cout << "2. Post Monte Carlo Overlaps Bead: " << glycoprotein_builder::CalculateOverlaps(glycosites, BEAD) << ". Atomic: " << glycoprotein_builder::CalculateOverlaps(glycosites, ATOMIC) << std::endl;
    resolve_overlaps::wiggle(glycosites, BEAD, cycles);
    std::cout << "3. Post Wiggle Overlaps Bead: " << glycoprotein_builder::CalculateOverlaps(glycosites, BEAD) << ". Atomic: " << glycoprotein_builder::CalculateOverlaps(glycosites, ATOMIC) << std::endl;

    resolve_overlaps::wiggleFirstLinkages(glycosites, BEAD, cycles);
    std::cout << "4. Post WiggleFirst Overlaps Bead: " << glycoprotein_builder::CalculateOverlaps(glycosites, BEAD) << ". Atomic: " << glycoprotein_builder::CalculateOverlaps(glycosites, ATOMIC) << std::endl;
    resolve_overlaps::weighted_protein_global_overlap_random_descent(glycosites, BEAD, cycles, use_monte_carlo);
    std::cout << "5. Post Monte Carlo Overlaps Bead: " << glycoprotein_builder::CalculateOverlaps(glycosites, BEAD) << ". Atomic: " << glycoprotein_builder::CalculateOverlaps(glycosites, ATOMIC) << std::endl;
    resolve_overlaps::wiggle(glycosites, BEAD, cycles);
    std::cout << "6. Post Wiggle Overlaps Bead: " << glycoprotein_builder::CalculateOverlaps(glycosites, BEAD) << ". Atomic: " << glycoprotein_builder::CalculateOverlaps(glycosites, ATOMIC) << std::endl;

    resolve_overlaps::wiggleFirstLinkages(glycosites, ATOMIC, cycles);
    std::cout << "7. Post WiggleFirst Overlaps Bead: " << glycoprotein_builder::CalculateOverlaps(glycosites, BEAD) << ". Atomic: " << glycoprotein_builder::CalculateOverlaps(glycosites, ATOMIC) << std::endl;
    resolve_overlaps::weighted_protein_global_overlap_random_descent(glycosites, ATOMIC, cycles, use_monte_carlo);
    std::cout << "8. Post Monte Carlo Overlaps Bead: " << glycoprotein_builder::CalculateOverlaps(glycosites, BEAD) << ". Atomic: " << glycoprotein_builder::CalculateOverlaps(glycosites, ATOMIC) << std::endl;
    resolve_overlaps::wiggle(glycosites, ATOMIC, cycles);
    std::cout << "9. Post Wiggle Overlaps Bead: " << glycoprotein_builder::CalculateOverlaps(glycosites, BEAD) << ". Atomic: " << glycoprotein_builder::CalculateOverlaps(glycosites, ATOMIC) << std::endl;

    glycoprotein_builder::CalculateOverlaps(glycosites);
    std::cout << "Pre remove beads overlap: " << glycoprotein_builder::GetGlobalOverlap(glycosites) << "\n";
    beads::Remove_Beads(glycoprotein); //Remove beads and write a final PDB & PRMTOP

    outputPdbFileGlycoProteinAll = glycoprotein.BuildPdbFileStructureFromAssembly(-1,0,-1,false);
    outputPdbFileGlycoProteinAll->Write(working_Directory + "/GlycoProtein_All_Resolved.pdb");
    double loose_overlap_tolerance = 1.0;
    glycoprotein_builder::DeleteSitesIterativelyWithAtomicOverlapAboveTolerance(glycosites, loose_overlap_tolerance);


//    // Testing algorithms:
////    glycoprotein_builder::SetRandomChi1Chi2Values(glycosites);
////    int max_cycles = 100;
////    for(int i=0; i < 5; ++i)
////    {
////        resolve_overlaps::weighted_protein_global_overlap_random_descent(glycosites, max_cycles);
////        resolve_overlaps::weighted_protein_global_overlap_monte_carlo(glycosites, max_cycles);
////    }

    std::cout << "Atomic overlap is " << glycoprotein_builder::CalculateOverlaps(glycosites, ATOMIC) << "\n";
    std::cout << "Global overlap is " << glycoprotein_builder::GetGlobalOverlap(glycosites) << "\n";

    glycoprotein_builder::PrintDihedralAnglesAndOverlapOfGlycosites(glycosites);
    outputPdbFileGlycoProteinAll = glycoprotein.BuildPdbFileStructureFromAssembly(-1,0,-1,false);
    outputPdbFileGlycoProteinAll->Write(working_Directory + "/GlycoProtein_Resolved.pdb");

    std::cout << "Program got to end ok" << std::endl;
    return 0;
}




