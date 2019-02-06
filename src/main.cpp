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
    if (argc == 2)
    {
        working_Directory = argv[1];
    }
    std::cout << "Working directory is " << working_Directory << "\n"; // Assumes folder with glycans inside is present.

    //************************************************//
    // Read input file                                //
    //************************************************//

    //std::string installation_Directory = Find_Program_Installation_Directory();
    GlycosylationSiteVector glycosites;
    std::string proteinPDB, glycanDirectory;
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
    PdbFileSpace::PdbFile *outputPdbFileGlycoProteinAll = glycoprotein.BuildPdbFileStructureFromAssembly(-1,0);
    outputPdbFileGlycoProteinAll->Write(working_Directory + "/GlycoProtein_Initial.pdb");
    //glycoprotein_builder::SetResonableChi1Chi2DihedralAngles(glycosites);



    //************************************************//
    // Prep for Resolve Overlaps                      //
    //************************************************//

    std::cout << "BuildGlycoproteinStructure"  << std::endl;
   // PdbFileSpace::PdbFile *outputPdbFileGlycoProteinAll = glycoprotein.BuildPdbFileStructureFromAssembly(-1,0);
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

    // Wiggle
    bool use_monte_carlo = true;
    glycoprotein_builder::SetDefaultDihedralAnglesUsingMetadata(glycosites); // Reset to reasonable starting points 
    resolve_overlaps::wiggleFirstLinkages(glycosites, 100);
    std::cout << "Real atomic overlaps is " << glycoprotein_builder::CalculateAtomicOverlaps(glycosites) << std::endl;
    resolve_overlaps::weighted_protein_global_overlap_random_descent(glycosites, 100, use_monte_carlo);
    std::cout << "Real atomic overlaps is " << glycoprotein_builder::CalculateAtomicOverlaps(glycosites) << std::endl;
    resolve_overlaps::wiggle(glycosites, 100);
    std::cout << "Real atomic overlaps is " << glycoprotein_builder::CalculateAtomicOverlaps(glycosites) << std::endl;

    resolve_overlaps::wiggleFirstLinkages(glycosites, 100);
    std::cout << "Real atomic overlaps is " << glycoprotein_builder::CalculateAtomicOverlaps(glycosites) << std::endl;
    resolve_overlaps::weighted_protein_global_overlap_random_descent(glycosites, 100, use_monte_carlo);
    std::cout << "Real atomic overlaps is " << glycoprotein_builder::CalculateAtomicOverlaps(glycosites) << std::endl;
    resolve_overlaps::wiggle(glycosites, 100);
    std::cout << "Real atomic overlaps is " << glycoprotein_builder::CalculateAtomicOverlaps(glycosites) << std::endl;

    resolve_overlaps::wiggleFirstLinkages(glycosites, 100);
    std::cout << "Real atomic overlaps is " << glycoprotein_builder::CalculateAtomicOverlaps(glycosites) << std::endl;
    resolve_overlaps::weighted_protein_global_overlap_random_descent(glycosites, 100, use_monte_carlo);
    std::cout << "Real atomic overlaps is " << glycoprotein_builder::CalculateAtomicOverlaps(glycosites) << std::endl;
    resolve_overlaps::wiggle(glycosites, 100);
    std::cout << "Real atomic overlaps is " << glycoprotein_builder::CalculateAtomicOverlaps(glycosites) << std::endl;

    resolve_overlaps::wiggle(glycosites, 300);
    std::cout << "Real atomic overlaps is " << glycoprotein_builder::CalculateAtomicOverlaps(glycosites) << std::endl;

  //  double loose_overlap_tolerance = 1.0;
    outputPdbFileGlycoProteinAll->Write(working_Directory + "/GlycoProtein_All_Resolved.pdb");
    //glycoprotein_builder::DeleteSitesIterativelyWithOverlapAboveTolerance(glycosites, loose_overlap_tolerance);


//    // Testing algorithms:
////    glycoprotein_builder::SetRandomChi1Chi2Values(glycosites);
////    int max_cycles = 100;
////    for(int i=0; i < 5; ++i)
////    {
////        resolve_overlaps::weighted_protein_global_overlap_random_descent(glycosites, max_cycles);
////        resolve_overlaps::weighted_protein_global_overlap_monte_carlo(glycosites, max_cycles);
////    }

    std::cout << "Global overlap is " << glycoprotein_builder::GetGlobalOverlap(glycosites) << "\n";
    beads::Remove_Beads(glycoprotein); //Remove beads and write a final PDB & PRMTOP

    glycoprotein_builder::PrintDihedralAnglesAndOverlapOfGlycosites(glycosites);
    outputPdbFileGlycoProteinAll = glycoprotein.BuildPdbFileStructureFromAssembly(-1,0);
    outputPdbFileGlycoProteinAll->Write(working_Directory + "/GlycoProtein_Resolved.pdb");

    std::cout << "Program got to end ok" << std::endl;
    return 0;
}




