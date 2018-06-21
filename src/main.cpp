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
#include "../includes/genetic_algorithm.h"
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
    GlycosylationSiteVector glycoSites;
    std::string proteinPDB, glycanDirectory;
    std::cout << "Read_Input_File\n";
    glycoprotein_builder::Read_Input_File(glycoSites, proteinPDB, glycanDirectory, working_Directory);

    //************************************************//
    // Load Protein PDB file                          //
    //************************************************//
    std::cout << "Build glycoprotein Structure By Distance" << std::endl;
    Assembly glycoprotein( (working_Directory + "/" + proteinPDB), gmml::InputFileType::PDB);
    glycoprotein.BuildStructureByDistance(4, 1.91); // 4 threads, 1.91 cutoff to allow C-S in Cys and Met to be bonded.

    //************************************************//
    // Load Glycans and Attach to glycosites          //
    //************************************************//

    std::cout << "AttachGlycansToGlycosites"  << std::endl;
    glycoprotein_builder::AttachGlycansToGlycosites(glycoprotein, glycoSites, glycanDirectory);

    //************************************************//
    // Resolve Overlaps                               //
    //************************************************//

    std::cout << "BuildGlycoproteinStructure"  << std::endl;
    PdbFileSpace::PdbFile *outputPdbFileGlycoProteinAll = glycoprotein.BuildPdbFileStructureFromAssembly(-1,0);
    outputPdbFileGlycoProteinAll->Write(working_Directory + "/GlycoProtein_Initial.pdb");
    // Add beads. They make the overlap calculation faster.
    std::cout << "Add_Beads"  << std::endl;

    Add_Beads(glycoprotein, glycoSites);
    std::cout << "protein_first_monte_carlo" << std::endl ;
    resolve_overlaps::protein_first_monte_carlo(glycoSites);
    Remove_Beads(glycoprotein); //Remove beads and write a final PDB & PRMTOP

//    std::cout << "In main, the following sites are in the glycoSite vector:\n";
//    for(GlycosylationSiteVector::iterator current_glycosite = glycoSites.begin(); current_glycosite != glycoSites.end(); ++current_glycosite)
//    {
//        std::cout << current_glycosite->GetResidueNumber() << "\n";
//    }

    outputPdbFileGlycoProteinAll = glycoprotein.BuildPdbFileStructureFromAssembly(-1,0);
    outputPdbFileGlycoProteinAll->Write(working_Directory + "/GlycoProtein_Resolved.pdb");

    std::cout << "Program got to end ok" << std::endl;
    return 0;
}




