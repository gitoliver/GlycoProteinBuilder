#include <fstream>
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

constexpr auto PI = 3.14159265358979323846;

using namespace MolecularModeling;

/*******************************************/
/* Function Declarations                   */
/*******************************************/
void Read_Input_File(GlycosylationSiteVector &glycoSites, std::string &proteinPDB, std::string &glycanDirectory, const std::string working_Directory);
void AttachGlycansToGlycosites(Assembly &glycoprotein, GlycosylationSiteVector &glycoSites, const std::string glycanDirectory);

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
    Read_Input_File(glycoSites, proteinPDB, glycanDirectory, working_Directory);

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
    AttachGlycansToGlycosites(glycoprotein, glycoSites, glycanDirectory);

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
    //resolve_overlaps::genetic_algorithm(&glycoprotein, &glycoSites);
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

/*******************************************/
/* Functions                               */
/*******************************************/

void AttachGlycansToGlycosites(Assembly &glycoprotein, GlycosylationSiteVector &glycoSites, const std::string glycanDirectory)
{
    // Find protein residues in Glycoprotein that will get a glycan added. Set Residue in Glycosite.
    // Often there are multiple chains with the same residue number. If the input file repeats a residue number, it will go on the next instance (next chain)
    // of that residue number. This is not a good long term solution, as users will want to select chains.
    ResidueVector protein_residues = glycoprotein.GetResidues();
    for(GlycosylationSiteVector::iterator glycosite = glycoSites.begin(); glycosite != glycoSites.end(); ++glycosite)
    {
        bool stop = false;
        std::string glycosite_number = glycosite->GetResidueNumber();
        ResidueVector::iterator it2 = protein_residues.begin();
        while ( (!stop) && (it2 != protein_residues.end()) ) {
           // for (ResidueVector::iterator it2 = protein_residues.begin(); it2 != protein_residues.end(); ++it2)
            {
                Residue *protein_residue = *it2;
                std::string id = protein_residue->GetId();
                std::string formatted_glycosite_number = "_" + glycosite_number + "_";
                if( id.compare(5, formatted_glycosite_number.size(), formatted_glycosite_number) == 0)
                {
                    //std::cout << "glycosite: " << glycosite_number << std::endl;
                    std::cout << "glycosite id:" << id << std::endl;
                    glycosite->SetResidue(protein_residue);
                    stop = true;
                    // Remove residue from list, so if user has listed same residue name twice, it will go on next instance (i.e. on next chain) of residue number
                    protein_residues.erase(std::remove(protein_residues.begin(), protein_residues.end(), *it2), protein_residues.end()); // Note need #include <algorithm>
                }
                ++it2;
            }
        }
    }
//    // Find protein residues in Glycoprotein that will get a glycan added. Set Residue in Glycosite.
//    ResidueVector protein_residues = glycoprotein.GetResidues();
//    for (ResidueVector::iterator it2 = protein_residues.begin(); it2 != protein_residues.end(); ++it2)
//    {
//        Residue *protein_residue = *it2;
//        std::string id = protein_residue->GetId();
//        for(GlycosylationSiteVector::iterator glycosite = glycoSites.begin(); glycosite != glycoSites.end(); ++glycosite)
//        {
//            std::string glycosite_number = glycosite->GetResidueNumber();
//            std::string formatted_glycosite_number = "_" + glycosite_number + "_";
//            if( id.compare(5, formatted_glycosite_number.size(), formatted_glycosite_number) == 0)
//            {
//                //std::cout << "glycosite: " << glycosite_number << std::endl;
//                std::cout << "glycosite id:" << id << std::endl;
//                glycosite->SetResidue(protein_residue);
//            }
//        }
//    }
    // Load glycan files from directory             
    //std::cout << "Glycan directory: " << glycanDirectory << std::endl;
    std::string filepath;
    DIR *dp; // A directory stream
    struct dirent *dirp; // Contains file serial number and name (char d_name[])
    struct stat filestat; // Contains info about file, such as device ID, user ID, access time etc

    dp = opendir( glycanDirectory.c_str() ); //.c_str adds a null character to the end.
    if (dp == NULL)
    {
        std::cout << "Error(" << errno << ") opening " << glycanDirectory << std::endl;
        return;
    }
    while ((dirp = readdir ( dp )))
    {
        filepath = glycanDirectory + "/" + dirp->d_name;
        // If the file is a directory (or is in some way invalid) we'll skip it
        if (stat( filepath.c_str(), &filestat )) continue; // Is it a valid file?
        if (S_ISDIR( filestat.st_mode ))         continue; // Is it a directory?
        for (GlycosylationSiteVector::iterator glycosite = glycoSites.begin(); glycosite != glycoSites.end(); ++glycosite)
        {
            //std::cout << "Glycan is " << glycosite->GetGlycanName() << ". d_name is " << dirp->d_name << std::endl;
            if (glycosite->GetGlycanName().compare(0, glycosite->GetGlycanName().size(), dirp->d_name, 0, glycosite->GetGlycanName().size()) == 0 )
            {
                Assembly input_glycan(filepath, gmml::InputFileType::PDB);               
                input_glycan.BuildStructureByDistance();
                glycosite->AttachGlycan(input_glycan, glycoprotein);
                //std::cout << "Added " << glycosite->GetGlycanName() << " to " << glycosite->GetResidueNumber();
            }
        }
    }
    closedir( dp ); 
}


void Read_Input_File(GlycosylationSiteVector &glycoSites, std::string &proteinPDB, std::string &glycanDirectory, const std::string working_Directory)
{
    std::string buffer;
    std::ifstream inf (working_Directory + "/input.txt");
    if (!inf)
    {
        std::cerr << "Uh oh, input file could not be opened for reading!" << std::endl;
        std::exit(1);
    }
    while (inf) // While there's still stuff left to read
    {
        std::string strInput;
        getline(inf, strInput);
        if(strInput == "Protein:")
            getline(inf, proteinPDB);
        if(strInput == "Glycans:")
        {
            getline(inf, glycanDirectory);
            glycanDirectory = working_Directory + "/" + glycanDirectory;
        }
        if(strInput == "Protein Residue, Glycan Name:")
        {
            getline(inf, buffer);
            while(buffer != "END")
            {
                StringVector splitLine = split(buffer, ',');
                glycoSites.emplace_back(splitLine.at(1), splitLine.at(0)); // Creates GlycosylationSite instance on the vector. Love it.
                getline(inf, buffer);
            }
        }
    }
}

/*void FindConnectedAtoms(Atom *atom, Assembly::AtomVector &visitedAtoms){

    visitedAtoms.push_back(atom);

    //std::cout << "AtomNodes must be set by an e.g. Assembly::BuildStructureByDistance() call or a SegFault will happen now" << std::endl;
    Assembly::AtomVector neighbors = atom->GetNode()->GetNodeNeighbors();
    //std::cout << "Okay, AtomNodes were set" << std::endl;
    bool alreadyVisited = false;


   // std::cout << "Current atom is " << atom->GetName() << "." << atom->GetId() << std::endl;

    for (auto &neighbor : neighbors){
        alreadyVisited = false; // reset for each neighbor
        for (auto &visitedAtom : visitedAtoms) { // check against each visited atom
            if ( neighbor->GetId() == visitedAtom->GetId() )
                alreadyVisited = true;
        }
        if (!alreadyVisited) {
            std::cout << "Found unvisited neighbor, Going to " << neighbor->GetId() << " from " << atom->GetId() << std::endl;
            FindConnectedAtoms(neighbor, visitedAtoms); // recursive function call
        }
    }
    return;
}
*/




