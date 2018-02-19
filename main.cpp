#include <fstream>
#include <iostream>
#include <cstdlib> // for exit()

//#include <QApplication>
// Includes for directory reading
#include <string>
#include <dirent.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <stdlib.h>     /* getenv */
#include <fstream>      // std::ifstream

#include "glycosylationsite.h"
#include "io.h"
#include "resolve_overlaps.h"
#include "bead_residues.h"

constexpr auto PI = 3.14159265358979323846;

using namespace MolecularModeling;
typedef std::vector<Atom*> AtomVector;
typedef std::vector<GeometryTopology::Coordinate*> CoordinateVector;
typedef std::vector<Residue*> ResidueVector;
typedef std::vector<Assembly*> AssemblyVector;
typedef GeometryTopology::Coordinate Vector;
typedef std::vector<GlycosylationSite> GlycosylationSiteVector;
typedef std::vector<std::string> StringVector;
//typedef std::vector<AttachedGlycan*> AttachedGlycanVector;

/*******************************************/
/* Function Declarations                   */
/*******************************************/
//void Find_and_Prepare_Protein_Residues_for_Glycosylation(GlycoSiteVector *glycosites, ResidueVector *protein_residues, const std::vector<std::string> *glycositeResidueList);
void AttachGlycans(Assembly *glycoprotein, GlycosylationSiteVector *glycoSites);


int main()
{
    std::string working_Directory = Find_Program_Working_Directory();
    //std::string installation_Directory = Find_Program_Installation_Directory();

    //************************************************//
    // Reading input file                             //
    //************************************************//

    GlycosylationSiteVector glycoSites;
    std::string proteinPDB, glycanDirectory, buffer;
    StringVector glycositeResidueList, listOfGlycans;
    std::ifstream inf (working_Directory + "/inputs/" + "input.txt");
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
            getline(inf, glycanDirectory);
        if(strInput == "Protein Residue, Glycan Name:")
        {
            getline(inf, buffer);
            while(buffer != "END")
            {
                StringVector splitLine = split(buffer, ',');
                glycoSites.emplace_back(splitLine.at(0), splitLine.at(1));
                getline(inf, buffer);
            }
        }
        if(strInput == "Protein Residue list:") // Reads lines until it finds a line with "END"
        {
            getline(inf, buffer);
            while(buffer != "END")
            {
                glycositeResidueList.push_back(buffer);
                std::cout << "Found glycosite:" << buffer << std::endl;
                getline(inf, buffer);
            }
        }
        buffer = "Flushed";
        if(strInput == "Glycan id list:")
        {
            getline(inf, buffer);
            while(buffer != "END")
            {
                listOfGlycans.push_back(buffer);
                //glycoSites.push_back(new GlycosylationSite(buffer));
                std::cout << "Found glycan:" << buffer << std::endl;
                getline(inf, buffer);          
            }
        }
    }
   // GlycosylationSiteVector glycoSites;
    //glycoSites.resize(listOfGlycans.size()); // Efficient

    //************************************************//
    // Load Protein PDB file                          //
    //************************************************//

    Assembly glycoprotein ((working_Directory + "/inputs/" + proteinPDB), gmml::InputFileType::PDB);
    glycoprotein.BuildStructureByDistance();

    AttachGlycans(glycoprotein.GetResidues(), glycositeResidueList, &glycoSites, listOfGlycans);


  //  ResidueVector glycoprotein_residues = ;
  //  Find_and_Prepare_Protein_Residues_for_Glycosylation(&glycoSites, &glycoprotein_residues, &glycositeResidueList);

    //************************************************//
    // Load glycan files from directory               //
    //************************************************//

    std::string directory = working_Directory + "/inputs/" + glycanDirectory ;
    std::cout << "directory: " << directory << std::endl;
    std::string filepath;
    DIR *dp; // A directory stream
    struct dirent *dirp; // Contains file serial number and name (char d_name[])
    struct stat filestat; // Contains info about file, such as device ID, user ID, access time etc

    dp = opendir( directory.c_str() ); //.c_str adds a null character to the end.
    if (dp == NULL)
    {
        std::cout << "Error(" << errno << ") opening " << directory << std::endl;
        return errno;
    }
    while ((dirp = readdir ( dp )))
    {
        filepath = directory + "/" + dirp->d_name;
        // If the file is a directory (or is in some way invalid) we'll skip it
        if (stat( filepath.c_str(), &filestat )) continue; // Is it a valid file?
        if (S_ISDIR( filestat.st_mode ))         continue; // Is it a directory?
        for (GlycosylationSiteVector::iterator it = glycoSites.begin(); it != glycoSites.end(); ++it)
        {
            GlycosylationSite* glycosite = &(*it);
            if (glycosite->GetGlycanName().compare(0, glycosite->GetGlycanName().size(), dirp->d_name, 0, glycosite->GetGlycanName().size()) == 0 )
            {
                Assembly input_glycan(filepath, gmml::InputFileType::PDB);
                input_glycan.BuildStructureByDistance();
                glycosite->AttachGlycan(input_glycan, &glycoprotein);
            }
        }
    }
    closedir( dp );

    
    PdbFileSpace::PdbFile *outputPdbFileGlycoProteinAll = glycoprotein.BuildPdbFileStructureFromAssembly(-1,0);
    outputPdbFileGlycoProteinAll->Write(working_Directory + "/outputs/GlycoProtein.pdb");

    //************************************************//
    // Overlaps                                       //
    //************************************************//
   // Assembly *glycan_atoms = glycoSites.at(0)->GetAttachedGlycan();
   // ResidueVector residues = glycan_atoms->GetResidues();
   // for(ResidueVector::iterator it1 = residues.begin(); it1 != residues.end(); ++it1)
   // {
  //      Residue *residue = *it1;
//std::cout << residue->GetId() << std::endl;
  //  }
    //glycan_atoms->Print();



    //resolve_overlaps::monte_carlo(glycoprotein, glycoSites);
    Add_Beads(&glycoprotein, &glycoSites);
    resolve_overlaps::example_for_Gordon(glycoprotein, glycoSites);

    std::cout << "Program got to end ok" << std::endl;
    return 0;
}

void AttachGlycans(Assembly *glycoprotein, GlycosylationSiteVector *glycoSites)
{
    int i = 0;


}
/*void Find_and_Prepare_Protein_Residues_for_Glycosylation(GlycoSiteVector *glycosites, ResidueVector *protein_residues, const std::vector<std::string> *glycositeResidueList)
{
    int i = 0;
    for(std::vector<std::string>::const_iterator it1 = glycositeResidueList->begin(); it1 != glycositeResidueList->end(); ++it1)
    {
        std::string glycosite = (*it1);
        std::cout << "Preparing glycosite: " << glycosite << std::endl;
        for(ResidueVector::iterator it2 = protein_residues->begin(); it2 != protein_residues->end(); ++it2)
        {
            // Comparing strings is easier:
            Residue *residue = (*it2);
            std::string id = residue->GetId();
            std::string formatted_glycosite = "_" + glycosite + "_";
            // If the residue number in the input file is equal to the current residue number
            //std::cout << "Comparing " << id << " with " << formatted_glycosite << std::endl;
            if( id.compare(5, formatted_glycosite.size(), formatted_glycosite) == 0)
            {
                std::cout << "glycosite: " << glycosite << std::endl;
                std::cout << "glycosite id:" << id << std::endl;
                //glycosites->push_back(residue);
                glycosites->at(i)->SetResidue(residue);
                ++i;
            }
        }
    }
}
*/





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




