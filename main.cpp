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

constexpr auto PI = 3.14159265358979323846;

using namespace MolecularModeling;
typedef std::vector<Atom*> AtomVector;
typedef std::vector<GeometryTopology::Coordinate*> CoordinateVector;
typedef std::vector<Residue*> ResidueVector;
typedef std::vector<Assembly*> AssemblyVector;
typedef GeometryTopology::Coordinate Vector;
typedef std::vector<GlycosylationSite*> GlycoSiteVector;
//typedef std::vector<AttachedGlycan*> AttachedGlycanVector;

/*******************************************/
/* Function Declarations                   */
/*******************************************/

void GetCenterOfGeometry(Assembly *assembly, GeometryTopology::Coordinate *center);
double GetDistanceToAtom(Atom *A, Atom *otherAtom);
double CalculateAtomicOverlaps(Assembly *assemblyA, Assembly *assemblyB);

void GetResidueRingCenter (Residue *residue, GeometryTopology::Coordinate *center);
void Find_and_Prepare_Protein_Residues_for_Glycosylation(GlycoSiteVector *glycosites, ResidueVector *protein_residues, const std::vector<std::__cxx11::string> *glycositeList);

//int main(int argc, char *argv[])
int main()
{
    std::string working_Directory = Find_Program_Working_Directory();
    std::string installation_Directory = Find_Program_Installation_Directory();

    //************************************************//
    // Reading input file                             //
    //************************************************//

    GlycoSiteVector glycoSites;
    std::string proteinPDB, glycanDirectory, buffer;
    std::vector<std::string> glycositeList, listOfGlycans;
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
        if(strInput == "Protein Residue list:") // Reads lines until it finds a line with "END"
        {
            getline(inf, buffer);
            while(buffer != "END")
            {
                glycositeList.push_back(buffer);
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
                glycoSites.push_back(new GlycosylationSite(buffer));
                std::cout << "Found glycan:" << buffer << std::endl;
                getline(inf, buffer);          
            }
        }
    }

    //************************************************//
    // Load Protein PDB file                          //
    //************************************************//

    Assembly glycoprotein ((working_Directory + "/inputs/" + proteinPDB), gmml::InputFileType::PDB);
    glycoprotein.BuildStructureByDistance();

    ResidueVector protein_residues = glycoprotein.GetResidues();
    Find_and_Prepare_Protein_Residues_for_Glycosylation(&glycoSites, &protein_residues, &glycositeList);

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
        for (GlycoSiteVector::iterator it = glycoSites.begin(); it != glycoSites.end(); ++it)
        {
            GlycosylationSite* glycosite = *it;
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

    resolve_overlaps::monte_carlo(glycoprotein, glycoSites);

    std::cout << "Program got to end ok" << std::endl;
    return 0;
}

void Find_and_Prepare_Protein_Residues_for_Glycosylation(GlycoSiteVector *glycosites, ResidueVector *protein_residues, const std::vector<std::string> *glycositeList)
{
    int i = 0;
    for(std::vector<std::string>::const_iterator it = glycositeList->begin(); it != glycositeList->end(); ++it)
    {
        std::string glycosite = (*it);
        std::cout << "Preparing glycosite: " << glycosite << std::endl;
        for(ResidueVector::iterator itt = protein_residues->begin(); itt != protein_residues->end(); ++itt)
        {
            // Comparing strings is easier:
            Residue *residue = (*itt);
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
                //prepare_Glycans_For_Superimposition_To_Particular_Residue(working_Directory, residue->GetName(), &glycan); //OG NOT HAVING glycan SET WILL BREAK IT LATER!!!!!!!!!!!!!!!!!!!!!!!!!
              /*  if (residue->GetName().compare("SER")==0)
                    residue->SetName("OLS");
                if (residue->GetName().compare("THR")==0)
                    residue->SetName("OLT");
                if (residue->GetName().compare("ASN")==0)
                    residue->SetName("NLN");
                if (residue->GetName().compare("TYR")==0)
                    residue->SetName("OLY");
                    */
            }
        }
    }
}

double GetDistanceToAtom(Atom *A, Atom *otherAtom)
{
    double x = ( A->GetCoordinates().at(0)->GetX() - otherAtom->GetCoordinates().at(0)->GetX() );
    double y = ( A->GetCoordinates().at(0)->GetY() - otherAtom->GetCoordinates().at(0)->GetY() );
    double z = ( A->GetCoordinates().at(0)->GetZ() - otherAtom->GetCoordinates().at(0)->GetZ() );

    return sqrt( (x*x) + (y*y) + (z*z) );
}

double CalculateAtomicOverlaps(Assembly *assemblyA, Assembly *assemblyB)
{
    AtomVector assemblyAAtoms = assemblyA->GetAllAtomsOfAssembly();
    AtomVector assemblyBAtoms = assemblyB->GetAllAtomsOfAssembly();

    double rA = 0.0, rB = 0.0, distance = 0.0, totalOverlap = 0.0;

    for(AtomVector::iterator atomA = assemblyAAtoms.begin(); atomA != assemblyAAtoms.end(); atomA++)
    {
        for(AtomVector::iterator atomB = assemblyBAtoms.begin(); atomB != assemblyBAtoms.end(); atomB++)
        {
            distance = GetDistanceToAtom((*atomA), (*atomB));
            if ( ( distance < 3.6 ) && ( distance > 0.0 ) ) //Close enough to overlap, but not the same atom
            {
                // element info not set, so I look at first letter of atom name.
                if ((*atomA)->GetName().at(0) == 'C') rA = 1.70;
                if ((*atomA)->GetName().at(0) == 'O') rA = 1.52;
                if ((*atomA)->GetName().at(0) == 'N') rA = 1.55;
                if ((*atomA)->GetName().at(0) == 'S') rA = 1.80;

                if ((*atomB)->GetName().at(0) == 'C') rB = 1.70;
                if ((*atomB)->GetName().at(0) == 'O') rB = 1.52;
                if ((*atomB)->GetName().at(0) == 'N') rB = 1.55;
                if ((*atomB)->GetName().at(0) == 'S') rB = 1.80;
       //         std::cout << "Distance=" << distance << " rA=" << rA << " rB=" << rB << std::endl;
                if (rA + rB > distance + 0.6){ // 0.6 overlap is deemed acceptable. (Copying chimera:)

                    totalOverlap += ( 2 * PI * rA* ( rA - distance / 2 - ( ( (rA*rA) - (rB*rB) ) / (2 * distance) ) ) );

                    //std::cout << "Overlap=" << totalOverlap << std::endl;
                }
            }
        }
    }
    return totalOverlap;
}

void GetCenterOfGeometry(Assembly *assembly, GeometryTopology::Coordinate *center)
{
    double sumX=0.0;
    double sumY=0.0;
    double sumZ=0.0;

    CoordinateVector all_coords = assembly->GetAllCoordinates();
    for(CoordinateVector::iterator it = all_coords.begin(); it != all_coords.end(); it++)
    {

        GeometryTopology::Coordinate coord = *it;
        sumX += coord.GetX();
        sumY += coord.GetY();
        sumZ += coord.GetZ();
    }

    center->SetX( (sumX / all_coords.size()) );
    center->SetY( (sumY / all_coords.size()) );
    center->SetZ( (sumZ / all_coords.size()) );

    return;
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

// Change Coordinate so I can do coordinate.subtract(coordinate);
GeometryTopology::Coordinate subtract_coordinates(GeometryTopology::Coordinate minuaend, GeometryTopology::Coordinate subtrahend)
{
    GeometryTopology::Coordinate new_coordinate ( (minuaend.GetX()-subtrahend.GetX()), (minuaend.GetY()-subtrahend.GetY()), (minuaend.GetZ()-subtrahend.GetZ()) );
    return new_coordinate;
}

void GetResidueRingCenter (Residue *residue, GeometryTopology::Coordinate *center)
{
    double sumX = 0.0, sumY = 0.0, sumZ = 0.0;
    int numberOfRingAtoms = 0;
    AtomVector atoms = residue->GetAtoms();

    for(Assembly::AtomVector::iterator atom = atoms.begin(); atom != atoms.end(); atom++)
    {

        if ( (*atom)->GetIsRing() )
        {
            numberOfRingAtoms++;
            std::cout << "Atom is ring: " << (*atom)->GetName() << std::endl;
            sumX += (*atom)->GetCoordinates().at(0)->GetX();
            sumY += (*atom)->GetCoordinates().at(0)->GetY();
            sumZ += (*atom)->GetCoordinates().at(0)->GetZ();
        }
    }
    center->SetX( sumX / numberOfRingAtoms  );
    center->SetY( sumY / numberOfRingAtoms  );
    center->SetZ( sumZ / numberOfRingAtoms  );
    return;
}

