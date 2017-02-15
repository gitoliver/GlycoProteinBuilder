#include <fstream>
#include <iostream>
#include <cstdlib> // for exit()
#include "glycoproteinbuilder.h"
#include <QApplication>
// Includes for directory reading
#include <string>
#include <dirent.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "/home/oliver/Programs/gems/gmml/includes/gmml.hpp"

#include <Eigen/Geometry>



constexpr auto PI = 3.14159265358979323846;

using namespace MolecularModeling;
typedef std::vector<Atom*> AtomVector;
typedef std::vector<GeometryTopology::Coordinate*> CoordinateVector;
typedef std::vector<Residue*> ResidueVector;
typedef std::vector<Assembly*> AssemblyVector;
typedef GeometryTopology::Coordinate Vector;

//void FindConnectedAtoms(Atom *atom, Assembly::AtomVector &visitedAtoms);

/*******************************************/
/* Function Declarations                   */
/*******************************************/

void GetCenterOfGeometry(Assembly *assembly, GeometryTopology::Coordinate *center);
double GetDistanceToAtom(Atom *A, Atom *otherAtom);
double CalculateAtomicOverlaps(Assembly *assemblyA, Assembly *assemblyB);
// Superimposition Functions
void TestFind3DAffineTransform();
Eigen::Affine3d Find3DAffineTransform(Eigen::Matrix3Xd in, Eigen::Matrix3Xd out);
void GenerateMatrixFromAssembyCoordinates(Assembly *assembly, Eigen::Matrix3Xd *matrix);
void GenerateMatrixFromAtomVectorCoordinates(AtomVector *atoms, Eigen::Matrix3Xd *matrix);
void ReplaceAssemblyCoordinatesFromMatrix(Assembly *assembly, Eigen::Matrix3Xd *matrix);
void Superimpose(AtomVector *moving, AtomVector *target);
void Superimpose(Assembly *moving, Assembly *target);
void Superimpose(Assembly *moving, Assembly *target, Assembly *alsoMoving);
void Superimpose(Assembly *moving, Assembly *target, AssemblyVector *alsoMoving);
void GetResidueRingCenter (Residue *residue, GeometryTopology::Coordinate *center);
GeometryTopology::Coordinate subtract_coordinates(GeometryTopology::Coordinate minuaend, GeometryTopology::Coordinate subtrahend);
GeometryTopology::Coordinate get_cartesian_point_from_internal_coords(GeometryTopology::Coordinate a, GeometryTopology::Coordinate b, GeometryTopology::Coordinate c, double theta_Degrees, double phi_Degrees, double distance_Angstrom);
GeometryTopology::Coordinate get_cartesian_point_from_internal_coords(Atom *a, Atom *b, Atom *c, double theta_Degrees, double phi_Degrees, double distance_Angstrom);
void prepare_Glycans_For_Superimposition_To_Any_Residue(std::string working_Directory, std::string amino_acid_name, Assembly *glycan);


int main(int argc, char *argv[])
{
    QApplication abss(argc, argv);
    glycoproteinBuilder w;
    w.show();

    std::string working_Directory = "/home/oliver/work/zq.Sarah_Flowers_OLink_GlycoProteinBuilder/3.GMML_GlycoProteinBuilder/";
    std::string installtionDirectory;

    //************************************************//
    // Details for loading in a PDB file              //
    //************************************************//

    std::vector<std::string> amino_libs, glycam_libs, other_libs, prep;
    amino_libs.push_back("/home/oliver/Programs/gems/gmml/dat/CurrentParams/leaprc.ff12SB_2014-04-24/amino12.lib");
    amino_libs.push_back("/home/oliver/Programs/gems/gmml/dat/CurrentParams/leaprc.ff12SB_2014-04-24/aminoct12.lib");
    amino_libs.push_back("/home/oliver/Programs/gems/gmml/dat/CurrentParams/leaprc.ff12SB_2014-04-24/aminont12.lib");

    glycam_libs.push_back("/home/oliver/Programs/gems/gmml/dat/CurrentParams/leaprc_GLYCAM_06j-1_2014-03-14/GLYCAM_amino_06j_12SB.lib");
    glycam_libs.push_back("/home/oliver/Programs/gems/gmml/dat/CurrentParams/leaprc_GLYCAM_06j-1_2014-03-14/GLYCAM_aminoct_06j_12SB.lib");
    glycam_libs.push_back("/home/oliver/Programs/gems/gmml/dat/CurrentParams/leaprc_GLYCAM_06j-1_2014-03-14/GLYCAM_aminont_06j_12SB.lib");

    other_libs.push_back("/home/oliver/Programs/gems/gmml/dat/CurrentParams/leaprc.ff12SB_2014-04-24/nucleic12.lib");
    other_libs.push_back("/home/oliver/Programs/gems/gmml/dat/CurrentParams/leaprc.ff12SB_2014-04-24/nucleic12.lib");
    other_libs.push_back("/home/oliver/Programs/gems/gmml/dat/CurrentParams/other/solvents.lib");

    prep.push_back("/home/oliver/Programs/gems/gmml/dat/CurrentParams/leaprc_GLYCAM_06j-1_2014-03-14/GLYCAM_06j-1.prep");

    std::string pdb_file_path = working_Directory + "test.pdb";
    std::string parameter_file_path = "/home/oliver/Programs/gems/gmml/dat/CurrentParams/leaprc_GLYCAM_06j-1_2014-03-14/GLYCAM_06j.dat";
    std::string ion_parameter_file_path = "/home/oliver/Programs/gems/gmml/dat/CurrentParams/other/atomic_ions.lib";

    //************************************************//
    // Build PDB file from condensed sequence         //
    //************************************************//
/*
    Assembly assemblyC;
    assemblyC.BuildAssemblyFromCondensedSequence("DGlcpb1-4DNeup5Aca2-4DGlcpb1-OME", prep.at(0) ,parameter_file_path, true);
    assemblyC.BuildStructureByDistance();
    PdbFileSpace::PdbFile *outputPdbFile1 = assemblyC.BuildPdbFileStructureFromAssembly();
    outputPdbFile1->Write(working_Directory + "test-buildFromStructure.pdb");

/*
    TopologyFileSpace::TopologyFile *outputTopologyFile1 = assemblyC.BuildTopologyFileStructureFromAssembly(parameter_file_path, ion_parameter_file_path);
    outputTopologyFile1->Write(working_Directory + "test-BuildFromSequence.parm7");

    CoordinateFileSpace::CoordinateFile *outputCoordinateFile = assemblyC.BuildCoordinateFileStructureFromAssembly();
    outputCoordinateFile->Write(working_Directory + "test-BuildFromSequence.rst7");
*/

    //Assembly assemblyD;
    //assemblyD.BuildAssemblyFromPdbFile(pdb_file_path, amino_libs, glycam_libs, other_libs, prep, parameter_file_path );
    //assemblyD.BuildStructureByDistance();

    //************************************************//
    // Reading an input file                          //
    //************************************************//

    std::string proteinPDB, glycanDirectory, buffer;
    std::vector<std::string> glycositeList, listOfGlycans;

    std::ifstream inf ("/home/oliver/work/zq.Sarah_Flowers_OLink_GlycoProteinBuilder/3.GMML_GlycoProteinBuilder/input.txt");
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
        if(strInput == "Protein Residue list:")
        {
            getline(inf, buffer);
            while(buffer != "END")
            {
                glycositeList.push_back(buffer);
                std::cout << "residue: " << buffer << std::endl;
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
                std::cout << "glycan: " << buffer << std::endl;
                getline(inf, buffer);
            }
        }
        //std::cout << strInput << std::endl;
    }
    //std::cout << proteinPDB << ", " << glycanDirectory << std::endl;



    //************************************************//
    // Create list containing pairs:                  //
    // protein residue - Glycan                       //
    //************************************************//

    /*
     std::list < std::pair<std::string, std::string> > protein_glycanList;
     protein_glycanList.push_back( std::pair<std::string, std::string>(glycositeList.at(i), listOfGlycans.at(i)) );
     */

    //************************************************//
    // Load Files From a Directory                    //
    //************************************************//

    std::ifstream fin;
    std::string directory = "/home/oliver/work/zq.Sarah_Flowers_OLink_GlycoProteinBuilder/3.GMML_GlycoProteinBuilder/" + glycanDirectory ;
    std::cout << "directory: " << directory << std::endl;
    std::string filepath;
    DIR *dp; // A directory stream
    struct dirent *dirp; // Contains file serial number and name (char d_name[])
    struct stat filestat; // Contains info about file, such as device ID, user ID, access time etc
    int num;

    dp = opendir( directory.c_str() ); //.c_str adds a null character to the end.
    if (dp == NULL)
    {
        std::cout << "Error(" << errno << ") opening " << directory << std::endl;
        return errno;
    }
    while ((dirp = readdir ( dp )))
    {
        // std::cout << "Here" << std::endl;
        filepath = directory + "/" + dirp->d_name;

        // If the file is a directory (or is in some way invalid) we'll skip it
        if (stat( filepath.c_str(), &filestat )) continue; // Is it a valid file?
        if (S_ISDIR( filestat.st_mode ))         continue; // Is it a directory?

        std::string temp = listOfGlycans.at(0);
        if ( temp.compare(0, temp.size(), dirp->d_name, 0, temp.size()) == 0 )
        {   // Does glycan filename start with appropriate code?
            std::cout << "BOOM " << temp << " " << dirp->d_name << std::endl;
           // assembly.BuildAssemblyFromPdbFile(filepath, amino_libs, glycam_libs, other_libs, prep, parameter_file_path);
           // assembly.BuildStructureByDistance();
        }
       // Endeavor to read a single number from the file and display it
       // fin.open( filepath.c_str() );
       // if (fin >> num)
       //     std::cout << filepath << ": " << num << std::endl;
       // fin.close();

    }
    closedir( dp );

    //************************************************//
    // Building the necessary N-link Torsions         //
    //************************************************//


    //Yohanna example
/*
    Assembly yohAssembly;
    pdb_file_path = "/home/oliver/Downloads/1hzh.pdb";
    yohAssembly.BuildAssemblyFromPdbFile(pdb_file_path, amino_libs, glycam_libs, other_libs, prep, parameter_file_path);
    yohAssembly.BuildStructureByDistance();

    ResidueVector yohResidues = yohAssembly.GetAllResiduesOfAssembly();

    GeometryTopology::Coordinate *centerOfRing = new GeometryTopology::Coordinate();
    for(ResidueVector::iterator it = yohResidues.begin(); it != yohResidues.end(); ++it)
    {
        Residue* current_residue = *it;
        if(current_residue->GetName().compare("TYR")==0 || current_residue->GetName().compare("PHE")==0 ||
           current_residue->GetName().compare("TRP")==0 || current_residue->GetName().compare("HIS")==0)

        {
            // Want it to look like this:
           //centerOfRing = current_residue->GetResidueRingCenter;

           GetResidueRingCenter(current_residue, centerOfRing);
           std::cout << "Residue is " << current_residue->GetId() << std::endl;
           centerOfRing->Print();
           std::cout << std::endl;
        }
    }
    */

    Assembly glycan;
    //pdb_file_path = working_Directory + listOfGlycans.at(0);
    pdb_file_path = working_Directory + "/glycans/SiaCore1_-g.pdb";
    glycan.BuildAssemblyFromPdbFile(pdb_file_path, amino_libs, glycam_libs, other_libs, prep, parameter_file_path);
    glycan.BuildStructureByDistance();



    //prepare_Glycans_For_Superimposition_To_Any_Residue(working_Directory, "THR", &glycan);
   // prepare_Glycans_For_Superimposition_To_Any_Residue(working_Directory, "ASN", &glycan);

    //************************************************//
    // COM translating                                //
    //************************************************//
    /*
    GeometryTopology::Coordinate center_of_geometry;
    GetCenterOfGeometry(&assembly, &center_of_geometry);

    center_of_geometry.Print();

    CoordinateVector assemblyCoordinates = assembly.GetAllCoordinates();
    for(CoordinateVector::iterator it = assemblyCoordinates.begin(); it != assemblyCoordinates.end(); it++)
    {
        (*it)->Translate( -(center_of_geometry.GetX()), -(center_of_geometry.GetY()), -(center_of_geometry.GetZ()) );
    }

    PdbFileSpace::PdbFile *outputPdbFile = assembly.BuildPdbFileStructureFromAssembly();
    outputPdbFile->Write(working_Directory + "test-Translate.pdb");
    */
    //************************************************//
    // Adjusting a torsion                            //
    //************************************************//

    /*
    Atom *atom1, *atom2, *atom3, *atom4;

    Assembly::AtomVector assAtoms = assembly.GetAllAtomsOfAssembly();

    Residue* residue;
    std::string residueId, residueName, atomName;

    for (auto &atom : assAtoms){
        atomName = atom->GetName();
        residue = atom->GetResidue();
        residueId = residue->GetId();
        if ((atomName.compare("O5") == 0) && (residueId.compare(6,1,"4") == 0)) {
            atom1 = atom;
            std::cout << residue->GetId() << std::endl;
        //    std::cout << atomName << std::endl;
        }
        if ((atomName.compare("C1") == 0) && (residueId.compare(6,1,"4") == 0))
            atom2 = atom;
        if ((atomName.compare("O4") == 0) && (residueId.compare(6,1,"3") == 0))
            atom3 = atom;
        if ((atomName.compare("C4") == 0) && (residueId.compare(6,1,"3") == 0))
            atom4 = atom;
    }
    std::cout << atom1->GetName() << ", " << atom2->GetName() << ", " << atom3->GetName() << ", " << atom4->GetName() << std::endl;

    double torsion1 = 80.0;
    assembly.SetDihedral(atom1, atom2, atom3, atom4, torsion1);

    outputPdbFile->Write(working_Directory + "test-Dihedral.pdb");
    */
    //************************************************//
    // Atomic Overlaps                                //
    //************************************************//

    //std::cout << "Distance is " << GetDistanceToAtom(atom2, atom4) << std::endl;


    // Write a PDB file containing the moved co-ordinates
   // PdbFileSpace::PdbFile *outputPdbFileA = assemblyMoving.BuildPdbFileStructureFromAssembly();
  //  outputPdbFileA->Write(working_Directory + "Superimposed-moved.pdb");


  //  std::cout << "Overlap is " << CalculateAtomicOverlaps(&assemblyA, &assemblyB) << std::endl;

    //assemblyA.AddAssembly(&assemblyB);
   // outputPdbFile->Write(working_Directory + "test-MergeAssemblies.pdb");

    //************************************************//
    // Superimposition                                //
    //************************************************//
    std::cout << "Superimposition" << std::endl;

    Assembly protein;
    protein.BuildAssemblyFromPdbFile( (working_Directory + proteinPDB), amino_libs, glycam_libs, other_libs, prep, parameter_file_path );
    protein.BuildStructureByDistance();

    ResidueVector *glycosites = new ResidueVector();

    // Get list of Residues that we will superimpose the glycans onto
    for(std::vector<std::string>::iterator it = glycositeList.begin(); it < glycositeList.end(); ++it)
    {
        std::string glycosite=(*it);
        ResidueVector residues = protein.GetResidues();
        for(ResidueVector::iterator itt = residues.begin(); itt != residues.end(); ++itt)
        {
            // Comparing strings is easier:
            Residue *residue = (*itt);
            std::string id = residue->GetId();
            std::string formatted_glycosite = "_" + glycosite + "_";
            // If the residue number set in the input file is equal to the current residue number
            if( id.compare(5,formatted_glycosite.size(),formatted_glycosite) == 0)
            {
                std::cout << "glycosite: " << glycosite << std::endl;
                std::cout << "glycosite id:" << id << std::endl;
                glycosites->push_back(residue);
                prepare_Glycans_For_Superimposition_To_Any_Residue(working_Directory, residue->GetName(), &glycan); //OG NOT HAVING glycan SET WILL BREAK IT LATER!!!!!!!!!!!!!!!!!!!!!!!!!
                if (residue->GetName().compare("SER")==0)
                    residue->SetName("OLS");
                if (residue->GetName().compare("THR")==0)
                    residue->SetName("OLT");
                if (residue->GetName().compare("ASN")==0)
                    residue->SetName("NLN");
                if (residue->GetName().compare("TYR")==0)
                    residue->SetName("OLY");
            }
        }
    }
    // Want: AtomVector = assembly.Selection(@CA,N,ND2,CG,CB)
    // Superimpose onto each glycosite
   // std::vector<Assembly> *addedGlycans;
    Assembly *addedGlycans = new Assembly[glycosites->size()];
    int i=0;

    for(ResidueVector::iterator it = glycosites->begin(); it != glycosites->end(); ++it)
    {
        Residue *glycosite=(*it);
        std::cout << "Glycosite is " << glycosite->GetId() << std::endl;
        AtomVector atoms = glycosite->GetAtoms();
        Assembly* assemblyTarget = new Assembly();

        Residue* residueTarget = new Residue();
        for(AtomVector::iterator itt = atoms.begin(); itt != atoms.end(); ++itt)
        {
            Atom *atom = (*itt);
            std::cout << "Checking: " << atom->GetName() << std::endl;
            if (atom->GetName() == "CB")
                residueTarget->AddAtom(atom);
            if (atom->GetName() == "CA")
                residueTarget->AddAtom(atom);
            if (atom->GetName() == "N")
                residueTarget->AddAtom(atom);
        }
        residueTarget->SetAssembly(assemblyTarget);
        assemblyTarget->AddResidue(residueTarget);
        assemblyTarget->BuildStructureByDistance();

        pdb_file_path = working_Directory + "/glycans/SiaCore1_-g.pdb";
        addedGlycans[i].BuildAssemblyFromPdbFile(pdb_file_path, amino_libs, glycam_libs, other_libs, prep, parameter_file_path);
        addedGlycans[i].BuildStructureByDistance();

        pdb_file_path = working_Directory + glycosite->GetName() + "_AlignedToGlycan.pdb";
        std::cout << "Attempting to open " << pdb_file_path << std::endl;
        Assembly *residueAminoAcid = new Assembly();
        residueAminoAcid->BuildAssemblyFromPdbFile(pdb_file_path, amino_libs, glycam_libs, other_libs, prep, parameter_file_path);
        residueAminoAcid->BuildStructureByDistance();

        AssemblyVector *atomsToAlsoMove = new AssemblyVector; // using AddAssembly for the NLN part crashes the PDB writer.
        atomsToAlsoMove->push_back((&addedGlycans[i]));
        atomsToAlsoMove->push_back(residueAminoAcid);

        pdb_file_path = working_Directory + glycosite->GetName() + "_superimposition_atoms.pdb";
        std::cout << "Attempting to open " << pdb_file_path << std::endl;

        Assembly *assemblyMoving = new Assembly();
        assemblyMoving->BuildAssemblyFromPdbFile(pdb_file_path, amino_libs, glycam_libs, other_libs, prep, parameter_file_path );

        //Superimpose(assemblyMoving, assemblyTarget, &addedGlycans[i]);
        Superimpose(assemblyMoving, assemblyTarget, atomsToAlsoMove);

        // Want: atom.SetX()
        // Want: atom.ReplaceCoords()
        // Replace Co-ords of sidechain
        // Want: assembly.Select(:residue@atom)
        AtomVector newSideChainAtoms = residueAminoAcid->GetAllAtomsOfAssembly();
        for(AtomVector::iterator itt = atoms.begin(); itt != atoms.end(); ++itt)
        {
            Atom *atom = (*itt);
            std::cout << "Checking: " << atom->GetName() << std::endl;
            if ( (atom->GetName() != "N") && (atom->GetName() != "CA") && (atom->GetName() != "CB") )
            {
                std::cout << "Replacing co-ords of protein " << atom->GetName() << std::endl;
                for(AtomVector::iterator ittt = newSideChainAtoms.begin(); ittt != newSideChainAtoms.end(); ++ittt)
                {
                    Atom *atom1 = (*ittt);
                    std::cout << "Comparing with " << atom1->GetName() << std::endl;
                    if (atom->GetName() == atom1->GetName() )
                    {

                        //    std::cout << "Before X=" << atom->GetCoordinates().at(0)->GetX() << std::endl;
                            atom->GetCoordinates().at(0)->SetX( atom1->GetCoordinates().at(0)->GetX() );
                            atom->GetCoordinates().at(0)->SetY( atom1->GetCoordinates().at(0)->GetY() );
                            atom->GetCoordinates().at(0)->SetZ( atom1->GetCoordinates().at(0)->GetZ() );
                         //   std::cout << "After X=" << atom->GetCoordinates().at(0)->GetX() << std::endl;
                    }
                }
            }
        }

        // Write out each glycan to check superimposition
        std::stringstream ss;
        ss << working_Directory + "GlycoProtein" << i << ".pdb";
        PdbFileSpace::PdbFile *outputPdbFileGlycoProtein = addedGlycans[i].BuildPdbFileStructureFromAssembly(-1,0);
        outputPdbFileGlycoProtein->Write(ss.str());
        protein.AddAssembly(&addedGlycans[i]);
        ++i;
    }
    std::cout << "Got this far" << std::endl;
    // Add superimposed glycan (assemblyAlsoMoving) to protein assembly
    i = 1;

    /*
    for(std::vector<Assembly>::iterator addedGlycan = addedGlycans->begin(); addedGlycan != addedGlycans->end(); addedGlycan++)
    {
        std::stringstream ss;
        ss << working_Directory + "GlycoProtein" << i << ".pdb";
        i++;
        PdbFileSpace::PdbFile *outputPdbFileGlycoProtein = (*addedGlycan).BuildPdbFileStructureFromAssembly(-1,0);
        outputPdbFileGlycoProtein->Write(ss.str());
        protein.AddAssembly(&(*addedGlycan));
    }
    */
    //protein.Print();

    std::cout << "Got this far" << std::endl;
    PdbFileSpace::PdbFile *outputPdbFileGlycoProteinAll = protein.BuildPdbFileStructureFromAssembly(-1,0);
    outputPdbFileGlycoProteinAll->Write(working_Directory + "GlycoProtein.pdb");
    delete[] addedGlycans;

      //************************************************//
      // Selections                                     //
      //************************************************//

      //std::vector<int> vec;
      // .. put in some values ..
     // int int_to_remove = n;
    //  vec.erase(std::remove(vec.begin(), vec.end(), int_to_remove), vec.end());


     /* AtomVector selection = assemblyA.GetAllAtomsOfAssembly();

      for (auto &atom : assemblyA.GetAllAtomsOfAssembly()){
          atomName = atom->GetName();
          if (atomName.compare("O5") == 0)  {
            //  selection.erase(std::remove(selection.begin(), selection.end(), atom), selection.end());
              assemblyA.GetAllAtomsOfAssembly().erase(std::remove(assemblyA.GetAllAtomsOfAssembly().begin(), assemblyA.GetAllAtomsOfAssembly().end(), atom), assemblyA.GetAllAtomsOfAssembly().end());
          }
       }
       */

   //   for (auto &atom : selection)
    //      std::cout << "Atom: " << atom->GetName() << std::endl;


  //     AtomVector selection2 = assemblyMoving.GetAllAtomsOfAssembly();
  //     for (auto &atom : selection2){
    //       atomName = atom->GetName();
   //        if (atomName.compare("O5") == 0)  {
   //            std::cout << "Boo found again" << std::endl;
   //        }
 //     }
 //      outputPdbFileA->Write(working_Directory + "test-Delete-Function.pdb");

    std::cout << "Program got to end ok" << std::endl;
    return abss.exec();
}


void GenerateMatrixFromAssembyCoordinates(Assembly *assembly, Eigen::Matrix3Xd *matrix)
{
    int col = 0; // Column index for matrix
    CoordinateVector assemblyCoordinates = assembly->GetAllCoordinates();
    for(CoordinateVector::iterator it = assemblyCoordinates.begin(); it != assemblyCoordinates.end(); it++)
    {
        (*matrix)(0, col) = (*it)->GetX();
        (*matrix)(1, col) = (*it)->GetY();
        (*matrix)(2, col) = (*it)->GetZ();
        col++;
    }
}

void ReplaceAssemblyCoordinatesFromMatrix(Assembly *assembly, Eigen::Matrix3Xd *matrix)
{
    int col = 0;
    CoordinateVector assemblyCoordinates = assembly->GetAllCoordinates();
    for(CoordinateVector::iterator it = assemblyCoordinates.begin(); it != assemblyCoordinates.end(); it++)
    {
        (*it)->SetX( (*matrix)(0, col) );
        (*it)->SetY( (*matrix)(1, col) );
        (*it)->SetZ( (*matrix)(2, col) );
        col++;
    }
}

//Atom Vector version, may be removed
void GenerateMatrixFromAtomVectorCoordinates(AtomVector *atoms, Eigen::Matrix3Xd *matrix)
{
    int col = 0; // Column index for matrix
    for(AtomVector::iterator it = atoms->begin(); it != atoms->end(); it++)
    {
        (*matrix)(0, col) = (*it)->GetCoordinates().at(0)->GetX();
        (*matrix)(1, col) = (*it)->GetCoordinates().at(0)->GetY();
        (*matrix)(2, col) = (*it)->GetCoordinates().at(0)->GetZ();
        col++;
    }
}

void ReplaceAtomVectorCoordinatesFromMatrix(AtomVector *atoms, Eigen::Matrix3Xd *matrix)
{
    int col = 0; // Column index for matrix
    for(AtomVector::iterator it = atoms->begin(); it != atoms->end(); it++)
    {
        (*it)->GetCoordinates().at(0)->SetX( (*matrix)(0, col) );
        (*it)->GetCoordinates().at(0)->SetY( (*matrix)(1, col) );
        (*it)->GetCoordinates().at(0)->SetZ( (*matrix)(2, col) );
        col++;
    }
}
//End atom vector version


void Superimpose(AtomVector *moving, AtomVector *target)
{
    Eigen::Matrix3Xd movingMatrix(3, moving->size()), targetMatrix(3, target->size()) ;

    //Create Matrices containing co-ordinates of moving and target
    GenerateMatrixFromAtomVectorCoordinates(moving, &movingMatrix);
    GenerateMatrixFromAtomVectorCoordinates(target, &targetMatrix);

    // Figure out how to move assembly moving onto target
    Eigen::Affine3d Affine = Find3DAffineTransform(movingMatrix, targetMatrix);

    // Create a matirx containing the moved co-ordinates of assembly moving.
    Eigen::Matrix3Xd movedMatrix = (Affine * movingMatrix);

    // Replace co-ordinates of moving with moved co-ordinates
    ReplaceAtomVectorCoordinatesFromMatrix(moving, &movedMatrix);
}

void Superimpose(Assembly *moving, Assembly *target)
{
    // Was thinking to do everything via AtomVector, but probably won't
    //AtomVector movingAtoms = moving->GetAllAtomsOfAssembly();
    //AtomVector targetAtoms = target->GetAllAtomsOfAssembly();
    //Superimpose(&movingAtoms, &targetAtoms);

    Eigen::Matrix3Xd movingMatrix(3, moving->GetAllCoordinates().size()), targetMatrix(3, target->GetAllCoordinates().size());

    // Create a matrices containing co-ordinates of assembly moving and target
    GenerateMatrixFromAssembyCoordinates(moving, &movingMatrix);
    GenerateMatrixFromAssembyCoordinates(target, &targetMatrix);

    // Figure out how to move assembly moving onto target
    Eigen::Affine3d Affine = Find3DAffineTransform(movingMatrix, targetMatrix);

    // Create a matirx containing the moved co-ordinates of assembly moving.
    Eigen::Matrix3Xd movedMatrix = (Affine * movingMatrix);

    // Replace co-ordinates of moving with moved co-ordinates
    ReplaceAssemblyCoordinatesFromMatrix(moving, &movedMatrix);

    // Write a PDB file containing the moved co-ordinates
   // PdbFileSpace::PdbFile *outputPdbFile = moving->BuildPdbFileStructureFromAssembly();
   // outputPdbFile->Write("/home/oliver/Desktop/test-Superimposed-Function.pdb");
}

void Superimpose(Assembly *moving, Assembly *target, Assembly *alsoMoving)
{
    // Was thinking to do everything via AtomVector, but probably won't
    //AtomVector movingAtoms = moving->GetAllAtomsOfAssembly();
    //AtomVector targetAtoms = target->GetAllAtomsOfAssembly();
    //Superimpose(&movingAtoms, &targetAtoms);

    Eigen::Matrix3Xd movingMatrix(3, moving->GetAllCoordinates().size()), targetMatrix(3, target->GetAllCoordinates().size());
    Eigen::Matrix3Xd alsoMovingMatrix(3, alsoMoving->GetAllCoordinates().size()); // separate from above line for clarity

    // Create a matrices containing co-ordinates of assembly moving and target
    GenerateMatrixFromAssembyCoordinates(moving, &movingMatrix);
    GenerateMatrixFromAssembyCoordinates(target, &targetMatrix);
    GenerateMatrixFromAssembyCoordinates(alsoMoving, &alsoMovingMatrix);

    // Figure out how to move assembly moving onto target
    Eigen::Affine3d Affine = Find3DAffineTransform(movingMatrix, targetMatrix);

    // Create a matrix containing the moved co-ordinates of assembly moving and alsoMoving.
    Eigen::Matrix3Xd movedMatrix = (Affine * movingMatrix);
    Eigen::Matrix3Xd alsoMovedMatrix = (Affine * alsoMovingMatrix);


    // Replace co-ordinates of moving with moved co-ordinates
    ReplaceAssemblyCoordinatesFromMatrix(moving, &movedMatrix);
    ReplaceAssemblyCoordinatesFromMatrix(alsoMoving, &alsoMovedMatrix);


    // Write a PDB file containing the moved co-ordinates
    PdbFileSpace::PdbFile *outputPdbFile = moving->BuildPdbFileStructureFromAssembly();
    outputPdbFile->Write("/home/oliver/Desktop/test-Superimposed-Function-moved.pdb");
    PdbFileSpace::PdbFile *outputPdbFile1 = alsoMoving->BuildPdbFileStructureFromAssembly();
    outputPdbFile1->Write("/home/oliver/Desktop/test-Superimposed-Function-AlsoMoved.pdb");
}

void Superimpose(Assembly *moving, Assembly *target, AssemblyVector *alsoMoving)
{
    // Was thinking to do everything via AtomVector, but probably won't
    //AtomVector movingAtoms = moving->GetAllAtomsOfAssembly();
    //AtomVector targetAtoms = target->GetAllAtomsOfAssembly();
    //Superimpose(&movingAtoms, &targetAtoms);

    Eigen::Matrix3Xd movingMatrix(3, moving->GetAllCoordinates().size()), targetMatrix(3, target->GetAllCoordinates().size());

    // Create a matrices containing co-ordinates of assembly moving and target
    GenerateMatrixFromAssembyCoordinates(moving, &movingMatrix);
    GenerateMatrixFromAssembyCoordinates(target, &targetMatrix);


    // Figure out how to move assembly moving onto target
    Eigen::Affine3d Affine = Find3DAffineTransform(movingMatrix, targetMatrix);

    // Create a matrix containing the moved co-ordinates of assembly moving and alsoMoving.
    Eigen::Matrix3Xd movedMatrix = (Affine * movingMatrix);

    // Replace co-ordinates of moving with moved co-ordinates
    ReplaceAssemblyCoordinatesFromMatrix(moving, &movedMatrix);

    // Write a PDB file containing the moved co-ordinates
    PdbFileSpace::PdbFile *outputPdbFile = moving->BuildPdbFileStructureFromAssembly();
    outputPdbFile->Write("/home/oliver/Desktop/output-Superimposed-Function-moved.pdb");

    // Also move every assembly in also moving
    int i=1; // To create unique output PDB file names
    for(AssemblyVector::iterator it = alsoMoving->begin(); it != alsoMoving->end(); it++)
    {
        Assembly *assembly = (*it); // I hate (*it)
        Eigen::Matrix3Xd alsoMovingMatrix(3, assembly->GetAllCoordinates().size());
        GenerateMatrixFromAssembyCoordinates(assembly, &alsoMovingMatrix);
        Eigen::Matrix3Xd alsoMovedMatrix = (Affine * alsoMovingMatrix);
        ReplaceAssemblyCoordinatesFromMatrix(assembly, &alsoMovedMatrix);
        PdbFileSpace::PdbFile *outputPdbFile1 = assembly->BuildPdbFileStructureFromAssembly();
        std::stringstream ss;
        ss << "/home/oliver/Desktop/output-Superimposed-Function-AlsoMoved" << i << ".pdb";
        outputPdbFile1->Write(ss.str());
        ++i;
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

Eigen::Affine3d Find3DAffineTransform(Eigen::Matrix3Xd in, Eigen::Matrix3Xd out) {

  // Default output
  Eigen::Affine3d A;
  A.linear() = Eigen::Matrix3d::Identity(3, 3);
  A.translation() = Eigen::Vector3d::Zero();

  if (in.cols() != out.cols())
    throw "Find3DAffineTransform(): input data mis-match";

  // First find the scale, by finding the ratio of sums of some distances,
  // then bring the datasets to the same scale.
  double dist_in = 0, dist_out = 0;
  for (int col = 0; col < in.cols()-1; col++) {
    dist_in  += (in.col(col+1) - in.col(col)).norm();
    dist_out += (out.col(col+1) - out.col(col)).norm();
  }
  if (dist_in <= 0 || dist_out <= 0)
    return A;
  double scale = dist_out/dist_in;
  out /= scale;

  // Find the centroids then shift to the origin
  Eigen::Vector3d in_ctr = Eigen::Vector3d::Zero();
  Eigen::Vector3d out_ctr = Eigen::Vector3d::Zero();
  for (int col = 0; col < in.cols(); col++) {
    in_ctr  += in.col(col);
    out_ctr += out.col(col);
  }
  in_ctr /= in.cols();
  out_ctr /= out.cols();
  for (int col = 0; col < in.cols(); col++) {
    in.col(col)  -= in_ctr;
    out.col(col) -= out_ctr;
  }

  // SVD
  Eigen::MatrixXd Cov = in * out.transpose();
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(Cov, Eigen::ComputeThinU | Eigen::ComputeThinV);

  // Find the rotation
  double d = (svd.matrixV() * svd.matrixU().transpose()).determinant();
  if (d > 0)
    d = 1.0;
  else
    d = -1.0;
  Eigen::Matrix3d I = Eigen::Matrix3d::Identity(3, 3);
  I(2, 2) = d;
  Eigen::Matrix3d R = svd.matrixV() * I * svd.matrixU().transpose();

  // The final transform
  A.linear() = scale * R;
  A.translation() = scale*(out_ctr - R*in_ctr);

  return A;
}

// A function to test Find3DAffineTransform()

void TestFind3DAffineTransform(){

  // Create datasets with known transform
  Eigen::Matrix3Xd in(3, 100), out(3, 100);
  Eigen::Quaternion<double> Q(1, 3, 5, 2);
  Q.normalize();
  Eigen::Matrix3d R = Q.toRotationMatrix();
  double scale = 2.0;
  for (int row = 0; row < in.rows(); row++) {
    for (int col = 0; col < in.cols(); col++) {
      in(row, col) = log(2*row + 10.0)/sqrt(1.0*col + 4.0) + sqrt(col*1.0)/(row + 1.0);
    }
  }
  Eigen::Vector3d S;
  S << -5, 6, -27;
  for (int col = 0; col < in.cols(); col++)
    out.col(col) = scale*R*in.col(col) + S;

  Eigen::Affine3d A = Find3DAffineTransform(in, out);

  // See if we got the transform we expected
  if ( (scale*R-A.linear()).cwiseAbs().maxCoeff() > 1e-13 ||
       (S-A.translation()).cwiseAbs().maxCoeff() > 1e-13)
    throw "Could not determine the affine transform accurately enough";
}

// Change Coordinate so I can do coordinate.subtract(coordinate);
GeometryTopology::Coordinate subtract_coordinates(GeometryTopology::Coordinate minuaend, GeometryTopology::Coordinate subtrahend)
{
    GeometryTopology::Coordinate new_coordinate ( (minuaend.GetX()-subtrahend.GetX()), (minuaend.GetY()-subtrahend.GetY()), (minuaend.GetZ()-subtrahend.GetZ()) );
    return new_coordinate;
}

GeometryTopology::Coordinate get_cartesian_point_from_internal_coords(GeometryTopology::Coordinate a, GeometryTopology::Coordinate b, GeometryTopology::Coordinate c, double theta_Degrees, double phi_Degrees, double distance_Angstrom)
{
    //Convert from Degrees to Radians
    theta_Degrees = ( (theta_Degrees * PI) / 180 );
    phi_Degrees = ( (phi_Degrees * PI) / 180 );

    // theta is the angle between 3 atoms. Phi is the torsion between 4 atoms.
    Vector lmn_x, lmn_y, lmn_z;
    double x_p, y_p, z_p;

    Vector cb = subtract_coordinates(b, c);
    Vector ba = subtract_coordinates(a, b);

    //cb.Print(); std::cout << "^cb" << std::endl;
    //ba.Print(); std::cout << "^ba" << std::endl;
    //std::cout << std::endl;

    lmn_y = ba;
    lmn_y.CrossProduct(cb);
    lmn_y.Normalize();

    lmn_z = cb;
    lmn_z.Normalize();

    lmn_x = lmn_z;
    lmn_x.CrossProduct(lmn_y);

    /*
    lmn_x.Print();
    std::cout << "^lmn_x" << std::endl;
    lmn_y.Print();
    std::cout << "^lmn_y" << std::endl;
    lmn_z.Print();
    std::cout << "^lmn_z" << std::endl;
    */

    /*
    lmn_y = normalize_vec(get_crossprod(ba, cb));
    lmn_z = normalize_vec(cb);
    lmn_x = get_crossprod(lmn_y, lmn_z);
    */

    x_p = distance_Angstrom *  sin(theta_Degrees) * cos(phi_Degrees);
    y_p = distance_Angstrom * sin(theta_Degrees) * sin(phi_Degrees);
    z_p = distance_Angstrom * cos(theta_Degrees);

    //std::cout << "x_p=" << x_p << "y_p=" << y_p << "z_p=" << z_p << std::endl;

    return (GeometryTopology::Coordinate)
    {
        lmn_x.GetX()*x_p + lmn_y.GetX()*y_p + lmn_z.GetX()*z_p + c.GetX(),
        lmn_x.GetY()*x_p + lmn_y.GetY()*y_p + lmn_z.GetY()*z_p + c.GetY(),
        lmn_x.GetZ()*x_p + lmn_y.GetZ()*y_p + lmn_z.GetZ()*z_p + c.GetZ()
    };

}

GeometryTopology::Coordinate get_cartesian_point_from_internal_coords(Atom *a, Atom *b, Atom *c, double theta_Degrees, double phi_Degrees, double distance_Angstrom)
{
    GeometryTopology::Coordinate new_Coord = get_cartesian_point_from_internal_coords(a->GetCoordinates().at(0), b->GetCoordinates().at(0), c->GetCoordinates().at(0), theta_Degrees, phi_Degrees, distance_Angstrom);
    return {new_Coord.GetX(), new_Coord.GetY(), new_Coord.GetZ()};
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

void prepare_Glycans_For_Superimposition_To_Any_Residue(std::string working_Directory, std::string amino_acid_name, Assembly *glycan)
{

    //Dear future self, the order that you add the atoms to the residue matters for superimposition ie N, CA, CB , not CB, CA, N.

    Residue* reducing_Residue = glycan->GetAllResiduesOfAssembly().at(1);
    std::cout << "Reducing residue is " << reducing_Residue->GetName() << std::endl;
    AtomVector reducing_Atoms = reducing_Residue->GetAtoms();

    Atom* atomC5; // = reducing_Residue.FindAtom("C5");
    Atom* atomO5;
    Atom* atomC1;

    // Want: residue.FindAtom(string Name);
    for(AtomVector::iterator it = reducing_Atoms.begin(); it != reducing_Atoms.end(); it++)
    {
       Atom* atom = *it;
       if(atom->GetName().compare("C5")==0)
           atomC5 = atom;
       if(atom->GetName().compare("O5")==0)
           atomO5 = atom;
       if(atom->GetName().compare("C1")==0)
           atomC1 = atom;
    }

    // Want: atom.Coordinate (the first fucking coordinate).
    // Want: residue.DeleteAtom(int index);
    // Want: assembly.AddAtom(int index);
   // std::cout << "Atoms are: " << atomC5->GetName() << ", " << atomO5->GetName() << ", " << atomC1->GetName() << std::endl;
   // std::cout << std::endl;

    Assembly* assembly = new Assembly();
    Residue* residue = new Residue();
    residue->SetAssembly(assembly);
    assembly->AddResidue(residue);
    PdbFileSpace::PdbFile *outputPdbFile = new PdbFileSpace::PdbFile();

    if (amino_acid_name.compare("ASN")==0)
    {
        residue->SetName("NLN");
        residue->SetId("NLN_ _1_ _ _1");

        Atom *atomND2 = new Atom(residue, "ND2", (get_cartesian_point_from_internal_coords(atomC5, atomO5, atomC1, 109.3, 180, 1.53)));
        Atom *atomCG = new Atom(residue, "CG", (get_cartesian_point_from_internal_coords(atomO5, atomC1, atomND2, 109.3, 261, 1.325)));
        Atom *atomOD1 = new Atom(residue, "OD1", (get_cartesian_point_from_internal_coords(atomC1, atomND2, atomCG, 126, 0, 1.22)));
        Atom *atomCB = new Atom(residue, "CB", (get_cartesian_point_from_internal_coords(atomC1, atomND2, atomCG, 114, 177.3, 1.53)));
        Atom *atomCA = new Atom(residue, "CA", (get_cartesian_point_from_internal_coords(atomND2, atomCG, atomCB, 111, 177.6, 1.53)));
        Atom *atomN = new Atom(residue, "N", (get_cartesian_point_from_internal_coords(atomCG, atomCB, atomCA, 111, 191.6, 1.453)));

        residue->AddAtom(atomN);
        residue->AddAtom(atomCA);
        residue->AddAtom(atomCB);
        assembly->BuildStructureByDistance();
        outputPdbFile = assembly->BuildPdbFileStructureFromAssembly();
        outputPdbFile->Write(working_Directory + "NLN_superimposition_atoms.pdb");

        residue->AddAtom(atomCG);
        residue->AddAtom(atomOD1);
        residue->AddAtom(atomND2);
        assembly->BuildStructureByDistance();
        outputPdbFile = assembly->BuildPdbFileStructureFromAssembly();
        outputPdbFile->Write(working_Directory + "NLN_AlignedToGlycan.pdb");

    }

    else if (amino_acid_name.compare("THR")==0 || amino_acid_name.compare("SER")==0)
    {
        residue->SetName("OLS");
        residue->SetId("OLS_ _1_ _ _1");

        Atom *atomOG1 = new Atom(residue, "OG", (get_cartesian_point_from_internal_coords(atomC5, atomO5, atomC1, 112, 68, 1.46)));
        Atom *atomCB = new Atom(residue, "CB", (get_cartesian_point_from_internal_coords(atomO5, atomC1, atomOG1, 109.3, 75, 1.53)));
        Atom *atomCA = new Atom(residue, "CA", (get_cartesian_point_from_internal_coords(atomC1, atomOG1, atomCB, 109.3, 125, 1.53)));
        Atom *atomN = new Atom(residue, "N", (get_cartesian_point_from_internal_coords(atomOG1, atomCB, atomCA, 109.3, 180, 1.53)));
        Atom *atomCG2 = new Atom(residue, "CG2", (get_cartesian_point_from_internal_coords(atomC1, atomOG1, atomCB, 109.3, -60, 1.53)));

        residue->AddAtom(atomN);
        residue->AddAtom(atomCA);
        residue->AddAtom(atomCB);
        assembly->BuildStructureByDistance();
        outputPdbFile = assembly->BuildPdbFileStructureFromAssembly();
        outputPdbFile->Write(working_Directory + "OLT_superimposition_atoms.pdb");
        outputPdbFile->Write(working_Directory + "OLS_superimposition_atoms.pdb");

        residue->AddAtom(atomOG1);
        assembly->BuildStructureByDistance();
        outputPdbFile = assembly->BuildPdbFileStructureFromAssembly();
        outputPdbFile->Write(working_Directory + "OLS_AlignedToGlycan.pdb");

        residue->SetName("OLT"); //Thr = Ser + CG2
        residue->SetId("OLT_ _1_ _ _1");
        residue->AddAtom(atomCG2);
        atomOG1->SetName("OG1"); // It's OG in Ser.
        assembly->BuildStructureByDistance();
        outputPdbFile = assembly->BuildPdbFileStructureFromAssembly();
        outputPdbFile->Write(working_Directory + "OLT_AlignedToGlycan.pdb");
    }

    else if (amino_acid_name.compare("TYR")==0)
    {
        residue->SetName("OLY");
        residue->SetId("OLY_ _1_ _ _1");

        Atom *atomOH = new Atom(residue, "OH", (get_cartesian_point_from_internal_coords(atomC5, atomO5, atomC1, 112, 68, 1.46)));
        Atom *atomCZ = new Atom(residue, "CZ", (get_cartesian_point_from_internal_coords(atomO5, atomC1, atomOH, 117, 60, 1.35)));
        Atom *atomCE1 = new Atom(residue, "CE1", (get_cartesian_point_from_internal_coords(atomC1, atomOH, atomCZ, 120, 180, 1.37)));
        Atom *atomCD1 = new Atom(residue, "CD1", (get_cartesian_point_from_internal_coords(atomOH, atomCZ, atomCE1, 120, 180, 1.37)));
        Atom *atomCE2 = new Atom(residue, "CE2", (get_cartesian_point_from_internal_coords(atomC1, atomOH, atomCZ, 120, 0, 1.37)));
        Atom *atomCD2 = new Atom(residue, "CD2", (get_cartesian_point_from_internal_coords(atomOH, atomCZ, atomCE2, 120, 180, 1.37)));
        Atom *atomCG = new Atom(residue, "CG", (get_cartesian_point_from_internal_coords(atomCZ, atomCE2, atomCD2, 120, 0, 1.37)));
        Atom *atomCB = new Atom(residue, "CB", (get_cartesian_point_from_internal_coords(atomCE2, atomCD2, atomCG, 122, 180, 1.51)));
        Atom *atomCA = new Atom(residue, "CA", (get_cartesian_point_from_internal_coords(atomCD2, atomCG, atomCB, 111, -107, 1.55)));
        Atom *atomN = new Atom(residue, "N", (get_cartesian_point_from_internal_coords(atomCG, atomCB, atomCA, 114, -170, 1.44)));

        residue->AddAtom(atomN);
        residue->AddAtom(atomCA);
        residue->AddAtom(atomCB);
        assembly->BuildStructureByDistance();
        outputPdbFile = assembly->BuildPdbFileStructureFromAssembly();
        outputPdbFile->Write(working_Directory + "OLY_superimposition_atoms.pdb");

        residue->AddAtom(atomOH);
        residue->AddAtom(atomCZ);
        residue->AddAtom(atomCE1);
        residue->AddAtom(atomCD1);
        residue->AddAtom(atomCE2);
        residue->AddAtom(atomCD2);
        residue->AddAtom(atomCG);
        assembly->BuildStructureByDistance();
        outputPdbFile = assembly->BuildPdbFileStructureFromAssembly();
        outputPdbFile->Write(working_Directory + "OLY_AlignedToGlycan.pdb");

    }
    return;
}
