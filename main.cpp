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

GeometryTopology::Coordinate subtract_coordinates(GeometryTopology::Coordinate minuaend, GeometryTopology::Coordinate subtrahend);

GeometryTopology::Coordinate get_cartesian_point_from_internal_coords(GeometryTopology::Coordinate a, GeometryTopology::Coordinate b, GeometryTopology::Coordinate c, double theta, double phi, double distance);

int main(int argc, char *argv[])
{
    QApplication abss(argc, argv);
    glycoproteinBuilder w;
    w.show();

    //************************************************//
    // Details for loading in a PDB file              //
    //************************************************//

    std::vector<std::string> amino_libs, glycam_libs, other_libs, prep;
    amino_libs.push_back("/home/oliver/Programs/gems-Nov2016/gmml/dat/CurrentParams/leaprc.ff12SB_2014-04-24/amino12.lib");
    amino_libs.push_back("/home/oliver/Programs/gems-Nov2016/gmml/dat/CurrentParams/leaprc.ff12SB_2014-04-24/aminoct12.lib");
    amino_libs.push_back("/home/oliver/Programs/gems-Nov2016/gmml/dat/CurrentParams/leaprc.ff12SB_2014-04-24/aminont12.lib");

    glycam_libs.push_back("/home/oliver/Programs/gems-Nov2016/gmml/dat/CurrentParams/leaprc_GLYCAM_06j-1_2014-03-14/GLYCAM_amino_06j_12SB.lib");
    glycam_libs.push_back("/home/oliver/Programs/gems-Nov2016/gmml/dat/CurrentParams/leaprc_GLYCAM_06j-1_2014-03-14/GLYCAM_aminoct_06j_12SB.lib");
    glycam_libs.push_back("/home/oliver/Programs/gems-Nov2016/gmml/dat/CurrentParams/leaprc_GLYCAM_06j-1_2014-03-14/GLYCAM_aminont_06j_12SB.lib");

    other_libs.push_back("/home/oliver/Programs/gems-Nov2016/gmml/dat/CurrentParams/leaprc.ff12SB_2014-04-24/nucleic12.lib");
    other_libs.push_back("/home/oliver/Programs/gems-Nov2016/gmml/dat/CurrentParams/leaprc.ff12SB_2014-04-24/nucleic12.lib");
    other_libs.push_back("/home/oliver/Programs/gems-Nov2016/gmml/dat/CurrentParams/other/solvents.lib");

    prep.push_back("/home/oliver/Programs/gems-Nov2016/gmml/dat/CurrentParams/leaprc_GLYCAM_06j-1_2014-03-14/GLYCAM_06j-1.prep");

    std::string pdb_file_path = "/home/oliver/Desktop/test.pdb";
    std::string parameter_file_path = "/home/oliver/Programs/gems-Nov2016/gmml/dat/CurrentParams/leaprc_GLYCAM_06j-1_2014-03-14/GLYCAM_06j.dat";
    std::string ion_parameter_file_path = "/home/oliver/Programs/gems-Nov2016/gmml/dat/CurrentParams/other/atomic_ions.lib";

    //************************************************//
    // Build PDB file from condensed sequence         //
    //************************************************//
/*
    Assembly assemblyC;
    assemblyC.BuildAssemblyFromCondensedSequence("DGlcpb1-4DNeup5Aca2-4DGlcpb1-OME", prep.at(0) ,parameter_file_path, true);
    assemblyC.BuildStructureByDistance();
    PdbFileSpace::PdbFile *outputPdbFile1 = assemblyC.BuildPdbFileStructureFromAssembly();
    outputPdbFile1->Write("/home/oliver/Desktop/test-BuildFromSequence.pdb");

    TopologyFileSpace::TopologyFile *outputTopologyFile1 = assemblyC.BuildTopologyFileStructureFromAssembly(parameter_file_path, ion_parameter_file_path);
    outputTopologyFile1->Write("/home/oliver/Desktop/test-BuildFromSequence.parm7");
*/

    //Assembly assemblyD;
    //assemblyD.BuildAssemblyFromPdbFile(pdb_file_path, amino_libs, glycam_libs, other_libs, prep, parameter_file_path );
    //assemblyD.BuildStructureByDistance();

    //************************************************//
    // Reading an input file                          //
    //************************************************//

    std::string proteinPDB, glycanDirectory, buffer;
    std::vector<std::string> glycositeList, listOfGlycans;

    std::ifstream inf ("/home/oliver/Desktop/input.txt");
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
            while(buffer != "END")
            {
                getline(inf, buffer);
                glycositeList.push_back(buffer);
                std::cout << "residue: " << buffer << std::endl;
            }
        }
        buffer = "Flushed";
        if(strInput == "Glycan id list:")
        {
            while(buffer != "END")
            {
                getline(inf, buffer);
                listOfGlycans.push_back(buffer);
                std::cout << "glycan: " << buffer << std::endl;
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

    Assembly assembly;
    /*
    std::ifstream fin;
    std::string directory = "/home/oliver/Desktop/" + glycanDirectory ;
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
            assembly.BuildAssemblyFromPdbFile(filepath, amino_libs, glycam_libs, other_libs, prep, parameter_file_path);
            assembly.BuildStructureByDistance();
        }
       // Endeavor to read a single number from the file and display it
       // fin.open( filepath.c_str() );
       // if (fin >> num)
       //     std::cout << filepath << ": " << num << std::endl;
       // fin.close();

    }
    closedir( dp );
    */
    //************************************************//
    // Building the necessary N-link Torsions         //
    //************************************************//

    pdb_file_path = "/home/oliver/Desktop/GAL-omega-gt.pdb";
    assembly.BuildAssemblyFromPdbFile(pdb_file_path, amino_libs, glycam_libs, other_libs, prep, parameter_file_path);
    assembly.BuildStructureByDistance();

    Residue* residue4YB = assembly.GetAllResiduesOfAssembly().at(1);
    std::cout << "residue4YB = " << residue4YB->GetName() << std::endl;

    AtomVector atoms4YB = residue4YB->GetAtoms();

    Atom* atomC5; // = residue4YB.FindAtom("C5");
    Atom* atomO5;
    Atom* atomC1;

    // Want: residue.FindAtom(string Name);
    for(AtomVector::iterator it = atoms4YB.begin(); it != atoms4YB.end(); it++)
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
    // Want: Constructor Atom atom(residue, name, coordinate).
    // Want: residue.DeleteAtom(int index);
    // Want: assembly.AddAtom(int index);
    std::cout << "Atoms are: " << atomC5->GetName() << ", " << atomO5->GetName() << ", " << atomC1->GetName() << std::endl;
    std::cout << std::endl;

    Assembly* assembly_NLN = new Assembly();
    Assembly* superimposition_atoms = new Assembly();
    Residue* residue_NLN = new Residue();
    residue_NLN->SetName("NLN");
    CoordinateVector Coords_Vector;

    Atom *atomND2 = new Atom;
  //  double distance =  1.42;
   // double angle    =  1.910629385; // (109.147 degrees)
    double distance = 1.53;
    double angle = 1.90765; // 109.3 degrees
    double torsion = 3.141592653; // (~180 degrees)
    GeometryTopology::Coordinate new_Coords = get_cartesian_point_from_internal_coords(atomC5->GetCoordinates().at(0), atomO5->GetCoordinates().at(0), atomC1->GetCoordinates().at(0), angle, torsion, distance);
    atomND2->AddCoordinate(new GeometryTopology::Coordinate(new_Coords.GetX(), new_Coords.GetY(), new_Coords.GetZ()) );
    atomND2->SetName("ND2");
    atomND2->SetDescription("Atom;");
    atomND2->SetId("ND2_1_NLN_ _1_ _ _1");
    atomND2->SetResidue(residue_NLN);

    //GeometryTopology::Coordinate *new_Coords_ptr = &new_Coords;
    //Coords_Vector.push_back(new_Coords_ptr);
    //Atom *atomND2 = new Atom(residue_NLN, "ND2", Coords_Vector);
    //Atom *atomND2 = new Atom(residue_NLN, "ND2",new_Coords);
  //  std::cout << "Id is " << atomND2->GetId() << std::endl;
    atomND2->SetId("ND2_1_NLN_ _1_ _ _1");

    Atom* atomCG = new Atom;
    distance = 1.325;
    angle = 1.90765; // 109.3 degrees
    torsion = 4.55531; // 261 degrees
    new_Coords = get_cartesian_point_from_internal_coords(atomO5->GetCoordinates().at(0), atomC1->GetCoordinates().at(0), atomND2->GetCoordinates().at(0), angle, torsion, distance);
    atomCG->AddCoordinate(new GeometryTopology::Coordinate(new_Coords.GetX(), new_Coords.GetY(), new_Coords.GetZ()) );
    atomCG->SetName("CG");
    atomCG->SetDescription("Atom;");
    atomCG->SetId("CG_2_NLN_ _1_ _ _1");
    atomCG->SetResidue(residue_NLN);

    Atom* atomOD1 = new Atom;
    distance = 1.22;
    angle = 2.19911; // 126 degrees
    torsion = 0; // 0 degrees
    new_Coords = get_cartesian_point_from_internal_coords(atomC1->GetCoordinates().at(0), atomND2->GetCoordinates().at(0), atomCG->GetCoordinates().at(0), angle, torsion, distance);
    atomOD1->AddCoordinate(new GeometryTopology::Coordinate(new_Coords.GetX(), new_Coords.GetY(), new_Coords.GetZ()) );
    atomOD1->SetName("OD1");
    atomOD1->SetDescription("Atom;");
    atomOD1->SetId("OD1_3_NLN_ _1_ _ _1");
    atomOD1->SetResidue(residue_NLN);

    Atom* atomCB = new Atom;
    distance = 1.53;
    angle = 1.98968; // 114 degrees
    torsion = 3.0945; // 177.3 degrees
    new_Coords = get_cartesian_point_from_internal_coords(atomC1->GetCoordinates().at(0), atomND2->GetCoordinates().at(0), atomCG->GetCoordinates().at(0), angle, torsion, distance);
    atomCB->AddCoordinate(new GeometryTopology::Coordinate(new_Coords.GetX(), new_Coords.GetY(), new_Coords.GetZ()) );
    atomCB->SetName("CB");
    atomCB->SetDescription("Atom;");
    atomCB->SetId("CB_4_NLN_ _1_ _ _1");
    atomCB->SetResidue(residue_NLN);

    Atom* atomCA = new Atom;
    distance = 1.53;
    angle = 1.93732; // 111 degrees
    torsion = 3.0997; // 177.6 degrees
    new_Coords = get_cartesian_point_from_internal_coords(atomND2->GetCoordinates().at(0), atomCG->GetCoordinates().at(0), atomCB->GetCoordinates().at(0), angle, torsion, distance);
    atomCA->AddCoordinate(new GeometryTopology::Coordinate(new_Coords.GetX(), new_Coords.GetY(), new_Coords.GetZ()) );
    atomCA->SetName("CA");
    atomCA->SetDescription("Atom;");
    atomCA->SetId("CA_5_NLN_ _1_ _ _1");
    atomCA->SetResidue(residue_NLN);

    Atom* atomN = new Atom;
    distance = 1.453;
    angle = 1.93732; // 111 degrees
    torsion = 3.3441; // 191.6 degrees
    new_Coords = get_cartesian_point_from_internal_coords(atomCG->GetCoordinates().at(0), atomCB->GetCoordinates().at(0), atomCA->GetCoordinates().at(0), angle, torsion, distance);
    atomN->AddCoordinate(new GeometryTopology::Coordinate(new_Coords.GetX(), new_Coords.GetY(), new_Coords.GetZ()) );
    atomN->SetName("N");
    atomN->SetDescription("Atom;");
    atomN->SetId("N_6_NLN_ _1_ _ _1");
    atomN->SetResidue(residue_NLN);


    residue_NLN->SetName("NLN");
    residue_NLN->AddAtom(atomN);
    residue_NLN->AddAtom(atomCA);
    residue_NLN->AddAtom(atomCB);

    residue_NLN->SetId("NLN_ _1_ _ _1");

    residue_NLN->SetAssembly(superimposition_atoms);
    superimposition_atoms->AddResidue(residue_NLN);
    superimposition_atoms->BuildStructureByDistance();
    PdbFileSpace::PdbFile *outputPdbFileSuperAtoms = superimposition_atoms->BuildPdbFileStructureFromAssembly(-1,0);
    outputPdbFileSuperAtoms->Write("/home/oliver/Desktop/superimposition_atoms.pdb");

    residue_NLN->AddAtom(atomCG);
    residue_NLN->AddAtom(atomOD1);
    residue_NLN->AddAtom(atomND2);

    residue_NLN->SetAssembly(assembly_NLN);
    assembly_NLN->AddResidue(residue_NLN);
    assembly_NLN->BuildStructureByDistance();
    //assembly_NLN->Print();
    PdbFileSpace::PdbFile *outputPdbFileNLN = assembly_NLN->BuildPdbFileStructureFromAssembly(-1,0);
    outputPdbFileNLN->Write("/home/oliver/Desktop/NLN_AlignedToGlycan.pdb");

/*
    torsion = 3.0944688; // 177.3
    new_Coords = get_cartesian_point_from_internal_coords(atomC1->GetCoordinates().at(0), atomND2.GetCoordinates().at(0), atomCG.GetCoordinates().at(0), angle, torsion, distance);
    Atom atomCB(&residue_NLN, "CB", Coords_Vector);

    torsion = 3.0997048; // 177.6
    new_Coords = get_cartesian_point_from_internal_coords(atomND2.GetCoordinates().at(0), atomCG.GetCoordinates().at(0), atomCB.GetCoordinates().at(0), angle, torsion, distance);
    Atom atomCA(&residue_NLN, "CA", Coords_Vector);

    torsion = 3.3440508; // 191.6
    new_Coords = get_cartesian_point_from_internal_coords(atomCG.GetCoordinates().at(0), atomCB.GetCoordinates().at(0), atomCA.GetCoordinates().at(0), angle, torsion, distance);
    Atom atomN(&residue_NLN, "N", Coords_Vector);

    torsion = 3.141592653; // (~180 degrees)
    new_Coords = get_cartesian_point_from_internal_coords(atomND2.GetCoordinates().at(0), atomCB.GetCoordinates().at(0), atomCG.GetCoordinates().at(0), angle, torsion, distance);
    Atom atomOD1(&residue_NLN, "OD1", Coords_Vector);
*/
   // assembly.AddAssembly(assembly_NLN);
 //   assembly_NLN->BuildStructureByDistance();
//    assembly.AddAssembly(assembly_NLN);
   // PdbFileSpace::PdbFile *outputPdbFileNLN = assembly_NLN->BuildPdbFileStructureFromAssembly();
    // outputPdbFileNLN->Write("/home/oliver/Desktop/test-NLN.pdb");

  /*
    GeometryTopology::Coordinate aO1(-0.847,   0.445,  -2.872);
    GeometryTopology::Coordinate aC1(0.535,   0.914,  -3.092);
    GeometryTopology::Coordinate aO5(1.398,   0.495,  -1.964);

    GeometryTopology::Coordinate expected_aN1(-1.251, -0.827, -3.599);

    double distance = 1.52, torsion = -1.415, angle = 2.03;

  //  distance = 1.52, torsion = -1.727, angle = 2.03;

    GeometryTopology::Coordinate actual_aN1 =  get_cartesian_point_from_internal_coords( aO5, aC1, aO1, angle, torsion, distance);

    expected_aN1.Print();
    std::cout << std::endl;
    actual_aN1.Print();
    std::cout << std::endl;
*/

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
    outputPdbFile->Write("/home/oliver/Desktop/test-Translate.pdb");
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

    outputPdbFile->Write("/home/oliver/Desktop/test-Dihedral.pdb");
    */
    //************************************************//
    // Atomic Overlaps                                //
    //************************************************//

    //std::cout << "Distance is " << GetDistanceToAtom(atom2, atom4) << std::endl;


    // Write a PDB file containing the moved co-ordinates
   // PdbFileSpace::PdbFile *outputPdbFileA = assemblyMoving.BuildPdbFileStructureFromAssembly();
  //  outputPdbFileA->Write("/home/oliver/Desktop/Superimposed-moved.pdb");


  //  std::cout << "Overlap is " << CalculateAtomicOverlaps(&assemblyA, &assemblyB) << std::endl;

    //assemblyA.AddAssembly(&assemblyB);
   // outputPdbFile->Write("/home/oliver/Desktop/test-MergeAssemblies.pdb");

    //************************************************//
    // Superimposition                                //
    //************************************************//

    Assembly protein;
    protein.BuildAssemblyFromPdbFile(proteinPDB, amino_libs, glycam_libs, other_libs, prep, parameter_file_path );
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
               // std::cout << "glycosite: " << glycosite << std::endl;
                std::cout << "glycosite id:" << id << std::endl;
                glycosites->push_back(residue);
                residue->SetName("NLN");
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

        pdb_file_path = "/home/oliver/Desktop/GAL-omega-gt.pdb";
        addedGlycans[i].BuildAssemblyFromPdbFile(pdb_file_path, amino_libs, glycam_libs, other_libs, prep, parameter_file_path);
        addedGlycans[i].BuildStructureByDistance();

        pdb_file_path = "/home/oliver/Desktop/NLN_AlignedToGlycan.pdb";
        Assembly *residueNLN = new Assembly();
        residueNLN->BuildAssemblyFromPdbFile(pdb_file_path, amino_libs, glycam_libs, other_libs, prep, parameter_file_path);
        residueNLN->BuildStructureByDistance();

        AssemblyVector *atomsToAlsoMove = new AssemblyVector; // using AddAssembly for the NLN part crashes the PDB writer.
        atomsToAlsoMove->push_back((&addedGlycans[i]));
        atomsToAlsoMove->push_back(residueNLN);

        pdb_file_path = "/home/oliver/Desktop/superimposition_atoms.pdb";
        Assembly *assemblyMoving = new Assembly();
        assemblyMoving->BuildAssemblyFromPdbFile(pdb_file_path, amino_libs, glycam_libs, other_libs, prep, parameter_file_path );

        //Superimpose(assemblyMoving, assemblyTarget, &addedGlycans[i]);
        Superimpose(assemblyMoving, assemblyTarget, atomsToAlsoMove);

        // Want: atom.SetX()
        // Want: atom.ReplaceCoords()
        // Replace Co-ords of sidechain
        // Want: assembly.Select(:residue@atom)
        AtomVector newSideChainAtoms = residueNLN->GetAllAtomsOfAssembly();
        for(AtomVector::iterator itt = atoms.begin(); itt != atoms.end(); ++itt)
        {
            Atom *atom = (*itt);
            std::cout << "Checking: " << atom->GetName() << std::endl;
            if ( (atom->GetName() == "ND2") || (atom->GetName() == "OD1") || (atom->GetName() == "CG") )
            {
                std::cout << "Replacing co-ords of protein " << atom->GetName() << std::endl;
                for(AtomVector::iterator ittt = newSideChainAtoms.begin(); ittt != newSideChainAtoms.end(); ++ittt)
                {
                    Atom *atom1 = (*ittt);
                    std::cout << "Comparing with " << atom1->GetName() << std::endl;
                    if (atom->GetName() == atom1->GetName() )
                    {

                            std::cout << "Before X=" << atom->GetCoordinates().at(0)->GetX() << std::endl;
                            atom->GetCoordinates().at(0)->SetX( atom1->GetCoordinates().at(0)->GetX() );
                            atom->GetCoordinates().at(0)->SetY( atom1->GetCoordinates().at(0)->GetY() );
                            atom->GetCoordinates().at(0)->SetZ( atom1->GetCoordinates().at(0)->GetZ() );
                            std::cout << "After X=" << atom->GetCoordinates().at(0)->GetX() << std::endl;
                    }
                }
            }
        }

        // Write out each glycan to check superimposition
        std::stringstream ss;
        ss << "/home/oliver/Desktop/GlycoProtein" << i << ".pdb";
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
        ss << "/home/oliver/Desktop/GlycoProtein" << i << ".pdb";
        i++;
        PdbFileSpace::PdbFile *outputPdbFileGlycoProtein = (*addedGlycan).BuildPdbFileStructureFromAssembly(-1,0);
        outputPdbFileGlycoProtein->Write(ss.str());
        protein.AddAssembly(&(*addedGlycan));
    }
    */
    //protein.Print();
    std::cout << "Got this far" << std::endl;
    PdbFileSpace::PdbFile *outputPdbFileGlycoProteinAll = protein.BuildPdbFileStructureFromAssembly(-1,0);
    outputPdbFileGlycoProteinAll->Write("/home/oliver/Desktop/GlycoProtein.pdb");
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
 //      outputPdbFileA->Write("/home/oliver/Desktop/test-Delete-Function.pdb");

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
    PdbFileSpace::PdbFile *outputPdbFile = moving->BuildPdbFileStructureFromAssembly();
    outputPdbFile->Write("/home/oliver/Desktop/test-Superimposed-Function.pdb");
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

GeometryTopology::Coordinate get_cartesian_point_from_internal_coords(GeometryTopology::Coordinate a, GeometryTopology::Coordinate b, GeometryTopology::Coordinate c, double theta, double phi, double distance)
{
// theta is the angle between 3 atoms. Phi is the torsion between 4 atoms.
Vector lmn_x, lmn_y, lmn_z;
double x_p, y_p, z_p;

Vector cb = subtract_coordinates(b, c);
Vector ba = subtract_coordinates(a, b);

cb.Print(); std::cout << "^cb" << std::endl;
ba.Print(); std::cout << "^ba" << std::endl;
std::cout << std::endl;

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

x_p = distance *  sin(theta) * cos(phi);
y_p = distance * sin(theta) * sin(phi);
z_p = distance * cos(theta);

//std::cout << "x_p=" << x_p << "y_p=" << y_p << "z_p=" << z_p << std::endl;

return (GeometryTopology::Coordinate)
    {
                lmn_x.GetX()*x_p + lmn_y.GetX()*y_p + lmn_z.GetX()*z_p + c.GetX(),
                lmn_x.GetY()*x_p + lmn_y.GetY()*y_p + lmn_z.GetY()*z_p + c.GetY(),
                lmn_x.GetZ()*x_p + lmn_y.GetZ()*y_p + lmn_z.GetZ()*z_p + c.GetZ()
    };

}
