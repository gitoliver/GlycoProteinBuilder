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

constexpr auto PI = 3.14159265358979323846;

using namespace MolecularModeling;
typedef std::vector<Atom*> AtomVector;
typedef std::vector<GeometryTopology::Coordinate*> CoordinateVector;
typedef std::vector<Residue*> ResidueVector;
typedef std::vector<Assembly*> AssemblyVector;
typedef GeometryTopology::Coordinate Vector;
typedef std::vector<GlycosylationSite*> GlycoSiteVector;
typedef std::vector<AttachedRotamer*> AttachedRotamerVector;

//void FindConnectedAtoms(Atom *atom, Assembly::AtomVector &visitedAtoms);

/*******************************************/
/* Function Declarations                   */
/*******************************************/

void GetCenterOfGeometry(Assembly *assembly, GeometryTopology::Coordinate *center);
double GetDistanceToAtom(Atom *A, Atom *otherAtom);
double CalculateAtomicOverlaps(Assembly *assemblyA, Assembly *assemblyB);

void GetResidueRingCenter (Residue *residue, GeometryTopology::Coordinate *center);
//GeometryTopology::Coordinate subtract_coordinates(GeometryTopology::Coordinate minuaend, GeometryTopology::Coordinate subtrahend);
//GeometryTopology::Coordinate get_cartesian_point_from_internal_coords(GeometryTopology::Coordinate a, GeometryTopology::Coordinate b, GeometryTopology::Coordinate c, double theta_Degrees, double phi_Degrees, double distance_Angstrom);
//GeometryTopology::Coordinate get_cartesian_point_from_internal_coords(Atom *a, Atom *b, Atom *c, double theta_Degrees, double phi_Degrees, double distance_Angstrom);
//void prepare_Glycans_For_Superimposition_To_Particular_Residue(std::string working_Directory, std::string amino_acid_name, Assembly *glycan);
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
    std::string proteinPDB, glycanDirectory, buffer, parameterDirectory;
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

    Assembly protein ((working_Directory + "/inputs/" + proteinPDB), gmml::InputFileType::PDB);
    protein.BuildStructureByDistance();

    ResidueVector protein_residues = protein.GetResidues();
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
                Assembly temp_assembly(filepath, gmml::InputFileType::PDB);
                temp_assembly.BuildStructureByDistance();
                glycosite->AddRotamer(new AttachedRotamer(temp_assembly));
            }
        }
    }
    closedir( dp );

    //************************************************//
    // Superimposition                                //
    //************************************************//
   // std::cout << "Superimposition" << std::endl;

    int j = 0, i = 0;
    Assembly glycoProtein;
    glycoProtein.AddAssembly(&protein); //Added assemblies get written into top of output PDB file.

    for (GlycoSiteVector::iterator it = glycoSites.begin(); it != glycoSites.end(); ++it)
    {
        GlycosylationSite* glycosite = *it;
        //std::cout << "Adding " << glycosite->GetGlycanName() << " to " << glycosite->GetResidue()->GetName() << std::endl;
        AttachedRotamerVector rotamers = glycosite->GetAttachedRotamers();
        for (AttachedRotamerVector::iterator itt = rotamers.begin(); itt != rotamers.end(); ++itt)
        {
            AttachedRotamer *rotamer = *itt;
           // std::cout << "Calling prepare glycans with " << glycosite->GetResidue()->GetName() << std::endl;
            rotamer->Prepare_Glycans_For_Superimposition_To_Particular_Residue(glycosite->GetResidue()->GetName());
            rotamer->Superimpose_Glycan_To_Glycosite(glycosite->GetResidue());
            //Write out a pdb file:
           // std::stringstream ss;
            //ss << working_Directory + "/outputs/addedGlycan_" << i << "_" << j << ".pdb";
            //PdbFileSpace::PdbFile *outputPdbFileGlycoProtein = rotamer->GetAttachedRotamer()->BuildPdbFileStructureFromAssembly(-1,0);
            //outputPdbFileGlycoProtein->Write(ss.str());

            // This is just temporary, I need to decide which rotamer to add
            if (j == 0)
            {
                glycoProtein.AddAssembly(rotamer->GetAttachedRotamer());
            }
            ++j;
        }
        j = 0; // reset rotamer counter. Used only for output file names.
        ++i; // increment residue counter. Used only for output file names.
    }

    PdbFileSpace::PdbFile *outputPdbFileGlycoProteinAll = glycoProtein.BuildPdbFileStructureFromAssembly(-1,0);
    outputPdbFileGlycoProteinAll->Write(working_Directory + "/outputs/GlycoProtein.pdb");

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
    outputPdbFile->Write(working_Directory + "/test-Translate.pdb");
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

    outputPdbFile->Write(working_Directory + "/test-Dihedral.pdb");
    */
    //************************************************//
    // Atomic Overlaps                                //
    //************************************************//

    //std::cout << "Distance is " << GetDistanceToAtom(atom2, atom4) << std::endl;


    // Write a PDB file containing the moved co-ordinates
   // PdbFileSpace::PdbFile *outputPdbFileA = assemblyMoving.BuildPdbFileStructureFromAssembly();
  //  outputPdbFileA->Write(working_Directory + "/Superimposed-moved.pdb");


   //std::cout << "Overlap is " << CalculateAtomicOverlaps(&assemblyA, &assemblyB) << std::endl;

    //assemblyA.AddAssembly(&assemblyB);
   // outputPdbFile->Write(working_Directory + "/test-MergeAssemblies.pdb");


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
 //      outputPdbFileA->Write(working_Directory + "/test-Delete-Function.pdb");



    //************************************************//
    // Build file from condensed sequence             //
    //************************************************//
/*
    if (fileExists(parameterDirectory + "/amino12.lib"))
        std::cout << "Using user provided parameters" << std::endl;
    else
        std::cout << "Using default parameters from installation folder" << std::endl;
        parameterDirectory = installation_Directory + "/CurrentParams";

    std::vector<std::string> amino_libs, glycam_libs, other_libs, prep;
    amino_libs.push_back(parameterDirectory + "/amino12.lib");
    amino_libs.push_back(parameterDirectory + "/aminoct12.lib");
    amino_libs.push_back(parameterDirectory + "/aminont12.lib");

    glycam_libs.push_back(parameterDirectory + "/GLYCAM_amino_06j_12SB.lib");
    glycam_libs.push_back(parameterDirectory + "/GLYCAM_aminoct_06j_12SB.lib");
    glycam_libs.push_back(parameterDirectory + "/GLYCAM_aminont_06j_12SB.lib");

    other_libs.push_back(parameterDirectory + "/nucleic12.lib");
    other_libs.push_back(parameterDirectory + "/nucleic12.lib");
    other_libs.push_back(parameterDirectory + "/solvents.lib");

    prep.push_back(parameterDirectory + "/GLYCAM_06j-1.prep");

    std::string pdb_file_path = working_Directory + "/outputs/" + "test.pdb";
    std::string parameter_file_path = parameterDirectory + "/GLYCAM_06j.dat";
    std::string ion_parameter_file_path = parameterDirectory + "/atomic_ions.lib";

    Assembly assemblyC;
    assemblyC.BuildAssemblyFromCondensedSequence("DGlcpb1-4[DManpa1-2]DNeup5Aca2-4DGlcpb1-OME", prep.at(0) ,parameter_file_path, true);
    assemblyC.BuildStructureByDistance();
    PdbFileSpace::PdbFile *outputPdbFile1 = assemblyC.BuildPdbFileStructureFromAssembly();
    outputPdbFile1->Write(working_Directory + "/test-buildFromStructure.pdb");

    TopologyFileSpace::TopologyFile *outputTopologyFile1 = assemblyC.BuildTopologyFileStructureFromAssembly(parameter_file_path, ion_parameter_file_path);
    outputTopologyFile1->Write(working_Directory + "/test-BuildFromSequence.parm7");

    CoordinateFileSpace::CoordinateFile *outputCoordinateFile = assemblyC.BuildCoordinateFileStructureFromAssembly();
    outputCoordinateFile->Write(working_Directory + "/test-BuildFromSequence.rst7");

*/
    //Assembly assemblyD;
    //assemblyD.BuildAssemblyFromPdbFile(pdb_file_path, amino_libs, glycam_libs, other_libs, prep, parameter_file_path );
    //assemblyD.BuildStructureByDistance();


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

/* Have moved this to GMML. Can delete if all is ok.
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

void TestFind3DAffineTransform()
{

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

/*
void prepare_Glycans_For_Superimposition_To_Particular_Residue(std::string working_Directory, std::string amino_acid_name, Assembly *glycan)
{
    //Dear future self, the order that you add the atoms to the residue matters for superimposition ie N, CA, CB , not CB, CA, N.

    Residue* reducing_Residue = glycan->GetAllResiduesOfAssembly().at(1);
    std::cout << "Reducing residue is " << reducing_Residue->GetName() << std::endl;
    AtomVector reducing_Atoms = reducing_Residue->GetAtoms();

    Atom* atomC5;
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
        residue->SetId("NLN_?_1_?_?_1");

        Atom *atomND2 = new Atom(residue, "ND2", (gmml::get_cartesian_point_from_internal_coords(atomC5, atomO5, atomC1, 109.3, 180, 1.53)));
        Atom *atomCG = new Atom(residue, "CG", (gmml::get_cartesian_point_from_internal_coords(atomO5, atomC1, atomND2, 109.3, 261, 1.325)));
        Atom *atomOD1 = new Atom(residue, "OD1", (gmml::get_cartesian_point_from_internal_coords(atomC1, atomND2, atomCG, 126, 0, 1.22)));
        Atom *atomCB = new Atom(residue, "CB", (gmml::get_cartesian_point_from_internal_coords(atomC1, atomND2, atomCG, 114, 177.3, 1.53)));
        Atom *atomCA = new Atom(residue, "CA", (gmml::get_cartesian_point_from_internal_coords(atomND2, atomCG, atomCB, 111, 177.6, 1.53)));
        Atom *atomN = new Atom(residue, "N", (gmml::get_cartesian_point_from_internal_coords(atomCG, atomCB, atomCA, 111, 191.6, 1.453)));

        residue->AddAtom(atomN);
        residue->AddAtom(atomCA);
        residue->AddAtom(atomCB);
        assembly->BuildStructureByDistance();
        outputPdbFile = assembly->BuildPdbFileStructureFromAssembly();
        outputPdbFile->Write(working_Directory + "/outputs/NLN_superimposition_atoms.pdb");

        residue->AddAtom(atomCG);
        residue->AddAtom(atomOD1);
        residue->AddAtom(atomND2);
        assembly->BuildStructureByDistance();
        outputPdbFile = assembly->BuildPdbFileStructureFromAssembly();
        outputPdbFile->Write(working_Directory + "/outputs/NLN_AlignedToGlycan.pdb");
    }
    else if (amino_acid_name.compare("THR")==0 || amino_acid_name.compare("SER")==0)
    {
        residue->SetName("OLS");
        residue->SetId("OLS_?_1_?_?_1");

        Atom *atomOG1 = new Atom(residue, "OG", (gmml::get_cartesian_point_from_internal_coords(atomC5, atomO5, atomC1, 112, 68, 1.46)));
        Atom *atomCB = new Atom(residue, "CB", (gmml::get_cartesian_point_from_internal_coords(atomO5, atomC1, atomOG1, 109.3, 75, 1.53)));
        Atom *atomCA = new Atom(residue, "CA", (gmml::get_cartesian_point_from_internal_coords(atomC1, atomOG1, atomCB, 109.3, 125, 1.53)));
        Atom *atomN = new Atom(residue, "N", (gmml::get_cartesian_point_from_internal_coords(atomOG1, atomCB, atomCA, 109.3, 180, 1.53)));
        Atom *atomCG2 = new Atom(residue, "CG2", (gmml::get_cartesian_point_from_internal_coords(atomC1, atomOG1, atomCB, 109.3, -60, 1.53)));

        residue->AddAtom(atomN);
        residue->AddAtom(atomCA);
        residue->AddAtom(atomCB);
        assembly->BuildStructureByDistance();
        outputPdbFile = assembly->BuildPdbFileStructureFromAssembly();
        outputPdbFile->Write(working_Directory + "/outputs/OLT_superimposition_atoms.pdb");
        outputPdbFile->Write(working_Directory + "/outputs/OLS_superimposition_atoms.pdb");

        residue->AddAtom(atomOG1);
        assembly->BuildStructureByDistance();
        outputPdbFile = assembly->BuildPdbFileStructureFromAssembly();
        outputPdbFile->Write(working_Directory + "/outputs/OLS_AlignedToGlycan.pdb");

        residue->SetName("OLT"); //Thr = Ser + CG2
        residue->SetId("OLT_?_1_?_?_1");
        residue->AddAtom(atomCG2);
        atomOG1->SetName("OG1"); // It's OG in Ser.
        assembly->BuildStructureByDistance();
        outputPdbFile = assembly->BuildPdbFileStructureFromAssembly();
        outputPdbFile->Write(working_Directory + "/outputs/OLT_AlignedToGlycan.pdb");
    }

    else if (amino_acid_name.compare("TYR")==0)
    {
        residue->SetName("OLY");
        residue->SetId("OLY_?_1_?_?_1");

        Atom *atomOH = new Atom(residue, "OH", (gmml::get_cartesian_point_from_internal_coords(atomC5, atomO5, atomC1, 112, 68, 1.46)));
        Atom *atomCZ = new Atom(residue, "CZ", (gmml::get_cartesian_point_from_internal_coords(atomO5, atomC1, atomOH, 117, 60, 1.35)));
        Atom *atomCE1 = new Atom(residue, "CE1", (gmml::get_cartesian_point_from_internal_coords(atomC1, atomOH, atomCZ, 120, 180, 1.37)));
        Atom *atomCD1 = new Atom(residue, "CD1", (gmml::get_cartesian_point_from_internal_coords(atomOH, atomCZ, atomCE1, 120, 180, 1.37)));
        Atom *atomCE2 = new Atom(residue, "CE2", (gmml::get_cartesian_point_from_internal_coords(atomC1, atomOH, atomCZ, 120, 0, 1.37)));
        Atom *atomCD2 = new Atom(residue, "CD2", (gmml::get_cartesian_point_from_internal_coords(atomOH, atomCZ, atomCE2, 120, 180, 1.37)));
        Atom *atomCG = new Atom(residue, "CG", (gmml::get_cartesian_point_from_internal_coords(atomCZ, atomCE2, atomCD2, 120, 0, 1.37)));
        Atom *atomCB = new Atom(residue, "CB", (gmml::get_cartesian_point_from_internal_coords(atomCE2, atomCD2, atomCG, 122, 180, 1.51)));
        Atom *atomCA = new Atom(residue, "CA", (gmml::get_cartesian_point_from_internal_coords(atomCD2, atomCG, atomCB, 111, -107, 1.55)));
        Atom *atomN = new Atom(residue, "N", (gmml::get_cartesian_point_from_internal_coords(atomCG, atomCB, atomCA, 114, -170, 1.44)));

        residue->AddAtom(atomN);
        residue->AddAtom(atomCA);
        residue->AddAtom(atomCB);
        assembly->BuildStructureByDistance();
        outputPdbFile = assembly->BuildPdbFileStructureFromAssembly();
        outputPdbFile->Write(working_Directory + "/outputs/OLY_superimposition_atoms.pdb");

        residue->AddAtom(atomOH);
        residue->AddAtom(atomCZ);
        residue->AddAtom(atomCE1);
        residue->AddAtom(atomCD1);
        residue->AddAtom(atomCE2);
        residue->AddAtom(atomCD2);
        residue->AddAtom(atomCG);
        assembly->BuildStructureByDistance();
        outputPdbFile = assembly->BuildPdbFileStructureFromAssembly();
        outputPdbFile->Write(working_Directory + "/outputs/OLY_AlignedToGlycan.pdb");

    }
    return;
}
*/
