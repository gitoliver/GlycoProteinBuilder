#include <iostream>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <string>

// Change this next line to your PATH:
//#include "/home/ubunter/software/gems/gmml/includes/MolecularModeling/assembly.hpp"
#include "/home/oliver/Programs/gems/gmml/includes/MolecularModeling/assembly.hpp"
#include "/home/oliver/Programs/gems/gmml/includes/MolecularModeling/overlaps.hpp"
#include "resolve_overlaps.h"

using namespace std;
using namespace MolecularModeling;
using namespace gmml;


void resolve_overlaps::monte_carlo(Assembly glycoprotein, GlycoSiteVector glycosites){
  //glycosites contains pointers to the residues in glycoprotein that have a glycan attached to them. GetResidue()
  // Each glycosite can have multiple "rotamers" aka "glycan shapes" that are attached. This is due to an old design plan.

 // Assembly glycoprotein("5rsa_cheat_solve.pdb", PDB);
 // glycoprotein.BuildStructureByDistance();
  ResidueVector residues = glycoprotein.GetAllResiduesOfAssembly();
  /////////////////// SEED THE RANDOMNESS N STUFF //////////////////////////////
  int seed = time(NULL);
  srand(seed);
  cout << "Using seed: " << seed << endl;

  ///////////// Get pointers to protein and glycan parts of assembly ///////////
  AtomVector protein = glycoprotein.GetAllAtomsOfAssemblyWithinProteinResidues();
  AtomVector glycans = glycoprotein.GetAllAtomsOfAssemblyNotWithinProteinResidues();

  /////////////////// Get a base score ////////////////////////////////
  //double base_score = glycoprotein.CalculateAtomicOverlaps( &glycoprotein );
  double base_score = gmml::CalculateAtomicOverlaps(protein, glycans); // I've made a function you can pass atomvectors to.
  cout << "INITIAL SCORE: " << base_score << endl;
  int cycle = 1;

 // OG: This next part is unnecessary now, as have the residues in glycosites class
  /////////////////// MAKE A VECTOR OF THE RESIDUES TO MOVE ////////////////////
  std::vector<Residue*> residues_with_glycans; // basically all the residues you CAN move
  for (ResidueVector::iterator iter = residues.begin() ; iter!=residues.end() ; ++iter){
    if ( (*iter)->GetName().compare("NLN")==0 ){
      residues_with_glycans.push_back(*iter);
    }
  }

  //Here's an example of accessing the data in Glycosites:
  for(GlycoSiteVector::iterator it1 = glycosites.begin(); it1 != glycosites.end(); ++it1)
  {
      GlycosylationSite *glycosite = *it1; // I always do this so I can autocomplete GlycosylationSite functions. There are better ways now with auto, but this is C++98 compliant.
      Assembly *glycan = glycosite->GetAttachedRotamers().at(0)->GetAttachedRotamer(); // In future there will be only one rotamer per site, this looks silly as it is now.
      //Assembly *glycan = glycosite->GetGlycan(); // This is what it will be.
      Residue *linkage_residue = glycosite->GetResidue();
      AtomVector atoms = linkage_residue->GetAtoms(); // Fingers crossed these will be pointers to atoms in glycoprotein.
      Atom *atom1, *atom2, *atom3, *atom4, *atom5;
      for(AtomVector::iterator atom_iter = atoms.begin(); atom_iter != atoms.end(); ++atom_iter)
      {
          if ( (*atom_iter)->GetName().compare("N")==0 )
          {
              atom1 = *atom_iter;
          } // Do this with brackets so I can find pairs easily. Also it's the law. Also four spaces for indentation. Law's the law.
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
              atom5 = *atom_iter; // So can do chi2 at same time.
          }
      }
      double random_dihedral = (rand() % 360) + 1 - 180;
      std::cout << "Trying dihedral: " << random_dihedral << std::endl;
      glycoprotein.SetDihedral(atom1, atom2, atom3, atom4, random_dihedral); // CHI1
      random_dihedral = (rand() % 360) + 1 - 180;
      glycoprotein.SetDihedral(atom2, atom3, atom4, atom5, random_dihedral); // CHI2

      std::cout << "SCORE: " << gmml::CalculateAtomicOverlaps(protein, glycan->GetAllAtomsOfAssembly()) << std::endl;
  }
  PdbFileSpace::PdbFile *outputPdbFile1 = glycoprotein.BuildPdbFileStructureFromAssembly(-1,0);
  outputPdbFile1->Write("outputs/Glycoprotein_OlyExample.pdb");
  // END OF OG EXAMPLE. I left your old code below.

  // cout << "OVERLAP: " << glycoprotein.CalculateAtomicOverlaps( &glycoprotein ) << endl << endl;

  /////////////////// CHOOSE WHICH GLYCORESIDUE TO MOVE ////////////////////////
  cout << "GlycoResidues found: " << residues_with_glycans.size() << "\n\n";

  int do_the_thing_this_many_times = 10000;
  for ( int repss = 0; repss!=do_the_thing_this_many_times; repss++)
  {
    std::vector<Residue*> move_these_guys;  // basically all the residues you WILL move
    random_shuffle(residues_with_glycans.begin(), residues_with_glycans.end());
    int random_number_of_residues = rand() % residues_with_glycans.size();
    move_these_guys.push_back( residues_with_glycans.at(0) ); // choose one for now (will make it random later).

    /////////////////// ROTATE THE DIHEDRALS /////////////////////////////////////
    int how_many_rotamers = 1;
    for ( int reps = 0; reps!=how_many_rotamers; reps++){
      for (ResidueVector::iterator iter = move_these_guys.begin() ; iter!=move_these_guys.end() ; ++iter){

        cout << "CURRENT CYCLE: " << cycle << endl;
        AtomVector atoms = (*iter)->GetAtoms();
        Atom *atom1, *atom2, *atom3, *atom4;

        cout << "Current Target: ASN-chi1" << endl;
        for(AtomVector::iterator atom_iter = atoms.begin(); atom_iter != atoms.end(); ++atom_iter){
          if ( (*atom_iter)->GetName().compare("N")==0 ){
            atom1 = *atom_iter;}
          if ( (*atom_iter)->GetName().compare("CA")==0 ){
            atom2 = *atom_iter;}
          if ( (*atom_iter)->GetName().compare("CB")==0 ){
            atom3 = *atom_iter;}
          if ( (*atom_iter)->GetName().compare("CG")==0 ){
            atom4 = *atom_iter;}
        }

        // cout << ConvertRadian2Degree(glycoprotein.CalculateTorsionAngleByAtoms(atom1, atom2, atom3, atom4)) << endl;
        double random_dihedral = (rand() % 360) + 1 - 180;
        cout << "Trying dihedral: " << random_dihedral;
        // cout << "Selected dihedral: " << "\n" << atom1->GetId()
        //                               << "\n" << atom2->GetId()
        //                               << "\n" << atom3->GetId()
        //                               << "\n" << atom4->GetId() << "\n";
        glycoprotein.SetDihedral(atom1, atom2, atom3, atom4, random_dihedral);
        // cout << ConvertRadian2Degree(glycoprotein.CalculateTorsionAngleByAtoms(atom1, atom2, atom3, atom4)) << endl;

        cout << "\nCurrent Target: ASN-chi2" << endl;
        for(AtomVector::iterator atom_iter = atoms.begin(); atom_iter != atoms.end(); ++atom_iter){
          if ( (*atom_iter)->GetName().compare("CA")==0 ){
            atom1 = *atom_iter;}
          if ( (*atom_iter)->GetName().compare("CB")==0 ){
            atom2 = *atom_iter;}
          if ( (*atom_iter)->GetName().compare("CG")==0 ){
            atom3 = *atom_iter;}
          if ( (*atom_iter)->GetName().compare("OD1")==0 ){
            atom4 = *atom_iter;}
        }

        // cout << ConvertRadian2Degree(glycoprotein.CalculateTorsionAngleByAtoms(atom1, atom2, atom3, atom4)) << endl;
        random_dihedral = (rand() % 360) + 1 - 180;
        cout << "Trying dihedral: " << random_dihedral << "\n";
        // cout << "Selected dihedral: " << "\n" << atom1->GetId()
        //                               << "\n" << atom2->GetId()
        //                               << "\n" << atom3->GetId()
        //                               << "\n" << atom4->GetId() << "\n";
        glycoprotein.SetDihedral(atom1, atom2, atom3, atom4, random_dihedral);
        // cout << ConvertRadian2Degree(glycoprotein.CalculateTorsionAngleByAtoms(atom1, atom2, atom3, atom4)) << endl;
      }

      cout << endl;

      double current_score = gmml::CalculateAtomicOverlaps(protein, glycans);

      if ( (current_score - 20) < base_score ){
        cout << ("\noutput/test_test_"+to_string(cycle)+".pdb\n");
        PdbFileSpace::PdbFile *outputPdbFile = glycoprotein.BuildPdbFileStructureFromAssembly(-1,0);
        outputPdbFile->Write("output/test_test_"+to_string(cycle)+".pdb");
      }

      cout << "SCORE: " << current_score << endl;

      // despite the variable not being exactly right, i can MINIMIZE this
      cout << "---------------------\n"<< endl;
      cycle++;
    }
  }


  return;

}
