#ifndef RESIDUE_LINKAGE_H
#define RESIDUE_LINKAGE_H
/*
 * This class figures out the rotatable bonds between two residues
 * Starts/ends at the CA atoms in proteins. Looks for cycles (as they aren't rotatable).
 * Stores each rotatable bond as a rotatable_dihedral object.
 */
#include "rotatable_dihedral.h"
#include "../../../includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"

class Residue_linkage
{
public:
    //////////////////////////////////////////////////////////
    //                    TYPE DEFINITION                   //
    //////////////////////////////////////////////////////////

    typedef std::vector<Rotatable_dihedral> RotatableDihedralVector;
    typedef std::vector<Residue_linkage> ResidueLinkageVector;

    //////////////////////////////////////////////////////////
    //                       CONSTRUCTOR                    //
    //////////////////////////////////////////////////////////

    Residue_linkage();
    Residue_linkage(Residue *residue1, Residue *residue2);

    //////////////////////////////////////////////////////////
    //                       ACCESSOR                       //
    //////////////////////////////////////////////////////////

    ResidueVector GetResidues();
    RotatableDihedralVector GetRotatableDihedrals() const;
    int GetNumberOfRotatableBonds();

    //////////////////////////////////////////////////////////
    //                       MUTATOR                        //
    //////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////

    void SetDefaultDihedralAnglesUsingMetadata();
    void SetRandomDihedralAnglesUsingMetadata();
    void SetCustomDihedralAngles(std::vector <double> dihedral_angles);
    void SetPreviousDihedralAngles();
    void SetRandomDihedralAngles();
    void DetermineAtomsThatMove();

    //////////////////////////////////////////////////////////
    //                       DISPLAY FUNCTION               //
    //////////////////////////////////////////////////////////

    void Print();

    //////////////////////////////////////////////////////////
    //                       OPERATORS                      //
    //////////////////////////////////////////////////////////

private:

    //////////////////////////////////////////////////////////
    //                    PRIVATE FUNCTIONS                 //
    //////////////////////////////////////////////////////////

    void InitializeClass(Residue *residue1, Residue *residue2);
    RotatableDihedralVector FindRotatableBondsConnectingResidues(Atom *connection_atom1, Atom *connection_atom2);
    // Previous function generates a list of linearly connected atoms that define the rotatable bonds
    // This function splits that list into groups of 4 and creates rotatable_dihedral objects
    RotatableDihedralVector SplitAtomVectorIntoRotatableBonds(AtomVector atoms);
    gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector FindMetadata(Atom *connection_atom1, Atom *connection_atom2);
    void AddMetadataToRotatableDihedrals(gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector metadata);
    void SetResidues(Residue *residue1, Residue *residue2);
    void SetConnectionAtoms(Residue *residue1, Residue *residue2);

    //////////////////////////////////////////////////////////
    //                       ATTRIBUTES                     //
    //////////////////////////////////////////////////////////

    Residue* residue1_;
    Residue* residue2_;
    Atom* connection_atom1_;
    Atom* connection_atom2_;
    RotatableDihedralVector rotatable_bonds_;
    //gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector metadata_;
};

std::ostream& operator<<(std::ostream& os, const Residue_linkage&);

#endif // RESIDUE_LINKAGE_H
