#ifndef RESIDUE_LINKAGE_H
#define RESIDUE_LINKAGE_H
/*
 * This class figures out the rotatable bonds between two residues
 * Starts/ends at the CA atoms in proteins. Looks for cycles (as they aren't rotatable).
 * Stores each rotatable bond as a rotatable_dihedral object.
 */
#include "rotatable_dihedral.h"
#include "selections.h"

class Residue_linkage
{
public:
    //////////////////////////////////////////////////////////
    //                    TYPE DEFINITION                   //
    //////////////////////////////////////////////////////////

    typedef std::vector<Rotatable_dihedral> RotatableDihedralVector;

    //////////////////////////////////////////////////////////
    //                       CONSTRUCTOR                    //
    //////////////////////////////////////////////////////////

    Residue_linkage();
    Residue_linkage(Residue *residue1, Residue *residue2);
    Residue_linkage(Residue *residue1, Residue *residue2, std::vector <double> dihedral_angles);



    //////////////////////////////////////////////////////////
    //                       ACCESSOR                       //
    //////////////////////////////////////////////////////////

    ResidueVector GetResidues();
    RotatableDihedralVector GetRotatableDihedrals() const;

    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////

    RotatableDihedralVector FindRotatableBondsConnectingResidues(Residue *first_residue, Residue *second_residue);
    // Previous function generates a list of linearly connected atoms that define the rotatable bonds
    // This function splits that list into groups of 4 and creates rotatable_dihedral objects
    RotatableDihedralVector SplitAtomVectorIntoRotatableBonds(AtomVector atoms);
    // This is bad, but sets Chi1 and Chi2 to 180. These bonds are present in protein residues except Gly.
    void SetReasonableChi1Chi2DihedralAngles();

    //////////////////////////////////////////////////////////
    //                       MUTATOR                        //
    //////////////////////////////////////////////////////////

    void DetermineAtomsThatMove();
    void SetResidues(Residue *residue1, Residue *residue2);
    void SetDihedralAngles(std::vector <double> dihedral_angles);
    void ResetDihedralAngles();
    double RandomizeDihedralAnglesWithinTheirRanges();
    double RandomizeDihedralAngles();

    //////////////////////////////////////////////////////////
    //                       DISPLAY FUNCTION               //
    //////////////////////////////////////////////////////////

    void Print();

    //////////////////////////////////////////////////////////
    //                       OPERATORS                      //
    //////////////////////////////////////////////////////////

private:

    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////


    //////////////////////////////////////////////////////////
    //                       ATTRIBUTES                     //
    //////////////////////////////////////////////////////////

    Residue* residue1_;
    Residue* residue2_;
    RotatableDihedralVector rotatable_bonds_;

};

std::ostream& operator<<(std::ostream& os, const Residue_linkage&);

#endif // RESIDUE_LINKAGE_H
