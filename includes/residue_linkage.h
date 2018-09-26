#ifndef RESIDUE_LINKAGE_H
#define RESIDUE_LINKAGE_H



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
    RotatableDihedralVector GetRotatableDihedrals();

    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////

    RotatableDihedralVector FindRotatableBondsConnectingResidues(Residue *first_residue, Residue *second_residue);
    RotatableDihedralVector SplitAtomVectorIntoRotatableBonds(AtomVector atoms);
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

#endif // RESIDUE_LINKAGE_H
