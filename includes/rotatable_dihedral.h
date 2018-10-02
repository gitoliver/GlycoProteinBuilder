#ifndef ROTATABLE_DIHEDRAL_H
#define ROTATABLE_DIHEDRAL_H

//#include "../../../includes/gmml.hpp"
#include "gmml.hpp"

using namespace MolecularModeling;

class Rotatable_dihedral
{
public:
    //////////////////////////////////////////////////////////
    //                    TYPE DEFINITION                   //
    //////////////////////////////////////////////////////////

    typedef std::vector<Rotatable_dihedral> RotatableDihedralVector;

    //////////////////////////////////////////////////////////
    //                       CONSTRUCTOR                    //
    //////////////////////////////////////////////////////////

    Rotatable_dihedral();
    Rotatable_dihedral(Atom *atom1, Atom *atom2, Atom *atom3, Atom *atom4);
    Rotatable_dihedral(AtomVector atoms);
    Rotatable_dihedral(AtomVector atoms, AtomVector atoms_that_move);


    //////////////////////////////////////////////////////////
    //                       ACCESSOR                       //
    //////////////////////////////////////////////////////////

    double CalculateDihedralAngle() const;
    AtomVector GetAtoms() const;
    AtomVector GetAtomsThatMove();
    double GetPreviousDihedralAngle();

    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////


    //////////////////////////////////////////////////////////
    //                       MUTATOR                        //
    //////////////////////////////////////////////////////////

    void DetermineAtomsThatMove();
    void SetDihedralAngle(double dihedral_angle);
    void ResetDihedralAngle();
    double RandomizeDihedralAngleWithinRanges(std::vector<std::pair<double,double>> ranges);
    double RandomizeDihedralAngle();
    double RandomizeDihedralAngleWithinRange(double min, double max);

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

    void SetAtoms(AtomVector atoms);
    void SetAtomsThatMove(AtomVector atoms);
    void RecordPreviousDihedralAngle(double dihedral_angle);

    //////////////////////////////////////////////////////////
    //                       ATTRIBUTES                     //
    //////////////////////////////////////////////////////////

    Atom *atom1_;
    Atom *atom2_;
    Atom *atom3_;
    Atom *atom4_;
    AtomVector atoms_that_move_;
    double previous_dihedral_angle_;

};

std::ostream& operator<<(std::ostream& os, const Rotatable_dihedral&);

#endif // ROTATABLE_DIHEDRAL_H
