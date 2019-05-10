#ifndef RESIDUE_LINKAGE_H
#define RESIDUE_LINKAGE_H
/*
 * This class figures out the rotatable bonds between two residues
 * Starts/ends at the CA atoms in proteins. Looks for cycles (as they aren't rotatable).
 * Stores each rotatable bond as a rotatable_dihedral object.
 */
#include "gmml.hpp"
#include "rotatable_dihedral.h"
#include "selections.h"

// Seed with a real random value, if available
static pcg_extras::seed_seq_from<std::random_device> seed_source;
// Make a random number engine
static pcg32 rng(seed_source);

typedef std::vector<Rotatable_dihedral> RotatableDihedralVector;
class GlycosylationSite;
class Residue_linkage
{
public:
    //////////////////////////////////////////////////////////
    //                    TYPE DEFINITION                   //
    //////////////////////////////////////////////////////////

//    typedef std::vector<Rotatable_dihedral> RotatableDihedralVector;
//    typedef std::vector<Residue_linkage> ResidueLinkageVector;

    //////////////////////////////////////////////////////////
    //                       CONSTRUCTOR                    //
    //////////////////////////////////////////////////////////

    Residue_linkage();
    Residue_linkage(Residue *residue1, Residue *residue2, GlycosylationSite *glycosite);

    //////////////////////////////////////////////////////////
    //                       ACCESSOR                       //
    //////////////////////////////////////////////////////////

    ResidueVector GetResidues();
    RotatableDihedralVector GetRotatableDihedrals() const;
    RotatableDihedralVector GetRotatableDihedralsWithMultipleRotamers();
    int GetNumberOfRotatableDihedrals();
    int GetNumberOfShapes();
    Residue* GetFromThisResidue1();
    Residue* GetToThisResidue2();
    Atom* GetFromThisConnectionAtom1();
    Atom* GetToThisConnectionAtom2();
    bool CheckIfConformer();
    GlycosylationSite* GetAssociatedGlycosite();

    //////////////////////////////////////////////////////////
    //                       MUTATOR                        //
    //////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////

    void GenerateAllShapesUsingMetadata();
    void SetDefaultShapeUsingMetadata();
    void SetRandomShapeUsingMetadata(bool useRanges = true);
    void SetSpecificShapeUsingMetadata(int shapeNumber, bool useRanges = false);
    void SetCustomDihedralAngles(std::vector <double> dihedral_angles);
    void SetShapeToPrevious();
    void SetRandomDihedralAngles();
    void DetermineAtomsThatMove();
    void SetAssociatedGlycosylationSite(GlycosylationSite *glycosite);

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

    void InitializeClass(Residue *from_this_residue1, Residue *to_this_residue2, GlycosylationSite *glycosite);
    RotatableDihedralVector FindRotatableDihedralsConnectingResidues(Atom *from_this_connection_atom1, Atom *to_this_connection_atom2);
    // Previous function generates a list of linearly connected atoms that define the rotatable bonds
    // This function splits that list into groups of 4 and creates rotatable_dihedral objects
    RotatableDihedralVector SplitAtomVectorIntoRotatableDihedrals(AtomVector atoms);
    gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector FindMetadata(Atom *from_this_connection_atom1, Atom *to_this_connection_atom2);
    void AddMetadataToRotatableDihedrals(gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector metadata);
    void SetResidues(Residue *residue1, Residue *residue2);
    void SetConnectionAtoms(Residue *residue1, Residue *residue2);
    void SetConformerUsingMetadata(bool useRanges = false, int conformerNumber = 0);

    //////////////////////////////////////////////////////////
    //                       ATTRIBUTES                     //
    //////////////////////////////////////////////////////////

    GlycosylationSite* associatedGlycosite_;
    Residue* from_this_residue1_;
    Residue* to_this_residue2_;
    Atom* from_this_connection_atom1_;
    Atom* to_this_connection_atom2_;
    RotatableDihedralVector rotatable_dihedrals_;
    //gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector metadata_;
};

std::ostream& operator<<(std::ostream& os, const Residue_linkage&);

#endif // RESIDUE_LINKAGE_H
