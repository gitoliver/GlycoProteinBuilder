#ifndef GLYCOSYLATIONSITE_H
#define GLYCOSYLATIONSITE_H

#include "gmml.hpp"
#include "residue_linkage.h"
#include "selections.h"
#include <iomanip> // For setting precision and formating in std::cout 
#include <algorithm> //  std::erase, std::remove

using namespace MolecularModeling;

class GlycosylationSite
{
public:
    //////////////////////////////////////////////////////////
    //                    TYPE DEFINITION                   //
    //////////////////////////////////////////////////////////

    typedef std::vector<GlycosylationSite> GlycosylationSiteVector;
    typedef std::vector<GlycosylationSite*> GlycosylationSitePointerVector;
   // typedef std::vector<Rotatable_dihedral> RotatableDihedralVector;


    //////////////////////////////////////////////////////////
    //                       CONSTRUCTOR                    //
    //////////////////////////////////////////////////////////

    GlycosylationSite();
    GlycosylationSite(std::string glycan_name);
    GlycosylationSite(std::string glycan_name, std::string residue_number);
    //GlycosylationSite(std::string glycan_name, std::string residue_number_, Assembly glycan, Residue* residue);
    ~GlycosylationSite();
    //////////////////////////////////////////////////////////
    //                       ACCESSOR                       //
    //////////////////////////////////////////////////////////

    std::string GetGlycanName();
    std::string GetResidueNumber();
    Residue* GetResidue();
    Assembly* GetAttachedGlycan();
    Assembly* GetGlycoprotein();
    double GetOverlap();
    double GetWeightedOverlap(double glycan_weight, double protein_weight);
    double GetGlycanOverlap();
    double GetProteinOverlap();
    AtomVector GetSelfGlycanBeads();
    AtomVector GetProteinBeads();
    AtomVector GetOtherGlycanBeads();
    Residue_linkage GetRotatableBonds();



    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////
    void AttachGlycan(Assembly glycan, Assembly &glycoprotein);
    double Calculate_bead_overlaps(std::string overlap_type = "total", bool record = true);
    double Calculate_and_print_bead_overlaps();
    void UpdateAtomsThatMoveInLinkages();
    void Rename_Protein_Residue_From_GLYCAM_To_Standard();

   // void SetChiAtoms(Residue* residue);


    //////////////////////////////////////////////////////////
    //                       MUTATOR                        //
    //////////////////////////////////////////////////////////

    void SetGlycanName(std::string glycan_name);
    void SetResidueNumber(std::string residue_number);
    void SetResidue(Residue* residue);
    void SetGlycan(Assembly glycan);
    void SetGlycanOverlap(double overlap);
    void SetProteinOverlap(double overlap);
    void SetSelfGlycanBeads(AtomVector *beads);
    void SetProteinBeads(AtomVector *beads);
    void SetOtherGlycanBeads(AtomVector *beads);
    void RandomizeDihedralAngles();
    void SetResonableChi1Chi2DihedralAngles();
    void ResetDihedralAngles();

    //////////////////////////////////////////////////////////
    //                       DISPLAY FUNCTION               //
    //////////////////////////////////////////////////////////

    void Print_bead_overlaps();
    void Print(std::string type = "All");

    //////////////////////////////////////////////////////////
    //                       OPERATORS                      //
    //////////////////////////////////////////////////////////

    inline bool operator==(const GlycosylationSite &rhs) const
    {
        return rhs.residue_number_ == residue_number_;
    }

private:

    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////

    void Prepare_Glycans_For_Superimposition_To_Particular_Residue(std::string amino_acid_name);
    void Superimpose_Glycan_To_Glycosite(Residue *glycosite_residue);
    double CalculateTorsionAngle(AtomVector atoms);
    double Calculate_bead_overlaps(AtomVector &atomsA, AtomVector &atomsB);
    void SetRotatableBonds(Residue *residue1, Residue *residue2);
    void Rename_Protein_Residue_To_GLYCAM_Nomenclature();

    //////////////////////////////////////////////////////////
    //                       ATTRIBUTES                     //
    //////////////////////////////////////////////////////////

    std::string glycan_name_;
    std::string residue_number_;
    Residue* residue_;                                  /*!< A pointer back to the residue for this glycosite >*/
    Assembly glycan_;
    AtomVector superimposition_atoms_;               /*!< The 3 atoms used for superimposition of glycan to sidechain >*/
    double glycan_overlap_;
    double protein_overlap_;
    Residue_linkage rotatable_bonds_; // This should become a vector of residue_linkages?
    AtomVector self_glycan_beads_;
    AtomVector other_glycan_beads_;
    AtomVector protein_beads_;
};

#endif // GLYCOSYLATIONSITE_H
