#ifndef GLYCOSYLATIONSITE_H
#define GLYCOSYLATIONSITE_H


#include <iomanip> // For setting precision and formating in std::cout 
#include <algorithm> //  std::erase, std::remove

#include "gmml.hpp"
#include "residue_linkage.h"
#include "rotatable_dihedral.h"

typedef std::vector<Residue_linkage> ResidueLinkageVector;
typedef std::vector<Rotatable_dihedral> RotatableDihedralVector;


using namespace MolecularModeling;

class GlycosylationSite
{
public:
    //////////////////////////////////////////////////////////
    //                    TYPE DEFINITION                   //
    //////////////////////////////////////////////////////////

    typedef std::vector<GlycosylationSite> GlycosylationSiteVector;
    typedef std::vector<GlycosylationSite*> GlycosylationSitePointerVector;


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
    ResidueLinkageVector GetRotatableBonds();


    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////

    void AttachGlycan(Assembly glycan, Assembly &glycoprotein);
    double Calculate_bead_overlaps(std::string overlap_type = "total", bool record = true);
    double Calculate_and_print_bead_overlaps();
    double CalculateAtomicOverlaps();
    void UpdateAtomsThatMoveInLinkages();
    void Rename_Protein_Residue_From_GLYCAM_To_Standard();
    void Wiggle(int *output_pdb_id, double tolerance = 0.1, int interval = 5);
    void WiggleFirstLinkage(int *output_pdb_id, double tolerance = 0.1, int interval = 5);

    // Do not keep this here:
    void write_pdb_file(MolecularModeling::Assembly *glycoprotein, int cycle, std::string summary_filename, double overlap);


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
    void SetDefaultDihedralAnglesUsingMetadata();
    void SetRandomDihedralAnglesUsingMetadata();
    void SetRandomDihedralAnglesUsingMetadataForNthLinkage(int linkage_number);
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
    //                  PRIVATE FUNCTIONS                   //
    //////////////////////////////////////////////////////////

    void Prepare_Glycans_For_Superimposition_To_Particular_Residue(std::string amino_acid_name);
    void Superimpose_Glycan_To_Glycosite(Residue *glycosite_residue);
    double CalculateTorsionAngle(AtomVector atoms);
    double Calculate_bead_overlaps(AtomVector &atomsA, AtomVector &atomsB);
    //void SetRotatableBonds(Residue *residue1, Residue *residue2);
    void Rename_Protein_Residue_To_GLYCAM_Nomenclature();
    void FigureOutResidueLinkagesInGlycan(Residue *from_this_residue1, Residue *to_this_residue2, ResidueLinkageVector *residue_linkages);
    void RecursivelyGetAllNeighboringResidues(Atom* current_atom, ResidueVector* neighbors);
    Atom* GetConnectingProteinAtom(std::string residue_name);
    void WiggleOneLinkage(Residue_linkage &linkage, int *output_pdb_id, double tolerance = 0.1, int interval = 5);

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
    //Residue_linkage residue_linkage_; // This should become a vector of residue_linkages?
    ResidueLinkageVector all_residue_linkages_;
    AtomVector self_glycan_beads_;
    AtomVector other_glycan_beads_;
    AtomVector protein_beads_;
};

#endif // GLYCOSYLATIONSITE_H
