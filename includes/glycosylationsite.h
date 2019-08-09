#ifndef GLYCOSYLATIONSITE_H
#define GLYCOSYLATIONSITE_H


#include <iomanip> // For setting precision and formating in std::cout 
#include <algorithm> //  std::erase, std::remove

#include "gmml.hpp"
#include "residue_linkage.h"
#include "rotatable_dihedral.h"

typedef std::vector<Residue_linkage> ResidueLinkageVector;
typedef std::vector<Rotatable_dihedral> RotatableDihedralVector;

enum OverlapType
{
    BEAD,
    ATOMIC,
};

enum MoleculeType
{
    PROTEIN,
    GLYCAN,
    ALL
};

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

    Assembly* GetGlycoprotein();
    std::string GetGlycanName();
    std::string GetResidueNumber();
    Residue* GetResidue();
    Assembly* GetAttachedGlycan();
    double GetOverlap();
    double GetWeightedOverlap(double glycan_weight, double protein_weight);
    double GetGlycanOverlap();
    double GetProteinOverlap();
    AtomVector GetSelfGlycanBeads();
    AtomVector GetProteinBeads();
    AtomVector GetOtherGlycanBeads();
    ResidueLinkageVector GetRotatableBonds();
    ResidueLinkageVector GetFirstAnd1_6Linkages();
    ResidueLinkageVector GetFirstAnd2_XLinkages();



    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////

    void AttachGlycan(Assembly glycan, Assembly &glycoprotein);
    double CalculateOverlaps(OverlapType overlapType = BEAD, MoleculeType moleculeType = ALL, bool recordOverlap = true, bool printOverlap = false);
    void UpdateAtomsThatMoveInLinkages();
    void StashCoordinates();
    void SetStashedCoordinates();
    void Rename_Protein_Residue_From_GLYCAM_To_Standard();
    void Wiggle(int *output_pdb_id, double tolerance = 0.1, int interval = 5);
    void WiggleFirstLinkage(int *output_pdb_id, double tolerance = 0.1, int interval = 5);
    GlycosylationSiteVector GetXClosestSitesWithinOverlapDistanceY(GlycosylationSiteVector &glycosites, int maxNumberOfSitesToConsider);

    // Do not keep this here:
    void write_pdb_file(MolecularModeling::Assembly *glycoprotein, int cycle, std::string summary_filename, double overlap);


    //////////////////////////////////////////////////////////
    //                       MUTATOR                        //
    //////////////////////////////////////////////////////////

    void SetGlycoprotein(Assembly* glycoprotein);
    void SetGlycanName(std::string glycan_name);
    void SetResidueNumber(std::string residue_number);
    void SetResidue(Residue* residue);
    void SetGlycan(Assembly glycan);
    void SetOverlap(MoleculeType moleculeType, double overlap);
    void SetGlycanOverlap(double overlap);
    void SetProteinOverlap(double overlap);
    void SetSelfGlycanBeads(AtomVector *beads);
    void SetProteinBeads(AtomVector *beads);
    void SetOtherGlycanBeads(AtomVector *beads);
    void SetDefaultDihedralAnglesUsingMetadata();
    void SetRandomDihedralAnglesUsingMetadata();
    void SetRandomDihedralAnglesUsingMetadataForNthLinkage(int linkage_number);
    void ResetDihedralAngles();
    void SetOtherGlycosites(GlycosylationSiteVector &glycosites);

    //////////////////////////////////////////////////////////
    //                       DISPLAY FUNCTION               //
    //////////////////////////////////////////////////////////

    void PrintOverlaps();
    void Print(std::string type = "All");

    //////////////////////////////////////////////////////////
    //                       OPERATORS                      //
    //////////////////////////////////////////////////////////

    inline bool operator==(const GlycosylationSite &rhs) const
    {
        return rhs.residue_->GetId() == residue_->GetId();
    }

    inline bool operator!=(const GlycosylationSite &rhs) const
    {
        return residue_->GetId() == rhs.residue_->GetId();
    }

    inline bool operator<(const GlycosylationSite &rhs) const
    {
        return residue_->GetId() < rhs.residue_->GetId();
    }


    inline bool operator>(const GlycosylationSite &rhs) const
    {
        return residue_->GetId() > rhs.residue_->GetId();
    }

private:

    //////////////////////////////////////////////////////////
    //                  PRIVATE FUNCTIONS                   //
    //////////////////////////////////////////////////////////

    void Prepare_Glycans_For_Superimposition_To_Particular_Residue(std::string amino_acid_name);
    void Superimpose_Glycan_To_Glycosite(Residue *glycosite_residue);
    double CalculateTorsionAngle(AtomVector atoms);
    //void SetRotatableBonds(Residue *residue1, Residue *residue2);
    void Rename_Protein_Residue_To_GLYCAM_Nomenclature();
    void FigureOutResidueLinkagesInGlycan(Residue *from_this_residue1, Residue *to_this_residue2, ResidueLinkageVector *residue_linkages);
    void RecursivelyGetAllNeighboringResidues(Atom* current_atom, ResidueVector* neighbors);
    Atom* GetConnectingProteinAtom(std::string residue_name);
    void WiggleOneLinkage(Residue_linkage &linkage, int *output_pdb_id, double tolerance = 0.1, int interval = 5);
    double CalculateBeadOverlaps(MoleculeType moleculeType = ALL);
    //double Calculate_and_print_bead_overlaps();
    double CalculateAtomicOverlaps(MoleculeType moleculeType = ALL);
    double CalculateBeadOverlaps(AtomVector &atomsA, AtomVector &atomsB);

    //////////////////////////////////////////////////////////
    //                       ATTRIBUTES                     //
    //////////////////////////////////////////////////////////

    Assembly* glycoprotein_;
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
    GlycosylationSiteVector other_glycosites_;
};

#endif // GLYCOSYLATIONSITE_H
