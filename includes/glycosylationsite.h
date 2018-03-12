#ifndef GLYCOSYLATIONSITE_H
#define GLYCOSYLATIONSITE_H

#include "gmml.hpp"
#include "overlap_record.h"
#include <iomanip> // For setting precision and formating in std::cout 
#include <algorithm> //  std::erase, std::remove
using namespace MolecularModeling;
class GlycosylationSite
{
public:
    //////////////////////////////////////////////////////////
    //                    TYPE DEFINITION                   //
    //////////////////////////////////////////////////////////

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
    double GetOverlap();
    double GetGlycanOverlap();
    double GetProteinOverlap();
    double GetChi1Value();
    double GetChi2Value();
    AtomVector GetSelfGlycanBeads();
    AtomVector GetProteinBeads();
    AtomVector GetOtherGlycanBeads();
    Overlap_record GetBestOverlapRecord();
    Overlap_record GetBestProteinOverlapRecord();


    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////
    void AttachGlycan(Assembly glycan, Assembly *glycoprotein);
    double Calculate_bead_overlaps();
    double Calculate_protein_bead_overlaps();
    double Calculate_other_glycan_bead_overlaps();
    double Calculate_and_print_bead_overlaps();
    void SetChiAtoms(Residue* residue);

    //////////////////////////////////////////////////////////
    //                       MUTATOR                        //
    //////////////////////////////////////////////////////////

    void SetGlycanName(std::string glycan_name);
    void SetResidueNumber(std::string residue_number);
    void SetResidue(Residue* residue);
    void SetGlycan(Assembly glycan);
    void SetGlycanOverlap(double overlap);
    void SetProteinOverlap(double overlap);
    void SetChi1Value(double angle, Assembly *glycoprotein);
    void SetChi2Value(double angle, Assembly *glycoprotein);
    void SetSelfGlycanBeads(AtomVector *beads);
    void SetProteinBeads(AtomVector *beads);
    void SetOtherGlycanBeads(AtomVector *beads);
    void SetBestOverlapRecord(double overlap, double chi1, double chi2);
    void SetBestProteinOverlapRecord(double overlap, double chi1, double chi2);

    //////////////////////////////////////////////////////////
    //                       DISPLAY FUNCTION               //
    //////////////////////////////////////////////////////////

    void Print_bead_overlaps();

private:

    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////

    void Prepare_Glycans_For_Superimposition_To_Particular_Residue(std::string amino_acid_name);
    void Superimpose_Glycan_To_Glycosite(Residue *glycosite_residue);
    double CalculateTorsionAngle(AtomVector atoms);
    double Calculate_bead_overlaps(AtomVector &atomsA, AtomVector &atomsB);

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
    AtomVector chi1_;
    AtomVector chi2_;
    AtomVector self_glycan_beads_;
    AtomVector other_glycan_beads_;
    AtomVector protein_beads_;
    Overlap_record best_overlap_record_;
    Overlap_record best_protein_overlap_record_;
};

#endif // GLYCOSYLATIONSITE_H
