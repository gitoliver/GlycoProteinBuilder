#ifndef GLYCOSYLATIONSITE_H
#define GLYCOSYLATIONSITE_H

//#include <iostream>
//#include "/home/ubunter/software/gems/gmml/includes/gmml.hpp"
#include "../../gems/gmml/includes/gmml.hpp"

//#include "attachedglycan.h"



class GlycosylationSite
{
public:
    //////////////////////////////////////////////////////////
    //                    TYPE DEFINITION                   //
    //////////////////////////////////////////////////////////
    //typedef std::vector<Atom*> AtomVector;

    //////////////////////////////////////////////////////////
    //                       CONSTRUCTOR                    //
    //////////////////////////////////////////////////////////

    GlycosylationSite();
    GlycosylationSite(std::string glycan_name);
    GlycosylationSite(std::string glycan_name, Residue* residue, Assembly glycan);
    ~GlycosylationSite();
    //////////////////////////////////////////////////////////
    //                       ACCESSOR                       //
    //////////////////////////////////////////////////////////

    std::string GetGlycanName();
    Residue* GetResidue();
    Assembly* GetAttachedGlycan();
    double GetTotalOverlap();
    double GetGlycanOverlap();
    double GetProteinOverlap();
    double GetChi1Value();
    double GetChi2Value();
    AtomVector GetSelfGlycanBeads();
    AtomVector GetProteinBeads();
    AtomVector GetOtherGlycanBeads();

    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////
    void AttachGlycan(Assembly glycan, Assembly *glycoprotein);

    //double calculate_overlaps(Atomvector all_atoms);
    //void calculate_protein_overlap(Atomvector );
    double calculate_bead_overlaps();
    void SetChiAtoms(Residue* residue);

    //////////////////////////////////////////////////////////
    //                       MUTATOR                        //
    //////////////////////////////////////////////////////////

    void SetGlycanName(std::string glycan_name);
    void SetResidue(Residue* residue);
    void SetGlycan(Assembly glycan);
    void SetGlycanOverlap(double overlap);
    void SetProteinOverlap(double overlap);
    void SetChi1Value(double angle, Assembly *glycoprotein);
    void SetChi2Value(double angle, Assembly *glycoprotein);
    void SetSelfGlycanBeads(AtomVector *beads);
    void SetProteinBeads(AtomVector *beads);
    void SetOtherGlycanBeads(AtomVector *beads);


    //////////////////////////////////////////////////////////
    //                       DISPLAY FUNCTION               //
    //////////////////////////////////////////////////////////

    //void Print(std::ostream& out = std::cout);

private:

    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////

    void Prepare_Glycans_For_Superimposition_To_Particular_Residue(std::string amino_acid_name);
    void Superimpose_Glycan_To_Glycosite(Residue *glycosite_residue);

    //////////////////////////////////////////////////////////
    //                       ATTRIBUTES                     //
    //////////////////////////////////////////////////////////

    std::string glycan_name_;
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

    //Assembly alternate_sidechain_;

};

#endif // GLYCOSYLATIONSITE_H
