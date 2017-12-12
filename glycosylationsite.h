#ifndef GLYCOSYLATIONSITE_H
#define GLYCOSYLATIONSITE_H

//#include <iostream>
//#include "/home/ubunter/software/gems/gmml/includes/gmml.hpp"
#include "../../../gems/gmml/includes/gmml.hpp"
//#include "attachedglycan.h"

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
    GlycosylationSite(std::string glycan_name, Residue* residue, Assembly glycan);
    ~GlycosylationSite();
    //////////////////////////////////////////////////////////
    //                       ACCESSOR                       //
    //////////////////////////////////////////////////////////

    std::string GetGlycanName();
    Residue* GetResidue();
    Assembly* GetAttachedGlycan();

    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////
    void AttachGlycan(Assembly glycan, Assembly *glycoprotein);
    void Prepare_Glycans_For_Superimposition_To_Particular_Residue(std::string amino_acid_name);
    void Superimpose_Glycan_To_Glycosite(Residue *glycosite_residue);

    //////////////////////////////////////////////////////////
    //                       MUTATOR                        //
    //////////////////////////////////////////////////////////

    void SetGlycanName(std::string glycan_name);
    void SetResidue(Residue* residue);
    void SetGlycan(Assembly glycan);

    //////////////////////////////////////////////////////////
    //                       DISPLAY FUNCTION               //
    //////////////////////////////////////////////////////////

    //void Print(std::ostream& out = std::cout);

private:

    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////
    //                       ATTRIBUTES                     //
    //////////////////////////////////////////////////////////

    std::string glycan_name_;
    Residue* residue_;                                  /*!< A pointer back to the residue for this glycosite >*/
    Assembly glycan_;
    AtomVector superimposition_atoms_;               /*!< The 3 atoms used for superimposition of glycan to sidechain >*/
    //Assembly alternate_sidechain_;

};

#endif // GLYCOSYLATIONSITE_H
