#ifndef GLYCOSYLATIONSITE_H
#define GLYCOSYLATIONSITE_H

//#include <iostream>
//#include "/home/oliver/Programs/gems/gmml/includes/gmml.hpp"
#include "attachedrotamer.h"

class GlycosylationSite
{
public:
    //////////////////////////////////////////////////////////
    //                    TYPE DEFINITION                   //
    //////////////////////////////////////////////////////////

    typedef std::vector<AttachedRotamer*> AttachedRotamerVector;

    //////////////////////////////////////////////////////////
    //                       CONSTRUCTOR                    //
    //////////////////////////////////////////////////////////

    GlycosylationSite();
    GlycosylationSite(std::string glycan_name);
    GlycosylationSite(std::string glycan_name, Residue* residue, AttachedRotamerVector attached_rotamers);
    ~GlycosylationSite();
    //////////////////////////////////////////////////////////
    //                       ACCESSOR                       //
    //////////////////////////////////////////////////////////

    std::string GetGlycanName();
    Residue* GetResidue();
    AttachedRotamerVector GetAttachedRotamers();

    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////
    //                       MUTATOR                        //
    //////////////////////////////////////////////////////////

    void SetGlycanName(std::string glycan_name);
    void SetResidue(Residue* residue);
    void SetAttachedRotamers(AttachedRotamerVector attached_rotamers);
    void AddRotamer(AttachedRotamer *rotamer);

    //////////////////////////////////////////////////////////
    //                       DISPLAY FUNCTION               //
    //////////////////////////////////////////////////////////

    //void Print(std::ostream& out = std::cout);

private:
    //////////////////////////////////////////////////////////
    //                       ATTRIBUTES                     //
    //////////////////////////////////////////////////////////

    std::string glycan_name_;
    Residue* residue_;                                  /*!< A pointer back to the residue for this glycosite >*/
    AttachedRotamerVector attached_rotamers_;

};

#endif // GLYCOSYLATIONSITE_H
