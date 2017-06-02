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
    GlycosylationSite(Residue* residue, AttachedRotamerVector attached_rotamers);
    ~GlycosylationSite();
    //////////////////////////////////////////////////////////
    //                       ACCESSOR                       //
    //////////////////////////////////////////////////////////

    Residue* GetResidue();
    AttachedRotamerVector GetAttachedRotamers();

    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////
    //                       MUTATOR                        //
    //////////////////////////////////////////////////////////

    void SetResidue(Residue* residue);
    void SetAttachedRotamers(AttachedRotamerVector attached_rotamers);

    //////////////////////////////////////////////////////////
    //                       DISPLAY FUNCTION               //
    //////////////////////////////////////////////////////////

    //void Print(std::ostream& out = std::cout);

private:
    //////////////////////////////////////////////////////////
    //                       ATTRIBUTES                     //
    //////////////////////////////////////////////////////////

    Residue* residue_;                                  /*!< A pointer back to the residue for this glycosite >*/
    AttachedRotamerVector attached_rotamers_;

};

#endif // GLYCOSYLATIONSITE_H
