#include "glycosylationsite.h"

typedef std::vector<AttachedRotamer*> AttachedRotamerVector; // Had to define it here to use it as return type.

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
GlycosylationSite::GlycosylationSite()
{
    residue_ = NULL;
    attached_rotamers_ = AttachedRotamerVector();
}

GlycosylationSite::GlycosylationSite(Residue* residue, AttachedRotamerVector attached_rotamers)
{
    residue_ = residue;
    attached_rotamers_ = attached_rotamers;
}


//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////

Residue* GlycosylationSite::GetResidue()
{
    return residue_;
}

AttachedRotamerVector GlycosylationSite::GetAttachedRotamers()
{
    return attached_rotamers_;
}

//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void GlycosylationSite::SetResidue(Residue* residue)
{
    residue_ = residue;
}

void GlycosylationSite::SetAttachedRotamers(AttachedRotamerVector attached_rotamers)
{
    attached_rotamers_ = attached_rotamers;
}

//////////////////////////////////////////////////////////
//                       DISPLAY FUNCTION               //
//////////////////////////////////////////////////////////

/*void GlycosylationSite::Print(ostream &out)
{
    std::out << "Residue ID: " << residue_->GetId() << endl;
    //out << "Glycan sequence: " << this->GetGlycanSequence() << endl;
}
*/
