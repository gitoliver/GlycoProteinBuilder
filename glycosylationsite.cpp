#include "glycosylationsite.h"

typedef std::vector<AttachedRotamer*> AttachedRotamerVector; // Had to define it here to use it as return type.

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
GlycosylationSite::GlycosylationSite()
{
    glycan_name_ = "";
    residue_ = NULL;
    attached_rotamers_ = AttachedRotamerVector();
}

GlycosylationSite::GlycosylationSite(std::string glycan_name)
{
    glycan_name_ = glycan_name;
}

GlycosylationSite::GlycosylationSite(std::string glycan_name, Residue* residue, AttachedRotamerVector attached_rotamers)
{
    glycan_name_ = glycan_name;
    residue_ = residue;
    attached_rotamers_ = attached_rotamers;
}


//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////

std::string GlycosylationSite::GetGlycanName()
{
    return glycan_name_;
}

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

void GlycosylationSite::SetGlycanName(std::string glycan_name)
{
    glycan_name_ = glycan_name;
}

void GlycosylationSite::SetResidue(Residue* residue)
{
    residue_ = residue;
}

void GlycosylationSite::SetAttachedRotamers(AttachedRotamerVector attached_rotamers)
{
    attached_rotamers_ = attached_rotamers;
}

void GlycosylationSite::AddRotamer(AttachedRotamer *rotamer)
{
    attached_rotamers_.push_back(rotamer);
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
