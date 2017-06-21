#include "overlap.h"

//////////////////////////////////////////////////////////
//                    TYPE DEFINITION                   //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////

Overlap::Overlap()
{
    overlap_ = 0.0;
}

Overlap::Overlap(Residue *residue1, Residue *residue2, double overlap)
{
    residue1_ = residue1;
    residue2_ = residue2;
    overlap_ = overlap;
}

//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////

double Overlap::GetOverlap()
{
    return overlap_;
}

Residue* Overlap::GetResidue1()
{
    return residue1_;
}

Residue* Overlap::GetResidue2()
{
    return residue2_;
}

//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void Overlap::SetOverlap(double overlap)
{
    overlap_ = overlap;
}

void Overlap::SetResidue1(Residue *residue1)
{
    residue1_ = residue1;
}

void Overlap::SetResidue2(Residue *residue2)
{
    residue2_ = residue2;
}

//////////////////////////////////////////////////////////
//                       DISPLAY FUNCTION               //
//////////////////////////////////////////////////////////

void Overlap::Print(std::ostream& out)
{
    out << residue1_->GetId();
    out << " has a " << overlap_;
    out << " overlap with " << residue2_->GetId() << std::endl;
}




