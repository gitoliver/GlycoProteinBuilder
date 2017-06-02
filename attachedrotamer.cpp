#include "attachedrotamer.h"


//////////////////////////////////////////////////////////
//                    TYPE DEFINITION                   //
//////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////

AttachedRotamer::AttachedRotamer()
{
    glycan_rotamer_ = NULL;
}

AttachedRotamer::AttachedRotamer(Assembly *glycan_rotamer)
{
    glycan_rotamer_ = glycan_rotamer;
}

//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////

Assembly* AttachedRotamer::GetAttachedRotamer()
{
    return glycan_rotamer_;
}

std::string AttachedRotamer::GetGlycanPDBPath()
{
    return glycan_pdb_path_;
}

double AttachedRotamer::GetTotalOverlap() // Recalculate and return double each time, do not store.
{
    double total_overlap = 0.0;
    for (OverlapVector::iterator it = overlaps_.begin(); it != overlaps_.end(); ++it)
    {
        Overlap *overlap = (*it);
        total_overlap += overlap->GetOverlap();
    }
    return total_overlap;
}

//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void AttachedRotamer::SetGlycanPDBPath(std::string glycan_pdb_path)
{
    glycan_pdb_path_ = glycan_pdb_path;
}

//////////////////////////////////////////////////////////
//                       DISPLAY FUNCTION               //
//////////////////////////////////////////////////////////

//void Print(std::ostream& out = std::cout) {
//    std::out << "Residue ID: " << total_overlap_ << endl;
//}



