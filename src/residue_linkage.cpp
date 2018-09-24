#include "../includes/residue_linkage.h"


//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////

Residue_linkage::Residue_linkage()
{

}

Residue_linkage::Residue_linkage(Residue *residue1, Residue *residue2)
{
    this->SetResidues(residue1, residue2);
    this->RandomizeDihedralAngles();
}


Residue_linkage::Residue_linkage(Residue *residue1, Residue *residue2, std::vector <double> dihedral_angles)
{
    this->SetResidues(residue1, residue2);
    this->SetDihedralAngles(dihedral_angles);
}


//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////

ResidueVector Residue_linkage::GetResidues()
{
    ResidueVector residues {residue1_, residue2_};
    return residues;
}

RotatableDihedralVector Residue_linkage::GetRotatableDihedrals()
{
    return rotatable_bonds_;
}

//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void Residue_linkage::SetResidues(Residue *residue1, Residue *residue2)
{
    residue1_ = residue1;
    residue2_ = residue2;
}

void Residue_linkage::SetDihedralAngles(std::vector <double> dihedral_angles)
{
    if(dihedral_angles.size() == rotatable_bonds_.size())
    {
        std::vector <double>::iterator it0 = dihedral_angles.begin();
        for(RotatableDihedralVector::iterator it1 = rotatable_bonds_.begin(); it1 != rotatable_bonds_.end(); ++it1)
        {
            double *dihedral_angle = *it0;
            Rotatable_dihedral *rotatable_bond = *it1;
            rotatable_bond->SetDihedralAngle(dihedral_angle);
            ++it0;
        }
    }
    else
    {
        // Really need to figure out this throwing exceptions lark.
        std::cout << "ERROR; attempted to set dihedral angles for set of dihedrals but with mismatching number of bonds to angles\n" << std::endl;
    }
}

void Residue_linkage::ResetDihedralAngles()
{

}

double Residue_linkage::RandomizeDihedralAnglesWithinTheirRanges()
{

}

double Residue_linkage::RandomizeDihedralAngles()
{

}

//////////////////////////////////////////////////////////
//                       DISPLAY FUNCTION               //
//////////////////////////////////////////////////////////

//  void Residue_linkage::Print();

//////////////////////////////////////////////////////////
//                       OPERATORS                      //
//////////////////////////////////////////////////////////



