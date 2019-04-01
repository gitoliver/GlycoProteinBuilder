#include "../includes/residue_linkage.h"

//////////////////////////////////////////////////////////
//                    TYPE DEFINITION                   //
//////////////////////////////////////////////////////////

typedef std::vector<Rotatable_dihedral> RotatableDihedralVector;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////

Residue_linkage::Residue_linkage() {} // Do nothin

Residue_linkage::Residue_linkage(Residue *residue1, Residue *residue2)
{
    this->InitializeClass(residue1, residue2);
}

//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////

ResidueVector Residue_linkage::GetResidues()
{
    ResidueVector residues {from_this_residue1_, to_this_residue2_};
    return residues;
}

RotatableDihedralVector Residue_linkage::GetRotatableDihedrals() const
{
    return rotatable_dihedrals_;
}

int Residue_linkage::GetNumberOfRotatableDihedrals()
{
    return rotatable_dihedrals_.size();
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////

void Residue_linkage::SetDefaultDihedralAnglesUsingMetadata()
{
    for (auto &entry : rotatable_dihedrals_)
    {
        entry.SetDihedralAngleUsingMetadata();
    }
}

// Range should be inherent to each dihedral. Should add that to the class.
void Residue_linkage::SetRandomDihedralAnglesUsingMetadata()
{
    for (auto &entry : rotatable_dihedrals_)
    {
        bool use_ranges = true;
        entry.SetDihedralAngleUsingMetadata(use_ranges);
    }
}

void Residue_linkage::SetCustomDihedralAngles(std::vector <double> dihedral_angles)
{
    if(dihedral_angles.size() == rotatable_dihedrals_.size())
    {
        std::vector <double>::iterator dihedral_angle_iterator = dihedral_angles.begin();
        for (auto &rotatable_dihedral : rotatable_dihedrals_)
        {
            rotatable_dihedral.SetDihedralAngle(*dihedral_angle_iterator);
            ++dihedral_angle_iterator;
        }
    }
    else
    {   // Really need to figure out this throwing exceptions lark.
        std::cout << "ERROR; attempted to set dihedral angles for set of dihedrals but with mismatching number of bonds to angles\n" << std::endl;
    }
}

void Residue_linkage::SetDihedralAnglesToPrevious()
{
    for(RotatableDihedralVector::iterator rotatable_dihedral = rotatable_dihedrals_.begin(); rotatable_dihedral != rotatable_dihedrals_.end(); ++rotatable_dihedral)
    {
        rotatable_dihedral->SetDihedralAngleToPrevious();
    }
}


void Residue_linkage::SetRandomDihedralAngles()
{
    for(RotatableDihedralVector::iterator rotatable_dihedral = rotatable_dihedrals_.begin(); rotatable_dihedral != rotatable_dihedrals_.end(); ++rotatable_dihedral)
    {
        rotatable_dihedral->RandomizeDihedralAngle();
    }
}

void Residue_linkage::DetermineAtomsThatMove()
{
    for(RotatableDihedralVector::iterator rotatable_dihedral = rotatable_dihedrals_.begin(); rotatable_dihedral != rotatable_dihedrals_.end(); ++rotatable_dihedral)
    {
        rotatable_dihedral->DetermineAtomsThatMove();
    }
}

//////////////////////////////////////////////////////////
//                       DISPLAY FUNCTION               //
//////////////////////////////////////////////////////////

void Residue_linkage::Print()
{
    for(RotatableDihedralVector::iterator rotatable_dihedral = rotatable_dihedrals_.begin(); rotatable_dihedral != rotatable_dihedrals_.end(); ++rotatable_dihedral)
    {
        rotatable_dihedral->Print();
    }
}

//////////////////////////////////////////////////////////
//                       OPERATORS                      //
//////////////////////////////////////////////////////////

std::ostream& operator<<(std::ostream& os, const Residue_linkage& residue_linkage)
{
    RotatableDihedralVector rotatable_dihedrals = residue_linkage.GetRotatableDihedrals();
    for(RotatableDihedralVector::iterator rotatable_dihedral = rotatable_dihedrals.begin(); rotatable_dihedral != rotatable_dihedrals.end(); ++rotatable_dihedral)
    {
        os << (*rotatable_dihedral);
    }
    return os;
} // operator<<


//////////////////////////////////////////////////////////
//                    PRIVATE FUNCTIONS                 //
//////////////////////////////////////////////////////////

void Residue_linkage::InitializeClass(Residue *from_this_residue1, Residue *to_this_residue2)
{
    this->SetResidues(from_this_residue1, to_this_residue2);
    //std::cout << "Finding connection between " << from_this_residue1->GetId() << " :: " << to_this_residue2->GetId() << std::endl;
    this->SetConnectionAtoms(from_this_residue1_, to_this_residue2_);
    rotatable_dihedrals_ = this->FindRotatableDihedralsConnectingResidues(from_this_connection_atom1_, to_this_connection_atom2_);
    gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector metadata = this->FindMetadata(from_this_connection_atom1_, to_this_connection_atom2_);
    this->AddMetadataToRotatableDihedrals(metadata);
}

RotatableDihedralVector Residue_linkage::FindRotatableDihedralsConnectingResidues(Atom *from_this_connection_atom1, Atom *to_this_connection_atom2)
{
    // Going to ignore tags etc.
    // Given two residues that are connected. Find connecting atoms.
    // Search neighbors other than connected atom. Ie search out in both directions, but remain within same residue.
    // Warning, residue may have fused cycles!
    // Will fail for non-protein residues without cycles. As don't have a non-rotatable bond to anchor from. Can code that later (and deal with branches from these residues).
    //std::cout << "Finding rot bonds for " << from_this_connection_atom1->GetResidue()->GetId() << " and " << to_this_connection_atom2->GetResidue()->GetId() << "\n";

    AtomVector from_this_residue1_cycle_points = selection::FindCyclePoints(from_this_connection_atom1);
   // std::cout << "Moving onto second residue.\n";
    AtomVector to_this_residue2_cycle_points = selection::FindCyclePoints(to_this_connection_atom2);
    // Need to reverse one of these, so when concatenated, they are ordered ok. This might not be ok.
//    std::reverse(to_this_residue2_cycle_points.begin(), to_this_residue2_cycle_points.end());
    std::reverse(from_this_residue1_cycle_points.begin(), from_this_residue1_cycle_points.end());
    // Now concatenate:
    from_this_residue1_cycle_points.insert( from_this_residue1_cycle_points.end(), to_this_residue2_cycle_points.begin(), to_this_residue2_cycle_points.end() );
    // Now that have a list of rotation points. Split into pairs and find rotatable bonds between them
    bool found = false;
    AtomVector connecting_atoms = {from_this_connection_atom1, to_this_connection_atom2};
    for(int i = 0; i < from_this_residue1_cycle_points.size(); i = i+2)
    {
        Atom *cycle_point1 = from_this_residue1_cycle_points.at(i);
        Atom *cycle_point2 = from_this_residue1_cycle_points.at(i+1);

        found = false;
        connecting_atoms.clear();
        //std::cout << "Finding Path between " << rotation_point1->GetId() << " and " << rotation_point2->GetId() << ".\n";
        selection::FindPathBetweenTwoAtoms(cycle_point1, cycle_point2, &connecting_atoms, &found);
        selection::ClearAtomDescriptions(cycle_point1->GetResidue());
        selection::ClearAtomDescriptions(cycle_point2->GetResidue());
        // Find neighboring atoms needed to define dihedral. Pass in connecting atoms so don't find any of those.
        Atom *neighbor1 =  selection::FindCyclePointNeighbor(connecting_atoms, cycle_point1);
        Atom *neighbor2 =  selection::FindCyclePointNeighbor(connecting_atoms, cycle_point2);
        // Insert these neighbors into list of connecting atoms, at beginning and end of vector.
        // connecting_atoms gets populated as it falls out, so list is reveresed from what you'd expect
        std::reverse(connecting_atoms.begin(), connecting_atoms.end());
        connecting_atoms.insert(connecting_atoms.begin(), neighbor1);
        connecting_atoms.push_back(neighbor2);

//        std::cout << "Updated Path between " << cycle_point1->GetId() << " and " << cycle_point2->GetId() << ":\n";
//        for(AtomVector::iterator it1 = connecting_atoms.begin(); it1 != connecting_atoms.end(); ++it1)
//        {
//            Atom *atom = *it1;
//            std::cout << atom->GetId() << "\n";
//        }
//        for (const auto& atom : connecting_atoms)
//        {
//            std::cout << atom->GetId() << "\n";
//        }
    }
    RotatableDihedralVector rotatable_dihedrals = this->SplitAtomVectorIntoRotatableDihedrals(connecting_atoms);
    return rotatable_dihedrals;
}

RotatableDihedralVector Residue_linkage::SplitAtomVectorIntoRotatableDihedrals(AtomVector atoms)
{
    //Ok looking for sets of four atoms, but shifting along vector by one atom for each dihedral.
    // So four atoms will make one rotatable bond, five will make two bonds, six will make three etc.
    RotatableDihedralVector rotatable_dihedrals_generated;
    if(atoms.size() < 4)
    {
        std::cout << "ERROR; in Residue_linkage::SplitAtomVectorIntoRotatableDihedrals, not enough atoms in atom vector: " << atoms.size() << std::endl;
    }
    else
    {
        for(AtomVector::iterator it1 = atoms.begin(); it1 != (atoms.end()-3); ++it1)
        {
            Atom *atom1 = *it1;
            Atom *atom2 = *(it1+1);
            Atom *atom3 = *(it1+2);
            Atom *atom4 = *(it1+3);
            rotatable_dihedrals_generated.emplace_back(atom1, atom2, atom3, atom4);
        }
    }
    return rotatable_dihedrals_generated;
}

gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector Residue_linkage::FindMetadata(Atom *from_this_connection_atom1, Atom *to_this_connection_atom2)
{
    gmml::MolecularMetadata::GLYCAM::DihedralAngleDataContainer DihedralAngleMetadata;
    gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector matching_entries = DihedralAngleMetadata.GetEntriesForLinkage(from_this_connection_atom1, to_this_connection_atom2);
    std::cout << "Found these " << matching_entries.size() << " entries:\n";
    for (const auto& entry : matching_entries)
    {
        std::cout << entry.index_ << " : " << entry.atom1_ << ", " << entry.atom2_ << ", " << entry.atom3_ << ", " << entry.atom4_ << ", " << entry.default_angle_value_ << "\n";
    }
    return matching_entries;
}

void Residue_linkage::AddMetadataToRotatableDihedrals(gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector metadata)
{
    for (const auto& entry : metadata)
    {
//        int bond_number = int (entry.number_of_bonds_from_anomeric_carbon_); // typecast to an int
        int vector_position = (entry.number_of_bonds_from_anomeric_carbon_ - 1); // vectors start at 0.
        std::cout << "Adding to position: "<< vector_position << " in vector of size: " << rotatable_dihedrals_.size() << std::endl;
        if (vector_position <= rotatable_dihedrals_.size())
        {
            // I think that here I need to check for conformers? No wait that's being handled already? Hmm..
            rotatable_dihedrals_.at(vector_position).AddMetadata(entry);
            std::cout << "Added " << entry.index_ << " = " << entry.default_angle_value_ << " to: \n";
            rotatable_dihedrals_.at(vector_position).Print();
        }
        else
        {
            std::cout << "Huge problem in residue_linkage.cpp AddMetadataToRotatableDihedrals. Tried to add metadata to a rotatable bond that does not exist.\n"
                         "Check both dihedralangledata metadata and Residue_linkage::FindRotatableDihedralsConnectingResidues." << std::endl;
        }
    }
}


void Residue_linkage::SetResidues(Residue *residue1, Residue *residue2)
{
    from_this_residue1_ = residue1;
    to_this_residue2_ = residue2;
}

void Residue_linkage::SetConnectionAtoms(Residue *residue1, Residue *residue2)
{
    AtomVector connecting_atoms;
    bool found = false;
    selection::FindAtomsConnectingResidues(residue1->GetAtoms().at(0), residue2, &connecting_atoms, &found);
    from_this_connection_atom1_ = connecting_atoms.at(0);
    to_this_connection_atom2_ = connecting_atoms.at(1);
}
