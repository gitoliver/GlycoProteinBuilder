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
    ResidueVector residues {residue1_, residue2_};
    return residues;
}

RotatableDihedralVector Residue_linkage::GetRotatableDihedrals() const
{
    return rotatable_bonds_;
}

int Residue_linkage::GetNumberOfRotatableBonds()
{
    return rotatable_bonds_.size();
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////

void Residue_linkage::SetDefaultDihedralAnglesUsingMetadata()
{
    for (auto entry : rotatable_bonds_)
    {
        entry.SetDihedralAngleUsingMetadata();
    }
}

// Range should be inherent to each dihedral. Should add that to the class.
void Residue_linkage::SetRandomDihedralAnglesUsingMetadata()
{
    for (auto entry : rotatable_bonds_)
    {
        bool use_ranges = true;
        entry.SetDihedralAngleUsingMetadata(use_ranges);
    }
}

void Residue_linkage::SetCustomDihedralAngles(std::vector <double> dihedral_angles)
{
    if(dihedral_angles.size() == rotatable_bonds_.size())
    {
        std::vector <double>::iterator dihedral_angle = dihedral_angles.begin();
        for(RotatableDihedralVector::iterator rotatable_bond = rotatable_bonds_.begin(); rotatable_bond != rotatable_bonds_.end(); ++rotatable_bond)
        {
            rotatable_bond->SetDihedralAngle(*dihedral_angle);
            ++dihedral_angle;
        }
    }
    else
    {   // Really need to figure out this throwing exceptions lark.
        std::cout << "ERROR; attempted to set dihedral angles for set of dihedrals but with mismatching number of bonds to angles\n" << std::endl;
    }
}

void Residue_linkage::SetPreviousDihedralAngles()
{
    for(RotatableDihedralVector::iterator rotatable_bond = rotatable_bonds_.begin(); rotatable_bond != rotatable_bonds_.end(); ++rotatable_bond)
    {
        rotatable_bond->SetPreviousDihedralAngle();
    }
}



void Residue_linkage::SetRandomDihedralAngles()
{
    for(RotatableDihedralVector::iterator rotatable_bond = rotatable_bonds_.begin(); rotatable_bond != rotatable_bonds_.end(); ++rotatable_bond)
    {
        rotatable_bond->RandomizeDihedralAngle();
    }
}

void Residue_linkage::DetermineAtomsThatMove()
{
    for(RotatableDihedralVector::iterator rotatable_bond = rotatable_bonds_.begin(); rotatable_bond != rotatable_bonds_.end(); ++rotatable_bond)
    {
        rotatable_bond->DetermineAtomsThatMove();
    }
}

//////////////////////////////////////////////////////////
//                       DISPLAY FUNCTION               //
//////////////////////////////////////////////////////////

void Residue_linkage::Print()
{
    for(RotatableDihedralVector::iterator rotatable_bond = rotatable_bonds_.begin(); rotatable_bond != rotatable_bonds_.end(); ++rotatable_bond)
    {
        rotatable_bond->Print();
    }
}

//////////////////////////////////////////////////////////
//                       OPERATORS                      //
//////////////////////////////////////////////////////////

std::ostream& operator<<(std::ostream& os, const Residue_linkage& residue_linkage)
{
    RotatableDihedralVector rotatable_bonds = residue_linkage.GetRotatableDihedrals();
    for(RotatableDihedralVector::iterator rotatable_bond = rotatable_bonds.begin(); rotatable_bond != rotatable_bonds.end(); ++rotatable_bond)
    {
        os << (*rotatable_bond);
    }
    return os;
} // operator<<


//////////////////////////////////////////////////////////
//                    PRIVATE FUNCTIONS                 //
//////////////////////////////////////////////////////////

void Residue_linkage::InitializeClass(Residue *residue1, Residue *residue2)
{
    this->SetResidues(residue1, residue2);
    this->SetConnectionAtoms(residue1_, residue2_);
    rotatable_bonds_ = this->FindRotatableBondsConnectingResidues(connection_atom1_, connection_atom2_);
    gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector metadata = this->FindMetadata(connection_atom1_, connection_atom2_);
    this->AddMetadataToRotatableDihedrals(metadata);
}

RotatableDihedralVector Residue_linkage::FindRotatableBondsConnectingResidues(Atom *connection_atom1, Atom *connection_atom2)
{
    // Going to ignore tags etc.
    // Given two residues that are connected. Find connecting atoms.
    // Search neighbors other than connected atom. Ie search out in both directions, but remain within same residue.
    // Warning, residue may have fused cycles!
    // Will fail for non-protein residues without cycles. As don't have a non-rotatable bond to anchor from. Can code that later (and deal with branches from these residues).
    std::cout << "Finding rot bonds for " << connection_atom1->GetResidue()->GetId() << " and " << connection_atom2->GetResidue()->GetId() << "\n";

    AtomVector rotation_points = selection::FindRotationPoints(connection_atom1);
    AtomVector rotation_points2 = selection::FindRotationPoints(connection_atom2);
    // Need to reverse one of these, so when concatenated, they are ordered ok. This might not be ok.
    std::reverse(rotation_points2.begin(), rotation_points2.end());
    // Now concatenate:
    rotation_points.insert( rotation_points.begin(), rotation_points2.begin(), rotation_points2.end() );
    // Now that have a list of rotation points. Split into pairs and find rotatable bonds between them
    // BUT remember that first and last atom in list are just there to define dihedrals. Skip those
    bool found = false;
    AtomVector connecting_atoms = {connection_atom1, connection_atom2};
    for(int i = 0; i < rotation_points.size(); i = i+2)
    {
        Atom *rotation_point1 = rotation_points.at(i);
        Atom *rotation_point2 = rotation_points.at(i+1);

        found = false;
        connecting_atoms.clear();
        // std::cout << "Finding Path between " << rotation_point1->GetId() << " and " << rotation_point2->GetId() << ".\n";
        selection::FindPathBetweenTwoAtoms(rotation_point1, rotation_point2, &connecting_atoms, &found);
        selection::ClearAtomDescriptions(rotation_point1->GetResidue());
        selection::ClearAtomDescriptions(rotation_point2->GetResidue());
        // Find neighboring atoms needed to define dihedral. Pass in connecting atoms so don't find any of those.
        Atom *neighbor1 =  selection::FindCyclePointNeighbor(connecting_atoms, rotation_point1);
        Atom *neighbor2 =  selection::FindCyclePointNeighbor(connecting_atoms, rotation_point2);
        // Insert these neighbors into list of connecting atoms, at beginning and end of vector.
        // connecting_atoms gets populated as it falls out, so list is reveresed from what you'd expect
        connecting_atoms.insert(connecting_atoms.begin(), neighbor2);
        connecting_atoms.push_back(neighbor1);

        std::cout << "Updated Path between " << rotation_point1->GetId() << " and " << rotation_point2->GetId() << ":\n";
        for(AtomVector::iterator it1 = connecting_atoms.begin(); it1 != connecting_atoms.end(); ++it1)
        {
            Atom *atom = *it1;
            std::cout << atom->GetId() << "\n";
        }
    }
    RotatableDihedralVector rotatable_bonds = this->SplitAtomVectorIntoRotatableBonds(connecting_atoms);
    return rotatable_bonds;
}

RotatableDihedralVector Residue_linkage::SplitAtomVectorIntoRotatableBonds(AtomVector atoms)
{
    //Ok looking for sets of four atoms, but shifting along vector by one atom for each dihedral.
    // So four atoms will make one rotatable bond, five will make two bonds, six will make three etc.
    RotatableDihedralVector rotatable_bonds_generated;
    if(atoms.size() < 4)
    {
        std::cout << "ERROR; in Residue_linkage::SplitAtomVectorIntoRotatableBonds, not enough atoms in atom vector: " << atoms.size() << std::endl;
    }
    else
    {
        for(AtomVector::iterator it1 = atoms.begin(); it1 != (atoms.end()-3); ++it1)
        {
            Atom *atom1 = *it1;
            Atom *atom2 = *(it1+1);
            Atom *atom3 = *(it1+2);
            Atom *atom4 = *(it1+3);
            rotatable_bonds_generated.emplace_back(atom1, atom2, atom3, atom4);
        }
    }
    return rotatable_bonds_generated;
}

gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector Residue_linkage::FindMetadata(Atom *connection_atom1, Atom *connection_atom2)
{
    gmml::MolecularMetadata::GLYCAM::DihedralAngleDataContainer DihedralAngleMetadata;
    gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector matching_entries = DihedralAngleMetadata.GetEntriesForLinkage(connection_atom1, connection_atom2);
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
        int bond_number = int (entry.index_); // typecast to an int
        int vector_position = (bond_number - 1); // vectors start at 0.
        std::cout << "Adding to position: "<< vector_position << " in vector of size: " << rotatable_bonds_.size() << std::endl;
        if (vector_position <= rotatable_bonds_.size())
        {
            rotatable_bonds_.at(vector_position).AddMetadata(entry);
            std::cout << "Added " << entry.index_ << " = " << entry.default_angle_value_ << " to: \n";
            rotatable_bonds_.at(vector_position).Print();
        }
        else
        {
            std::cout << "Huge problem in residue_linkage.cpp AddMetadataToRotatableDihedrals. Tried to add metadata to a rotatable bond that does not exist.\n"
                         "Check both dihedralangledata metadata and Residue_linkage::FindRotatableBondsConnectingResidues." << std::endl;
        }
    }
}


void Residue_linkage::SetResidues(Residue *residue1, Residue *residue2)
{
    residue1_ = residue1;
    residue2_ = residue2;
}

void Residue_linkage::SetConnectionAtoms(Residue *residue1, Residue *residue2)
{
    AtomVector connecting_atoms;
    bool found = false;
    selection::FindAtomsConnectingResidues(residue1->GetAtoms().at(0), residue2, &connecting_atoms, &found);
    connection_atom1_ = connecting_atoms.at(0);
    connection_atom2_ = connecting_atoms.at(1);
}
