#include "../includes/rotatable_dihedral.h"
#include <random>

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
rotatable_dihedral::rotatable_dihedral()
{
    this->angle_ = gmml::dNotSet;
}

rotatable_dihedral::rotatable_dihedral(Atom *atom1, Atom *atom2, Atom *atom3, Atom *atom4)
{
    atom1_ = atom1;
    atom2_ = atom2;
    atom3_ = atom3;
    atom4_ = atom4;
    AtomVector atoms_that_move;
    atoms_that_move.push_back(atom2_);
    atom3_->FindConnectedAtoms(atoms_that_move);
    atoms_that_move_ = atoms_that_move;
    angle_ = this->CalculateDihedral();
}

rotatable_dihedral::rotatable_dihedral(AtomVector atoms)
{
    atom1_ = atoms.at(0);
    atom2_ = atoms.at(1);
    atom3_ = atoms.at(2);
    atom4_ = atoms.at(3);
    AtomVector atoms_that_move;
    atoms_that_move.push_back(atom2_);
    atom3_->FindConnectedAtoms(atoms_that_move);
    atoms_that_move_ = atoms_that_move;
    angle_ = this->CalculateDihedral();
}

//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////

double rotatable_dihedral::CalculateAngle()
{
    GeometryTopology::Coordinate* a1 = atom1_->GetCoordinate();
    GeometryTopology::Coordinate* a2 = atom2_->GetCoordinate();
    GeometryTopology::Coordinate* a3 = atom3_->GetCoordinate();
    GeometryTopology::Coordinate* a4 = atom4_->GetCoordinate();

    GeometryTopology::Coordinate b1 = a2;
    b1.operator -(*a1);
    GeometryTopology::Coordinate b2 = a3;
    b2.operator -(*a2);
    GeometryTopology::Coordinate b3 = a4;
    b3.operator -(*a3);
    GeometryTopology::Coordinate b4 = b2;
    b4.operator *(-1);

    GeometryTopology::Coordinate b2xb3 = b2;
    b2xb3.CrossProduct(b3);

    GeometryTopology::Coordinate b1_m_b2n = b1;
    b1_m_b2n.operator *(b2.length());

    GeometryTopology::Coordinate b1xb2 = b1;
    b1xb2.CrossProduct(b2);

    double current_dihedral = atan2(b1_m_b2n.DotProduct(b2xb3), b1xb2.DotProduct(b2xb3));

    return current_dihedral;
}

AtomVector rotatable_dihedral::GetAtoms()
{

}
AtomVector rotatable_dihedral::GetAtomsThatMove()
{

}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void rotatable_dihedral::SetAtoms(AtomVector atoms)
{

}

void rotatable_dihedral::SetAtomsThatMove(AtomVector atoms)
{

}

double rotatable_dihedral::CalculateAngle()
{

}


double rotatable_dihedral::RandomizeAngle()
{
    double random_angle = rotatable_dihedral::RandomizeAngleWithinRange(0.0, 360.0);
    this->SetAngle(random_angle);
    return random_angle;
    //return (rand() % 360) + 1 - 180; // Can get same one everytime for testing
}

double rotatable_dihedral::RandomizeAngleWithinRange(double min, double max)
{
    std::random_device rd1; // obtain a random number from hardware
    std::mt19937 eng1(rd1()); // seed the generator
    std::uniform_real_distribution<> angle_distribution(min, max); // define the range

    return angle_distribution(eng1);
    //return rand() % (max + 1 - min) + min; // Can get same one everytime for testing
}

double rotatable_dihedral::RandomizeAngleWithinRanges(std::vector<std::pair<double,double>> ranges)
{
    // For usage, can do ranges.emplace_back(min, max);
    // Pass in a vector of pairs of ranges.
    // First select one of those ranges.
    // Then create an angle within the selected range.

    // Rando stuff from slack overflow:
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 eng(rd()); // seed the generator
    std::uniform_int_distribution<> distr(0, (ranges.size() - 1)); // define the range

    // Select one of the ranges
    int range_selection = distr(eng);
    std::pair *selected_range = ranges.at(range_selection);

    // create an angle within the selected range
    return rotatable_dihedral::RandomizeAngleWithinRange(selected_range->first, selected_range->second);
}

void rotatable_dihedral::SetDihedral(double torsion, AtomVector &atomsToRotate)
{
    GeometryTopology::Coordinate* a1 = atom1_->GetCoordinate();
    GeometryTopology::Coordinate* a2 = atom2_->GetCoordinate();
    GeometryTopology::Coordinate* a3 = atom3_->GetCoordinate();
    GeometryTopology::Coordinate* a4 = atom4_->GetCoordinate();

    GeometryTopology::Coordinate b1 = a2;
    b1.operator -(*a1);
    GeometryTopology::Coordinate b2 = a3;
    b2.operator -(*a2);
    GeometryTopology::Coordinate b3 = a4;
    b3.operator -(*a3);
    GeometryTopology::Coordinate b4 = b2;
    b4.operator *(-1);

    GeometryTopology::Coordinate b2xb3 = b2;
    b2xb3.CrossProduct(b3);

    GeometryTopology::Coordinate b1_m_b2n = b1;
    b1_m_b2n.operator *(b2.length());

    GeometryTopology::Coordinate b1xb2 = b1;
    b1xb2.CrossProduct(b2);

    double current_dihedral = atan2(b1_m_b2n.DotProduct(b2xb3), b1xb2.DotProduct(b2xb3));

    double** torsion_matrix = gmml::GenerateRotationMatrix(&b4, a2, current_dihedral - gmml::ConvertDegree2Radian(torsion));

//    AtomVector atomsToRotate = AtomVector();
//    atomsToRotate.push_back(atom2);
//    atom3->FindConnectedAtoms(atomsToRotate);
 //   std::cout << "Moving: ";
    for(AtomVector::iterator it = atomsToRotate.begin(); it != atomsToRotate.end(); it++)
    {
        Atom *atom = *it;
        GeometryTopology::Coordinate* atom_coordinate = atom->GetCoordinate();
        GeometryTopology::Coordinate result;
        result.SetX(torsion_matrix[0][0] * atom_coordinate->GetX() + torsion_matrix[0][1] * atom_coordinate->GetY() +
                torsion_matrix[0][2] * atom_coordinate->GetZ() + torsion_matrix[0][3]);
        result.SetY(torsion_matrix[1][0] * atom_coordinate->GetX() + torsion_matrix[1][1] * atom_coordinate->GetY() +
                torsion_matrix[1][2] * atom_coordinate->GetZ() + torsion_matrix[1][3]);
        result.SetZ(torsion_matrix[2][0] * atom_coordinate->GetX() + torsion_matrix[2][1] * atom_coordinate->GetY() +
                torsion_matrix[2][2] * atom_coordinate->GetZ() + torsion_matrix[2][3]);

        atom->GetCoordinate()->SetX(result.GetX());
        atom->GetCoordinate()->SetY(result.GetY());
        atom->GetCoordinate()->SetZ(result.GetZ());
    }
//    std::cout << "\n";
    return;
}
