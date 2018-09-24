#include "../includes/rotatable_dihedral.h"
#include <random>

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
Rotatable_dihedral::Rotatable_dihedral()
{
    std::cout << "This is a bad idea, use one of the other constructors\n";
}

// Of the next two forms, I'll probably use only one and delete the other later
Rotatable_dihedral::Rotatable_dihedral(Atom *atom1, Atom *atom2, Atom *atom3, Atom *atom4)
{
    atom1_ = atom1;
    atom2_ = atom2;
    atom3_ = atom3;
    atom4_ = atom4;
    AtomVector atoms_that_move;
    atoms_that_move.push_back(atom2_);
    atom3_->FindConnectedAtoms(atoms_that_move);
    atoms_that_move_ = atoms_that_move;
}

Rotatable_dihedral::Rotatable_dihedral(AtomVector atoms)
{
    this->SetAtoms(atoms);
    AtomVector atoms_that_move;
    atoms_that_move.push_back(atom2_);
    atom3_->FindConnectedAtoms(atoms_that_move);
    atoms_that_move_ = atoms_that_move;
}

// Not sure if I'll want to trust connection atoms or not, maybe pass it in like this:
Rotatable_dihedral::Rotatable_dihedral(AtomVector atoms, AtomVector atoms_that_move)
{
    this->SetAtoms(atoms);
    atoms_that_move_ = atoms_that_move;
}

//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////

double Rotatable_dihedral::CalculateDihedralAngle()
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

    double current_dihedral_angle = atan2(b1_m_b2n.DotProduct(b2xb3), b1xb2.DotProduct(b2xb3));

    return current_dihedral_angle;
}

AtomVector Rotatable_dihedral::GetAtoms()
{
    AtomVector atoms = {atom1_, atom2_, atom3_, atom4_};
    return atoms;
}
AtomVector Rotatable_dihedral::GetAtomsThatMove()
{
    return atoms_that_move_;
}

double Rotatable_dihedral::GetPreviousDihedralAngle()
{
    return previous_dihedral_angle_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void Rotatable_dihedral::SetAtoms(AtomVector atoms)
{
    atom1_ = atoms.at(0);
    atom2_ = atoms.at(1);
    atom3_ = atoms.at(2);
    atom4_ = atoms.at(3);
}

void Rotatable_dihedral::SetAtomsThatMove(AtomVector atoms)
{
    atoms_that_move_ = atoms;
}

double Rotatable_dihedral::RandomizeDihedralAngle()
{
    return Rotatable_dihedral::RandomizeDihedralAngleWithinRange(0.0, 360.0);
    //return (rand() % 360) + 1 - 180; // Can get same one everytime for testing
}

double Rotatable_dihedral::RandomizeDihedralAngleWithinRange(double min, double max)
{
    std::random_device rd1; // obtain a random number from hardware
    std::mt19937 eng1(rd1()); // seed the generator
    std::uniform_real_distribution<> angle_distribution(min, max); // define the range

    double random_angle = angle_distribution(eng1);

    /*******************************************/
    /*               IMPORTANT                 */
    /*******************************************/

    this->SetDihedralAngle(random_angle); // THIS IS IMPORTANT!!! THIS SHOULD BE SEPARATED?!?!

    /*******************************************/
    /*               IMPORTANT                 */
    /*******************************************/

    return random_angle;
    //return rand() % (max + 1 - min) + min; // Can get same one everytime for testing
}

double Rotatable_dihedral::RandomizeDihedralAngleWithinRanges(std::vector<std::pair<double,double>> ranges)
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

    // create an angle within the selected range
    return Rotatable_dihedral::RandomizeDihedralAngleWithinRange(ranges.at(range_selection).first, ranges.at(range_selection).second);
}

void Rotatable_dihedral::RecordPreviousDihedralAngle(double dihedral_angle)
{
    previous_dihedral_angle_=dihedral_angle;
}

void Rotatable_dihedral::ResetDihedralAngle()
{
    this->SetDihedralAngle(this->GetPreviousDihedralAngle());
}

void Rotatable_dihedral::SetDihedralAngle(double dihedral_angle)
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
    double** dihedral_angle_matrix = gmml::GenerateRotationMatrix(&b4, a2, current_dihedral - gmml::ConvertDegree2Radian(dihedral_angle));
    this->RecordPreviousDihedralAngle(current_dihedral);

//    AtomVector atomsToRotate = AtomVector();
//    atomsToRotate.push_back(atom2);
//    atom3->FindConnectedAtoms(atomsToRotate);
 //   std::cout << "Moving: ";
    // Yo you should add something here that checks if atoms_that_move_ is set. Yeah you.
    for(AtomVector::iterator it = atoms_that_move_.begin(); it != atoms_that_move_.end(); it++)
    {
        Atom *atom = *it;
        GeometryTopology::Coordinate* atom_coordinate = atom->GetCoordinate();
        GeometryTopology::Coordinate result;
        result.SetX(dihedral_angle_matrix[0][0] * atom_coordinate->GetX() + dihedral_angle_matrix[0][1] * atom_coordinate->GetY() +
                dihedral_angle_matrix[0][2] * atom_coordinate->GetZ() + dihedral_angle_matrix[0][3]);
        result.SetY(dihedral_angle_matrix[1][0] * atom_coordinate->GetX() + dihedral_angle_matrix[1][1] * atom_coordinate->GetY() +
                dihedral_angle_matrix[1][2] * atom_coordinate->GetZ() + dihedral_angle_matrix[1][3]);
        result.SetZ(dihedral_angle_matrix[2][0] * atom_coordinate->GetX() + dihedral_angle_matrix[2][1] * atom_coordinate->GetY() +
                dihedral_angle_matrix[2][2] * atom_coordinate->GetZ() + dihedral_angle_matrix[2][3]);

        atom->GetCoordinate()->SetX(result.GetX());
        atom->GetCoordinate()->SetY(result.GetY());
        atom->GetCoordinate()->SetZ(result.GetZ());
    }
//    std::cout << "\n";
    return;
}

//////////////////////////////////////////////////////////
//                       DISPLAY FUNCTION               //
//////////////////////////////////////////////////////////

void Rotatable_dihedral::Print()
{
    std::cout << atom1_->GetName() << ", " << atom2_->GetName() << ", " << atom3_->GetName() << ", " << atom4_->GetName() << ": " << this->CalculateDihedralAngle() << ".";
}
