#include "../includes/glycosylationsite.h"

constexpr auto PI = 3.14159265358979323846;

//typedef std::vector<Overlap_record> OverlapRecordVector;


//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
GlycosylationSite::GlycosylationSite()
{
    SetGlycanName("");
    SetGlycanOverlap(0.0);
    SetProteinOverlap(0.0);
}

GlycosylationSite::GlycosylationSite(std::string glycan_name)
{
    SetGlycanName(glycan_name);
    SetGlycanOverlap(0.0);
    SetProteinOverlap(0.0);
}

GlycosylationSite::GlycosylationSite(std::string glycan_name, std::string residue_number)
{
    SetGlycanName(glycan_name);
    SetResidueNumber(residue_number);
    SetGlycanOverlap(0.0);
    SetProteinOverlap(0.0);
}

/*GlycosylationSite::GlycosylationSite(std::string glycan_name, Assembly glycan, Residue* residue)
{
    SetGlycanName(glycan_name);
    SetGlycan(glycan);
    SetResidue(residue);
    SetGlycanOverlap(0.0);
    SetProteinOverlap(0.0);
}
*/

GlycosylationSite::~GlycosylationSite()
{

}

//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////

Assembly* GlycosylationSite::GetGlycoprotein()
{
    return glycoprotein_;
}

std::string GlycosylationSite::GetGlycanName()
{
    return glycan_name_;
}

std::string GlycosylationSite::GetResidueNumber()
{
    return residue_number_;
}

Residue* GlycosylationSite::GetResidue()
{
    return residue_;
}

Assembly* GlycosylationSite::GetAttachedGlycan()
{
    return &glycan_;
}

double GlycosylationSite::GetOverlap()
{
    return (glycan_overlap_ + protein_overlap_);
}

double GlycosylationSite::GetWeightedOverlap(double glycan_weight, double protein_weight)
{
    return ( (glycan_overlap_ * glycan_weight) + (protein_overlap_ * protein_weight) );
}


double GlycosylationSite::GetGlycanOverlap()
{
    return glycan_overlap_;
}

double GlycosylationSite::GetProteinOverlap()
{
    return protein_overlap_;
}

AtomVector GlycosylationSite::GetSelfGlycanBeads()
{
    return self_glycan_beads_;
}

AtomVector GlycosylationSite::GetProteinBeads()
{
    return protein_beads_;
}

AtomVector GlycosylationSite::GetOtherGlycanBeads()
{
    return other_glycan_beads_;
}

//Residue_linkage GlycosylationSite::GetRotatableBonds()
//{
//    return residue_linkage_;
//}

ResidueLinkageVector GlycosylationSite::GetRotatableBonds()
{
    return all_residue_linkages_;
}

ResidueLinkageVector GlycosylationSite::GetFirstAnd1_6Linkages()
{
    ResidueLinkageVector returnedVector;
    returnedVector.push_back(all_residue_linkages_.front()); // Gives you the Protein-Glycan linkage

    for (auto &residue_linkage : all_residue_linkages_)
    {
        std::string residue2Name = residue_linkage.GetToThisResidue2()->GetName();
        std::string residue2AtomName = residue_linkage.GetToThisConnectionAtom2()->GetName();
        std::string residue1AtomName = residue_linkage.GetFromThisConnectionAtom1()->GetName();
        if( (residue2Name.compare(1,2,"MB")==0) && (residue2AtomName.compare("O6")==0) && (residue1AtomName.compare("C1")==0) )
        {
            std::cout << "Have found the inner core 1-6 linkage\n";
            returnedVector.push_back(residue_linkage);
        }
    }
    return returnedVector;
}

ResidueLinkageVector GlycosylationSite::GetFirstAnd2_XLinkages()
{
    ResidueLinkageVector returnedVector;
    returnedVector.push_back(all_residue_linkages_.front()); // Gives you the Protein-Glycan linkage

    for (auto &residue_linkage : all_residue_linkages_)
    {
        std::string residue1AtomName = residue_linkage.GetFromThisConnectionAtom1()->GetName();
        if( (residue1AtomName.compare("C2")==0) )
        {
            std::cout << "Have found the 2-X linkage\n";
            returnedVector.push_back(residue_linkage);
        }
    }
    return returnedVector;
}

std::vector<GlycosylationSite*> GlycosylationSite::GetOtherGlycosites()
{
    return other_glycosites_;
}


//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////

// Only need glycoprotein so can merge assemblies and set bonding for the connecting atom
// Bond by distance wouldn't work as may have overlaps after superimposition.
void GlycosylationSite::AttachGlycan(Assembly glycan, Assembly &glycoprotein)
{
    this->SetGlycan(glycan);
    this->SetGlycoprotein(&glycoprotein);
  //  std::cout << "\nGlycan set" << std::endl;
    this->Prepare_Glycans_For_Superimposition_To_Particular_Residue(residue_->GetName());
 //   std::cout << "Superimpose prep done" << std::endl;
    this->Superimpose_Glycan_To_Glycosite(residue_);
  //  std::cout << "Suerpimposed" << std::endl;
    this->Rename_Protein_Residue_To_GLYCAM_Nomenclature();
    glycoprotein.MergeAssembly(&glycan_); // Add glycan to glycoprotein assembly, allows SetDihedral later. May not be necessary anymore with new Rotatable Dihedral class.
 //   std::cout << "Merge done" << std::endl;
    //this->SetRotatableBonds(glycan_.GetResidues().at(0), residue_);
  // ResidueLinkageVector temp;
   // all_residue_linkages_.emplace_back(glycan_.GetResidues().at(0), residue_, this);
    all_residue_linkages_.emplace_back(glycan_.GetResidues().at(0), residue_);
    this->FigureOutResidueLinkagesInGlycan(glycan_.GetResidues().at(0), glycan_.GetResidues().at(0), &all_residue_linkages_);
  //  all_residue_linkages_ = temp;
}

/*
 * This function prepares the glycan molecule in the glycan_ assembly for superimpostion onto an amino acid in the protein
 * It does this by "growing" the atoms of the amino acid side chain (e.g. Asn, Thr or Ser) out from the glycan reducing terminal
 * Another function will use these additional atoms to superimpose the glycan onto residue_
 * This function assumes that the glycan_ assembly for this glycosylation site has already been set
*/
void GlycosylationSite::Prepare_Glycans_For_Superimposition_To_Particular_Residue(std::string amino_acid_name)
{
    //Dear future self, the order that you add the atoms to the residue matters for superimposition ie N, CA, CB , not CB, CA, N.
    // Want: assembly.FindResidueByTag("reducing-residue");
    Residue* reducing_Residue = glycan_.GetAllResiduesOfAssembly().at(1); // I assume I assumed something stupid here.
    // Want: residue.FindAtomByTag("anomeric-carbon"); The below is risky as it uses atoms names, i.e. would break for Sialic acid.
    Atom *atomC5 = reducing_Residue->GetAtom("C5");
    Atom *atomO5 = reducing_Residue->GetAtom("O5");
    Atom *atomC1 = reducing_Residue->GetAtom("C1");

    // Delete aglycon atoms from glycan.
    Residue * aglycon = glycan_.GetAllResiduesOfAssembly().at(0); // Oh jeez these assumptions are really building up.
    AtomVector aglycon_Atoms = aglycon->GetAtoms();
    for(AtomVector::iterator it = aglycon_Atoms.begin(); it != aglycon_Atoms.end(); ++it)
    {
       Atom* atom = *it;
       aglycon->RemoveAtom(atom); // Note only removes from residue. Atoms still exist and can be found through AtomNodes. They aren't written out in a PDB file though.
    }

    // Ok so going to set it so that the new "superimposition residue" is the old aglycon residue i.e. .at(0)
    // This avoids having to delete the algycon residue object from assembly and adding the super residue to assembly.
    // Deleting the residue is actually hard as the atoms still exist and are referenced from other places.
    Residue* superimposition_residue = aglycon; // "renaming" so the below reads better.
    superimposition_residue->SetName("SUP");
    superimposition_residue->SetId("SUP_?_1_?_?_1");

    // I put both the regular name and the O/N-linked glycam name here, as I'm not sure when it will be renamed.
    if ( (amino_acid_name.compare("ASN")==0) || (amino_acid_name.compare("NLN")==0) )
    {
        Atom *atomND2 = new Atom(superimposition_residue, "ND2", (GeometryTopology::get_cartesian_point_from_internal_coords(atomC5, atomO5, atomC1, 109.3, 180, 1.53)));
        Atom *atomCG = new Atom(superimposition_residue, "CG", (GeometryTopology::get_cartesian_point_from_internal_coords(atomO5, atomC1, atomND2, 109.3, 261, 1.325)));
        Atom *atomOD1 = new Atom(superimposition_residue, "OD1", (GeometryTopology::get_cartesian_point_from_internal_coords(atomC1, atomND2, atomCG, 126, 0, 1.22)));
        superimposition_residue->AddAtom(atomCG);
        superimposition_residue->AddAtom(atomOD1);
        superimposition_residue->AddAtom(atomND2);
        superimposition_atoms_ = superimposition_residue->GetAtoms();
    }
    else if ( (amino_acid_name.compare("THR")==0) || (amino_acid_name.compare("SER")==0) || (amino_acid_name.compare("OLT")==0) || (amino_acid_name.compare("OLS")==0) )
    {
        Atom *atomOG1 = new Atom(superimposition_residue, "OG", (GeometryTopology::get_cartesian_point_from_internal_coords(atomC5, atomO5, atomC1, 112, 68, 1.46)));
        Atom *atomCB = new Atom(superimposition_residue, "CB", (GeometryTopology::get_cartesian_point_from_internal_coords(atomO5, atomC1, atomOG1, 109.3, 75, 1.53)));
        Atom *atomCA = new Atom(superimposition_residue, "CA", (GeometryTopology::get_cartesian_point_from_internal_coords(atomC1, atomOG1, atomCB, 109.3, 125, 1.53)));
        superimposition_residue->AddAtom(atomCA);
        superimposition_residue->AddAtom(atomCB);
        superimposition_residue->AddAtom(atomOG1);
        superimposition_atoms_ = superimposition_residue->GetAtoms();
        std::cout << "Ding: " << amino_acid_name <<  std::endl;
        if ( (amino_acid_name.compare("THR")==0) || (amino_acid_name.compare("OLT")==0) )
        {
            atomOG1->SetName("OG1"); // It's OG in Ser.
            std::cout << "Ding" << std::endl;
        }
        atomOG1->GetCoordinate()->Print();
    }
    else if ( (amino_acid_name.compare("TYR")==0) || (amino_acid_name.compare("OLY")==0) )
    {
        Atom *atomOH = new Atom(superimposition_residue, "OH", (GeometryTopology::get_cartesian_point_from_internal_coords(atomC5, atomO5, atomC1, 112, 68, 1.46)));
        Atom *atomCZ = new Atom(superimposition_residue, "CZ", (GeometryTopology::get_cartesian_point_from_internal_coords(atomO5, atomC1, atomOH, 117, 60, 1.35)));
        Atom *atomCE1 = new Atom(superimposition_residue, "CE1", (GeometryTopology::get_cartesian_point_from_internal_coords(atomC1, atomOH, atomCZ, 120, 180, 1.37)));
        superimposition_residue->AddAtom(atomCE1);
        superimposition_residue->AddAtom(atomCZ);
        superimposition_residue->AddAtom(atomOH);
        superimposition_atoms_ = superimposition_residue->GetAtoms();
    }
    else
    {
        // OK I DON'T KNOW HOW TO HANDLE EXCEPTIONS. This will never happen though I promise.
        std::cout << "Problem in glycosylationsite::Prepare_Glycans_For_Superimposition_To_Particular_Residue(). Expect Segfault soon." << std::endl;
    }
    return;
}

void GlycosylationSite::Superimpose_Glycan_To_Glycosite(Residue *glycosite_residue)
{
    // Get the 3 target atoms from protein residue.
    //AtomVector protein_atoms = glycosite_residue->GetAtoms();
    AtomVector target_atoms;
    //AtomVector super_atoms = superimposition_atoms_.GetAllAtomsOfAssembly();
  //  Assembly* assemblyTarget = new Assembly();
   // Residue* residueTarget = new Residue();
    //Atom* last_protein_atom;
    // superimposition_atoms_ points to three atoms that were added to the glycan. Based on their names e.g. CG, ND2, we will superimpose them onto
    // the correspoinding "target" atoms in the protein residue (glycosite_residue).
    for (auto &superimposition_atom : superimposition_atoms_)
    {
        for(auto &protein_atom : glycosite_residue->GetAtoms())
        {
            if (protein_atom->GetName() == superimposition_atom->GetName())
            {
                target_atoms.push_back(protein_atom);
             //   std::cout << "Superimposition target acquired is " << protein_atom->GetName() << std::endl;
            }
        }
    }
//    for(AtomVector::iterator it1 = protein_atoms.begin(); it1 != protein_atoms.end(); ++it1)
//    {
//        Atom *protein_atom = (*it1);
//        for(AtomVector::iterator it2 = superimposition_atoms_.begin(); it2 != superimposition_atoms_.end(); ++it2)
//        {
//            Atom *super_atom = (*it2);
//            if (protein_atom->GetName() == super_atom->GetName())
//            {
//                target_atoms.push_back(protein_atom);
//                if (super_atom == superimposition_atoms_.back()) // Convention is that the last atom will be the one connecting to the glycan
//                {
//                    last_protein_atom = protein_atom; // To set the bonding correctly later.
//                }
//            }
//        }
//    }

    AtomVector glycan_atoms = glycan_.GetAllAtomsOfAssembly();

    gmml::Superimpose(superimposition_atoms_, target_atoms, glycan_atoms);

    Residue* reducing_Residue = glycan_.GetResidues().at(1); // I assume I assumed something stupid here.
  //  std::cout << "Reducing residue is " << reducing_Residue->GetName() << std::endl;
    AtomVector reducing_Atoms = reducing_Residue->GetAtoms();
    Atom* atomC1;
    for(AtomVector::iterator it = reducing_Atoms.begin(); it != reducing_Atoms.end(); it++)
    {
       Atom* atom = *it;
       if(atom->GetName().compare("C1")==0)
       {
           atomC1 = atom;
       }
    }
    //Connect the glycan and protein atoms to each other.
    Atom *protein_connection_atom = this->GetConnectingProteinAtom(glycosite_residue->GetName());

    protein_connection_atom->GetNode()->AddNodeNeighbor(atomC1);

    atomC1->GetNode()->AddNodeNeighbor(protein_connection_atom);
    //Delete the atoms used to superimpose the glycan onto the protein. Remove the residue.

    Residue *superimposition_residue = glycan_.GetAllResiduesOfAssembly().at(0);
    glycan_.RemoveResidue(superimposition_residue);
//    std::cout << "glycan_ now contains: ";
//    ResidueVector residues = glycan_.GetResidues();
//    for (ResidueVector::iterator it = residues.begin(); it != residues.end(); ++it)
//    {
//        std::cout << (*it)->GetName() << ", ";
//    }
//    std::cout << std::endl;
}

void GlycosylationSite::Rename_Protein_Residue_To_GLYCAM_Nomenclature()
{
    std::string amino_acid_name = this->GetResidue()->GetName();
    if (amino_acid_name.compare("ASN")==0) {this->GetResidue()->SetName("NLN");}
    if (amino_acid_name.compare("SER")==0) {this->GetResidue()->SetName("OLS");}
    if (amino_acid_name.compare("THR")==0) {this->GetResidue()->SetName("OLT");}
    if (amino_acid_name.compare("TYR")==0) {this->GetResidue()->SetName("OLY");}
}

void GlycosylationSite::Rename_Protein_Residue_From_GLYCAM_To_Standard()
{
    std::string amino_acid_name = this->GetResidue()->GetName();
    if (amino_acid_name.compare("NLN")==0) {this->GetResidue()->SetName("ASN");}
    if (amino_acid_name.compare("OLS")==0) {this->GetResidue()->SetName("SER");}
    if (amino_acid_name.compare("OLT")==0) {this->GetResidue()->SetName("THR");}
    if (amino_acid_name.compare("OLY")==0) {this->GetResidue()->SetName("TYR");}
}

double GlycosylationSite::CalculateOverlaps(OverlapType overlapType, MoleculeType moleculeType, bool recordOverlap, bool printOverlap)
{
    double overlap = 0.0;
    if (moleculeType == ALL) // glycan and protein overlap are calculated and stored sepatately. They get combined only when reporting total overlaps. Combined value is never stored.
    {
        double proteinOverlap = this->CalculateOverlaps(overlapType, PROTEIN, recordOverlap, printOverlap);
        double glycanOverlap = this->CalculateOverlaps(overlapType, GLYCAN, recordOverlap, printOverlap);
        overlap = proteinOverlap + glycanOverlap;
    }
    else
    {
        if(overlapType == ATOMIC)
        {
            overlap = this->CalculateAtomicOverlaps(moleculeType);
        }
        else if (overlapType == BEAD)
        {
            overlap = this->CalculateBeadOverlaps(moleculeType);
        }
        if (recordOverlap)
        {
            this->SetOverlap(moleculeType, overlap);
        }
    }
    if (printOverlap)
    {
        this->PrintOverlaps();
    }
    return overlap;
}

double GlycosylationSite::CalculateAtomicOverlaps(MoleculeType moleculeType)
{
    double overlap = 0.0;
    if(moleculeType == ALL)
    {
        overlap = ( this->CalculateAtomicOverlaps(PROTEIN) + this->CalculateAtomicOverlaps(GLYCAN) );
    }
    if(moleculeType == PROTEIN)
    {
        // WARNING. protein_atoms and glycan_atoms contains the bead atoms, but this function in GMML ignores them:
        AtomVector allProteinAtoms = glycoprotein_->GetAllAtomsOfAssemblyWithinProteinResidues();
        AtomVector attachmentResidueAtoms = residue_->GetAtoms();
        AtomVector proteinAtomsNotInAttachmentResidue = selection::GetAtomsin_a_Notin_b_AtomVectors(allProteinAtoms, attachmentResidueAtoms);
        overlap = gmml::CalculateAtomicOverlaps(proteinAtomsNotInAttachmentResidue, glycan_.GetAllAtomsOfAssembly());
      //  std::cout << residue_->GetId() << "-Protein: " << overlap << "\n";
    }
    if(moleculeType == GLYCAN)
    {
        for(auto &other_glycosite : other_glycosites_)
        {
        //    double current_overlap = gmml::CalculateAtomicOverlaps(other_glycosite->GetAttachedGlycan()->GetAllAtomsOfAssembly(), glycan_.GetAllAtomsOfAssembly());
        //    overlap += current_overlap;
        //    std::cout << residue_->GetId() << "-" << other_glycosite->GetResidue()->GetId() << ": " << current_overlap << "\n";
        overlap += gmml::CalculateAtomicOverlaps(other_glycosite->GetAttachedGlycan()->GetAllAtomsOfAssembly(), glycan_.GetAllAtomsOfAssembly());
        }
    }
    return overlap;
}

void GlycosylationSite::PrintOverlaps()
{
    std::cout << std::fixed; // Formating ouput
    std::cout << std::setprecision(2); // Formating ouput
    std::cout
        << std::setw(17) << this->GetResidue()->GetId() << " | "
        << std::setw(6)  << this->GetOverlap()     << " |  "
        << std::setw(6)  << this->GetProteinOverlap()   << " | "
        << std::setw(6)  << this->GetGlycanOverlap()    <<
    std::endl;
}

double GlycosylationSite::CalculateBeadOverlaps(MoleculeType moleculeType)
{
    double overlap = 0.0;
    if(moleculeType == ALL)
    {
        overlap = ( this->CalculateBeadOverlaps(PROTEIN) + this->CalculateBeadOverlaps(GLYCAN) );
    }
    if(moleculeType == PROTEIN)
    {
        overlap = this->CalculateBeadOverlaps(self_glycan_beads_, protein_beads_);

    }
    if(moleculeType == GLYCAN)
    {
        overlap = this->CalculateBeadOverlaps(self_glycan_beads_, other_glycan_beads_);
    }
    return overlap;
}


double GlycosylationSite::CalculateBeadOverlaps(AtomVector &atomsA, AtomVector &atomsB)
{
    double radius = 3.0; //Using same radius for all beads.
    double distance = 0.0, overlap = 0.0, current_overlap = 0.0;
 //   std::cout << "Glycosite" << this->GetResidueNumber() << "\n";
  //  std::cout << "About to check " << atomsA.size() << " atoms vs " << atomsB.size() << " atoms" << std::endl;
    for(AtomVector::iterator it1 = atomsA.begin(); it1 != atomsA.end(); ++it1)
    {
        Atom *atomA = *it1;
        for(AtomVector::iterator it2 = atomsB.begin(); it2 != atomsB.end(); ++it2)
        {
            Atom *atomB = *it2;
            if ( (atomA->GetCoordinates().at(0)->GetX() - atomB->GetCoordinates().at(0)->GetX()) < (radius * 2) ) // This is faster than calulating distance, and rules out tons of atom pairs.
            {
                distance = atomA->GetDistanceToAtom(atomB);
                if ( ( distance < (radius + radius) ) && ( distance > 0.0 ) ) //Close enough to overlap, but not the same atom
                {
                    current_overlap = gmml::CalculateAtomicOverlaps(atomA, atomB, radius, radius); // This calls the version with radius values
                    overlap += current_overlap;
                    //std::cout << atomA->GetResidue()->GetId() << " overlaping with " << atomB->GetResidue()->GetId() << ": " << current_overlap << "\n";
                }
            }
        }
    }
    //return (overlap );
    return (overlap / gmml::CARBON_SURFACE_AREA); //Normalise to area of a buried carbon
}


// Copied this from gmml Assembly. Oly updated to fix memory leaks
double GlycosylationSite::CalculateTorsionAngle(AtomVector atoms)
{

    double current_dihedral = 0.0;

    GeometryTopology::Coordinate *atom1_crd = atoms.at(0)->GetCoordinates().at(0);
    GeometryTopology::Coordinate *atom2_crd = atoms.at(1)->GetCoordinates().at(0);
    GeometryTopology::Coordinate *atom3_crd = atoms.at(2)->GetCoordinates().at(0);
    GeometryTopology::Coordinate *atom4_crd = atoms.at(3)->GetCoordinates().at(0);

    // Oliver updates to solve memory leaks
    GeometryTopology::Coordinate b1 = *atom2_crd; // deep copy
    GeometryTopology::Coordinate b2 = *atom3_crd;
    GeometryTopology::Coordinate b3 = *atom4_crd;
    GeometryTopology::Coordinate b4 = b2;
    b1.operator -(*atom1_crd);
    b2.operator -(*atom2_crd);
    b3.operator -(*atom3_crd);
    b4.operator *(-1);

    GeometryTopology::Coordinate b2xb3 = b2; // deep copy
    b2xb3.CrossProduct(b3);

    GeometryTopology::Coordinate b1_m_b2n = b1;
    b1_m_b2n.operator *(b2.length());

    GeometryTopology::Coordinate b1xb2 = b1; // deep copy
    b1xb2.CrossProduct(b2);

    current_dihedral = atan2(b1_m_b2n.DotProduct(b2xb3), b1xb2.DotProduct(b2xb3));
    return (current_dihedral * (180 / PI ) ); // Convert to DEGREES

}

void GlycosylationSite::WiggleFirstLinkage(int *output_pdb_id, double tolerance, int interval)
{
    //Residue_linkage *linkage = all_residue_linkages_.front();
    this->WiggleOneLinkage(all_residue_linkages_.front(), output_pdb_id, tolerance, interval);
    return;
}

void GlycosylationSite::Wiggle(int *output_pdb_id, double tolerance, int interval)
{ // I want to find the lowest overlap as close to each bonds default as possibe. So code is a bit more complicated.
    for(auto &linkage : all_residue_linkages_)
    {
        this->WiggleOneLinkage(linkage, output_pdb_id, tolerance, interval);
    }
    return;
}

void GlycosylationSite::write_pdb_file(Assembly *glycoprotein, int cycle, std::string summary_filename, double overlap)
{
    std::stringstream ss;
    ss << this->GetResidueNumber() << "_cyc_" << cycle << "ovrlp_" << overlap << ".pdb";
   // ss << cycle << "_cycle_" << ".pdb";

    PdbFileSpace::PdbFile *outputPdbFile = glycoprotein->BuildPdbFileStructureFromAssembly(-1,0);
    outputPdbFile->Write(ss.str());

}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void GlycosylationSite::SetGlycoprotein(Assembly* glycoprotein)
{
    glycoprotein_ = glycoprotein;
}

void GlycosylationSite::SetGlycanName(std::string glycan_name)
{
    glycan_name_ = glycan_name;
}

void GlycosylationSite::SetResidueNumber(std::string residue_number)
{
    residue_number_ = residue_number;
}

void GlycosylationSite::SetResidue(Residue* residue)
{
    residue_ = residue;
}

void GlycosylationSite::SetGlycan(Assembly glycan)
{
    glycan_ = glycan;
}

void GlycosylationSite::SetOverlap(MoleculeType moleculeType, double overlap)
{
    switch (moleculeType)
    {
    case PROTEIN:
        this->SetProteinOverlap(overlap);
        break;
    case GLYCAN:
        this->SetGlycanOverlap(overlap);
        break;
    default:
        std::cerr << "ERROR, No type specified in GlycosylationSite::SetOverlap. Required.\n";
        break;
    }
}

void GlycosylationSite::SetGlycanOverlap(double overlap)
{
    glycan_overlap_ = overlap;
}

void GlycosylationSite::SetProteinOverlap(double overlap)
{
    protein_overlap_ = overlap;
}

//void GlycosylationSite::SetChi1Value(double angle)
//{
//    Atom *atom1 = chi1_.at(0); // horrific, fix later.
//    Atom *atom2 = chi1_.at(1);
//    Atom *atom3 = chi1_.at(2);
//    Atom *atom4 = chi1_.at(3);

////    std::cout << std::fixed;
////    std::cout << std::setprecision(10);
////    std::cout << angle << std::endl;
////    std::cout << "Setting dihedral for " << atom1->GetName() << ", " << atom2->GetName() << ", " << atom3->GetName() << ", " << atom4->GetName() << "\n";
//    this->GetGlycoprotein()->SetDihedral(atom1, atom2, atom3, atom4, angle);
////    std::cout << this->GetChi1Value() << std::endl;

//}


void GlycosylationSite::SetSelfGlycanBeads(AtomVector *beads)
{
    self_glycan_beads_ = *beads;
}

void GlycosylationSite::SetProteinBeads(AtomVector *beads)
{
    protein_beads_ = *beads; // This creates a copy
    //Remove beads from attachment point residue. Don't want to count overlaps between glycan and residue it is attached to.
    for(AtomVector::iterator it1 = protein_beads_.begin(); it1 != protein_beads_.end(); /* Not incrementing here as erasing increments*/)
    {
        Atom *atom = *it1;
        if (atom->GetResidue() == this->GetResidue()) // this->GetResidue returns the glycosites protein residue.
        {
            protein_beads_.erase(std::remove(protein_beads_.begin(), protein_beads_.end(), *it1), protein_beads_.end());
        }
        else
        {
            ++it1; // "erase" increments it1, so if no erase happens this does an increment instead.
        }
    }
}

void GlycosylationSite::SetOtherGlycanBeads(AtomVector *beads)
{
    other_glycan_beads_ = *beads;
}

void GlycosylationSite::SetDefaultDihedralAnglesUsingMetadata()
{
    for(auto &linkage : all_residue_linkages_)
    {
        linkage.SetDefaultShapeUsingMetadata();
    }
    //residue_linkage_.SetDefaultDihedralAnglesUsingMetadata();
}

void GlycosylationSite::SetRandomDihedralAnglesUsingMetadata()
{
    for(auto &linkage : all_residue_linkages_)
    {
        linkage.SetRandomShapeUsingMetadata();
    }
    //residue_linkage_.SetRandomDihedralAnglesUsingMetadata();
}

void GlycosylationSite::SetRandomDihedralAnglesUsingMetadataForNthLinkage(int linkage_number)
{
    if(linkage_number < all_residue_linkages_.size())
    {
        all_residue_linkages_.at(linkage_number).SetRandomShapeUsingMetadata();
    }
}

void GlycosylationSite::ResetDihedralAngles()
{
    for(auto &linkage : all_residue_linkages_)
    {
        linkage.SetShapeToPrevious();
    }
    //residue_linkage_.SetPreviousDihedralAngles();
}

void GlycosylationSite::UpdateAtomsThatMoveInLinkages()
{
    for(auto &linkage : all_residue_linkages_)
    {
        linkage.DetermineAtomsThatMove();
    }
    //residue_linkage_.DetermineAtomsThatMove();
}

void GlycosylationSite::StashCoordinates()
{
    std::cout << "Stashing coordinates for glycosite " << this->GetResidueNumber() << std::endl;
    AtomVector atoms = this->GetAttachedGlycan()->GetAllAtomsOfAssembly();
    AtomVector sidechain_atoms = this->GetResidue()->GetAtoms();
    atoms.insert(atoms.end(), sidechain_atoms.begin(), sidechain_atoms.end() );
    for(auto &atom : atoms)
    {
        atom->AddCoordinate(new GeometryTopology::Coordinate(atom->GetCoordinate())); // push back the currect coordinates onto the end of the coordinates vector.
    }
}

void GlycosylationSite::SetStashedCoordinates() // When a lower overlap is found, the coords are pushed back. Last one pushed back is the lowest overlap.
{
    std::cout << "Setting stashed coordinates for site " << this->GetResidueNumber() << "\n";
    AtomVector atoms = this->GetAttachedGlycan()->GetAllAtomsOfAssembly();
    AtomVector sidechain_atoms = this->GetResidue()->GetAtoms();
    atoms.insert(atoms.end(), sidechain_atoms.begin(), sidechain_atoms.end() );
    for(auto &atom : atoms)
    {
        GeometryTopology::Coordinate *first_coords = atom->GetCoordinate();
        GeometryTopology::Coordinate *last_coords = atom->GetCoordinates().back();
        *first_coords = *last_coords;
    }
    //std::cout << "Finished Setting stashed coordinates" << std::endl;
}

// Very interesting that I have to specify std::vector<GlycosylationSite>, but can use the typedef in the function.
// A HUGE clue. MASSIVE.
std::vector<GlycosylationSite> GlycosylationSite::GetXClosestSitesWithinOverlapDistanceY(GlycosylationSiteVector &glycosites, int maxNumberOfSitesToConsider)
{
    typedef std::pair<double, GlycosylationSite> siteDistancePair;
    typedef std::vector<siteDistancePair> siteDistancePairVector;
    siteDistancePairVector sitesWithinOverlapRange;
    GlycosylationSiteVector sitesToReturn;
    double maxLengthOfThisGlycan = (selection::GetMaxDistanceBetweenAtoms(this->GetAttachedGlycan()->GetAllAtomsOfAssembly()) + 5);
    Atom *this_cb_atom = this->GetResidue()->GetAtom("CB");
    for(auto &glycosite : glycosites)
    {
        std::cout << "    distance to " << glycosite.GetResidue()->GetId() << ": ";
        if (this->GetResidue()->GetId().compare(glycosite.GetResidue()->GetId())==0) // if not this site. Oly you should overload the = operator?
        {
            std::cout << "Skipping\n";
            continue; // I think this skips to next item in for loop.
        }
        double maxLengthOfThatGlycan = (selection::GetMaxDistanceBetweenAtoms(glycosite.GetAttachedGlycan()->GetAllAtomsOfAssembly()) + 5);
        Atom *that_cb_atom = glycosite.GetResidue()->GetAtom("CB");
        double distanceBetweenCBAtoms = this_cb_atom->GetDistanceToAtom(that_cb_atom);
        std::cout << distanceBetweenCBAtoms;
        if ( ( distanceBetweenCBAtoms ) <= (maxLengthOfThisGlycan + maxLengthOfThatGlycan) )
        {
            std::cout << " withing range! (" << (maxLengthOfThisGlycan + maxLengthOfThatGlycan) << ")" << std::endl;
            // push back onto vector site and distance pairs
            sitesWithinOverlapRange.emplace_back(distanceBetweenCBAtoms, glycosite);
        }
        else
            std::cout << "\n";
    }
    // Sort list by distance. Remember that you might have more or less sites than the maxNumberOfSitesToConsider
    std::sort(sitesWithinOverlapRange.begin(), sitesWithinOverlapRange.end());
    // Return up to maxNumberOfSitesToConsider
    for(auto &pairInstance : sitesWithinOverlapRange)
    {
        std::cout << "Sorted distance: " << pairInstance.first << "\n";
        if (sitesToReturn.size() < maxNumberOfSitesToConsider)
        {
            sitesToReturn.push_back(pairInstance.second);
        }
    }
    return sitesToReturn;
}

void GlycosylationSite::SetOtherGlycosites(GlycosylationSiteVector &glycosites)
{
    for(auto &glycosite : glycosites)
    {
        if(this->GetResidue()->GetId().compare(glycosite.GetResidue()->GetId())!=0)
            other_glycosites_.push_back(&glycosite);
    }
    return;
}

//////////////////////////////////////////////////////////
//                       DISPLAY FUNCTION               //
//////////////////////////////////////////////////////////

void GlycosylationSite::Print(std::string type)
{
    if (type.compare("All")==0)
    {
        std::cout << "Residue ID: " << this->GetResidue()->GetId()
                  << ", overlap: " << this->GetOverlap()
                  << ", Gly--Pro " << this->GetGlycanOverlap() << "--" << this->GetProteinOverlap()
                  << std::endl;
        //std::cout << "Dihedrals: \n";
        //this->GetRotatableBonds().Print();
        //rotatable_bonds_.Print();
        //std::cout << "\n";
    }
}


//////////////////////////////////////////////////////////
//                   PRIVATE FUNCTIONS                 //
//////////////////////////////////////////////////////////

//void GlycosylationSite::SetRotatableBonds(Residue *residue1, Residue *residue2)
//{

//    Residue_linkage rotatable_bonds(residue1, residue2);
//    residue_linkage_ = rotatable_bonds;
//    // Copy is ok for now. I set up Residue_linkage to be consructed, but would need Glycosidic linkage to be constructed all at once too. Need to look into that.
//    // I.e. Construct everything all at once via "constructors" <-- lolyes.

//}

// Below doesn't work as ResidueNode is not set.
//void GlycosylationSite::FigureOutResidueLinkagesInGlycanOld(Residue *residue1, Residue *residue2, ResidueLinkageVector *residue_linkages)
//{
//    std::cout << "Brutus?" << std::endl;
//    for(auto &node : residue2->GetNode()->GetResidueNodeNeighbors())
//    {
//        Residue *next_residue = node->GetResidue();
//        if ( next_residue->GetId().compare(residue1->GetId())!=0) // if not the previous residue
//        {
//            std::cout << "Et tu," << std::endl;
//            Residue_linkage steve(residue2, next_residue);
//            std::cout << "Buttalis?" << std::endl;
//            //residue_linkages->push_back(steve);
//            FigureOutResidueLinkagesInGlycanOld(residue2, next_residue, residue_linkages);
//        }
//    }
//    return;
//}

// These functions do not fit here. Also, as ResidueNodes are not set, they are overly complex (nested recursive functions anyone?).
// They traverse the residues. I want to ignore the algycon for now.
void GlycosylationSite::FigureOutResidueLinkagesInGlycan(Residue *from_this_residue1, Residue *to_this_residue2, ResidueLinkageVector *residue_linkages)
{
    Atom *start_atom = to_this_residue2->GetAtoms().at(0); // Need to start somewhere.
    ResidueVector neighbors;
    RecursivelyGetAllNeighboringResidues(start_atom, &neighbors);
    ResidueVector glycan_residues = glycan_.GetResidues();
    for(auto &neighbor : neighbors)
    {
        if( (neighbor->GetIndex() != from_this_residue1->GetIndex()) && (std::find(glycan_residues.begin(), glycan_residues.end(), neighbor) != glycan_residues.end())  )
        {
            residue_linkages->emplace_back(neighbor, to_this_residue2);
        }
    }
    for(auto &neighbor : neighbors)
    {
        if( (neighbor->GetIndex() != from_this_residue1->GetIndex()) && (std::find(glycan_residues.begin(), glycan_residues.end(), neighbor) != glycan_residues.end())  )
        {
            this->FigureOutResidueLinkagesInGlycan(to_this_residue2, neighbor, residue_linkages);
        }
    }
    return;
}

void GlycosylationSite::RecursivelyGetAllNeighboringResidues(Atom* current_atom, ResidueVector* neighbors)
{
    current_atom->SetDescription("VisitedByRecursivelyGetAllNeighboringResidues");
    for(auto &neighboring_atom : current_atom->GetNode()->GetNodeNeighbors())
    {
        unsigned long long neighbor_index = neighboring_atom->GetResidue()->GetIndex();
        if(neighbor_index != current_atom->GetResidue()->GetIndex())
        {
            neighbors->push_back(neighboring_atom->GetResidue());
//            std::cout << "Foreign neighbor of " << current_atom->GetId() << " is " << neighboring_atom->GetId() << "\n";
//            std::cout << "size is now " << neighbors->size() << "\n";
        }
        else if (neighboring_atom->GetDescription().compare("VisitedByRecursivelyGetAllNeighboringResidues")!=0)
        { // If in current residue, and not visited already:
            this->RecursivelyGetAllNeighboringResidues(neighboring_atom, neighbors);
        }
    }
}

Atom* GlycosylationSite::GetConnectingProteinAtom(std::string residue_name)
{
    if(residue_name == "NLN" || residue_name == "ASN")
    {
        return residue_->GetAtom("ND2");
    }
    else if(residue_name == "OLT" || residue_name == "THR")
    {
        return residue_->GetAtom("OG1");
    }
    else if(residue_name == "OLS" || residue_name == "SER")
    {
        return residue_->GetAtom("OG");
    }
    else if(residue_name == "OLY" || residue_name == "TYR")
    {
        return residue_->GetAtom("OH");
    }
    else
    {
        std::cout << "In GlycosylationSite::GetConectinProteinAtom you passed in a residue (" << residue_name << ") that isn't ON THE LIST. Exiting early: " << std::endl;
        exit(1);
    }
}


// OG re-reading. I'm pretty sure this belongs in Residue_linkage. linkage.wiggle(output_pdb_id, tolerance, interval
// OG re-re-reading. It's because you need to calculate overlaps you nonce.
void GlycosylationSite::WiggleOneLinkage(Residue_linkage &linkage, int *output_pdb_id, double tolerance, int interval)
{
    double current_overlap = this->CalculateBeadOverlaps();
    double lowest_overlap = current_overlap;
    // Reverse as convention is Glc1-4Gal and I want to wiggle in opposite direction i.e. from first rotatable bond in Asn outwards
    RotatableDihedralVector reversed_rotatable_bond_vector = linkage.GetRotatableDihedrals();
    std::reverse(reversed_rotatable_bond_vector.begin(), reversed_rotatable_bond_vector.end());
    for(auto &rotatable_dihedral : reversed_rotatable_bond_vector)
    {
        double best_dihedral_angle = rotatable_dihedral.CalculateDihedralAngle();
       // std::cout << "Starting new linkage with best angle as " << best_dihedral_angle << "\n";
        gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector metadata_entries = rotatable_dihedral.GetMetadata();
        for(auto &metadata : metadata_entries)
        {
            double lower_bound = (metadata.default_angle_value_ - metadata.lower_deviation_);
            double upper_bound = (metadata.default_angle_value_ + metadata.upper_deviation_);
            double current_dihedral = lower_bound;
            while(current_dihedral <= upper_bound )
            {
                rotatable_dihedral.SetDihedralAngle(current_dihedral);
                //GlycosylationSite::write_pdb_file(this->GetGlycoprotein(), *output_pdb_id, "wiggle", lowest_overlap);
                //++(*output_pdb_id);
                current_overlap = this->CalculateBeadOverlaps();
              //  std::cout << this->GetResidueNumber() << ": current dihedral : overlap " << current_dihedral << " : " << current_overlap << ". Best dihedral : overlap: " << best_dihedral_angle << " : "<< lowest_overlap << "\n";
                if (lowest_overlap >= (current_overlap + 0.01)) // 0.01 otherwise rounding errors
                {
                  //  std::cout << "Setting id " << *output_pdb_id << " index: " << metadata.index_ << ": ";
                    //rotatable_dihedral.Print();
                    lowest_overlap = current_overlap;
                    //                        glycoprotein_builder::write_pdb_file(this->GetGlycoprotein(), *output_pdb_id, "wiggle", lowest_overlap);
                    //                        ++(*output_pdb_id);
                    best_dihedral_angle = current_dihedral;
                //    std::cout << "Best angle is now " << best_dihedral_angle << "\n";
                //    GlycosylationSite::write_pdb_file(this->GetGlycoprotein(), *output_pdb_id, "wiggle", lowest_overlap);
                }
                // Perfer angles closer to default.
                else if ( (lowest_overlap == current_overlap) &&
                          (abs(metadata.default_angle_value_ - best_dihedral_angle ) > abs(metadata.default_angle_value_ - current_dihedral)) )
                {
                    best_dihedral_angle = current_dihedral;
                }
                current_dihedral += interval; // increment
            }
        }
        //std::cout << "Setting best angle as " << best_dihedral_angle << "\n";
        rotatable_dihedral.SetDihedralAngle(best_dihedral_angle);
        // std::cout << "\n";
        if(lowest_overlap <= tolerance) return;
    }
    return; // Note possibility of earlier return above
}




