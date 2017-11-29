#include "glycosylationsite.h"

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
GlycosylationSite::GlycosylationSite()
{
    glycan_name_ = "";
    residue_ = NULL;
}

GlycosylationSite::GlycosylationSite(std::string glycan_name)
{
    glycan_name_ = glycan_name;
}

GlycosylationSite::GlycosylationSite(std::string glycan_name, Residue* residue, Assembly glycan)
{
    glycan_name_ = glycan_name;
    residue_ = residue;
    glycan_ = glycan;
}

GlycosylationSite::~GlycosylationSite()
{

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

Assembly* GlycosylationSite::GetAttachedGlycan()
{
    return &glycan_;
}

//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////

// Only need glycoprotein so can merge assemblies and set bonding for the connecting atom
// Bond by distance wouldn't work as may have overlaps after superimposition.
void GlycosylationSite::AttachGlycan(Assembly glycan, Assembly *glycoprotein)
{
    this->SetGlycan(glycan);
    std::cout << "Here" << std::endl;
    glycoprotein->MergeAssembly(&glycan_); // Add glycan to glycoprotein assembly, allows SetDihedral later.
    this->Prepare_Glycans_For_Superimposition_To_Particular_Residue(residue_->GetName());
    this->Superimpose_Glycan_To_Glycosite(residue_);
}

void GlycosylationSite::Prepare_Glycans_For_Superimposition_To_Particular_Residue(std::string amino_acid_name)
{
    //Dear future self, the order that you add the atoms to the residue matters for superimposition ie N, CA, CB , not CB, CA, N.

    Residue* reducing_Residue = glycan_.GetAllResiduesOfAssembly().at(1); // I assume I assumed something stupid here.
    //std::cout << "Reducing residue is " << reducing_Residue->GetName() << std::endl;
    AtomVector reducing_Atoms = reducing_Residue->GetAtoms();

    Atom* atomC5;
    Atom* atomO5;
    Atom* atomC1;

    // Want: residue.FindAtom(string Name);
    for(AtomVector::iterator it = reducing_Atoms.begin(); it != reducing_Atoms.end(); it++)
    {
       Atom* atom = *it;
       if(atom->GetName().compare("C5")==0)
           atomC5 = atom;
       if(atom->GetName().compare("O5")==0)
           atomO5 = atom;
       if(atom->GetName().compare("C1")==0)
           atomC1 = atom;
    }
    // Delete aglycon atoms from glycan.
    Residue * aglycon = glycan_.GetAllResiduesOfAssembly().at(0);
    AtomVector aglycon_Atoms = aglycon->GetAtoms();
    for(AtomVector::iterator it = aglycon_Atoms.begin(); it != aglycon_Atoms.end(); it++)
    {
       Atom* atom = *it;
//       std::cout << "Removing " << atom->GetName() << std::endl;
       aglycon->RemoveAtom(atom);
    }

   // Assembly* assembly = new Assembly();
    Residue* super_residue = new Residue();
    Residue* alt_residue = new Residue();
    super_residue->SetAssembly(&superimposition_atoms_);
    superimposition_atoms_.AddResidue(super_residue);

   // PdbFileSpace::PdbFile *outputPdbFile = new PdbFileSpace::PdbFile();

    if (amino_acid_name.compare("ASN")==0)
    {
        super_residue->SetName("NLN");
        super_residue->SetId("NLN_?_1_?_?_1");

        Atom *atomND2 = new Atom(super_residue, "ND2", (gmml::get_cartesian_point_from_internal_coords(atomC5, atomO5, atomC1, 109.3, 180, 1.53)));
        Atom *atomCG = new Atom(super_residue, "CG", (gmml::get_cartesian_point_from_internal_coords(atomO5, atomC1, atomND2, 109.3, 261, 1.325)));
        Atom *atomOD1 = new Atom(super_residue, "OD1", (gmml::get_cartesian_point_from_internal_coords(atomC1, atomND2, atomCG, 126, 0, 1.22)));
        Atom *atomCB = new Atom(super_residue, "CB", (gmml::get_cartesian_point_from_internal_coords(atomC1, atomND2, atomCG, 114, 177.3, 1.53)));
        Atom *atomCA = new Atom(super_residue, "CA", (gmml::get_cartesian_point_from_internal_coords(atomND2, atomCG, atomCB, 111, 177.6, 1.53)));
        Atom *atomN = new Atom(super_residue, "N", (gmml::get_cartesian_point_from_internal_coords(atomCG, atomCB, atomCA, 111, 191.6, 1.453)));

        super_residue->AddAtom(atomCG);
        super_residue->AddAtom(atomOD1);
        super_residue->AddAtom(atomND2);
        superimposition_atoms_.BuildStructureByDistance();

        alt_residue->SetName("NLN");
        alt_residue->SetId("NLN_?_1_?_?_1");
        alt_residue->SetAssembly(&alternate_sidechain_);
        alternate_sidechain_.AddResidue(alt_residue);
        alt_residue->AddAtom(atomN);
        alt_residue->AddAtom(atomCA);
        alt_residue->AddAtom(atomCB);
        alt_residue->AddAtom(atomCG);
        alt_residue->AddAtom(atomOD1);
        alt_residue->AddAtom(atomND2);
        alternate_sidechain_.BuildStructureByDistance();

       // outputPdbFile = alternate_sidechain_.BuildPdbFileStructureFromAssembly();
       // outputPdbFile->Write("/home/oliver/Programs/Cplusplus/GlycoproteinBuilder/glycoproteinBuilder/outputs/NLN_AlignedToGlycan.pdb");
    }
    else if (amino_acid_name.compare("THR")==0 || amino_acid_name.compare("SER")==0)
    {
        super_residue->SetName("OLS");
        super_residue->SetId("OLS_?_1_?_?_1");

        Atom *atomOG1 = new Atom(super_residue, "OG", (gmml::get_cartesian_point_from_internal_coords(atomC5, atomO5, atomC1, 112, 68, 1.46)));
        Atom *atomCB = new Atom(super_residue, "CB", (gmml::get_cartesian_point_from_internal_coords(atomO5, atomC1, atomOG1, 109.3, 75, 1.53)));
        Atom *atomCA = new Atom(super_residue, "CA", (gmml::get_cartesian_point_from_internal_coords(atomC1, atomOG1, atomCB, 109.3, 125, 1.53)));
        Atom *atomN = new Atom(super_residue, "N", (gmml::get_cartesian_point_from_internal_coords(atomOG1, atomCB, atomCA, 109.3, 180, 1.53)));
        Atom *atomCG2 = new Atom(super_residue, "CG2", (gmml::get_cartesian_point_from_internal_coords(atomC1, atomOG1, atomCB, 109.3, -60, 1.53)));

        super_residue->AddAtom(atomCA);
        super_residue->AddAtom(atomCB);
        super_residue->AddAtom(atomOG1);
        superimposition_atoms_.BuildStructureByDistance();

        alt_residue->SetName("OLS");
        alt_residue->SetId("OLS_?_1_?_?_1");
        alt_residue->SetAssembly(&alternate_sidechain_);
        alternate_sidechain_.AddResidue(alt_residue);
        alt_residue->AddAtom(atomN);
        alt_residue->AddAtom(atomCA);
        alt_residue->AddAtom(atomCB);
        alt_residue->AddAtom(atomOG1);

        if (amino_acid_name.compare("SER")==0)
        {
            alternate_sidechain_.BuildStructureByDistance();
        }
        else
        {
            alt_residue->SetName("OLT"); //Thr = Ser + CG2
            alt_residue->SetId("OLT_?_1_?_?_1");
            alt_residue->AddAtom(atomCG2);
            atomOG1->SetName("OG1"); // It's OG in Ser.
            alternate_sidechain_.BuildStructureByDistance();
        }
    }

    else if (amino_acid_name.compare("TYR")==0)
    {
        super_residue->SetName("OLY");
        super_residue->SetId("OLY_?_1_?_?_1");

        Atom *atomOH = new Atom(super_residue, "OH", (gmml::get_cartesian_point_from_internal_coords(atomC5, atomO5, atomC1, 112, 68, 1.46)));
        Atom *atomCZ = new Atom(super_residue, "CZ", (gmml::get_cartesian_point_from_internal_coords(atomO5, atomC1, atomOH, 117, 60, 1.35)));
        Atom *atomCE1 = new Atom(super_residue, "CE1", (gmml::get_cartesian_point_from_internal_coords(atomC1, atomOH, atomCZ, 120, 180, 1.37)));
        Atom *atomCD1 = new Atom(super_residue, "CD1", (gmml::get_cartesian_point_from_internal_coords(atomOH, atomCZ, atomCE1, 120, 180, 1.37)));
        Atom *atomCE2 = new Atom(super_residue, "CE2", (gmml::get_cartesian_point_from_internal_coords(atomC1, atomOH, atomCZ, 120, 0, 1.37)));
        Atom *atomCD2 = new Atom(super_residue, "CD2", (gmml::get_cartesian_point_from_internal_coords(atomOH, atomCZ, atomCE2, 120, 180, 1.37)));
        Atom *atomCG = new Atom(super_residue, "CG", (gmml::get_cartesian_point_from_internal_coords(atomCZ, atomCE2, atomCD2, 120, 0, 1.37)));
        Atom *atomCB = new Atom(super_residue, "CB", (gmml::get_cartesian_point_from_internal_coords(atomCE2, atomCD2, atomCG, 122, 180, 1.51)));
        Atom *atomCA = new Atom(super_residue, "CA", (gmml::get_cartesian_point_from_internal_coords(atomCD2, atomCG, atomCB, 111, -107, 1.55)));
        Atom *atomN = new Atom(super_residue, "N", (gmml::get_cartesian_point_from_internal_coords(atomCG, atomCB, atomCA, 114, -170, 1.44)));

        super_residue->AddAtom(atomCE1);
        super_residue->AddAtom(atomCZ);
        super_residue->AddAtom(atomOH);
        superimposition_atoms_.BuildStructureByDistance();

        alt_residue->SetName("OLY");
        alt_residue->SetId("OLY_?_1_?_?_1");
        alt_residue->SetAssembly(&alternate_sidechain_);
        alternate_sidechain_.AddResidue(alt_residue);
        alt_residue->AddAtom(atomN);
        alt_residue->AddAtom(atomCA);
        alt_residue->AddAtom(atomCB);
        alt_residue->AddAtom(atomOH);
        alt_residue->AddAtom(atomCZ);
        alt_residue->AddAtom(atomCE1);
        alt_residue->AddAtom(atomCD1);
        alt_residue->AddAtom(atomCE2);
        alt_residue->AddAtom(atomCD2);
        alt_residue->AddAtom(atomCG);
        alternate_sidechain_.BuildStructureByDistance();
    }
    return;
}

void GlycosylationSite::Superimpose_Glycan_To_Glycosite(Residue *glycosite_residue)
{
    // Get the 3 target atoms from protein residue.
    std::cout << "Superimposing to " << glycosite_residue->GetName() << std::endl;
    AtomVector target_atoms = glycosite_residue->GetAtoms();
    AtomVector super_atoms = superimposition_atoms_.GetAllAtomsOfAssembly();
    Assembly* assemblyTarget = new Assembly();
    Residue* residueTarget = new Residue();
    Atom* last_protein_atom;
    for(AtomVector::iterator it1 = target_atoms.begin(); it1 != target_atoms.end(); ++it1)
    {
        Atom *target_atom = (*it1);
        for(AtomVector::iterator it2 = super_atoms.begin(); it2 != super_atoms.end(); ++it2)
        {
            Atom *super_atom = (*it2);
            if (target_atom->GetName() == super_atom->GetName())
            {
                residueTarget->AddAtom(target_atom);
                std::cout << "Here1" << std::endl;
                if (super_atom == super_atoms.back()) // This is just to set the bonding correctly.
                {
                    std::cout << "here2" << std::endl;
                    last_protein_atom = target_atom;
                }
            }
        }
    }
    residueTarget->SetAssembly(assemblyTarget);
    assemblyTarget->AddResidue(residueTarget);
    assemblyTarget->BuildStructureByDistance();

    //std::cout << "Here now" << std::endl;
   //superimposition_atoms_.Print();
   // std::cout << "How now" << std::endl;
   // glycan_.Print();

    gmml::Superimpose(&superimposition_atoms_, assemblyTarget, &glycan_);

    Residue* reducing_Residue = glycan_.GetAllResiduesOfAssembly().at(1); // I assume I assumed something stupid here.
    //std::cout << "Reducing residue is " << reducing_Residue->GetName() << std::endl;
    AtomVector reducing_Atoms = reducing_Residue->GetAtoms();
    Atom* atomC1;
    for(AtomVector::iterator it = reducing_Atoms.begin(); it != reducing_Atoms.end(); it++)
    {
       Atom* atom = *it;
       if(atom->GetName().compare("C1")==0)
           atomC1 = atom;
    }

    last_protein_atom->GetNode()->AddNodeNeighbor(atomC1);
    atomC1->GetNode()->AddNodeNeighbor(last_protein_atom);

   // std::cout << "Brown cow" << std::endl;

    //AtomVector newSideChainAtoms = alternate_sidechain_.GetAllAtomsOfAssembly();
    //glycosite_residue->ReplaceAtomCoordinates(&newSideChainAtoms);
}

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

void GlycosylationSite::SetGlycan(Assembly glycan)
{
    glycan_ = glycan;
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
