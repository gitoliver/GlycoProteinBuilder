#ifndef ATTACHEDROTAMER_H
#define ATTACHEDROTAMER_H

#include "overlap.h"


class AttachedRotamer
{

public:
    //////////////////////////////////////////////////////////
    //                    TYPE DEFINITION                   //
    //////////////////////////////////////////////////////////

   // typedef std::vector < std::tuple<Residue*, Residue*, double> > Residue_overlap;
    typedef std::vector<Overlap*> OverlapVector;

    //////////////////////////////////////////////////////////
    //                       CONSTRUCTOR                    //
    //////////////////////////////////////////////////////////

    AttachedRotamer();
    AttachedRotamer(Assembly glycan_rotamer);
    ~AttachedRotamer(); // destructor

    //////////////////////////////////////////////////////////
    //                       ACCESSOR                       //
    //////////////////////////////////////////////////////////

    Assembly *GetAttachedRotamer();
    Assembly* GetSuperimpositionAtoms();
    Assembly* GetAlternateSidechain();
    double GetTotalOverlap(); // Recalculate and return double each time, do not store.

    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////

    AttachedRotamer(const AttachedRotamer &L);             // copy constructor
    AttachedRotamer & operator=(const AttachedRotamer &L); // assignment

    void Prepare_Glycans_For_Superimposition_To_Particular_Residue(std::string amino_acid_name);
    void Superimpose_Glycan_To_Glycosite(Residue *glycosite_residue);
    //////////////////////////////////////////////////////////
    //                       MUTATOR                        //
    //////////////////////////////////////////////////////////

    void SetAttachedRotamer(Assembly *glycan_assembly);
    void SetGlycanPDBPath(std::string glycan_pdb_path);
    void AddOverlap(Overlap *overlap);
    void AddOverlap(Residue *residue1, Residue *residue2, double overlap);

    //////////////////////////////////////////////////////////
    //                       DISPLAY FUNCTION               //
    //////////////////////////////////////////////////////////

    //void Print(std::ostream& out = std::cout);
     void Print(std::ostream &output) const; // print the list to output

private:
    //////////////////////////////////////////////////////////
    //                       ATTRIBUTES                     //
    //////////////////////////////////////////////////////////

    Assembly glycan_rotamer_;                                            /*!< An assembly for this shape >*/
    Assembly superimposition_atoms_;                                     /*!< The 3 atoms used for superimposition of glycan to sidechain >*/
    Assembly alternate_sidechain_;                                       /*!< Not sure. I think I was going to keep the original sidechain position and only alter the alternative >*/
    OverlapVector overlaps_;
};

#endif // ATTACHEDROTAMER_H
