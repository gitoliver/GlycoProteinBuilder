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
    AttachedRotamer(Assembly *glycan_rotamer);
    ~AttachedRotamer(); // destructor

    //////////////////////////////////////////////////////////
    //                       ACCESSOR                       //
    //////////////////////////////////////////////////////////

    Assembly *GetAttachedRotamer();
    std::string GetGlycanPDBPath();
    double GetTotalOverlap(); // Recalculate and return double each time, do not store.

    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////

    AttachedRotamer(const AttachedRotamer &L);             // copy constructor
    AttachedRotamer & operator=(const AttachedRotamer &L); // assignment

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

    Assembly* glycan_rotamer_;                                        /*!< A pointer to the assembly for this shape >*/
    // double total_overlap_;                                            /*!< Total vdW overlap with other molecules >*
    std::string glycan_pdb_path_;
    OverlapVector overlaps_;
};

#endif // ATTACHEDROTAMER_H
