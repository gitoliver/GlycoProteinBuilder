

#include "resolve_overlaps.h"
#include "glycosylationsite.h"


typedef std::vector<GlycosylationSite*> GlycosylationSitePointerVector;

using namespace gmml;

// random number generator; allows full range rotation
double RandomAngle_360range()
{
    return (rand() % 360) + 1 - 180;
}

// random number generator; specify a maximum step size relative to a start point
double RandomAngle_PlusMinusX(double start_point, int max_step_size)
{
    return start_point + (rand() % (max_step_size * 2)) - max_step_size;
}

void write_pdb_file(Assembly *glycoprotein, int cycle, std::string summary_filename, double score)
{
    std::string pdb_filename = "outputs/pose_" + std::to_string(cycle) + ".pdb";
    PdbFileSpace::PdbFile *outputPdbFile = glycoprotein->BuildPdbFileStructureFromAssembly(-1,0);
    outputPdbFile->Write(pdb_filename);
    std::ofstream summary;   // write a file that describes the best conformations found
    summary.open(summary_filename, std::ios::out | std::ios::app);
    summary << score << "\t" << "pose_" << cycle << ".pdb\n";
    summary.close();
}
GlycosylationSitePointerVector DetermineSitesWithOverlap(GlycosylationSiteVector *glycosites, double tolerance)
{
    GlycosylationSitePointerVector sites_with_overlaps;
    std::cout << "      Site        | Total | Protein | Glycan " << std::endl;
    for (GlycosylationSiteVector::iterator current_glycosite = glycosites->begin(); current_glycosite != glycosites->end(); ++current_glycosite)
    {
        if (current_glycosite->Calculate_and_print_bead_overlaps() > tolerance)
        {
           // std::cout << current_glycosite->GetResidue()->GetId() << " added to sites_with_overlaps." << std::endl;
            sites_with_overlaps.push_back(&(*current_glycosite));
        }
    }
    return sites_with_overlaps;
}

/* New algorithm idea:
 * For each site that has overlap:
 1) Try to resolve all protein overlaps, ignoring glycan-glycan. Delete glycan from any site with unresolvable protein overlaps, report to user.
 2) Resolve glycan overlaps. Either
    a) Loop through and check every site each cycle, and move those with overlaps

    b) Keep list of glycans with overlaps, move those.
       Check total_overlaps after move, remove any sites with no overlap from list of overlapping sites.
       After moving all sites (X times?), check glycan overlaps for those sites that didn't have overlaps before move. Add any newly overlapping sites to list.
       Update no_overlaps list with sites that now have no overlap.

    c) Create a tree structure for sites that overlap with each other.
        Could use assembly index/id of bead atom? Can work on individual trees and report nice info to users.
 3)

*/
void resolve_overlaps::dumb_monte_carlo(Assembly *glycoprotein, GlycosylationSiteVector *glycosites)
{
    /* Algorithm:
     * Determine which sites have overlaps
     * For each site with overlaps:
     *        Randomly change all chi1 and chi2 values
     */
    double tolerance = 0.1;
    int cycle = 1, max_cycles = 10;
    GlycosylationSitePointerVector sites_with_overlaps = DetermineSitesWithOverlap(glycosites, tolerance);
    bool stop = false;

    while ( (cycle <= max_cycles) && (stop == false) )
    {
        std::cout << "Cycle " << cycle << " of " << max_cycles << std::endl;
        for(GlycosylationSitePointerVector::iterator it1 = sites_with_overlaps.begin(); it1 != sites_with_overlaps.end(); ++it1)
        {
            GlycosylationSite *current_glycosite = (*it1);
            current_glycosite->SetChi1Value(RandomAngle_360range(), glycoprotein);
            current_glycosite->SetChi2Value(RandomAngle_360range(), glycoprotein);
            //double percent_overlap = ((current_glycosite->GetTotalOverlap() / (current_glycosite->GetAttachedGlycan()->GetAllAtomsOfAssembly().size()) ) + 0.01);
            //  new_dihedral_value = RandomAngle_PlusMinusX(current_glycosite->GetChi1Value(), (180 * percent_overlap) ); // scaled to degree of overlap
            //        sites_with_overlaps.erase(std::remove(sites_with_overlaps.begin(), sites_with_overlaps.end(), *it1), sites_with_overlaps.end());
        }
        std::cout << "Updating list of sites with overlaps." << std::endl;
        sites_with_overlaps = DetermineSitesWithOverlap(glycosites, tolerance); // Moved glycans may clash with other glycans. Need to check.
        if (sites_with_overlaps.size() == 0)
        {
            std::cout << "Stapping" << std::endl;
            stop = true;
        }
        ++cycle;
    }
    write_pdb_file(glycoprotein, 1, "./outputs/summary", 0.0);
}

void resolve_overlaps::example_for_Gordon(Assembly *glycoprotein, GlycosylationSiteVector *glycosites)
{
    double site_total_overlap = 0.0, site_glycan_overlap = 0.0, site_protein_overlap = 0.0, new_dihedral_value = 0.0, total_system_overlap = 0.0;
    std::cout << "      Site        | Total | Protein | Glycan " << std::endl;

    for (GlycosylationSiteVector::iterator current_glycosite = glycosites->begin(); current_glycosite != glycosites->end(); ++current_glycosite)
    {
        site_total_overlap = current_glycosite->Calculate_bead_overlaps(); // Must repeat after rotating chi1, chi2.
        site_glycan_overlap = current_glycosite->GetGlycanOverlap(); // If you wish to have this level of detail
        site_protein_overlap = current_glycosite->GetProteinOverlap(); // If you wish to have this level of detail
        current_glycosite->Print_bead_overlaps();
        new_dihedral_value = RandomAngle_PlusMinusX(current_glycosite->GetChi1Value(), 6);
        //std::cout << "Changing chi1 from " << current_glycosite->GetChi1Value() << " to " << new_dihedral_value << std::endl;
        current_glycosite->SetChi1Value(new_dihedral_value, glycoprotein);
        new_dihedral_value = RandomAngle_PlusMinusX(current_glycosite->GetChi2Value(), 6);
        current_glycosite->SetChi2Value(new_dihedral_value, glycoprotein);
        current_glycosite->Calculate_and_print_bead_overlaps();
        total_system_overlap += site_total_overlap;
    }
    write_pdb_file(glycoprotein, 1, "./outputs/summary", total_system_overlap);
}

void resolve_overlaps::genetic_algorithm(Assembly *glycoprotein, GlycosylationSiteVector *glycosites)
{
    // Gordon.
}
