#include "../includes/resolve_overlaps.h"

typedef std::vector<GlycosylationSite*> GlycosylationSitePointerVector;

using namespace gmml;

/* New algorithm idea:
 * For each site that has overlap:
 1) Try to resolve all protein overlaps, ignoring glycan-glycan. Delete glycan from any site with unresolvable protein overlaps, report to user.
 2) Resolve glycan overlaps. Either
    a) Loop through and check every site each cycle, and move those with overlaps

    c) Create a tree structure for sites that overlap with each other.
        Could use assembly index/id of bead atom? Can work on individual trees and report nice info to users.


Note: Chi1 is 180,-60,60 +/- 30 degrees. Want a function that keeps the values within these states.
Note: Need to save best structure.
*/


void resolve_overlaps::protein_first_monte_carlo(GlycosylationSiteVector &glycosites)
{
    /* Algorithm plan:
     *  Resolve all protein overlap first, reject sites that cannot be resolved. Record cycles to resolve?
     *  Calculate protein overlaps, for each site with overlaps
     *        Randomly change all chi1 and chi2 values
     *  Once finished;
     *  Calculate glycan overlaps, for each site with overlaps
     *        Change chi1 and chi2, scaled by degree of overlap
     *        Reject changes causing protein overlaps
     */
    double strict_tolerance = 0.1, loose_tolerance = 1.0; // Aim for <0.1 when resolving, but keep any less than 1 when culling.
    int max_cycles = 500;
    GlycosylationSitePointerVector sites_with_protein_overlaps = DetermineSitesWithOverlap(glycosites, strict_tolerance, "protein");
    Monte_Carlo_Torsions(sites_with_protein_overlaps, glycosites, max_cycles, strict_tolerance, "protein");
    GlycosylationSitePointerVector temp_for_Gordon = DetermineSitesWithOverlap(glycosites, strict_tolerance, "total");
    DeleteSitesWithOverlapRecordsAboveTolerance(glycosites, loose_tolerance, "protein");
    sites_with_protein_overlaps = DetermineSitesWithOverlap(glycosites, strict_tolerance, "protein");
    SetBestChi1Chi2(sites_with_protein_overlaps, "protein");
    GlycosylationSitePointerVector remaining_sites = DetermineSitesWithOverlap(glycosites, strict_tolerance, "total");
    Monte_Carlo_Torsions(remaining_sites, glycosites, max_cycles, strict_tolerance, "total");
    SetBestChi1Chi2(remaining_sites);
    //GlycosylationSitePointerVector resolved_sites = DeleteSitesWithOverlaps(glycosites, loose_tolerance, "total");
    //PrintOverlaps(glycosites);
}

void resolve_overlaps::dumb_monte_carlo(GlycosylationSiteVector &glycosites)
{
    /* Algorithm:
     * Determine which sites have overlaps greater than tolerance. Stop if zero sites.
     * For each site with overlaps:
     *        Randomly change all chi1 and chi2 values
     */
    double tolerance = 0.1;
    int cycle = 0, max_cycles = 5;
    GlycosylationSitePointerVector sites_with_overlaps = DetermineSitesWithOverlap(glycosites, tolerance);
    bool stop = false;

    while ( (cycle < max_cycles) && (stop == false) )
    {
        ++cycle;
        std::cout << "Cycle " << cycle << " of " << max_cycles << std::endl;
        for(GlycosylationSitePointerVector::iterator it1 = sites_with_overlaps.begin(); it1 != sites_with_overlaps.end(); ++it1)
        {
            GlycosylationSite *current_glycosite = (*it1);
            current_glycosite->SetChi1Value(RandomAngle_360range());
            current_glycosite->SetChi2Value(RandomAngle_360range());
            //double percent_overlap = ((current_glycosite->GetTotalOverlap() / (current_glycosite->GetAttachedGlycan()->GetAllAtomsOfAssembly().size()) ) + 0.01);
            //  new_dihedral_value = RandomAngle_PlusMinusX(current_glycosite->GetChi1Value(), (180 * percent_overlap) ); // scaled to degree of overlap
            //        sites_with_overlaps.erase(std::remove(sites_with_overlaps.begin(), sites_with_overlaps.end(), *it1), sites_with_overlaps.end());
        }
        //std::cout << "Updating list of sites with overlaps." << std::endl;
        sites_with_overlaps = DetermineSitesWithOverlap(glycosites, tolerance); // Moved glycans may clash with other glycans. Need to check.
        if (sites_with_overlaps.size() == 0)
        {
            std::cout << "Stopping with all overlaps resolved.\n";
            stop = true;
        }
    }
}


// random number generator; allows full range rotation
double RandomAngle_360range()
{
    return (rand() % 360) + 1 - 180;
}

// random number generator; specify a maximum step size relative to a start point
double RandomAngle_PlusMinusX(double start_point, int max_step_size)
{
    return start_point + (rand() % (max_step_size * 2) + 1) - max_step_size;
}

double GetNewAngleScaledToDegreeOfOverlap(double current_angle, double overlap, int number_of_atoms)
{
    int max_step_size = 1 + std::round( 180 * ( overlap / number_of_atoms ) ); // Always allow at least 1 degrees of movement
    return RandomAngle_PlusMinusX(current_angle, max_step_size);
}

void write_pdb_file(Assembly &glycoprotein, int cycle, std::string summary_filename, double score)
{
    std::string pdb_filename = "outputs/pose_" + std::to_string(cycle) + ".pdb";
    PdbFileSpace::PdbFile *outputPdbFile = glycoprotein.BuildPdbFileStructureFromAssembly(-1,0);
    outputPdbFile->Write(pdb_filename);
    std::ofstream summary;   // write a file that describes the best conformations found
    summary.open(summary_filename, std::ios::out | std::ios::app);
    summary << score << "\t" << "pose_" << cycle << ".pdb\n";
    summary.close();
}

void PrintOverlaps(GlycosylationSiteVector &glycosites)
{
    for (GlycosylationSiteVector::iterator current_glycosite = glycosites.begin(); current_glycosite != glycosites.end(); ++current_glycosite)
    {
        current_glycosite->Print_bead_overlaps();
    }
}

void PrintOverlaps(GlycosylationSitePointerVector &glycosites)
{
    for (GlycosylationSitePointerVector::iterator it = glycosites.begin(); it != glycosites.end(); ++it)
    {
        GlycosylationSite *current_glycosite = (*it);
        current_glycosite->Print_bead_overlaps();
    }
}

void SetBestChi1Chi2(GlycosylationSitePointerVector &glycosites, std::string overlap_type)
{
    /****
     * Ok, the problem here is that two sites that are close may independantly find a best set that when set together causes them to overlap. This can happen
     * when just looking at protein overlaps, as the sfat atom on the NLN is moving around. It's worse when you consider glycan overlaps.
     * When searching for the best overlaps, the monte_carlo function must always treat new low/equal overlaps as better. If it found zero for a site,
     * overlaps are reintroduced by another site moving, then when it finds a new zero for this site, the new one should take precedence. Going to check
     * that now.
     */
    for (GlycosylationSitePointerVector::iterator it = glycosites.begin(); it != glycosites.end(); ++it)
    {
        GlycosylationSite *current_glycosite = (*it);
        current_glycosite->SetChi1Value(current_glycosite->GetBestOverlapRecord(overlap_type).GetChi1());
        current_glycosite->SetChi2Value(current_glycosite->GetBestOverlapRecord(overlap_type).GetChi2());
        current_glycosite->Calculate_bead_overlaps();
    }

}

GlycosylationSitePointerVector DetermineSitesWithOverlap(GlycosylationSiteVector &glycosites, double tolerance, std::string overlap_type)
{
    GlycosylationSitePointerVector sites_to_return;
    double overlap = 0.0;
    std::cout << "      Site        |  Total | Protein | Glycan " << std::endl;
    for (GlycosylationSiteVector::iterator current_glycosite = glycosites.begin(); current_glycosite != glycosites.end(); ++current_glycosite)
    {
        overlap = current_glycosite->Calculate_bead_overlaps(overlap_type);
        if ( overlap > tolerance)
        {
            current_glycosite->Print_bead_overlaps();
            sites_to_return.push_back(&(*current_glycosite));
        }
    }
    return sites_to_return;
}


GlycosylationSitePointerVector DeleteSitesWithOverlaps(GlycosylationSiteVector &glycosites, double tolerance, std::string overlap_type)
{
    GlycosylationSitePointerVector sites_to_return;
    double overlap = 0.0;
    for (GlycosylationSiteVector::iterator current_glycosite = glycosites.begin(); current_glycosite != glycosites.end();)
    {
        overlap = current_glycosite->Calculate_bead_overlaps(overlap_type);
        // Delete site from list if overlap is greater than the tolerance value
        std::cout << "Site " << current_glycosite->GetResidueNumber() << ": " << overlap << " :";
        if ( overlap > tolerance)
        {
            std::cout << "Removed\n";
            ResidueVector glycan_residues = current_glycosite->GetAttachedGlycan()->GetResidues();
            for(ResidueVector::iterator it = glycan_residues.begin(); it != glycan_residues.end(); ++it)
            {
                current_glycosite->GetGlycoprotein()->RemoveResidue(*it);
                //glycoprotein.RemoveResidue(*it);
            }
            glycosites.erase(std::remove(glycosites.begin(), glycosites.end(), *current_glycosite), glycosites.end()); // Note need #include <algorithm>
        }
        else
        {
            std::cout << "Retained\n";
            sites_to_return.push_back(&(*current_glycosite));
            ++current_glycosite; // This will get you. Erase/Remove advances current_glycosite.
        }
    }
    return sites_to_return;
}

void DeleteSitesWithOverlapRecordsAboveTolerance(GlycosylationSiteVector &glycosites, double tolerance, std::string overlap_type)
{
    for (GlycosylationSiteVector::iterator current_glycosite = glycosites.begin(); current_glycosite != glycosites.end();)
    {
        double overlap = current_glycosite->GetBestOverlapRecord(overlap_type).GetOverlap();
        std::cout << "Site " << current_glycosite->GetResidueNumber() << ": " << overlap << " :";
        if ( overlap > tolerance)
        {
            std::cout << "Removed\n";
            ResidueVector glycan_residues = current_glycosite->GetAttachedGlycan()->GetResidues();
            for(ResidueVector::iterator it = glycan_residues.begin(); it != glycan_residues.end(); ++it)
            {
                current_glycosite->GetGlycoprotein()->RemoveResidue(*it);
            }
            glycosites.erase(std::remove(glycosites.begin(), glycosites.end(), *current_glycosite), glycosites.end()); // Note need #include <algorithm>
        }
        else
        {
            std::cout << "Retained\n";
            ++current_glycosite; // This will get you. Erase/Remove advances current_glycosite.
        }
    }
}

void Monte_Carlo_Torsions(GlycosylationSitePointerVector &sites, GlycosylationSiteVector &glycosites, int max_cycles, double tolerance, std::string overlap_type)
{
    int cycle = 0;
    bool stop = false;
    while ( (cycle < max_cycles) && (stop == false) )
    {
        ++cycle;
        std::cout << "Cycle " << cycle << " of " << max_cycles << ".\n";
        RandomizeTorsions(sites);
        std::cout << "Updating list of sites with protein overlaps.\n";
        sites = DetermineSitesWithOverlap(glycosites, tolerance, overlap_type);
        if (sites.size() == 0)
        {
            std::cout << "No more protein overlaps\n";
            stop = true;
        }
    }
    return;
}

void RandomizeTorsions(GlycosylationSitePointerVector &sites)
{
    double new_dihedral_angle = 0.0;
    for(GlycosylationSitePointerVector::iterator it1 = sites.begin(); it1 != sites.end(); ++it1)
    {
        GlycosylationSite *current_glycosite = (*it1);
        new_dihedral_angle = GetNewAngleScaledToDegreeOfOverlap(current_glycosite->GetChi1Value(), current_glycosite->GetProteinOverlap(), current_glycosite->GetAttachedGlycan()->GetAllAtomsOfAssembly().size());
        current_glycosite->SetChi1Value(new_dihedral_angle);
        new_dihedral_angle = GetNewAngleScaledToDegreeOfOverlap(current_glycosite->GetChi2Value(), current_glycosite->GetProteinOverlap(), current_glycosite->GetAttachedGlycan()->GetAllAtomsOfAssembly().size());
        current_glycosite->SetChi2Value(new_dihedral_angle);
    }
    return;
}
