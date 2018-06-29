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

void resolve_overlaps::weighted_protein_global_overlap_monte_carlo(GlycosylationSiteVector &glycosites)
{
    /* First attempt to resolve protein overlaps with scaled random walk.
     * Track global overlap for all sites.
     * Weigh protein overlaps 5x over glycan overlaps.
     * Move a site and monte carlo accept movement based on global overlap.
     */

    //Resolve_Overlaps_Random_Walk_Scaled_To_Overlap(glycosites, "protein", 500);

    double strict_tolerance = 0.1, loose_tolerance = 1.0;

    GlycosylationSitePointerVector sites_with_overlaps = DetermineSitesWithOverlap(glycosites, strict_tolerance, "total");

    std::cout << "Initial torsions and overlaps:\n";
    //glycoprotein_builder::CalculateOverlaps(glycosites);
    glycoprotein_builder::PrintDihedralAnglesAndOverlapOfGlycosites(glycosites);

    int cycle = 0;
    int max_cycles = 50;
    bool stop = false;
    double previous_chi1;
    double previous_chi2;
//    double previous_overlap, new_overlap;
    double previous_glycan_overlap, new_glycan_overlap, previous_protein_overlap, new_protein_overlap;
//    double lowest_global_overlap = glycoprotein_builder::GetGlobalOverlap(glycosites);
//    double new_global_overlap;
    bool accept_change;

    while ( (cycle < max_cycles) && (stop == false) )
    {
        ++cycle;
        std::cout << "Cycle " << cycle << " of " << max_cycles << std::endl;
        std::random_shuffle (sites_with_overlaps.begin(), sites_with_overlap.end());
        for(GlycosylationSitePointerVector::iterator it1 = sites_with_overlaps.begin(); it1 != sites_with_overlaps.end(); ++it1)
        {
            GlycosylationSite *current_glycosite = (*it1);
           // std::cout << "Checking " << current_glycosite->GetResidue()->GetId() << "\n";
            //previous_overlap = current_glycosite->GetWeightedOverlap(1.0, 5.0);
            previous_glycan_overlap = current_glycosite->GetGlycanOverlap();
            previous_protein_overlap = current_glycosite->GetProteinOverlap();
            previous_chi1 = current_glycosite->GetChi1Value();
            current_glycosite->SetChi1Value(RandomAngle_360range());
            previous_chi2 = current_glycosite->GetChi2Value();
            current_glycosite->SetChi2Value(RandomAngle_360range());
            new_glycan_overlap = current_glycosite->Calculate_bead_overlaps_noRecord_noSet("glycan");
            new_protein_overlap = current_glycosite->Calculate_bead_overlaps_noRecord_noSet("protein");

//            previous_overlap = current_glycosite->GetOverlap();
//            previous_chi1 = current_glycosite->GetChi1Value();
//            previous_chi2 = current_glycosite->GetChi2Value();
//            current_glycosite->SetChi1Value(RandomAngle_360range());
//            current_glycosite->SetChi2Value(RandomAngle_360range());
//            new_overlap = current_glycosite->Calculate_bead_overlaps_noRecord_noSet("total");

            //new_overlap = current_glycosite->GetWeightedOverlap(1.0, 5.0);
//            accept_change = monte_carlo::accept_via_metropolis_criterion((new_glycan_overlap + (new_protein_overlap*2)) - (previous_glycan_overlap + (previous_protein_overlap*2)));
//            if (!accept_change)
//            if (new_overlap >= previous_overlap)
            // Added weights to emphasis protein overlap as more important to relieve
            if ((new_glycan_overlap + (new_protein_overlap*2)) >= (previous_glycan_overlap + (previous_protein_overlap*2)))
            {
                current_glycosite->SetChi1Value(previous_chi1);
                current_glycosite->SetChi2Value(previous_chi2);
            }
//            else
//            {
//                std::cout << "Accepted a change of " << ((new_overlap) - (previous_overlap)) << "\n";
//                current_glycosite->Calculate_bead_overlaps();
//            }
        }
        //write_pdb_file(glycosites.at(0).GetGlycoprotein(), cycle, "test", 5.0);
        // Only need to write out when at lowest when doing metropolis where I might end up higher in overlap.
        //CalculateOverlaps(glycosites);
//        new_global_overlap = glycoprotein_builder::GetGlobalOverlap(glycosites);
//        std::cout << "Lowest: " << lowest_global_overlap << ", Current: " << new_global_overlap << "\n";
//        if ( lowest_global_overlap > new_global_overlap)
//        {
//            //Assembly *glycoprotein = glycosites.at(0).GetGlycoprotein();
//            write_pdb_file(glycosites.at(0).GetGlycoprotein(), cycle, "best", new_global_overlap);
//            lowest_global_overlap = new_global_overlap;
//        }

        //std::cout << "Updating list of sites with overlaps." << std::endl;
        sites_with_overlaps = DetermineSitesWithOverlap(glycosites, strict_tolerance); // Moved glycans may clash with other glycans. Need to check.
        if (sites_with_overlaps.size() == 0)
        {
            std::cout << "Stopping with all overlaps resolved.\n";
            stop = true;
        }
    }
    std::cout << "Global overlap before deleting sites is " << glycoprotein_builder::GetGlobalOverlap(glycosites) << "\n";
    std::cout << "Finished torsions and overlaps:\n";
    glycoprotein_builder::PrintDihedralAnglesAndOverlapOfGlycosites(glycosites);
    DeleteSitesIterativelyWithOverlapAboveTolerance(glycosites, loose_tolerance);
}

void resolve_overlaps::Resolve_Overlaps_Random_Walk_Scaled_To_Overlap(GlycosylationSiteVector &glycosites, std::string type, int max_cycles, double strict_tolerance, double loose_tolerance)
{
    //double strict_tolerance = 0.1, loose_tolerance = 1.0; // Aim for <0.1 when resolving, but keep any less than 1 when culling.
    GlycosylationSitePointerVector sites_with_overlaps = DetermineSitesWithOverlap(glycosites, strict_tolerance, type);
    Overlap_Weighted_Adjust_Torsions_For_X_Cycles(sites_with_overlaps, glycosites, max_cycles, strict_tolerance, type);
    DeleteSitesWithOverlapRecordsAboveTolerance(glycosites, loose_tolerance, type);
    sites_with_overlaps = DetermineSitesWithOverlap(glycosites, strict_tolerance, type);
    SetBestChi1Chi2(sites_with_overlaps, type);
}

void resolve_overlaps::protein_first_random_walk_scaled_to_overlap(GlycosylationSiteVector &glycosites)
{
    // NOTE!!! Upon testing, this algorithm fails with densely clustered glycans. It remembers good spots for individual glycans, but they can be in the same place.

    /* Algorithm overview:
     *  Resolve all protein overlap first, reject sites that cannot be resolved.
     *  Calculate protein overlaps, for each site with overlaps
     *        Change chi1 and chi2 values, record values that reduce overlaps
     *  Once finished delete any sites that could not be resolved
     *  Calculate total (protein+glycan) overlaps, for each site with overlaps
     *        Change chi1 and chi2 values, record values that reduce overlaps
     *  Once finished delete any sites that could not be resolved
     */

//    std::cout << "Initial overlaps for all sites\n";
//    std::cout << "      Site        |  Total | Protein | Glycan \n";
//    for(GlycosylationSiteVector::iterator it1 = glycosites.begin(); it1 != glycosites.end(); ++it1)
//    {
//        GlycosylationSite current_glycosite = *it1;
//        current_glycosite.Calculate_and_print_bead_overlaps();
//    }

     Resolve_Overlaps_Random_Walk_Scaled_To_Overlap(glycosites, "protein");

//    std::cout << "Overlaps for all sites after protein resolution and deletion\n";
//    std::cout << "      Site        |  Total | Protein | Glycan \n";
//    for(GlycosylationSiteVector::iterator it1 = glycosites.begin(); it1 != glycosites.end(); ++it1)
//    {
//        GlycosylationSite current_glycosite = *it1;
//        current_glycosite.Calculate_and_print_bead_overlaps();
//    }

     Resolve_Overlaps_Random_Walk_Scaled_To_Overlap(glycosites, "total");
    // glycoprotein_builder::PrintDihedralAnglesOfGlycosites(glycosites);
    //PrintOverlaps(glycosites);
}


bool resolve_overlaps::dumb_random_walk(GlycosylationSiteVector &glycosites)
{
    /* Algorithm:
     * Determine which sites have overlaps greater than tolerance. Stop if zero sites.
     * For each site with overlaps:
     *        Randomly change all chi1 and chi2 values
     */
    double tolerance = 0.1;
    int cycle = 0, max_cycles = 10;
    GlycosylationSitePointerVector sites_with_overlaps = DetermineSitesWithOverlap(glycosites, tolerance);
    bool resolved = false;

    while ( (cycle < max_cycles) && (resolved == false) )
    {
        ++cycle;
        std::cout << "Cycle " << cycle << " of " << max_cycles << std::endl;
        for(GlycosylationSitePointerVector::iterator it1 = sites_with_overlaps.begin(); it1 != sites_with_overlaps.end(); ++it1)
        {
            GlycosylationSite *current_glycosite = (*it1);
            current_glycosite->SetChi1Value(RandomAngle_360range());
            current_glycosite->SetChi2Value(RandomAngle_360range());
            //double percent_overlap = ((current_glycosite->GetTotalOverlap() / (current_glycosite->GetAttachedGlycan()->GetAllAtomsOfAssembly().size()) ) + 0.01);
            //  new_dihedral_value = Angle_PlusMinusX(current_glycosite->GetChi1Value(), (180 * percent_overlap) ); // scaled to degree of overlap
            //        sites_with_overlaps.erase(std::remove(sites_with_overlaps.begin(), sites_with_overlaps.end(), *it1), sites_with_overlaps.end());
        }
        //std::cout << "Updating list of sites with overlaps." << std::endl;
        sites_with_overlaps = DetermineSitesWithOverlap(glycosites, tolerance); // Moved glycans may clash with other glycans. Need to check.
        if (sites_with_overlaps.size() == 0)
        {
            std::cout << "Stopping with all overlaps resolved.\n";
            resolved = true;
        }
    }
    return resolved;
}

// random number generator; allows full range rotation
double RandomAngle_360range()
{
    return (rand() % 360) + 1 - 180;
}

double RandomAngle_range(int min, int max)
{
   // double angle = rand() % (max + 1 - min) + min;
    //std::cout << "Angle in range " << min << " - " << max << " is " << angle << "\n";
    // return angle;
    return rand() % (max + 1 - min) + min;
}

// random number generator; specify a maximum step size relative to a start point
double Angle_PlusMinusX(double start_point, int max_step_size)
{
    return start_point + (rand() % (max_step_size * 2) + 1) - max_step_size;
}

double GetNewAngleScaledToDegreeOfOverlap(double current_angle, double overlap, int number_of_atoms)
{
    int max_step_size = 1 + std::round( 180 * ( overlap / number_of_atoms ) ); // Always allow at least 1 degrees of movement
    return Angle_PlusMinusX(current_angle, max_step_size);
}

void write_pdb_file(Assembly *glycoprotein, int cycle, std::string summary_filename, double score)
{
    std::stringstream ss;
    ss << summary_filename << "_cycle_" << cycle << "overlap_" << score << ".pdb";
    PdbFileSpace::PdbFile *outputPdbFile = glycoprotein->BuildPdbFileStructureFromAssembly(-1,0);
    outputPdbFile->Write(ss.str());
//    std::ofstream summary;   // write a file that describes the best conformations found
//    summary.open(summary_filename, std::ios::out | std::ios::app);
//    summary << score << "\t" << "cycle_" << cycle << ".pdb\n";
//    summary.close();
}

void PrintOverlaps(GlycosylationSiteVector &glycosites)
{
    std::cout << "      Site        |  Total | Protein | Glycan \n";
    for (GlycosylationSiteVector::iterator current_glycosite = glycosites.begin(); current_glycosite != glycosites.end(); ++current_glycosite)
    {
        current_glycosite->Print_bead_overlaps();
    }
}

void PrintOverlaps(GlycosylationSitePointerVector &glycosites)
{
    std::cout << "      Site        |  Total | Protein | Glycan \n";
    for (GlycosylationSitePointerVector::iterator current_glycosite = glycosites.begin(); current_glycosite != glycosites.end(); ++current_glycosite)
    {
        (*current_glycosite)->Print_bead_overlaps(); // Good old pointers to pointers.
    }
}



void CalculateAndPrintOverlaps(GlycosylationSiteVector &glycosites)
{
    std::cout << "      Site        |  Total | Protein | Glycan \n";
    for (GlycosylationSiteVector::iterator current_glycosite = glycosites.begin(); current_glycosite != glycosites.end(); ++current_glycosite)
    {
        current_glycosite->Calculate_and_print_bead_overlaps();
    }
}

void CalculateOverlaps(GlycosylationSiteVector &glycosites)
{
    for (GlycosylationSiteVector::iterator current_glycosite = glycosites.begin(); current_glycosite != glycosites.end(); ++current_glycosite)
    {
        current_glycosite->Calculate_bead_overlaps();
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
//            std::cout << "Site " << current_glycosite->GetResidue()->GetId() << " is over tolerance with " << overlap << "\n";
            current_glycosite->Print_bead_overlaps();
            sites_to_return.push_back(&(*current_glycosite));
        }
    }
    return sites_to_return;
}

GlycosylationSitePointerVector GetSitesWithOverlap(GlycosylationSiteVector &glycosites, double tolerance)
{
    GlycosylationSitePointerVector sites_to_return;
    double overlap = 0.0;
    std::cout << "      Site        |  Total | Protein | Glycan " << std::endl;
    for (GlycosylationSiteVector::iterator current_glycosite = glycosites.begin(); current_glycosite != glycosites.end(); ++current_glycosite)
    {
        overlap = current_glycosite->GetOverlap();
        if ( overlap > tolerance)
        {
//            std::cout << "Site " << current_glycosite->GetResidue()->GetId() << " is over tolerance with " << overlap << "\n";
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

void DeleteSitesIterativelyWithOverlapAboveTolerance(GlycosylationSiteVector &glycosites, double tolerance)
{
    bool continue_deleting = true;
    // While overlap for any site is > tolerance
    // Delete site with highest overlap.
    // Re-calculate overlaps.
    while (continue_deleting)
    {
        GlycosylationSite *worst_site = glycosites.data(); // Pointer to the first glycosite. Remember an erase/remove "advances"
        for (GlycosylationSiteVector::iterator current_glycosite = glycosites.begin(); current_glycosite != glycosites.end(); ++current_glycosite)
        {
            if  ( current_glycosite->GetOverlap() > worst_site->GetOverlap())
            {
                worst_site = &(*current_glycosite); // The C is strong with this one.
            }
        }
        if (worst_site->GetOverlap() > tolerance)
        {
            continue_deleting = true;
            std::cout << "Site " << worst_site->GetResidueNumber() << ": " << worst_site->GetOverlap() << " :" << "Removed\n";
            ResidueVector glycan_residues = worst_site->GetAttachedGlycan()->GetResidues();
            for(ResidueVector::iterator it = glycan_residues.begin(); it != glycan_residues.end(); ++it)
            {
                worst_site->GetGlycoprotein()->RemoveResidue(*it);
            }
            glycosites.erase(std::remove(glycosites.begin(), glycosites.end(), *worst_site), glycosites.end()); // Note need #include <algorithm>
            Set_Other_Glycan_Beads(glycosites); // After "erasing", the actual atoms still exist and pointers to them are valid. Need to reset what beads are part of "other".
            glycoprotein_builder::CalculateOverlaps(glycosites); // After deleting, other sites will have lower values
        }
        else
        {
            continue_deleting = false;
        }
    }
    return;
}

void Overlap_Weighted_Adjust_Torsions_For_X_Cycles(GlycosylationSitePointerVector &sites, GlycosylationSiteVector &glycosites, int max_cycles, double tolerance, std::string overlap_type)
{
    int cycle = 0;
    bool stop = false;
    while ( (cycle < max_cycles) && (stop == false) )
    {
        ++cycle;
        std::cout << "Cycle " << cycle << " of " << max_cycles << ".\n";
        Overlap_Weighted_Adjust_Torsions(sites);
        std::cout << "Updating list of sites with overlaps.\n";
        sites = DetermineSitesWithOverlap(glycosites, tolerance, overlap_type);
        if (sites.size() == 0)
        {
            std::cout << "No more overlaps\n";
            stop = true;
        }
    }
    return;
}

void Overlap_Weighted_Adjust_Torsions(GlycosylationSitePointerVector &sites)
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
