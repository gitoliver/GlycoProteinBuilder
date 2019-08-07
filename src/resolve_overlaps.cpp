#include "../includes/resolve_overlaps.h"

typedef std::vector<GlycosylationSite*> GlycosylationSitePointerVector;

using namespace gmml;
using namespace glycoprotein_builder;
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

// Lol if you're reading this good luck. I haven't figured out a better way than this mess:
//void resolve_overlaps::generatePermutationsWithinLinkageRecursively(double &lowest_overlap, RotatableDihedralVector::iterator currentRotatableBond, RotatableDihedralVector::iterator rotEnd, ResidueLinkageVector::iterator linkage, ResidueLinkageVector::iterator end, GlycosylationSiteVector &glycosites)
//{
//    for(int rotamerNumber = 0; rotamerNumber < currentRotatableBond->GetNumberOfRotamers(); ++rotamerNumber)
//    {
//        currentRotatableBond->SetSpecificAngleEntryUsingMetadata(false, rotamerNumber);
//     //   std::cout << linkage->GetFromThisResidue1()->GetId() << "-" << linkage->GetToThisResidue2()->GetId() << ": " << rotamerNumber <<" \n";
//        if(std::next(linkage) != end)
//        {
//            resolve_overlaps::generateLinkagePermutationsRecursively(lowest_overlap, std::next(linkage), end, glycosites);
//        }
//        if(std::next(currentRotatableBond) != rotEnd)
//        {
//            resolve_overlaps::generatePermutationsWithinLinkageRecursively(lowest_overlap, std::next(currentRotatableBond), rotEnd, linkage, end, glycosites);
//        }
//        if((std::next(linkage) == end) && (std::next(currentRotatableBond) == rotEnd) )
//        {
//            // CALCULATE OVERLAPS
//            glycoprotein_builder::CalculateOverlaps(glycosites);
//           // std::cout << "Lowest: " << lowest_overlap << ". Global: " << glycoprotein_builder::GetGlobalOverlap(glycosites) << "\n";
//          //  glycoprotein_builder::write_pdb_file(glycosites.at(0).GetGlycoprotein(), rotamerNumber, "shapes", glycoprotein_builder::GetGlobalOverlap(glycosites) );
//            if (lowest_overlap >= (glycoprotein_builder::GetGlobalOverlap(glycosites) + 0.1) )
//            {
//                std::cout << "Saving " << glycoprotein_builder::GetGlobalOverlap(glycosites) << "\n";
//                glycoprotein_builder::StashCoordinates(glycosites);
//                lowest_overlap = glycoprotein_builder::GetGlobalOverlap(glycosites);
//            }
//        }
//    }
//}

//void resolve_overlaps::generateLinkagePermutationsRecursively(double &lowest_overlap, ResidueLinkageVector::iterator linkage, ResidueLinkageVector::iterator end, GlycosylationSiteVector &glycosites)
//{
//    if(linkage->CheckIfConformer())
//    {
//        for(int shapeNumber = 0; shapeNumber < linkage->GetNumberOfShapes(); ++shapeNumber)
//        {
//            linkage->SetSpecificShapeUsingMetadata(shapeNumber);
//         //   std::cout << linkage->GetFromThisResidue1()->GetId() << "-" << linkage->GetToThisResidue2()->GetId() << ": " << shapeNumber << "\n";
//            if(std::next(linkage) != end)
//            {
//                resolve_overlaps::generateLinkagePermutationsRecursively(lowest_overlap, std::next(linkage), end, glycosites);
//            }
//            else // At the end
//            {
//                // CALCULATE OVERLAPS
//                glycoprotein_builder::CalculateOverlaps(glycosites);
//               // std::cout << "Lowest: " << lowest_overlap << ". Global: " << glycoprotein_builder::GetGlobalOverlap(glycosites) << "\n";
//             //   glycoprotein_builder::write_pdb_file(glycosites.at(0).GetGlycoprotein(), shapeNumber, "shapes", glycoprotein_builder::GetGlobalOverlap(glycosites) );
//                if (lowest_overlap >= (glycoprotein_builder::GetGlobalOverlap(glycosites) + 0.1) )
//                {
//                    std::cout << "Saving " << glycoprotein_builder::GetGlobalOverlap(glycosites) << "\n";
//                    glycoprotein_builder::StashCoordinates(glycosites);
//                    lowest_overlap = glycoprotein_builder::GetGlobalOverlap(glycosites);
//                }
//            }
//        }
//    }
//    else
//    {
//        RotatableDihedralVector rotatableDihedrals = linkage->GetRotatableDihedralsWithMultipleRotamers();
//        resolve_overlaps::generatePermutationsWithinLinkageRecursively(lowest_overlap, rotatableDihedrals.begin(), rotatableDihedrals.end(), linkage, end, glycosites);
//    }
//}

void resolve_overlaps::generateLinkagePermutationsRecursively(double &lowest_overlap, ResidueLinkageVector::iterator linkage, ResidueLinkageVector::iterator end, GlycosylationSiteVector &glycosites)
{
    for(int shapeNumber = 0; shapeNumber < linkage->GetNumberOfShapes(); ++shapeNumber)
    {
        linkage->SetSpecificShapeUsingMetadata(shapeNumber);
        //std::cout << linkage->GetFromThisResidue1()->GetId() << "-" << linkage->GetToThisResidue2()->GetId() << ": " << (shapeNumber + 1) << " of " << linkage->GetNumberOfShapes() <<  "\n";
        if(std::next(linkage) != end)
        {
            resolve_overlaps::generateLinkagePermutationsRecursively(lowest_overlap, std::next(linkage), end, glycosites);
        }
        else // At the end
        {
            // CALCULATE OVERLAPS
            glycoprotein_builder::CalculateOverlaps(glycosites);
            // std::cout << "Lowest: " << lowest_overlap << ". Global: " << glycoprotein_builder::GetGlobalOverlap(glycosites) << "\n";
            //   glycoprotein_builder::write_pdb_file(glycosites.at(0).GetGlycoprotein(), shapeNumber, "shapes", glycoprotein_builder::GetGlobalOverlap(glycosites) );
            if (lowest_overlap >= (glycoprotein_builder::GetGlobalOverlap(glycosites) + 0.01) )
            {
                std::cout << "Saving " << glycoprotein_builder::GetGlobalOverlap(glycosites) << "\n";
                glycoprotein_builder::StashCoordinates(glycosites);
                lowest_overlap = glycoprotein_builder::GetGlobalOverlap(glycosites);
            }
        }
    }
}

void resolve_overlaps::generateLinkagePermutationsRecursivelyWritePDB(GlycosylationSite glycosite, ResidueLinkageVector::iterator linkage, ResidueLinkageVector::iterator end, int previousShapeNumber)
{
    for(int shapeNumber = 0; shapeNumber < linkage->GetNumberOfShapes(); ++shapeNumber)
    {
        linkage->SetSpecificShapeUsingMetadata(shapeNumber);

        //std::cout << linkage->GetFromThisResidue1()->GetId() << "-" << linkage->GetToThisResidue2()->GetId() << ": " << (shapeNumber + 1) << " of " << linkage->GetNumberOfShapes() <<  "\n";
        if(std::next(linkage) != end)
        {
            resolve_overlaps::generateLinkagePermutationsRecursivelyWritePDB(glycosite, std::next(linkage), end, shapeNumber);
        }
        else // At the end
        {
            PdbFileSpace::PdbFile *outputPdbFileGlycoProteinAll = glycosite.GetGlycoprotein()->BuildPdbFileStructureFromAssembly(-1,0);
            std::stringstream outname;
            outname << "Rotamer_" << (previousShapeNumber + 1) << "_" << glycosite.GetResidue()->GetId() << linkage->GetFromThisResidue1()->GetId() << "-" << linkage->GetToThisResidue2()->GetId() << "_" << (shapeNumber + 1) << "of" << linkage->GetNumberOfShapes() << ".pdb";
            std::cout << "Ding: " << outname.str() << "\n";
            outputPdbFileGlycoProteinAll->Write(outname.str());

//            // CALCULATE OVERLAPS
//            glycoprotein_builder::CalculateOverlaps(glycosites);
//            // std::cout << "Lowest: " << lowest_overlap << ". Global: " << glycoprotein_builder::GetGlobalOverlap(glycosites) << "\n";
//            //   glycoprotein_builder::write_pdb_file(glycosites.at(0).GetGlycoprotein(), shapeNumber, "shapes", glycoprotein_builder::GetGlobalOverlap(glycosites) );
//            if (lowest_overlap >= (glycoprotein_builder::GetGlobalOverlap(glycosites) + 0.01) )
//            {
//                std::cout << "Saving " << glycoprotein_builder::GetGlobalOverlap(glycosites) << "\n";
//                glycoprotein_builder::StashCoordinates(glycosites);
//                lowest_overlap = glycoprotein_builder::GetGlobalOverlap(glycosites);
//            }
        }
    }
}

void resolve_overlaps::rotamer_permutator(GlycosylationSiteVector &glycosites)
{
    std::cout << "Rotamer Permutator\n";
    glycoprotein_builder::CalculateOverlaps(glycosites);
    double lowest_global_overlap = glycoprotein_builder::GetGlobalOverlap(glycosites);
    std::cout << "Initial overlap: " << lowest_global_overlap << "\n";
    // Go through each glycosite and select the first and inner 1-6 linkages
    ResidueLinkageVector allSelectedLinkages = glycoprotein_builder::GetAllFirstAnd1_6Linkages(glycosites);
   // ResidueLinkageVector allSelectedLinkages = glycoprotein_builder::GetAllFirstAnd2_XLinkages(glycosites); // Just for testing
    ResidueLinkageVector linkagesOrderedForPermutation = glycoprotein_builder::SplitLinkagesIntoPermutants(allSelectedLinkages);
    resolve_overlaps::generateLinkagePermutationsRecursively(lowest_global_overlap, linkagesOrderedForPermutation.begin(), linkagesOrderedForPermutation.end(), glycosites);
    glycoprotein_builder::SetStashedCoordinatesWithLowestOverlap(glycosites);
    glycoprotein_builder::CalculateOverlaps(glycosites);
}

void resolve_overlaps::writePDB_rotamer_permutator(GlycosylationSiteVector &glycosites)
{
    for(auto &glycosite : glycosites)
    {
        GlycosylationSiteVector fakeVector;
        fakeVector.push_back(glycosite);
        // Go through each glycosite and select the first and inner 1-6 linkages
        ResidueLinkageVector allSelectedLinkages = glycoprotein_builder::GetAllFirstAnd1_6Linkages(fakeVector);
        // ResidueLinkageVector allSelectedLinkages = glycoprotein_builder::GetAllFirstAnd2_XLinkages(glycosites); // Just for testing
        ResidueLinkageVector linkagesOrderedForPermutation = glycoprotein_builder::SplitLinkagesIntoPermutants(allSelectedLinkages);
        resolve_overlaps::generateLinkagePermutationsRecursivelyWritePDB(glycosite, linkagesOrderedForPermutation.begin(), linkagesOrderedForPermutation.end(), 0);
    }
}

// Generate each rotamer combo for two sites and check if they overlap
void resolve_overlaps::rotamerPermutatorReasonableLimits(GlycosylationSiteVector &glycosites, int maxNumberOfSitesToConsider)
{

    for(auto &glycosite : glycosites)
    {
        glycoprotein_builder::CalculateOverlaps(glycosites);
        std::cout << "Per perm Initial overlap: " << glycoprotein_builder::GetGlobalOverlap(glycosites) << "\n";
        std::cout << "Finding partners of " << glycosite.GetResidue()->GetId() << ":" << std::endl;
        GlycosylationSiteVector closestSites = glycosite.GetXClosestSitesWithinOverlapDistanceY(glycosites, maxNumberOfSitesToConsider);
//        for(auto partner : closestSites)
//        {
//            std::cout << "    " << partner.GetResidue()->GetId() << "\n";
//        }
        // Is it worth it to now generate each rotamer combo for two sites and check if they overlap?
        // GlycosylationSiteVector relevantSites = glycosite->GetSitesWithRotamersThatOverlaps(closestSites);
        std::cout << "    Checking all possible rotamers with partners\n";
        resolve_overlaps::rotamer_permutator(closestSites);
        glycoprotein_builder::CalculateOverlaps(glycosites);
        std::cout << "Post perm Current overlap: " << glycoprotein_builder::GetGlobalOverlap(glycosites) << "\n";
    }

}

void resolve_overlaps::weighted_protein_global_overlap_random_descent(GlycosylationSiteVector &glycosites, OverlapType overlapType, int max_cycles, bool monte_carlo)
{
    std::cout << "Weighted Protein, Global Overlap with Random Decent for " << max_cycles << " cycles and monte carlo set as " << std::boolalpha << monte_carlo << ".\n";
    int cycle = 1;
    bool stop = false;
    bool record_overlap = false;
    double previous_glycan_overlap, new_glycan_overlap, previous_protein_overlap, new_protein_overlap;
    double lowest_global_overlap = glycoprotein_builder::GetGlobalOverlap(glycosites);
    double new_global_overlap;
    double overlap_difference = 0.0;
    double strict_tolerance = 0.1, loose_tolerance = 1.0;

    GlycosylationSitePointerVector sites_with_overlaps = DetermineSitesWithOverlap(glycosites, strict_tolerance, overlapType);

    std::cout << "Initial torsions and overlaps:\n";
    //glycoprotein_builder::PrintDihedralAnglesAndOverlapOfGlycosites(glycosites);
    while ( (cycle < max_cycles) && (stop == false) )
    {
        std::cout << "Cycle " << cycle << "/" << max_cycles << "\n";
        ++cycle;
        std::random_shuffle (sites_with_overlaps.begin(), sites_with_overlaps.end());
        for(auto &current_glycosite : sites_with_overlaps)
        {
            // std::cout << "Checking " << current_glycosite->GetResidue()->GetId() << "\n";
            previous_glycan_overlap = current_glycosite->GetGlycanOverlap();
            previous_protein_overlap = current_glycosite->GetProteinOverlap();
            current_glycosite->SetRandomDihedralAnglesUsingMetadata();
            std::cout << "Site: " << current_glycosite->GetResidueNumber() << "\n";
            new_glycan_overlap = current_glycosite->CalculateOverlaps(overlapType, GLYCAN, record_overlap);
            new_protein_overlap = current_glycosite->CalculateOverlaps(overlapType, PROTEIN, record_overlap);
            overlap_difference = (new_glycan_overlap + (new_protein_overlap*5)) - (previous_glycan_overlap + (previous_protein_overlap*5));
            if (overlap_difference >= 0.0) // if the change made it worse
            {
                current_glycosite->ResetDihedralAngles();
            }
            else if ( (monte_carlo) && (! monte_carlo::accept_via_metropolis_criterion(overlap_difference)) )
            {
                current_glycosite->ResetDihedralAngles();
            }
            else
            {
                std::cout << "Accepted a change of " << overlap_difference << "\n";
            }
        }
        //std::cout << "Updating list of sites with overlaps." << std::endl;
        new_global_overlap = glycoprotein_builder::GetGlobalOverlap(glycosites);
        sites_with_overlaps = DetermineSitesWithOverlap(glycosites, strict_tolerance, overlapType); // Moved glycans may clash with other glycans. Need to check.
        if (sites_with_overlaps.size() == 0)
        {
            std::cout << "Stopping with all overlaps resolved and global overlap is " << new_global_overlap <<  "\n";
            stop = true;
        }
        //write_pdb_file(glycosites.at(0).GetGlycoprotein(), cycle, "current", new_global_overlap);
        if ( lowest_global_overlap > new_global_overlap + 1 )
        {
            //   std::cout << "Lowest: " << lowest_global_overlap << ", Current: " << new_global_overlap << "\n";
           // write_pdb_file(glycosites.at(0).GetGlycoprotein(), cycle, "best", new_global_overlap);
            lowest_global_overlap = new_global_overlap;
        }
    }
  //  glycoprotein_builder::DeleteSitesIterativelyWithOverlapAboveTolerance(glycosites, loose_tolerance);
    return;
}

void resolve_overlaps::wiggle(GlycosylationSiteVector &glycosites, OverlapType overlapType, int max_cycles)
{
    std::cout << "Wiggling all linkages for " << max_cycles << " cycles.\n";
    int output_pdbfile_id = 0;
    double strict_tolerance = 0.1, loose_tolerance = 1.0;
    GlycosylationSitePointerVector sites_with_overlaps = DetermineSitesWithOverlap(glycosites, strict_tolerance, overlapType);
    int cycle = 0;
    bool stop = false;
    while ( (cycle < max_cycles) && (stop == false) )
    {
        ++cycle;
        std::cout << "Cycle " << cycle << "/" << max_cycles << "\n";
        std::random_shuffle (sites_with_overlaps.begin(), sites_with_overlaps.end());
        for(auto &glycosite : sites_with_overlaps)
        {
           // std::cout << "Wiggling " << glycosite->GetResidueNumber() << " with id at " << output_pdbfile_id << "\n";
            glycosite->Wiggle(&output_pdbfile_id, strict_tolerance);
        }
        // Check which sites still have overlaps, stop if none.
        sites_with_overlaps = DetermineSitesWithOverlap(glycosites, strict_tolerance, overlapType); // Moved glycans may clash with other glycans. Need to check.
        if (sites_with_overlaps.size() == 0)
        {
            std::cout << "Stopping with all overlaps resolved.\n";
            stop = true;
        }
    }
 //   glycoprotein_builder::DeleteSitesIterativelyWithOverlapAboveTolerance(glycosites, loose_tolerance);
    return;
}

void resolve_overlaps::wiggleFirstLinkages(GlycosylationSiteVector &glycosites, OverlapType overlapType, int max_cycles)
{
    std::cout << "Wiggling first linkages for " << max_cycles << " cycles.\n";
    int output_pdbfile_id = 0;
    double strict_tolerance = 0.1, loose_tolerance = 1.0;
    GlycosylationSitePointerVector sites_with_overlaps = DetermineSitesWithOverlap(glycosites, strict_tolerance, overlapType);
    int cycle = 0;
    bool stop = false;
    while ( (cycle < max_cycles) && (stop == false) )
    {
        std::cout << "Cycle " << cycle << "/" << max_cycles << "\n";
        ++cycle;
        std::random_shuffle (sites_with_overlaps.begin(), sites_with_overlaps.end());
        for(auto &glycosite : sites_with_overlaps)
        {
            glycosite->WiggleFirstLinkage(&output_pdbfile_id, strict_tolerance);
        }
        // Check which sites still have overlaps, stop if none.
        sites_with_overlaps = DetermineSitesWithOverlap(glycosites, strict_tolerance, overlapType); // Moved glycans may clash with other glycans. Need to check.
        if (sites_with_overlaps.size() == 0)
        {
            std::cout << "Stopping with all overlaps resolved.\n";
            stop = true;
        }
    }
    return;
}

//void resolve_overlaps::weighted_protein_global_overlap_monte_carlo(GlycosylationSiteVector &glycosites, int max_cycles)
//{
//    bool accept_change;
//    int cycle = 1;
//    bool stop = false;
//    double previous_chi1;
//    double previous_chi2;
//    double previous_glycan_overlap, new_glycan_overlap, previous_protein_overlap, new_protein_overlap;
//    double lowest_global_overlap = glycoprotein_builder::GetGlobalOverlap(glycosites);
//    double new_global_overlap;
//    double strict_tolerance = 0.1, loose_tolerance = 1.0;

//    GlycosylationSitePointerVector sites_with_overlaps = DetermineSitesWithOverlap(glycosites, strict_tolerance, "total");

//    std::cout << "Initial torsions and overlaps:\n";
//    glycoprotein_builder::PrintDihedralAnglesAndOverlapOfGlycosites(glycosites);
//    while ( (cycle < max_cycles) && (stop == false) )
//    {
//        ++cycle;
//        std::cout << "Cycle " << cycle << " of " << max_cycles << std::endl;
//        std::random_shuffle (sites_with_overlaps.begin(), sites_with_overlaps.end());
//        for(GlycosylationSitePointerVector::iterator it1 = sites_with_overlaps.begin(); it1 != sites_with_overlaps.end(); ++it1)
//        {
//            GlycosylationSite *current_glycosite = (*it1);
//           // std::cout << "Checking " << current_glycosite->GetResidue()->GetId() << "\n";
//            previous_glycan_overlap = current_glycosite->GetGlycanOverlap();
//            previous_protein_overlap = current_glycosite->GetProteinOverlap();
//            previous_chi1 = current_glycosite->GetChi1Value();
//            current_glycosite->SetChi1Value(RandomAngle_360range());
//            previous_chi2 = current_glycosite->GetChi2Value();
//            current_glycosite->SetChi2Value(RandomAngle_360range());
//            new_glycan_overlap = current_glycosite->Calculate_bead_overlaps_noRecord_noSet("glycan");
//            new_protein_overlap = current_glycosite->Calculate_bead_overlaps_noRecord_noSet("protein");


//            accept_change = monte_carlo::accept_via_metropolis_criterion((new_glycan_overlap + (new_protein_overlap*5)) - (previous_glycan_overlap + (previous_protein_overlap*5)));
//            if (!accept_change)
//            {
//                current_glycosite->SetChi1Value(previous_chi1);
//                current_glycosite->SetChi2Value(previous_chi2);
//            }
////            else
////            {
////              //  std::cout << "Accepted a change of " << ((new_overlap) - (previous_overlap)) << "\n";
////            }
//        }
//        //std::cout << "Updating list of sites with overlaps." << std::endl;
//        sites_with_overlaps = DetermineSitesWithOverlap(glycosites, strict_tolerance); // Moved glycans may clash with other glycans. Need to check.
//        if (sites_with_overlaps.size() == 0)
//        {
//            std::cout << "Stopping with all overlaps resolved.\n";
//            stop = true;
//        }

//        new_global_overlap = glycoprotein_builder::GetGlobalOverlap(glycosites);
////        std::cout << "Lowest: " << lowest_global_overlap << ", Current: " << new_global_overlap << "\n";
//        if ( lowest_global_overlap > new_global_overlap + 1 )
//        {
//            write_pdb_file(glycosites.at(0).GetGlycoprotein(), cycle, "best", new_global_overlap);
//            lowest_global_overlap = new_global_overlap;
//        }
//    }
////    std::cout << "Global overlap before deleting sites is " << glycoprotein_builder::GetGlobalOverlap(glycosites) << "\n";
////    std::cout << "Finished torsions and overlaps:\n";
////    glycoprotein_builder::PrintDihedralAnglesAndOverlapOfGlycosites(glycosites);
//    //DeleteSitesIterativelyWithOverlapAboveTolerance(glycosites, loose_tolerance);
//}

//void resolve_overlaps::Resolve_Overlaps_Random_Walk_Scaled_To_Overlap(GlycosylationSiteVector &glycosites, std::string type, int max_cycles, double strict_tolerance, double loose_tolerance)
//{
//    //double strict_tolerance = 0.1, loose_tolerance = 1.0; // Aim for <0.1 when resolving, but keep any less than 1 when culling.
//    GlycosylationSitePointerVector sites_with_overlaps = DetermineSitesWithOverlap(glycosites, strict_tolerance, type);
//    Overlap_Weighted_Adjust_Torsions_For_X_Cycles(sites_with_overlaps, glycosites, max_cycles, strict_tolerance, type);
//    DeleteSitesWithOverlapRecordsAboveTolerance(glycosites, loose_tolerance, type);
//    sites_with_overlaps = DetermineSitesWithOverlap(glycosites, strict_tolerance, type);
//    SetBestChi1Chi2(sites_with_overlaps, type);
//}

//void resolve_overlaps::protein_first_random_walk_scaled_to_overlap(GlycosylationSiteVector &glycosites)
//{
//    // NOTE!!! Upon testing, this algorithm fails with densely clustered glycans. It remembers good spots for individual glycans, but they can be in the same place.

//    /* Algorithm overview:
//     *  Resolve all protein overlap first, reject sites that cannot be resolved.
//     *  Calculate protein overlaps, for each site with overlaps
//     *        Change chi1 and chi2 values, record values that reduce overlaps
//     *  Once finished delete any sites that could not be resolved
//     *  Calculate total (protein+glycan) overlaps, for each site with overlaps
//     *        Change chi1 and chi2 values, record values that reduce overlaps
//     *  Once finished delete any sites that could not be resolved
//     */

////    std::cout << "Initial overlaps for all sites\n";
////    std::cout << "      Site        |  Total | Protein | Glycan \n";
////    for(GlycosylationSiteVector::iterator it1 = glycosites.begin(); it1 != glycosites.end(); ++it1)
////    {
////        GlycosylationSite current_glycosite = *it1;
////        current_glycosite.Calculate_and_print_bead_overlaps();
////    }

//     Resolve_Overlaps_Random_Walk_Scaled_To_Overlap(glycosites, "protein");

////    std::cout << "Overlaps for all sites after protein resolution and deletion\n";
////    std::cout << "      Site        |  Total | Protein | Glycan \n";
////    for(GlycosylationSiteVector::iterator it1 = glycosites.begin(); it1 != glycosites.end(); ++it1)
////    {
////        GlycosylationSite current_glycosite = *it1;
////        current_glycosite.Calculate_and_print_bead_overlaps();
////    }

//     Resolve_Overlaps_Random_Walk_Scaled_To_Overlap(glycosites, "total");
//    // glycoprotein_builder::PrintDihedralAnglesOfGlycosites(glycosites);
//    //PrintOverlaps(glycosites);
//}


bool resolve_overlaps::dumb_random_walk(GlycosylationSiteVector &glycosites, OverlapType overlapType)
{
    /* Algorithm:
     * Determine which sites have overlaps greater than tolerance. Stop if zero sites.
     * For each site with overlaps:
     *        Randomly change all chi1 and chi2 values
     */
    double tolerance = 0.1;
    int cycle = 0, max_cycles = 10;
    GlycosylationSitePointerVector sites_with_overlaps = DetermineSitesWithOverlap(glycosites, tolerance, overlapType);
    bool resolved = false;

    while ( (cycle < max_cycles) && (resolved == false) )
    {
        ++cycle;
        std::cout << "Cycle " << cycle << " of " << max_cycles << std::endl;
        for(GlycosylationSitePointerVector::iterator it1 = sites_with_overlaps.begin(); it1 != sites_with_overlaps.end(); ++it1)
        {
            GlycosylationSite *current_glycosite = (*it1);
            current_glycosite->SetRandomDihedralAnglesUsingMetadata();
        }
        //std::cout << "Updating list of sites with overlaps." << std::endl;
        sites_with_overlaps = DetermineSitesWithOverlap(glycosites, tolerance, overlapType); // Moved glycans may clash with other glycans. Need to check.
        if (sites_with_overlaps.size() == 0)
        {
            std::cout << "Stopping with all overlaps resolved.\n";
            resolved = true;
        }
    }
    return resolved;
}


