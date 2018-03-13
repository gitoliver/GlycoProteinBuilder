#include "../includes/resolve_overlaps.h"
#include "../includes/glycosylationsite.h"


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
    return start_point + (rand() % (max_step_size * 2) + 1) - max_step_size;
}

double GetNewAngleScaledToDegreeOfOverlap(double current_angle, double overlap, int number_of_atoms)
{
    int max_step_size = 1 + std::round( 180 * ( overlap / number_of_atoms ) ); // Always allow at least 1 degrees of movement
    return RandomAngle_PlusMinusX(current_angle, max_step_size);
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

void PrintOverlaps(GlycosylationSiteVector *glycosites)
{
    for (GlycosylationSiteVector::iterator current_glycosite = glycosites->begin(); current_glycosite != glycosites->end(); ++current_glycosite)
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

void SetBestProteinChi1Chi2(GlycosylationSitePointerVector &glycosites, Assembly *glycoprotein)
{
    for (GlycosylationSitePointerVector::iterator it = glycosites.begin(); it != glycosites.end(); ++it)
    {
        GlycosylationSite *current_glycosite = (*it);
        current_glycosite->SetChi1Value(current_glycosite->GetBestProteinOverlapRecord().GetChi1(), glycoprotein);
        current_glycosite->SetChi2Value(current_glycosite->GetBestProteinOverlapRecord().GetChi2(), glycoprotein);
    }
}

GlycosylationSitePointerVector DetermineSitesWithOverlap(GlycosylationSiteVector *glycosites, double tolerance, std::string type = "total")
{
    GlycosylationSitePointerVector sites_with_overlaps;
    double overlap = 0.0;
    std::cout << "      Site        |  Total | Protein | Glycan " << std::endl;
    for (GlycosylationSiteVector::iterator current_glycosite = glycosites->begin(); current_glycosite != glycosites->end(); ++current_glycosite)
    {
        if(type.compare("total")==0)
        {
            overlap = current_glycosite->Calculate_bead_overlaps();
        }
        else if (type.compare("protein")==0)
        {
            overlap = current_glycosite->Calculate_protein_bead_overlaps();
        }
        else if (type.compare("glycan")==0)
        {
            overlap = current_glycosite->Calculate_other_glycan_bead_overlaps();
        }
        if ( overlap > tolerance)
        {
            current_glycosite->Print_bead_overlaps();
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

    c) Create a tree structure for sites that overlap with each other.
        Could use assembly index/id of bead atom? Can work on individual trees and report nice info to users.


Note: Chi1 is 180,-60,60 +/- 30 degrees. Want a function that keeps the values within these states.
Note: Need to save best structure.
Note: Remember to delete the SUP atoms at some stage.
*/

void resolve_overlaps::monte_carlo(Assembly *glycoprotein, GlycosylationSiteVector *glycosites)
{
    /* Algorithm:
     *  Resolve all protein overlap first, reject sites that cannot be resolved. Record cycles to resolve?
     *  Calculate protein overlaps, for each site with overlaps
     *        Randomly change all chi1 and chi2 values
     *  Once finished;
     *  Calculate glycan overlaps, for each site with overlaps
     *        Change chi1 and chi2, scaled by degree of overlap
     *        Reject changes causing protein overlaps
     *
     */
    double tolerance = 0.1;
    double new_dihedral_angle = 0.0;
    int cycle = 0, max_cycles = 5;

    while (cycle < max_cycles)
    {
        ++cycle;
        for(GlycosylationSiteVector::iterator current_glycosite = glycosites->begin(); current_glycosite != glycosites->end(); ++current_glycosite)
        {
            current_glycosite->Calculate_protein_bead_overlaps();
            new_dihedral_angle = GetNewAngleScaledToDegreeOfOverlap(current_glycosite->GetChi1Value(), current_glycosite->GetProteinOverlap(), current_glycosite->GetAttachedGlycan()->GetAllAtomsOfAssembly().size());
            current_glycosite->SetChi1Value(new_dihedral_angle, glycoprotein);
            new_dihedral_angle = GetNewAngleScaledToDegreeOfOverlap(current_glycosite->GetChi2Value(), current_glycosite->GetProteinOverlap(), current_glycosite->GetAttachedGlycan()->GetAllAtomsOfAssembly().size());
            current_glycosite->SetChi2Value(new_dihedral_angle, glycoprotein);
            current_glycosite->Calculate_protein_bead_overlaps();

        }
    }

//    GlycosylationSitePointerVector sites_with_protein_overlaps = DetermineSitesWithOverlap(glycosites, tolerance, "protein");
//    bool stop = false;

//    while ( (cycle < max_cycles) && (stop == false) )
//    {
//        ++cycle;
//        std::cout << "Cycle " << cycle << " of " << max_cycles << std::endl;
//        for(GlycosylationSitePointerVector::iterator it1 = sites_with_protein_overlaps.begin(); it1 != sites_with_protein_overlaps.end(); ++it1)
//        {
//            GlycosylationSite *current_glycosite = (*it1);
//            current_glycosite->Calculate_protein_bead_overlaps();

//            new_dihedral_angle = GetNewAngleScaledToDegreeOfOverlap(current_glycosite->GetChi1Value(), current_glycosite->GetProteinOverlap(), current_glycosite->GetAttachedGlycan()->GetAllAtomsOfAssembly().size());
//           // std::cout << "Changing chi1 from " << current_glycosite->GetChi1Value() << " to " << new_dihedral_angle << " based on protein overlap of " << current_glycosite->GetProteinOverlap() << std::endl;
//            current_glycosite->SetChi1Value(new_dihedral_angle, glycoprotein);
//          //  std::cout << "New value: " << current_glycosite->GetChi1Value() << std::endl;
//            new_dihedral_angle = GetNewAngleScaledToDegreeOfOverlap(current_glycosite->GetChi2Value(), current_glycosite->GetProteinOverlap(), current_glycosite->GetAttachedGlycan()->GetAllAtomsOfAssembly().size());
//         //   std::cout << "Changing chi2 from " << current_glycosite->GetChi2Value() << " to " << new_dihedral_angle << " based on protein overlap of " << current_glycosite->GetProteinOverlap() << std::endl;
//            current_glycosite->SetChi2Value(new_dihedral_angle, glycoprotein);
//      //      std::cout << "New value: " << current_glycosite->GetChi2Value() << std::endl;

////            current_glycosite->Calculate_protein_bead_overlaps();
////            current_glycosite->Print_bead_overlaps();
////            current_glycosite->SetChi1Value(-73.4055522836, glycoprotein);
////            current_glycosite->SetChi2Value(-132.8354561298, glycoprotein);
////            current_glycosite->Calculate_protein_bead_overlaps();
////            current_glycosite->Print_bead_overlaps();
////            current_glycosite->SetChi1Value(-73.4055522836, glycoprotein);
////            current_glycosite->SetChi2Value(-132.8354561298, glycoprotein);
////            current_glycosite->Calculate_protein_bead_overlaps();
////            current_glycosite->Print_bead_overlaps();
////            current_glycosite->SetChi1Value(0.00, glycoprotein);
////            current_glycosite->SetChi2Value(0.00, glycoprotein);
////            current_glycosite->Calculate_protein_bead_overlaps();
////            current_glycosite->Print_bead_overlaps();
////            current_glycosite->SetChi1Value(-73.4055522836, glycoprotein);
////            current_glycosite->SetChi2Value(-132.8354561298, glycoprotein);
//               current_glycosite->Calculate_protein_bead_overlaps();
//         //   current_glycosite->Print_bead_overlaps();

//            //double percent_overlap = ((current_glycosite->GetTotalOverlap() / (current_glycosite->GetAttachedGlycan()->GetAllAtomsOfAssembly().size()) ) + 0.01);
//            //  new_dihedral_value = RandomAngle_PlusMinusX(current_glycosite->GetChi1Value(), (180 * percent_overlap) ); // scaled to degree of overlap
//            //        sites_with_overlaps.erase(std::remove(sites_with_overlaps.begin(), sites_with_overlaps.end(), *it1), sites_with_overlaps.end());
//        }
//        //std::cout << "Updating list of sites with overlaps." << std::endl;
//       // sites_with_protein_overlaps = DetermineSitesWithOverlap(glycosites, tolerance, "protein"); // Moved glycans may clash with other glycans. Need to check.
//        if (sites_with_protein_overlaps.size() == 0)
//        {
//            std::cout << "No more protein overlaps" << std::endl;
//            stop = true;
//        }
//    }
    //PrintOverlaps(glycosites);
  //  write_pdb_file(glycoprotein, 1, "./outputs/summary", 0.0);
    // Remove sites that cannot be resolved.
 //   std::cout << "Setting best chi1 and chi2 found so far" << std::endl;
  //  SetBestProteinChi1Chi2(sites_with_protein_overlaps, glycoprotein);
 //   std::cout << "Could not resolve protein overlaps for these sites: " << std::endl;
   // sites_with_protein_overlaps = DetermineSitesWithOverlap(glycosites, tolerance, "protein");
    //PrintOverlaps(sites_with_protein_overlaps);
}

void resolve_overlaps::dumb_monte_carlo(Assembly *glycoprotein, GlycosylationSiteVector *glycosites)
{
    /* Algorithm:
     * Determine which sites have overlaps greater than tolerance. Stop if zero sites.
     * For each site with overlaps:
     *        Randomly change all chi1 and chi2 values
     */
    double tolerance = 0.1;
    int cycle = 0, max_cycles = 2;
    GlycosylationSitePointerVector sites_with_overlaps = DetermineSitesWithOverlap(glycosites, tolerance);
    bool stop = false;

    while ( (cycle < max_cycles) && (stop == false) )
    {
        ++cycle;
        std::cout << "Cycle " << cycle << " of " << max_cycles << std::endl;
        for(GlycosylationSitePointerVector::iterator it1 = sites_with_overlaps.begin(); it1 != sites_with_overlaps.end(); ++it1)
        {
            GlycosylationSite *current_glycosite = (*it1);

            current_glycosite->SetChi1Value(RandomAngle_360range(), glycoprotein);

           // current_glycosite->SetChi2Value(RandomAngle_360range(), glycoprotein);
            //double percent_overlap = ((current_glycosite->GetTotalOverlap() / (current_glycosite->GetAttachedGlycan()->GetAllAtomsOfAssembly().size()) ) + 0.01);
            //  new_dihedral_value = RandomAngle_PlusMinusX(current_glycosite->GetChi1Value(), (180 * percent_overlap) ); // scaled to degree of overlap
            //        sites_with_overlaps.erase(std::remove(sites_with_overlaps.begin(), sites_with_overlaps.end(), *it1), sites_with_overlaps.end());
        }
        write_pdb_file(glycoprotein, cycle, "./outputs/summary", 0.0);
        //std::cout << "Updating list of sites with overlaps." << std::endl;
        sites_with_overlaps = DetermineSitesWithOverlap(glycosites, tolerance); // Moved glycans may clash with other glycans. Need to check.
        if (sites_with_overlaps.size() == 0)
        {
            std::cout << "Stopping with all overlaps resolved.\n";
            stop = true;
        }
    }
    write_pdb_file(glycoprotein, (cycle + 1), "./outputs/summary", 0.0);
}



void resolve_overlaps::example_for_Gordon(Assembly *glycoprotein, GlycosylationSiteVector *glycosites)
{
    double site_total_overlap = 0.0, site_glycan_overlap = 0.0, site_protein_overlap = 0.0, new_dihedral_value = 0.0, total_system_overlap = 0.0;
    std::cout << "      Site        |  Total | Protein | Glycan " << std::endl;

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
