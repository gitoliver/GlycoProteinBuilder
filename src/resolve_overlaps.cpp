

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
GlycosylationSitePointerVector DetermineSitesWithOverlap(GlycosylationSiteVector *glycosites)
{
    GlycosylationSitePointerVector sites_with_overlaps;
    for (GlycosylationSiteVector::iterator current_glycosite = glycosites->begin(); current_glycosite != glycosites->end(); ++current_glycosite)
    {
        if (current_glycosite->Calculate_bead_overlaps() > tolerance)
        {
            sites_with_overlaps.push_back(&(*current_glycosite));
        }
    }
    return sites_with_overlaps;
}

void resolve_overlaps::monte_carlo(Assembly *glycoprotein, GlycosylationSiteVector *glycosites)
{
    double new_dihedral_value = 0.0, worst_site_overlap = 0.0, tolerance = 0.1;
    int cycle = 0, max_cycles = 1000;
    GlycosylationSitePointerVector sites_with_overlaps =  DetermineSitesWithOverlap(glycosites);

    bool stop = false;
    while ( (cycle < max_cycles) && (stop == false) )
    {
        std::cout << "Cycle " << cycle << " of " << max_cycles << std::endl;
        std::cout << "      Site        | Total | Protein | Glycan " << std::endl;
        worst_site_overlap = 0.0; // reset.
        for(GlycosylationSitePointerVector::iterator it1 = sites_with_overlaps.begin(); it1 != sites_with_overlaps.end(); /* Nothing here as need to erase() from vector */)
        {
            GlycosylationSite *current_glycosite = (*it1);
            //double percent_overlap = ((current_glycosite->GetTotalOverlap() / (current_glycosite->GetAttachedGlycan()->GetAllAtomsOfAssembly().size()) ) + 0.01);
          //  new_dihedral_value = RandomAngle_PlusMinusX(current_glycosite->GetChi1Value(), (180 * percent_overlap) ); // scaled to degree of overlap
            new_dihedral_value = RandomAngle_360range();
           // std::cout << "Changing chi1 from " << current_glycosite->GetChi1Value() << " to " << new_dihedral_value << std::endl;
            current_glycosite->SetChi1Value(new_dihedral_value, glycoprotein);
            //new_dihedral_value = RandomAngle_PlusMinusX(current_glycosite->GetChi2Value(), (180 * percent_overlap));
            new_dihedral_value = RandomAngle_360range();
            current_glycosite->SetChi2Value(new_dihedral_value, glycoprotein);
            current_glycosite->Calculate_and_print_bead_overlaps();
            if (current_glycosite->GetTotalOverlap() > worst_site_overlap)
            {
                worst_site_overlap = current_glycosite->GetTotalOverlap();
            }
            if (current_glycosite->GetTotalOverlap() <= tolerance) // Remove site fromlist of sites with overlaps
            {
              //  std::cout << "segfaulting?" << std::endl;
                sites_with_overlaps.erase(std::remove(sites_with_overlaps.begin(), sites_with_overlaps.end(), *it1), sites_with_overlaps.end());
            }
            else
            {
                ++it1;
            }
        }
      //  std::cout << "Worst overlapping site has overlap of " << worst_site_overlap << std::endl;
        if (worst_site_overlap <= tolerance)
        {
            stop = true;
            std::cout << "We are like, totally STOPPING?" << std::endl;
        }
        ++cycle;
    }
    write_pdb_file(glycoprotein, 1, "./outputs/summary", worst_site_overlap);
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
