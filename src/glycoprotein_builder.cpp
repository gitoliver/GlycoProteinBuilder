#include "../includes/glycoprotein_builder.h"

/*******************************************/
/* Functions                               */
/*******************************************/

void glycoprotein_builder::AttachGlycansToGlycosites(MolecularModeling::Assembly &glycoprotein, GlycosylationSiteVector &glycoSites, const std::string glycanDirectory)
{
    // Find protein residues in Glycoprotein that will get a glycan added. Set Residue in Glycosite.
    // Often there are multiple chains with the same residue number. If the input file repeats a residue number, it will go on the next instance (next chain)
    // of that residue number. This is not a good long term solution, as users will want to select chains.
    ResidueVector protein_residues = glycoprotein.GetResidues();
    for(GlycosylationSiteVector::iterator glycosite = glycoSites.begin(); glycosite != glycoSites.end(); ++glycosite)
    {
        bool stop = false;
        std::string glycosite_number = glycosite->GetResidueNumber();
        ResidueVector::iterator it2 = protein_residues.begin();
        while ( (!stop) && (it2 != protein_residues.end()) ) {
           // for (ResidueVector::iterator it2 = protein_residues.begin(); it2 != protein_residues.end(); ++it2)
            {
                MolecularModeling::Residue *protein_residue = *it2;
                std::string id = protein_residue->GetId();
                std::string formatted_glycosite_number = "_" + glycosite_number + "_";
                if( id.compare(5, formatted_glycosite_number.size(), formatted_glycosite_number) == 0)
                {
                    //std::cout << "glycosite: " << glycosite_number << std::endl;
                    std::cout << "glycosite id:" << id << std::endl;
                    glycosite->SetResidue(protein_residue);
                    stop = true;
                    // Remove residue from list, so if user has listed same residue name twice, it will go on next instance (i.e. on next chain) of residue number
                    protein_residues.erase(std::remove(protein_residues.begin(), protein_residues.end(), *it2), protein_residues.end()); // Note need #include <algorithm>
                }
                ++it2;
            }
        }
    }
//    // Find protein residues in Glycoprotein that will get a glycan added. Set Residue in Glycosite.
//    ResidueVector protein_residues = glycoprotein.GetResidues();
//    for (ResidueVector::iterator it2 = protein_residues.begin(); it2 != protein_residues.end(); ++it2)
//    {
//        Residue *protein_residue = *it2;
//        std::string id = protein_residue->GetId();
//        for(GlycosylationSiteVector::iterator glycosite = glycoSites.begin(); glycosite != glycoSites.end(); ++glycosite)
//        {
//            std::string glycosite_number = glycosite->GetResidueNumber();
//            std::string formatted_glycosite_number = "_" + glycosite_number + "_";
//            if( id.compare(5, formatted_glycosite_number.size(), formatted_glycosite_number) == 0)
//            {
//                //std::cout << "glycosite: " << glycosite_number << std::endl;
//                std::cout << "glycosite id:" << id << std::endl;
//                glycosite->SetResidue(protein_residue);
//            }
//        }
//    }
    // Load glycan files from directory
    //std::cout << "Glycan directory: " << glycanDirectory << std::endl;
    std::string filepath;
    DIR *dp; // A directory stream
    struct dirent *dirp; // Contains file serial number and name (char d_name[])
    struct stat filestat; // Contains info about file, such as device ID, user ID, access time etc

    dp = opendir( glycanDirectory.c_str() ); //.c_str adds a null character to the end.
    if (dp == NULL)
    {
        std::cout << "Error(" << errno << ") opening " << glycanDirectory << std::endl;
        return;
    }
    while ((dirp = readdir ( dp )))
    {
        filepath = glycanDirectory + "/" + dirp->d_name;
        // If the file is a directory (or is in some way invalid) we'll skip it
        if (stat( filepath.c_str(), &filestat )) continue; // Is it a valid file?
        if (S_ISDIR( filestat.st_mode ))         continue; // Is it a directory?
        for (GlycosylationSiteVector::iterator glycosite = glycoSites.begin(); glycosite != glycoSites.end(); ++glycosite)
        {
            //std::cout << "Glycan is " << glycosite->GetGlycanName() << ". d_name is " << dirp->d_name << std::endl;
            if (glycosite->GetGlycanName().compare(0, glycosite->GetGlycanName().size(), dirp->d_name, 0, glycosite->GetGlycanName().size()) == 0 )
            {
                MolecularModeling::Assembly input_glycan(filepath, gmml::InputFileType::PDB);
                input_glycan.BuildStructureByDistance();
                glycosite->AttachGlycan(input_glycan, glycoprotein);
                //std::cout << "Added " << glycosite->GetGlycanName() << " to " << glycosite->GetResidueNumber();
            }
        }
    }
    closedir( dp );
}

void glycoprotein_builder::Read_Input_File(GlycosylationSiteVector &glycoSites, std::string &proteinPDB, std::string &glycanDirectory, const std::string working_Directory)
{
    std::string buffer;
    std::ifstream inf (working_Directory + "/input.txt");
    if (!inf)
    {
        std::cerr << "Uh oh, input file could not be opened for reading!" << std::endl;
        std::exit(1);
    }
    while (inf) // While there's still stuff left to read
    {
        std::string strInput;
        getline(inf, strInput);
        if(strInput == "Protein:")
            getline(inf, proteinPDB);
        if(strInput == "Glycans:")
        {
            getline(inf, glycanDirectory);
            glycanDirectory = working_Directory + "/" + glycanDirectory;
        }
        if(strInput == "Protein Residue, Glycan Name:")
        {
            getline(inf, buffer);
            while(buffer != "END")
            {
                StringVector splitLine = split(buffer, ',');
                glycoSites.emplace_back(splitLine.at(1), splitLine.at(0)); // Creates GlycosylationSite instance on the vector. Love it.
                getline(inf, buffer);
            }
        }
    }
}
