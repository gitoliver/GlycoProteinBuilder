# GlycoProteinBuilder
Uses GEMS/GMML to add and adapt 3D structures of N-glycans and O-glycans onto glycoproteins.
Protein residues currently handled: Asn, Ser, Thr and Tyr

# Compilation in a bash shell
export GEMSHOME=<Your Path To Gems > # eg: export GEMSHOME=/home/oliver/Programs/gems
qmake
make

# Setup
Edit or create an input.txt file and place in a folder called inputs. See inputs/input.txt for an example.
If running outside of the program directory, create a directory called outputs/
Providing Parameters is optional, everything else is required.
In Glycan id list, list

You must provide:
    a protein 3D structure
    glycan 3D structure(s)
    input.txt, which contains:
        parameter files (optional)
        protein file name
        name of the folder containing the glycans e.g. glycans
        the protein residue numbers you want to attach to (no automatic detection of sequons)
        the name of the glycan you want to attach. Be careful that the name matches the start of the name of the glycan file.
        e.g. m9 will match to m9* in the glycan folder, while m9-gtgt will cause just that rotamer to be added.
        Note that for each protein residue number provided a glycan must be detailed. Allowing you to add different glycans to different sites

