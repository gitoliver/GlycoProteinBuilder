#########################
###   FUNCTIONS        ##
#########################

void print_coords(Atom* atom)
		Simply used to print the xyz coordinates of an atom

void Implant_FatAtoms(Assembly glycoprotein, GlycoSiteVector glycosites)
		Adds fat atoms to the assembly objects
		Goes through every residue in assembly with CheckIfProtein(); if the residue is an amino acid then place a fat atom at CA; if not (assumed to be glycan) then place fat atom at ring center (GetRingCenter() currently is not working)
		All fat atoms are named like "3fatom" where the 3 is assigned radius

void Sacrifice_FatAtoms(Assembly glycoprotein)
		Goes through every residue in the assembly and delete any that have the word "fatom" in its name

AtomVector Filter_FatAtoms(AtomVector atoms)
		Goes through every residue in the assembly and make a new atom vector containing any residue with "fatom" in its name

double Atomwise_CalculateAtomicOverlaps(Atom *atomA, Atom *atomB, double radiusA, double radiusB)
		Mostly copied from GMML's CalculateAtomicOverlaps()
		Intended to only be used when scoring fat atoms; tweak as needed	

double modified_CalculateAtomicOverlaps(AtomVector atomsA, AtomVector atomsB)
		Mostly copied from GMML's CalculateAtomicOverlaps()
		Intended to only be used when scoring fat atoms; tweak as needed

ResidueVector ResiFilter_ScoreFatAtomOverlap(Assembly* glycoprotein, GlycoSiteVector* glycosites, double* overlap_score, double threshold)
		Scores a glycoprotein assembly object by each glycan chain using FAT ATOMS; provides a rough estimate and is pretty fast
		Pass in a threshold overlap score, all glycan chains that have a higher clash score than the set threshold will be output as a ResidueVector of the amino acid residue that holds the glycan chain (in other words, the output ResidueVector will hold residue connected to glycan chains that scored higher than the set threshold)
		The glycans attached to the output ResidueVector will be move by a seperate function (ResiRotor)
		Writes the aggregate ovelap score all of glycan chains to the a double

ResidueVector ResiFilter_ScoreTrueOverlap(Assembly* glycoprotein, GlycoSiteVector* glycosites, double* overlap_score, double threshold)
		Scores a glycoprotein assembly object by each glycan chain using all atoms; provides a more exact score is ia horribly slow
		Pass in a threshold overlap score, all glycan chains that have a higher clash score than the set threshold will be output as a ResidueVector of the amino acid residue that holds the glycan chain (in other words, the output ResidueVector will hold residue connected to glycan chains that scored higher than the set threshold)
		The glycans attached to the output ResidueVector will be move by a seperate function (ResiRotor)
		Writes the aggregate ovelap score all of glycan chains to the a double

ResidueVector ResiFilter_Aggregate(Assembly* glycoprotein, GlycoSiteVector* glycosites)
		Outputs all residues that have a glycan chain attached to them
		Intended to be used with ResiRotor functions to shuffle all glycan chains in the whole glycoprotein

double RandomAngle_360range()
		Generates a random number between 1 and 360

double RandomAngle_PlusMinusX(double start_point, int max_step_size)
		Generates a random number within a given step size from the given start point
		Example: RandomAngle_PlusMinusX(10, 3) would yield numbers between 7-13

void ResiRotor_FullRange(Assembly* glycoprotein, ResidueVector* move_these_guys)
		Rotates the dihedrals of the residues in the given ResidueVector
		Rotation is done with RandomAngle_360range() so full range rotation is allowed
		Will move ASN(X1 & X2), TYR(X1 & X2), THR(X1), SER(X1)
		== Alternative naming is also supported for glycosylated residues (NLN, OLY, OLT, OLS)

void ResiRotor_BabyStep(Assembly* glycoprotein, ResidueVector* move_these_guys)
		Rotates the dihedrals of the residues in the given ResidueVector
		Rotation is done with RandomAngle_PlusMinusX() so move range is limited
		Will move ASN(X1 & X2), TYR(X1 & X2), THR(X1), SER(X1)
		== Alternative naming is also supported for glycosylated residues (NLN, OLY, OLT, OLS)

void write_pdb_file(Assembly glycoprotein, int cycle, string summary_filename, double score)
		Outputs a pdb file and adds an entry into the log file

#########################
###   MAIN FUNCTION    ##
#########################

The input glycoprotein is immediately scored for overlaps using fat atoms.

The lowest/ minimum score found will recorded each cycle.

Whenever an overlap score (done with fat atoms) reaches a minimum, it gets scored atom-wise for its true overlap score. Then, the conformation will be output as a pdb file.

Glycan chains found to be above the threshold during the score function (ResiFilter) will be moved with ResiRotor

The number of cycles since the last minima found is also recorded, after a certain number of cycles without progress, ALL glycans chains will be shuffled.




