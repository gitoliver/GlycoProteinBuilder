#ifndef GENETIC_ALGORITHM_H
#define GENETIC_ALGORITHM_H

#include <iostream>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <string>
#include "bead_residues.h"

# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

typedef std::vector<GlycosylationSite*> GlycosylationSitePointerVector;

//int main ( );
void crossover ( Assembly *glycoprotein, GlycosylationSiteVector *glycosites);
void elitist ( );
void evaluate ( Assembly *glycoprotein, GlycosylationSiteVector *glycosites);
int i4_uniform_ab ( int a, int b, int &seed );
void initialize ( std::string filename, int &seed );
void keep_the_best ( );
void mutate (Assembly *glycoprotein, GlycosylationSiteVector *glycosites );
double r8_uniform_ab ( double a, double b, int &seed );
void report (  int generation, Assembly *glycoprotein, GlycosylationSiteVector *glycosites );
void selector ( int &seed );
void timestamp ( );
void Xover ( int one, int two, int &seed, Assembly *glycoprotein, GlycosylationSiteVector *glycosites );

namespace resolve_overlaps
{
    void genetic_algorithm(MolecularModeling::Assembly *glycoprotein, GlycosylationSiteVector *glycosites);
}

#endif // GENETIC_ALGORITHM_H
