#ifndef RESOLVE_OVERLAPS_H
#define RESOLVE_OVERLAPS_H

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

//int main ( );
void crossover ( int &seed );
void elitist ( );
void evaluate ( Assembly *glycoprotein, GlycosylationSiteVector *glycosites);
int i4_uniform_ab ( int a, int b, int &seed );
void initialize ( std::string filename, int &seed );
void keep_the_best ( );
void mutate ( int &seed );
double r8_uniform_ab ( double a, double b, int &seed );
void report ( int generation );
void selector ( int &seed );
void timestamp ( );
void Xover ( int one, int two, int &seed );

namespace resolve_overlaps
{
    void dumb_monte_carlo(MolecularModeling::Assembly *glycoprotein, GlycosylationSiteVector *glycosites);
    void example_for_Gordon(MolecularModeling::Assembly *glycoprotein, GlycosylationSiteVector *glycosites);
    void genetic_algorithm(MolecularModeling::Assembly *glycoprotein, GlycosylationSiteVector *glycosites);
}

#endif // RESOLVE_OVERLAPS_H
