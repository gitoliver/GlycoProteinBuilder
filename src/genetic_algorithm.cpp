// genetic algorithm 

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <cstring>
#include <ostream>

#include <string>
#include <vector> 
#include <stdio.h>
#include <stdlib.h>
#include <sstream> 
#include <math.h>
#include <ctime>
#include <cstring>
#include <time.h>

#include "../includes/genetic_algorithm.h"
#include "../includes/resolve_overlaps.h"

#define POPSIZE 100
#define MAXGENS 100
#define NVARS 52
#define PXOVER 0.9
#define PMUTATION 0.4
#define parameter 60
#define minimum_fitness 1.0
#define angle_change 5.0
#define elite_factor .10
#define NUMBER_CROSSOVERS 1

double test=0;
double value=0;

double site_total_overlap = 0.0; 
double site_glycan_overlap = 0.0;  
double site_protein_overlap = 0.0; 
double new_dihedral_value = 0.0; 
double total_system_overlap = 0.0;

struct genotype
{
  double gene[NVARS];
  double fitness;
  double upper[NVARS];
  double lower[NVARS];
  double rfitness;
  double cfitness;
};


struct genotype population[POPSIZE];
struct genotype newpopulation[POPSIZE]; 
struct genotype population_reordered[POPSIZE];
struct genotype previous_population[(int) (POPSIZE*elite_factor)];
struct genotype previous_best[1];

struct genotype population2[POPSIZE];
struct genotype temp_population_reordered;


time_t t;
time_t timer; 
time_t timezero=time(0);

using namespace gmml;

//
//  Each GENOTYPE is a member of the population, with
//  gene: a string of variables,
//  fitness: the fitness
//  upper: the variable upper bounds,
//  lower: the variable lower bounds,
//  rfitness: the relative fitness,
//  cfitness: the cumulative fitness.
//

void resolve_overlaps::genetic_algorithm(Assembly *glycoprotein, GlycosylationSiteVector *glycosites)
{


previous_best[0].fitness=-1000000;

  double site_total_overlap = 0.0; 
  double site_glycan_overlap = 0.0;  
  double site_protein_overlap = 0.0; 
  double new_dihedral_value = 0.0; 
  double total_system_overlap = 0.0;
  
  std::cout << "      Site        |  Total | Protein | Glycan " << std::endl;

   /*
    for (GlycosylationSiteVector::iterator current_glycosite = glycosites->begin(); current_glycosite != glycosites->end(); ++current_glycosite)
    {
        site_total_overlap = current_glycosite->Calculate_bead_overlaps(); // Must repeat after rotating chi1, chi2.
        site_glycan_overlap = current_glycosite->GetGlycanOverlap(); // If you wish to have this level of detail
        site_protein_overlap = current_glycosite->GetProteinOverlap(); // If you wish to have this level of detail
        current_glycosite->Print_bead_overlaps();
        new_dihedral_value = RandomAngle_PlusMinusX(current_glycosite->GetChi1Value(), 6);
        //std::cout << "Changing chi1 from " << current_glycosite->GetChi1Value() << " to " << new_dihedral_value << std::endl;
        current_glycosite->SetChi1Value(new_dihedral_value);
        new_dihedral_value = RandomAngle_PlusMinusX(current_glycosite->GetChi2Value(), 6);
        current_glycosite->SetChi2Value(new_dihedral_value);
        current_glycosite->Calculate_and_print_bead_overlaps();
        total_system_overlap += site_total_overlap;
    }
   */

  srand(time(NULL));

  std::string filename = "simple_ga_input.txt";
  int generation;
  int i;
  int seed;

  timestamp ( );

  if ( NVARS < 2 )
  {
    cout << "\n";
    cout << "  The crossover modification will not be available,\n";
    cout << "  since it requires 2 <= NVARS.\n";
  }

  seed = 123456789;

  initialize ( filename, seed );

  cout << "population size " << POPSIZE << endl;
  cout << "max iterations " << MAXGENS << endl;

  for ( generation = 0; generation < MAXGENS; generation++ ){
    cout << endl << "generation " << generation << endl << endl; 


    cout << endl;
    
    //    cout << "population fitness before operations" << endl;

    //    for(int i=0; i<POPSIZE; i++){
    //       cout << population[i].fitness << endl;
    //    }


    //    cout << "selector " << endl;
    //    selector ( seed );

    cout << "crossover " << endl;
    crossover (glycoprotein,glycosites);

    cout << "mutate " << endl;
    mutate (glycoprotein, glycosites );
    
    cout << "evaluate " << endl << endl;
    evaluate(glycoprotein,glycosites);

    cout << "keep the best" << endl; 
    keep_the_best();
	
    cout << "elitist" << endl;
    elitist();

    cout << "report " << endl;
    report(generation,glycoprotein,glycosites);

    //    cout << population[POPSIZE-1].fitness << endl;

  }

  ofstream geneticoutput;

  geneticoutput.open("genetic_algorithm_output.txt",ios::out|ios::app);
  geneticoutput << endl;  
  geneticoutput << "  Best member after " << MAXGENS << " generations:\n";

  for(int i = 0; i < NVARS; i++ )
      {
      cout << "best :  var(" << i << ") = " << population[POPSIZE-1].gene[i] << "\n";
      geneticoutput <<  "  var(" << i << ") = " << population[POPSIZE-1].gene[i] << "\n";
  }

  cout << "\n";
  cout << "  Best fitness = " << population[POPSIZE-1].fitness << "\n";

  geneticoutput << "  Best fitness = " << population[POPSIZE-1].fitness << "\n";

  cout << "\n";
  cout << "SIMPLE_GA:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  geneticoutput << "\n  Normal end of execution.\n";
  geneticoutput.close();

//**********

// output of pdb_file

//    write_pdb_file(glycoprotein, 1, "./outputs/summary", total_system_overlap);
}


//******************************************************************

void crossover (Assembly *glycoprotein, GlycosylationSiteVector *glycosites){

// 
//  Purpose:
//
//    CROSSOVER selects two parents for the single point crossover.
//

  const double a = 0.0;
  const double b = 1.0;
  int mem=0;
  int first = 0;
  double x;
  int dont_crossover=0;
  int one=0;
  int two=0;
  int member=0;
  int total=0;
  double total_system_overlap=0;
  double site_total_overlap=0;
	
  int i;
  int point;
  double t;
  
// printout of before crossover
  
//  cout << endl;
//  cout << "population fitness before crossover" << endl << endl;
//  for(int i=0; i<POPSIZE; i++){
//     cout << population[i].fitness << endl;
//  }
//  cout << endl;

  
// multiple single point crossovers set by NUMBER_CROSSOVERS
  
for(int index_crossover=0; index_crossover<NUMBER_CROSSOVERS; index_crossover++){

// for each chromosome 

  for(int mem = 0; mem < POPSIZE*(1-elite_factor); ++mem ){

// set default

      dont_crossover=0;

	  
// just a test for crossover eliminated from elitist part of population 
	  
     for(int test=0; test<POPSIZE*elite_factor; test++){

// printout of elite part of population 

//       cout << previous_population[test].fitness << endl;

        if(population[mem].fitness==previous_population[test].fitness){ 
          dont_crossover=1; 
        }
     } // end of test
     cout << endl;

	 
// dont_crossover=0;

// if dont_crossover=0, proceed with crossover	 
	 
	 one=mem;
	 
     if(dont_crossover==0){

	    x=(double) rand()/RAND_MAX;

        if(x < PXOVER ){  // crossover
		
// the second chromosome is a random choice for crossover
			
	  two= (round) ((double) (POPSIZE-1)*rand()/RAND_MAX); 
  
//
//  Select the crossover point.
//

            point = (round) ((double) NVARS/2 * rand()/RAND_MAX);

//
//  Swap genes in positions 0 through POINT-1.
//

  double test=(double) rand()/RAND_MAX;

  //  cout << test << " " << point << endl;

  
// swap 2 different ways depending upon probability
  
    if(test<=.5){
      for(i = 0; i < 2*point; i++){
	     t= population[one].gene[i];
         population[one].gene[i] = population[two].gene[i];
         population[two].gene[i] = t;
      }
    }

    if(test>.5){
      for(i=2*point; i<NVARS; i++){
	     t=population[one].gene[i];
         population[one].gene[i]=population[two].gene[i];
         population[two].gene[i]=t;
      }
    } 
	
	
	  } // end x<PXOVER
   } // end dont_crossover check
  
 
  } // member 
} // index crossover
	
	
// print out the population fitnesses of the crossovered population
// calculate the new fitnesses of the changed population
	
      total_system_overlap=0;
	  
	 for(int member=0; member<POPSIZE; member++){
	  
      site_total_overlap=0;
      total=0;
	  
      for(GlycosylationSiteVector::iterator current_glycosite = glycosites->begin(); current_glycosite!=glycosites->end(); ++current_glycosite){

        current_glycosite->SetChi1Value(population[member].gene[total]);
        current_glycosite->SetChi2Value(population[member].gene[total+1]);

        total=total+2;

        site_total_overlap = current_glycosite->Calculate_bead_overlaps();
        total_system_overlap = total_system_overlap + site_total_overlap;
      }
      population[member].fitness=-total_system_overlap;
  
     } // member 
	
// printout
	
//  cout << "population fitness after crossover" << endl;

//  for(int i=0; i<POPSIZE; i++){
//    cout << population[i].fitness << endl;
//  }
//  cout << endl;


  return;
}
//*****************************************************************

void elitist ( ){

//  Purpose:
//
//    ELITIST stores the best member of the previous generation.
//
//  Discussion:
//
//    The best member of the previous generation is stored as 
//    the last in the array. If the best member of the current 
//    generation is worse then the best member of the previous 
//    generation, the latter one would replace the worst member 
//    of the current population.
//
//  Local parameters:
//
//    Local, double BEST, the best fitness value.
//
//    Local, double WORST, the worst fitness value.
//

/*

  int worst_mem, best_mem;
  int i;  

  double best = population[0].fitness;
  double worst = population[0].fitness;

// double best=-10000;
// double worst=-100000000;
  
  for(i = 0; i < POPSIZE-1 ; ++i ){
     if(population[i+1].fitness < population[i].fitness ){
       if(best <= population[i].fitness ){
         best = population[i].fitness;
         best_mem = i;
       }
      if(population[i+1].fitness <= worst ){
        worst = population[i+1].fitness;
        worst_mem = i + 1;
      }
     }
     else{
         if(population[i].fitness <= worst){
           worst = population[i].fitness;
           worst_mem = i;
         }
         if(best <= population[i+1].fitness){
           best = population[i+1].fitness;
           best_mem = i + 1;
         }
     }
  }

// 
//  If the best individual from the new population is better than 
//  the best individual from the previous population, then 
//  copy the best from the new population; else replace the 
//  worst individual from the current population with the 
//  best one from the previous generation                     
//

  if(population[POPSIZE-1].fitness <= best){
    for(i = 0; i < NVARS; i++){
       population[POPSIZE-1].gene[i] = population[best_mem].gene[i];
    }
    population[POPSIZE-1].fitness = population[best_mem].fitness;
  }
  else
  {
   for(i = 0; i < NVARS; i++){
      population[worst_mem].gene[i] = population[POPSIZE-1].gene[i];
   }
   population[worst_mem].fitness = population[POPSIZE-1].fitness;
  } 

*/

  if(previous_best[0].fitness>population[POPSIZE-1].fitness){
    population[POPSIZE-1]=previous_best[0];
    population[POPSIZE-1].fitness=previous_best[0].fitness;
  }

  if(previous_best[0].fitness<population[POPSIZE-1].fitness){
    previous_best[0]=population[POPSIZE-1];
    previous_best[0].fitness=population[POPSIZE-1].fitness;
  }

  cout << endl;

  cout << " elitist " << endl;
  
  for(int i=0; i<POPSIZE; i++){
    cout << population[i].fitness << endl;
  }


  return;
}


//********************************************************************
//
// evaluate the population

void evaluate (Assembly *glycoprotein, GlycosylationSiteVector *glycosites){

    //
    //    EVALUATE implements the user-defined valuation function
    //
    //    Uses Oliver's code.
    //

    double site_total_overlap = 0.0;
    site_glycan_overlap = 0.0;
    site_protein_overlap = 0.0;
    new_dihedral_value = 0.0;
    total_system_overlap = 0.0;

    int total=0;

    ofstream geneticoutput;

    geneticoutput.open("genetic_algorithm_output.txt",ios::out|ios::app);


    // first calculate the fitnesses, again, place them in the .fitness part of the
    //  structure

    for(int member = 0; member < POPSIZE; member++){

        total=0;

        total_system_overlap=0;

        for(GlycosylationSiteVector::iterator current_glycosite = glycosites->begin(); current_glycosite != glycosites->end(); ++current_glycosite){

            std::cout << "setting population[" << member << "].gene[" << total << "]:" << population[member].gene[total] << ", " << population[member].gene[total+1] << "\n";
            current_glycosite->SetChi1Value(population[member].gene[total]);
            current_glycosite->SetChi2Value(population[member].gene[total+1]);

            total=total+2;

            site_total_overlap = current_glycosite->Calculate_bead_overlaps();
            total_system_overlap = total_system_overlap + site_total_overlap;
        }
        population[member].fitness=-total_system_overlap;
    }




    // next reorder the population in accordance with fitness - bad to best

    for(int i=0; i<POPSIZE; i++)
    {
        cout << "Member " << i << " fitness " << population[i].fitness << " gene[1]: " << population[i].gene[1] << endl;
        population_reordered[i]=population[i];
        population_reordered[i].fitness=population[i].fitness;
    }

    double temp=0;

    for(int i=0;i<=(int) POPSIZE;i++)
    {
        for(int j=0;j<=(int) POPSIZE-i;j++)
        {
            if(population_reordered[j].fitness>population_reordered[j+1].fitness)
            {
                temp=population_reordered[j].fitness;
                temp_population_reordered=population_reordered[j];

                population_reordered[j].fitness=population_reordered[j+1].fitness;
                population_reordered[j]=population_reordered[j+1];

                population_reordered[j+1].fitness=temp;
                population_reordered[j+1]=temp_population_reordered;

            }
        }
    }

    // population is reordered

    for(int i=0; i<POPSIZE; i++){
        population[i]=population_reordered[i];
        population[i].fitness=population_reordered[i].fitness;
        cout << "Member " << i << " fitness " << population[i].fitness << " gene[1]: " << population[i].gene[1] << endl;
    }




    // check

    // this should agree with reordered member=0 fitness

    // does not

    total=0;

    double total_overlap=0;
    total_system_overlap = 0;

    int mem=0;

    for(GlycosylationSiteVector::iterator current_glycosite = glycosites->begin(); current_glycosite != glycosites->end(); current_glycosite++){
        std::cout << "setting population[" << mem << "].gene[" << total << "]:" << population[mem].gene[total] << ", " << population[mem].gene[total+1] << "\n";
        //std::cout << "settting " << population[mem].gene[total] << ", " << population[mem].gene[total+1] << "\n";
        current_glycosite->SetChi1Value(population[mem].gene[total]);
        current_glycosite->SetChi2Value(population[mem].gene[total+1]);

        total=total+2;

        site_total_overlap = current_glycosite->Calculate_bead_overlaps();
        total_system_overlap = total_system_overlap + site_total_overlap;
    }

    cout << endl << "check " << -total_system_overlap << endl << endl;



    // printout reordered population according to fitness

    cout << "population reordered " << endl << endl;

    for(int i=0; i<POPSIZE; i++){
        cout << "member " << i << " fitness " << population_reordered[i].fitness << " " << population[i].fitness << endl;
        cout << "member " << i << " gene " << population_reordered[i].gene[0] << " " << population[i].gene[0] << endl;
    }
    cout << endl << endl;




    site_total_overlap=0;

    // calculate fitness

    for(int mem=0; mem<POPSIZE; mem++){

        total=0;
        geneticoutput << "member " << mem << endl << endl;

        double total_overlap=0;

        // does not agree with previous fitness

        for(GlycosylationSiteVector::iterator current_glycosite = glycosites->begin(); current_glycosite != glycosites->end(); ++current_glycosite){

            current_glycosite->SetChi1Value(population[mem].gene[total]);
            current_glycosite->SetChi2Value(population[mem].gene[total+1]);

            total=total+2;

            site_total_overlap = current_glycosite->Calculate_bead_overlaps();
            total_overlap=total_overlap+site_total_overlap;

            geneticoutput << "glycan " << total/2 << " : " << - site_total_overlap << endl;

        }

        population[mem].fitness=-total_overlap;
        geneticoutput << endl << "total overlap : " << population[mem].fitness << endl;
        geneticoutput << endl << endl;

        cout << mem << " " << population[mem].fitness << endl;

    }

    // OG added:
    for(int i=0; i<POPSIZE; i++)
    {
        cout << "member " << i << " " << population_reordered[i].fitness << " " << population[i].fitness << endl;
    }
    cout << endl << endl;

    geneticoutput.close();

    return;
}

//****************

void initialize ( std::string filename, int &seed ){

//
//  Purpose:
//
//    INITIALIZE initializes the genes within the variables bounds.
//
//  Discussion:
//
//    It also initializes (to zero) all fitness values for each
//    member of the population. It reads upper and lower bounds
//    of each variable from the input file `gadata.txt'. It
//    randomly generates values between these bounds for each
//    gene of each genotype in the population. The format of
//    the input file `gadata.txt' is
//
//      var1_lower_bound var1_upper bound
//      var2_lower_bound var2_upper bound ...

//  Parameters:
//
//    Input, string FILENAME, the name of the input file.
//
//    Input/output, int &SEED, a seed for the random number generator.
//

  int i;
  ifstream input;
  int j;
  double lbound=0;
  double ubound=360;

  srand(time(NULL));

    for(j = 0; j < POPSIZE; j++){
        for(int i=0; i<NVARS/2; i++){
           population[j].fitness = 0;
           population[j].rfitness = 0;
           population[j].cfitness = 0;
           // population[j].lower[2*i-1] = lbound;
           // population[j].upper[2*i-1]= ubound;

          double x=(double) rand()/RAND_MAX;
		  
          if(x<=(double) 1/3){
            population[j].gene[2*i] =180-30+(double) angle_change*round((double) 60/angle_change*rand()/RAND_MAX);
          }

          if(x<=(double) 2/3&& x>(double) 1/3){
            population[j].gene[2*i] =60-30+(double) angle_change*round((double) 60/angle_change*rand()/RAND_MAX);
          }

          if(x>2/3){
            population[j].gene[2*i] =-60-30+(double) angle_change*round((double) 60/angle_change*rand()/RAND_MAX);
          }

      }
    }

    // chi2

    for(j = 0; j < POPSIZE; j++){
        for(int i=0; i<NVARS/2; i++){
           population[j].fitness = 0;
           population[j].rfitness = 0;
           population[j].cfitness = 0;
           population[j].lower[2*i+1] = lbound;
           population[j].upper[2*i+1]= ubound;
           population[j].gene[2*i+1] =round((double) rand()/RAND_MAX * 360) - 180;

       }
    }
	
    //    for(int i=0; i<POPSIZE; i++){
    //       for(int i2=0; i2<NVARS; i2++){
    //          cout << "initial mem " << i << " var " << i2 << " : " << population[i].gene[i2] << endl;
    //       }
    //    }   
    //    cout << endl;
	
    return; 

}

	 
//***************************************	 
	 
void keep_the_best ( ){

//
//  Purpose:
//
//    KEEP_THE_BEST keeps track of the best member of the population.
//

  double cur_best=0;

  int cur_best_individual=0;
  cur_best = -100000000;
  
  for(int mem = 0; mem < POPSIZE; mem++){
    if(population[mem].fitness > cur_best){
      cur_best_individual=mem;
      cur_best = population[mem].fitness;
      population[POPSIZE-1].fitness = cur_best;
    }
  }
  //  cout << endl;
  //  cout << cur_best << " " << cur_best_individual << endl << endl;

  //  for(int i=0; i<POPSIZE; i++){
  //    cout << population[i].fitness << endl;
  // }
  
//
//  Once the best member in the population is found, copy the genes 
//   into the last member.
//

  for(int i2 = 0; i2 < NVARS; i2++){
     population[POPSIZE-1].gene[i2] = population[cur_best_individual].gene[i2];
  }

  return;
}


//********************************************************************
// mutation function

void mutate (Assembly *glycoprotein, GlycosylationSiteVector *glycosites ){

//
//  Purpose:
//
//    MUTATE performs a random uniform mutation.
//
//
//  Input/output, int &SEED, a seed for the random number generator.
//

  const double a = 0.0;
  const double b = 1.0;
  int i=0;
  int j=0;
  double lbound=0;
  double ubound=0;
  double x=0;
  double site_total_overlap=0;
  int total=0;
  
  double first_chi1_angle,first_chi2_angle;
  double temp_chi1_angle,temp_chi2_angle;
  double best_chi1,best_chi2;
  double min_overlap,temp_overlap;
  
  
  //  cout << "population fitness before mutation" << endl << endl;

  //  for(int i=0; i<POPSIZE; i++){
  //   cout << population[i].fitness << endl;
  // }
  // cout << endl;
  
  
  // mutation

  for ( int i = 0; i < (int) POPSIZE*(1-elite_factor); i++ ){

       total=0;

       for (GlycosylationSiteVector::iterator current_glycosite = glycosites->begin(); current_glycosite != glycosites->end(); ++current_glycosite){

     current_glycosite->SetChi1Value(population[i].gene[2*total]);
     current_glycosite->SetChi2Value(population[i].gene[2*total+1]);

	 x=(double) rand()/RAND_MAX;

	 if(x<=1/3){
	   population[j].gene[2*total] =180-30+angle_change*round((double) 60/angle_change*rand()/RAND_MAX);
	 }

	 if(x<=2/3&& x>1/3){
	   population[j].gene[2*total] =60-30+angle_change*round((double) 60/angle_change*rand()/RAND_MAX);
	 }

	 if(x>2/3){
	   population[j].gene[2*total] =-60-30+angle_change*round((double) 60/angle_change*rand()/RAND_MAX);
	 }

	 population[i].gene[2*total+1] =  ((double) rand()/RAND_MAX * 360) + 1 - 180;
	 
       } // end glycan
     } // end i 

     
/*
     cout << "pop member " << i << endl;
 
     for(GlycosylationSiteVector::iterator current_glycosite = glycosites->begin(); current_glycosite != glycosites->end(); ++current_glycosite){

//        current_glycosite->SetChi1Value(population[i].gene[total]);
		
//        current_glycosite->SetChi2Value(population[i].gene[total+1]);
	
	min_overlap=1000000;
		
        site_total_overlap = current_glycosite->Calculate_bead_overlaps();

//	   first_chi1_angle=population[i].gene[total];
//	   first_chi2_angle=population[i].gene[total+1];
			
        x = (float) rand()/RAND_MAX;
	   
      if(x < PMUTATION){
        if(site_total_overlap>minimum_fitness){

// test to find the minimal fitness placement of mutation of angle
// for a glycosite, including both chi1 and chi2
		
	  for(temp_chi1_angle=150; temp_chi1_angle<=210; temp_chi1_angle=temp_chi1_angle+angle_change){		  
   	     for(temp_chi2_angle=-180; temp_chi2_angle<=180; temp_chi2_angle=temp_chi2_angle+angle_change){

		 population[i].gene[total]=temp_chi1_angle;
		 population[i].gene[total+1]=temp_chi2_angle;
		  
         current_glycosite->SetChi1Value(population[i].gene[total]);
         current_glycosite->SetChi2Value(population[i].gene[total+1]);
		  
		 temp_overlap=current_glycosite->Calculate_bead_overlaps();
			 
		 if(temp_overlap<min_overlap){
		   best_chi1=temp_chi1_angle;
		   best_chi2=temp_chi2_angle;
		   min_overlap=temp_overlap;
		 }	 
	      }
	  }
		
	  for(temp_chi1_angle=-90; temp_chi1_angle<=-30; temp_chi1_angle=temp_chi1_angle+angle_change){		  
   	     for(temp_chi2_angle=-180; temp_chi2_angle<=180; temp_chi2_angle=temp_chi2_angle+angle_change){

		 population[i].gene[total]=temp_chi1_angle;
		 population[i].gene[total+1]=temp_chi2_angle;
		  
         current_glycosite->SetChi1Value(population[i].gene[total]);
         current_glycosite->SetChi2Value(population[i].gene[total+1]);
		  
		 temp_overlap=current_glycosite->Calculate_bead_overlaps();
			 
		 if(temp_overlap<min_overlap){
		   best_chi1=temp_chi1_angle;
		   best_chi2=temp_chi2_angle;
		   min_overlap=temp_overlap;
		 }	 
	     }
	  }
		
	  for(temp_chi1_angle=30; temp_chi1_angle<=90; temp_chi1_angle=temp_chi1_angle+angle_change){		  
   	     for(temp_chi2_angle=-180; temp_chi2_angle<=180; temp_chi2_angle=temp_chi2_angle+angle_change){

		 population[i].gene[total]=temp_chi1_angle;
		 population[i].gene[total+1]=temp_chi2_angle;
		  
         current_glycosite->SetChi1Value(population[i].gene[total]);
         current_glycosite->SetChi2Value(population[i].gene[total+1]);
		  
		 temp_overlap=current_glycosite->Calculate_bead_overlaps();
			 
		 if(temp_overlap<min_overlap){
		   best_chi1=temp_chi1_angle;
		   best_chi2=temp_chi2_angle;
		   min_overlap=temp_overlap;
		 }	 
	     }
	  }

       } // mutation test 
    }

    cout << "done " << endl;
    cout << "total " << total << endl;
    cout << "member " << i << endl;

    population[i].gene[total]=best_chi1;
    population[i].gene[total+1]=best_chi2;
	  
    current_glycosite->SetChi1Value(population[i].gene[total]);
    current_glycosite->SetChi2Value(population[i].gene[total+1]);
	
    total=total+2;
   
       } 
    }

*/

	
//  cout << "population fitness after mutation" << endl << endl;

//  for(int i=0; i<POPSIZE; i++){
//    cout << population[i].fitness << endl;
//  }
//  cout << endl;

  return;
}


//**********************************************************

void report ( int generation, Assembly *glycoprotein, GlycosylationSiteVector *glycosites ){

//
//  Purpose:
//
//    REPORT reports progress of the simulation.
//
//  Local parameters:
//
//    Local, double avg, the average population fitness.
//
//    Local, best_val, the best population fitness.
//
//    Local, double square_sum, square of sum for std calc.
//
//    Local, double stddev, standard deviation of population fitness.
//
//    Local, double sum, the total population fitness.
//
//    Local, double sum_square, sum of squares for std calc.
//


// population has changed to unordered


//  cout << "populaton in report" << endl;

  for(int i=0; i<POPSIZE; i++){
//    population[i]=population_reordered[i];
    cout << population[i].fitness << endl;
  }
  cout << endl;

  cout << population[POPSIZE-1].fitness << " " << previous_best[0].fitness << endl;

  

  double best_fitness=population[POPSIZE-1].fitness;
  double avg=0;
  double best_val=0;
  int i;
  double square_sum=0;
  double stddev=0;
  double test_sum=0;
  double test_sum_square=0;
  double current=0;

  if(generation == 0){
    cout << "\n";
    cout << "  Generation       Best            Average       Standard       time \n";
    cout << "  number           value           fitness       deviation      calculation \n";
    cout << "\n";
  }

  test_sum = 0.0;
  test_sum_square = 0.0;

  for(i = 0; i < POPSIZE; i++){

     test_sum = test_sum + population[i].fitness;
     test_sum_square = test_sum_square + population[i].fitness * population[i].fitness;

  }

  avg = (float) test_sum / POPSIZE;
  square_sum = avg * avg * POPSIZE;
  stddev = sqrt ( ( test_sum_square - square_sum ) / ( POPSIZE - 1 ) );

  t=time(0)-timezero;

  ofstream geneticoutput;

  geneticoutput.open("genetic_algorithm_output.txt",ios::out|ios::app);

  cout << "  " << setw(8) << generation
       << "  " << setw(14) << best_fitness
       << "  " << setw(14) << (float) test_sum/POPSIZE
       << "  " << setw(14) << stddev
       << "  " << t << "\n";
  
  cout << endl;

  geneticoutput << "  " << setw(8) << generation
       << "  " << setw(14) << best_fitness
       << "  " << setw(14) << (float) test_sum/POPSIZE
       << "  " << setw(14) << stddev
       << "  " << t << "\n";

  geneticoutput.close();


  //  cout << "previous_population in report" << endl;

  for(int i = 0; i < POPSIZE*elite_factor; i++){
    previous_population[i] = population[(int) (POPSIZE*(1-elite_factor)+i)];
  }
  
  //  cout << endl;
  //  cout << previous_population[0].fitness << "  " << previous_population[1].fitness << endl;
  //  cout << endl;

  return;
}
//*****************************************************

void selector ( int &seed ){

//
//  Purpose:
//
//    SELECTOR is the selection function.
//
//  Discussion:
//
//    Standard proportional selection for maximization problems incorporating
//    the elitist model.  This makes sure that the best member always survives.
//  Parameters:
//
//    Input/output, int &SEED, a seed for the random number generator.
//

  const double a = 0.0;
  const double b = 1.0;
  int i;
  int j;
  int mem;
  double p;
  double sum;


  for(int i = 0; i < POPSIZE*elite_factor; i++){
     previous_population[i] = population[(int) (POPSIZE*(1-elite_factor)+i-1)];
  }

//
//  Find the total fitness of the population.
//
  sum = 0.0;
  for(mem = 0; mem < POPSIZE; mem++){
     sum = sum + population[mem].fitness;
  }
//
//  Calculate the relative fitness of each member.
//
  for(mem = 0; mem < POPSIZE; mem++){
     population[mem].rfitness = population[mem].fitness / sum;
  }
//
//  Calculate the cumulative fitness.
//
  population[0].cfitness = population[0].rfitness;
  for(mem = 1; mem < POPSIZE; mem++){
     population[mem].cfitness = population[mem-1].cfitness +
     population[mem].rfitness;
  }
//
//  Select survivors using cumulative fitness.
//
  for(i = 0; i < POPSIZE; i++){

  //     p = r8_uniform_ab ( a, b, seed);

   p=(double) rand()/RAND_MAX;

     if(p < population[0].cfitness){
       newpopulation[i] = population[0];
       }
     else{
         for(j = 0; j < POPSIZE; j++){
            if(population[j].cfitness <= p && p < population[j+1].cfitness){
              newpopulation[i] = population[j+1];
            }
         }
      }
   }
//
//  Overwrite the old population with the new one.
//
  for(i = 0; i < POPSIZE; i++){
     population[i] = newpopulation[i];
     population[i].fitness=newpopulation[i].fitness;
  }

  return;
}

//******************************************************************

void timestamp ( ){

//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    May 31 2001 09:45:54 AM
//
//   

# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  //  cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
