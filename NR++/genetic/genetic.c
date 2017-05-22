/*#include "strategy.h"*/
#include "genetic.h"
#include "sort2.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/*
    Allocates one chromosome 
*/

Chromosome GA_allocate_chromosome (long length)
{
  Chromosome c = NULL;
  
  c = (GA_Allele *) malloc (length * sizeof (GA_Allele));
  
  return c;
}

/*  
    Allocates individuals of old and new pops
    Initializes old population
    Returns pointer to GA
*/

GA * GA_initialization (short     selection, 
                        long      size, 
                        long      length, 
                        double    prob_cross, 
                        double    prob_mut)
{
  long   i, j, cont;

  GA * g = (GA *) malloc (sizeof (GA));

  size = (size % 2 == 0) ? size : size +1;

  /* Type of selection */
  g->selection = selection;

  /* Other default options */
  g->crossover_type = 0; /*uniform-point*/
  g->elitism = 0;        /*no elitism used*/

  /* Seting of probabilities of crossover and mutation  */
  g->prob_cross = prob_cross;
  g->prob_mut = prob_mut;

  /* Initialization of number of crossovers and mutations  */
  g->n_cross = 0;
  g->n_mut = 0;
  g->gen_number = 0;

  /* Seting of length and size  */
  g->length = length;
  g->pop_size = size;

  /* Old and new populations  */
  g->old = (Population *) malloc (sizeof (Population));
  g->nueva = (Population *) malloc (sizeof (Population));

  g->old->pop_size = size;
  g->nueva->pop_size = size;

  /* Allocation of individuals  */
  g->old->ind = (Individual *) malloc (size * sizeof (Individual));
  g->nueva->ind = (Individual *) malloc (size * sizeof (Individual));

  for (i=0; i<size; i++)
  {
    g->old->ind [i].chrom = GA_allocate_chromosome (length);
    g->nueva->ind [i].chrom = GA_allocate_chromosome (length);
  }

  /* Initialization of individuals  */
  for (i=0; i<size; i++)
  {
    g->old->ind [i].length = length;
    g->nueva->ind [i].length = length;
  }

  /* Semilla  */
  /* srand (time (NULL) / 2); */

  /* Generation of diagonal initial population */
  for (i=0; i<size; i++)
  {
    for (j=0; j<length; j++)
    {
      /*  cont es el numero de individuos modulo la longitud  */
      cont = i;
      while (cont >= length)
        cont -= length;

      g->old->ind [i].chrom [j] = flip(0.5);/*flip(0.5);*//*(cont == j) ? 1 : 0;*/
      /*  g->old->ind [i].chrom [j] = (short) random_integer (1);  */
#ifdef WRITE_INITIAL_POPULATION
      printf ("%d ", g->old->ind [i].chrom [j]);
#endif
    }
#ifdef WRITE_INITIAL_POPULATION
      printf ("\n");
#endif
  }
  
  g->gen_without_improve = 0;
  
  /* Calculation of fitness */
  GA_calculate_population_fitness (g, g->old);

  /* Calculate best individual  */
  GA_calculate_best_individual (g);
   
  /* Variance  */
  GA_fitness_variance (g);
 
  /* Calculation of selection probabilities */
  GA_calculate_probabilities (g->old);
  
  /* Report results 
  GA_report (g);*/

  return g;
}
 
/*  
    Frees 'g'
*/

void GA_free (GA *g)
{
  long   i;

  for (i=0; i<g->old->pop_size; i++)
  {
    free(g->old->ind [i].chrom);
    free(g->nueva->ind [i].chrom);
  }

  free(g->old->ind);
  free(g->nueva->ind);

  free(g->old);
  free(g->nueva);

  free(g);
  
}
 
GA * GA_initialization_from_file (char * name)
{
  long   i, j;
  FILE   * fp = NULL;

  GA * g = (GA *) malloc (sizeof (GA));
  fp = fopen (name, "r");

  /* Type of selection  */
  fscanf (fp, "%hd ", &(g->selection));

  /* Seting of length and size  */
  fscanf (fp, "%ld ", &(g->pop_size));
  fscanf (fp, "%ld ", &(g->length));

  /* Seting of probabilities of crossover and mutation  */
  fscanf (fp, "%lf ", &(g->prob_cross));
  fscanf (fp, "%lf ", &(g->prob_mut));

  /* Other default options */
  g->crossover_type = 1; /*one-point*/
  g->elitism = 0;        /*no elitism used*/

  /* Initialization of number of crossovers and mutations  */
  g->n_cross = 0;
  g->n_mut = 0;
  g->gen_number = 0;

  /* Old and new populations  */
  g->old = (Population *) malloc (sizeof (Population));
  g->nueva = (Population *) malloc (sizeof (Population));

  g->old->pop_size = g->pop_size;
  g->nueva->pop_size = g->pop_size;

  /* Allocation of individuals  */
  g->old->ind = (Individual *) malloc (g->pop_size * sizeof (Individual));
  g->nueva->ind = (Individual *) malloc (g->pop_size * sizeof (Individual));

  for (i=0; i<g->pop_size; i++)
  {
    g->old->ind [i].chrom = GA_allocate_chromosome (g->length);
    g->nueva->ind [i].chrom = GA_allocate_chromosome (g->length);
  }

  /* Initialization of individuals  */
  for (i=0; i<g->pop_size; i++)
  {
    g->old->ind [i].length = g->length;
    g->nueva->ind [i].length = g->length;
  }

  /* Reading initial population */
  for (i=0; i<g->pop_size; i++)
    for (j=0; j<g->length; j++)
      fscanf (fp, "%ld ", &(g->old->ind [i].chrom [j]));

  fclose (fp);

  /* Calculation of fitness */
  GA_calculate_population_fitness (g, g->old);

  /* Calculate best individual  */
  GA_calculate_best_individual (g);
   
  /* Variance  */
  GA_fitness_variance (g);
 
  /* Calculation of selection probabilities */
  GA_calculate_probabilities (g->old);
  
  /* Report results
  GA_report (g);*/

  return g;
}


void GA_fitness (GA * g, Individual * ind)
{
  ind->objective = GA_objective (g, ind->chrom, ind->length);
  ind->fitness = 1.0 / (1.0 + ind->objective);
}

/*  Calculates fitness and average value  */

void GA_calculate_population_fitness (GA * x, Population * pop)
{
  long i;

  pop->avg_fitness = 0.0;
  pop->avg_objective = 0.0;

  for (i=0; i<pop->pop_size; i++)
  {
    GA_fitness (x, &(pop->ind [i]));
    pop->avg_fitness += pop->ind [i].fitness;
    pop->avg_objective += pop->ind [i].objective;
  }
  
  pop->avg_fitness /= x->pop_size;
  pop->avg_objective /= x->pop_size;
}

/* Determines which are the best and worst individuals  */

void GA_calculate_best_individual (GA * g)
{
  long     i;
  long   * ind;
  double * fitness;
  /*static*/ double best;

  best = g->best_fitness;

  ind = (long *) malloc (g->pop_size * sizeof (long));
  fitness = (double *) malloc (g->pop_size * sizeof (double));

  for (i=0; i<g->pop_size; i++)
  {
    ind [i] = i;
    fitness [i] = g->old->ind [i].fitness;
  }

  sort2 (g->pop_size, &(fitness[-1]), &(ind[-1]));

  g->best = ind [g->pop_size-1];
  g->worst = ind [0];
  g->best_fitness = g->old->ind [g->best].fitness;

  g->old->best_fitness = g->old->ind [g->best].fitness;
  g->old->worst_fitness = g->old->ind [g->worst].fitness;
  
  /*  Best fitness is saved in best  */
  if (g->old->best_fitness <= best)
    g->gen_without_improve++;
  else
  {
    g->gen_without_improve=0;
    best = g->old->best_fitness;
  }
  
  free (ind);
  free (fitness);
}

/* Calculates fitness variance  */

void GA_fitness_variance (GA * g)
{
  long     i;

  g->old->var = 0.0;
  for (i=0; i<g->pop_size; i++)
    g->old->var += g->old->ind [i].fitness * g->old->ind [i].fitness;
  g->old->var /= g->pop_size;
  g->old->var -= g->old->avg_fitness * g->old->avg_fitness;
}

/* Returns a random number in [0, 1)  */

double random_double ()
{
  return (1.0 * rand () / (RAND_MAX + 1.0));
}

/* Returns a random integer in [0,n]  */

long random_integer (long n)
{
  return (long) floor ( (n+1) * random_double () );
}

/* Returns 1 with probability prob  */

long flip (double prob)
{
  double pos = random_double ();

  if (pos < prob) return 1;

  return 0;
}

/* 
   Calculates probabilities for 
   proportional selection 
*/

void GA_calculate_probabilities (Population * pop)
{
  long i;

  pop->sum_fitness = 0.0;
  for (i=0; i<pop->pop_size; i++)
    pop->sum_fitness += pop->ind [i].fitness;

  for (i=0; i<pop->pop_size; i++)
    pop->ind [i].prob_proportional = pop->ind [i].fitness / pop->sum_fitness;
}

/*  
   Returns the index of an individual
   selected by proportional selection
*/

long GA_proportional_selection (Population * pop)
{
  long   i;
  double pos, val;

  pos = (1.0 * rand ()) / (RAND_MAX + 1.0);

  val = 0.0;
  i = 0;
  while (val<pos && i<pop->pop_size)
  {
    val += pop->ind [i].prob_proportional;
    i++;
  }

  if (i>0) i--;
  return i;
}

/*
   Tournament selection
*/

long GA_tournament_selection (Population * pop, long num)
{
  long   i, pos, best;
  double max = 0.0;

  for (i=0; i<num; i++)
  {
    pos = random_integer (pop->pop_size - 1);
    if (pop->ind [pos].fitness > max)
    {
      max = pop->ind [pos].fitness;
      best = pos;
    }
  }

  return best;
}

/*  
   Mutation of each bit with probability prob_mut  
*/

void GA_mutation (GA * g, Individual * x)
{
  long   i;
  
  for (i=0; i<x->length; i++)
    if (flip (g->prob_mut) == 1)
    {
      x->chrom [i] = 1 - x->chrom [i];
      g->n_mut++;
    }
 
}

/*
   Generational genetic algorithm
*/

void GA_generational_step (GA * x)
{
  long       i, j, k, mate1, mate2, child_incr;

/*  printf ("Generation %d: \n", x->gen_number);*/

  /* Selection + crossover + mutation j */
  child_incr = /*x->crossover_type == 0 ? 1 :*/ 2;
  for (j=0; j<x->pop_size; j+=child_incr)
  {
    /* Selection of parents */
    if (x->selection == 0)       // proportional
    {
      mate1 = GA_proportional_selection (x->old);
      mate2 = GA_proportional_selection (x->old);
    }
    else if (x->selection == 1)  // tournament
    {
      mate1 = GA_tournament_selection (x->old, x->tournament);
      mate2 = GA_tournament_selection (x->old, x->tournament);
    }

    /* One point crossover of parents mate1 and mate2  */
    if (flip (x->prob_cross) == 1)
    {
      if (x->crossover_type == 0)
      {
        GA_uniform_crossover (x, mate1, mate2, &(x->nueva->ind [j]), &(x->nueva->ind [j+1]));
      }
      else if (x->crossover_type == 1) 
      {
        GA_one_point_crossover (x, mate1, mate2, &(x->nueva->ind [j]), &(x->nueva->ind [j+1]));
      }
/*      else if (x->crossover_type == 2)
        GA_two_point_crossover (x, mate1, mate2, &(x->nueva->ind [j]), &(x->nueva->ind [j+1]));*/
    }
    else  // clonation
    {
      for (i=0; i<x->length; i++)
      {
        x->nueva->ind [j].chrom [i] = x->old->ind [mate1].chrom [i];
        if (child_incr == 2) 
          x->nueva->ind [j+1].chrom [i] = x->old->ind [mate2].chrom [i];
      }
      x->nueva->ind [j].parent1 = mate1;
      x->nueva->ind [j].parent2 = mate1;
      x->nueva->ind [j].xpoint = 0;
      if (child_incr == 2) 
      {
        x->nueva->ind [j+1].parent1 = mate2;
        x->nueva->ind [j+1].parent2 = mate2;
        x->nueva->ind [j+1].xpoint = 0;
      }
    }

    /* Mutation of children */
    for (k=0; k<child_incr; k++)
      GA_mutation (x, &(x->nueva->ind [j+k]));
  }
  
  /* Elistism */
  if (x->elitism == 1)
  {
     for (i=0; i<x->length; i++)
        x->nueva->ind [0].chrom [i] = x->old->ind [x->best].chrom [i];
    x->nueva->ind [0].xpoint = 0;
    x->nueva->ind [0].parent1 = x->best;
    x->nueva->ind [0].parent2 = x->best;
     
    for (i=0; i<x->length; i++)
        x->nueva->ind [1].chrom [i] = x->old->ind [x->best].chrom [i];
    x->nueva->ind [1].xpoint = 0;
    x->nueva->ind [1].parent1 = x->best;
    x->nueva->ind [1].parent2 = x->best;
  }

  /* Generational replacement of populations  */
  for (j=0; j<x->pop_size; j++)
  {
    for (i=0; i<x->length; i++)
        x->old->ind [j].chrom [i] = x->nueva->ind [j].chrom [i];
    x->old->ind [j].xpoint = x->nueva->ind [j].xpoint;
    x->old->ind [j].parent1 = x->nueva->ind [j].parent1;
    x->old->ind [j].parent2 = x->nueva->ind [j].parent2;
  }

  /* Calculation of fitness */
  GA_calculate_population_fitness (x, x->old);

  /* Calculate best individual  */
  GA_calculate_best_individual (x);
   
  /* Variance  */
  GA_fitness_variance (x);
 
  /* Calculation of selection probabilities */
  GA_calculate_probabilities (x->old);
  
  /* Report results 
  GA_report (x);*/

  /* Updating of generation number  */
  x->gen_number++;
}

/*  One point crossover function  */

void GA_one_point_crossover (GA * x, long mate1, long mate2, 
                             Individual * child1, Individual * child2)
{
  long  i, xpoint;
/*  static long cont=0;*/

  /*  xpoint belongs to [0, x->length-2]  */
  xpoint = random_integer (x->length-2);
  x->n_cross++;

  /*  Set parents and crossover points  */
  child1->parent1 = mate1;
  child1->parent2 = mate2;
  child1->xpoint = xpoint;

  child2->parent1 = mate1;
  child2->parent2 = mate2;
  child2->xpoint = xpoint;

  for (i=0; i<=xpoint; i++)
  {
    child1->chrom [i] = x->old->ind [mate1].chrom [i];
    child2->chrom [i] = x->old->ind [mate2].chrom [i];
  }

  for (i=xpoint+1; i<x->length; i++)
  {
    child1->chrom [i] = x->old->ind [mate2].chrom [i];
    child2->chrom [i] = x->old->ind [mate1].chrom [i];
  }
}

/*  Uniform crossover function  */

void GA_uniform_crossover      (GA         * x, 
                                long       mate1, 
                                long       mate2, 
                                Individual * child1,
                                Individual * child2)
{
  long  i;

  x->n_cross++;

  /*  Set parents and crossover points  */
  child1->parent1 = mate1;
  child1->parent2 = mate2;
  child1->xpoint = -1;
  child2->parent1 = mate1;
  child2->parent2 = mate2;
  child2->xpoint = -1;

  for (i=0; i<x->length; i++)
  {
    if (flip(0.5)==1) 
    {
       child1->chrom [i] = x->old->ind [mate1].chrom [i];
       child2->chrom [i] = x->old->ind [mate2].chrom [i];
    }
    else
    {
       child1->chrom [i] = x->old->ind [mate2].chrom [i];
       child2->chrom [i] = x->old->ind [mate1].chrom [i];
    }
  }
}

/*  
    Reports about individuals, 
    their fitness, parents, cross points 
    and strings  
*/

void GA_report (GA * g)
{
  //static int k=0, m=1;
  long     i, j;
  FILE    * fp = NULL, * fit = NULL;

return;

  //    Open files
  fp = fopen ("ga.out", "a+");
  fit = fopen ("fit.out", "a+");

  fprintf (fp, "--------------------------------------------------\n");
  fprintf (fp, "Generation:             %ld\n", g->gen_number);
  fprintf (fp, "Fitness variance:    %g\n", g->old->var);
  fprintf (fp, "Average fitness:     %g\n", g->old->avg_fitness);
  fprintf (fp, "Best Fitness         %g\n", g->old->best_fitness);
  fprintf (fp, "Worst Fitness        %g\n", g->old->worst_fitness);
  fprintf (fp, "Best objective:      %g\n", g->old->ind [g->best].objective);
  fprintf (fp, "Number of mutations:      %ld\n", g->n_mut);

  /* Print best chromosome  */
  for (j=0; j<g->length; j++){
//    if (k%m==0) printf ("%d ", g->old->ind [g->best].chrom [j]);
    fprintf (fp, "%ld ", g->old->ind [g->best].chrom [j]);
  }
  if (g->gen_number > 0) 
  {
    fprintf (fp, "(%3ld,", g->old->ind [g->best].parent1);
    fprintf (fp, "%3ld) ", g->old->ind [g->best].parent2);
    fprintf (fp, "-> %3ld <- ", g->old->ind [g->best].xpoint);
  }
  fprintf (fp, "%g ", g->old->ind [g->best].fitness);
  fprintf (fp, "\n");
//  if (k%m==0) printf (" <--  %g\n", g->old->ind [g->best].fitness);
  //k++;


  /* Print average and best fitness and objective  */
  fprintf (fit, "%g %g %g %g ", 
      g->old->avg_fitness, 
      g->old->ind [g->best].fitness,
      g->old->avg_objective, 
      g->old->ind [g->best].objective);

#ifdef REPORT_LONG
  for (i=0; i<g->pop_size; i++)
  {
    fprintf (fp, "%3ld : ", i);
   /* for (j=0; j<g->length; j++)
      fprintf (fp, "%d ", g->old->ind [i].chrom [j]);
    fprintf (fp, "(%3ld,", g->old->ind [i].parent1);
    fprintf (fp, "%3ld) ", g->old->ind [i].parent2);
    fprintf (fp, "-> %3ld <- ", g->old->ind [i].xpoint);
*/    fprintf (fp, "%.6g ", g->old->ind [i].objective);
    fprintf (fp, "%.6g ", g->old->ind [i].fitness);
    fprintf (fp, "%.6g ", g->old->ind [i].prob_proportional);
    fprintf (fp, "%3ld ", (long ) floor (g->old->ind [i].prob_proportional * g->pop_size));
    if (i == g->best)
      fprintf (fp, "<-- best");
    fprintf (fp, "\n");
  }
#endif
  fprintf (fp, "--------------------------------------------------\n");

  fclose (fp);
  fclose (fit);
}

/*  
    Whole genetic algorithm
    Returning best fitness
*/

double GA_evolution (short selection, long size, long length, double prob_cross, double prob_mut)
{
  GA * g;

  g = GA_initialization (selection, size, length, prob_cross, prob_mut);

  do 
  {
    GA_generational_step (g);
  } while (g->gen_without_improve < 10);

  return g->best_fitness;
}

/*  Hill climbing  */

void GA_hill_climbing (GA * g)
{
  long   c [500];
  long   i;
  double fitness, best_fitness;

  /*  best_fintess is the fitness of the best individual  */
  best_fitness = g->old->best_fitness;

  for (i=0; i<g->length; i++)
    c [i] = g->old->ind [g->best].chrom [i];

  for (i=0; i<g->length; i++)
  {
    if (c [i] == 1)
    {
      c [i] = 0;
      fitness = GA_objective (g, c, g->length);
      if (fitness < best_fitness)
      {
        best_fitness = fitness;
      }
      else 
        c [i] = 1;
    }
  }
  
  for (i=0; i<g->length; i++)
    g->old->ind [g->best].chrom [i] = c[i];
}


