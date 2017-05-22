#include "genetic.h"
#include <stdio.h>
#include <stdlib.h>

#define TAM 32
#define EVOLS 500

double GA_objective                    (GA *ag, Chromosome x, long length);

int main()
{
  GA *g;
  int i, j, gen, elit;
  long best[TAM];
  double fitness = 0.0;
  double cross, mutan;

//  printf("Prob cross = "); scanf("%lf", &cross);
//  printf("Prob mutan = "); scanf("%lf", &mutan);
  printf("Elitismo = "); scanf("%d", &elit);

  g = GA_initialization(0, TAM<10 ? 20 : TAM, TAM, 0.65, 1.0/TAM);
  g->elitism = elit;
  g->crossover_type = 0;

  i=0;
  do {
    i++;
    GA_generational_step(g);
    if (fitness<g->best_fitness) {
      for (j=0; j<g->length; j++)
        best[j] = g->nueva->ind [g->best].chrom [j];
      gen = i;
      fitness = g->best_fitness;
    }
    if (i%50==0) {
/*      printf ("G%05d:\n", i);
      for (j=0; j<g->length; j++)
        printf ("%ld ", g->nueva->ind [g->best].chrom [j]);
      printf (" %lf\n", g->nueva->ind [g->best].fitness);
  */  }
  } while (i<EVOLS);// && g->gen_without_improve < 1000);
 
  printf ("B%05d: ", gen);
  for (j=0; j<g->length; j++)
    printf ("%ld ", best [j]);
  printf (" %lf\n", fitness);

  free(g); 

  return 0;
}

/*double GA_objective                    (GA *ag, Chromosome x, long length)
{
  int i;
  double fitness = 0.0;

  for (i=0;i<length;i++) 
    if (1==x[i]) fitness += 1.0/TAM;

  return 1.0-fitness;
}*/
/*double GA_objective                    (GA *ag, Chromosome x, long length)
{
  int i;
  double fitness = 0.0;

  for (i=0;i<length;i++) 
    if (i%2==x[i]) fitness += 1.0/TAM;

  return 1.0-fitness;
}*/
double GA_objective                    (GA *ag, Chromosome x, long length)
{
  int i;
  double fitness = 1.0/TAM;

  for (i=1;i<length;i++) 
    if (x[i]!=x[i-1]) fitness += 1.0/TAM;

  return 1.0-fitness;
}

