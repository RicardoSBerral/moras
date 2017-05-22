/*#define REPORT_LONG*/

typedef long GA_Allele;

typedef GA_Allele * Chromosome;

typedef struct {
  Chromosome chrom;
  long       length;
  double     fitness;
  double     objective;
  double     prob_proportional;
  long       parent1, parent2, xpoint;
} Individual;

typedef struct {
  Individual   * ind;
  long         pop_size;
  double       sum_fitness;
  double       avg_fitness;
  double       avg_objective;
  double       best_fitness;
  double       worst_fitness;
  double       var;
} Population;

typedef struct {
  Population   * old, * nueva;
  long         length;
  long         pop_size;
  long         best;
  long         worst;
  double       best_fitness;
  double       prob_cross;
  long         n_cross;
  double       prob_mut;
  long         n_mut;    
  long         gen_number;
  short        selection;
  long         tournament;
  int          elitism;
  long         n_parameters;
  double       min;
  double       max;
  long         gen_without_improve;
  long         crossover_type; /* 0=uniform; 1=one point; 2=two point*/
} GA;
 
/* Initial population */
GA * GA_initialization         (short      selection,
                                long       size, 
                                long       length, 
                                double     prob_cross, 
                                double     prob_mut);

void GA_free(GA *g); 

GA * GA_initialization_from_file 
                               (char       * name);

double GA_evolution            (short      selection,
                                long       size, 
                                long       length, 
                                double     prob_cross, 
                                double     prob_mut);

/* Random number generation */
long   random_integer          (long       n);
double random_double           ();
long   flip                    (double     prob);

/* Selection function */
void GA_calculate_probabilities(Population * pop);
long GA_proportional_selection (Population * pop);
long GA_tournament_selection   (Population * pop, long size);

/* Crossover operators */
void GA_one_point_crossover    (GA         * x, 
                                long       mate1, 
                                long       mate2, 
                                Individual * child1, 
                                Individual * child2);
void GA_uniform_crossover      (GA         * x, 
                                long       mate1, 
                                long       mate2, 
                                Individual * child1, 
                                Individual * child2);
 
/* Mutation operators  */
void GA_mutation               (GA         * g, 
                                Individual * x);

/* Replacement strategies */
void GA_generational_step      (GA         * x);
 
/* Allocation */
Chromosome GA_allocate_chromosome (long      length);

/* Fitness function */
double GA_objective                    (GA *, Chromosome x, long length);
void   GA_fitness                      (GA *, Individual * ind);
void   GA_calculate_population_fitness (GA *, Population * pop);
void   GA_calculate_best_individual    (GA *);
void   GA_fitness_variance             (GA * g);
void   GA_hill_climbing                (GA * g);

/* Reports  */   
void   GA_report                (GA * g);

