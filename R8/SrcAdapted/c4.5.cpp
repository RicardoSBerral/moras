/*************************************************************************/
/*									 */
/*	Main routine, c4.5						 */
/*	------------------						 */
/*									 */
/*************************************************************************/


#include "c4.5.h"

#include "getnames.h"
#include "getdata.h"
#include "trees.h"
#include "getopt.h"
#include "header.h"
#include "besttree.h"
#include "build.h"
#include "st-thresh.h"
#include "stats.h"


#include <stdlib.h>


    /*  External data, described in extern.i  */
namespace c45
{

short		MaxAtt, MaxClass, MaxDiscrVal = 2;

ItemNo		MaxItem;

Description	*Item = 0;

DiscrValue	*MaxAttVal = 0;

char		*SpecialStatus = 0;

String		*ClassName = 0,
		*AttName = 0,
		**AttValName = 0,
		FileName = "DF",
                WeightsFileName = "";

short		VERBOSITY = 0,
		TRIALS    = 10;

Boolean		GAINRATIO  = true,
		SUBSET     = false,
		BATCH      = true,
		UNSEENS    = false,
		PROBTHRESH = false,
		MDLPENALTY = true;

ItemNo		MINOBJS   = 2,
		WINDOW    = 0,
		INCREMENT = 0;

float		CF = 0.25;

Tree		*Pruned = 0;

Boolean		AllKnown = true;

void SetGlobalOpt(int opt, char *arg)
{
	switch (opt)
	{
	case 'f':   FileName = arg;
/*		    printf("\tFile stem <%s>\n", FileName);*/
		    break;
	case 'W':   WeightsFileName = arg;
/*		    printf("\tFile stem <%s>\n", FileName);*/
		    break;
	case 'b':   BATCH = true;
		    printf("\tWindowing disabled (now the default)\n");
		    break;
	case 'u':   UNSEENS = true;
		    printf("\tTrees evaluated on unseen cases\n");
		    break;
	case 'd':   MDLPENALTY = MDLPENALTY ? false : true;
		    printf("\tMDL-Penalty is %s\n", MDLPENALTY ? "ON" : "OFF");
		    break;
	case 'p':   PROBTHRESH = true;
		    printf("\tProbability thresholds used\n");
		    break;
	case 'v':   VERBOSITY = atoi(arg);
		    printf("\tVerbosity level %d\n", VERBOSITY);
		    break;
	case 't':   TRIALS = atoi(arg);
		    printf("\tWindowing enabled with %d trials\n", TRIALS);
		    Check(TRIALS, 1, 10000);
		    BATCH = false;
		    break;
	case 'w':   WINDOW = atoi(arg);
		    printf("\tInitial window size of %d items\n", WINDOW);
		    Check(WINDOW, 1, 1000000);
		    BATCH = false;
		    break;
	case 'i':   INCREMENT = atoi(arg);
		    printf("\tMaximum window increment of %d items\n",
			   INCREMENT);
		    Check(INCREMENT, 1, 1000000);
		    BATCH = false;
		    break;
	case 'g':   GAINRATIO = false;
		    printf("\tGain criterion used\n");
		    break;
	case 's':   SUBSET = true;
		    printf("\tTests on discrete attribute groups\n");
		    break;
	case 'm':   MINOBJS = atoi(arg);
		    printf("\tSensible test requires 2 branches with >=%d cases\n",
			    MINOBJS);
		    Check(MINOBJS, 1, 1000000);
		    break;
	case 'c':   CF = atof(arg);
		    printf("\tPruning confidence level %g%%\n", CF);
		    Check(CF, Epsilon, 100);
		    CF /= 100;
                    ResetCoeff();
		    break;
	case '?':   printf("unrecognised option\n");
		    exit(1);
	}
}
void InitParams(int Argc, char *Argv[])
{
    int o;
    extern  char *optarg;
/*    extern  int c45::optind;*/
    Boolean FirstTime=true;

    while ( (o = getopt(Argc, Argv, "f:budpv:t:w:i:gsm:c:")) != EOF )
    {
	if ( FirstTime )
	{
	    printf("\n    Options:\n");
	    FirstTime = false;
	}
        SetGlobalOpt(o, optarg);
    }

}
int    run(int Argc, char *Argv[])
/*  ----  
    int Argc;
    char *Argv[];*/
{
    short Best/*, BestTree()*/;

    PrintHeader("decision tree generator");

    /*  Process options */ 
    InitParams(Argc, Argv);
    
    /*  Initialise  */
    GetNames();
    GetData(".data");
    printf("\nRead %d cases (%d attributes) from %s.data\n",
	   MaxItem+1, MaxAtt+1, FileName);

    /*  Build decision trees  */
    if ( BATCH )
    {
	TRIALS = 1;
	OneTree();
	Best = 0;
    }
    else
    {
	Best = BestTree();
    }

    /*  Soften thresholds in best tree  */
    if ( PROBTHRESH )
    {
	printf("Softening thresholds");
	if ( ! BATCH ) printf(" for best tree from trial %d", Best);
	printf("\n");
	SoftenThresh(Pruned[Best]);
	printf("\n");
	PrintTree(Pruned[Best]);
    }

    /*  Save best tree  */
    if ( BATCH || TRIALS == 1 )
    {
	printf("\nTree saved\n");
    }
    else
    {
	printf("\nBest tree from trial %d saved\n", Best);
    }
    SaveTree(Pruned[Best], ".tree");

    /*  Evaluation  */
    printf("\n\nEvaluation on training data (%d items):\n", MaxItem+1);
    Evaluate(false, Best);

    if ( UNSEENS )
    {   
        GetData(".test");
        printf("\nEvaluation on test data (%d items):\n", MaxItem+1);
        Evaluate(true, Best);
    }

    exit(0);

    return 0;
}
}/*namespace c45*/

#ifdef RUNRUN
int    main(int Argc, char *Argv[])
{
  return c45::run(Argc, Argv);
}
#endif

