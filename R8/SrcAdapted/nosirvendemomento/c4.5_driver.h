/*************************************************************************/
/*									 */
/*	Main routine, c4.5						 */
/*	------------------						 */
/*									 */
/*************************************************************************/


#include "defns.i"
#include "types.i"
#include "extern.i"

    /*  External data, described in extern.i  */

short		MaxAtt, MaxClass, MaxDiscrVal = 2;

ItemNo		MaxItem;

Description	*Item;

DiscrValue	*MaxAttVal;

char		*SpecialStatus;

String		*ClassName,
		*AttName,
		**AttValName,
		FileName = "DF";

short		VERBOSITY = 0,
		TRIALS    = 10;

Boolean		GAINRATIO  = true,
		SUBSET     = false,
		BATCH      = true,
		UNSEENS    = false,
		PROBTHRESH = false;

ItemNo		MINOBJS   = 2,
		WINDOW    = 0,
		INCREMENT = 0;

float		CF = 0.25;

Tree		*Pruned;

Boolean		AllKnown = true;


int initialize_with_params(int argc, char *argv[]);
double error(int ini, int fin, Tree tree);
void OneTree();
void GetNames();
void GetData(String Extension);
void SaveTree(Tree T, String Extension);
