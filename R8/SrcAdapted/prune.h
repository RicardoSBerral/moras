/*************************************************************************/
/*									 */
/*	Prune a decision tree and predict its error rate		 */
/*	------------------------------------------------		 */
/*									 */
/*************************************************************************/


#include "defns.i"
#include "types.i"
#include "extern.i"

#define	LocalVerbosity(x)	if (Sh >= 0 && VERBOSITY >= x)
#define	Intab(x)		Indent(x, "| ")

namespace c45
{

Boolean Prune(Tree T);
float EstimateErrors(Tree T, ItemNo Fp, ItemNo Lp, short Sh, Boolean UpdateTree);
void    CheckPossibleValues(Tree T);

}//namespace c45
