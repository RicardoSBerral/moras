/*************************************************************************/
/*									 */
/*	Routines for displaying, building, saving and restoring trees	 */
/*	-------------------------------------------------------------	 */
/*									 */
/*************************************************************************/


#include "defns.i"
#include "types.i"
#include "extern.i"

#define	Tab		"|   "
#define	TabSize		4
#define	Width		80	/* approx max width of printed trees */

namespace c45
{

	/*  If lines look like getting too long while a tree is being
	    printed, subtrees are broken off and printed separately after
	    the main tree is finished	 */

void    PrintTree(Tree T);
void    Show(Tree T, short Sh);
void    ShowBranch(short Sh, Tree T, DiscrValue v);
short MaxLine(Tree St);
void    Indent(short Sh, char *Mark);
void    SaveTree(Tree T, String Extension);
void    OutTree(Tree T);
Tree GetTree(String Extension);
Tree InTree(void);
void    StreamOut(String s, int n);
void    StreamIn(String s, int n);
void    ReleaseTree(Tree Node);
Tree Leaf(ItemCount *ClassFreq, ClassNo NodeClass, ItemCount Cases, ItemCount Errors);
void    Sprout(Tree Node, DiscrValue Branches);
int    TreeSize(Tree Node);
Tree CopyTree(Tree T);
void    SaveDiscreteNames();
void    RecoverDiscreteNames();

}//namespace c45
