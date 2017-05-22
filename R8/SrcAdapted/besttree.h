/*************************************************************************/
/*									 */
/*	Routines to manage tree growth, pruning and evaluation		 */
/*	------------------------------------------------------		 */
/*									 */
/*************************************************************************/


#include "defns.i"
#include "types.i"
#include "extern.i"


#include "buildex.i"

namespace c45
{

void    OneTree(void);
short BestTree(void);
void    FormTarget(ItemNo Size);
void    FormInitialWindow(void);
void    Shuffle(void);
Tree Iterate(ItemNo Window, ItemNo IncExceptions);
void    Evaluate(Boolean CMInfo, short Saved);

}//namespace c45
