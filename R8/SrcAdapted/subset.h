/*************************************************************************/
/*									 */
/*      Evaluation of the subsetting of a discrete attribute		 */
/*      ----------------------------------------------------		 */
/*									 */
/*************************************************************************/


#include "buildex.i"

namespace c45
{

void    EvalSubset(Attribute Att, ItemNo Fp, ItemNo Lp, ItemCount Items);
void    Combine(DiscrValue x, DiscrValue y, DiscrValue Last);
void    Uncombine(DiscrValue x, DiscrValue y);
void    PrintSubset(Attribute Att, Set Ss);
void    SubsetTest(Tree Node, Attribute Att);

}//namespace c45
