/*************************************************************************/
/*									 */
/*	Evaluation of a test on a discrete valued attribute		 */
/*      ---------------------------------------------------		 */
/*									 */
/*************************************************************************/


#include "buildex.i"

namespace c45
{

void    EvalDiscreteAtt(Attribute Att, ItemNo Fp, ItemNo Lp, ItemCount Items);
void    ComputeFrequencies(Attribute Att, ItemNo Fp, ItemNo Lp);
float DiscrKnownBaseInfo(ItemCount KnownItems, DiscrValue MaxVal);
void    DiscreteTest(Tree Node, Attribute Att);

}//namespace c45
