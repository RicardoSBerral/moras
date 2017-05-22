/*************************************************************************/
/*                                                                	 */
/*	Evaluation of a test on a continuous valued attribute	  	 */
/*	-----------------------------------------------------	  	 */
/*								  	 */
/*************************************************************************/


#include "buildex.i"

namespace c45
{

void    EvalContinuousAtt(Attribute Att, ItemNo Fp, ItemNo Lp);
void    ContinTest(Tree Node, Attribute Att);
float GreatestValueBelow(Attribute Att, float t);

}//namespace c45
