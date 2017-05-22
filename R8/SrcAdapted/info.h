/*************************************************************************/
/*									 */
/*	Calculate information, information gain, and print dists	 */
/*	--------------------------------------------------------	 */
/*									 */
/*************************************************************************/


#include "buildex.i"

namespace c45
{
float Worth(float ThisInfo, float ThisGain, float MinGain);
void    ResetFreq(DiscrValue MaxVal);
float ComputeGain(float BaseInfo, float UnknFrac, DiscrValue MaxVal, ItemCount TotalItems);
float TotalInfo(ItemCount V[], DiscrValue MinVal, DiscrValue MaxVal);
void    PrintDistribution(Attribute Att, DiscrValue MaxVal, Boolean ShowNames);

}//namespace c45
