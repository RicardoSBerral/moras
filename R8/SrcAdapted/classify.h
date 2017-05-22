/*************************************************************************/
/*                                                              	 */
/*  Determine the class of a case description from a decision tree	 */
/*  --------------------------------------------------------------	 */
/*                                                              	 */
/*************************************************************************/
#ifndef __C45_CLASSIFY_H
#define __C45_CLASSIFY_H


#include "defns.i"
#include "types.i"
#include "extern.i"


namespace c45
{
extern float	*ClassSum;		/* ClassSum[c] = total weight of class c */
extern Tree	FinalNode;		

ClassNo Category(Description CaseDesc, Tree DecisionTree);
void    Classify(Description CaseDesc, Tree T, float Weight);

}//namespace c45

#endif
