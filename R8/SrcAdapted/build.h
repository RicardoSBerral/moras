/*************************************************************************/
/*								 	 */
/*    Central tree-forming algorithm incorporating all criteria  	 */
/*    ---------------------------------------------------------	 	 */
/*								 	 */
/*************************************************************************/


#include "defns.i"
#include "types.i"
#include "extern.i"

namespace c45
{

void      FreeTreeData(void);
void      InitialiseTreeData(void);
void      InitialiseWeights(void);
Tree      FormTree(ItemNo Fp, ItemNo Lp);
ItemNo    Group(DiscrValue V, ItemNo Fp, ItemNo Lp, Tree TestNode);
ItemCount CountItems(ItemNo Fp, ItemNo Lp);
void      Swap(ItemNo a, ItemNo b);

}//namespace c45
