/*************************************************************************/
/*									 */
/*	Get case descriptions from data file				 */
/*	--------------------------------------				 */
/*									 */
/*************************************************************************/


#include "defns.i"
#include "types.i"
#include "extern.i"

#define Inc 2048

namespace c45
{
void    ReleaseData(void);
void    GetData(String Extension);
Description GetDescription(FILE *Df);
}//namespace c45
