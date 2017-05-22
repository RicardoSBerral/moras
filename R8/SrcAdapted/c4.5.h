/*************************************************************************/
/*									 */
/*	Main routine, c4.5						 */
/*	------------------						 */
/*									 */
/*************************************************************************/


#include "defns.i"
#include "types.i"


    /*  External data, described in extern.i  */
namespace c45
{
void SetGlobalOpt(int opt, char *arg);
void InitParams(int Argc, char *Argv[]);
int    run(int Argc, char *Argv[]);
}
