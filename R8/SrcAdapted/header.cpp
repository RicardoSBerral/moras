/*************************************************************************/
/*									 */
/*	Print header for all C4.5 programs				 */
/*	----------------------------------				 */
/*									 */
/*************************************************************************/

#include "header.h"
#include <stdio.h>
#include <time.h>
#include <string.h>

namespace c45
{

void    PrintHeader(char *Title)
/*  -----------  
    char *Title;*/
{
    char /**ctime(),*/ TitleLine[80];
    long clock/*, time()*/;
    short Underline;

    clock = time(0);
    sprintf(TitleLine, "C4.5 [release %s] %s", RELEASE, Title);
    printf("\n%s\t%s", TitleLine, ctime(&clock));

    Underline = strlen(TitleLine);
    while ( Underline-- ) putchar('-');
    putchar('\n');
}


}//namespace c45
