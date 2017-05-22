/*************************************************************************/
/*									 */
/*	Get names of classes, attributes and attribute values		 */
/*	-----------------------------------------------------		 */
/*									 */
/*************************************************************************/

namespace c45 {

void    ReleaseNames(void);
Boolean ReadName(FILE *f, String s);
void    GetNames();
int Which(String Val, String List[], short First, short Last);
String CopyString(String x);
void    Error(short n, String s1, String s2);

}//namespace c45
