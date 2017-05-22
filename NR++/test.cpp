#include "geom_utils.h"
#include <stdlib.h>
#include <stdio.h>

int main(int argc,char**argv)
{
	double **pl = new double*[200];
srand(21);

	FILE*f=fopen("kk", "rt");
	for(int i=0;i<25;i++) pl[i] = new double[2];
	for(int i=0;i<25;i++) {
		fscanf(f, "%lf %lf", &pl[i][0], &pl[i][1]);
	}
fclose(f);	
//	pl = {{0,0},{1,1},{2,2},{3,3},{1,2},{2,1}};
/*pl[0][0]=141; pl[0][1]=72;
pl[1][0]=194; pl[1][1]=69;
pl[2][0]=195; pl[2][1]=106;
pl[3][0]=126; pl[3][1]=103;
pl[4][0]=56; pl[4][1]=160;
*/
	convexhull(pl, 25, 2);
/*	convexhull(pl, 4, 2);
	convexhull(pl, 2+rand()%19, 2);
	convexhull(pl, 2+rand()%19, 2);
	convexhull(pl, 2+rand()%19, 2);
	convexhull(pl, 5, 2);
	convexhull(pl, 4, 2);
	convexhull(pl, 2+rand()%19, 2);
	printf("=========================================================\n");
	convexhull(pl, 5, 2);
	convexhull(pl, 5, 2);
	convexhull(pl, 5, 2);
*/
	return 0;
}

