#include <stdio.h>
#include "hull.h"
#define NUMS 6
double pts[NUMS][2] = {{0,0},{1,1},{2,2},{3,3},{1,2},{2,1}};
FILE *DFILE;

long su_numero(site p)
{
  int j = (p-pts[0])/2;
//  printf("(%g,%g) is %d\n", p[0], p[1], j);
  return j;
}

void *primir(simplex *s, void *p) 
{
  int i; 
  printf("K: ");
  for (i=0;i<cdim;i++)
    printf("%d ", su_numero(s->neigh[i].vert));
  printf("\n");

  return 0;
}
site err_sitio_numero(int j)
{
  if(j>=NUMS || j<0) return 0;
  //printf("(%g,%g) IS %d\n", pts[j][0], pts[j][1], j);
  return pts[j];
}
site err_siguiente(void)
{
  static int k=0;
  return err_sitio_numero(k++);
}
int main()
{
  DFILE=stderr;

  simplex *root = build_convex_hull(err_siguiente, su_numero, 2, 0);
  
  visit_hull(root, primir);
	  
  return 0;
}
