#include <stdio.h>

int main(int argc, char**argv)
{

FILE*f;
char *c, *d;

//system("ls -S1F > kk");
if (argc<2) return 1;
f = fopen (argv[1], "rt");
if (!f) return 1;
c=(char*)malloc(100);
d=(char*)malloc(101);

while(1==fscanf(f, "%s", c)){
	int lo;
	lo = strlen(c);
	if (c[lo-1]=='/' || c[lo-1]=='@') continue;
	if (c[lo-1]=='*') c[lo-1] = '\0';
	printf("------------------------------\n");
	FILE *f1=fopen(c, "r");
	if (!f1) {strcpy(d, "No se abre");}
	else {
		int i = fread(d, 1, 100, f1);
		d[i]='\0';
		fclose(f1);
	}
	printf("%s: %s\n", c, d);
	//scanf("%s\n", d);
	
}
fclose(f);

return 0;
}
