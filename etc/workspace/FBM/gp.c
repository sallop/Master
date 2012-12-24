#include <stdio.h>

int main()
{
  int i;
  FILE *gp[]={
    popen("gnuplot","w"),
    popen("gnuplot","w"),
    popen("gnuplot","w"),
  };
#define GPBEGIN(n) fprintf(gp[(n)],
#define GPEND(n)   )

#define GPCMD (n, cmd ) fprintf(gp[(n)], cmd )
#define GPDATA(n, data) fprintf(gp[(n)], data)
  GPCMD (0,"plot '-' with line\n");
  GPDATA(0,"%d %d",0,1);
  GPDATA(0,"%d %d",1,1);
  GPDATA(0,"%d %d",2,2);
  GPDATA(0,"%d %d",3,3);
  
  
  GPBEGIN(1)
    GPCMD("plot '-' with lines\n")
    GPDATA("%d %d\n",0,4)
    GPDATA("%d %d\n",1,3)
    GPDATA("%d %d\n",2,2)
    GPDATA("%d %d\n",3,1)
    GPEND(1);
  
  for (i=0; i<sizeof(gp)/sizeof(gp[0]); i++)
    fclose(gp[i]);
  return 0;
}
