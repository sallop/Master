#include "FBM.hh"

#define BUFSZ 512
int main(int argc, char *argv[])
{
  char idir[32]="dir-dat",prefix[32]="L";
  int result;
  while ((result = getopt(argc,argv,"d:p:")) != -1){
    switch (result){
    case 'd': strcpy(idir  , optarg); break;
    case 'p': strcpy(prefix, optarg); break;
    case '?': std::cout << "Not implimented " << optarg << std::endl;
    }
  }

  Layer L;
  Probability gibbs_P;

  struct {std::ofstream ifs; char id[BUFSZ]} *idx,ftbl[] = {
    {std::ifstream, "-mat.dat"    },
    {std::ifstream, "-t_gibbs.dat"},
    {std::ifstream, "-cfg.dat"    },
    {NULL,NULL}
  };

  char fname[64];
  for (idx=&ftbl[0]; idx->id != NULL; idx++){
    sprintf(fname, "%s/%s%s",idir,prefix,idx->id);
    printf("%s",fname);
    //    idx->ifs.open(fname);
  }

  

//   char line[BUFSZ];
//   while (ftbl[0].os.getline(line,BUFSZ)){
//     int t;
//     char mat[BUFSZ];
//     sscanf(line, "%d %s", &t, &mat);
//     printf("time=%d,matrix=%s",t,mat);
//  }

  
  for (idx=&ftbl[0]; idx->id != NULL; idx++)
    idx->os.close();
  return 0;
}
