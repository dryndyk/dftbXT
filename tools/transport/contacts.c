#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SQRT3 1.732050808
#define SQRT6 2.449489743
#define LEN 250

int main(int argc,char *argv[]){

FILE *FPin;
char **atm,*filein, *tmpstr;
int i,j,k,*ind,nsp,natm,adv,npl1,npl2,istart1,iend1,istart2,iend2;
double *x,*y,*z,latt_v1,latt_v2;

if (argc<10)
{ 
  printf("Usage:  contacts filein.gen start1 end1 npl1 latt_vect1 start2 end2 npl2 latt_vect2 \n");
  printf("      npl: layer repetitions (if structure has only 1 contact PL then npl=2)\n");
  printf("latt_vect: scalar assuming transport along z-direction \n");
  printf("           mind sign reflecting contact direction \n");
  return;
}

nsp=20;
filein=malloc(60*sizeof(char));
tmpstr=malloc(LEN*sizeof(char));

strcpy(filein,argv[1]);
istart1=atoi(argv[2]);
iend1=atoi(argv[3]);
npl1=atoi(argv[4]);
latt_v1=atof(argv[5]);
istart2=atoi(argv[6]);
iend2=atoi(argv[7]);
npl2=atoi(argv[8]);
latt_v2=atof(argv[9]);


atm=malloc(nsp*sizeof(char *));
for (i=0;i<nsp;i++)
  {
    atm[i]=malloc(3*sizeof(char));
    //atm[i]=NULL;
  }
 
FPin=fopen(filein,"r");
if (!FPin)
  { 
    printf("File error or doesn't exist\n");
    free(filein);
    return;
  }

fgets(tmpstr,LEN,FPin);
sscanf(tmpstr," %d",&natm);

x = malloc((natm+1)*sizeof(double));
y = malloc((natm+1)*sizeof(double));
z = malloc((natm+1)*sizeof(double));
ind = malloc((natm+1)*sizeof(int));
	


fgets(tmpstr,LEN,FPin);
adv = strcspn(tmpstr,"\n");
tmpstr[adv]=' ';

adv=0;
i=0;
while(sscanf(tmpstr+adv," %2c",atm[i]) > 0)
 {
   adv=(int) ( strstr(tmpstr,atm[i])- tmpstr) +2;
   i++;
 }

// Start to write the new gen
printf("%d C \n",natm+(npl1-1)*(iend1-istart1+1)+(npl2-1)*(iend2-istart2+1));
nsp=i;
for(i=1; i<=nsp;i++){printf("%s ",atm[i-1]);}
printf("\n");

for(i=1;i<=natm;i++)
  {
    fgets(tmpstr,LEN,FPin);
    sscanf(tmpstr," %d %d %lg %lg %lg",&j,ind+i,x+i,y+i,z+i);
    //if(ind>nsp){printf("Error in file: at row %d, type %d > nsp %d\n",j,ind,nsp);return;}
  }


for(i=1;i<=natm;i++)
{
  if( (i>=istart1 && i<=iend1) || (i>=istart2 && i<=iend2)){continue;}

  printf("%d %d %f %f %f\n",i,ind[i],x[i],y[i],z[i]); 
}

for(k=0;k<npl1;k++)
{
  for(i=istart1;i<=iend1;i++)
  { 
    printf("%d %d %f %f %f\n",i,ind[i],x[i],y[i],z[i]+k*latt_v1);
  }
}

for(k=0;k<npl2;k++)
{
  for(i=istart2;i<=iend2;i++)
  { 
    printf("%d %d %f %f %f\n",i,ind[i],x[i],y[i],z[i]+k*latt_v2);
  }
}


fclose(FPin);

free(filein);
free(tmpstr);
for(i=0;i<nsp;i++){ free(atm[i]); }
free(atm);

}
	 	
    /* printf("HETATM%5d%3s%12d     %6.3f  %6.3f  %6.3f\n",i,atm[ind],i,x,y,z);*/	 	
