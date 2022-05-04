#include "Func.h"
int Input (double *Matrix,char *file_name,int n, double* Norm, int p, int number)
{
        FILE *In;
	*Norm=0;
	double NormShet, prob, eps=1;
	double* MultInd;

	while(eps>eps/2){
		eps=eps/2;
	}
	
	eps=eps*100;

        In = fopen (file_name,"r");
        if (In==NULL)
        {
                fclose(In);
                return -1;
        } 
	*Norm=0;
        MultInd=Matrix;
        for(int j=0; j<n; j++){
        	if(fscanf(In,"%lf",&MultInd[j])!=1)
                {
                	fclose(In);
                        return -2;
                }
                *Norm+=fabs(MultInd[j]);
        }
        

        for (int i = 1; i<n ;i++)
        {
		NormShet=0;
		MultInd=Matrix+i/p*n;
		for(int j=0; j<n; j++){
			if(i-p*(i/p)==number){
                		if(fscanf(In,"%lf",&MultInd[j])!=1)
                		{
                        		fclose(In);
                        		return -2;
                		}
				NormShet+=fabs(MultInd[j]);
			}else{
				if(fscanf(In,"%lf",&prob)!=1){
					fclose(In);
					return-2;
				}	
			}
		}
		if((fabs(NormShet-*Norm)>eps)&&(i-p*(i/p)==number)){*Norm=NormShet;}
        }

        fclose(In);
        return 0;
}  // чтение файла

double generate (int k, int n, int i, int j)
{

                if(k==1){
                	return (double)(n+1+(-1*abs(i-j)-abs(i+1)-abs(j+1))/2);
                }
                if(k==2){
			if(i>j){
				return i*1.0;
			}else{
				return j*1.0;
                        }
                }
                if(k==3){
                	return abs(i-j)*1.0;
                }
                if(k==4){
                        return (double)(1.0/(i+j+1));
                }
		return 0;
}
