#include "Func.h"
void POut(double *Matrix,int l,int n,int m, double* X)
{
        int c,d,number,p;
        MPI_Comm_size(MPI_COMM_WORLD, &p);
        MPI_Comm_rank(MPI_COMM_WORLD, &number);
        c = l+m-(abs(l-m)+l+m)/2;
        d = n+m-(abs(n-m)+n+m)/2;
        MPI_Status status;
        //double shet=1, shet2=0;
        double eps = 1.0;
        double* MultInd;
        while (1.0 + eps/ 2.0 > 1.0)
        {
	        eps /= 2.0;
        }
        eps=eps*1000;
        for(int i = 0;i<c;i++)
        {
	        MultInd=Matrix+(i/p)*n;
                for(int j=0; j<n; j++){
	                X[j]=MultInd[j];
                }
                if(i-(i/p)*p!=0){
        	        //printf("HI, number %d , i-p(i/p) %d 1 \n", number, i-p*(i/p));
	                if(number==0){
	        	        MPI_Recv(X, n, MPI_DOUBLE, i-(i/p)*p,  0, MPI_COMM_WORLD, &status);
                	}else{
	                	if(number == i-(i/p)*p){
	                        	MPI_Send(X, n, MPI_DOUBLE, 0,  0, MPI_COMM_WORLD);
	                        }
        	        }
                }
                if(number==0){
	                for(int j = 0;j<d;j++){
        	                printf("%e ",X[j]);
                        }
                        printf("\n");
                }
                MPI_Barrier(MPI_COMM_WORLD);
	}
}
 
void POut2(double *InvMatrix,int l,int n,int m)
{
	int c,d,number,p;
        MPI_Comm_size(MPI_COMM_WORLD, &p);
        MPI_Comm_rank(MPI_COMM_WORLD, &number);
        c = l+m-(abs(l-m)+l+m)/2;
        d = n+m-(abs(n-m)+n+m)/2;
        MPI_Status status;
	double a;

	for(int i=0; i<c; i++){
		for(int j=0; j<d; j++){
			MPI_Barrier(MPI_COMM_WORLD);
			a=InvMatrix[(j/p)*n+i];
			if((number==0)&&(0!=j-(j/p)*p)){
                                MPI_Recv(&a, 1, MPI_DOUBLE, j-(j/p)*p,  0, MPI_COMM_WORLD, &status);
                        }else{
                                if((number == j-(j/p)*p)&&(0!=j-(j/p)*p)){
                                        MPI_Send(&InvMatrix[(j/p)*n+i], 1, MPI_DOUBLE, 0,  0, MPI_COMM_WORLD);
                                }
                        }
			if(0==number){
				printf("%e ", a);
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if(number==0){printf("\n");}
		MPI_Barrier(MPI_COMM_WORLD);
	}
}






