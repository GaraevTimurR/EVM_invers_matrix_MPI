#include "Func.h"

double NormSol(double * Matrix, double* InvMatrix, int n, double* X)
{
        double shet, s=0;
	//double* f;
        int i, j,  k, p, number;
	double Norm=0, eps=1;

	while(eps>eps/2){
		eps=eps/2;
	}
	eps=eps*1000;

	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &number);

        for(i=0; i<n; i++){
		shet=0;
		for(j=0; j<n; j++){
			X[j]=Matrix[(i/p)*n+j];
		}
		MPI_Bcast(X, n, MPI_DOUBLE, i-(i/p)*p, MPI_COMM_WORLD);
		//for(j=0; j<n; j++){
		//	printf("%f", X[j]);
		//}
		//printf("\n");
                for(j=number; j<n; j+=p){
                        s=0;
                        if(i==j){
				for(k=0; k<n; k++){
                                	s+=X[k]*InvMatrix[(j/p)*n+k];
                        	}
                                s=s-1;
                        }else{
				for(k=0; k<n; k++){
                                	s+=X[k]*InvMatrix[(j/p)*n+k];
                        	}
			}
			shet+=abs(s);
                }
		Norm+=shet;
		//printf("%f\n", shet);
		//if(abs(Norm-shet)<eps){Norm=shet;}
        }
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allreduce(&Norm, &s, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	//printf("%e\n", s);
        return s;
}
