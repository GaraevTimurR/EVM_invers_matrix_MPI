#include "Func.h"

int Gauss(double* Matrix, double* InvMatrix, int n, int* swap, double* X, double* Y){
        double prob;
        double Del, BigDel;
	double max, eps=mach_eps();
	int maxJ, j;
        int k, ans=0;
        double* MultInd;
	//double* MultInd2;
//	double* MultInd3;
	int p, number;
	//MPI_Status Status;

	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &number);
		
	for(int i=0; i<n; i++){swap[i]=i;}
        for(int i=0; i<n; i++){
                MultInd= Matrix+(i/p)*n;
		max=MultInd[i];
		maxJ=i;
		if(i-(i/p)*p==number){
               		for(int j=i; j<n; j++){
                       		if(fabs(max)<fabs(MultInd[j])){
                               		max=MultInd[j];
                               		maxJ=j;
                       		}
		
       	        	}
		}
		
               	if((fabs(max)<eps)&&(i-(i/p)*p==number)){
			ans=-2;
			MPI_Bcast(&ans, 1, MPI_INT, i-(i/p)*p, MPI_COMM_WORLD);
                       	return -2;
               	}else{
			if(i-(i/p)*p==number){
				MPI_Bcast(&ans, 1, MPI_INT, i-(i/p)*p, MPI_COMM_WORLD);
				MPI_Bcast(&maxJ, 1, MPI_INT, i-(i/p)*p, MPI_COMM_WORLD);	
			}else{
				MPI_Bcast(&ans, 1, MPI_INT, i-(i/p)*p, MPI_COMM_WORLD);
				if(ans!=-2){
					MPI_Bcast(&maxJ, 1, MPI_INT, i-(i/p)*p, MPI_COMM_WORLD);
				}else{
					return -2;
				}
			}
		}


             	k=swap[i];
                swap[i]=swap[maxJ];
                swap[maxJ]=k;
		

		for(k=number; k<n; k+=p){
                        prob=Matrix[(k/p)*n+i];
                        Matrix[(k/p)*n+i]=Matrix[(k/p)*n+maxJ];
                        Matrix[(k/p)*n+maxJ]=prob;
		}

 	       	//меняем столбцы	
	        //меняем столбцы

		MultInd=Matrix+(i/p)*n;

		for(k=0; k<n; k++){
			X[k]=MultInd[k];
			Y[k]=Matrix[(k/p)*n+i];
		}

		MPI_Bcast(X, n, MPI_DOUBLE, i-(i/p)*p, MPI_COMM_WORLD);
		for(int k=0; k<n; k++){
			MPI_Bcast(&Y[k], 1, MPI_DOUBLE, k-(k/p)*p, MPI_COMM_WORLD);
		}

		BigDel=X[i];
	        // диагонализируем
		MPI_Barrier(MPI_COMM_WORLD);

		for(k=0; k<n; k+=p){
			for(int t=i+1; t<n; t++){
				InvMatrix[(k/p)*n+t]=InvMatrix[(k/p)*n+t]-InvMatrix[(k/p)*n+i]*Y[t]/Y[i];
			}
		}

		if((i/p)*p+number<i){k=(i/p)*p+p+ number;}else{k=(i/p)*p+ number;}

        	for(; k<n; k+=p){
	 	       if((k/p)*p+number!=i){
        		       Del = Matrix[n*(k/p)+i]/BigDel;
                	       for(int t=0; t<n; t++){
                		      Matrix[(k/p)*n+t]=Matrix[(k/p)*n+t]-Del*X[t];
				}
               		}
              	}

	        MPI_Barrier(MPI_COMM_WORLD);
		// диагонализируе
	}


	for(int k=0; k<n; k++){Y[k]=0;}

	for(int k=n-1; k>-1; k--){
		MultInd=Matrix+(k/p)*n;
		if((k/p)*p+number==k){
			for(int i=0; i<n; i++){
				X[i]=MultInd[i];
			}
		}
		MPI_Bcast(X, n, MPI_DOUBLE, k-(k/p)*p, MPI_COMM_WORLD);

		BigDel=X[k];
		for(int i=number; i<n; i+=p){
			Del=0;
			for(j=k+1; j<n; j++){
				Del+=X[j]*InvMatrix[(i/p)*n+j];
			}
			InvMatrix[(i/p)*n+k]=(InvMatrix[(i/p)*n+k]-Del)/BigDel;
			
		}

		MPI_Barrier(MPI_COMM_WORLD);


        }

	for(int i=0; i<n; i++){
		while(swap[i]!=i){
			int j=swap[i];	
			for(k=number; k<n; k+=p){
				Del=InvMatrix[(k/p)*n+i];
				InvMatrix[(k/p)*n+i]=InvMatrix[(k/p)*n+j];
				InvMatrix[(k/p)*n+j]=Del;
			}

			swap[i]=swap[j];
			swap[j]=j;
			MPI_Barrier(MPI_COMM_WORLD);
		}
	}
        return 0;
}

