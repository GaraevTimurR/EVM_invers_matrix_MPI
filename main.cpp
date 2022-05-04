#include "Func.h"

long double get_time()
{
	    struct timeval t;
	        gettimeofday(&t, 0);
		    return t.tv_sec + t.tv_usec/1000000.0;
}
double mach_eps(void)
{
	double eps = 1.0;

	while (1.0 + eps / 2.0 > 1.0)
	{
		eps /= 2.0;
	}
	return eps;
} // машинное эпсилон

int main(int argc, char* argv[])
{
	int n; 
	int m; 
	int prob; 
	int ans; 
	int k;
	double* Matrix;
	double* InvMatrix;
	int* swap;
	int number;
	double* X;
	double* Y;
	double* MultInd;
	double ForNorm,	 Norm=0, eps=mach_eps();
	double time=0, time1;
	int p;

	MPI_Init(&argc, &argv);

	n = atoi(argv[1]);
        m = atoi(argv[2]);
        prob=0;
        ans=0;
        k=atoi(argv[3]);
	
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	Matrix = (double*)malloc(n * (n/p+1) * (sizeof(double)));
	InvMatrix = (double*)malloc(n * (n/p+1) * (sizeof(double)));
	swap=(int*)malloc(n *  (sizeof(int)));
	X=(double*)malloc(n*(sizeof(double)));
	Y=(double*)malloc(n*(sizeof(double)));

	MPI_Comm_rank(MPI_COMM_WORLD, &number);
	//printf("p = %d\n", p);

	if(argc<4){
		printf("Can not initilaze mutex\n");
		free(X);
		free(Y);
		free(swap);
		free(InvMatrix);
		free(Matrix);
		return -1;
	}


	if (n == 0 || n < m) {
		printf("Incorrect input\n");
		free(X);
		free(Y);
		free(swap);
		free(InvMatrix);
		free(Matrix);
		return -1;
	}

	if ((k < 0) || (k > 4)) {
		printf("Incorrect input\n");
		free(X);
		free(Y);
		free(swap);
		free(InvMatrix);
		free(Matrix);
		return -1;
	}

	if (Matrix == NULL)
	{
		printf("No memory allocated\n");
		free(X);
		free(Y);
		free(swap);
		free(InvMatrix);
		free(Matrix);
		return -1;
	}
	if ((argc == 5) && (k == 0)) {
		
		prob = Input(Matrix, argv[4], n, &Norm, p, number);
		if (prob == -1) {
			printf("File not found\n");
			free(X);
			free(Y);
			free(swap);
			free(InvMatrix);
			free(Matrix);
			return -1;
		}
		if (prob == -2)
		{
			printf("Insufficient data \n");
			free(X);
			free(Y);
			free(swap);
			free(InvMatrix);
			free(Matrix);
			return -1;
		}
	}
	if((argc!=5)&&(k==0))
        {
                printf("File don't exist\n");
		free(Matrix);
		free(X);
		free(Y);
		free(swap);
		free(InvMatrix);
                return -1;
        }
        if(k>0){
                for(int i=number; i<n ; i+=p){
                        ForNorm=0;
			MultInd=Matrix+(i/p)*n;
                        for(int j=0; j<n; j++){
                                MultInd[j]=generate(k, n, i, j);
                                ForNorm+=fabs(generate(k, n, i, j));
                        }
                        if(ForNorm>Norm){Norm=ForNorm;}
                }
        }

	//      }
	//}


	if (InvMatrix == NULL)
	{
		printf("No memory allocated\n");
		free(Matrix);
		free(InvMatrix);
		free(X);
		free(Y);
		free(swap);
		return -1;
	}
	for(int i=number; i<n ; i+=p){
	        MultInd=InvMatrix+(i/p)*n;
	        for(int j=0; j<n; j++){
			if(i!=j){
		        	MultInd[j]=0;
			}else{
				MultInd[j]=1;
			}
		}
	}

	//инициализация и заполнение

	//решение
	if(number==0){
		printf("Matrix : \n");
	}
	POut(Matrix, n, n, m, X);
	if(number==0){
		printf("\n");
	}
	MPI_Allreduce(&Norm, &ForNorm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	if(ForNorm/eps<eps){
		printf("Degenerate matrix");
                free(Matrix);
                free(InvMatrix);
                free(X);
		free(Y);
                free(swap);
		return -1;
	};
	
	Norm=1/ForNorm;

	if(number==0){
		printf("Norm %e\n", ForNorm);
	}
	
	//Norm = 1 / Norm;

	for (int i = number; i < n; i+=p)
	{
		MultInd = Matrix + (i/p)* n;
		for (int j = 0; j < n; j++)
		{
			MultInd[j] = MultInd[j] * Norm;
		}
	}
	
	
	if (InvMatrix == NULL)
	{
                printf("No memory allocated\n");
                free(Matrix);
                free(InvMatrix);
		free(X);
		free(Y);
		free(swap);
                return -1;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	time=MPI_Wtime();
	ans = Gauss(Matrix, InvMatrix, n, swap, X, Y);
	time=MPI_Wtime()-time;

	if (ans == -2) {
		printf("Degenerate matrix");
	        free(Matrix);
	        free(InvMatrix);
		free(X);
		free(Y);
		free(swap);
		return -1;
	}
	
	time1=time;

	MPI_Allreduce(&time, &time1, 1, MPI_DOUBLE, MPI_MAX,  MPI_COMM_WORLD);
	if(number==0){
		printf("time; %f \n", time);
	}
	for (int i = number; i < n; i+=p)
	{
		MultInd = InvMatrix + (i/p) * n;
		for (int j = 0; j < n; j++)
		{
			MultInd[j] = MultInd[j] * Norm;
		}
	}
	//Gauss(Matrix, InvMatrix, n, swap, X);
	if(number==0){
		printf("Inverse Matrix: \n");
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(number==0){
		 printf("\n");
	}
	POut2(InvMatrix, n, n, m);

	MPI_Barrier(MPI_COMM_WORLD);

	if ((argc == 5) && (k == 0)) {
		prob = Input(Matrix, argv[5], n, &Norm, p, number);
		if (prob == -1) {
			printf("File not found\n");
		        free(Matrix);
			free(swap);
			free(X);
			free(Y);
		        free(InvMatrix);
			return -1;
		}
		if (prob == -2)
		{
			printf("Insufficient data \n");
		        free(Matrix);
		        free(swap);
			free(X);
			free(Y);
		        free(InvMatrix);
			return -1;
		}
	}
	
	if(k>0){
                for(int i=number; i<n ; i+=p){
			MultInd=Matrix+(i/p)*n;
                        for(int j=0; j<n; j++){
                                MultInd[j]=generate(k, n, i, j);
                        }
                }
        }
	//POut(Matrix, n, n, m, X);
	time1=MPI_Wtime();
	MPI_Barrier(MPI_COMM_WORLD);
	eps=NormSol(Matrix, InvMatrix, n, X);
	MPI_Barrier(MPI_COMM_WORLD);
	time1=MPI_Wtime()-time1;
	time=time1+time;
	time1=time;

	printf("%f\n", eps);

	MPI_Reduce(&time, &time1, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

	if(number==0){
		printf("%s : residual = %e elapsed = %f s = %d n = %d m = %d p = %d\n", argv[0], eps, time1, k, n, m, p);
	}
	free(Matrix);
	free(swap);
	free(X);
	free(Y);
	free(InvMatrix);
	MPI_Finalize();
	return 0;
}

