/* computing the possion equation */
/* Hanieh Soleimani */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void * _smalloc(size_t  count,  size_t  size)
{
        void * new = malloc (count * size);
        if (new == NULL)
        {
                printf("ERROR during memory allocation!\n");
                exit (7);
        }
        return  new;
}

#define  mem_free(var)                  do { free(var); var = NULL; } while (0)
#define  mem_alloc(var, cnt, typ) var = (typ *) _smalloc  ((cnt), sizeof(typ))


/* DGESV prototype */
extern void dgesv_( int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info );

/* DGEMV prototype */
extern void dgmv(char* trans, int* m, int* n, double* alpha, double* a, int* lda, double* x, int* incx, double* beta, double* y, int* incy);

/* Auxiliary routines prototypes */
extern void print_matrix( char* desc, int m, int n, double* a, int lda );
extern void print_int_vector( char* desc, int n, int* a );
extern void print_double_vector( char* desc, int n, double* a );
/* compute Kronecker Product, source from http://www.mymathlib.com/c_source/matrices/arithmetic/kronecker_product.c */
extern void Kronecker_Product(double *C, double *A, int nrows, int ncols, double *B, int mrows, int mcols);
/* Equation used */
extern double f(double x, double y);
extern double g(double x, double y);

int main() {
	/* Locals */
	int   _n[3] = {50, 100, 200};
    int _rhs[3] = {1,   1,    1}; 
    int nrhs;
    int n;                    //      n: number of approximation points 
	int n_p;				  //	n_p: dimension of the poisson matrix
    int ldpoisson, ldb;             //    ldx: leading dimensions of a and b
    int info;                 //   info: output parameter for SCALAPACK functions
	
	double min_bound;		//	Store the minimum of the square that we want to discretize
	double max_bound;		// Store the maximum of the square that we want to dicretize
	double delta;			// Step used to compute cooridnates of the points
	
	/* Local arrays */
   	int    *ipiv 	= NULL;   //   ipiv: permutation applied on the columns of a
	
	double *poisson 	= NULL;   //      a: matrix
	double *poisson_copy	= NULL;   // a_copy: copy of a
	
	double* solution 	= NULL; 	// solution: store the analytic solution, with boundaries
	double* solution_copy	 = NULL;	// solution_copy: copy of the solution 
	
	double *J 	= NULL;	//	matrix in the diagonal of A
	double *I_ 	= NULL;	//	-identity matrix for A
	double *I 	= NULL;	//	identity matrix for A
	
	double *b	= NULL;   //      b: Right-Hand-Sides
	double *b_copy	= NULL;   // b_copy: copy of b
	
	double *u 	= NULL;   //      u: solution of the system
	double *u1 	= NULL;   //     u1: solution of the system
	
	double *res 	= NULL;   //    res: residual
	
	int i, j, k, ndim;
	
    for (ndim = 0; ndim < 2; ndim++) {
        n    = _n[ndim];
		n_p = n*n;
        nrhs = _rhs[ndim];
        ldpoisson  = n*n;
        ldb  = n*n;
		min_bound = 0;
		max_bound = 1;
		delta = (max_bound - min_bound) / (n + 1);
		
		/* Allocate space for matrices */
		mem_alloc(poisson, n_p * n_p, double);
		
		mem_alloc(ipiv, n_p, int);
		
		mem_alloc(J, n * n, double);
		mem_alloc(I_, n * n, double);
		mem_alloc(I, n*n, double);
		
		mem_alloc(b, n_p, double);
		mem_alloc(b_copy, n_p, double);
		
		mem_alloc(u, n_p, double);
		mem_alloc(u1, n_p, double);
		mem_alloc(res, n_p, double);
		
		mem_alloc(solution, (n+2) * (n+2), double);
		mem_alloc(solution_copy, (n+2) * (n+2), double);
		
		printf("##### Computing Poisson equation using n: %d #####\n", n);
		
		/* Create Poisson matrix */
		/* Create J */
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				if (i == j) {
					J[i * n + j] = 4;
				} else if (i == j+1 || j == i+1) {
					J[i * n + j] = -1;
				}
			}
		}
		// print_matrix("J matrix", n, n, J, n);
		
		/* Create -I matrix */
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				if (i == j) {
					I_[i * n + j] = -1;
				}
			}
		}
		// print_matrix("-I matrix", n, n, I_, n);
		
		/* Create I matrix */
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				if (i == j) {
					I[i * n + j] = 1;
				}
			}
		}
		// print_matrix("I matrix", n, n, I, n);
		
		/* Create first part of poisson matrix */
		Kronecker_Product(poisson, I, n, n, J, n, n);
		
		/* Create matrix with two diagonals */
		double* I_double;
		mem_alloc(I_double, n*n, double);
		/* Fill the matrix */
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				if (i == j+1 || j == i+1) {
					I_double[i * n + j] = 1;
				}
			}
		}
		// print_matrix("I with two diagonals matrix", n, n, I_double, n);
		
		/* Create second part of the poisson matrix */
		double* poisson_2;
		mem_alloc(poisson_2, n_p * n_p, double);
		Kronecker_Product(poisson_2, I_double, n, n, I_, n, n);
		
		/* Sum the two matrices */
		for (i = 0; i < n_p; i++) {
			for  (j = 0; j < n_p; j++) {
				poisson[i * n_p + j] += poisson_2[i * n_p + j];
				poisson[i * n_p + j] /= delta * delta;
			}
		}
		
		/* Free the temporary created matrices */
		mem_free(poisson_2);
		mem_free(I_double);
		
		// print_matrix("Poisson matrix", n_p, n_p, poisson, n_p);
		
		/* Compute the matrix with the analytic solution */
		for (i = 0; i < n+2; i++) {
			for (j = 0; j < n+2; j++) {
				solution[i * (n+2) + j] = g(i * delta, j * delta);
				solution_copy[i * (n+2) + j] = g(i * delta, j * delta);
			}
		}
		
		// print_matrix("Analytic solution", n+2, n+2, solution, n+2);
		
		/* Create the right hand side vector */
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				b[i * n + j] = -f((i + 1) * delta, (j + 1) * delta);
			}
		}
		
		/* Add boundary conditions, rember that matrix are stored colum-wise */
		/* First row */
		for (i = 0; i < n; i++) {
			b[i*n] += solution[(i+1) * (n+2)] / (delta * delta); 
		}
		/* Last row */
		for (i = 0; i < n; i++) {
			b[i*n + (n-1)] += solution[(i+1) * (n+2) + (n+1)] / (delta * delta);
		}
		/* First colum */
		for (i = 0; i < n; i++) {
			b[i] += solution[i+1] / (delta * delta);
		}
		/* Last column */
		for (i = 0; i < n; i++) {
			b[(n-1) * n + i] += solution[(n+1)*(n+2) + i + 1] / (delta * delta);
		}
		// print_double_vector("Right hand side", n_p, b);
		
		/* Copy b */
		for (i = 0; i < n_p; i++)
			b_copy[i] = b[i];
		
		/* Solve system using dgesv routine */
		printf("Computing solution using gauss\n");
		dgesv_(&n_p, &nrhs, poisson, &ldpoisson, ipiv, b, &ldb, &info);
		/* b now contains the solution to our system */
		// print_double_vector("Solution Gauss", n_p, b);
		
		/* Put our solution into the center of solution_copy */
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				solution_copy[(i+1)*(n+2) + j+1] = b[i*n + j];
			}
		}
		


        for (i = 0; i < nrhs; i++) {
            /* Auxiliary variables */
            double alpha    =  1.0;
            double zero     =  0.0;
            double beta     = -1.0;
            double residual =  0.0;
            double normb    =  0.0;
            int    one_inc  =  1;
			double rel_res;

            normb = 0.0;
			
            for (j = 0; j < n_p; j++)
            {
                res[j] = b_copy[i*n_p + j];
                normb += b_copy[i*n_p + j] * b_copy[i*n_p + j];
            }

            normb = sqrt(normb);


            dgemv_("N", &n_p, &n_p, &alpha, poisson, &ldpoisson, b, &one_inc, &beta, res, &one_inc );

            for (j = 0; j < n_p; j++)
                residual += res[j] * res[j];

            residual = sqrt(residual);

            rel_res = residual/normb;
            printf( "[%d-rhs] Norm of || A x - b || / || b  ||  %17.5e\n", i, rel_res);
    	}
		
		/* Compute solution using Jacobi iteration */
		printf("Computing solution using Jacobi\n");
        for(k = 0; k < nrhs; k++) {
            double alpha    =  1.0;
            double zero     =  0.0;
            double beta     = -1.0;
            double residual =  0.0;
            double normb    =  0.0;
            int    one_inc  =  1;

            double rel_res = 1e100;

            printf( "\n");

            /* Initial guess x is zero */
            for(i = 0; i < n_p; i++ )
                u[i]  = 0;

			/* Repeat iteration until error is below a certain value */
            while (rel_res > 1e-4) {
				/* Loop over the columns */
                for(i = 0; i < n_p; i++ ) {
                    u1[i] = b_copy[i + k*n_p];
					/* Summ the row */
                    for(j = 0; j < n_p; j++ )
                    {
                        if (i != j )
                            u1[i] -= poisson[j * n_p + +i] * u[j];
                    }
					/* Now x1_i = (1/a_ii) * (b_i - sum_{i=0,...,n; i \neq j} x_j * a_ij); */
                    u1[i] /=  poisson[i*n_p + i];      
                }

				/* Updating the old u with the new u1 */
                for(i = 0; i < n_p; i++ )
                    u[i] = u1[i]; 

                normb = 0.0;

                for(j = 0; j < n_p; j++ ) {
                    res[j] = b_copy[k*n_p + j];
                    normb += b_copy[k*n_p + j] * b_copy[k*n_p + j];
                }

                normb = sqrt(normb);

                dgemv_("N", &n_p, &n_p, &alpha, poisson, &ldpoisson, u, &one_inc, &beta, res, &one_inc );

                residual = 0.0;

                for(j = 0; j < n_p; j++ )
                    residual += res[j]*res[j];

                residual = sqrt(residual);

                rel_res = residual/normb;
                printf( "[%d-rhs] Norm of || A x - b || / || b  ||  %17.5e\n", k, rel_res);

            }
        }
		
		mem_free(poisson);
		mem_free(ipiv);
		mem_free(J);
		mem_free(I_);
		mem_free(I);
		mem_free(res);
		mem_free(b);
		mem_free(b_copy);
		mem_free(solution);
		mem_free(solution_copy);
	}
	
	return 0;
}

/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, int m, int n, double* a, int lda ) {
	int i, j;
	printf( "\n %s\n", desc );
	for( i = 0; i < m; i++ ) {
		for( j = 0; j < n; j++ ) printf( " %6.2f", a[i+j*lda] );
		printf( "\n" );
	}
}

/* Auxiliary routine: printing a vector of integers */
void print_int_vector( char* desc, int n, int* a ) {
	int j;
	printf( "\n %s\n", desc );
	for( j = 0; j < n; j++ ) printf( " %6i", a[j] );
	printf( "\n" );
}

/* Auxiliary routine: printing a vector of double */
void print_double_vector( char* desc, int n, double* a ) {
	int j;
	printf( "\n %s\n", desc );
	for( j = 0; j < n; j++ ) printf( " %e", a[j] );
	printf( "\n" );
}

void Kronecker_Product(double *C, double *A, int nrows, int ncols, double *B, int mrows, int mcols) {
   int ccols, i, j, k, l;
   int block_increment;
   double *pB;
   double *pC, *p_C;
 
   ccols = ncols * mcols;
   block_increment = mrows * ccols;
   for (i = 0; i < nrows; C += block_increment, i++)
      for (p_C = C, j = 0; j < ncols; p_C += mcols, A++, j++) 
         for (pC = p_C, pB = B, k = 0; k < mrows; pC += ccols, k++)
            for (l = 0; l < mcols; pB++, l++) *(pC+l) = *A * *pB; 

}

double f(double x, double y) {
	return 2*(cos(x + y)-(1 + x)*sin(x + y));
}

double g(double x, double y) {
	return (x+1)*sin(x+y);
}
