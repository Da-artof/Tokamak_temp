#define  _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <mkl_lapacke.h>

long read_input(char *fname, long *Nt, long *Nz, double *tf, long *Im) {

        FILE* fptr=fopen(fname,"r");

        if (fptr==NULL) {
		printf("NULL");
		return 1;
	}
	
        if (4!=fscanf(fptr,"%ld %ld %lf %ld", Nt, Nz, tf, Im)) {
                return 1;
        }

        fclose(fptr);
        return 0;
}

long read_coefficients(char *fname, double Q11[], double Q22[], double Q12[], double S[], double R[]) {

	
	FILE * fp;
    	char * line = NULL;
	char * token = NULL;
	int i = 0, j = 0;
	size_t len = 0;
    	ssize_t read;

    	fp = fopen(fname, "r");
   	
	if (fp == NULL)
		return 1;


  	while ((read = getline(&line, &len, fp)) != -1) {
		token = strtok(line, " ");
   		while( token != NULL ) {

             		if (j==0) {
				Q11[i] = atof(token);
			}

             		if (j==1) {
				Q22[i+1] = atof(token);
			}

             		if (j==2) {
				Q12[i] = atof(token);
			}

             		if (j==3) {
				S[i] = atof(token);
			}

             		if (j==4) {
				R[i] = atof(token);
			}

                   	token = strtok(NULL, " ");
			j += 1;
              	}
		j = 0;
		i += 1;
    	}
    	fclose(fp);
    	if (line) {
        	free(line);
	}
	return 0;
}

long fold(long i, long N) {
	if ( (double) i < (double) N / 2) {
		return 2 * i;
	}
	else {
		return 2 * (N - i) - 1;
	}
}

long indx( long j, long p, long P) {
  	return j*P + p;
}

void write_output(double *Ti, double dth, double dz, long Nt, long Nz, double tf) {
	FILE *fptr;
	fptr = fopen("output.txt", "w");
	if(fptr == NULL) {
      		printf("File write error");   
      		exit(1);             
   	}
	
	fprintf(fptr, "0.0 0.0 \n");
	
	long i, j;
	double theta = 0, zeta = 0;

	for (i = 0; i < Nt; i++) {
		for (j=0; j < Nz; j++) {
			fprintf(fptr, "%lf %lf %lf %lf\n", tf, theta, zeta, Ti[indx(i, j, Nz)]);
			zeta += dz;
		}
		zeta = 0;
		theta += dth;
	}

	fclose(fptr);
}

struct band_mat{
	
 	long ncol;       
  	long nbrows;    
  	long nbands_up;
  	long nbands_low; 
  	double *array;  
  	long nbrows_inv;
  	double *array_inv;
  	int *ipiv;       

};

typedef struct band_mat band_mat;

int init_band_mat(band_mat *bmat, long nbands_lower, long nbands_upper, long n_columns) {
 	
	bmat->nbrows = nbands_lower + nbands_upper + 1;
 	bmat->ncol   = n_columns;
 	bmat->nbands_up = nbands_upper;
 	bmat->nbands_low= nbands_lower;
 	bmat->array      = (double *) malloc(sizeof(double)*bmat->nbrows*bmat->ncol);
 	bmat->nbrows_inv = bmat->nbands_up*2 + bmat->nbands_low + 1;
 	bmat->array_inv  = (double *) malloc(sizeof(double)*(bmat->nbrows+bmat->nbands_low)*bmat->ncol);
 	bmat->ipiv       = (int *) malloc(sizeof(int)*bmat->ncol);
 	
	if (bmat->array==NULL||bmat->array_inv==NULL) {
 		return 0;
	}  
  	
	long i;
  	
	for (i=0;i<bmat->nbrows*bmat->ncol;i++) {
    		bmat->array[i] = 0.0;
  	}
  	
	return 1;

}; 

double *getp(band_mat *bmat, long row, long column) {
  	int bandno = bmat->nbands_up + row - column;
  	if(row<0 || column<0 || row>=bmat->ncol || column>=bmat->ncol ) {
    		printf("Indexes out of bounds in getp: %ld %ld %ld \n",row,column,bmat->ncol);
    		exit(1);
  	} 	
	return &bmat->array[bmat->nbrows*column + bandno];
}

double getv(band_mat *bmat, long row, long column) {
  	return *getp(bmat,row,column);
}

double setv(band_mat *bmat, long row, long column, double val) {
  	*getp(bmat,row,column) = *getp(bmat,row,column) + val;
  	return val;
}

int solve_Ax_eq_b(band_mat *bmat, double *x, double *b) {
  	int i,bandno;
  	for(i=0;i<bmat->ncol;i++) { 
    		for (bandno=0;bandno<bmat->nbrows;bandno++) {
      			bmat->array_inv[bmat->nbrows_inv*i+(bandno+bmat->nbands_low)] = bmat->array[bmat->nbrows*i+bandno];
  		}
    		x[i] = b[i];
  	}

  	long nrhs = 1;
  	long ldab = bmat->nbands_low*2 + bmat->nbands_up + 1;
  	int info = LAPACKE_dgbsv( LAPACK_COL_MAJOR, bmat->ncol, bmat->nbands_low, bmat->nbands_up, nrhs, bmat->array_inv, ldab, bmat->ipiv, x, bmat->ncol);
  	return info;
}

int calculate_y_eq_Ax_plus_y(double **A, double *x, double *y, long ncol) {
	long i, j;
	for (i=0; i < ncol; i++) {
        	for (j=0; j < ncol; j++) {
            		y[i] += A[i][j] * x[j];
		}
    	}

	return 1;
}

void finalise_band_mat(band_mat *bmat) {
  	free(bmat->array);
  	free(bmat->array_inv);
  	free(bmat->ipiv);
}


long boundary_indx(long j, long k, long Nt, long Nz) {
	long j_new = j, k_new = k;
	if (j < 0) {j_new += Nt;}
	if (j >= Nt) {j_new -= Nt;}
	if (k < 0) {k_new += Nz;}
	if (k >= Nz) {k_new -= Nz;}
	return indx(j_new, k_new, Nz);
}


int main() {

	long Nt, Nz, Im; 
	double tf;

	if (read_input("input.txt", &Nt, &Nz, &tf, &Im)) {
		printf("File read error\n");
	}

	long L = Nt * Nz;

	double *Q11, *Q22, *Q12, *S, *R;
	Q11 = malloc(sizeof(double)*L);
	Q22 = malloc(sizeof(double)*L);
	Q12 = malloc(sizeof(double)*L);
	S = malloc(sizeof(double)*L);
	R = malloc(sizeof(double)*L);

	if (Q11 == NULL || Q22 == NULL || Q12 == NULL || S == NULL || R == NULL) {
		printf("Coefficent memory allocation error");
		return 1;
	}
	
	read_coefficients("coefficients.txt", Q11, Q22, Q12, S, R);

	double pi = 3.141592654;
	double dt = tf / (double)Im, dth = 2*pi/(double)Nt, dz = 2*pi/(double)Nz;

	/* Storage for Temperature at time i and RHS of equation*/
	
	double *Ti, *RHS;
		
	Ti = malloc(sizeof(double)*L);
	RHS = malloc(sizeof(double)*L);

	
	if (Ti == NULL || RHS == NULL) {
		printf("Temperature memory allocation error");
		return 1;
	}

	/*Initialise T1 = 0*/
	
	long i = 0;
	for (; i < L; i++) {
		Ti[i] = 0;
	}

	/*Storage for matrices*/
	
	band_mat A1;
	double **A2 = (double **)malloc(L * sizeof(double *));
   	for (i = 0; i < L; i++) {
        	A2[i] = (double *)malloc(L * sizeof(double));
    	}

	if (A2 == NULL) {
		printf("Matrix memory allocation error");
		return 1;
	}
	
	long nbands_low;
	
	if (Nt>Nz) {nbands_low = Nt-1;} else {nbands_low = Nz-1;};
	
	long nbands_up = nbands_low;
	
	init_band_mat(&A1, nbands_up, nbands_low, L);
	
	/*Set matrix values*/
	
	long j, k, unknown_indx, bound_indx;
	
	for (j = 0; j < Nt; j++) {	
		for (k = 0; k < Nz; k++) {
                        unknown_indx = indx(j,k,Nz);
			
			bound_indx = boundary_indx(j-1, k-1, Nt, Nz);
			setv(&A1, unknown_indx, bound_indx, -(Q12[boundary_indx(j-1, k, Nt, Nz)] + Q12[boundary_indx(j, k-1, Nt, Nz)])/(8*dth*dz));
                        A2[unknown_indx][bound_indx] = (Q12[boundary_indx(j-1, k, Nt, Nz)] + Q12[boundary_indx(j, k-1, Nt, Nz)])/(8*dth*dz);
       	                
		        bound_indx = boundary_indx(j-1, k+1, Nt, Nz);
                        setv(&A1, unknown_indx, bound_indx, -(Q12[boundary_indx(j-1, k, Nt, Nz)] + Q12[boundary_indx(j, k+1, Nt, Nz)])/(8*dth*dz));
                        A2[unknown_indx][bound_indx] =(Q12[boundary_indx(j-1, k, Nt, Nz)] + Q12[boundary_indx(j, k+1, Nt, Nz)])/(8*dth*dz);
                                
			bound_indx = boundary_indx(j+1, k-1, Nt, Nz);
                        setv(&A1, unknown_indx, bound_indx, -(Q12[boundary_indx(j+1, k, Nt, Nz)] + Q12[boundary_indx(j, k-1, Nt, Nz)])/(8*dth*dz));
                        A2[unknown_indx][bound_indx] = (Q12[boundary_indx(j+1, k, Nt, Nz)] + Q12[boundary_indx(j, k-1, Nt, Nz)])/(8*dth*dz);
                                
			bound_indx = boundary_indx(j+1, k+1, Nt, Nz);
                        setv(&A1, unknown_indx, bound_indx, -(Q12[boundary_indx(j+1, k, Nt, Nz)] + Q12[boundary_indx(j, k+1, Nt, Nz)])/(8*dth*dz));
                        A2[unknown_indx][bound_indx] = (Q12[boundary_indx(j+1, k, Nt, Nz)] + Q12[boundary_indx(j, k+1, Nt, Nz)])/(8*dth*dz);
                                
			bound_indx = boundary_indx(j-2, k, Nt, Nz);
                        setv(&A1, unknown_indx, bound_indx, -Q11[boundary_indx(j-1, k, Nt, Nz)]/(8*dth*dth));
                        A2[unknown_indx][bound_indx] = Q11[boundary_indx(j-1, k, Nt, Nz)]/(8*dth*dth);
                                
			bound_indx = boundary_indx(j+2, k, Nt, Nz);
                        setv(&A1, unknown_indx, bound_indx, -Q11[boundary_indx(j+1, k, Nt, Nz)]/(8*dth*dth) );
                        A2[unknown_indx][bound_indx] = Q11[boundary_indx(j+1, k, Nt, Nz)]/(8*dth*dth);
                                
			bound_indx = boundary_indx(j, k-2, Nt, Nz);
                        setv(&A1, unknown_indx, bound_indx, -Q22[boundary_indx(j, k-1, Nt, Nz)]/(8*dz*dz));
                        A2[unknown_indx][bound_indx] = Q22[boundary_indx(j, k-1, Nt, Nz)]/(8*dz*dz);
                                
			bound_indx = boundary_indx(j, k+2, Nt, Nz);
                        setv(&A1, unknown_indx, bound_indx, -Q22[boundary_indx(j, k+1, Nt, Nz)]/(8*dz*dz));
                        A2[unknown_indx][bound_indx] = Q22[boundary_indx(j, k+1, Nt, Nz)]/(8*dz*dz);
               	        
			setv(&A1, unknown_indx, unknown_indx,-((Q11[boundary_indx(j+1, k, Nt, Nz)]-Q11[boundary_indx(j-1, k, Nt, Nz)])/(8*dth*dth) + (Q22[boundary_indx(j, k+1, Nt, Nz)]-Q22[boundary_indx(j, k-1, Nt, Nz)])/(8*dz*dz) - R[unknown_indx]/2 + 1/dt));
               	        A2[unknown_indx][unknown_indx] = (Q11[boundary_indx(j+1, k, Nt, Nz)]-Q11[boundary_indx(j-1, k, Nt, Nz)])/(8*dth*dth) + (Q22[boundary_indx(j, k+1, Nt, Nz)]-Q22[boundary_indx(j, k-1, Nt, Nz)])/(8*dz*dz) - R[unknown_indx]/2 + 1/dt;
                }               
        }

	/*Solve A1 * Ti1 = A2 * Ti + S for i = 1, ... , Im */               
	
	for (i = 0; i < Im - 1; i++) {
		for (j = 0; j < L; j++) {
			RHS[j] = S[j];
		}
		calculate_y_eq_Ax_plus_y(A2, Ti, RHS, L);
		solve_Ax_eq_b(&A1, Ti, RHS);
	}

	write_output(Ti, dth, dz, Nt, Nz, tf);
	
	free(Q11);
	free(Q22);
	free(Q12);
	free(S);
	free(R);
	free(Ti);
	free(RHS);
	free(A2);
	finalise_band_mat(&A1);
}

