#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
typedef struct matcsr{
	double *val;
	int	*colindx, *rowptr;
	int	valsize, colsize, rowptrsize;	
}*MAT_CSR;

typedef struct matdense{
	double *val;
	int	row, col;
}*MAT_DENSE;

typedef struct vec{
	double *val;
	int	size;
}*VEC;

int Max(int array[],int size){
	int max;
	int i;
	max =array[0];

	for(i=1;i<size;i++){
		if(array[i] > max){
			max = array[i];
		}	
	}

	return max;
}

void VecView(VEC vec){
	int i;
	for(i = 0; i < vec->size;i++){
		printf("%f\n", vec->val[i]);
	}
}

void MatView(MAT_DENSE matdense){
	int i,j;

        for(i = 0; i< matdense->row; i++){
		printf("row %d: ", i);
                for(j = 0; j < matdense->col;j++){
            		printf("%f   ",  matdense->val[j+matdense->col*i]);
                }
		printf("\n");
        }
}

MAT_DENSE Csr2Dense(MAT_CSR mat){
	int i,j;
	int start, end;
	MAT_DENSE matdense;

	matdense = (MAT_DENSE) malloc(sizeof(struct matdense));
	matdense->row = mat->rowptrsize-1;
	matdense->col = 1 + Max(mat->colindx, mat->colsize);
	matdense->val = (double*)malloc(matdense->row*matdense->col*sizeof(double));
	
	for(i = 0; i< matdense->row; i++){
		for(j = 0; j < matdense->col;j++){
			matdense->val[j+matdense->col*i] = 0.0;		
		}
	}
        for(i = 0; i< matdense->row; i++){
                start = mat->rowptr[i];;
                end = mat->rowptr[i+1];
                for(j = start; j < end; j++){
                       matdense->val[i*matdense->col+mat->colindx[j]] = mat->val[j];
                }
        }
	return matdense;
}

VEC omp_spmv_csr(MAT_CSR mat, VEC vec){

	int i,j;
	int rowsize;
	int start, end;
	VEC sol;
	int tmp;

	sol = (VEC)malloc(sizeof(struct vec));
	rowsize = mat->rowptrsize - 1;
	sol->size = rowsize;
	sol->val = (double*)malloc(sol->size*sizeof(double));
	#pragma omp parallel for private(i,j,tmp)
	for(i = 0; i< rowsize; i++){
		sol->val[i] = 0.0; 
		for(j = mat->rowptr[i]; j < mat->rowptr[i+1]; j++){
			tmp  = mat->val[j]*vec->val[mat->colindx[j]];
			#pragma omp atomic
			sol->val[i]+=tmp;
		}	
	}

	return sol;
}

int main (int argc, char *argv[]) 
{
int tid, nthreads;

int a[9] = {1,7,2,8,5,3,9,6,4};
int b[9] = {0,1,1,2,0,2,3,1,3};
int c[5] = {0,2,4,7,9};
int d[4] = {1,1,1,1};

int i;

MAT_CSR mat;
VEC	vec;
VEC	sol;

mat = (MAT_CSR) malloc(sizeof(struct matcsr));
vec = (VEC)malloc(sizeof(struct vec));

mat->valsize = 9;
mat->colsize = 9;
mat->rowptrsize = 5;
vec->size = 4;

mat->val = (double*)malloc(mat->valsize*sizeof(double));
mat->colindx = (int*)malloc(mat->colsize*sizeof(int));
mat->rowptr = (int*)malloc(mat->rowptrsize*sizeof(int));
vec->val = (double*)malloc(vec->size*sizeof(double));

for(i = 0;i < mat->valsize;i++){
	mat->val[i] = a[i];
}
for(i = 0;i < mat->colsize;i++){
        mat->colindx[i] = b[i];
}
for(i = 0;i < mat->rowptrsize;i++){
        mat->rowptr[i] = c[i];
}
for(i = 0;i < vec->size;i++){
        vec->val[i] = d[i];
}

sol = omp_spmv_csr(mat, vec);
VecView(sol);

/*
MAT_DENSE mt;
mt = Csr2Dense(mat);
MatView(mt);
int j;
j = mt->col;
printf("j = %d\n",j);
*/

/* 
       for(i = 0; i< mt->row; i++){
                printf("row %d: ", i);
                for(j = 0; j < mt->col;j++){
                        printf("%f   ",  mt->val[i+mt->col*j]);
                }
                printf("\n");
        }
*/
/*
double *val;
int *colindx, *rowptr;

double *vec;
int i;

val = malloc(9*sizeof(double));
colindx = malloc(9*sizeof(int));
rowptr = malloc(5*sizeof(int));
vec = malloc(4*sizeof(double));


val[9] ={1,7,2,8,5,3,9,6,4};
colindx[9] = {0,1,1,2,0,2,3,1,3};
rowptr[5] = {0,2,4,7,9};
vec[4] = {1,1,1,1};
*/

/* Fork a team of threads giving them their own copies of variables */
#pragma omp parallel private(nthreads, tid)
  {

  /* Obtain thread number */
  tid = omp_get_thread_num();
  printf("Hello World from thread = %d\n", tid);

  /* Only master thread does this */
  if (tid == 0) 
    {
    nthreads = omp_get_num_threads();
    printf("Number of threads = %d\n", nthreads);
    }

  }  /* All threads join master thread and disband */

}

