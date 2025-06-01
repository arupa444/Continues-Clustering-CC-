/*
*   Matrix Market I/O example program
*
*   Read a real (non-complex) sparse matrix from a Matrix Market (v. 2.0) file.
*   and copies it to stdout.  This porgram does nothing useful, but
*   illustrates common usage of the Matrix Matrix I/O routines.
*   (See http://math.nist.gov/MatrixMarket for details.)
*
*   Usage:  a.out [filename] > output
*
*
*   NOTES:
*
*   1) Matrix Market files are always 1-based, i.e. the index of the first
*      element of a matrix is (1,1), not (0,0) as in C.  ADJUST THESE
*      OFFSETS ACCORDINGLY offsets accordingly when reading and writing
*      to files.
*
*   2) ANSI C requires one to use the "l" format modifier when reading
*      double precision floating point numbers in scanf() and
*      its variants.  For example, use "%lf", "%lg", or "%le"
*      when reading doubles, otherwise errors will occur.
*/

#include <stdio.h>
#include <stdlib.h>
#include<time.h>
#include<math.h>
#include<omp.h>
#include "mmio.h"

/* coo is a structure for COO format
*  structure members:
*  nnz: integer value that stores the number of non-zero elements in the sparse matrix
*  rows: integer array to store the row index of the non-zero entries
*  cols: integer array to store the column index of the non-zero entries
*  values: double array to store the non-zero values
*
*  row index of rows[i] and column index of cols[i] contains value values[i]
*
*/

struct coo {
    int nrows;
    int ncols;
    int nnz;
    int* rows;
    int* cols;
    double* values;
};

/* dia is structure for DIA or Diagonal format
*  structure members:
*  nrows: integer value that stores the number of rows of the sparse matrix
*  ncols: integer value that stores the number of columns of the sparse matrix
*  ndiags: integer value that stores the number of diagonals that contain
*  at-least one non-zero element.
*  offset: integer array to store the offset of each diagonal from the main diagonal.
*  By convention, the main diagonal corresponds to offset 0, while i>0 represents ith
*  super-diagonal and i<0 the ith sub-diagonal.
*  data: double 1D array to store the non-zero values. Basically data stores the
*  diagonals in a 2D array in a column-major order. We linearize it so that elements
*  within each diagonal are placed adjacently.
*/

struct dia {
    int nrows;
    int ncols;
    int nnz;
    int ndiags;
    int* offset;
    double* data;
    int data_len;

};

struct ccs {
    int nrows;
    int ncols;
    int nnz;
    int* cluster_rows;
    int* cluster_cols;
    int* freq;
    double* values;
};

struct coo* allocate_mem_coo(struct coo* my_sparse_matrix, int nrows, int ncols, int nnz) {
    my_sparse_matrix = (struct coo*)malloc(sizeof(struct coo));
    my_sparse_matrix->nrows = nrows;
    my_sparse_matrix->ncols = ncols;
    my_sparse_matrix->nnz = nnz;
    my_sparse_matrix->rows = (int*)malloc(sizeof(int) * nnz);
    my_sparse_matrix->cols = (int*)malloc(sizeof(int) * nnz);
    my_sparse_matrix->values = (double*)malloc(sizeof(double) * nnz);
    return my_sparse_matrix;
}

int compare(const void* elem1, const void* elem2)
{
    int f = *((int*)elem1);
    int s = *((int*)elem2);
    if (f > s) return  1;
    if (f < s) return -1;
    return 0;
}

struct ccs* coo_to_ccs(struct coo* coo_mat) {
    struct ccs* ccs_mat = (struct ccs*)malloc(sizeof(struct ccs));
    ccs_mat->nrows = coo_mat->nrows;
    ccs_mat->ncols = coo_mat->ncols;
    ccs_mat->nnz = coo_mat->nnz;

    /* Creating an array "diff" which stores the 
    *  value of [col_index - row_index] of every
    *  non-zero element present in the matrix. Same
    *  [diff] value will be in a single diaginal.
    */

    int* diff = (int*) malloc(sizeof(int) * coo_mat->nnz);

    for (int i = 0; i < coo_mat->nnz; i++) {        
        diff[i] = coo_mat->cols[i] - coo_mat->rows[i];                
    }

    /*
    *  Sorting the [diff] array so that the array is
    *  sorted according to the [diff] value. This will
    *  sort and cluster all the diagonal indices.
    */

    qsort(diff, coo_mat->nnz, sizeof(int), compare);

    /* Calculating the size of the cluster_rows and cluster_cols
    *  arrays i.e. calculating how many clusters are there.
    */

    for (int i = 0; i < coo_mat->nnz - 1; i++) {

    }

    return ccs_mat;
}

struct dia* coo_to_dia(struct coo* coo_mat) {
    struct dia* dia_mat = (struct dia*)malloc(sizeof(struct dia));
    dia_mat->nrows = coo_mat->nrows;
    dia_mat->ncols = coo_mat->ncols;
    dia_mat->nnz = coo_mat->nnz;
    int* duplicate_offsets = (int*)malloc(sizeof(int) * dia_mat->nnz);
    int ndiags = 1;

    for (int i = 0; i < coo_mat->nnz; i++) {
        duplicate_offsets[i] = coo_mat->cols[i] - coo_mat->rows[i];
    }

    qsort(duplicate_offsets, coo_mat->nnz, sizeof(int), compare);

    for (int i = 1; i < coo_mat->nnz; i++) {
        if (duplicate_offsets[i - 1] != duplicate_offsets[i]) {
            ndiags++;
        }
    }

    dia_mat->ndiags = ndiags;

    dia_mat->offset = (int*)malloc(sizeof(int) * ndiags);
    ndiags = 1;
    dia_mat->offset[0] = duplicate_offsets[0];
    for (int i = 1; i < coo_mat->nnz; i++) {
        if (duplicate_offsets[i - 1] != duplicate_offsets[i]) {
            dia_mat->offset[ndiags] = duplicate_offsets[i];
            ndiags++;
        }
    }

    free(duplicate_offsets);

    double** temp_data = (double*)calloc(max(dia_mat->nrows, dia_mat->ncols), sizeof(double*));
    for (int i = 0; i < max(dia_mat->nrows, dia_mat->ncols); i++) {
        temp_data[i] = (double*)calloc(ndiags, sizeof(double));
    }


    for (int i = 0; i < coo_mat->nnz; i++) {
        int row = coo_mat->rows[i];
        int col = coo_mat->cols[i];
        int val = coo_mat->values[i];
        int diff = col - row;
        int offset_ind = lin_search(dia_mat->offset, ndiags, diff);

        temp_data[row][offset_ind] = val;

    }

    /* counting the total number of entries in linear data array */

    /*int data_len = 0;
    for (int i = 0; i < ndiags; i++) {
        if (dia_mat->offset[i] < 0)
            data_len += max(dia_mat->nrows, dia_mat->ncols) + dia_mat->offset[i];
        else if (dia_mat->offset[i] > 0)
            data_len += max(dia_mat->nrows, dia_mat->ncols) - dia_mat->offset[i];
        else
            data_len += max(dia_mat->nrows, dia_mat->ncols);
    }
    dia_mat->data_len = data_len;*/
    dia_mat->data = (double*)malloc(sizeof(double) * ndiags * dia_mat->nrows);

    int k = 0;
    for (int i = 0; i < ndiags; ++i) {
        int start = 0, end = max(dia_mat->nrows, dia_mat->ncols);
        /*if (dia_mat->offset[i] < 0)
            start -= dia_mat->offset[i];
        else
            end -= dia_mat->offset[i];*/
        for (int j = 0; j < dia_mat->nrows; ++j) {
            dia_mat->data[k++] = temp_data[j][i];
        }
    }

    free(temp_data);


    return dia_mat;
}

int get_size_coo(int nnz) {
    return ((sizeof(int) * (2 * nnz + 1)) + (nnz * (sizeof(double))));
}

/* create a vector of length "size" with random double values
*  from "min" to "max". Returns a pointer to the vector.
*/

double* get_random_vector(int size, int min, int max) {
    double* vec = (double*)malloc(sizeof(double) * size);
    for (int i = 0; i < size; i++) {
        int rd_num = rand() % (max - min + 1) + min;
        vec[i] = rd_num;
    }

    return vec;
}

/* Does SpMV a COO matrix with a dense vector
*  and store the result in another dense vector
*/

void SpMV_COO(struct coo* A, double* x, double* y) {
    for (int i = 0; i < A->nnz; ++i) {
        int column = A->cols[i];
        int row = A->rows[i];
        y[row] += A->values[i] * x[column];
    }
}

/* Does SpMV a DIA matrix with a dense vector
*  and store the result in another dense vector
*/

void SpMV_DIA(struct dia* A, double* x, double* y) {
    int i = 0;
    int stride = fmax(A->nrows, A->ncols);
    int n = 0;
#pragma omp parallel for schedule(static, A->ndiags)
    for (i = 0; i < A->ndiags; i++) {
        int k = A->offset[i];
        int istart = fmax(0, -k);
        int jstart = max(0, k);

        int N = fmin(A->nrows - istart, A->ncols - jstart);

        for (n = 0; n < N; n++) {
            y[istart + n] += A->data[istart + i * stride + n] * x[jstart + n];
        }
    }


}

void SpMV_openmp(struct coo* A, double* x, double* y) {
    int i = 0;
#pragma omp parallel for
    for (i = 0; i < A->nnz; i++)
        y[A->rows[i]] += A->values[i] * x[A->cols[i]];

}

/* Helper function to get a minimum of two numbers */

int minimum(int a, int b) {
    return a < b ? a : b;
}

/* Utility finction to search an value in an array offset of length size */

int lin_search(int* offset, int size, int value) {
    int i = 0;
    for (i = 0; i < size; i++) {
        if (offset[i] == value)
            break;
    }

    return i;
}

int main(int argc, char* argv[])
{
    int ret_code;
    MM_typecode matcode;
    FILE* f;
    int M, N, nz;
    double* val;

    if (argc < 2)
    {

        fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
        exit(1);
    }
    else
    {
        if ((f = fopen(argv[1], "r")) == NULL)
            exit(1);
    }

    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }


    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */

    if (mm_is_complex(matcode) && mm_is_matrix(matcode) &&
        mm_is_sparse(matcode))
    {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }

    /* find out size of sparse matrix .... */

    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) != 0)
        exit(1);

    /* create a sparse matrix structure in coo format and reserve memory */

    struct coo* my_sparse_matrix = NULL;
    my_sparse_matrix = allocate_mem_coo(my_sparse_matrix, M, N, nz);
    //printf("The COO format takes %lu bytes space in memory\n", get_size_coo(nz));

    for (int i = 0; i < nz; i++)
    {
        fscanf(f, "%d %d %lf\n", &my_sparse_matrix->rows[i], &my_sparse_matrix->cols[i], &my_sparse_matrix->values[i]);
        my_sparse_matrix->rows[i]--;  /* adjust from 1-based to 0-based */
        my_sparse_matrix->cols[i]--;
    }

    if (f != stdin) fclose(f);

    /* Create a random vector of length nz */

    //double* x = get_random_vector(N, 5, 10);

    /* Create a random vector for storing the result of SpMV */

    //double* y = (double*)calloc(M, sizeof(double));

    /* SpMV using COO storage format in sequential mode */

    // clock_t start = clock();
    SpMV_COO(my_sparse_matrix, x, y);
    clock_t end = clock();
    double time_required = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Time for SpMV for matrix having %d non-zeros is %f\n", my_sparse_matrix->nnz, time_required);
    printf("GFlops for SpMV for matrix having %d non-zeros is %f\n", my_sparse_matrix->nnz, (2 * my_sparse_matrix->nnz) / (time_required * pow(10.0, 9.0)));*/

    /* SpMV using COO storage format in parallel mode */

    clock_t start = clock();
    SpMV_openmp(my_sparse_matrix, x, y);
    clock_t end = clock();
    double time_required = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Time for SpMV for matrix having %d non-zeros is %f\n", my_sparse_matrix->nnz, time_required);
    printf("GFlops for SpMV for matrix having %d non-zeros is %f\n", my_sparse_matrix->nnz, (2 * my_sparse_matrix->nnz) / (time_required * pow(10.0, 9.0)));*/

    //struct dia* dia_mat = coo_to_dia(my_sparse_matrix);
    struct ccs* ccs_mat = coo_to_ccs(my_sparse_matrix);

    free(my_sparse_matrix->rows);
    free(my_sparse_matrix->cols);
    free(my_sparse_matrix->values);
    free(my_sparse_matrix);

    /* printing the information of a dia matrix */

    printf("Number of rows: %d\n", dia_mat->nrows);
    printf("Number of columns: %d\n", dia_mat->ncols);
    printf("Number of non-zero elements: %d\n", dia_mat->nnz);
    printf("Number of diagonals: %d\n", dia_mat->ndiags);*/

    clock_t start = clock();
    SpMV_DIA(dia_mat, x, y);
    clock_t end = clock();
    double time_required = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Time for SpMV for matrix having %d non-zeros is %f\n", dia_mat->nnz, time_required);
    printf("GFlops for SpMV for matrix having %d non-zeros is %f\n", dia_mat->nnz, (2 * dia_mat->nnz) / (time_required * pow(10.0, 9.0)));*/

    free(x);
    free(y);
    free(ccs_mat->offset);
    free(dia_mat->data);
    free(ccs_mat);


    return 0;
}