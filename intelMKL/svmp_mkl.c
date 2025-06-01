#include <stdio.h>
#include <stdlib.h>
#include <mkl.h>

void sparseMatrixVectorMultiplication(int rows, int cols, int* csrRowIndex, int* colIndex, double* values, double* vector, double* result) {
    sparse_matrix_t sparseMatrix;
    struct matrix_descr descr;
    descr.type = SPARSE_MATRIX_TYPE_GENERAL;

    // Create the sparse matrix handle
    mkl_sparse_d_create_csr(&sparseMatrix, SPARSE_INDEX_BASE_ZERO,
                            rows, cols,
                            csrRowIndex, csrRowIndex + 1, colIndex, values);

    // Perform matrix-vector multiplication
    mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE,
                    1.0,                 // Alpha
                    sparseMatrix,        // Sparse matrix handle
                    descr,               // Matrix descriptor
                    vector,              // Input vector
                    0.0,                 // Beta
                    result);             // Output vector

    // Print the result
    printf("Result of Sparse Matrix-Vector Multiplication:\n");
    for (int i = 0; i < rows; i++) {
        printf("%f\n", result[i]);
    }

    // Destroy the sparse matrix handle
    mkl_sparse_destroy(sparseMatrix);
}

int main() {
    FILE* file = fopen("try.mtx", "r");
    if (file == NULL) {
        perror("Error opening file");
        return EXIT_FAILURE;
    }

    char line[256];
    // Skip comments and headers in Matrix Market file
    while (fgets(line, sizeof(line), file)) {
        if (line[0] != '%' && line[0] != '\n') break;
    }

    int rows, cols, nnz;
    sscanf(line, "%d %d %d", &rows, &cols, &nnz);

    // Allocate memory for COO format
    int* rowIndex = (int*)malloc(nnz * sizeof(int));
    int* colIndex = (int*)malloc(nnz * sizeof(int));
    double* values = (double*)malloc(nnz * sizeof(double));

    if (!rowIndex || !colIndex || !values) {
        fprintf(stderr, "Memory allocation failed\n");
        return EXIT_FAILURE;
    }

    // Read COO data from file
    int row, col, i = 0;
    double value;
    while (fscanf(file, "%d %d %lf", &row, &col, &value) == 3) {
        if (row >= 1 && row <= rows && col >= 1 && col <= cols) {
            rowIndex[i] = row - 1;  // Convert 1-based to 0-based index
            colIndex[i] = col - 1; // Convert 1-based to 0-based index
            values[i] = value;
            i++;
        } else {
            fprintf(stderr, "Warning: Index out of bounds (%d, %d)\n", row, col);
        }
    }
    fclose(file);

    if (i != nnz) {
        fprintf(stderr, "Error: Number of non-zero elements read (%d) does not match expected (%d)\n", i, nnz);
        return EXIT_FAILURE;
    }

    // Convert COO format to CSR format
    int* csrRowIndex = (int*)malloc((rows + 1) * sizeof(int));
    if (!csrRowIndex) {
        fprintf(stderr, "Memory allocation failed\n");
        return EXIT_FAILURE;
    }

    memset(csrRowIndex, 0, (rows + 1) * sizeof(int));
    for (int j = 0; j < nnz; j++) {
        csrRowIndex[rowIndex[j] + 1]++;
    }
    for (int j = 1; j <= rows; j++) {
        csrRowIndex[j] += csrRowIndex[j - 1];
    }

    // Allocate and initialize vector and result arrays
    double* vector = (double*)malloc(cols * sizeof(double));
    double* result = (double*)malloc(rows * sizeof(double));
    if (!vector || !result) {
        fprintf(stderr, "Memory allocation failed\n");
        return EXIT_FAILURE;
    }

    for (int j = 0; j < cols; j++) vector[j] = 1.0;  // Initialize vector with 1.0
    for (int j = 0; j < rows; j++) result[j] = 0.0;  // Initialize result with 0.0

    // Perform the sparse matrix-vector multiplication
    sparseMatrixVectorMultiplication(rows, cols, csrRowIndex, colIndex, values, vector, result);

    // Free dynamically allocated memory
    free(rowIndex);
    free(colIndex);
    free(values);
    free(csrRowIndex);
    free(vector);
    free(result);

    return 0;
}
