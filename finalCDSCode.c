#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h> // Include for timing

// Structure to hold the Compressed Diagonal Storage (CDS) representation
typedef struct
{
    int *LA;       // LA: Array of diagonal offsets. LA[i] is the offset of the i-th diagonal from the main diagonal.
    int **AD;      // AD: 2D array storing the values of the diagonals.  AD[i][j] is the element at row i of the j-th diagonal.
    int num_diags; // Number of diagonals stored in the CDS format.
    int N;         // Dimension of the original matrix (N x N).
} CDSMatrix;

// Function to initialize a CDSMatrix structure
CDSMatrix *createCDSMatrix(int N, int max_diags)
{
    // Allocate memory for the CDSMatrix structure itself
    CDSMatrix *cds = (CDSMatrix *)malloc(sizeof(CDSMatrix));
    if (cds == NULL)
    {
        fprintf(stderr, "Memory allocation failed for CDSMatrix.\n");
        exit(EXIT_FAILURE); // Exit program if allocation fails
    }

    // Allocate memory for the LA (diagonal offsets) array
    cds->LA = (int *)malloc(max_diags * sizeof(int));
    if (cds->LA == NULL)
    {
        fprintf(stderr, "Memory allocation failed for LA array.\n");
        free(cds);          // Free previously allocated memory
        exit(EXIT_FAILURE); // Exit program if allocation fails
    }

    // Allocate memory for the rows of the AD (diagonal values) array
    cds->AD = (int **)malloc(N * sizeof(int *));
    if (cds->AD == NULL)
    {
        fprintf(stderr, "Memory allocation failed for AD rows.\n");
        free(cds->LA);      // Free previously allocated memory
        free(cds);          // Free previously allocated memory
        exit(EXIT_FAILURE); // Exit program if allocation fails
    }

    // Allocate memory for each column of the AD array (for each row)
    for (int i = 0; i < N; i++)
    {
        cds->AD[i] = (int *)malloc(max_diags * sizeof(int));
        if (cds->AD[i] == NULL)
        {
            fprintf(stderr, "Memory allocation failed for AD column.\n");
            // Clean up allocated memory - important for preventing memory leaks
            for (int j = 0; j < i; j++)
            {
                free(cds->AD[j]);
            }
            free(cds->AD);
            free(cds->LA);
            free(cds);
            exit(EXIT_FAILURE); // Exit program if allocation fails
        }
    }

    // Initialize the members of the CDSMatrix structure
    cds->num_diags = 0; // Initially, no diagonals are stored
    cds->N = N;         // Store the matrix dimension
    return cds;         // Return the pointer to the newly created CDSMatrix
}

// Function to deallocate the CDSMatrix structure
void freeCDSMatrix(CDSMatrix *cds)
{
    // Check if the pointer is not NULL before freeing
    if (cds)
    {
        // Free the LA array if it's not NULL
        if (cds->LA)
            free(cds->LA);
        // Free the AD array (rows first, then the array of rows)
        if (cds->AD)
        {
            for (int i = 0; i < cds->N; i++)
            {
                if (cds->AD[i])
                    free(cds->AD[i]); // Free each row
            }
            free(cds->AD); // Free the array of row pointers
        }
        free(cds); // Free the CDSMatrix structure itself
    }
}

// Function to convert a dense matrix to CDS format
CDSMatrix *convertToCDS(int **A, int N, int max_diags)
{
    // Create a new CDSMatrix structure
    CDSMatrix *cds = createCDSMatrix(N, max_diags);

    // Iterate through all possible diagonal offsets
    for (int d = -(N - 1); d <= (N - 1); d++)
    {
        // Check if the diagonal has any non-zero values
        int has_values = 0;
        for (int i = 0; i < N; i++)
        {
            int j = i + d; // Calculate the column index for this diagonal
            if (j >= 0 && j < N && A[i][j] != 0)
            {
                has_values = 1; // Found a non-zero value on this diagonal
                break;          // No need to check further
            }
        }

        // If the diagonal has non-zero values, add it to the CDS format
        if (has_values)
        {
            // Check if the maximum number of diagonals has been exceeded
            if (cds->num_diags >= max_diags)
            {
                fprintf(stderr, "Exceeded maximum number of diagonals allowed.\n");
                freeCDSMatrix(cds); // Free the CDSMatrix structure
                return NULL;        // Return NULL to indicate failure
            }

            // Store the diagonal offset
            cds->LA[cds->num_diags] = d;

            // Store the values of the diagonal in the AD array
            for (int i = 0; i < N; i++)
            {
                int j = i + d; // Calculate the column index for this diagonal
                if (j >= 0 && j < N)
                    cds->AD[i][cds->num_diags] = A[i][j]; // Store the value from the dense matrix
                else
                    cds->AD[i][cds->num_diags] = 0; // If outside the matrix bounds, store 0 (padding)
            }
            cds->num_diags++; // Increment the number of stored diagonals
        }
    }

    return cds; // Return the pointer to the CDSMatrix structure
}

// Function to print a dense matrix (commented out)
// void printMatrix(int **mat, int rows, int cols) {
//     for (int i = 0; i < rows; i++) {
//         for (int j = 0; j < cols; j++) {
//             printf("%3d ", mat[i][j]);
//         }
//         printf("\n");
//     }
// }

// Function to print a 1D array (commented out)
// void printArray(int arr[], int size) {
//     printf("(");
//     for (int i = 0; i < size; i++) {
//         printf("%d", arr[i]);
//         if (i < size - 1) printf(", ");
//     }
//     printf(")\n");
// }

// Function to print a CDS matrix (commented out)
// void printCDSMatrix(CDSMatrix* cds) {
//     if (!cds) return;

//     printf("CDS Matrix AD:\n");
//     printMatrix(cds->AD, cds->N, cds->num_diags);

//     printf("\nDiagonal Offsets (LA): ");
//     printArray(cds->LA, cds->num_diags);
// }

// Function to read matrix from .mtx file (Matrix Market format)
int **readMatrixFromMTX(const char *filename, int *N)
{
    FILE *fp = fopen(filename, "r"); // Open the file in read mode
    if (!fp)
    {
        fprintf(stderr, "Error opening file: %s\n", filename);
        return NULL; // Return NULL if file opening fails
    }

    // Skip header lines (lines starting with '%')
    char line[256]; // Buffer to store each line
    while (fgets(line, sizeof(line), fp))
    {
        if (line[0] != '%')
        {
            break; // Exit loop when a non-comment line is found
        }
    }

    // Read dimensions (rows, cols) and number of non-zero entries (nnz)
    int rows, cols, nnz;
    if (sscanf(line, "%d %d %d", &rows, &cols, &nnz) != 3)
    {
        fprintf(stderr, "Error reading matrix dimensions and nnz.\n");
        fclose(fp);  // Close the file
        return NULL; // Return NULL if reading fails
    }

    // Check if the matrix is square
    if (rows != cols)
    {
        fprintf(stderr, "Matrix must be square.\n");
        fclose(fp);  // Close the file
        return NULL; // Return NULL if matrix is not square
    }

    *N = rows; // Set the matrix dimension N

    // Allocate memory for the matrix
    int **A = (int **)malloc(rows * sizeof(int *));
    if (!A)
    {
        fprintf(stderr, "Memory allocation failed.\n");
        fclose(fp);  // Close the file
        return NULL; // Return NULL if allocation fails
    }
    for (int i = 0; i < rows; i++)
    {
        A[i] = (int *)malloc(cols * sizeof(int));
        if (!A[i])
        {
            fprintf(stderr, "Memory allocation failed.\n");
            for (int j = 0; j < i; j++)
                free(A[j]); // Free previously allocated memory
            free(A);
            fclose(fp);  // Close the file
            return NULL; // Return NULL if allocation fails
        }
        // Initialize the matrix elements to 0
        for (int j = 0; j < cols; j++)
        {
            A[i][j] = 0;
        }
    }

    // Read non-zero entries from the file
    int row, col, value;
    for (int i = 0; i < nnz; i++)
    {
        if (fgets(line, sizeof(line), fp) == NULL)
        {
            fprintf(stderr, "Error reading matrix entry.\n");
            for (int j = 0; j < rows; j++)
                free(A[j]); // Free previously allocated memory
            free(A);
            fclose(fp);  // Close the file
            return NULL; // Return NULL if reading fails
        }
        if (sscanf(line, "%d %d %d", &row, &col, &value) != 3)
        {
            fprintf(stderr, "Error parsing matrix entry.\n");
            for (int j = 0; j < rows; j++)
                free(A[j]); // Free previously allocated memory
            free(A);
            fclose(fp);  // Close the file
            return NULL; // Return NULL if parsing fails
        }
        A[row - 1][col - 1] = value; // MTX format is 1-indexed, so subtract 1
    }

    fclose(fp); // Close the file
    return A;   // Return the pointer to the allocated matrix
}

// Function to perform SpMV (Sparse Matrix-Vector Multiplication) using CDS format
void spmvCDS(CDSMatrix *cds, double *y)
{
    // Check if the CDSMatrix pointer is valid
    if (!cds)
        return;

    int N = cds->N;                 // Matrix dimension
    int num_diags = cds->num_diags; // Number of diagonals

    // Initialize the output vector y to zero
    for (int i = 0; i < N; i++)
    {
        y[i] = 0.0;
    }

    // Iterate over the rows of the matrix
    for (int i = 0; i < N; i++)
    {
        // Iterate over the diagonals
        for (int diag_idx = 0; diag_idx < num_diags; diag_idx++)
        {
            int d = cds->LA[diag_idx]; // Get the diagonal offset
            int j = i + d;             // Calculate the column index for this diagonal

            // Check if the column index is within the matrix bounds
            if (j >= 0 && j < N)
            {
                // Accumulate the result (y[i] += A[i][j] * x[j]) - NOTE: x is no longer used.
                y[i] += (double)cds->AD[i][diag_idx]; // Accumulate values from diagonals of A
            }
        }
    }
}

int main(int argc, char *argv[])
{
    int N;                                        // Matrix dimension
    int **A = NULL;                               // Pointer to the dense matrix
    char *filename = "matrixDataset/bfw398a.mtx"; // Default filename

    // Check for command-line arguments
    if (argc > 1)
    {
        filename = argv[1]; // Use filename provided as a command-line argument
    }

    // Read the matrix from the .mtx file
    A = readMatrixFromMTX(filename, &N);

    // Check if matrix reading was successful
    if (!A)
    {
        return 1; // Exit with an error code
    }

    int MAX_DIAGS = 2 * N - 1; // Maximum possible diagonals in an N x N matrix

    // Convert the dense matrix to CDS format
    CDSMatrix *cds_matrix = convertToCDS(A, N, MAX_DIAGS);

    // (Commented out) Print the original matrix
    // printf("Original Matrix A (from %s):\n", filename);
    // printMatrix(A, N, N);

    // Check if the conversion to CDS format was successful
    if (cds_matrix)
    {
        // (Commented out) Print the CDS matrix
        // printf("\nCompressed Diagonal Storage (CDS) Matrix:\n");
        // printCDSMatrix(cds_matrix);

        // Create an output vector y
        double *y = (double *)malloc(N * sizeof(double));
        if (!y)
        {
            fprintf(stderr, "Memory allocation failed for vector y.\n");
            freeCDSMatrix(cds_matrix); // Free the CDS matrix
            for (int i = 0; i < N; i++)
            {
                free(A[i]); // Free the rows of the original matrix
            }
            free(A);  // Free the original matrix
            return 1; // Exit with an error code
        }

        // Time the SpMV operation
        clock_t start_time = clock();

        // Perform SpMV - No longer passing x
        spmvCDS(cds_matrix, y);

        clock_t end_time = clock();
        double elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

        // (Commented out) Print the result vector y
        // printf("Resulting vector y (A * x): (");
        for (int i = 0; i < N; i++)
        {
            // printf("%.6f", y[i]);
            // if(i < N - 1) printf(" ");
        }
        // printf(")\n");

        printf("Time taken for SpMV: %f seconds\n", elapsed_time);

        free(y);                   // Free the output vector
        freeCDSMatrix(cds_matrix); // Deallocate the memory associated with CDSMatrix
    }

    // Deallocate the dynamically allocated matrix A
    for (int i = 0; i < N; i++)
    {
        free(A[i]); // Free each row of the matrix
    }
    free(A); // Free the matrix itself

    return 0; // Exit with success code
}