#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

int offsetFunct(int rows, int cols, double** matrix, int** offsets, int x)
{
    // Allocate memory for offsets
    *offsets = (int *)malloc((rows + cols - 1) * sizeof(int));
    if (*offsets == NULL)
    {
        printf("Memory allocation failed!\n");
        return -1; // Return an error code if memory allocation fails
    }

    int off_count = 0; // Counter for non-zero diagonals

    // Initialize offsets to -1 (indicating no valid offset found yet)
    for (int i = 0; i < rows + cols - 1; i++)
    {
        (*offsets)[i] = x;
    }

    // Fill offsets based on non-zero entries in the matrix
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            if (matrix[i][j] != 0)
            {
                // Calculate the diagonal index
                int diag_index = j - i + (rows - 1);
                (*offsets)[diag_index] = j - i; // Store the offset
                off_count++;
            }
        }
    }
    return off_count; // Return the count of non-zero offsets
}

int compressedDIA(double** matrix, int rows, int cols, int x)
{

    clock_t start, end;
    double cpu_time_used;
    long int space=0;
    double flops =0;
    double* vector = (double*)malloc((cols) * sizeof(double));

    for (int i = 0; i < cols; i++) {
        vector[i] = 1.0; 
    }

    double* result = (double*)malloc((rows) * sizeof(double));

    for (int i = 0; i < rows; i++) {
        result[i] = 0.0; 
    }
    int *offsets = NULL;
    int count = offsetFunct(rows, cols, matrix, &offsets, x);

    space += sizeof(int)*count;

    if (count != -1)
    {


        start = clock();
        for (int i = 0; i < (rows + cols - 1); i++)
        {
            if (offsets[i] != x)
            { // Check for valid offsets
                // Print diagonal values based on offsets
                if (offsets[i] < 0)
                {                             // If the offset is negative
                    int row = 0 - offsets[i]; // Calculate starting row
                    int col = 0;
                    // Only check if row is within bounds since col starts at 0
                    while (row < rows)
                    {
                        if(matrix[row][col] != 0){
                            result[row] += matrix[row][col] * vector[col];
                        }
                        row++;
                        col++;
                    }
                }
                else
                {                         // If the offset is non-negative
                    int col = offsets[i]; // Calculate starting column
                    int row = 0;
                    // Only check if col is within bounds since row starts at 0
                    while (col < cols)
                    {
                        if(matrix[row][col] != 0){
                            result[row] += matrix[row][col] * vector[col];
                        }
                        col++;
                        row++;
                    }
                }
            }
        }
        end = clock();
    }

    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;

    if (count != -1)
    {


        start = clock();
        for (int i = 0; i < (rows + cols - 1); i++)
        {
            if (offsets[i] != x)
            { // Check for valid offsets
                // Print diagonal values based on offsets
                if (offsets[i] < 0)
                {                             // If the offset is negative
                    int row = 0 - offsets[i]; // Calculate starting row
                    int col = 0;
                    // Only check if row is within bounds since col starts at 0
                    while (row < rows)
                    {
                        if(matrix[row][col] != 0){

                        }
                        row++;
                        col++;
                    }
                }
                else
                {                         // If the offset is non-negative
                    int col = offsets[i]; // Calculate starting column
                    int row = 0;
                    // Only check if col is within bounds since row starts at 0
                    while (col < cols)
                    {
                        if(matrix[row][col] != 0){
                            
                        }
                        col++;
                        row++;
                    }
                }
            }
        }
        end = clock();
    }

    double cpu_time_used_without;
    

    cpu_time_used_without = ((double)(end - start)) / CLOCKS_PER_SEC;

    int intOperations=0;

    if (count != -1)
    {
        for (int i = 0; i < (rows + cols - 1); i++)
        {
            if (offsets[i] != x)
            { // Check for valid offsets
                // Print diagonal values based on offsets
                if (offsets[i] < 0)
                {                             // If the offset is negative
                    int row = 0 - offsets[i]; // Calculate starting row
                    int col = 0;
                    // Only check if row is within bounds since col starts at 0
                    while (row < rows)
                    {
                        space += sizeof(double);
                        if(matrix[row][col] != 0)
                            flops++;
                        row++;
                        col++;
                        intOperations++;
                    }
                }
                else
                {                         // If the offset is non-negative
                    int col = offsets[i]; // Calculate starting column
                    int row = 0;
                    // Only check if col is within bounds since row starts at 0
                    while (col < cols)
                    {
                        space += sizeof(double);
                        if(matrix[row][col] != 0)
                            flops++;
                        col++;
                        row++;
                        intOperations++;
                    }
                }
            }
        }
        // Free allocated memory
        free(offsets);
    }
    flops*=2;


    printf("\n\n\nTime used: %f seconds", cpu_time_used);
    printf("\n\n\nTime used: %f seconds", cpu_time_used_without);

    double gflops;

    gflops = (flops/cpu_time_used)/1e9;

    printf("\n\nFlops : %lf",flops);
        
    printf("\n\n\n%d int operations",intOperations + (rows + cols - 1));

    printf("\n\n\n%lf flops with the loop and all",gflops);

    gflops = (flops/(cpu_time_used - cpu_time_used_without))/1e9;

    printf("\n\n\n%lf flops without the loop and all only the svmp",gflops);
    
    printf("\n\n\n%d byte\n\n\n",space);

    for (int i = 0; i < 10; i++) {
        printf("%lf\t",result[i]);
    }


    free(vector);
    free(result);
}

int main(){
    FILE* file = fopen("matrixDataset/af23560.mtx", "r"); // file fetching
    if (file == NULL) {
        perror("Error opening file");
        return EXIT_FAILURE; // show a failure status if the directory is not there...
    }

    char line[256]; // a char array of max size to store each line
    while (fgets(line, sizeof(line), file) != NULL) {
        if (line[0] != '%' && line[0] != '\n') {
            break;
        }
    }

    int rows, cols, nonzeros;
    sscanf(line, "%d %d %d", &rows, &cols, &nonzeros);
    // Allocate the matrix dynamically



    double** matrix = (double**)malloc(rows * sizeof(double*));
    if (matrix == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        return EXIT_FAILURE;
    }
    for (int i = 0; i < rows; i++) {
        matrix[i] = (double*)malloc(cols * sizeof(double));
    	memset(matrix[i], 0, cols * sizeof(double));
        if (matrix[i] == NULL) {
            fprintf(stderr, "Memory allocation failed\n");
            return EXIT_FAILURE;
        }
    }

    printf("The no. of Rows: %d \nThe no. of Columns: %d \nThe No. of terms in this matrix: %d\n", rows, cols, nonzeros);

    int row, col;
    double value;
    while (fscanf(file, "%d %d %lf", &row, &col, &value) == 3) {
        if (row >= 1 && row <= rows && col >= 1 && col <= cols ) {
            matrix[row - 1][col - 1] = value;
        } else {
            fprintf(stderr, "Warning: Index out of bounds (%d, %d)\n", row, col);
        }
    }

    fclose(file);
    
    int x = rows + 5;


    compressedDIA(matrix, rows, cols, x);

    // Free the allocated memory
    for (int i = 0; i < cols; i++) {
        free(matrix[i]);
    }
    free(matrix);

    return 0;
}