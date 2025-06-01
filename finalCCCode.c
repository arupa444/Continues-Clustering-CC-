#include <stdio.h>     // Standard input/output library for functions like printf and file operations
#include <stdlib.h>    // Standard library for general utilities like memory allocation (malloc, free)
#include <time.h>

// Function to swap two integers (used in sorting algorithms)
// Parameters:
//   - *a: Pointer to the first integer
//   - *b: Pointer to the second integer
// Return Value: void (modifies the integers directly through the pointers)
void swap(int *a, int *b) {
    int temp = *a; // Store the value of *a in temp
    *a = *b;       // Assign the value of *b to *a
    *b = temp;     // Assign the value of temp (original *a) to *b
}

// Function to swap two doubles (used in sorting algorithms)
// Parameters:
//   - *a: Pointer to the first double
//   - *b: Pointer to the second double
// Return Value: void (modifies the doubles directly through the pointers)
void swapDoubllleee(double *a, double *b) {
    double temp = *a; // Store the value of *a in temp
    *a = *b;       // Assign the value of *b to *a
    *b = temp;     // Assign the value of temp (original *a) to *b
}

// Structure to store matrix data during diagonal traversal and sorting
struct StoreDiag {
    int* storeOffsets;  // Array to store the (col - row) offset for each element
    int* storeRow;      // Array to store the row indices of elements
    int* storeCol;      // Array to store the column indices of elements
    double* storeValues; // Array to store the values of elements
    int nonZeros;       // Number of non-zero elements
};

// Structure to store matrix data in the Clustered Column (CC) format
struct ourMethodStr {
    int rows;           // Number of rows in the matrix
    int cols;           // Number of columns in the matrix
    int nonZeros;       // Total number of non-zero elements
    int* clusterSizes;  // Array to store the size of each cluster
    int numOfClusters;  // Number of clusters
    int* startRowClus;  // Array to store the starting row index of each cluster
    int* startColClus;  // Array to store the starting column index of each cluster
    double* storeValues; // Array to store the values of elements, grouped by cluster
};

// Function to allocate memory for a StoreDiag structure
// Parameters:
//   - diagStorage: A pointer to a StoreDiag struct
//   - nonZeros: The number of non-zero elements to allocate space for
// Return Value:
//  - A pointer to the StoreDiag struct with allocated memory
struct StoreDiag* allocate_mem_diag(struct StoreDiag* diagStorage, int nonZeros) {
    diagStorage->storeOffsets = (int*)malloc(nonZeros * sizeof(int)); // Allocate memory for offsets
    diagStorage->storeRow = (int*)malloc(nonZeros * sizeof(int));     // Allocate memory for row indices
    diagStorage->storeCol = (int*)malloc(nonZeros * sizeof(int));     // Allocate memory for column indices
    diagStorage->storeValues = (double*)malloc(nonZeros * sizeof(double)); // Allocate memory for values
    diagStorage->nonZeros = nonZeros;                                // Store the number of non-zero elements
    return diagStorage; // Return the pointer to the modified struct
}

// Function to allocate memory for an ourMethodStr structure
// Parameters:
//   - CCStorage: A pointer to an ourMethodStr struct
//   - rows: The number of rows in the matrix
//   - cols: The number of columns in the matrix
//   - nonZeros: The number of non-zero elements to allocate space for
// Return Value:
//  - A pointer to the ourMethodStr struct with allocated memory
struct ourMethodStr* allocate_mem_cc(struct ourMethodStr* CCStorage, int rows, int cols, int nonZeros) {
    CCStorage->clusterSizes = (int*)malloc(1 * sizeof(int));         // Allocate memory for cluster sizes
    CCStorage->startRowClus = (int*)malloc(1 * sizeof(int));        // Allocate memory for starting row indices
    CCStorage->startColClus = (int*)malloc(1 * sizeof(int));        // Allocate memory for starting column indices
    CCStorage->storeValues = (double*)malloc(nonZeros * sizeof(double));  // Allocate memory for values
    CCStorage->cols = cols;                                        // Store the number of columns
    CCStorage->rows = rows;                                        // Store the number of rows
    CCStorage->nonZeros = nonZeros;                                 // Store the number of non-zero elements
    return CCStorage;                                              // Return the pointer to the modified struct
}

// Function to merge two sorted subarrays during merge sort
// Parameters:
//  - CCStorage: A pointer to the StoreDiag struct
//  - l: The left index of the subarray to be merged
//  - m: The middle index of the subarray to be merged
//  - r: The right index of the subarray to be merged
// Return Value: void (modifies the StoreDiag struct directly)
void merge(struct StoreDiag* CCStorage, int l, int m, int r) {
    int i, j, k;
    int n1 = m - l + 1; // Size of the left subarray
    int n2 = r - m;     // Size of the right subarray

    // Temporary arrays for left and right subarrays
    int* LOffsets = (int*)malloc(n1 * sizeof(int));
    int* LRows = (int*)malloc(n1 * sizeof(int));
    int* LCols = (int*)malloc(n1 * sizeof(int));
    double* LValues = (double*)malloc(n1 * sizeof(double));

    int* ROffsets = (int*)malloc(n2 * sizeof(int));
    int* RRows = (int*)malloc(n2 * sizeof(int));
    int* RCols = (int*)malloc(n2 * sizeof(int));
    double* RValues = (double*)malloc(n2 * sizeof(double));

    if (!LOffsets || !LRows || !LCols || !LValues || !ROffsets || !RRows || !RCols || !RValues) {
        fprintf(stderr, "Memory allocation failed in merge\n");
        exit(1); // Handle memory allocation failure
    }

    // Copy data to temporary arrays
    for (i = 0; i < n1; i++) {
        LOffsets[i] = CCStorage->storeOffsets[l + i];
        LRows[i] = CCStorage->storeRow[l + i];
        LCols[i] = CCStorage->storeCol[l + i];
        LValues[i] = CCStorage->storeValues[l + i];
    }
    for (j = 0; j < n2; j++) {
        ROffsets[j] = CCStorage->storeOffsets[m + 1 + j];
        RRows[j] = CCStorage->storeRow[m + 1 + j];
        RCols[j] = CCStorage->storeCol[m + 1 + j];
        RValues[j] = CCStorage->storeValues[m + 1 + j];
    }

    // Merge the temporary arrays back into the original arrays in CCStorage
    i = 0;
    j = 0;
    k = l;
    while (i < n1 && j < n2) {
        if (LOffsets[i] <= ROffsets[j]) {
            CCStorage->storeOffsets[k] = LOffsets[i];
            CCStorage->storeRow[k] = LRows[i];
            CCStorage->storeCol[k] = LCols[i];
            CCStorage->storeValues[k] = LValues[i];
            i++;
        }
        else {
            CCStorage->storeOffsets[k] = ROffsets[j];
            CCStorage->storeRow[k] = RRows[j];
            CCStorage->storeCol[k] = RCols[j];
            CCStorage->storeValues[k] = RValues[j];
            j++;
        }
        k++;
    }

    // Copy remaining elements of LOffsets[], if any
    while (i < n1) {
        CCStorage->storeOffsets[k] = LOffsets[i];
        CCStorage->storeRow[k] = LRows[i];
        CCStorage->storeCol[k] = LCols[i];
        CCStorage->storeValues[k] = LValues[i];
        i++;
        k++;
    }

    // Copy remaining elements of ROffsets[], if any
    while (j < n2) {
        CCStorage->storeOffsets[k] = ROffsets[j];
        CCStorage->storeRow[k] = RRows[j];
        CCStorage->storeCol[k] = RCols[j];
        CCStorage->storeValues[k] = RValues[j];
        j++;
        k++;
    }

    free(LOffsets); free(LRows); free(LCols); free(LValues); // Free temporary arrays for left subarray
    free(ROffsets); free(RRows); free(RCols); free(RValues); // Free temporary arrays for right subarray
}

// Function to perform merge sort on the StoreDiag structure
// Parameters:
//  - CCStorage: A pointer to the StoreDiag struct
//  - l: The left index of the subarray to be sorted
//  - r: The right index of the subarray to be sorted
// Return Value: void (modifies the StoreDiag struct directly)
void mergeSort(struct StoreDiag* CCStorage, int l, int r) {
    if (l < r) {
        int m = l + (r - l) / 2; // Calculate middle index

        mergeSort(CCStorage, l, m);      // Recursively sort the left subarray
        mergeSort(CCStorage, m + 1, r);  // Recursively sort the right subarray

        merge(CCStorage, l, m, r);       // Merge the sorted subarrays
    }
}

// Function to sort the StoreDiag structure using merge sort
// Parameters:
//   - diagStorage: A pointer to a StoreDiag struct
// Return Value: void (modifies the StoreDiag struct directly)
void sortOurMethodStr(struct StoreDiag* diagStorage) {
    mergeSort(diagStorage, 0, diagStorage->nonZeros - 1); // Sort based on offsets
}

// Function to print an integer array to the console
// Parameters:
//   - arr: The integer array to be printed
//   - size: The size of the array
// Return Value: void (prints to the console)
void printArray(int arr[], int size) {
    for (int i = 0; i < size; i++) {
        printf("%d ", arr[i]); // Print each element followed by a space
    }
    printf("\n"); // Print a newline character at the end
}

// Function to print a double array to the console
// Parameters:
//   - arr: The double array to be printed
//   - size: The size of the array
// Return Value: void (prints to the console)
void printArrayValDouble(double arr[], int size) {
    for (int i = 0; i < size; i++) {
        printf("%lf ", arr[i]);  // Print each element followed by a space
    }
    printf("\n");  // Print a newline character at the end
}

// Function to perform clustering and sparse matrix-vector multiplication
// Parameters:
//   - diagStorage: A pointer to the StoreDiag struct (original data)
//   - CC: A pointer to the ourMethodStr struct (where the clustered data will be stored)
// Return Value: void (modifies the ourMethodStr struct and prints results)
void cc(struct StoreDiag* diagStorage, struct ourMethodStr* CC){
    CC->numOfClusters =1; // Initialize the number of clusters to 1
    for(int i=1; i< diagStorage->nonZeros; i++ ){  // Iterate through the non-zero elements
         //Check to see if a new cluster starts or not based on row and col values
        if((diagStorage->storeRow[i-1] != diagStorage->storeRow[i]-1) || (diagStorage->storeCol[i]-1 != diagStorage->storeCol[i-1])){
            CC->numOfClusters++; // Increment cluster count if a new cluster is found
        }
    }
    
    CC->clusterSizes = (int*)malloc(CC->numOfClusters * sizeof(int));  // Allocate memory for cluster sizes
    CC->startRowClus = (int*)malloc(CC->numOfClusters * sizeof(int));   // Allocate memory for starting row indices
    CC->startColClus = (int*)malloc(CC->numOfClusters * sizeof(int));   // Allocate memory for starting column indices
    
     if (!CC->clusterSizes || !CC->startRowClus || !CC->startColClus) {
         fprintf(stderr, "Memory allocation failed in cc for cluster info.\n"); // Memory check
        exit(1);
    }
    CC->storeValues = (double*)malloc(CC->nonZeros * sizeof(double));    // Allocate memory for matrix values
    if (!CC->storeValues) {
         fprintf(stderr, "Memory allocation failed in cc for storeValues.\n");
        exit(1);
    }

    CC->clusterSizes[0] = 1;           // Initialize size of first cluster
    CC->startRowClus[0] = diagStorage->storeRow[0]; // Initialize starting row of first cluster
    CC->startColClus[0] = diagStorage->storeCol[0]; // Initialize starting column of first cluster
    int jump =0;
    for(int i=1; i< diagStorage->nonZeros; i++ ){    // Iterate through elements to assign to clusters
        if((diagStorage->storeRow[i-1] != diagStorage->storeRow[i]-1) || (diagStorage->storeCol[i]-1 != diagStorage->storeCol[i-1])){
            ++jump; // Increment jump variable for new clusters
            CC->clusterSizes[jump] = 1;             // Initialize size of new cluster
            CC->startRowClus[jump] = diagStorage->storeRow[i]; // Initialize starting row of new cluster
            CC->startColClus[jump] = diagStorage->storeCol[i]; // Initialize starting column of new cluster
        }else{
            CC->clusterSizes[jump]++; // Increment size of current cluster if next element belongs to the same cluster
        }
    }
     for(int i = 0; i< CC->nonZeros; i++){ // Copy the values to the CC struct
            CC->storeValues[i] = diagStorage->storeValues[i];
        }

    // spmv (Sparse Matrix-Vector Multiplication)
    double* vector = (double*)malloc((CC->cols) * sizeof(double)); // Allocate memory for the input vector
      if (!vector){
          fprintf(stderr, "Memory allocation failed in cc for vector.\n");
          exit(1);
    }
    for (int i = 0; i < CC->cols; i++) {
        vector[i] = 1.0; // Initialize input vector with 1.0
    }

    double* result = (double*)malloc((CC->rows) * sizeof(double)); // Allocate memory for result vector
      if (!result){
        fprintf(stderr, "Memory allocation failed in cc for result.\n");
          exit(1);
    }
    for (int i = 0; i < CC->rows; i++) {
        result[i] = 0.0;  // Initialize result vector with 0.0
    }
    clock_t start = clock();
    int drive=0;
    for(int i=0; i<CC->numOfClusters;i++){
        for (int j = 0; j < CC->clusterSizes[i]; j++){
            result[CC->startRowClus[i]+j-1] += CC->storeValues[drive] * vector[CC->startColClus[i]+j-1]; // Perform multiplication and addition
            drive++;
        }
    }
    clock_t end = clock();
 
    double time_required = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Time taken for SPMV_CC is %lf\n",time_required);
    printf("Num of clusters : %d\n",CC->numOfClusters);
    double space = (2 * sizeof(int) * CC->numOfClusters) + (sizeof(double) * CC->nonZeros);
    printf("Space occupied : %lfMB",space/(1024*1024));
    //printArrayValDouble(result, CC->rows); // Print the result vector

    // printArrayValDouble(result, CC->rows);

    free(vector);  // Free the memory used by the vector
    free(result); // Free the memory used by the result vector
}

// Main function: Entry point of the program
// Parameters: void
// Return Value: int (0 for success, non-zero for failure)
int main() {
    FILE* file = fopen("matrixDataset/crystk02.mtx", "r"); // Open matrix data file
     if (file == NULL) { 
        perror("Error opening file");  // Print an error if the file cannot be opened
        return EXIT_FAILURE; // Return failure status
    }
    char line[256]; // Character array to store each line

    // Skip comments and empty lines in the input file to retrieve the matrix dimensions
    while (fgets(line, sizeof(line), file) != NULL) {
        if (line[0] != '%' && line[0] != '\n') {
            break;
        }
    }

    int rows, cols, nonZeros;  // Variables to store matrix dimensions and number of non-zero elements
    sscanf(line, "%d %d %d", &rows, &cols, &nonZeros);  // Read dimensions and non-zero count from the line

     struct StoreDiag* diagStorage = (struct StoreDiag*)malloc(sizeof(struct StoreDiag)); // Allocate memory for storing diagonal data
     if(diagStorage == NULL){
          fprintf(stderr, "Memory allocation for diagStorage failed\n");
        fclose(file);
        return EXIT_FAILURE;
    }
    diagStorage = allocate_mem_diag(diagStorage, nonZeros); // Allocate memory for storing diagonal data

    printf("The no. of Rows: %d \nThe no. of Columns: %d \nThe No. of NonZeros in this matrix: %d\n", rows, cols, nonZeros); // Print dimensions


    int row, col;  // Variable for row index and column index
    double value;  // Variable for the value at the specified index
    int i=0;
    // Read the matrix data from file and store it in the diagStorage struct
    while (fscanf(file, "%d %d %lf", &row, &col, &value) == 3) {
        if (row >= 1 && row <= rows && col >= 1 && col <= cols ) { // Check index to ensure within bounds
            diagStorage->storeRow[i] = row;    // Store the row index in storeRow array
            diagStorage->storeCol[i] = col;    // Store the column index in storeCol array
            diagStorage->storeValues[i] = value; // Store the value in storeValues array
            diagStorage->storeOffsets[i++] = col-row;    // Store the offset in storeOffsets array and increment the counter
        } else {
            fprintf(stderr, "Warning: Index out of bounds (%d, %d)\n", row, col);
        }
    }

    fclose(file); // Close the matrix file

    //   printf("Original array: \n"); // Print unsorted data
    // printArray(diagStorage->storeOffsets, diagStorage->nonZeros);
    // printArray(diagStorage->storeCol, diagStorage->nonZeros);
    // printArray(diagStorage->storeRow, diagStorage->nonZeros);
    // printArrayValDouble(diagStorage->storeValues, diagStorage->nonZeros);

    sortOurMethodStr(diagStorage); // Call merge sort instead to sort the data based on offsets

    //   printf("Sorted array with all sorted: \n"); // Print sorted data
    // printArray(diagStorage->storeOffsets, diagStorage->nonZeros);
    // printArray(diagStorage->storeCol, diagStorage->nonZeros);
    // printArray(diagStorage->storeRow, diagStorage->nonZeros);
    // printArrayValDouble(diagStorage->storeValues, diagStorage->nonZeros);
     // Allocate memory for the ourMethodStr structure
    struct ourMethodStr* CCStorage = (struct ourMethodStr*)malloc(sizeof(struct ourMethodStr));
        if (CCStorage == NULL) {
        fprintf(stderr, "Memory allocation for CCStorage failed\n");
          free(diagStorage->storeOffsets);
        free(diagStorage->storeCol);
        free(diagStorage->storeRow);
        free(diagStorage->storeValues);
         free(diagStorage);
        return EXIT_FAILURE;
    }
    CCStorage = allocate_mem_cc(CCStorage,rows, cols, nonZeros);  // Allocate memory

     CCStorage->rows = rows;     // Store the number of rows
     CCStorage->cols = cols;     // Store the number of columns
    CCStorage->nonZeros = nonZeros; // Store the number of non-zero elements
    cc(diagStorage,CCStorage); // Call the cc function for clustering and spmv

    // printf("\nCluster Information:\n"); // Print cluster information
    // for (int i = 0; i < CCStorage->numOfClusters; i++) {
    //     printf("  Cluster %d:\n", i + 1); // Print cluster number
    //     printf("    Size: %d\n", CCStorage->clusterSizes[i]);  // Print the size of the cluster
    //     printf("    Start Row: %d\n", CCStorage->startRowClus[i]);  // Print the starting row of the cluster
    //     printf("    Start Col: %d\n", CCStorage->startColClus[i]); // Print the starting column of the cluster
    // }

    // Free the allocated memory
    free(diagStorage->storeOffsets);
    free(diagStorage->storeCol);
    free(diagStorage->storeRow);
    free(diagStorage->storeValues);
    free(diagStorage);
    free(CCStorage->clusterSizes);
    free(CCStorage->startRowClus);
    free(CCStorage->startColClus);
    free(CCStorage->storeValues);
    free(CCStorage);
    // return....
    return 0; // Return success
}