#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#define numThread 1


void ourmethod(double** matrix, int rows, int cols) {
    clock_t start, end;
    double cpu_time_used;
    int nonZeroCount = 0;
    double gflops =0;

    int startRow = rows-1;
    int startCol = 0;
    int counterEle = 0;
    int i=0;
    int conju=0;
    int space = 0;
    omp_set_num_threads(numThread);
    
    while(startRow >= 0 && startCol < cols) {// 4 0
        i = 0;
        while(i <= counterEle) {

            if(matrix[startRow + i][startCol + i] != 0 && conju==0){
                nonZeroCount++;
                conju=1;
            }
            if(conju!=0){
                if(matrix[startRow + i][startCol + i] != 0){
                }else{
                    conju=0;
                }
            }

            
            i++;
        }
        if(conju){
            conju=0;
        }

        if (startRow > 0) {
            startRow--;
            counterEle++;
        } else {
            startCol++;
            counterEle--;
        }
    }

    // malloc to the first **pointer

    double** storeMatrix = (double**)malloc(nonZeroCount * sizeof(double*));
    int* rowData = (int*)malloc(nonZeroCount * sizeof(int));
    int* colData = (int*)malloc(nonZeroCount * sizeof(int));
    int* clusterLen = (int*)malloc(nonZeroCount * sizeof(int));



    space += 2*(nonZeroCount * sizeof(int));


    // malloc to the *pointer 
    startRow = rows-1;
    startCol = 0;
    counterEle = 0;
    i=0;
    conju=0;
    int nonZeroCount1=0;
    int valueCountPerCluster=0;
    
    while(startRow >= 0 && startCol < cols) {
        i = 0;
        while(i <= counterEle) {

            if(matrix[startRow + i][startCol + i] != 0 && conju==0){
                nonZeroCount1++;
                conju=1;
            }
            if(conju!=0){
                if(matrix[startRow + i][startCol + i] != 0){
                    valueCountPerCluster++;
                }else{
                    conju=0;
                    storeMatrix[nonZeroCount1-1] = (double*)malloc(valueCountPerCluster * sizeof(double));
                    space += valueCountPerCluster * sizeof(double);
                    clusterLen[nonZeroCount1-1] = valueCountPerCluster;
                    valueCountPerCluster=0;
                }
            }

            
            i++;
        }
        if(conju){
            conju=0;
            storeMatrix[nonZeroCount1-1] = (double*)malloc(valueCountPerCluster * sizeof(double));
            space += valueCountPerCluster * sizeof(double);
            clusterLen[nonZeroCount1-1] = valueCountPerCluster;
            valueCountPerCluster=0;
        }

        if (startRow > 0) {
            startRow--;
            counterEle++;
        } else {
            startCol++;
            counterEle--;
        }
    }

    // store the values inside the pointers which we have created
    nonZeroCount1=0;
    startRow = rows-1;
    startCol = 0;
    counterEle = 0;
    i=0;
    conju=0;
    valueCountPerCluster=0;

    while(startRow >= 0 && startCol < cols) {
        i = 0;
        while(i <= counterEle) {

            if(matrix[startRow + i][startCol + i] != 0 && conju==0){
                nonZeroCount1++;
                rowData[nonZeroCount1-1] = startRow+i+1;
                colData[nonZeroCount1-1] = startCol+i+1;
                // fprintf(storeDia,"%d %10.7e \n",(((i*rows)+j)+1),matrix[i][j]);
                conju=1;
            }
            if(conju!=0){
                if(matrix[startRow + i][startCol + i] != 0){
                    storeMatrix[nonZeroCount1-1][valueCountPerCluster]=matrix[startRow + i][startCol + i];
                    valueCountPerCluster++;
                }else{
                    conju=0;
                    valueCountPerCluster=0;
                }
            }


            
            i++;
        }
        if(conju){
            conju=0;
            valueCountPerCluster=0;
        }

        if (startRow > 0) {
            startRow--;
            counterEle++;
        } else {
            startCol++;
            counterEle--;
        }
    }


    // Sparse Matrix-Vector Multiplication (SpMV)
    double* vector = (double*)malloc((cols) * sizeof(double));

    for (int i = 0; i < cols; i++) {
        vector[i] = 1.0; 
    }



    double* result = (double*)malloc((rows) * sizeof(double));

    for (int i = 0; i < rows; i++) {
        result[i] = 0.0; 
    }

    start = clock();
    //clock start
    for(int i=0; i<nonZeroCount1;i++){
        for(int j=0; j<clusterLen[i]; j++){
            result[rowData[i]+j-1] += storeMatrix[i][j] * vector[colData[i]+j-1];
        }
    }
    end = clock();


    // flop count
    for(int i=0; i<nonZeroCount1;i++){
        for(int j=0; j<clusterLen[i]; j++){
            gflops++;
        }
    }
    gflops *= 2;

    printf("\n\nResult Matrix : \n\n");

    for (int i = 0; i < 5; i++) {
        printf("%f\n",result[i]);
    }
    

    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("\n\n\nTime used: %f seconds", cpu_time_used);
        
    printf("\n\n\n%lf flop Arithmatic operation count",gflops);

    gflops = (gflops/cpu_time_used)/1e9;

    printf("\n\n\n%lf gflops",gflops);
    
    printf("\n\n\n%d byte",space);


    // Parallelized SpMV using OpenMP

    start = clock();
    // #pragma omp parallel for
    // for (int i = 0; i < nonZeroCount1; i++) {
    //     for (int j = 0; j < clusterLen[i]; j++) {
    //         double x = result[rowData[i] + j - 1];
    //         double y = storeMatrix[i][j];
    //         double z = vector[colData[i] + j - 1];
    //         result[rowData[i] + j - 1] = x + (y * z);
    //     }
    // }
    #pragma omp parallel for
    for(int i=0; i<nonZeroCount1;i++){
        for(int j=0; j<clusterLen[i]; j++){
            result[rowData[i]+j-1] += storeMatrix[i][j] * vector[colData[i]+j-1];
        }
    }
    // int j=0;

    // int threadsGo = nonZeroCount1/numThread;

    // #pragma omp parallel for
    // for(i=0; i<nonZeroCount1;i++){
    //     for(j=0; j<clusterLen[i]; j++){
    //         int ff = colData[i];
    //         double stor = storeMatrix[i][j];
    //         double vect = vector[ff+j-1];
    //         result[rowData[i]+j-1] += stor * vect;
    //     }
    // }

    // #pragma omp parallel for schedule(static, nonZeroCount1)
    // for(i=0; i<nonZeroCount1;i++){
    //     for(j=0; j<clusterLen[i]; j++){
    //         int ff = colData[i];
    //         double stor = storeMatrix[i][j];
    //         double vect = vector[ff+j-1];
    //         result[rowData[i]+j-1] += stor * vect;
    //     }
    // }
    end = clock();
    
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("\n\n\nParallel computing time used: %f seconds", cpu_time_used);
    printf("\n\nResult Matrix : \n\n");

    for (int i = 0; i < 5; i++) {
        printf("%f\n",result[i]);
    }
    free(vector);
    free(result);
    i=0;
    while(i<nonZeroCount){
        free(storeMatrix[i]);
        i++;
    }
    free(rowData);
    free(colData);
    free(clusterLen);
}

int main() {
    

    FILE* file = fopen("/Users/apple/Documents/PROGRAMMING/c/matrixDataset/af23560.mtx", "r"); // file fetching
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


    ourmethod(matrix, rows, cols);

    // Free the allocated memory
    for (int i = 0; i < cols; i++) {
        free(matrix[i]);
    }
    free(matrix);

    return 0;
}

