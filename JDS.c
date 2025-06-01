#include <stdio.h>
#include <stdlib.h>
#include <time.h>

typedef struct {
    int row, col;
    double value;
} Entry;

// Function to compare entries based on diagonal index (col - row)
int compare(const void *a, const void *b) {
    Entry *e1 = (Entry *)a;
    Entry *e2 = (Entry *)b;
    int diag1 = e1->col - e1->row;
    int diag2 = e2->col - e2->row;
    return (diag1 == diag2) ? e1->row - e2->row : diag1 - diag2;
}

// Function to perform Sparse Matrix-Vector Multiplication
void spmv(int rows, int nonZeros, Entry *entries, double *x, double *y) {
    for (int i = 0; i < rows; i++) {
        y[i] = 0.0; // Initialize result vector
    }

    for (int i = 0; i < nonZeros; i++) {
        int row = entries[i].row;
        int col = entries[i].col;
        double value = entries[i].value;
        y[row] += value * x[col];
    }
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        printf("Usage: %s <file.mtx>\n", argv[0]);
        return 1;
    }

    FILE *file = fopen(argv[1], "r");
    if (!file) {
        perror("Error opening file");
        return 1;
    }

    // Skip comments in the header
    char line[256];
    do {
        fgets(line, sizeof(line), file);
    } while (line[0] == '%');

    // Read matrix dimensions and non-zero count
    int rows, cols, nonZeros;
    sscanf(line, "%d %d %d", &rows, &cols, &nonZeros);

    // Allocate memory for non-zero entries
    Entry *entries = (Entry *)malloc(nonZeros * sizeof(Entry));
    if (!entries) {
        perror("Memory allocation failed");
        fclose(file);
        return 1;
    }

    // Read non-zero entries
    for (int i = 0; i < nonZeros; i++) {
        int row, col;
        double value;
        fscanf(file, "%d %d %lf", &row, &col, &value);
        entries[i].row = row - 1; // Convert to 0-based index
        entries[i].col = col - 1;
        entries[i].value = value;
    }
    fclose(file);

    // Sort entries by diagonal (col - row) and then row
    qsort(entries, nonZeros, sizeof(Entry), compare);

    // Space calculation
    size_t matrixStorage = nonZeros * sizeof(Entry);
    size_t inputVectorStorage = cols * sizeof(double);
    size_t outputVectorStorage = rows * sizeof(double);
    size_t totalStorage = matrixStorage + inputVectorStorage + outputVectorStorage;

    printf("\nTotal Storage: %.6f MB\n", totalStorage / (1024.0 * 1024.0));
    
    // Clean up
    free(entries);
    return 0;
}
