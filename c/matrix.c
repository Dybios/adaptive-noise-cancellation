#include "matrix.h"
//#define DEBUG

// Function to add two matrices
void add(double *mat1, double *mat2, double *result, int rows, int cols) {
    for(int i=0; i<rows; i++) {
        for(int j=0; j<cols; j++) {
            *(result + i*cols + j) = *(mat1 + i*cols + j) + *(mat2 + i*cols + j);
        }
    }
}

// Function to subtract two matrices
void subtract(double *mat1, double *mat2, double *result, int rows, int cols) {
    for(int i=0; i<rows; i++) {
        for(int j=0; j<cols; j++) {
            *(result + i*cols + j) = *(mat1 + i*cols + j) - *(mat2 + i*cols + j);
        }
    }
}

// Function to multiply two matrices
void multiply(double *mat1, double *mat2, double *result, int rows1, int cols1, int cols2) {
    for(int i=0; i<rows1; i++) {
        for(int j=0; j<cols2; j++) {
            *(result + i*cols2 + j) = 0;
            for(int k=0; k<cols1; k++) {
                *(result + i*cols2 + j) += (*(mat1 + i*cols1 + k)) * (*(mat2 + k*cols2 + j));
            }
        }
    }
}

// Function to print a matrix
void print_matrix(double *mat, int rows, int cols) {
    for(int i=0; i<rows; i++) {
        for(int j=0; j<cols; j++) {
            printf("%f ", *(mat + i*cols + j));
        }
        printf("\n");
    }
}

// Function to get determinant of matrix
double determinant(double *mat, int n)
{
#ifdef DEBUG
    printf("Incoming Matrix To Determinant\n");
    print_matrix(mat, n, n);
#endif

    double det = 1.0;
    int sign = 1;

    // Create a local copy of the mat
    double temp[n][n];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            temp[i][j] = *(mat + i*n + j);
        }
    }

    // Perform Gaussian elimination to convert the matrix to upper triangular form
    for (int i = 0; i < n; i++) {
        // Find the pivot element and swap rows if necessary
        if (temp[i][i] == 0.0) {
            int k = i + 1;
            while (k < n && temp[k][i] == 0.0) {
                k++;
            }
            if (k == n) {
                return 0.0;  // Matrix is singular, determinant is 0

            }

            for (int j = i; j < n; j++) {
                double val = temp[i][j];
                temp[i][j] = temp[k][j];
                temp[k][j] = val;
            }

            sign = -sign;
        }

        // Perform row operations to make the elements below the pivot element zero
        for (int j = i + 1; j < n; j++) {
            double ratio = temp[j][i] / temp[i][i];
            for (int k = i; k < n; k++) {
                temp[j][k] -= ratio * temp[i][k];
            }
        }
    }

    // Calculate the determinant as the product of diagonal elements of the upper triangular matrix
    for (int i = 0; i < n; i++) {
        det *= temp[i][i];
    }

    return sign * det;
}

// Function to get transpose of a matrix
void transpose(double *matrix, double *result, int rows, int cols)
{
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            *(result + j*rows + i) = *(matrix + i*cols + j);
        }
    }
}

// Function to calculate the adjoint matrix
void adjoint(double *matrix, double *adjointMatrix, int n) {
#ifdef DEBUG
    printf("Incoming original to adjoint\n");
    print_matrix(matrix, n, n);
#endif

    double temp[n-1][n-1];
    int sign = 1;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            // Get the submatrix of matrix[i][j]
            int sub_i = 0, sub_j = 0;
            for (int row = 0; row < n; row++) {
                for (int col = 0; col < n; col++) {
                    if (row != i && col != j) {
                        temp[sub_i][sub_j++] = *(matrix + row*n + col);
                        if (sub_j == n - 1) {
                            sub_j = 0;
                            sub_i++;
                        }
                    }
                }
            }

            // Sign of adj[j][i] positive if sum of row and column indexes is even

            sign = ((i + j) % 2 == 0) ? 1 : -1;

            // Calculate the determinant of the submatrix

            *(adjointMatrix + j*n + i) = sign * determinant((double *)temp, n-1);
        }
    }

#ifdef DEBUG
    printf("Adjoint Matrix\n");
    print_matrix(adjointMatrix, n, n);
#endif
}

// Function to calculate the inverse of a matrix
void inverse(double *matrix, double *inverseMatrix, int n) {
#ifdef DEBUG
    printf("Incoming Matrix\n");
    print_matrix(matrix, n, n);
#endif

    double det = determinant(matrix, n);
#ifdef DEBUG
    printf("Determinant = %f", det);
#endif

    double adjointMatrix[n][n];

    if (det == 0) {
        return;
    }

    // adjoint
    adjoint(matrix, (double *)adjointMatrix, n);

    // Invert
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            *(inverseMatrix + i*n + j) = adjointMatrix[i][j] / det;
        }
    }
}

// Function to calculate the covariance of a matrix
void covariance(double *matrix, double *covariance_mat, int rows, int cols) {
    double covariance[rows][cols];
    double mean[rows] = {0};

    // Calculate column-wise mean
    for (int j = 0; j < cols; j++) {
        for (int i = 0; i < rows; i++) {
            mean[j] += *(matrix + i*cols + j);
        }
        mean[j] /= rows;
    }

    // Calculate covariance
    for (int j = 0; j < cols; j++) {
        for (int k = 0; k < cols; k++) {
            double cov = 0;
            for (int i = 0; i < rows; i++) {
                cov += (*(matrix + i*cols + j) - mean[j]) * (*(matrix + i*cols + k) - mean[k]);
            }
            *(covariance_mat + j*cols + k) = cov / (rows - 1);
        }
    }
}


// Function to calculate the mean of a matrix; 1 = column mean, 2 = row mean
void mean(double *matrix, double *result, int rows, int cols, int flag) {
    if (flag == 1) {
        // Calculate the sum of the elements in a row
        for (int row_count = 0; row_count < rows; row_count++) {
            double sum = 0;
            for (int col_count = 0; col_count < cols; col_count++) {
                sum += *(matrix + row_count * cols + col_count);
            }

            // Calculate the mean
            *(result + row_count) = sum / cols;
        }
    }
    else if (flag == 2) {
        // Calculate the sum of the elements in the column
        for (int col_count = 0; col_count < cols; col_count++) {
            double sum = 0;
            for (int row_count = 0; row_count < rows; row_count++) {
                sum += *(matrix + row_count * cols + col_count);
            }

            // Calculate the mean
            *(result + col_count) = sum / rows;
        }
    }
}
