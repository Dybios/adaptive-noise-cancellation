#ifndef MATRIX_H
#define MATRIX_H
#endif

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

// Function to print a matrix
void print_matrix(double *mat, int rows, int cols);

// Function to add two matrices
void add(double *mat1, double *mat2, double *result, int rows, int cols);

// Function to subtract two matrices
void subtract(double *mat1, double *mat2, double *result, int rows, int cols);

// Function to multiply two matrices
void multiply(double *mat1, double *mat2, double *result, int rows1, int cols1, int cols2);

// Function to find the determinant of a square matrix
double determinant(double *matrix, int size);

// Function to find transpose of a matrix
void transpose(double *matrix, double *result, int rows, int cols);

// Function to find inverse of a matrix
void inverse(double *matrix, double *inverted_matrix, int size);

// Function to find the covariance of a matrix
void covariance(double *matrix, double *covariance_mat, int rows, int cols);

// Function to find the mean of a matrix (row/column)
void mean(double *matrix, double *result, int rows, int cols, int flag);

