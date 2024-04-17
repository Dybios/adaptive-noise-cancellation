#ifndef MATRIX_H
#define MATRIX_H
#endif

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

// Function to print a matrix
void print_matrix(float *mat, int rows, int cols);

// Function to add two matrices
void add(float *mat1, float *mat2, float *result, int rows, int cols);

// Function to subtract two matrices
void subtract(float *mat1, float *mat2, float *result, int rows, int cols);

// Function to multiply two matrices
void multiply(float *mat1, float *mat2, float *result, int rows1, int cols1, int cols2);

// Function to find the determinant of a square matrix
float determinant(float *matrix, int size);

// Function to find transpose of a matrix
void transpose(float *matrix, float *result, int rows, int cols);

// Function to find inverse of a matrix
void inverse(float *matrix, float *inverted_matrix, int size);

// Function to find the covariance of a matrix
void covariance(float *matrix, float *covariance_mat, int rows, int cols);

// Function to find the mean of a matrix (row/column)
void mean(float *matrix, float *result, int rows, int cols, int flag);

