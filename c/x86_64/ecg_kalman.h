#ifndef ECG_KALMAN_H
#define ECG_KALMAN_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#define N 5 // Averaging for state vector
#define ECG_LEADS 12

int preprocess_ecg_data(double **data, int rows, double **preprocessed_output, int *ecg_complex_length);
int process_kalman(double **data, int rows, int cols, int ecg_complex_length, double *output);

#endif // ECG_KALMAN_H
