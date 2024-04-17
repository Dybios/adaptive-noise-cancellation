#ifndef ECG_KALMAN_H
#define ECG_KALMAN_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#define N 5 // Averaging for state vector
#define ECG_LEADS 12

int preprocess_ecg_data(float **data, int rows, float **preprocessed_output, int *ecg_complex_length);
int process_kalman(float **data, int rows, int cols, int ecg_complex_length, float *output);

#ifdef __cplusplus
}
#endif

#endif // ECG_KALMAN_H
