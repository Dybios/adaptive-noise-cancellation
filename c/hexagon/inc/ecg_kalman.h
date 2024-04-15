#ifndef ECG_KALMAN_H
#define ECG_KALMAN_H

#ifdef __cplusplus
extern "C" {
#endif

int ecg_kalman_main(int nErr, double **data, int rows, int cols, double *output);

#ifdef __cplusplus
}
#endif

#endif // ECG_KALMAN_H
