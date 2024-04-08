#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"

#define T 1024 // Length between Q complexes TODO: Use pan-tompkins to find this instead

#define N 10 // Averaging for state vector
#define ECG_LEADS 12

double x_hat[T][1];     // State vector
double Y[T][ECG_LEADS];     // measurements without coeffs
double P;           // State covariance
double Q;           // Process noise covariance
double R;           // Measurement noise covariance
double K[1][T];     // Kalman gain
double w_hat[T][1];     // Measurement vector

int main(int argc, char *argv[]) {
    const int ch = 1;
    const int frame = 128;
    const int shift = 32;

    // Initialize Q
    Q = 0.005;

    // Keep residual history across the streams
    double model_residual_avg[T][1];

    // Input files from args or use defaults
    FILE *input = NULL;
    FILE *output = NULL;
    if (argc == 1) {
       // Use defaults
       input = fopen("data.csv", "rb");
       output = fopen("data/clean_output.wav", "wb");
    }
    else if (argc > 1 && argc < 3) {
       printf("Not enough arguments.\n");
       exit(1);
    }
    else if (argc > 3) {
       printf("Too many arguments.\n");
       exit(1);
    }
    else {
       input = fopen(argv[1], "rb");
       output = fopen(argv[3], "wb");
    }

    if (input1 == NULL || input2 == NULL || output == NULL) {
           printf("Cannot open file.\n");
    }	
	
    // Get the file size
    fseek(input1, 0, SEEK_END);
    int file_size = ftell(input1);

    // get aligned buffer loop counter
    int rem = data_size % T;
    int loop_count = (file_size + T - rem) / T;
    printf("Loop count = %d\n", loop_count);

    // read the 12 lead ECG data to 12 input channels
    double *in1, *in2, *in3, *in4, *in5, *in6, *in7, *in8, *in9, *in10, *in11, *in12;
    // TODO: read 12 ecg lead data into each of in1-in12
    fread(in_short1, data_size, 2, input);
    fclose(input1);

   // Initial state estimate and covariance
   double x0[T][N];
   for (int row = 0; row < T; row++) {
       for (int col = 0; col < N; col++) {
           x0[row][col] = (double) in1[row + (col * T)];
       }
   }
   mean((double *)x0, (double *)x_hat, T, N, 1);

   double P_mat[1][1];
   covariance((double *)x_hat, (double *)P_mat, T, 1);
   P = P_mat[0][0];

   double y[T][1], y1[T][1], y2[T][1], y3[T][1],
          y4[T][1], y5[T][1], y6[T][1], y7[T][1],
		  y8[T][1], y9[T][1], y10[T][1], y11[T][1];

   // Run the Adaptive Kalman filter
   for (int count = 0; count < loop_count; count++) {
       // Initial measurement matrix
       // Put T samples to yk buffer
       for (int i = 0; i < T; i++) {
          y[i][0] = in1[i + (count * T)];
          y1[i][0] = in2[i + (count * T)];
          y2[i][0] = in3[i + (count * T)];
          y3[i][0] = in4[i + (count * T)];
          y4[i][0] = in5[i + (count * T)];
          y5[i][0] = in6[i + (count * T)];
          y6[i][0] = in7[i + (count * T)];
          y7[i][0] = in8[i + (count * T)];
          y8[i][0] = in9[i + (count * T)];
          y9[i][0] = in10[i + (count * T)];
          y10[i][0] = in11[i + (count * T)];
          y11[i][0] = in12[i + (count * T)];
       }

       // eq4
       double y_inv[ECG_LEADS][ECG_LEADS], Y_transpose[ECG_LEADS][T], Y_res[ECG_LEADS][ECG_LEADS], Y_res2[1][T], y_coeff_hat[ECG_LEADS][1];
       memcpy(Y, y1[0], sizeof(y1));
       memcpy(Y + 1, y2[0], sizeof(y2));
       memcpy(Y + 2, y3[0], sizeof(y3));
       memcpy(Y + 3, y4[0], sizeof(y4));
       memcpy(Y + 4, y5[0], sizeof(y5));
       memcpy(Y + 5, y6[0], sizeof(y6));
       memcpy(Y + 6, y7[0], sizeof(y7));
       memcpy(Y + 7, y8[0], sizeof(y8));
       memcpy(Y + 8, y9[0], sizeof(y9));
       memcpy(Y + 9, y10[0], sizeof(y10));
       memcpy(Y + 10, y11[0], sizeof(y11));
       memcpy(Y + 11, y12[0], sizeof(y12));

       transpose((double *)Y, (double *)Y_transpose, T, 1);
       multiply((double *)Y_transpose, (double *)Y, (double *)Y_res, ECG_LEADS, T, ECG_LEADS);
       inverse((double *)Y_res, (double *)y_inv, ECG_LEADS);
       for (int rows = 0; rows < ECG_LEADS; rows++) {
           for (int cols = 0; cols < T; cols++) {
               Y_res2[rows][cols] = y_inv[rows][rows] * Y_transpose[rows][cols];
           }
	   }
       multiply((double *) Y_res2, (double *)y, (double *)y_coeff_hat, ECG_LEADS, T, 1);

       double y_i_hat[T][1];
	   multiply((double *)Y, (double *)y_coeff_hat, (double *)y_i_hat, T, ECG_LEADS, 1);
       subtract((double *)y, (double *) y_i_hat, (double *)w_hat, T, 1); // eq3

       double R_mat[1][1];
       covariance((double *)w_hat, (double*)R_mat, T, 1);
       R = R_mat[0][0];


       double K_scalar = (P + Q) / (R + P + Q); // eq9
       // eq7
       double temp_sub[T][1], temp_mul[T][1], temp_x_hat[T][1];
       subtract((double *)y, (double *)x_hat, (double *)temp_sub, T, 1);
       for (int i = 0; i < T; i++) {
           temp_mul[i][0] = temp_sub[i][0] * K_scalar;
       }
       add((double *) x_hat, (double *)temp_mul, (double *)x_hat, T, 1);
       P = P + Q - (K_scalar * (P + Q)); // eq8

       // update process noise covariance for each incoming input
       double model_residual[T][1], residual_mean[T][1], residual_transpose[1][T], Q_coeff[1][1];
       subtract((double *)y, (double *)x_hat, (double *) model_residual, T, 1); // eq11
       add((double *)model_residual, (double *)model_residual_avg, (double *)residual_mean, T, 1);
       if (loop_count % N == 0) {
           for (int i = 0; i < T; i++) {
               residual_mean[i][0] = residual_mean[i][0] / N;
           }
       }
       memcpy(model_residual_avg, residual_mean, T);
       transpose((double *)residual_mean, (double *)residual_transpose, T, 1);
       multiply((double *)residual_transpose, (double *)residual_mean, (double *)Q_coeff, 1, T, 1);
       Q = fmax(0, (Q_coeff[0][0] / T) - P - R); // eq14

       // Append the output buffer to output file
       fseek(output, 0, SEEK_END);
       fwrite(x_hat, sizeof(double), 2, output);
   }

   free(in1);
   free(in2);
   free(in3);
   free(in4);
   free(in5);
   free(in6);
   free(in7);
   free(in8);
   free(in9);
   free(in10);
   free(in11);
   free(in12);
   fclose(output);
   return 0;
}
