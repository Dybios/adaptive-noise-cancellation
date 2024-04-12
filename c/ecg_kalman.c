#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "matrix.h"
#include "qrs_detector/c89/qrs_detection.h"
#include "filter-c/filter.h"

#define N 5 // Averaging for state vector
#define ECG_LEADS 12

static int T;

int main(int argc, char *argv[]) {
    const int ch = 1;
    const int frame = 128;
    const int shift = 32;
    int rows = 0, cols = 0;
    char line[ECG_LEADS * 100]; // Buffer for a line, assuming max 100 characters per field

    // Input files from args or use defaults
    FILE *input = NULL;
    FILE *output_preprocess = NULL;
    FILE *output = NULL;
    if (argc == 1) {
       // Use defaults
       input = fopen("../data/data_synthesized_0dB.csv", "r");
       output_preprocess = fopen("../data/output/preprocessed.csv", "w");
       output = fopen("../data/output/clean_output.csv", "w");
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

    if (input == NULL || output_preprocess == NULL || output == NULL) {
           printf("Cannot open file.\n");
    }

    // Count rows and columns
    while (fgets(line, sizeof(line), input) != NULL) {
        cols = 0;
        char *token = strtok(line, ",");
        while (token != NULL) {
            cols++;
            token = strtok(NULL, ",");
        }
        rows++; // data in each channel length
    }

    fseek(input, 0, SEEK_SET); // Move back to beginning of file

    // Allocate memory for the 2D array
    double **data = (double **)malloc(rows * sizeof(double *));
    for (int i = 0; i < rows; i++) {
        data[i] = (double *)malloc(cols * sizeof(double));
    }

    // Read the CSV data into the array
    int i = 0, j = 0;
    while (fgets(line, sizeof(line), input) != NULL) {
        j = 0;
        char *token = strtok(line, ",");
        while (token != NULL) {
            data[i][j++] = atof(token);
            token = strtok(NULL, ",");
        }
        i++;
    }
    fclose(input);

    // Example: Print the array contents
    printf("rows (data/channel) = %d\tcols (channels) = %d\n", rows, cols);

    // read the 12 lead ECG data to 12 input channels
    double in1[rows], in2[rows], in3[rows], in4[rows],
           in5[rows], in6[rows], in7[rows], in8[rows],
           in9[rows], in10[rows], in11[rows], in12[rows];
    for (int i = 0; i < rows; i++) {
        in1[i] = data[i][0];
        in2[i] = data[i][1];
        in3[i] = data[i][2];
        in4[i] = data[i][3];
        in5[i] = data[i][4];
        in6[i] = data[i][5];
        in7[i] = data[i][6];
        in8[i] = data[i][7];
        in9[i] = data[i][8];
        in10[i] = data[i][9];
        in11[i] = data[i][10];
        in12[i] = data[i][11];
    }

    // Preprocess the noisy ECG signal in1
    /* High pass filter with 0.5Hz cutoff */
    BWHighPass* hp_filter = create_bw_high_pass_filter(4, 360, 0.5);
    for (int i = 0; i < rows; i++) {
        in1[i] = bw_high_pass(hp_filter, in1[i]);
    }
    free_bw_high_pass(hp_filter);

    /* Low pass filter with 0.5Hz cutoff */
    BWLowPass* lp_filter = create_bw_low_pass_filter(4, 360, 90);
    for (int i = 0; i < rows; i++) {
        in1[i] = bw_low_pass(lp_filter, in1[i]);
    }
    free_bw_low_pass(lp_filter);

    /* Notch filter centered around 50Hz for powerline/baseline wander */
    BWBandPass* bp_filter = create_bw_band_pass_filter(4, 360, 4, 61);
    for (int i = 0; i < rows; i++) {
        in1[i] = bw_band_pass(bp_filter, in1[i]);
    }
    free_bw_band_pass(bp_filter);

    // Output for preprocessing stage
    for (int i = 0; i < rows; i++) {
       fprintf(output_preprocess, "%f%c", in1[i], '\n');
    }
    fprintf(output_preprocess, "\n");  // Add newline at the end of each row

    // Detect QRS to get the value of T
    int buff_size = 600;
    char result[buff_size];
    double sampling_rate = 360.00;
    T = DetectQrsPeaks(in1, rows, result, sampling_rate);
    T -= 0.02 * T;
    printf("T = %d\n", T);

    double x_hat[T][1];     // State vector
    double Y[T][ECG_LEADS];     // measurements without coeffs
    double P;           // State covariance
    double Q;           // Process noise covariance
    double R;           // Measurement noise covariance
    double K[1][T];     // Kalman gain
    double w_hat[T][1];     // Measurement vector

    // Initialize Q
    Q = 0.005;

    // Keep residual history across the streams
    double model_residual_avg[T][1];

    // get aligned buffer loop counter
    int rem = rows % T;
    int loop_count = (rows + T - rem) / T;
    printf("Loop count = %d\n", loop_count);

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
       for (int row = 0; row < T; row++) {
           Y[row][0] = y1[row][0];
           Y[row][1] = y1[row][1];
           Y[row][2] = y1[row][2];
           Y[row][3] = y1[row][3];
           Y[row][4] = y1[row][4];
           Y[row][5] = y1[row][5];
           Y[row][6] = y1[row][6];
           Y[row][7] = y1[row][7];
           Y[row][8] = y1[row][8];
           Y[row][9] = y1[row][9];
           Y[row][10] = y1[row][10];
           Y[row][11] = y1[row][11];
       }

       transpose((double *)Y, (double *)Y_transpose, T, ECG_LEADS);
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

//       printf("P=%f,Q=%f,R=%f\n", P, Q, R);

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

//       print_matrix((double*)x_hat, T, 1);
       // Append the output buffer to output file
       for (int i = 0; i < T; i++) {
           fprintf(output, "%f%c", x_hat[i][0], '\n');
       }
       fprintf(output, "\n");  // Add newline at the end of each row
   }


   // Free the allocated memory
   for (int i = 0; i < rows; i++) {
       free(data[i]);
   }
   free(data);

   fclose(output);
   return 0;
}
