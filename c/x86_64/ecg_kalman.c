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
    double sampling_rate = 360.00;
    int rows = 0, cols = 0;
    char line[ECG_LEADS * 100];

    /** This section reads the files into the program memory.
     *  Input can be provided by either cmdline arguments or if run as is, defaults are taken.
     */
    FILE *input = NULL;
    FILE *output_preprocess = NULL;
    FILE *output = NULL;
    if (argc == 1) {
       // Use defaults
       input = fopen("../../data/data_synthesized_0dB_64k.csv", "r");
       output_preprocess = fopen("../../data/output/preprocessed.csv", "w");
       output = fopen("../../data/output/clean_output.csv", "w");
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
       input = fopen(argv[1], "r");
       output = fopen(argv[2], "w");
       output_preprocess = fopen("../../data/output/preprocessed.csv", "w");
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

    /** Once the data is stored in array, sort the data by channel onto individal arrays.
     *  These are later used as inputs for filter operation.
     */
    // read the 12 lead ECG data to 12 input channels
    double *in1, *in2, *in3, *in4, *in5, *in6,
           *in7, *in8, *in9, *in10, *in11, *in12;
    in1 = (double *) malloc(sizeof(double[rows]));
    in2 = (double *) malloc(sizeof(double[rows]));
    in3 = (double *) malloc(sizeof(double[rows]));
    in4 = (double *) malloc(sizeof(double[rows]));
    in5 = (double *) malloc(sizeof(double[rows]));
    in6 = (double *) malloc(sizeof(double[rows]));
    in7 = (double *) malloc(sizeof(double[rows]));
    in8 = (double *) malloc(sizeof(double[rows]));
    in9 = (double *) malloc(sizeof(double[rows]));
    in10 = (double *) malloc(sizeof(double[rows]));
    in11 = (double *) malloc(sizeof(double[rows]));
    in12 = (double *) malloc(sizeof(double[rows]));

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

    /** Preprocessing Stage (adis300)
     *  The input ECG data is passed through a series of filters to remove noise, artifacts and baseline wander.
     *  Filters implemented here:
     *       High Pass Filter - 90Hz
     *       Low Pass Filter - 0.5Hz
     *       Band Pass Filter - 4 to 61Hz
     *       # NOTE: The paper only states the usage of a notch filter centered around 50/60Hz to remove powerline
     *         noise with a width relatively small, but wide enough to account for fluctuations in grid. In this program,
     *         a good performance was observed with a bandpass filter with a low half-pass of 4Hz and a high of 61Hz.
     */

    /* High pass filter with 0.5Hz cutoff */
    BWHighPass* hp_filter = create_bw_high_pass_filter(4, (int)sampling_rate, 0.5);
    for (int i = 0; i < rows; i++) {
        in1[i] = bw_high_pass(hp_filter, in1[i]);
    }
    free_bw_high_pass(hp_filter);

    /* Low pass filter with 0.5Hz cutoff */
    BWLowPass* lp_filter = create_bw_low_pass_filter(4, (int)sampling_rate, 90);
    for (int i = 0; i < rows; i++) {
        in1[i] = bw_low_pass(lp_filter, in1[i]);
    }
    free_bw_low_pass(lp_filter);

    /* Notch filter centered around 50Hz for powerline/baseline wander */
    BWBandPass* bp_filter = create_bw_band_pass_filter(4, (int)sampling_rate, 4, 61);
    for (int i = 0; i < rows; i++) {
        in1[i] = bw_band_pass(bp_filter, in1[i]);
    }
    free_bw_band_pass(bp_filter);

    // Output for preprocessing stage
    for (int i = 0; i < rows; i++) {
       fprintf(output_preprocess, "%f%c", in1[i], '\n');
    }
    fprintf(output_preprocess, "\n");  // Add newline at the end of each row

    /** QRS length detection for finding the value of T. (kosachevds)
     *  In this paper, demarcation of ECG complex was defined as per 120% of the mean interval between consecutive
     *  consecutive heartbeats found using QRS complex detection. While the paper has used PCA for QRS detection,
     *  we found that the Pan-Tompkins algorithm is also sufficiently accurate in detecting the R waves with an error
     *  of about 2%. For this project, the length between all R peaks is averaged and a 2% error is factored into it
     *  for allowing overlap between ECG signals.
     */
    char result[rows]; // Dummy as we do not need this.
    T = DetectQrsPeaks(in1, rows, result, sampling_rate);
    T -= 0.02 * T;
    printf("T = %d\n", T);

    /* Initialize all Kalman filter parameters here. */
    double x_hat[T][1];     // State vector
    double Y[T][ECG_LEADS - 1];     // measurements without coeffs
    double P;           // State covariance
    double Q;           // Process noise covariance
    double R;           // Measurement noise covariance
    double K_scalar;     // Kalman gain
    double w_hat[T][1];     // Measurement vector

    // Initialize Q (Not explicitely needed)
    Q = 0.001;

    // Keep residual history across the looping filter
    double model_residual_avg[T][1];

    // get aligned buffer loop counter
    int rem = rows % T;
    int loop_count = (rows + T - rem) / T;
    printf("Loop count = %d\n", loop_count);

    // Initial state estimate (x_hat) and state covariance (P/psi) by taking average of the first N ECG complexes.
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

    /** Adaptive Kalman Filter operations based on Vullings, Vries & Bergmans' equations begin from here.
     *  Equations used here are directly based of the flow sequence described by Vullings, Vries & Bergmans and
     *  mentioned in the supported project document. I nstead of using element-wise operations, all
     *  equations are directly operated on in their matrix form and outputs are handled as scalars internally if needed.
     */
    double y[T][1], y1[T][1], y2[T][1], y3[T][1],
          y4[T][1], y5[T][1], y6[T][1], y7[T][1],
          y8[T][1], y9[T][1], y10[T][1], y11[T][1];

    for (int count = 0; count < loop_count; count++) {
       /* Measurement matrix (ŵ) and Measurement covariance (R/Σ) */
       // Get T samples into separate yk buffers
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

       /** EQ4: Coefficients of estimate γ̂ = [Yt * Y]^−1 * Yt * y(i)
        *           where t denotes tranpose.
        *  Subscript -i ommitted for clarity.
        */
       double y_inv[ECG_LEADS][ECG_LEADS], Y_transpose[ECG_LEADS][T], Y_res[ECG_LEADS][ECG_LEADS], Y_res2[1][T], y_coeff_hat[ECG_LEADS][1];
       for (int row = 0; row < T; row++) {
           Y[row][0] = y1[row][0];
           Y[row][1] = y2[row][0];
           Y[row][2] = y3[row][0];
           Y[row][3] = y4[row][0];
           Y[row][4] = y5[row][0];
           Y[row][5] = y6[row][0];
           Y[row][6] = y7[row][0];
           Y[row][7] = y8[row][0];
           Y[row][8] = y9[row][0];
           Y[row][9] = y10[row][0];
           Y[row][10] = y11[row][0];
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

       /** EQ3: Measurement noise vector for ith signal, ŵ(i) = y(i) − ŷ(i)
        *           where, y = ith measured value of ECG complex
        *                  ŷ = ith estimated value of ECG complex
        */
       double y_i_hat[T][1];
       multiply((double *)Y, (double *)y_coeff_hat, (double *)y_i_hat, T, ECG_LEADS, 1);
       subtract((double *)y, (double *) y_i_hat, (double *)w_hat, T, 1);

       /* Updating measurement noise covariance R/sigma */
       double R_mat[1][1];
       covariance((double *)w_hat, (double*)R_mat, T, 1);
       R = R_mat[0][0];

       /** EQ9: Kalman Gain for (k+1)th reading, K(k + 1) = [Ψ(k) + Λ(k)] / [Σ(k +1) + Ψ(k) + Λ(k)]
        *  The above equation assumes an initial value of K available at k = 0. The first loop of the filter operation
        *  provides the initial values for w_hat, R/Σ and K for k=0. This can then be used as inputs to the next interation
        *  which can also be referred to as the update sequence of the filter.
        *
        *  NOTE: K is considered as a scalar, since P/Ψ, Q/Ψ and R/Σ essentially become scalar matrices. This is an implcit assumption
        *  that process and measurement noise are spatially uncorrelated. More details provided in the original paper.
        */
       K_scalar = (P + Q) / (R + P + Q);

       /** EQ7: The next (k+1)th state estimate, x̂(k + 1) = x̂(k) + [K(k + 1) * [y(k + 1) − x̂(k)]]
        *  This updates the state estimate using the noisy input and the Kalman Gain previously calculated.
        */
       double temp_sub[T][1], temp_mul[T][1];
       subtract((double *)y, (double *)x_hat, (double *)temp_sub, T, 1);
       for (int i = 0; i < T; i++) {
           temp_mul[i][0] = temp_sub[i][0] * K_scalar;
       }
       add((double *) x_hat, (double *)temp_mul, (double *)x_hat, T, 1);

       /** EQ8: The (k+1)th state estimate covariance, Ψ(k + 1) = [Ψ(k) + Λ(k)] − [K(k + 1) * [Ψ(k) + Λ(k)]]
        */
       P = P + Q - (K_scalar * (P + Q));

       /** EQ11: Model residual for (k+1)th heartbeat is defined to be:
        *             ρ(k + 1) = y(k + 1) − E[y(k +1) | y(k), Λ(k), Σ(k)]
        *         or, ρ(k + 1) = y(k + 1) − x̂(k)
        */
       double model_residual[T][1], residual_mean[T][1], residual_transpose[1][T], Q_coeff[1][1];
       subtract((double *)y, (double *)x_hat, (double *) model_residual, T, 1);
       add((double *)model_residual, (double *)model_residual_avg, (double *)residual_mean, T, 1);
       if (loop_count % N == 0) {
           for (int i = 0; i < T; i++) {
               residual_mean[i][0] = residual_mean[i][0] / N;
           }
       }
       memcpy(model_residual_avg, residual_mean, T);
       transpose((double *)residual_mean, (double *)residual_transpose, T, 1);
       multiply((double *)residual_transpose, (double *)residual_mean, (double *)Q_coeff, 1, T, 1);

       /** EQ14: The process noise covariance, λ̂ in its scalar form is given as,
        *           λ̂(k) = [(1/T) * ρt(k + 1) * ρ(k + 1)] − ψ(k) − σ(k + 1)  ; if positive,
        *     and   λ̂(k) = 0                                                 ; otherwise.
        */
       Q = fmax(0, (Q_coeff[0][0] / T) - P - R);

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
