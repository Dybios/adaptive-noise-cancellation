#include "matrix.h"
#include "qrs_detector/c89/qrs_detection.h"
#include "filter-c/filter.h"

#include "ecg_kalman.h"

/** Preprocessing Stage (adis300)
 *  The input ECG data is passed through a series of filters to remove noise, artifacts and baseline wander.
 *  Filters implemented here:
 *       High Pass Filter - 90Hz
 *       Low Pass Filter - 0.5Hz
 *       Band Stop Filter - 49 to 71Hz
 *       # NOTE: The paper only states the usage of a notch filter centered around 50/60Hz to remove powerline
 *         noise with a width relatively small, but wide enough to account for fluctuations in grid. In this program,
 *         a good performance was observed with a bandstop/notch filter with a low half-pass of 49Hz and a high of 71Hz.
 */
int preprocess_ecg_data(double **data, int rows, double **preprocessed_output, int *ecg_complex_length) {
    double sampling_rate = 360.00;

    /* High pass filter with 0.5Hz cutoff */
    BWHighPass* hp_filter = create_bw_high_pass_filter(4, (int)sampling_rate, 0.5);
    for (int i = 0; i < rows; i++) {
        preprocessed_output[i][0] = bw_high_pass(hp_filter, data[i][0]);
#if 1
        preprocessed_output[i][1] = data[i][1];
        preprocessed_output[i][2] = data[i][2];
        preprocessed_output[i][3] = data[i][3];
        preprocessed_output[i][4] = data[i][4];
        preprocessed_output[i][5] = data[i][5];
        preprocessed_output[i][6] = data[i][6];
        preprocessed_output[i][7] = data[i][7];
        preprocessed_output[i][8] = data[i][8];
        preprocessed_output[i][9] = data[i][9];
        preprocessed_output[i][10] = data[i][10];
        preprocessed_output[i][11] = data[i][11];
#endif

#if 0 // Disable the above and use this when real inputs used for rest of the ECG leads
        preprocessed_output[i][1] = bw_high_pass(hp_filter, data[i][1]);
        preprocessed_output[i][2] = bw_high_pass(hp_filter, data[i][2]);
        preprocessed_output[i][3] = bw_high_pass(hp_filter, data[i][3]);
        preprocessed_output[i][4] = bw_high_pass(hp_filter, data[i][4]);
        preprocessed_output[i][5] = bw_high_pass(hp_filter, data[i][5]);
        preprocessed_output[i][6] = bw_high_pass(hp_filter, data[i][6]);
        preprocessed_output[i][7] = bw_high_pass(hp_filter, data[i][7]);
        preprocessed_output[i][8] = bw_high_pass(hp_filter, data[i][8]);
        preprocessed_output[i][9] = bw_high_pass(hp_filter, data[i][9]);
        preprocessed_output[i][10] = bw_high_pass(hp_filter, data[i][10]);
        preprocessed_output[i][11] = bw_high_pass(hp_filter, data[i][11]);
#endif
    }
    free_bw_high_pass(hp_filter);

    /* Low pass filter with 0.5Hz cutoff */
    BWLowPass* lp_filter = create_bw_low_pass_filter(4, (int)sampling_rate, 90);
    for (int i = 0; i < rows; i++) {
        preprocessed_output[i][0] = bw_low_pass(lp_filter, preprocessed_output[i][0]);

#if 0 // Enable when real inputs used for rest of the ECG leads
        preprocessed_output[i][1] = bw_low_pass(lp_filter, preprocessed_output[i][1]);
        preprocessed_output[i][2] = bw_low_pass(lp_filter, preprocessed_output[i][2]);
        preprocessed_output[i][3] = bw_low_pass(lp_filter, preprocessed_output[i][3]);
        preprocessed_output[i][4] = bw_low_pass(lp_filter, preprocessed_output[i][4]);
        preprocessed_output[i][5] = bw_low_pass(lp_filter, preprocessed_output[i][5]);
        preprocessed_output[i][6] = bw_low_pass(lp_filter, preprocessed_output[i][6]);
        preprocessed_output[i][7] = bw_low_pass(lp_filter, preprocessed_output[i][7]);
        preprocessed_output[i][8] = bw_low_pass(lp_filter, preprocessed_output[i][8]);
        preprocessed_output[i][9] = bw_low_pass(lp_filter, preprocessed_output[i][9]);
        preprocessed_output[i][10] = bw_low_pass(lp_filter, preprocessed_output[i][10]);
        preprocessed_output[i][11] = bw_low_pass(lp_filter, preprocessed_output[i][11]);
#endif
    }
    free_bw_low_pass(lp_filter);

    /* Notch filter centered around 50Hz for powerline/baseline wander */
    BWBandStop* bs_filter = create_bw_band_stop_filter(4, (int)sampling_rate, 49, 71);
    for (int i = 0; i < rows; i++) {
        preprocessed_output[i][0] = bw_band_stop(bs_filter, preprocessed_output[i][0]);

#if 0 // Enable when real inputs used for rest of the ECG leads
        preprocessed_output[i][1] = bw_band_stop(bs_filter, preprocessed_output[i][1]);
        preprocessed_output[i][2] = bw_band_stop(bs_filter, preprocessed_output[i][2]);
        preprocessed_output[i][3] = bw_band_stop(bs_filter, preprocessed_output[i][3]);
        preprocessed_output[i][4] = bw_band_stop(bs_filter, preprocessed_output[i][4]);
        preprocessed_output[i][5] = bw_band_stop(bs_filter, preprocessed_output[i][5]);
        preprocessed_output[i][6] = bw_band_stop(bs_filter, preprocessed_output[i][6]);
        preprocessed_output[i][7] = bw_band_stop(bs_filter, preprocessed_output[i][7]);
        preprocessed_output[i][8] = bw_band_stop(bs_filter, preprocessed_output[i][8]);
        preprocessed_output[i][9] = bw_band_stop(bs_filter, preprocessed_output[i][9]);
        preprocessed_output[i][10] = bw_band_stop(bs_filter, preprocessed_output[i][10]);
        preprocessed_output[i][11] = bw_band_stop(bs_filter, preprocessed_output[i][11]);
#endif
    }
    free_bw_band_stop(bs_filter);

    /** QRS length detection for finding the value of T. (kosachevds)
     *  In this paper, demarcation of ECG complex was defined as per 120% of the mean interval between consecutive
     *  consecutive heartbeats found using QRS complex detection. While the paper has used PCA for QRS detection,
     *  we found that the Pan-Tompkins algorithm is also sufficiently accurate in detecting the R waves with an error
     *  of about 2%. For this project, the length between all R peaks is averaged and a 2% error is factored into it
     *  for allowing overlap between ECG signals.
     */
    char result[rows]; // Dummy as we do not need this.
    double *in = (double *) malloc(sizeof(double[rows]));
    for (int i = 0; i < rows; i++) {
        in[i] = preprocessed_output[i][0];
    }
    *ecg_complex_length = DetectQrsPeaks(in, rows, result, sampling_rate);
    *ecg_complex_length -= 0.02 * (*ecg_complex_length);

    return 0;
}

/** Kalman Filtering Stage:
 *  Adaptive Kalman Filter operations based on Vullings, Vries & Bergmans' equations begin from here.
 *  Core implentation logic and mathematical equations used here are directly based of the flow sequence
 *  described by Vullings, Vries & Bergmans and mentioned in the supported project document. Instead of
 *  using element-wise operations, all equations are directly operated on in their matrix form and
 *  outputs are handled as scalars internally if needed.
 */
int process_kalman(double **data, int rows, int cols, int ecg_complex_length, double *output) {
    const int T = ecg_complex_length;

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
    int loop_count = ((rows + T - rem) / T) - 1;

    // Initial state estimate (x_hat) and state covariance (P/psi) by taking average of the first N ECG complexes.
    double x0[T][N];
    for (int row = 0; row < T; row++) {
        for (int col = 0; col < N; col++) {
            x0[row][col] = (double) data[row + (col * T)][0];
        }
    }
    mean((double *)x0, (double *)x_hat, T, N, 1);

    double P_mat[1][1];
    covariance((double *)x_hat, (double *)P_mat, T, 1);
    P = P_mat[0][0];

    double y[T][1];

    for (int count = 0; count < loop_count; count++) {
       /* Measurement matrix (ŵ) and Measurement covariance (R/Σ) */
       /** EQ4: Coefficients of estimate γ̂ = [Yt * Y]^−1 * Yt * y(i)
        *           where t denotes tranpose.
        *  Subscript -i ommitted for clarity.
        */
       double y_inv[ECG_LEADS][ECG_LEADS], Y_transpose[ECG_LEADS][T], Y_res[ECG_LEADS][ECG_LEADS], Y_res2[1][T], y_coeff_hat[ECG_LEADS][1];
       for (int row = 0; row < T; row++) {
           // Update y
           y[row][0] = data[row + (count * T)][0];

           // Update the measurement matrix
           Y[row][0] = data[row + (count * T)][1];
           Y[row][1] = data[row + (count * T)][2];
           Y[row][2] = data[row + (count * T)][3];
           Y[row][3] = data[row + (count * T)][4];
           Y[row][4] = data[row + (count * T)][5];
           Y[row][5] = data[row + (count * T)][6];
           Y[row][6] = data[row + (count * T)][7];
           Y[row][7] = data[row + (count * T)][8];
           Y[row][8] = data[row + (count * T)][9];
           Y[row][9] = data[row + (count * T)][10];
           Y[row][10] = data[row + (count * T)][11];
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
           output[i + (count * T)] = x_hat[i][0];
       }
   }

   return 0;
}
