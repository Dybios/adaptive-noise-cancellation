# Adaptive Kalman Filter for Noise Cancellation

This is the implementation of the Adaptive Kalman filter as proposed by Vullings, R. Vullings, B. de Vries and J. W. M. Bergmans (https://ieeexplore.ieee.org/document/5667049). This implementation utilizes Bayesian estimation techniques to estimate next state using at least four, simulataneous recorded inputs. The ECG data is first preprocessed with lowpass filter of 0.5Hz, highpass filter of 90Hz and a bandpass filter centered around 60Hz (measurements taken in America). After preprocessing, the data is run through a QRS detector using the Pan-Tompkins algorithm to find the length of QRS complexes. Following which, a residual analysis is performed to obtain the updated process noise covariance from the recorded input signal. We then iteratively update the Kalman Gain to prefer the estimated signal value instead of the recorded signal, in case the residue (difference between recorded and estimated signal) diverges too much.

LP, HP and BP filter operation logic by adis300: https://github.com/adis300/filter-c

QRS complex detection using Pan-Tompkin's algorithm by kosachevds: https://github.com/kosachevds/qrs_detector

## Adaptive Noise Supression (Input Noise Filtering)
The adaptive nature of this filter also allows for different applications such as adaptive noise cancellation from input audio signal. The sinewave example in MATLAB has been repurposed (adaptive_kalman_audio) to use a fixed buffer of 10k instead of length of each sine period, which then tries to estimate the voice signal from the recorded voice + noise input. The estimated voice signal is then subtracted from the recorded to obtain the noise vector. This noise vector can then be used to remove the noise from the recorded signal.
