#include "TransientDetector.h"
#include <cmath>
#include <complex>
#include <cstddef>
#include <vector>

/*
implements transient detection algorithm based on the paper:
"A Transient Detection Algorithm for Audio Using Iterative Analysis of STFT"
by Balaji Thoshkahna, Francois Xavier Nsabimana and  K.R.Ramakrishnan
*/

TransientDetector::TransientDetector(std::vector<double> audio, int samplerate,
                                     Parameters params)
    : audio(audio),
      samplerate(samplerate) {
  Parameters defaults = Parameters::getDefaults(samplerate);

  // Apply defaults for any unset parameters
  this->params.window_size =
      (params.window_size > 0) ? params.window_size : defaults.window_size;
  this->params.hop_size =
      (params.hop_size > 0) ? params.hop_size : defaults.hop_size;
  this->params.nu = (params.nu >= 0) ? params.nu : defaults.nu;
  this->params.beta = (params.beta > 0) ? params.beta : defaults.beta;
  this->params.tau = (params.tau >= 0) ? params.tau : defaults.tau;
  this->params.delta = (params.delta > 0) ? params.delta : defaults.delta;
  this->params.iterations =
      (params.iterations > 0) ? params.iterations : defaults.iterations;
  this->params.lambda_thr_fraction = (params.lambda_thr_fraction > 0)
                                         ? params.lambda_thr_fraction
                                         : defaults.lambda_thr_fraction;
}

std::vector<double> TransientDetector::blackmanHarrisWindow(std::size_t N) {
  std::vector<double> window(N);
  const double a0 = 0.35875;
  const double a1 = 0.48829;
  const double a2 = 0.14128;
  const double a3 = 0.01168;
  for (int n = 0; n < N; n++) {
    double w = a0 - a1 * cos((2 * M_PI * n) / (N - 1)) +
               a2 * cos((4 * M_PI * n) / (N - 1)) -
               a3 * cos((6 * M_PI * n) / (N - 1));
  }
  return window;
}

std::vector<std::vector<std::complex<double>>>
TransientDetector::computeSTFT() {
  // equation (2) from paper
  std::vector<double> w = blackmanHarrisWindow(params.window_size);

  // calculate number of frames
  // stop when cant fit a full frame
  int num_frames = (audio.size() - params.window_size) / params.hop_size + 1;

  // spectrogram container: X[frame_index][freq_bin]
  std::vector<std::vector<std::complex<double>>> spectrogram;
  spectrogram.reserve(num_frames);

  // loop over frames
  for (int i = 0; i < num_frames; ++i) {

    int start_sample = i * params.hop_size;
    std::vector<std::complex<double>> frame_spectrum;

    // compute frequncy bins to N/2+1
    int num_bins = params.window_size / 2 + 1;
    frame_spectrum.reserve(num_bins);

    // loop over frequency bins
    for (int k = 0; k < num_bins; ++k) {
      std::complex<double> sum(0.0, 0.0);

      for (int n = 0; n < params.window_size; ++n) {

        double x_n = audio[start_sample + n];
        double windowed_signal = x_n * w[n];

        double angle = (2.0 * M_PI * n * k) / params.window_size;
        std::complex<double> basis(cos(angle), -sin(angle));

        sum += windowed_signal * basis;
      }

      frame_spectrum.push_back(sum);
    }
    spectrogram.push_back(frame_spectrum);
  }

  return spectrogram;
}

std::vector<std::vector<double>> TransientDetector::computeTMinus(
    const std::vector<std::vector<std::complex<double>>> &X) {
  // equation (3) from paper

  int num_frames = X.size();
  if (num_frames == 0)
    return {};

  int num_bins = X[0].size();

  // initialise T- with same dimensions as X
  std::vector<std::vector<double>> T_minus(num_frames,
                                           std::vector<double>(num_bins, 0.0));

  // start from frame i=1
  for (int i = 1; i < num_frames; ++i) {
    for (int k = 0; k < num_bins; ++k) {
      double mag_current = std::abs(X[i][k]);
      double mag_previous = std::abs(X[i - 1][k]);
      T_minus[i][k] = mag_current - mag_previous;
    }
  }

  return T_minus;
}
std::vector<std::vector<double>> TransientDetector::computeTPlus(
    const std::vector<std::vector<std::complex<double>>> &X) {
  // equation (4) from paper

  int num_frames = X.size();
  if (num_frames == 0)
    return {};

  int num_bins = X[0].size();

  // initialise T- with same dimensions as X
  std::vector<std::vector<double>> T_minus(num_frames,
                                           std::vector<double>(num_bins, 0.0));

  // end at frame size-1
  for (int i = 0; i < num_frames - 1; ++i) {
    for (int k = 0; k < num_bins; ++k) {
      double mag_current = std::abs(X[i][k]);
      double mag_previous = std::abs(X[i + 1][k]);
      T_minus[i][k] = mag_current - mag_previous;
    }
  }

  return T_minus;
}

std::vector<std::vector<double>>
TransientDetector::computeF(const std::vector<std::vector<double>> &T_minus,
                            const std::vector<std::vector<double>> &T_plus) {
  // equation (5) from paper

  int num_frames = T_minus.size();
  if (num_frames == 0)
    return {};

  int num_bins = T_minus[0].size();

  std::vector<std::vector<double>> F(num_frames,
                                     std::vector<double>(num_bins, 0.0));

  for (int i = 0; i < num_frames; ++i) {
    for (int j = 0; j < num_bins; ++j) {
      double sum = 0.0;

      // sum over vertical neighbours: [j-ν, j+ν]
      int k_start = std::max(0, j - params.nu);
      int k_end = std::min(num_bins - 1, j + params.nu);

      for (int k = k_start; k <= k_end; ++k) {
        // half-wave rectification: (1 + sgn(x)) * x = x if x>=0, else 0
        double t_minus_contrib = (T_minus[i][k] >= 0) ? T_minus[i][k] : 0.0;
        double t_plus_contrib = (T_plus[i][k] >= 0) ? T_plus[i][k] : 0.0;

        sum += t_minus_contrib + t_plus_contrib;
      }

      F[i][j] = 0.5 * sum;
    }
  }

  return F;
}

// Equation (7): Compute dynamic threshold λ(i,j)
// λ(i,j) = β × mean(F(l,j)) for l in [i-τ, i+τ]
std::vector<std::vector<double>>
TransientDetector::computeLambda(const std::vector<std::vector<double>> &F) {

  int num_frames = F.size();
  if (num_frames == 0)
    return {};

  int num_bins = F[0].size();

  std::vector<std::vector<double>> lambda(num_frames,
                                          std::vector<double>(num_bins, 0.0));

  for (int i = 0; i < num_frames; ++i) {
    for (int j = 0; j < num_bins; ++j) {
      // Compute temporal neighbourhood [i-τ, i+τ]
      int l_start = std::max(0, i - params.tau);
      int l_end = std::min(num_frames - 1, i + params.tau);

      // Compute mean of F in the temporal window
      double sum = 0.0;
      int count = 0;
      for (int l = l_start; l <= l_end; ++l) {
        sum += F[l][j];
        count++;
      }

      double mean = (count > 0) ? sum / count : 0.0;
      lambda[i][j] = params.beta * mean;
    }
  }

  return lambda;
}

// Equation (8): Compute flag function Γ(i,j)
// Γ(i,j) = 1 if F(i,j) > λ(i,j), else 0
std::vector<std::vector<int>> TransientDetector::computeGamma(
    const std::vector<std::vector<double>> &F,
    const std::vector<std::vector<double>> &lambda) {

  int num_frames = F.size();
  if (num_frames == 0)
    return {};

  int num_bins = F[0].size();

  std::vector<std::vector<int>> Gamma(num_frames,
                                      std::vector<int>(num_bins, 0));

  for (int i = 0; i < num_frames; ++i) {
    for (int j = 0; j < num_bins; ++j) {
      Gamma[i][j] = (F[i][j] > lambda[i][j]) ? 1 : 0;
    }
  }

  return Gamma;
}

std::vector<int> TransientDetector::computeSigmaGamma(
    const std::vector<std::vector<int>> &Gamma) {

  // Equation (9): Compute Σ_Γ(i) - sum of flags across frequency bins

  int num_frames = Gamma.size();
  if (num_frames == 0)
    return {};

  int num_bins = Gamma[0].size();

  std::vector<int> Sigma_Gamma(num_frames, 0);

  for (int i = 0; i < num_frames; ++i) {
    int sum = 0;
    for (int j = 0; j < num_bins; ++j) {
      sum += Gamma[i][j];
    }
    Sigma_Gamma[i] = sum;
  }

  return Sigma_Gamma;
}

void TransientDetector::updateTransientAndMagnitude(
    std::vector<std::vector<std::complex<double>>> &X,
    std::vector<std::vector<double>> &P, const std::vector<int> &Sigma_Gamma,
    double lambda_thr) {
  // Equations (10-11): Update transient signal P and magnitude spectrum X

  int num_frames = X.size();
  if (num_frames == 0)
    return;

  int num_bins = X[0].size();

  for (int i = 0; i < num_frames; ++i) {
    // Check if this frame has a transient
    if (Sigma_Gamma[i] >= lambda_thr) {
      // This frame is declared as transient
      for (int j = 0; j < num_bins; ++j) {
        double mag = std::abs(X[i][j]);

        // Equation (10): Add fraction delta to transient signal
        P[i][j] += params.delta * mag;

        // Equation (11): Reduce magnitude in X by fraction delta
        X[i][j] *= (1.0 - params.delta);
      }
    }
  }
}

std::vector<double> TransientDetector::inverseSTFT(
    const std::vector<std::vector<std::complex<double>>> &P_complex,
    const std::vector<std::vector<std::complex<double>>> &phase_reference) {

  // Inverse STFT to convert spectrogram back to time domain
  int num_frames = P_complex.size();
  if (num_frames == 0)
    return {};

  int num_bins = P_complex[0].size();
  int N = params.window_size;

  // Output signal length
  int signal_length = (num_frames - 1) * params.hop_size + N;
  std::vector<double> output(signal_length, 0.0);
  std::vector<double> window_sum(signal_length, 0.0);

  std::vector<double> w = blackmanHarrisWindow(N);

  for (int i = 0; i < num_frames; ++i) {
    // Reconstruct full spectrum (including negative frequencies)
    std::vector<std::complex<double>> full_spectrum(N);

    // Positive frequencies (including DC and Nyquist)
    for (int k = 0; k < num_bins; ++k) {
      // Use magnitude from P and phase from original signal
      double mag = std::abs(P_complex[i][k]);
      double phase = std::arg(phase_reference[i][k]);
      full_spectrum[k] = std::polar(mag, phase);
    }

    // Negative frequencies (conjugate symmetry for real signals)
    for (int k = num_bins; k < N; ++k) {
      int mirror_k = N - k;
      full_spectrum[k] = std::conj(full_spectrum[mirror_k]);
    }

    // IDFT
    std::vector<double> frame(N, 0.0);
    for (int n = 0; n < N; ++n) {
      std::complex<double> sum(0.0, 0.0);
      for (int k = 0; k < N; ++k) {
        double angle = (2.0 * M_PI * n * k) / N;
        std::complex<double> basis(cos(angle), sin(angle));
        sum += full_spectrum[k] * basis;
      }
      frame[n] = sum.real() / N;
    }

    // Overlap-add with windowing
    int start_sample = i * params.hop_size;
    for (int n = 0; n < N; ++n) {
      output[start_sample + n] += frame[n] * w[n];
      window_sum[start_sample + n] += w[n] * w[n];
    }
  }

  // Normalise by window overlap
  for (int n = 0; n < signal_length; ++n) {
    if (window_sum[n] > 1e-10) {
      output[n] /= window_sum[n];
    }
  }

  return output;
}

std::vector<double> TransientDetector::detectTransients() {

  // algirthm 1 from paper
  // compute initial stft
  std::vector<std::vector<std::complex<double>>> X = computeSTFT();
  std::vector<std::vector<std::complex<double>>> X_original =
      X; // Keep for phase

  int num_frames = X.size();
  int num_bins = X[0].size();

  // initialise p to zeros (equation above algorithm 1)
  std::vector<std::vector<double>> P(num_frames,
                                     std::vector<double>(num_bins, 0.0));

  // calculate λ_thr from window size
  double lambda_thr = params.lambda_thr_fraction * params.window_size;

  // iterative algorithm (algorithm 1)
  for (int iteration = 0; iteration < params.iterations; ++iteration) {
    // step 3: compute detection functions from current x
    auto T_minus = computeTMinus(X);
    auto T_plus = computeTPlus(X);
    auto F = computeF(T_minus, T_plus);

    // Step 1: Compute dynamic thresholds
    auto lambda = computeLambda(F);
    auto Gamma = computeGamma(F, lambda);
    auto Sigma_Gamma = computeSigmaGamma(Gamma);

    // Step 2: Extract transients and update X
    updateTransientAndMagnitude(X, P, Sigma_Gamma, lambda_thr);
  }

  // convert p magnitude to complex with original phase
  std::vector<std::vector<std::complex<double>>> P_complex(num_frames);
  for (int i = 0; i < num_frames; ++i) {
    P_complex[i].resize(num_bins);
    for (int j = 0; j < num_bins; ++j) {
      double phase = std::arg(X_original[i][j]);
      P_complex[i][j] = std::polar(P[i][j], phase);
    }
  }

  // convert back to time domain using inverse stft
  std::vector<double> transient_signal = inverseSTFT(P_complex, X_original);

  return transient_signal;
}
