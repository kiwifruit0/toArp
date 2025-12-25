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
    : audio(audio), samplerate(samplerate) {
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
