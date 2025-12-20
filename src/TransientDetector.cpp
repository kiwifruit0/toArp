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
  this->params.window_size =
      (params.window_size > 0) ? params.window_size : defaults.window_size;
  this->params.hop_size =
      (params.hop_size > 0) ? params.hop_size : defaults.hop_size;
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

std::vector<std::vector<std::complex<double>>> TransientDetector::computeSTFT() {
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
