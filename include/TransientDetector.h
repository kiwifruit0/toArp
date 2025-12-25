#pragma once

#include <vector>

struct Parameters {
  // stft parameters
  int window_size = 0;
  int hop_size = 0;

  // algorithm parameters from paper
  int nu = 3;          // spectral spread/vertical neighbours (eq. 5)
  double beta = 2.0;   // threshold strength (eq. 7)
  int tau = 3;         // temporal context for threshold (eq. 7)
  double delta = 0.1;  // extraction fraction (eq. 10, 11)
  int iterations = 20; // number of iterations
  double lambda_thr_fraction = 1.0 / 6.0; // fraction of N (eq. 8-10)

  static Parameters getDefaults(int samplerate) {
    Parameters p;
    // paper uses 40ms window with 30ms overlap
    // this gives 40ms frame and 10ms hop
    p.window_size = static_cast<int>(0.040 * samplerate);
    p.hop_size = static_cast<int>(0.010 * samplerate);

    // paper's optimised values (table 2)
    p.nu = 3;
    p.beta = 2.0;
    p.tau = 3;
    p.delta = 0.1;
    p.iterations = 20;
    p.lambda_thr_fraction = 1.0 / 6.0; // N/6 from paper

    return p;
  }
};

class TransientDetector {
public:
  TransientDetector(std::vector<double> audio, int samplerate,
                    Parameters params);

private:
  std::vector<double> audio;
  int samplerate;
  Parameters params;

  std::vector<double> blackmanHarrisWindow(std::size_t n);

  std::vector<std::vector<std::complex<double>>> computeSTFT();

  std::vector<std::vector<double>>
  computeTMinus(const std::vector<std::vector<std::complex<double>>> &X);

  std::vector<std::vector<double>>
  computeTPlus(const std::vector<std::vector<std::complex<double>>> &X);

  std::vector<std::vector<double>>
  computeF(const std::vector<std::vector<double>> &T_minus,
           const std::vector<std::vector<double>> &T_plus);
};
