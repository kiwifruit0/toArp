#pragma once

#include <vector>

struct Parameters {
  int window_size = 0;
  int hop_size = 0;
  int nu = 3;

  static Parameters getDefaults(int samplerate) {
    Parameters p;
    p.window_size = static_cast<int>(0.040 * samplerate);
    p.hop_size = static_cast<int>(0.010 * samplerate);

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
