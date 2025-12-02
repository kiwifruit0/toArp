#include "find_peaks.hpp"
#include <iostream>
#include <sndfile.h>
#include <vector>

void findpeaks(std::vector<double> signal);

int main() {
  SF_INFO sfinfo = {};
  SNDFILE *sf = sf_open("assets/ball_bounce.wav", SFM_READ, &sfinfo);
  if (!sf) {
    std::cerr << "Error: " << sf_strerror(NULL);
    return 1;
  }

  if (sfinfo.frames / sfinfo.channels > 158760000) {
    std::cerr << "file too large, exiting";
    return 1;
  }

  // read audio to vector
  std::vector<double> samples(sfinfo.frames * sfinfo.channels);
  sf_readf_double(sf, samples.data(), sfinfo.frames);
  sf_close(sf);

  // convert to mono signal
  std::vector<double> mono;
  if (sfinfo.channels == 2) {
    mono.resize(sfinfo.frames);
    for (sf_count_t i = 0; i < sfinfo.frames; i++) {
      mono[i] = 0.5f * (samples[i * 2] + samples[i * 2 + 1]);
    }
  } else {
    mono = samples;
  }
  std::cout << sfinfo.frames << "\n";
  std::cout << sfinfo.samplerate;

  for (size_t i = 0; i < mono.size(); i++) {
    if (mono[i] > 0.5) {
      std::cout << mono[i] << " ";
    }
  }

  findpeaks(mono);

  return 0;
}

void findpeaks(std::vector<double> signal) {

  findPeaks::PeakConditions conditions;
  conditions.set_height(0);      // Minimum height of 2.0
  conditions.set_prominence(0.1); // Minimum prominence of 1.0
  conditions.set_distance(2);      // At least 2 samples between peaks
  conditions.set_width(1.0, 3.0);  // Width between 1.0 and 3.0 samples
  conditions.set_rel_height(0.3);  // Measure width at 70% of peak height

  std::vector<findPeaks::peak_result_t> peaks =
      findPeaks::find_peaks(signal, conditions);

  for (const auto &peak : peaks) {
    std::cout << "Peak at position " << peak.peak << ":\n";
    std::cout << "  Height: " << peak.peak_height << "\n";
    std::cout << "  Prominence: " << peak.prominence.prominence << "\n";
    std::cout << "  Width: " << peak.width.width << " at height "
              << peak.width.width_height << "\n";
    std::cout << "  Left threshold: " << peak.threshold.left_threshold << "\n";
    std::cout << "  Right threshold: " << peak.threshold.right_threshold
              << "\n";

    // Plateau information
    if (peak.plateau.plateau_size > 1) {
      std::cout << "  Plateau size: " << peak.plateau.plateau_size << "\n";
    }
  }
}
