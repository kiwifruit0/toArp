#include "MusicController.h"
#include "Theory.h"
#include <iostream>
#include <sndfile.h>
#include <vector>

int main() {
  SF_INFO sfinfo = {};
  SNDFILE *sf = sf_open("../assets/ball_bounce.wav", SFM_READ, &sfinfo);
  if (!sf) {
    std::cerr << "Error: " << sf_strerror(NULL);
    throw std::runtime_error("Failed to open audio file");
  }

  // limit file size to ~30 minutes at 44.1kHz
  if (sfinfo.frames / sfinfo.channels > 158760000) {
    std::cerr << "file too large, exiting";
    throw std::runtime_error("File too large");
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

  // music controller test
  Scale scale = {Notes::C, Mode::Major};

  MusicController mc(scale, sfinfo.samplerate);
  int notelength = 4000;
  for (int i = 0; i < 40; i++) {
    mc.addNextNote(notelength * i, notelength * i + notelength);
  }

  mc.write();

  return 0;
}
