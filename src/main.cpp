#include <iostream>
#include <sndfile.h>
#include <vector>
#include <cmath>

enum class Note { C, Cs, D, Ds, E, F, Fs, G, Gs, A, As, B };

enum class Mode { Major, Minor };

struct Scale {
  Mode mode;
  Note note;
};

class Synth {
public:
  Scale scale;
  int samplerate;
  Synth(Scale scale, int samplerate) {
    this->scale = scale;
    this->samplerate = samplerate;
  }
  double getFrequency(Note note, int octave) {
    int n = static_cast<int>(note) + (octave - 4) * 12;
    return 440.0 * pow(2.0, (n - 9) / 12.0);
  }
};

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

  // create synth
  Scale scale = {Mode::Major, Note::C};
  Synth synth(scale, sfinfo.samplerate);
  double freq = synth.getFrequency(Note::C, 4);

  std::cout << "frequency of note C4: " << freq;


  return 0;
}
