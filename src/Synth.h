#pragma once

#include <vector>
#include <cstddef>

enum class Notes { C, Cs, D, Ds, E, F, Fs, G, Gs, A, As, B };

struct Note {
  Notes note;
  int octave;
};

enum class Mode { Major, Minor };

struct Scale {
  Notes root;
  Mode mode;
};

class Synth {
public:
  Scale scale;
  std::vector<double> audio;
  int samplerate;

  Synth(Scale scale, int samplerate);

  void addNote(struct Note note, std::size_t start, std::size_t end);

  void create_wav();

private:
  void normaliseAudio();

  double getFrequency(Note note);

  std::vector<double> genSineWave(Note note, std::size_t num_samples);
};
