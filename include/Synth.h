#pragma once

#include "Theory.h"
#include <cstddef>
#include <vector>

class Synth {
public:
  std::vector<double> audio;
  int samplerate;

  Synth(int samplerate);

  void addNote(Note note, std::size_t start, std::size_t end);

  void create_wav();

private:
  void normaliseAudio();

  double getFrequency(Note note);

  std::vector<double> genSineWave(Note note, std::size_t num_samples);
};
