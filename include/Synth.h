#pragma once

#include "Theory.h"
#include <cstddef>
#include <vector>

struct Envelope {
  double attack_ms;  // ms
  double decay_ms;   // ms
  double sustain_level; // level (0.0 to 1.0)
  double release_ms; // ms
};

class Synth {
public:
  std::vector<double> audio;
  int samplerate;
  Envelope env;

  Synth(int samplerate, Envelope env);

  void addNote(Note note, std::size_t start, std::size_t end);

  void create_wav(const char *file_name);

  void envToSamples(Envelope env);

private:
  int attack_smp;
  int decay_smp;
  int release_smp;

  void normaliseAudio();

  double getFrequency(Note note);

  std::vector<double> genSineWave(Note note, std::size_t start, std::size_t end);

  int msToSamples(double ms);
};
