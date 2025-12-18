#pragma once

#include "Synth.h"
#include "Theory.h"
#include <cstddef>

class MusicController {
public:
  const Scale scale;
  const int samplerate;
  const ScaleIntervals intervals;
  const int num_octaves;
  const int start_octave;

  int position;
  bool moving_up;
  Synth synth;

  MusicController(Scale scale, int samplerate, int start_octave=3, int num_octaves = 3);
  void addNextNote(std::size_t start, std::size_t end);
  
  void write();

private:
  int getTotalNotes() const;
  Note positionToNote(int pos) const;
};
