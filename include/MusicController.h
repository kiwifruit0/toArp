#pragma once

#include "Synth.h"
#include "Theory.h"
#include <cstddef>

class MusicController {
public:
    const Scale scale;
    const int samplerate;
    const ScaleIntervals intervals;
    int last_note;
    bool movingUp;
    Synth synth;
  MusicController(Scale scale, int samplerate);

  void addNextNote(std::size_t start, std::size_t end);
};
