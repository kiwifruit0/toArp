#pragma once

#include "Synth.h"
#include "Theory.h"
#include <cstddef>

class MusicController {
public:
  Scale scale;
  Notes last_note;
  bool movingUp = true;
  int samplerate;
  Synth synth;

  MusicController(Scale scale, int samplerate);

  void addNextNote(Notes note, std::size_t start, std::size_t end);
};
