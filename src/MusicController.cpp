#include "MusicController.h"

MusicController::MusicController(Scale scale, int samplerate)
    : scale(scale), samplerate(samplerate),
      intervals(getScaleIntervals(scale.mode)), last_note(intervals.back()),
      movingUp(true), synth(samplerate) {}

void MusicController::addNextNote(std::size_t start, std::size_t end) {}
