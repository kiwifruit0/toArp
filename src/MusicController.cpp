#include "MusicController.h"

MusicController::MusicController(Scale scale, int samplerate)
    : scale(scale), samplerate(samplerate), synth(samplerate), last_note(), movingUp() {}

void MusicController::addNextNote(Notes note, std::size_t start, std::size_t end) { }
