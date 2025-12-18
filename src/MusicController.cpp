#include "MusicController.h"

MusicController::MusicController(Scale scale, int samplerate, int start_octave,
                                 int num_octaves)
    : scale(scale),
      samplerate(samplerate),
      intervals(getScaleIntervals(scale.mode)),
      position(0),
      start_octave(start_octave),
      num_octaves(num_octaves),
      moving_up(true),
      synth(samplerate) {}

int MusicController::getTotalNotes() const {
  return intervals.size() * num_octaves;
}

void MusicController::addNextNote(std::size_t start, std::size_t end) {
  // add note to audio
  Note note = positionToNote(this->position);
  synth.addNote(note, start, end);

  // increment position
  int total = getTotalNotes();

  if (moving_up) {
    position++;
    if (position >= total) {
      moving_up = false;
      // reverse without repeating top note
      position = total - 2;
    }
  } else {
    position--;
    if (position < 0) {
      moving_up = true;
      // reverse without repeating bottom note
      position = 1;
    }
  }
}

void MusicController::write() { synth.create_wav("output.wav"); }

Note MusicController::positionToNote(int pos) const {
  // find note in position in scale and octave
  int scale_degree = pos % intervals.size();
  int octave_offset = pos / intervals.size();

  int semitone_in_scale = intervals[scale_degree];

  // calculate total semitones from the root note
  int total_semitones =
      static_cast<int>(scale.root) + semitone_in_scale + (12 * octave_offset);

  // convert to note and octave
  // assumes base octave of 4
  Notes note = static_cast<Notes>(total_semitones % 12);
  int octave = start_octave + (total_semitones / 12);

  return Note{note, octave};
}
