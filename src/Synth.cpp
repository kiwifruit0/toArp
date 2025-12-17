#include "Synth.h"
#include <cmath>
#include <sndfile.h>

Synth::Synth(Scale scale, int samplerate)
    : scale(scale), samplerate(samplerate) {}

void Synth::addNote(struct Note note, std::size_t start, std::size_t end) {
  // if end sample is beyond current audio size, resize audio
  if (end >= audio.size()) {
    this->audio.resize(end, 0.0);
  }
  std::vector<double> sine = genSineWave(note, end - start);
  // add beep to audio at specified position
  for (std::size_t i = start; i < end; i++) {
    this->audio[i] += sine[i - start];
  }
}

void Synth::create_wav() {
  // setting up file info
  SF_INFO sfinfo = {};
  sfinfo.format = SF_FORMAT_WAV | SF_FORMAT_PCM_16;
  sfinfo.samplerate = this->samplerate;
  sfinfo.channels = 1;

  // normalise audio before writing
  normaliseAudio();

  // create and write to file
  SNDFILE *outfile = sf_open("output.wav", SFM_WRITE, &sfinfo);
  sf_write_double(outfile, this->audio.data(), this->audio.size());
  sf_close(outfile);
}

void Synth::normaliseAudio() {
  double max_amp = 0.0;
  for (double s : audio)
    max_amp = std::max(max_amp, std::abs(s));

  if (max_amp > 1.0) {
    for (double &s : audio)
      s /= max_amp;
  }
}

double Synth::getFrequency(Note note) {
  int n = static_cast<int>(note.note) + (note.octave - 4) * 12;
  return 440.0 * pow(2.0, (n - 9) / 12.0);
}

std::vector<double> Synth::genSineWave(Note note, std::size_t num_samples) {
  // initialise beep vector
  std::vector<double> wave(num_samples, 0.0);
  double freq = getFrequency(note);
  // find phase increment
  const double phase_incr = (2 * M_PI * freq) / this->samplerate;
  // fill beep vector with sine wave samples
  for (std::size_t i = 0; i < wave.size(); i++) {
    wave[i] = sin(phase_incr * i);
  }
  return wave;
}
